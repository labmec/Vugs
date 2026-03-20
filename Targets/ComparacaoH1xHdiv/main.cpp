#include <iostream>
#include <filesystem>
#include <math.h>
#include "pzcmesh.h"
#include "TPZElementMatrixT.h"
#include "pzmanvector.h"
#include "TPZGeoMeshTools.h"
#include "TPZCompMeshTools.h"
#include "TPZRefPattern.h"
#include "TPZGenGrid2D.h"
#include "TPZVTKGeoMesh.h"
#include "pzvec.h"
#include "TPZVTKGenerator.h"
#include <fstream>
#include "TPZMaterial.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFractureFlow.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZNullMaterialCS.h"
#include "TPZNullMaterial.h"
#include "TPZAnalysis.h"
#include "pzstepsolver.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"
#include <set>
#include "TPZAnalyticSolution.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TSFMixedDarcy.h"
#include "TPZHDivApproxCreator.h"

TPZCompMesh* HdivMesh(TPZGeoMesh *);
TPZCompMesh* Pressuremesh(TPZGeoMesh *, int order);
void GetAtomicIds(TPZGeoMesh *geomesh, std::set<int> &volId, std::set<int> &bcId);
void insertAtomicMaterials(TPZCompMesh *cmesh, std::set<int> matIdsVol, std::set<int> matIdsBcs);

TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void findElDim(TPZStack<TPZGeoElSide> &allneigh, int dim, TPZStack<TPZGeoElSide> &allneighdim);

int H1Vugs();
int H1Fracs();
int Hdiv_Fract();
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine){
            
    // Creating gmsh reader
    TPZGmshReader  GeometryFine;
    TPZGeoMesh *gmeshFine;
    REAL l = 1.0;
    GeometryFine.SetCharacteristiclength(l);
    
    // Reading mesh
    GeometryFine.SetDimNamePhysical(dim_name_and_physical_tagFine);
    gmeshFine = GeometryFine.GeometricGmshMesh(filename,nullptr,false);
    return gmeshFine;
}
int H1Fracs(){
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
    dim_name_and_physical_tagCoarse[2]["k11"] = 1;
    dim_name_and_physical_tagCoarse[1]["inlet"] = 2;
    dim_name_and_physical_tagCoarse[1]["outlet"] = 3;
    dim_name_and_physical_tagCoarse[1]["noflux"] = 4;
    dim_name_and_physical_tagCoarse[1]["SmallFract"] = 5;


    std::string filename="/Users/victorvillegassalabarria/Downloads/SingleFracture.msh";

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);
    std::ofstream file3("SingleFractureMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file3);
    TPZCompMesh *cmesh =  new TPZCompMesh(gmesh);
  
    int matId=1;
    int dim2d = 2;
    int matIdFract=5;
    int dim1d=1;
  
    TPZDarcyFlow *matDarcy = new TPZDarcyFlow(matId, dim2d);

    TPZDarcyFlow *matDarcyFract= new TPZDarcyFlow(matIdFract,dim1d);
    matDarcy->SetConstantPermeability(0.01);
    matDarcyFract->SetConstantPermeability(1e7);

    cmesh->InsertMaterialObject(matDarcy);
    int bc_id=2;
    int bc_typeN = 1;
    int bc_typeD = 0;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);
    
    int bcinletId = 2;
    int bcOutletId = 3;
    int bcNoFlux = 4;
    TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,bcNoFlux,bc_typeN,val1,val2);
    cmesh->InsertMaterialObject(face2);
    
    val2[0]=100; // Valor a ser impuesto como presión en la entrada
    TPZBndCond * face = matDarcy->CreateBC(matDarcy,bcinletId,bc_typeD,val1,val2);
    cmesh->InsertMaterialObject(face);
    
    val2[0]=10; // Valor a ser impuesto como presión en la salida
    TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,bcOutletId,bc_typeD,val1,val2);
    cmesh->InsertMaterialObject(face1);

    cmesh->InsertMaterialObject(matDarcyFract);
 
   
    cmesh->AutoBuild();
    //Esto hace que el espacio de aproxiación sea H1
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    //Inicializa el tamaño del vector solución
    cmesh->ExpandSolution();
    
  
    //CreateAnalisys
    TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh);
    bool mustOptimizeBandwidth = false;
   
    //Carga la solución a la malla computacional
    Analisys->LoadSolution();
    
   // Selecciona el método numérico para resolver el problema algebraico
    TPZStepSolver<STATE> step;
    
//    TPZSSpStructMatrix<STATE> matrix(cmesh);
    step.SetDirect(ELDLt);
  
//    Analisys->SetStructuralMatrix(matrix);
    

    Analisys->SetSolver(step);
    
    //Ensamblaje de la matriz de rigidez y vector de carga
    Analisys->Assemble();

    //Resolución del sistema algebraico
    Analisys->Solve();

    //Definición de variables escalares y vectoriales a posprocesar
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    
    //Configuración del posprocesamiento
    int ref =0; // Permite refinar la malla con la solucion obtenida
    std::string file_reservoir("DarcyFract_H1.vtk");
    Analisys->DefineGraphMesh(dim2d,scalnames,vecnames,file_reservoir);
    //Posprocesamiento
    Analisys->PostProcess(ref, dim2d);
 
    
    return 0;
}
int H1Vugs(){
      TPZGeoMesh *gmesh = new TPZGeoMesh;
      TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
      dim_name_and_physical_tagCoarse[2]["k11"] = 1;
      dim_name_and_physical_tagCoarse[2]["Vugs"] = 6;
      dim_name_and_physical_tagCoarse[1]["inlet"] = 2;
      dim_name_and_physical_tagCoarse[1]["outlet"] = 3;
      dim_name_and_physical_tagCoarse[1]["noflux"] = 4;



      
    std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/Meshes/testskelSLICE77SP.msh";
    //std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/Meshes/FewVugsMesh.msh";
    //std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/FastMesh.msh";
    //std::string filename="/home/itopo/Stokes-Darcy_Research/Vugs/testskelSLICE77SP.msh";

      gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);
     
      std::ofstream file3("TestGeoMesh2Dskel.vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file3);
      TPZCompMesh *cmesh =  new TPZCompMesh(gmesh);
    
        int matId=1;
        int dim2d = 2;
        int matIdsmallFract=2;
        int matIBigFract=2;
        int dim1d=1;
    
        TPZDarcyFlow *matDarcy = new TPZDarcyFlow(matId, dim2d);

        TPZDarcyFlow *matDarcySmallVug= new TPZDarcyFlow(6,dim2d);
  
    
        matDarcy->SetConstantPermeability(0.01);
        matDarcySmallVug->SetConstantPermeability(1e9);

    
        cmesh->InsertMaterialObject(matDarcy);
        int bc_id=2;
        int bc_typeN = 1;
        int bc_typeD = 0;
        TPZFMatrix<STATE> val1(1,1,0.0);
        TPZVec<STATE> val2(1,0.0);
        
        int bcinletId = 2;
        int bcOutletId = 3;
        int bcNoFlux = 4;
        TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,bcNoFlux,bc_typeN,val1,val2);
        cmesh->InsertMaterialObject(face2);
        
        val2[0]=100; // Valor a ser impuesto como presión en la entrada
        TPZBndCond * face = matDarcy->CreateBC(matDarcy,bcinletId,bc_typeD,val1,val2);
        cmesh->InsertMaterialObject(face);
        
        val2[0]=10; // Valor a ser impuesto como presión en la salida
        TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,bcOutletId,bc_typeD,val1,val2);
        cmesh->InsertMaterialObject(face1);

        cmesh->InsertMaterialObject(matDarcySmallVug);
     
       
        cmesh->AutoBuild();
        //Esto hace que el espacio de aproxiación sea H1
        cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
        
        //Inicializa el tamaño del vector solución
        cmesh->ExpandSolution();
        
      
        //CreateAnalisys
        TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh);
        bool mustOptimizeBandwidth = false;
       
        //Carga la solución a la malla computacional
        Analisys->LoadSolution();
        
       // Selecciona el método numérico para resolver el problema algebraico
        TPZStepSolver<STATE> step;
        
    //    TPZSSpStructMatrix<STATE> matrix(cmesh);
        step.SetDirect(ELDLt);
      
    //    Analisys->SetStructuralMatrix(matrix);
        

        Analisys->SetSolver(step);
        
        //Ensamblaje de la matriz de rigidez y vector de carga
        Analisys->Assemble();

        //Resolución del sistema algebraico
        Analisys->Solve();

        //Definición de variables escalares y vectoriales a posprocesar
        TPZStack<std::string,10> scalnames, vecnames;
        vecnames.Push("Flux");
        scalnames.Push("Pressure");
        
        //Configuración del posprocesamiento
        int ref =0; // Permite refinar la malla con la solucion obtenida
        std::string file_reservoir("Darcy_H1.vtk");
        Analisys->DefineGraphMesh(dim2d,scalnames,vecnames,file_reservoir);
        //Posprocesamiento
        Analisys->PostProcess(ref, dim2d);
     
        return 0;
}
int Hdiv_Fract(){
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
    dim_name_and_physical_tagCoarse[2]["k11"] = 1;
    dim_name_and_physical_tagCoarse[1]["inlet"] = 2;
    dim_name_and_physical_tagCoarse[1]["outlet"] = 3;
    dim_name_and_physical_tagCoarse[1]["noflux"] = 4;
    dim_name_and_physical_tagCoarse[1]["SmallFract"] = 5;


    std::string filename="/Users/victorvillegassalabarria/Downloads/SingleFracture.msh";

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);
    std::ofstream file3("SingleFractureMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file3);
    
    
    int order =1;
    TPZHDivApproxCreator hdivCreator(gmesh);
    hdivCreator.ProbType() = ProblemType::EDarcy;
    hdivCreator.SetDefaultOrder(order);
    hdivCreator.SetShouldCondense(false);
    
    
    // Add materials (weak formulation)
    TPZMixedDarcyFlow *matDarcy = new TPZMixedDarcyFlow(1,2);
    TPZMixedDarcyFractureFlow *matDarcyFracts= new TPZMixedDarcyFractureFlow(5,1);
    
    
    matDarcy->SetConstantPermeability(0.01);
    hdivCreator.InsertMaterialObject(matDarcy);
    
    matDarcyFracts->SetConstantPermeability(1e7);


    //hdivCreator.InsertMaterialObject(matDarcyFracts);
    
    int bc_id=2;
    int bc_typeN = 1;
    int bc_typeD = 0;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);
    
    int bcinletId = 2;
    int bcOutletId = 3;
    int bcNoFlux = 4;

    TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,bcNoFlux,bc_typeN,val1,val2);

    hdivCreator.InsertMaterialObject(face2);

    
    val2[0]=100; // Valor a ser impuesto como presión en la entrada
    TPZBndCond * face = matDarcy->CreateBC(matDarcy,bcinletId,bc_typeD,val1,val2);

    hdivCreator.InsertMaterialObject(face);

    val2[0]=14; // Valor a ser impuesto como presión en la salida
    TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,bcOutletId,bc_typeD,val1,val2);

    hdivCreator.InsertMaterialObject(face1);

    
    int lagmultilevel = 1;
    TPZManVector<TPZCompMesh *, 7> meshvec(2);
    
    auto flux_mesh=meshvec[0];
    auto pressure_mesh=meshvec[1];
    
    flux_mesh->InsertMaterialObject(matDarcyFracts);
    pressure_mesh->InsertMaterialObject(matDarcyFracts);
    hdivCreator.CreateAtomicMeshes(meshvec, lagmultilevel); // This method increments the lagmultilevel
    TPZMultiphysicsCompMesh *cmesh = nullptr;

    hdivCreator.CreateMultiPhysicsMesh(meshvec, lagmultilevel, cmesh);

    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();

    for (int iel = 0; iel < cmesh->NElements(); iel++) {
        auto cel = cmesh->Element(iel);
        if (!cel) continue;

        if (cel->Material()->Id() == 3) {
            std::cout << "BC element with ndof = "
                      << cel->NConnects() << std::endl;
        }
    }

    TPZLinearAnalysis anMixed(cmesh, RenumType::EMetis);
  #ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> matMixed(cmesh);
  #else
    TPZFStructMatrix<STATE> matMixed(cmesh);
  #endif
    matMixed.SetNumThreads(0);
    anMixed.SetStructuralMatrix(matMixed);
    TPZStepSolver<STATE> stepMixed;
    stepMixed.SetDirect(ELDLt);
    anMixed.SetSolver(stepMixed);
    anMixed.Run();

    // ---- Plotting ---

    {
      const std::string plotfile = "darcyFract_mixed";
      constexpr int vtkRes{0};
      TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
      auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
      vtk.Do();
    }

    // --- Clean up ---
    delete cmesh;
    return 0;
}
int Hdiv_MixedCT(){
    
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
    dim_name_and_physical_tagCoarse[2]["k11"] = 1;
    dim_name_and_physical_tagCoarse[1]["inlet"] = 2;
    dim_name_and_physical_tagCoarse[1]["outlet"] = 3;
    dim_name_and_physical_tagCoarse[1]["noflux"] = 4;
    dim_name_and_physical_tagCoarse[2]["Vugs"] = 6;

    
    std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/Meshes/testskelSLICE77SP.msh";
    //std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/FewVugsMesh.msh";
    //std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/FastMesh.msh";
    //std::string filename="/home/itopo/Stokes-Darcy_Research/Vugs/testskelSLICE77SP.msh";

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);

    
 
    std::ofstream file20("TestGeoMesh2D.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file20);
    
    
    int order =1;
    TPZHDivApproxCreator hdivCreator(gmesh);
    hdivCreator.ProbType() = ProblemType::EDarcy;
    hdivCreator.SetDefaultOrder(order);
    hdivCreator.SetShouldCondense(false);
    
    
    // Add materials (weak formulation)
    TPZMixedDarcyFlow *matDarcy = new TPZMixedDarcyFlow(1,2);
    TPZMixedDarcyFlow *matDarcyVugs= new TPZMixedDarcyFlow(6,1);


    matDarcy->SetConstantPermeability(0.01);
    hdivCreator.InsertMaterialObject(matDarcy);
    
    matDarcyVugs->SetConstantPermeability(1e9);


    hdivCreator.InsertMaterialObject(matDarcyVugs);
    
    int bc_id=2;
    int bc_typeN = 1;
    int bc_typeD = 0;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);
    
    int bcinletId = 2;
    int bcOutletId = 3;
    int bcNoFlux = 4;

    TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,bcNoFlux,bc_typeN,val1,val2);

    hdivCreator.InsertMaterialObject(face2);

    
    val2[0]=100; // Valor a ser impuesto como presión en la entrada
    TPZBndCond * face = matDarcy->CreateBC(matDarcy,bcinletId,bc_typeD,val1,val2);

    hdivCreator.InsertMaterialObject(face);

    val2[0]=14; // Valor a ser impuesto como presión en la salida
    TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,bcOutletId,bc_typeD,val1,val2);

    hdivCreator.InsertMaterialObject(face1);

    
    int lagmultilevel = 1;
    TPZManVector<TPZCompMesh *, 7> meshvec(2);
    TPZGeoEl *gel = gmesh->Element(10);
    hdivCreator.CreateAtomicMeshes(meshvec, lagmultilevel); // This method increments the lagmultilevel
    TPZMultiphysicsCompMesh *cmesh = nullptr;

    hdivCreator.CreateMultiPhysicsMesh(meshvec, lagmultilevel, cmesh);
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();

    for (int iel = 0; iel < cmesh->NElements(); iel++) {
        auto cel = cmesh->Element(iel);
        if (!cel) continue;

        if (cel->Material()->Id() == 3) {
            std::cout << "BC element with ndof = "
                      << cel->NConnects() << std::endl;
        }
    }

    TPZLinearAnalysis anMixed(cmesh, RenumType::EMetis);
  #ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> matMixed(cmesh);
  #else
    TPZFStructMatrix<STATE> matMixed(cmesh);
  #endif
    matMixed.SetNumThreads(0);
    anMixed.SetStructuralMatrix(matMixed);
    TPZStepSolver<STATE> stepMixed;
    stepMixed.SetDirect(ELDLt);
    anMixed.SetSolver(stepMixed);
    anMixed.Run();

    // ---- Plotting ---

    {
      const std::string plotfile = "darcy_mixed";
      constexpr int vtkRes{0};
      TPZManVector<std::string, 2> fields = {"Flux", "Pressure", "GradFluxX", "GradFluxY"};
      auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
      vtk.Do();
    }

    // --- Clean up ---
    delete cmesh;
    return 0;
}
void findElDim(TPZStack<TPZGeoElSide> &allneigh, int dim, TPZStack<TPZGeoElSide> &allneighdim){
    int nels = allneigh.size();
    for (int iel =0; iel<nels; iel++) {
        if (allneigh[iel].Element()->Dimension()==dim) {
            allneighdim.push_back(allneigh[iel]);
        }
    }
}


TPZCompMesh *HdivMesh(TPZGeoMesh *gmesh){
    int matId=1;
    int dim2d = 2;
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
 
    std::set<int> volId, bcId;
    GetAtomicIds(gmesh, volId, bcId);
    insertAtomicMaterials(cmesh, volId, bcId);
    
    //int pOrder=1;
    int pOrder=1;
    cmesh->SetDefaultOrder(pOrder);
    int meshdim = gmesh->Dimension();

    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    std::cout<<cmesh->NEquations() <<std::endl;
    return cmesh;


    //devuelve una malla Hdiv
}
TPZCompMesh *Pressuremesh(TPZGeoMesh *gmesh,int order){
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    std::set<int> volId, bcId;
    GetAtomicIds(gmesh, volId, bcId);
    insertAtomicMaterials(cmesh, volId, bcId);

   
    cmesh->AutoBuild();
    
    if(order>0){
        cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    else {
        cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    
    if(1 > 0){
           int64_t ncon = cmesh->NConnects();
           for(int64_t i=0; i<ncon; i++){
               TPZConnect &newnod = cmesh->ConnectVec()[i];
               newnod.SetLagrangeMultiplier(1);
           }
       }
    
    
    cmesh->SetDefaultOrder(order);

    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    
    return cmesh;

    //devuelve una malla L2
}
int main (){
    //H1Vugs();
    Hdiv_MixedCT();
    return 0;
}
void insertAtomicMaterials(TPZCompMesh *cmesh, std::set<int> matIdsVol, std::set<int> matIdsBcs){
    
    int dim = cmesh->Dimension();
    
    for (auto iD:matIdsVol) {
        TPZNullMaterial <STATE> *matDarcy = new TPZNullMaterial(iD, dim);
        cmesh->InsertMaterialObject(matDarcy);
      
    }
    for (auto iD:matIdsBcs) {
        
        TPZNullMaterial<STATE> * face2 = new TPZNullMaterial(iD, dim-1);
        cmesh->InsertMaterialObject(face2);

    }
}
void GetAtomicIds(TPZGeoMesh *geomesh, std::set<int> &volId, std::set<int> &bcId){
    int dim = geomesh->Dimension();
    for (auto gel: geomesh->ElementVec()) {
        if (! gel) {
            continue;
        }
        int matId = gel->MaterialId();
        int geldim = gel->Dimension();
        if (geldim == dim) {
            volId.insert(matId);
        }
        if (geldim == dim-1) {
            bcId.insert(matId);
        }
    }
    
    
}
