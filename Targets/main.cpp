#include <iostream>
#include <filesystem>
#include <math.h>
#include "pzcmesh.h"
#include "TPZElementMatrixT.h"

//#include <opencv2/opencv.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
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
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZNullMaterialCS.h"
#include "TPZNullMaterial.h"
#include "TPZAnalysis.h"
//#include "TPZCreateMultiphysicsSpace.h"
#include "pzstepsolver.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"
//#include "TPZStepSolver.h"
//using std::cout;
//using std::endl;
//using std::cin;
#include <set>
#include "TPZAnalyticSolution.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TSFMixedDarcy.h"
#include "TPZHDivApproxCreator.h"
//hola
int mainDarcy2d();
int mainDarcy3D();
TPZCompMesh* HdivMesh(TPZGeoMesh *);
TPZCompMesh* Pressuremesh(TPZGeoMesh *, int order);
void GetAtomicIds(TPZGeoMesh *geomesh, std::set<int> &volId, std::set<int> &bcId);
void insertAtomicMaterials(TPZCompMesh *cmesh, std::set<int> matIdsVol, std::set<int> matIdsBcs);

//using namespace cv;
//using namespace std;
//Function to generate a mesh using gmsh library
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void findElDim(TPZStack<TPZGeoElSide> &allneigh, int dim, TPZStack<TPZGeoElSide> &allneighdim);


int main3D();
int main2D();
int main2DFracVug();
int mainDarcy3D ();

//int main(){
//
//    return mainDarcy3D();
//}


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

int main2DFracVug(){
      TPZGeoMesh *gmesh = new TPZGeoMesh;
      TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
      dim_name_and_physical_tagCoarse[2]["k11"] = 1;
      //dim_name_and_physical_tagCoarse[2]["SmallVug"] = 7;
      //dim_name_and_physical_tagCoarse[2]["BigVug"] = 8;
      dim_name_and_physical_tagCoarse[2]["Vugs"] = 6;
      dim_name_and_physical_tagCoarse[1]["inlet"] = 2;
      dim_name_and_physical_tagCoarse[1]["outlet"] = 3;
      dim_name_and_physical_tagCoarse[1]["noflux"] = 4;
      //dim_name_and_physical_tagCoarse[1]["SmallFract"] = 5;
      //dim_name_and_physical_tagCoarse[1]["BigFract"] = 6;


      
      //std::string filename="/Users/victorvillegassalabarria/python-test/testskel4.msh";
      //std::string filename="/Users/victorvillegassalabarria/python-test/testskel30sp.msh";
      std::string filename="/Users/victorvillegassalabarria/python-test/testskelSLICE77SP.msh";

      gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);
     
      std::ofstream file3("TestGeoMesh2Dskel.vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file3);
      //Create CompMesh
      TPZCompMesh *cmesh =  new TPZCompMesh(gmesh);
    
        //Create Materials
        int matId=1;
        int dim2d = 2;
        int matIdsmallFract=2;
        int matIBigFract=2;
        int dim1d=1;
    
        TPZDarcyFlow *matDarcy = new TPZDarcyFlow(matId, dim2d);
//        TPZDarcyFlow *matDarcySmallFract= new TPZDarcyFlow(5,dim1d);
        //TPZDarcyFlow *matDarcyBigFract= new TPZDarcyFlow(6,dim1d);
        //TPZDarcyFlow *matDarcySmallVug= new TPZDarcyFlow(7,dim2d);
        TPZDarcyFlow *matDarcySmallVug= new TPZDarcyFlow(6,dim2d);
        //TPZDarcyFlow *matDarcyBigVug= new TPZDarcyFlow(8,dim2d);
    
        matDarcy->SetConstantPermeability(0.01);
        matDarcySmallVug->SetConstantPermeability(1e9);
        //matDarcyBigVug->SetConstantPermeability(1.0e9);
//        matDarcySmallFract->SetConstantPermeability(1e9);
        //matDarcyBigFract->SetConstantPermeability(1.0e9);
    int x, y;

//    // Definir una función de permeabilidad como una lambda
//    PermeabilityFunctionType perm_function = [](const TPZVec<REAL>& coord) -> STATE {
//        if (coord[0]>100 and coord[1]>100){
//            return 1000;
//        }
//        else{
//            return 1;
//
//        };
//    };
    PermeabilityFunctionType perm_function = [](const TPZVec<REAL>& coord) -> STATE {
        REAL x = coord[0];
        REAL y = coord[1];
        REAL arg = 2 * M_PI * x + 2 * M_PI * y;
        REAL cos_arg = cos(arg);
        REAL exp_term = exp(2.3 * cos_arg);
        
        return exp_term;
    };
//        // Ejemplo: Permeabilidad depende de x (coord[0])
//        return coord[0] * 1e-3;
//
//
    //matDarcy->SetPermeabilityFunction(perm_function);
    //Conseguir permeabilidade em um ponto da malha coord(x,y);
    
    TPZVec<REAL> coord(2);
    coord[0]=90.5;
    coord[1]=650;
    auto Perm=matDarcy->GetPermeability(coord);
    std::cout<<Perm<<std::endl;
    //Conseguir permeabilidade em um ponto da malha coord(x,y);
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
        //cmesh->InsertMaterialObject(matDarcySmallFract);
        //cmesh->InsertMaterialObject(matDarcyBigFract);
        cmesh->InsertMaterialObject(matDarcySmallVug);
        //cmesh->InsertMaterialObject(matDarcyBigVug);
       
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
//TPZ MixedDarcy Flow
int mainMixed(){
    
    //
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
    dim_name_and_physical_tagCoarse[2]["k11"] = 1;
    dim_name_and_physical_tagCoarse[1]["inlet"] = 2;
    dim_name_and_physical_tagCoarse[1]["outlet"] = 3;
    dim_name_and_physical_tagCoarse[1]["noflux"] = 4;
    

    
    //std::string filename="/Users/victorvillegassalabarria/python-test/testskel4.msh";
    //std::string filename="/Users/victorvillegassalabarria/python-test/testskel30sp.msh";
    std::string filename="/Users/victorvillegassalabarria/Downloads/MallaTriangles.msh";

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);
   
    std::ofstream file3("TestGeoMesh2Dskel.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file3);
    //Create CompMesh
    TPZCompMesh *cmesh_flux =  HdivMesh(gmesh);
    //TPZCompMesh *cmesh =  Pressuremesh(gmesh,1);
    TPZCompMesh *cmesh =  Pressuremesh(gmesh,0);
  //  TPZCompMesh *cmesh_paverage =  Pressuremesh(gmesh,0);
   // TPZCompMesh *cmesh_faverage =  Pressuremesh(gmesh,0);

    
    //MULTIFISICA
    TPZMultiphysicsCompMesh *cmesh_mult= new TPZMultiphysicsCompMesh(gmesh);
    TPZVec<TPZCompMesh *> meshvec(2);
    meshvec[0]= cmesh_flux;
    meshvec[1]= cmesh;
    //meshvec[2]= cmesh_paverage;
    //meshvec[3]= cmesh_faverage;
   //Config
    //cmesh_mult->BuildMultiphysicsSpace(meshvec);
    
  
      //Create Materials
      int matId=1;
      int dim2d = 2;
      int matIdsmallFract=2;
      int matIBigFract=2;
      int dim1d=1;
  
      TPZMixedDarcyFlow *matDarcy = new TPZMixedDarcyFlow(matId, dim2d);
 
  
  //    matDarcy->SetConstantPermeability(1);

  
      cmesh_mult->InsertMaterialObject(matDarcy);
      int bc_id=2;
      int bc_typeN = 1;
      int bc_typeD = 0;
      TPZFMatrix<STATE> val1(1,1,0.0);
      TPZVec<STATE> val2(1,0.0);
      
      int bcinletId = 2;
      int bcOutletId = 3;
      int bcNoFlux = 4;
      TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,bcNoFlux,bc_typeN,val1,val2);
      cmesh_mult->InsertMaterialObject(face2);
      
      val2[0]=100; // Valor a ser impuesto como presión en la entrada
      TPZBndCond * face = matDarcy->CreateBC(matDarcy,bcinletId,bc_typeD,val1,val2);
      cmesh_mult->InsertMaterialObject(face);
      
      val2[0]=10; // Valor a ser impuesto como presión en la salida
      TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,bcOutletId,bc_typeD,val1,val2);
      cmesh_mult->InsertMaterialObject(face1);
 
     
    //  cmesh_mult->AutoBuild();
      //Esto hace que el espacio de aproxiación sea H1
    //cmesh_mult->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh_mult->ApproxSpace().SetAllCreateFunctionsContinuous();

      //Inicializa el tamaño del vector solución
      cmesh_mult->ExpandSolution();
    cmesh_mult->ApproxSpace().Style()= TPZCreateApproximationSpace::EMultiphysics;
    cmesh_mult->BuildMultiphysicsSpace(meshvec);
    //TPZManVector<int, 2> active_approx_spaces(2, 1);
    //cmesh_mult->BuildMultiphysicsSpace(active_approx_spaces, meshvec);
    cmesh_mult->InitializeBlock();
    std::cout<<cmesh_mult->Element(1)<<std::endl;
    
    bool mustOp = false;

      //CreateAnalisys
    TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh_mult);
    
    //new TPZLinearAnalysis(cmesh_mult);
    

      //TPZAnalysis *Analisys = new TPZAnalysis(cmesh_mult,true);
//      bool mustOptimizeBandwidth = false;
     
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
      //Analisys->Solve();
    TPZElementMatrixT<double> mat, vec;
    std::ofstream file("matrixel.txt");
    int nels =cmesh_mult->NElements();
    for (int i=1;i<nels;i++){
        auto cel=cmesh_mult->Element(i);
        auto gel=cel->Reference();
        if(gel->Dimension()==2){
            cel->CalcStiff(mat,vec);
            mat.fMat.Print(file);
        }
    }
      //Definición de variables escalares y vectoriales a posprocesar
      TPZStack<std::string,10> scalnames, vecnames;
      vecnames.Push("Flux");
      scalnames.Push("Pressure");
      
      //Configuración del posprocesamiento
      int ref =0; // Permite refinar la malla con la solucion obtenida
      std::string file_reservoir("SolVictorCTmesh.vtk");
      Analisys->DefineGraphMesh(dim2d,scalnames,vecnames,file_reservoir);
      //Posprocesamiento
      Analisys->PostProcess(ref, dim2d);
   
      return 0;
    //
   
}
int mainMixedCT(){
    
    //
    TPZVec<int> nx(2, 10);
    //What does it mean the line 46?
    
    nx[0]=10;
    nx[1]=10;
    const TPZVec<REAL> x0(3, 0.);
    const TPZVec<REAL> x1(3, 0.);
    x1[0]=1.0;
    x1[1]=1.0;
    auto msh = TPZGenGrid2D(nx, x0, x1);
    
//    TPZGeoMesh *gmesh = new TPZGeoMesh;
//
//    msh.Read(gmesh);
//    msh.SetElementType(MMeshType::EQuadrilateral);
//
//    msh.SetBC(gmesh, 7, 2);
//    msh.SetBC(gmesh, 5, 3);
//    msh.SetBC(gmesh, 4, 4);
//    msh.SetBC(gmesh, 6, 4);
//    std::ofstream file20("TestGeoMesh2D.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file20);
    //
    //TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
    dim_name_and_physical_tagCoarse[2]["k11"] = 1;
    dim_name_and_physical_tagCoarse[1]["inlet"] = 2;
    dim_name_and_physical_tagCoarse[1]["outlet"] = 3;
    dim_name_and_physical_tagCoarse[1]["noflux"] = 4;
    dim_name_and_physical_tagCoarse[2]["Vugs"] = 6;

    
    //std::string filename="/Users/victorvillegassalabarria/python-test/testskel4.msh";
    //std::string filename="/Users/victorvillegassalabarria/python-test/testskel30sp.msh";
    //std::string filename="/Users/victorvillegassalabarria/Downloads/MallaTriangles.msh";
    std::string filename="/Users/victorvillegassalabarria/python-test/testskelSLICE77SP.msh";

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);

    
  
    // inserta contornos mat 5 / mat 6
//    int ncreated = 0;
//    int nels = gmesh->NElements();
//
//
//    for (int iel = 0; iel< nels; iel++) {
//        TPZGeoEl *gel = gmesh->Element(iel);
//        if (!gel){
//            continue;
//        }
//        if (gel->Dimension() != 2) {
//            continue;
//        }
//        int nsides= gel->NSides();
//        int ncorners= gel->NCornerNodes();
//        int firstside= nsides-ncorners-1;
//
//
//        for (int iside = firstside; iside<nsides; iside++) {
//            TPZGeoElSide gelside(gel, iside);
//            int matid = gelside.Element()->MaterialId();
//            TPZStack<TPZGeoElSide> allneigh;
//            gelside.AllNeighbours(allneigh);
////            std::cout<<allneigh[0].Element()<<std::endl;
//            int nneighs = allneigh.size();
//            //verify Dimension
//            int verify =0;
//
//            for (int ineigh=0; ineigh<nneighs; ineigh++) {
//                TPZGeoEl *gelneigh = allneigh[ineigh].Element();
//                int dimen = gelneigh->Dimension();
//
//                if (dimen == 1) {
//                    verify = 1;
//                }
//            }
//            if(verify == 1){
//                continue;
//            }
//
//            for (int ineigh=0; ineigh<nneighs; ineigh++) {
//                TPZGeoEl *gelneigh = allneigh[ineigh].Element();
//                int matNeigh = gelneigh->MaterialId();
//                if (matNeigh != matid && (gel->Dimension() == gelneigh->Dimension()) ) {
//                   gelside.Element()->CreateBCGeoEl(iside, 100);
//                   ncreated++;
//                }
//            }
//        }
//    }
//
//    std::cout<< "se crearon: " << ncreated << " elements"<<std::endl;
//
//    gmesh->BuildConnectivity();
//    int nels2 = gmesh->NElements();
//    TPZVec<int> verificador(nels2, 0);
//    // creador de contornos por ids
//    int mat=100;
//    for (int iel =nels-1; iel<nels2; iel++) {
//
//        TPZGeoEl * gel = gmesh->Element(iel);
//        if (gel->MaterialId() ==100) {
//            std::cout<<"ok "<<std::endl;
//        }
//        if (!gel) {
//            continue;
//        }
//        if (gel->Dimension() != 1) {
//            continue;
//        }
//        if (verificador[iel]==1) {
//            continue;
//        }
//        if (gel->MaterialId() != 100) {
//            continue;
//        }
//        int side = 1;
//
//        TPZGeoElSide gelside(gel, side);
//
//        TPZStack<TPZGeoElSide> allneigh;
//        gelside.AllNeighbours(allneigh);
//        TPZStack<TPZGeoElSide> allneighdim;
//        findElDim(allneigh, 1, allneighdim);
//        int ntest = allneighdim.size();
//        TPZGeoElSide gelneigh = allneighdim[0];
//        gel->SetMaterialId(mat);
//        while (gel != gelneigh.Element()) {
//            TPZStack<TPZGeoElSide> allneigh;
//            int sidetest = gelneigh.Side();
//            if (sidetest==0) {
//                gelneigh.SetSide(1);
//            }
//            else{
//                gelneigh.SetSide(0);
//            }
//            gelneigh.AllNeighbours(allneigh);
//            TPZStack<TPZGeoElSide> allneighdim;
//            findElDim(allneigh, 1, allneighdim);
//            int indexneig = gelneigh.Element()->Index();
//            verificador[indexneig] =1;
//            gelneigh.Element()->SetMaterialId(mat);
//            gelneigh =allneighdim[0];
//        }
//
//
//        mat++;
//        int ok=0;
//    }
    
    //
    std::ofstream file20("TestGeoMesh2D.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file20);
    
    
    int order =1;
    TPZHDivApproxCreator hdivCreator(gmesh);
    hdivCreator.ProbType() = ProblemType::EDarcy;
    hdivCreator.SetDefaultOrder(order);
    hdivCreator.SetShouldCondense(false);
    
    
    // Add materials (weak formulation)
    TPZMixedDarcyFlow *matDarcy = new TPZMixedDarcyFlow(1,2);
    TPZMixedDarcyFlow *matDarcyVugs= new TPZMixedDarcyFlow(6,2);

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

    //val2[0]=0;
    TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,bcNoFlux,bc_typeN,val1,val2);
    //auto face2v = matDarcyVugs->CreateBC(matDarcyVugs,bcNoFlux,bc_typeN,val1,val2);

    hdivCreator.InsertMaterialObject(face2);
    //hdivCreator.InsertMaterialObject(face2v);

    
    val2[0]=100; // Valor a ser impuesto como presión en la entrada
    TPZBndCond * face = matDarcy->CreateBC(matDarcy,bcinletId,bc_typeD,val1,val2);
    //TPZBndCond * facev = matDarcyVugs->CreateBC(matDarcyVugs,bcinletId,bc_typeD,val1,val2);

    hdivCreator.InsertMaterialObject(face);
    //hdivCreator.InsertMaterialObject(facev);

    val2[0]=14; // Valor a ser impuesto como presión en la salida
    TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,bcOutletId,bc_typeD,val1,val2);
    //TPZBndCond * face1v = matDarcyVugs->CreateBC(matDarcyVugs,bcOutletId,bc_typeD,val1,val2);

    hdivCreator.InsertMaterialObject(face1);
    //hdivCreator.InsertMaterialObject(face1v);

    //
    //
//    for (int p=100; p<=190; p++) {
//        val2[0]=50;
//        TPZBndCond *contorno=matDarcy->CreateBC(matDarcy,p,bc_typeD,val1,val2);
//        hdivCreator.InsertMaterialObject(contorno);
//    }
    //
    //
    
    int lagmultilevel = 1;
    TPZManVector<TPZCompMesh *, 7> meshvec(hdivCreator.NumMeshes());
    hdivCreator.CreateAtomicMeshes(meshvec, lagmultilevel); // This method increments the lagmultilevel
    TPZMultiphysicsCompMesh *cmesh = nullptr;

    hdivCreator.CreateMultiPhysicsMesh(meshvec, lagmultilevel, cmesh);
//    cmesh->CleanUpUnconnectedNodes();
//    cmesh->ExpandSolution();
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
      TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
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
    main2DFracVug();
    //mainDarcy3D();
    //mainMixed();
    mainMixedCT();
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
