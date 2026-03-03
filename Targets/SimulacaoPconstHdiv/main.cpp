#include "sources.h"


void Hdiv_MixedCT(){
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZGeoMesh *gmeshp = new TPZGeoMesh;

    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
    dim_name_and_physical_tagCoarse[2]["k11"] = EMatId;
    dim_name_and_physical_tagCoarse[1]["inlet"] = EbcInletId;
    dim_name_and_physical_tagCoarse[1]["outlet"] = EbcOutletId;
    dim_name_and_physical_tagCoarse[1]["noflux"] = EbcNoFlux;
    dim_name_and_physical_tagCoarse[2]["Vugs"] = 6;
    
    //std::string filename="/Users/victorvillegassalabarria/python-test/testskelSLICE77SP.msh";
    //std::string filename="/Users/victorvillegassalabarria/Documents/Github/Vugs/FastMesh.msh";
    //std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/FastMesh.msh"; 
    std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/FewVugsMesh.msh";

    std::ofstream file20("TestGeoMesh2D.vtk");
    // std::ofstream file21("Test_cmeshFlux.vtk");
    // std::ofstream file22("Test_cmeshPressure.vtk");
    // std::ofstream file23("Test_cmeshPressure.txt");
    // std::ofstream file24("Test_cmeshFlux.txt"); 
    // std::ofstream file25("Test_cmeshMulti.txt");

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);
    
    MeshWithSegmentVugs(gmesh);

    std::set<int> volId, bcId;
    GetAtomicIds(gmesh, volId, bcId);

    TPZCompMesh *Flux_cmesh=CreateFluxMesh(gmesh,volId,bcId);
    
    TPZCompMesh *Pressure_cmesh=CreatePressureMesh(gmesh,volId,bcId,0);
    
    TPZMultiphysicsCompMesh *cmesh_mult= new TPZMultiphysicsCompMesh(gmesh);
    cmesh_mult->SetName("MultiMesh");
    TPZVec<TPZCompMesh *> meshvec(2);
    
    meshvec[0]= Flux_cmesh;
    meshvec[1]= Pressure_cmesh;
    
    // Add materials (weak formulation)
    TPZMixedDarcyFlow *matDarcy = new TPZMixedDarcyFlow(EMatId,2);
    //TPZMixedDarcyFlow *matDarcyVugs= new TPZMixedDarcyFlow(EVugId,2); 

    //TODO VERIFICAR SE É ISSO fazer para cada vug
    for(auto vugId: vugIds) {
        TPZMixedDarcyFlow *matDarcyVugs= new TPZMixedDarcyFlow(vugId,2); 
        cmesh_mult->InsertMaterialObject(matDarcyVugs);
    }

    cmesh_mult->InsertMaterialObject(matDarcy);
    //cmesh_mult->InsertMaterialObject(matDarcyVugs);
    matDarcy->SetConstantPermeability(0.01);

    int bc_id=2;
    int bc_typeN = 1;
    int bc_typeD = 0;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);
    int dim2d=2;

    TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,EbcNoFlux,bc_typeN,val1,val2);
    cmesh_mult->InsertMaterialObject(face2);
    
    val2[0]=100; // Valor a ser impuesto como presión en la entrada
    TPZBndCond * face = matDarcy->CreateBC(matDarcy,EbcInletId,bc_typeD,val1,val2);
    cmesh_mult->InsertMaterialObject(face);

    val2[0]=14; // Valor a ser impuesto como presión en la salida
    TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,EbcOutletId,bc_typeD,val1,val2);
    cmesh_mult->InsertMaterialObject(face1);

    //TODO Create comp elements of Vug Boundary
    for(auto bcId: vugBcIds) {
        TPZBndCond *faceVug = matDarcy->CreateBC(matDarcy,bcId,bc_typeD,val1,val2);
        cmesh_mult->InsertMaterialObject(faceVug);
    }

    cmesh_mult->ExpandSolution();
    cmesh_mult->ApproxSpace().Style()= TPZCreateApproximationSpace::EMultiphysics;
    cmesh_mult->BuildMultiphysicsSpace(meshvec);

    CreateInterfaceGeoEls(gmesh);
    InsertInterfaceEls(cmesh_mult, gmesh);

    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file20);

    cmesh_mult->InitializeBlock();
    std::cout<<cmesh_mult->Element(1)<<std::endl;
    bool mustOp = false;

    //Show Shape
    TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(Pressure_cmesh);

    PrintCompMesh(Flux_cmesh);
    PrintCompMesh(Pressure_cmesh);
    PrintCompMesh(cmesh_mult);

    const std::string strShape = "Shape.vtk";
    TPZVec<int64_t> eqIndices(4, 0);
    eqIndices[0] = 293;
    eqIndices[1] = 294;
    eqIndices[2] = 295;
    eqIndices[3] = 296;


    Analisys->ShowShape(strShape, eqIndices);
    

    //TODO ANALYSIS
    //TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh_mult);
        

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
}


int main (){
    Hdiv_MixedCT();
    return 0;
}