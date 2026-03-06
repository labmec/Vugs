#include "sources.h"
#include "pzlog.h"
#include "Projection/TPZL2ProjectionCS.h"

void Hdiv_MixedCT(){
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZGeoMesh *gmeshp = new TPZGeoMesh;

    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
    dim_name_and_physical_tagCoarse[2]["k11"] = EMatId;
    dim_name_and_physical_tagCoarse[2]["Vugs"] = EVugId;
    dim_name_and_physical_tagCoarse[1]["inlet"] = EbcInletId;
    dim_name_and_physical_tagCoarse[1]["outlet"] = EbcOutletId;
    dim_name_and_physical_tagCoarse[1]["noflux"] = EbcNoFlux;

    
    //std::string filename="/Users/victorvillegassalabarria/python-test/testskelSLICE77SP.msh";
    //std::string filename="/Users/victorvillegassalabarria/Documents/Github/Vugs/FastMesh.msh";
    std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/FastMesh.msh";
    //std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/FewVugsMesh.msh";
    std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/FastMesh.msh";
    //std::string filename="/home/itopo/Stokes-Darcy_Research/Vugs/testskelSLICE77SP.msh";

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);

    
    std::ofstream file20("TestGeoMesh2D.vtk");
    // std::ofstream file21("Test_cmeshFlux.vtk");
    // std::ofstream file22("Test_cmeshPressure.vtk");
    // std::ofstream file23("Test_cmeshPressure.txt");
    // std::ofstream file24("Test_cmeshFlux.txt"); 
    // std::ofstream file25("Test_cmeshMulti.txt");

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);
    
    MeshWithSegmentVugs(gmesh); //!FAZER COM QUE ESSA FUNÇÃO RETORNE OS SETS THE IDS

    std::set<int> volId, bcId;
    GetAtomicIds(gmesh, volId, bcId);
    TPZCompMesh *Flux_cmesh = CreateFluxMesh(gmesh,volId,bcId);
    TPZCompMesh *Pressure_cmesh = CreatePressureMesh(gmesh,volId,bcId,0);
    
    TPZMultiphysicsCompMesh *cmesh_mult= new TPZMultiphysicsCompMesh(gmesh);
    TPZVec<TPZCompMesh *> meshvec(2);
    
    meshvec[0]= Flux_cmesh;
    meshvec[1]= Pressure_cmesh;
    
    
    // Add materials (weak formulation)
    TPZMixedDarcyFlow *matDarcy = new TPZMixedDarcyFlow(EMatId,2);
    // TPZMixedDarcyFlow *matDarcyVugs= new TPZMixedDarcyFlow(EVugId,2); 

   for(auto vugId: vugIds) { //TODO Vug elements
        //auto *matDarcyVugs= new TPZNullMaterialCS<REAL> (vugId,2,1);
        //auto *matDarcyVugs= new TPZDarcyFlow (vugId,2);
        auto *matDarcyVugs= new TPZL2ProjectionCS<REAL> (vugId,2,1);
        matDarcyVugs->SetScaleFactor(0);
        cmesh_mult->InsertMaterialObject(matDarcyVugs);
   }

    cmesh_mult->InsertMaterialObject(matDarcy);
    // cmesh_mult->InsertMaterialObject(matDarcyVugs);
    matDarcy->SetConstantPermeability(0.01);
    
    int bc_typeN = 1;
    int bc_typeD = 0;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);
    int dim2d=2;
    int dim2d=2;

    val2[0]=1;
    TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,EbcNoFlux,bc_typeD,val1,val2);
    cmesh_mult->InsertMaterialObject(face2);
                
    val2[0]=1; // Valor a ser impuesto como presión en la entrada
    TPZBndCond * face = matDarcy->CreateBC(matDarcy,EbcInletId,bc_typeD,val1,val2);
    cmesh_mult->InsertMaterialObject(face);

    val2[0]=1; // Valor a ser impuesto como presión en la salida
    TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,EbcOutletId,bc_typeD,val1,val2);
    cmesh_mult->InsertMaterialObject(face1);


    //TODO Create comp elements of Vug Boundary
    val2[0] = 0;
    for(auto bcId: vugBcIds) {
        TPZBndCond *faceVug = matDarcy->CreateBC(matDarcy,bcId,bc_typeD,val1,val2);
        cmesh_mult->InsertMaterialObject(faceVug);
    }

    cmesh_mult->ExpandSolution();
    cmesh_mult->ApproxSpace().Style()= TPZCreateApproximationSpace::EMultiphysics;
    cmesh_mult->BuildMultiphysicsSpace(meshvec);
    // cmesh_mult->CleanUpUnconnectedNodes();

    CreateInterfaceGeoEls(gmesh); //TODO Vug elements
    InsertInterfaceEls(cmesh_mult, gmesh); //TODO Vug elements
    SideOrientation(Flux_cmesh);

    PrintCompMesh(Pressure_cmesh);
    PrintCompMesh(Flux_cmesh);

    cmesh_mult->InitializeBlock();
    std::cout<<cmesh_mult->Element(1)<<std::endl;
    
    CreateInterfaceGeoEls(gmesh);
    InsertInterfaceEls(cmesh_mult, gmesh); //TODO
        
    bool mustOp = false;

    //Show Shape

    PrintCompMesh(Flux_cmesh);
    PrintCompMesh(Pressure_cmesh);
    PrintCompMesh(cmesh_mult);

    const std::string strShape = "Shape.vtk";
    TPZVec<int64_t> eqIndices(4, 0);
    eqIndices[0] = 293;
    eqIndices[1] = 294;
    eqIndices[2] = 295;
    eqIndices[3] = 296;


//    Analisys->ShowShape(strShape, eqIndices);


    cmesh_mult->Reference()->ResetReference();
    cmesh_mult->LoadReferences();
    //TPZLinearAnalysis anMixed(cmesh_mult,RenumType::EMetis);
    TPZLinearAnalysis anMixed(cmesh_mult,RenumType::EMetis);

    //anMixed->ShowShape(strShape, eqIndices);//new TPZLinearAnalysis(cmesh_mult);
    #ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> matMixed(cmesh_mult);
    #else
    TPZFStructMatrix<STATE> matMixed(cmesh_mult);
    #endif
    matMixed.SetNumThreads(0);
    anMixed.SetStructuralMatrix(matMixed);
    TPZStepSolver<STATE> stepMixed;
    stepMixed.SetDirect(ELDLt);
    anMixed.SetSolver(stepMixed);
    anMixed.Run();

        {
          const std::string plotfile = "Darcy_mixed";
          constexpr int vtkRes{0};
          TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
          auto vtk = TPZVTKGenerator(cmesh_mult, fields, plotfile, vtkRes);
          vtk.Do();
        }
    PrintCompMesh(cmesh_mult);
        // --- Clean up ---
        delete cmesh_mult;
}


int main (){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    Hdiv_MixedCT();
    return 0;
}