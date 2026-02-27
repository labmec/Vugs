#include "sources.h"
#include "DarcyFlow/TPZMixedDarcyFractureFlow.h"

//#include "sources.cpp"

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
    std::string filename="/Users/victorvillegassalabarria/Documents/Github/Vugs/FewVugsMesh.msh";
    
    std::ofstream file20("TestGeoMesh2D.vtk");
    // std::ofstream file21("Test_cmeshFlux.vtk");
    // std::ofstream file22("Test_cmeshPressure.vtk");
    // std::ofstream file23("Test_cmeshPressure.txt");
    // std::ofstream file24("Test_cmeshFlux.txt"); 
    // std::ofstream file25("Test_cmeshMulti.txt");

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);
    std::set<int> vugBcIds;
    std::set<int> vugIds;

    MeshWithSegmentVugs(gmesh, vugBcIds,vugIds);
    std::cout<<"Vug bc index final2: "<<vugBcIds.size()<<std::endl;

    //MeshWithSegmentVugs(gmesh, );

    std::set<int> volId, bcId;
    GetAtomicIds(gmesh, volId, bcId);
    int orderp=1;
    TPZCompMesh *cmesh =  new TPZCompMesh(gmesh);

    TPZCompMesh *Flux_cmesh=CreateFluxMesh(gmesh,volId,bcId,orderp);
    
    TPZCompMesh *Pressure_cmesh=CreatePressureMesh(gmesh,volId,bcId,orderp);
    
    TPZMultiphysicsCompMesh *cmesh_mult= new TPZMultiphysicsCompMesh(gmesh);
    cmesh_mult->SetName("MultiMesh");
    TPZVec<TPZCompMesh *> meshvec(2);
    
    meshvec[0]= Flux_cmesh;
    meshvec[1]= Pressure_cmesh;
    SideOrientation(Flux_cmesh);

    // Add materials (weak formulation)
    TPZMixedDarcyFlow *matDarcy = new TPZMixedDarcyFlow(EMatId,2);
    //TPZMixedDarcyFlow *matDarcyVugs= new TPZMixedDarcyFlow(EVugId,2);
    //TPZMixedDarcyFlow *matDarcyVugs= new TPZMixedDarcyFlow(600,2);
    //TODO VERIFICAR SE É ISSO fazer para cada vug
//    for(auto vugId: vugIds) {
//        std::cout<<"Vug index: "<<vugId<<std::endl;
//
//        TPZMixedDarcyFlow *matDarcyVugs= new TPZMixedDarcyFlow(vugId,2);
//        cmesh_mult->InsertMaterialObject(matDarcyVugs);
//    }
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

    val2[0]=10; // Valor a ser impuesto como presión en la salida
    TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,EbcOutletId,bc_typeD,val1,val2);
    cmesh_mult->InsertMaterialObject(face1);
    //TODO Create comp elements of Vug Boundary
    //
    //
    //
    //
//    for(auto bcId: vugBcIds) {
//        std::cout<<"Bc index: "<<bcId<<bcId>10<<std::endl;
//        int PContornoVug=30;
//        val2[0]=PContornoVug;
//        TPZBndCond *faceVug = matDarcy->CreateBC(matDarcy,bcId,bc_typeD,val1,val2);
//        cmesh_mult->InsertMaterialObject(faceVug);
//
//    }
    int Nvugs=vugBcIds.size();
    for(int i=0;i<Nvugs-46;i++){
        int bcId=0;
        int PContornoVug=0;
        if(i%2==0){
             bcId=100+(2*i);
             PContornoVug=90;
        }
        else{
             bcId=100+(2*i)+1;
             PContornoVug=30;

        }
    //int PContornoVug=50;
//        std::cout<<"Bc index: "<<bcId<<std::endl;
        val2[0]=PContornoVug;
        TPZBndCond *faceVug = matDarcy->CreateBC(matDarcy,bcId,bc_typeD,val1,val2);
        cmesh_mult->InsertMaterialObject(faceVug);

    }
    //CreateInterfaceGeoEls(gmesh,vugIds);
    //InsertInterfaceEls(cmesh_mult, gmesh,vugIds, vugBcIds);
    //
    //
    //
    //
//    int PContornoVug=50;
//    val2[0]=PContornoVug;
//    TPZBndCond *faceVug = matDarcy->CreateBC(matDarcy,300,bc_typeD,val1,val2);
//    cmesh_mult->InsertMaterialObject(faceVug);

    cmesh_mult->ExpandSolution();
    cmesh_mult->ApproxSpace().Style()= TPZCreateApproximationSpace::EMultiphysics;
    //SetUniqueVugConnect(gmesh, cmesh_mult);

    cmesh_mult->BuildMultiphysicsSpace(meshvec);
    cmesh_mult->CleanUpUnconnectedNodes();

    //SideOrientation(cmesh_mult);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file20);

    cmesh_mult->InitializeBlock();
    //std::cout<<cmesh_mult->Element(1)<<std::endl;
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
          const std::string plotfile = "darcy_mixed";
          constexpr int vtkRes{0};
          TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
          auto vtk = TPZVTKGenerator(cmesh_mult, fields, plotfile, vtkRes);
          vtk.Do();
        }

        // --- Clean up ---
        delete cmesh_mult;

    
}
void Hdiv_Fract(){
    TPZGeoMesh *gmesh = new TPZGeoMesh;
//    TPZGeoMesh *gmeshp = new TPZGeoMesh;

    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
    dim_name_and_physical_tagCoarse[2]["k11"] = EMatId;
    dim_name_and_physical_tagCoarse[1]["inlet"] = 2;
    dim_name_and_physical_tagCoarse[1]["outlet"] = 3;
    dim_name_and_physical_tagCoarse[1]["noflux"] = 4;
    dim_name_and_physical_tagCoarse[2]["Vugs"] = EVugId;

    
    //std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/testskelSLICE77SP.msh";
    std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/FastMesh.msh";
    //std::string filename="/home/itopo/Stokes-Darcy_Research/Vugs/testskelSLICE77SP.msh";

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);

    
    std::ofstream file20("TestGeoMesh2D.vtk");
    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);

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
    int dim2d=2;

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
      TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
      auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
      vtk.Do();
    }

    // --- Clean up ---
    delete cmesh;
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
}


TPZCompMesh* CreateCompMeshV(TPZGeoMesh *gmesh){
    
    TPZCompMesh *cmesh_v = new TPZCompMesh(gmesh);

    //cmesh_v->ApproxSpace().SetHDivFamily(HDivFamily::EHDivStandard);
    cmesh_v->SetAllCreateFunctionsHDiv();

    TPZNullMaterial<STATE> *nullMat = new TPZNullMaterial(matID_flux, gmesh->Dimension());  
    cmesh_v->InsertMaterialObject(nullMat);



    return cmesh_v;
}

TPZCompMesh *CreateCompMeshP(TPZGeoMesh *gmesh,int order){
    TPZCompMesh *cmesh_p = new TPZCompMesh(gmesh);
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
    
    // if(1 > 0){ //? What?
    //        int64_t ncon = cmesh->NConnects();
    //        for(int64_t i=0; i<ncon; i++){
    //            TPZConnect &newnod = cmesh->ConnectVec()[i];
    //            newnod.SetLagrangeMultiplier(1);
    //        }
    //    }
    
    
    cmesh->SetDefaultOrder(order);

    TPZNullMaterial<STATE> *nullMat = new TPZNullMaterial(matID_pressure, gmesh->Dimension());  
    cmesh_p->InsertMaterialObject(nullMat);

    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    
    return cmesh;
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



int main (){
    //main2DFracVug();
    //mainDarcy3D();
    //mainMixed();
    //mainMixedCT();
    return 0;
}
