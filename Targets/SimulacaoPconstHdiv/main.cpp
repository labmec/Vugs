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

    
    std::ofstream file20("TestGeoMesh2D.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file20);
    
    int nels = gmesh->NElements();
    int nElVugBound = 0;
    TPZVec<int64_t> els_cont1d(nels,-1); // vector to store the mat ids of the 1d elements (contours), the rest will be -1
    int nVugs = CreateBoundaryElements(gmesh, els_cont1d); //MeshWithSegmentVugs(gmesh);

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

    cmesh_mult->InitializeBlock();
    std::cout<<cmesh_mult->Element(1)<<std::endl;

    InsertInterfaceEls(cmesh_mult, gmesh); //TODO
        
    bool mustOp = false;


    TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh_mult);

    Analisys->LoadSolution();

    TPZStepSolver<STATE> step;
          
    step.SetDirect(ELDLt);
        
    Analisys->SetSolver(step);
          
    //Ensamblaje de la matriz de rigidez y vector de carga
    Analisys->Assemble();

    //TODO Não entendi
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

    int ref = 0;
    std::string file_reservoir("SolVictorCTmesh.vtk");
    Analisys->DefineGraphMesh(dim2d,scalnames,vecnames,file_reservoir);

    Analisys->PostProcess(ref, dim2d);
}

//------------------------------------------------------------------------------------------------------------------

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
    // GetAtomicIds(gmesh, volId, bcId);
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

void CreateCompMeshMP(TPZGeoMesh *gmesh, TPZManVector<TPZCompMesh *,2> &cmeshes){

    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDimModel(gmesh->Dimension());
    // cmesh->SetDefaultOrder(1);
    cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;

    TPZMixedDarcyFlow *mat = new TPZMixedDarcyFlow(EMatId, gmesh->Dimension());  
    cmesh->InsertMaterialObject(mat);


}
