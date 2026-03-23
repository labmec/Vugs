#ifndef SOURCES_H
#define SOURCES_H

#include "sources.h"

std::set<int> vugBcIds;
std::set<int> vugIds;

std::set<int> fracBcIds;
std::set<int> fracIds;

TPZGeoMesh* generateGMeshWithPhysTagVec(ReadJson inputData, std::string filename, std::string meshName){

    TPZGmshReader GeometryFine;
    TPZGeoMesh *gmeshFine;
    REAL l = 1.0;
    GeometryFine.SetCharacteristiclength(l);

    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag(4);

    std::map<std::string, int> domainData = inputData.DomainData();
    std::map<std::string, int> vugData = inputData.VugData();
    std::map<std::string, int> fracData = inputData.FracData();
    std::vector<BcData> bcData = inputData.BCInput();

    for (auto bc: bcData){
        dim_name_and_physical_tag[1][bc.name] = bc.matId;
    }

    for (auto info: domainData){
            dim_name_and_physical_tag[2][info.first] = info.second;
    }
    if (!vugData.empty()) {
        for (auto info: vugData){
            dim_name_and_physical_tag[2][info.first] = info.second;
        }
    }
    if (!fracData.empty()) {
        for (auto info: fracData){
            dim_name_and_physical_tag[1][info.first] = info.second;
        }
    }
    
    // Reading mesh
    GeometryFine.SetDimNamePhysical(dim_name_and_physical_tag);
    gmeshFine = GeometryFine.GeometricGmshMesh(filename,nullptr,false);
    gmeshFine->SetName(meshName);
    return gmeshFine;
}

void MeshWithSegmentPhil(ReadJson inputData, TPZGeoMesh *gmesh){
    int ncreated = 0;
    int nels = gmesh->NElements();
    int matid_vug=6;
    int meshdim = gmesh->Dimension();
    int EfractureId = 6; //CHANGE IT

    int VugBcMatId = EFracBcId;

    TPZVec<int64_t> verificador(nels,-1);
  
    for(int el = 0; el < nels; el++){
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) continue;
        int matId = gel->MaterialId();

        if(matId != matid_vug && matId != EfractureId) continue;

        if(verificador[el] != -1) continue;

        TPZStack<int64_t> tocheck;
        tocheck.Push(el);

        while(tocheck.size()){

            int64_t elcheck = tocheck.Pop();

            if(verificador[elcheck] != -1) continue;

            TPZGeoEl *gelcheck = gmesh->Element(elcheck);
            int gelDim = gelcheck->Dimension();
            verificador[elcheck] = VugBcMatId;
            gelcheck->SetMaterialId(VugBcMatId+500);  
            //vugBcIds.insert(mat); //TODO MELHORAR ISSO
            vugIds.insert(VugBcMatId+500);
            vugBcIds.insert(VugBcMatId);

            int lastside   = gelcheck->NSides();
            if (gelDim == 2) lastside--;

            int ncorners = gelcheck->NCornerNodes();
            int firstside = gelcheck->FirstSide(1);

            for(int iside = firstside; iside < lastside; iside++){
                TPZGeoElSide gelside(gelcheck, iside);
                bool hasDarcyNeigh = gelside.HasNeighbour(EMatId);
                bool hasVugBound   = gelside.HasNeighbour(EFracBcId);
                if (gelDim == 1 && !hasDarcyNeigh) DebugStop();
                if(gelDim == 2 && hasDarcyNeigh && !hasVugBound){
                    gelside.Element()->CreateBCGeoEl(iside, VugBcMatId);
                    //vugBcIds.insert(mat);
                    std::cout<<"Vug bc index: "<<vugBcIds.size()<<std::endl;
                    hasVugBound = true;
                }
                if(gelDim == 1 && hasVugBound) DebugStop();
                if(gelDim == 1){
                    TPZGeoElSide neighbour = gelside.Neighbour();
                    while(neighbour != gelside){
                        int neighmatId = neighbour.Element()->MaterialId();
                        if(neighmatId == EMatId){
                            TPZGeoElBC gbc(neighbour,VugBcMatId);
                            std::cout << "Creating neighbor for fracture " << gelcheck->Index() << "\n"; 
                        }
                        neighbour = neighbour.Neighbour();
                    }
                }
                TPZGeoElSide neighbour = gelside.Neighbour();
                TPZGeoEl *neighgel = neighbour.Element();

                if(gelDim == 2 && !hasDarcyNeigh){
                    int64_t neighindex = neighgel->Index();
                    if(verificador[neighindex] == -1)
                        tocheck.Push(neighindex);
                } 
                if(gelDim == 1){
                    for(int side = 0; side < 2; side++){
                        TPZGeoElSide gelside(gelcheck, side);
                        TPZGeoElSide neighFrac = gelside.HasNeighbour(EfractureId);
                        if(neighFrac){
                            TPZGeoEl* gelFrac = neighFrac.Element();
                            int elIndex = gelFrac->Index();
                            if(verificador[elIndex] == -1){
                                tocheck.Push(elIndex);
                            }
                        }
                    }
                }
            }
        }
      // incrementamos o valor de mat, para o seguinte vug.
        VugBcMatId++;
        ncreated++;
    }
    std::cout<<"Vug bc index final: "<<vugBcIds.size()<<std::endl;

}


//TODO Arrumar para fazer Vugs e Frac
void MeshWithSegment(ReadJson inputData, TPZGeoMesh *gmesh){
    int ncreated = 0;
    int nels = gmesh->NElements();
    int matid_vug = 0;
    int matid_frac = 0;
    int matFrac = EFracBcId;
    int mat = EVugBcId;

    std::map<std::string, int> domainData = inputData.DomainData();
    std::map<std::string, int> vugData = inputData.VugData();
    std::map<std::string, int> fracData = inputData.FracData();

    if (!vugData.empty()) {
        for (auto info: vugData){
            matid_vug = info.second;
        }
    }
    if (!fracData.empty()) {
        for (auto info: fracData){
            matid_frac = info.second;
        }
    }

    if (!matid_vug && !matid_frac) DebugStop();
    
    TPZVec<int64_t> verificador(nels,-1);
  
    for(int el = 0; el < nels; el++){
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) continue;

        if(gel->MaterialId() != matid_frac && gel->MaterialId() != matid_vug) continue; 

        if(verificador[el] != -1) continue;

        TPZStack<int64_t> tocheck;
        tocheck.Push(el); 

        while(tocheck.size()){

            int64_t elcheck = tocheck.Pop(); 

            if(verificador[elcheck] != -1) continue; 

            TPZGeoEl *gelcheck = gmesh->Element(elcheck); 
            verificador[elcheck] = mat+500; 
            gelcheck->SetMaterialId(mat+500); 

            // fracBcIds.insert(mat); 
            // fracIds.insert(mat+500); 

            vugBcIds.insert(mat); 
            vugIds.insert(mat+500); 

            int nsides   = gelcheck->NSides();
            int ncorners = gelcheck->NCornerNodes();
            int firstside = nsides - ncorners - 1;

            for(int iside = firstside; iside < nsides; iside++){ 
                
                
                TPZGeoElSide gelside(gelcheck, iside);
                TPZStack<TPZGeoElSide> allneigh;

                gelside.AllNeighbours(allneigh); 
                for(auto neigh: allneigh){
                    TPZGeoEl* gelNeigh = neigh.Element();
                    if(gelNeigh->MaterialId() == EMatId && gelside.Dimension() == gmesh->Dimension()-1){ 
                       gelside.Element()->CreateBCGeoEl(iside, mat);
                    }
                    if(gelNeigh->MaterialId() == matid_vug){ 
                        gelNeigh->SetMaterialId(mat+500); 
                        tocheck.Push(gelNeigh->Index());
                    }
                }  
            }
        }
        mat++;
        ncreated++;
    }
}


//TODO Melhorar
void SetUniqueVugConnect(TPZGeoMesh *gmesh, TPZCompMesh *cmesh){
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    int count = 0;
    int nels = gmesh->NElements();
    TPZVec<int64_t> gelIndex(nels, -1);
    std::map<int, int> matId_connect;
    std::map<int, int> origCon_newCon;

    for (int64_t el = 0; el < cmesh->NElements(); el++){

        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int meshDim = gmesh->Dimension();
        int gelMatId = gel->MaterialId();

        for(int side = 0; side < gel->NCornerNodes(); side++){ // sides associated with vertices
            TPZGeoElSide gelside(gel, side); // node i of gel
            TPZCompElSide celside = gelside.Reference();
            if(origCon_newCon.find(celside.ConnectIndex()) != origCon_newCon.end()){
                count++;
                cel->SetConnectIndex(side, origCon_newCon.at(celside.ConnectIndex()));
            }
        }

        if (fracIds.find(gelMatId) == fracIds.end() && vugIds.find(gelMatId) == vugIds.end() && vugBcIds.find(gelMatId) == vugBcIds.end()) continue; 

        int connIndex = 0;
        auto it = matId_connect.find(gel->MaterialId()); // check if the material ID of the vug already has an associated connect index
        if(it != matId_connect.end()){ // if it has, use the same connect index for all vug elements
            connIndex = matId_connect.at(gel->MaterialId());
        }
        else{ // if it doesn't, create a new connect index and associate it with the material ID of the vug
            connIndex = cel->ConnectIndex(0);
            matId_connect.insert({gel->MaterialId(), connIndex});
        }
        int nsides = gel->NSides();
        int nVertex = gel->NCornerNodes();

        for(int side = 0; side < nVertex; side++){ // sides associated with vertices
            TPZGeoElSide gelside(gel, side); // node i of gel
            TPZCompElSide celside = gelside.Reference();
            int originalConn = celside.ConnectIndex();
            origCon_newCon.insert({originalConn, connIndex});
            cel->SetConnectIndex(side, connIndex);
        }
    }
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();
    std::cout << "Map size " << origCon_newCon.size() << "\n";
}


void insertAtomicMaterials(TPZCompMesh *cmesh, std::set<int> matIdsEls, std::set<int> matIdsBcs, int typeMesh, ReadJson inputData){

    std::vector<BcData> bcData = inputData.BCInput();
    int bcType = 0;
    int bcId = 0;
    TPZManVector<double, 3> bcValue;

    int dim = cmesh->Dimension();
    if (typeMesh == 0){ // if flux mesh
        for (auto iD: matIdsEls) {
            if(iD == EMatId){
                TPZNullMaterial <STATE> *matDarcy = new TPZNullMaterial(iD, dim);
                cmesh->InsertMaterialObject(matDarcy);
                
                TPZFMatrix<STATE> val1(1,1,0.0);
                TPZVec<STATE> val2(1,0.0);

                for (auto bc: bcData){
                    bcId = bc.matId;
                    bcType = bc.type;
                    bcValue = bc.value;

                    TPZBndCond * face = matDarcy->CreateBC(matDarcy, bcId, bcType, val1, val2);
                    cmesh->InsertMaterialObject(face);
                }
                
                for(auto bcId: fracBcIds) {
                    TPZBndCond *faceVug = matDarcy->CreateBC(matDarcy, bcId, 0, val1, val2);
                    cmesh->InsertMaterialObject(faceVug);
                }
                for(auto bcId: vugBcIds) {
                    TPZBndCond *faceVug = matDarcy->CreateBC(matDarcy, bcId, 0, val1, val2);
                    cmesh->InsertMaterialObject(faceVug);
                }
            }            
        }
    }
    else if (typeMesh == 1){ //if pressure mesh
        int elDim;
        for (auto iD: matIdsEls) {
            if (fracIds.find(iD) != fracIds.end()) elDim = 1;
            else elDim = dim;
            TPZNullMaterial <STATE> *matDarcy = new TPZNullMaterial(iD, elDim); 
            cmesh->InsertMaterialObject(matDarcy);
        }
    }

    cmesh->LoadReferences();
}


void GetAtomicIds(TPZGeoMesh *geomesh, std::set<int> &elsId, std::set<int> &bcId){
    int dim = geomesh->Dimension();
    for (auto gel: geomesh->ElementVec()) {
        if (! gel) {
            continue;
        }
        int matId = gel->MaterialId();
        int geldim = gel->Dimension();
        if (geldim == dim || (fracIds.find(matId) != fracIds.end())) {
            elsId.insert(matId);
        }
        if (geldim == dim-1 && (fracIds.find(matId) == fracIds.end())) {
            bcId.insert(matId);
        }
    }
}


TPZCompMesh *CreateFluxMesh(TPZGeoMesh *gmesh, std::set<int> &elsId, std::set<int> &bcId, int &orderp, ReadJson inputData){
    int typeMesh = 0; //Malha de fluxo
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetName("FluxMesh");
     
    insertAtomicMaterials(cmesh, elsId, bcId, typeMesh, inputData);
        
    int pOrder=1;
    cmesh->SetDefaultOrder(pOrder);
    int meshdim = gmesh->Dimension();
    if (orderp==0){
        cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
    }

    //!ASK
    // gmesh->ResetReference();
    // for (auto cel : cmesh->ElementVec()) {
    //         if (!cel) {
    //             continue;
    //         }
    //         if (cel->Reference()->MaterialId() == 300) {
    //             cel->LoadElementReference();
    //             DebugStop();
    //         }
    // }

    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    cmesh->AutoBuild();

    DuplicateConnectFracture(gmesh, cmesh); 
    SideOrientation(cmesh);

    // cmesh->ExpandSolution();
    // cmesh->CleanUpUnconnectedNodes();
    // cmesh->InitializeBlock();

    return cmesh;
}


TPZCompMesh *CreatePressureMesh(TPZGeoMesh *gmesh, std::set<int> &elsId, std::set<int> &bcId,int order, ReadJson inputData){
    int TypeMesh = 1;
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetName("PressureMesh");
    //cmesh->SetDimModel(gmesh->Dimension()-1);
    cmesh->SetDefaultOrder(order);

    if(order>0){
            cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
        }
    else {
            cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
            // cmesh->ApproxSpace().CreateDisconnectedElements(true);
        }
    
    insertAtomicMaterials(cmesh, elsId, bcId, TypeMesh, inputData);
    
    
    cmesh->AutoBuild();

    if(1 > 0){
        int64_t ncon = cmesh->NConnects();
        for(int64_t i=0; i<ncon; i++){
            TPZConnect &newnod = cmesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }
    }

    PrintCompMesh(cmesh);
    SetUniqueVugConnect(gmesh, cmesh);
    
    cmesh->InitializeBlock();
 
    return cmesh;
}


TPZMultiphysicsCompMesh *CreateMultiMesh(TPZGeoMesh* gmesh, TPZVec<TPZCompMesh *> meshvec, ReadJson inputData){

    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetName("MultiMesh");
    if (meshvec.size() != 2) {
        std::cout << "It is expected 2 meshes.\n";
        DebugStop();
    }

    std::vector<BcData> bcData = inputData.BCInput();
    int bcType = 0;
    int bcId = 0;
    TPZManVector<double, 3> bcValue;
    
    // Add materials (weak formulation)
    TPZMixedDarcyFlow *matDarcy = new TPZMixedDarcyFlow(EMatId, inputData.dim());
    cmesh->InsertMaterialObject(matDarcy);
    matDarcy->SetConstantPermeability(inputData.perm());
    
    for(auto fracId: fracIds) { // Frac elements
        auto *matDarcyFrac = new TPZMixedDarcyFlow(fracId, 1);
        cmesh->InsertMaterialObject(matDarcyFrac);
    }

    for(auto vugId: vugIds) { // Vug elements
        auto *matDarcyVug = new TPZMixedDarcyFlow(vugId, inputData.dim());
        cmesh->InsertMaterialObject(matDarcyVug);
    }

    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);

    for (auto bc: bcData){
        bcId = bc.matId;
        bcType = bc.type;
        bcValue = bc.value;

        val2[0] = bcValue[0];
        TPZBndCond * face = matDarcy->CreateBC(matDarcy, bcId, bcType, val1, val2);
        cmesh->InsertMaterialObject(face);
    }

    val2[0] = 0;
    for(auto bcId: fracBcIds) { // Frac boundary elements
        TPZBndCond *faceFrac = matDarcy->CreateBC(matDarcy, bcId, 0, val1, val2);
        cmesh->InsertMaterialObject(faceFrac);
    }

    for(auto bcId: vugBcIds) { // Vug boundary elements 
        TPZBndCond *faceVug = matDarcy->CreateBC(matDarcy, bcId, 0, val1, val2);
        cmesh->InsertMaterialObject(faceVug);
    }

    cmesh->ExpandSolution();
    cmesh->ApproxSpace().Style()= TPZCreateApproximationSpace::EMultiphysics;
    cmesh->BuildMultiphysicsSpace(meshvec);
    cmesh->CleanUpUnconnectedNodes();

    CreateInterfaceGeoEls(gmesh);
    InsertInterfaceEls(cmesh, gmesh); 

    cmesh->InitializeBlock();

    return cmesh;
}


TPZCompMesh *CreateMesh(TPZGeoMesh* gmesh, ReadJson inputData){

    TPZCompMesh *cmesh =  new TPZCompMesh(gmesh);
    cmesh->SetName("CompMesh");

    std::vector<BcData> bcData = inputData.BCInput();
    int bcType = 0;
    int bcId = 0;
    TPZManVector<double, 3> bcValue;
    
    // Add materials (weak formulation)
    TPZDarcyFlow *matDarcy = new TPZDarcyFlow(EMatId, inputData.dim());
    matDarcy->SetConstantPermeability(inputData.perm());
    cmesh->InsertMaterialObject(matDarcy);

    TPZDarcyFlow *matDarcyFrac = nullptr;
    for(auto Id: fracIds) {
        matDarcyFrac = new TPZDarcyFlow(Id, 1);
        if(inputData.problemType() == 0){
            matDarcyFrac->SetConstantPermeability(inputData.permFrac());
        }
        cmesh->InsertMaterialObject(matDarcyFrac);
    }

    TPZDarcyFlow *matDarcyVug = nullptr;
    for(auto Id: vugIds) {
        matDarcyVug = new TPZDarcyFlow(Id, inputData.dim());
        if(inputData.problemType() == 0) {
            matDarcyVug->SetConstantPermeability(inputData.permVug());
        }
        cmesh->InsertMaterialObject(matDarcyVug);
    }

    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);

    for (auto bc: bcData){
        bcId = bc.matId;
        bcType = bc.type;
        bcValue = bc.value;

        val2[0] = bcValue[0];
        TPZBndCond * face = matDarcy->CreateBC(matDarcy, bcId, bcType, val1, val2);
        cmesh->InsertMaterialObject(face);
    }

    // val2[0] = 0;
    // for(auto bcId: fracBcIds) { // Frac boundary elements
    //     TPZBndCond *faceFrac = matDarcy->CreateBC(matDarcy, bcId, 0, val1, val2);
    //     cmesh->InsertMaterialObject(faceFrac);
    // }

    // for(auto bcId: vugBcIds) { // Vug boundary elements 
    //     TPZBndCond *faceVug = matDarcy->CreateBC(matDarcy, bcId, 0, val1, val2);
    //     cmesh->InsertMaterialObject(faceVug);
    // }

    cmesh->AutoBuild();

    PrintCompMesh(cmesh);

    if(inputData.problemType() == 1) SetUniqueVugConnect(gmesh, cmesh);

    cmesh->CleanUpUnconnectedNodes();

    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();

    cmesh->ExpandSolution();

    return cmesh;
}


void CreateInterfaceGeoEls(TPZGeoMesh *gmesh){

    int nEl = gmesh->NElements();
    for(int el = 0; el < nEl; el++){

        TPZGeoEl *gel = gmesh->Element(el);

        if(fracBcIds.find(gel->MaterialId()) == fracBcIds.end() && vugBcIds.find(gel->MaterialId()) == vugBcIds.end()) continue;

        int nSides = gel->NSides();
        TPZGeoElSide gelSide(gel, nSides - 1);

        TPZGeoElSide neighSide = gelSide.HasNeighbour(EVugId);
        if(vugBcIds.find(gelSide.Element()->MaterialId()) != vugBcIds.end()){
            for (auto id: vugIds){
                neighSide = gelSide.HasNeighbour(id); 
                if(neighSide) break;
            }
        }

        if(!neighSide) neighSide = gelSide.HasNeighbour(EFracId);
        if(fracBcIds.find(gelSide.Element()->MaterialId()) != fracBcIds.end()){
            for (auto id: fracIds){
                neighSide = gelSide.HasNeighbour(id); 
                if(neighSide) break;
            }
        }

        if(!neighSide) DebugStop();

        TPZGeoElBC gelInterface(neighSide, ELagrange);
    }
}


void InsertInterfaceEls(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh *gmesh){

    TPZLagrangeMultiplierCS<STATE> *matInterface = new TPZLagrangeMultiplierCS<STATE>(ELagrange, gmesh->Dimension()-1, 1);
    //matInterface->SetMultiplier(-1);
    //matInterface->SetMultiplier(-1);
    cmesh->InsertMaterialObject(matInterface);

    gmesh->ResetReference(); 
    cmesh->LoadReferences();

    int nEl = gmesh->NElements();

    std::set<int> neighIndices; // set to store already analyzed frac boundary elements

    for(int el = 0; el < nEl; el++){

        TPZGeoEl *gel = gmesh->Element(el);
        //TPZCompEl *cel = gel->Reference();
        if(!gel || gel->MaterialId() != ELagrange) continue;

        int nSides = gel->NSides();
        TPZGeoElSide gelSide(gel, nSides - 1);
        TPZStack<TPZGeoElSide> allneigh;

        TPZCompElSide neighPressure = gelSide.HasNeighbour(EVugId).Reference(); //or EFracId
        TPZCompElSide neighHdiv = gelSide.HasNeighbour(EVugBcId).Reference(); //or EFracBcId

        for (auto id: vugIds){
            neighPressure = gelSide.HasNeighbour(id).Reference();
            if(neighPressure) break;
        }

        for (auto id: vugBcIds){
            neighHdiv = gelSide.HasNeighbour(id).Reference();
            if(neighHdiv) break;
        }

        for (auto id: fracIds){
            neighPressure = gelSide.HasNeighbour(id).Reference();
            if(neighPressure) break;
        }

        TPZGeoElSide neighbour = gelSide.Neighbour();
        int neighIndex = 0;
        while(neighbour != gelSide){
            if(fracBcIds.find(neighbour.Element()->MaterialId()) != fracBcIds.end()){
                neighIndex = neighbour.Reference().Element()->Index();
                auto it = neighIndices.find(neighIndex);
                if (it == neighIndices.end()){
                    neighHdiv = neighbour.Reference();
                    neighIndices.insert(neighIndex);
                }
            }
            if(neighHdiv) break;
            neighbour = neighbour.Neighbour();
        }


        if(!neighPressure || !neighHdiv) DebugStop();

        TPZMultiphysicsInterfaceElement *interface = new TPZMultiphysicsInterfaceElement(*cmesh, gel, neighHdiv, neighPressure);
    }
}
void PrintCompMesh(TPZCompMesh *cmesh)
{
    std::cout << "\nPrinting comp mesh in .txt and .vtk formats...\n";

    std::ofstream VTKCompMeshFile(cmesh->Name() + ".vtk");
    std::ofstream TextCompMeshFile(cmesh->Name() + ".txt");

    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, VTKCompMeshFile);
    cmesh->Print(TextCompMeshFile);
}

void PrintGeoMesh(TPZGeoMesh *gmesh)
{
    std::cout << "\nPrinting geo mesh in .txt and .vtk formats...\n";

    std::ofstream VTKGeoMeshFile(gmesh->Name() + ".vtk");
    std::ofstream TextGeoMeshFile(gmesh->Name() + ".txt");

    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, VTKGeoMeshFile);
    gmesh->Print(TextGeoMeshFile);
}


void SideOrientation(TPZCompMesh *cmesh){ //CheckSideOrientation(TPZCompMesh *cmesh, TPZInterpolationSpace *intEl);

    for(int el = 0; el < cmesh->NElements(); el++){
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        
        //if(gel->MaterialId() < EVugBcId || gel->MaterialId() >= EVugId) continue;
        if(gel->MaterialId() != EMatId) continue;
        
        int nSides = gel->NSides(); 
        int nNodeSides = gel->NCornerNodes(); 

        for(int side = nNodeSides; side < nSides-1; side++){
            TPZGeoElSide gelSide(gel, side);
            TPZGeoElSide neigh = gelSide.HasNeighbour(vugBcIds);
            if(!neigh) neigh = gelSide.HasNeighbour(fracBcIds);
            if(neigh){                
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                int orientation = intel->GetSideOrient(side);
                std::cout << orientation << "\n";
                //int orientation = gel->NormalOrientation(side);
                intel->SetSideOrient(side, 1.);
            }
        }
    }
}  


void DuplicateConnectFracture(TPZGeoMesh *gmesh, TPZCompMesh *cmesh){
    int nels = gmesh->NElements();
    std::set<int> neighIndices; // set to store el indeces that has already been analyzed
    for (int64_t el = 0; el < cmesh->NElements(); el++){
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int meshDim = gmesh->Dimension();
        int elDim = gel->Dimension(); 
        int matId = gel->MaterialId();


        if(elDim != meshDim - 1 || matId < EFracBcId) continue; // Searching for frac boundary id = 200
        if(neighIndices.find(gel->Index()) != neighIndices.end()) continue;

        int iside = gel->NSides() - 1;

        TPZGeoElSide gelside(gel, iside);
        TPZCompElSide celside = gelside.Reference();

        TPZGeoElSide neighBCside = gelside.HasNeighbour(EFracBcId); 
        TPZGeoElSide neighDarcyside = gelside.HasNeighbour(EMatId); 
        TPZGeoElSide neighDarcy2side = neighDarcyside.HasNeighbour(EMatId);

        TPZGeoElSide neighbour = gelside.Neighbour();
        int neighIndex = 0;
        while(neighbour != gelside){
            if(fracBcIds.find(neighbour.Element()->MaterialId()) != fracBcIds.end()){
                neighIndex = neighbour.Element()->Index();
                auto it = neighIndices.find(neighIndex);
                if (it == neighIndices.end()){
                    neighBCside = neighbour;
                    neighIndices.insert(neighIndex);
                }
            }
            neighbour = neighbour.Neighbour();
        }


        if (!neighBCside || !neighDarcyside || !neighDarcy2side) DebugStop();


        TPZCompElSide celBCneigh = neighBCside.Reference(); 
        TPZCompElSide celDarcyneigh = neighDarcyside.Reference();
        TPZCompElSide celDarcy2neigh = neighDarcy2side.Reference();

        int connRight = celDarcyneigh.ConnectIndex();
        int connLeft = celDarcy2neigh.ConnectIndex();

        celDarcyneigh.SplitConnect(celDarcy2neigh);

        connRight = celDarcyneigh.ConnectIndex();
        connLeft = celDarcy2neigh.ConnectIndex();

        celside.Element()->SetConnectIndex(0, connRight);
        celBCneigh.Element()->SetConnectIndex(0, connLeft);

        cmesh->ExpandSolution();
        cmesh->ComputeNodElCon();
    }
}


void Solve(TPZLinearAnalysis* an, TPZCompMesh* cmesh, ReadJson inputData){

    // an->LoadSolution();

    #ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE> matMixed(cmesh);
    #else
    TPZFStructMatrix<STATE> matMixed(cmesh);
    #endif
    matMixed.SetNumThreads(0);
    an->SetStructuralMatrix(matMixed);

    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an->SetSolver(step);
    
    an->Assemble();
    an->Solve();
}


void PostProcess(ReadJson inputData){
    return;
}


ReadJson::ReadJson(std::string fileName){
    
    std::ifstream filejson(fileName);

    fInputFile = json::parse(filejson,nullptr,true,true,true); 

    fMeshName = fInputFile["MeshName"];
    
    fMeshDirectory = fInputFile["MeshDirectory"];

    if(fInputFile["DomainData"].find("NameDomain") == fInputFile["DomainData"].end()) DebugStop();

    fDomainData[fInputFile["DomainData"]["NameDomain"]] = fInputFile["DomainData"]["matId"];

    if(fInputFile["DomainData"].find("NameVug") != fInputFile["DomainData"].end()) fVugData[fInputFile["DomainData"]["NameVug"]] = fInputFile["DomainData"]["matIdVug"];

    if(fInputFile["DomainData"].find("NameFrac") != fInputFile["DomainData"].end()) fFracData[fInputFile["DomainData"]["NameFrac"]] = fInputFile["DomainData"]["matIdFrac"];
    
    fApproxType = fInputFile["Simulation"]["ApproximationType"];

    fProblemType = fInputFile["Simulation"]["ProblemType"];
    
    fPresspOrder = fInputFile["PressurepOrder"];

    fDim = fInputFile["Dimension"];
    
    fResolution = fInputFile["Resolution"];
    
    fPerm = fInputFile["ProblemData"]["Permeability"];

    if(fInputFile["ProblemData"].find("VugPermeability") == fInputFile["ProblemData"].end() && fProblemType == 0) DebugStop();

    if(fInputFile["ProblemData"].find("VugPermeability") != fInputFile["ProblemData"].end()) fVugPerm = fInputFile["ProblemData"]["VugPermeability"];

    if(fInputFile["ProblemData"].find("FracPermeability") != fInputFile["ProblemData"].end()) fFracPerm = fInputFile["ProblemData"]["FracPermeability"];

    fVisc = fInputFile["ProblemData"]["FluidViscosity"];

    BcData bcInput;

    for(auto& bcjson : fInputFile["BCs"]){
        if(bcjson.find("Name") == bcjson.end()) DebugStop(); // check if the information exists

        bcInput.name = bcjson["Name"];
        bcInput.matId = bcjson["matId"];
        bcInput.type = bcjson["Type"];
        bcInput.value[0] = bcjson["Value"];
        
        fBcDataVec.push_back(bcInput);
    }
}
    
std::string ReadJson::MeshName(){

    return fMeshName;
}

std::string ReadJson::MeshFile(){

    std::string meshPath = fMeshDirectory + fMeshName + ".msh";
    return meshPath;
}

std::map<std::string, int> ReadJson::DomainData(){
    return fDomainData;
}

std::map<std::string, int> ReadJson::VugData(){
    return fVugData;
}

std::map<std::string, int> ReadJson::FracData(){
    return fFracData;
}

int ReadJson::approxType(){
    return fApproxType;
}

int ReadJson::problemType(){
    return fProblemType;
}

int ReadJson::pressOrder(){
    return fPresspOrder;
}

int ReadJson::dim(){
    return fDim;
}

int ReadJson::resolution(){
    return fResolution;
}

double ReadJson::perm(){
    return fPerm;
}

double ReadJson::permVug(){
    return fVugPerm;
}

double ReadJson::permFrac(){
    return fFracPerm;
}
    
double ReadJson::visc(){
    return fVisc;
}

std::vector<BcData> ReadJson::BCInput(){
    return fBcDataVec;
}

#endif
