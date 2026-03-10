#ifndef SOURCES_H
#define SOURCES_H

#include "sources.h"
#include "pzintel.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZHDivApproxCreator.h"
#include "TPZLagrangeMultiplierCS.h"

TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine){

    TPZGmshReader  GeometryFine;
    TPZGeoMesh *gmeshFine;
    REAL l = 1.0;
    GeometryFine.SetCharacteristiclength(l);
    
    // Reading mesh
    GeometryFine.SetDimNamePhysical(dim_name_and_physical_tagFine);
    gmeshFine = GeometryFine.GeometricGmshMesh(filename,nullptr,false);
    return gmeshFine;
}


void MeshWithSegmentVugs(TPZGeoMesh *gmesh, std::set<int> &vugBcIds,std::set<int> &vugIds){
    int ncreated = 0;
    int nels = gmesh->NElements();
    int matid_vug=6;
    int matDarcytag=1;
    int meshdim=gmesh->Dimension();
    int matid_Vugbound=100;

    int mat = 100;

    TPZVec<int64_t> verificador(nels,-1);
  
    for(int el = 0; el < nels; el++){
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) continue;

        if(gel->MaterialId() != matid_vug) continue;

        if(verificador[el] != -1) continue;

        TPZStack<int64_t> tocheck;
        tocheck.Push(el);

        while(tocheck.size()){

            int64_t elcheck = tocheck.Pop();

            if(verificador[elcheck] != -1) continue;

            TPZGeoEl *gelcheck = gmesh->Element(elcheck);
            verificador[elcheck] = mat;
            gelcheck->SetMaterialId(mat+500);  
            //vugBcIds.insert(mat); //TODO MELHORAR ISSO
            vugIds.insert(mat+500);
            vugBcIds.insert(mat);

            int nsides   = gelcheck->NSides();
            int ncorners = gelcheck->NCornerNodes();
            int firstside = nsides - ncorners - 1;

            for(int iside = firstside; iside < nsides; iside++){
                TPZGeoElSide gelside(gelcheck, iside);
                bool hasDarcyNeigh = gelside.HasNeighbour(matDarcytag);
                bool hasVugBound   = gelside.HasNeighbour(matid_Vugbound);
                if(hasDarcyNeigh && !hasVugBound){
                    gelside.Element()->CreateBCGeoEl(iside, mat);
                    //vugBcIds.insert(mat);
                    std::cout<<"Vug bc index: "<<vugBcIds.size()<<std::endl;
                    hasVugBound = true;
                }
                TPZGeoElSide neighbour = gelside.Neighbour();
                TPZGeoEl *neighgel = neighbour.Element();

                if(hasVugBound){
                    neighgel->SetMaterialId(mat);
                }

                else if(!hasDarcyNeigh){
                    int64_t neighindex = neighgel->Index();
                    if(verificador[neighindex] == -1)
                        tocheck.Push(neighindex);
                } 
            }
        }
      // incrementamos o valor de mat, para o seguinte vug.
        mat++;
        ncreated++;
    }
    std::cout<<"Vug bc index final: "<<vugBcIds.size()<<std::endl;

}

void SetUniqueVugConnect(TPZGeoMesh *gmesh, TPZCompMesh *cmesh){
    int nels = gmesh->NElements();
    //TPZVec<int64_t> vugIndex(nVugs, -1);
    TPZVec<int64_t> gelIndex(nels, -1);
    std::map<int, int> matId_connect;

    for (int64_t el = 0; el < cmesh->NElements(); el++){

        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int meshDim = gmesh->Dimension();

        //if (gel->Dimension() != meshDim-1) continue; // only boundary elements (Dim-1 elems)
        if (gel->MaterialId() < EVugId) continue;

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
            TPZStack<TPZGeoElSide> allneigh; // all node neighbors of node i
            gelside.AllNeighbours(allneigh);
            int nneighs = allneigh.size();
            cel->SetConnectIndex(side, connIndex);
            // for(int neigh = 0; neigh < nneighs; neigh++){
            //     TPZGeoElSide neighside = allneigh[neigh];
            //     TPZGeoEl *gelneigh = neighside.Element();
            //     if (gelneigh->MaterialId() == EMatId || gelneigh->MaterialId() == EVugBcId) continue; // ignore neighbors that are not vug elements
            //     TPZCompEl *celneigh = gelneigh->Reference();
            //     celneigh->SetConnectIndex(neighside.Side(), connIndex);
            //}
        }
    }
}

void insertAtomicMaterials(TPZCompMesh *cmesh, std::set<int> matIdsVol, std::set<int> matIdsBcs, int typeMesh){
    int dim = cmesh->Dimension();
    if (typeMesh==0){//Se for malha de fluxo não insertar vugs
        for (auto iD:matIdsVol) {
            if(iD<499){
            TPZNullMaterial <STATE> *matDarcy = new TPZNullMaterial(iD, dim);
            cmesh->InsertMaterialObject(matDarcy);
                int bc_id=2;
                int bc_typeN = 1;
                int bc_typeD = 0;
                TPZFMatrix<STATE> val1(1,1,0.0);
                TPZVec<STATE> val2(1,0.0);
                int dim2d=2;

                TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,EbcNoFlux,bc_typeN,val1,val2);
                cmesh->InsertMaterialObject(face2);
                
                val2[0]=100; // Valor a ser impuesto como presión en la entrada
                TPZBndCond * face = matDarcy->CreateBC(matDarcy,EbcInletId,bc_typeD,val1,val2);
                cmesh->InsertMaterialObject(face);

                val2[0]=10; // Valor a ser impuesto como presión en la salida
                TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,EbcOutletId,bc_typeD,val1,val2);
                cmesh->InsertMaterialObject(face1);
                //val2[0]=30;
//                for(auto bcId: vugBcIds) {
//                    val2[0]=30;
//                    TPZBndCond *faceVug = matDarcy->CreateBC(matDarcy,bcId,bc_typeD,val1,val2);
//                    cmesh->InsertMaterialObject(faceVug);
//                }
                //int PContornoVug=-50;
                //val2[0]=PContornoVug;
                //TPZBndCond *faceVug = matDarcy->CreateBC(matDarcy,100,bc_typeD,val1,val2);
                //cmesh->InsertMaterialObject(faceVug);
                //TPZBndCond *faceVug1 = matDarcy->CreateBC(matDarcy,101,bc_typeD,val1,val2);
                //cmesh->InsertMaterialObject(faceVug1);
                //TPZBndCond *faceVug2 = matDarcy->CreateBC(matDarcy,102,bc_typeD,val1,val2);
                //cmesh->InsertMaterialObject(faceVug2);
                //TPZBndCond *faceVug3 = matDarcy->CreateBC(matDarcy,103,bc_typeD,val1,val2);
                //cmesh->InsertMaterialObject(faceVug3);
                
                
            }
            else if(iD>499)continue;

            //
                        
        }
        
        for (auto iD:matIdsBcs) {
//            TPZNullMaterial <STATE> *matDarcy = new TPZNullMaterial(iD, dim);
//            cmesh->InsertMaterialObject(matDarcy);
            if (!iD){
            TPZNullMaterial<STATE> * face2 = new TPZNullMaterial(iD, dim-1);
            cmesh->InsertMaterialObject(face2);
            std::cout<<"Bc Flux: "<<iD<<std::endl;
            }

        }
    }
    else if (typeMesh==1){//Se for malha de pressão
        for (auto iD:matIdsVol) {
            if(iD<599){
            std::cout<<"Material: "<<iD<<std::endl;
            TPZNullMaterial <STATE> *matDarcy = new TPZNullMaterial(iD, dim);
            cmesh->InsertMaterialObject(matDarcy);
            }
        }
//        for (auto iD:matIdsBcs) {
//            if(iD<99){
//            std::cout<<"Material: "<<iD<<std::endl;
//
//            TPZNullMaterial<STATE> * face2 = new TPZNullMaterial(iD, dim-1);
//            cmesh->InsertMaterialObject(face2);
//            }
//            else if(iD>99)continue;
//
//        }
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

//TODO VERIFY
TPZCompMesh *CreateFluxMesh(TPZGeoMesh *gmesh, std::set<int> &volId, std::set<int> &bcId, int &orderp){
    int dim2d = 2;
    int typeMesh=0;//Malha de fluxo
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetName("FluxMesh");
     
    //GetAtomicIds(gmesh, volId, bcId);
    insertAtomicMaterialsf(cmesh, volId, bcId);
    
    //int pOrder=1;
    int pOrder=1;
    cmesh->SetDefaultOrder(pOrder);
    int meshdim = gmesh->Dimension();
    if (orderp==0){
        cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
    }
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    cmesh->AutoBuild();
    //SetUniqueVugConnect(gmesh, cmesh);
    cmesh->InitializeBlock();
    //std::cout<<cmesh->NEquations() <<std::endl;
    return cmesh;
}

//TODO VERIFY
TPZCompMesh *CreatePressureMesh(TPZGeoMesh *gmesh, std::set<int> &volId, std::set<int> &bcId,int order){
    int TypeMesh=1;//Malha de pressão
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
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
        }
    
    insertAtomicMaterialsp(cmesh, volId, bcId);
    
    
    cmesh->AutoBuild();

    if(1 > 0){
        int64_t ncon = cmesh->NConnects();
        for(int64_t i=0; i<ncon; i++){
            TPZConnect &newnod = cmesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(2);
        }
    }

    //SetUniqueVugConnect(gmesh, cmesh);
    cmesh->InitializeBlock();
 
    return cmesh;
}

void CreateInterfaceGeoEls(TPZGeoMesh *gmesh){
 
    int nEl = gmesh->NElements();
    for(int el = 0; el < nEl; el++){

        TPZGeoEl *gel = gmesh->Element(el);

        if(!gel || gel->MaterialId() < EVugBcId || gel->MaterialId() >= EVugId) continue;

        int nSides = gel->NSides();
        TPZGeoElSide gelSide(gel, nSides - 1);

        TPZGeoElSide neighSide = gelSide.HasNeighbour(EVugId);

        for (auto id: vugIds){
            neighSide = gelSide.HasNeighbour(id); //TODO VERIFY
            if(neighSide) break;
        }

        if(!neighSide) DebugStop();

        TPZGeoElBC gelInterface(neighSide, ELagrange);
    }
}

void InsertInterfaceEls(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh *gmesh){

    TPZLagrangeMultiplierCS<STATE> *matInterface = new TPZLagrangeMultiplierCS<STATE>(ELagrange, gmesh->Dimension()-1, 1);
    cmesh->InsertMaterialObject(matInterface);

    gmesh->ResetReference(); //! ASK
    cmesh->LoadReferences();

    int nEl = gmesh->NElements();

    for(int el = 0; el < nEl; el++){

        TPZGeoEl *gel = gmesh->Element(el);
        //TPZCompEl *cel = gel->Reference();
        if(!gel || gel->MaterialId() != ELagrange) continue;

        int nSides = gel->NSides();
        TPZGeoElSide gelSide(gel, nSides - 1);

        TPZCompElSide neighVug = gelSide.HasNeighbour(EVugId).Reference();
        TPZCompElSide neighHdiv = gelSide.HasNeighbour(EVugBcId).Reference();

        for (auto id: vugIds){
            neighVug = gelSide.HasNeighbour(id).Reference();
            if(neighVug) break;
        }

        for (auto id: vugBcIds){
            neighHdiv = gelSide.HasNeighbour(id).Reference();
            if(neighHdiv) break;
        }

        if(!neighVug || !neighHdiv) DebugStop();

        TPZMultiphysicsInterfaceElement *interface = new TPZMultiphysicsInterfaceElement(*cmesh, gel, neighHdiv, neighVug);
    }
}

void PrintCompMesh(TPZCompMesh *cmesh)
{
    std::cout << "\nPrinting multiphysics mesh in .txt and .vtk formats...\n";

    std::ofstream VTKCompMeshFile(cmesh->Name() + ".vtk");
    std::ofstream TextCompMeshFile(cmesh->Name() + ".txt");

    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, VTKCompMeshFile);
    cmesh->Print(TextCompMeshFile);
}
void insertAtomicMaterialsf(TPZCompMesh *cmesh, std::set<int> matIdsVol, std::set<int> matIdsBcs){
    
    int dim = cmesh->Dimension();
    
    for (auto iD:matIdsVol) {
        if(iD==1){
        TPZNullMaterial <STATE> *matDarcy = new TPZNullMaterial(iD, dim);
        cmesh->InsertMaterialObject(matDarcy);
        
        int bc_id=2;
        int bc_typeN = 1;
        int bc_typeD = 0;
        TPZFMatrix<STATE> val1(1,1,0.0);
        TPZVec<STATE> val2(1,0.0);
        int dim2d=2;

        TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,EbcNoFlux,bc_typeN,val1,val2);
        cmesh->InsertMaterialObject(face2);
        
        val2[0]=100; // Valor a ser impuesto como presión en la entrada
        
        TPZBndCond * face = matDarcy->CreateBC(matDarcy,EbcInletId,bc_typeD,val1,val2);
        cmesh->InsertMaterialObject(face);
      
        val2[0]=10; // Valor a ser impuesto como presión en la salida
        TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,EbcOutletId,bc_typeD,val1,val2);

        cmesh->InsertMaterialObject(face1);
//            int PContornoVug=50;
//            val2[0]=PContornoVug;
//            TPZBndCond *faceVug = matDarcy->CreateBC(matDarcy,151,bc_typeD,val1,val2);
//            cmesh->InsertMaterialObject(faceVug);
//            for(auto bcId: vugBcIds) {
//                if(99<bcId<150){
//                std::cout<<"Bc index: "<<bcId<<std::endl;
//                int PContornoVug=30;
//                val2[0]=PContornoVug;
//                TPZBndCond *faceVug = matDarcy->CreateBC(matDarcy,bcId,bc_typeD,val1,val2);
//                cmesh->InsertMaterialObject(faceVug);
//                }
//            }
            for(int i=0;i<44;i++){
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
                std::cout<<"Bc index: "<<bcId<<std::endl;
                val2[0]=PContornoVug;
                TPZBndCond *faceVug = matDarcy->CreateBC(matDarcy,bcId,bc_typeD,val1,val2);
                cmesh->InsertMaterialObject(faceVug);
               
            }
    
        }
      
    }
    for (auto iD:matIdsBcs) {
        
        TPZNullMaterial<STATE> * face2 = new TPZNullMaterial(iD, dim-1);
        //cmesh->InsertMaterialObject(face2);

    }
}
void insertAtomicMaterialsp(TPZCompMesh *cmesh, std::set<int> matIdsVol, std::set<int> matIdsBcs){
    
    int dim = cmesh->Dimension();
    for (auto iD:matIdsVol) {
        TPZNullMaterial <STATE> *matDarcy = new TPZNullMaterial(iD, dim);
        cmesh->InsertMaterialObject(matDarcy);
        
    }
    for (auto iD:matIdsBcs) {
        
        TPZNullMaterial<STATE> * face2 = new TPZNullMaterial(iD, dim-1);
        //cmesh->InsertMaterialObject(face2);

    }
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
            if(neigh){
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                int orientation = intel->GetSideOrient(side);
                //std::cout << orientation << "\n";
                //int orientation = gel->NormalOrientation(side);
                intel->SetSideOrient(side, 1.);
                //std::cout << intel->GetSideOrient(side) << "\n";
            }
        }
    }
}

#endif
