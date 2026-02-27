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
#include "TPZLagrangeMultiplierCS.h"
#include "pzstepsolver.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"
#include <set>
#include "TPZAnalyticSolution.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TSFMixedDarcy.h"
#include "TPZHDivApproxCreator.h"
//hola

enum MatID{
    EMatId = 1,
    EVugId = 100,
    EVugBcId = 500,
    EbcInletId = 2,
    EbcOutletId = 3,
    EbcNoFlux = 4,
    ELagrange = 5,
    EInterface = 6,
    EInterfaceR = 7
};

TPZCompMesh *CreateCompMeshFlux(TPZGeoMesh *gmesh, int pOrder, int nVugs);
TPZCompMesh *CreateCompMeshPressure(TPZGeoMesh *gmesh, int pOrder);
void CreateCompMeshMP(TPZGeoMesh *gmesh, TPZManVector<TPZCompMesh *,2> &cmeshes);

int CreateBoundaryElements(TPZGeoMesh *gmesh, TPZVec<int64_t> &els_cont1d);
void CreateInterfaceGeoEls(TPZGeoMesh *gmesh);
void InsertInterfaceEls(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh *gmesh);
void SetMaterialIdVug(TPZGeoMesh *gmesh);
void SetUniqueVugConnect(TPZGeoMesh *gmesh, TPZCompMesh *cmesh);

//Function to generate a mesh using gmsh library
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void findElDim(TPZStack<TPZGeoElSide> &allneigh, int dim, TPZStack<TPZGeoElSide> &allneighdim);

void Hdiv_MixedCT();

//---------------------------MAIN-----------------------------------
int main (){

    Hdiv_MixedCT();
    return 0;
}
//-------------------------------------------------------------------

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

int CreateBoundaryElements(TPZGeoMesh *gmesh, TPZVec<int64_t> &els_cont1d){

    int nels = gmesh->NElements();
    int nElVugBound = 0;
    int nVug = 0; 

    for (int iel = 0; iel< nels; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        if (!gel){
            continue;
        }
        if (gel->Dimension() != 2) {
            continue;
        }
        int nsides= gel->NSides();
        int ncorners= gel->NCornerNodes();
        int firstside= nsides-ncorners-1;

        for (int iside = firstside; iside<nsides; iside++) {
            TPZGeoElSide gelside(gel, iside);
            int matid = gelside.Element()->MaterialId();
            TPZStack<TPZGeoElSide> allneigh;
            gelside.AllNeighbours(allneigh);
            int nneighs = allneigh.size();
            //verify Dimension
            int verify =0;

            for (int ineigh=0; ineigh<nneighs; ineigh++) {
                TPZGeoEl *gelneigh = allneigh[ineigh].Element();
                int dimen = gelneigh->Dimension();

                if (dimen == 1) {
                    verify = 1;
                }
            }
            if(verify == 1){
                continue;
            }

            for (int ineigh=0; ineigh<nneighs; ineigh++) {
                TPZGeoEl *gelneigh = allneigh[ineigh].Element();
                int matNeigh = gelneigh->MaterialId();
                if (matNeigh != matid && (gel->Dimension() == gelneigh->Dimension()) ) {
                    gelside.Element()->CreateBCGeoEl(iside, EVugBcId);
                    nElVugBound++;
                }
            }
        }
    }

    std::cout<< "Created " << nElVugBound << " boundary elements." << std::endl;

    gmesh->BuildConnectivity();
    const int nels2 = gmesh->NElements();
    //TPZVec<int64_t> els_cont1d(nels2,-1); // vector to store the mat ids of the 1d elements (contours), the rest will be -1
    els_cont1d.Resize(nels2, -1);
    TPZVec<int> verificador(nels2, 0);
    int mat = EVugBcId; 

    for (int iel =nels-1; iel<nels2; iel++) {

        TPZGeoEl * gel = gmesh->Element(iel);
        if (gel->MaterialId() == EVugBcId) {
            std::cout<<"ok "<<std::endl;
        }
        if (!gel) {
            continue;
        }
        if (gel->Dimension() != 1) {
            continue;
        }
        if (verificador[iel]==1) {
            continue;
        }
        if (gel->MaterialId() != EVugBcId) {
            continue;
        }
        int side = 1;

        TPZGeoElSide gelside(gel, side);
        TPZStack<TPZGeoElSide> allneigh;
        gelside.AllNeighbours(allneigh);
        TPZStack<TPZGeoElSide> allneighdim;
        findElDim(allneigh, 1, allneighdim);
        int ntest = allneighdim.size();
        TPZGeoElSide gelneigh = allneighdim[0];
        gel->SetMaterialId(mat);
        els_cont1d[iel]=mat; 
        while (gel != gelneigh.Element()) {
            TPZStack<TPZGeoElSide> allneigh;
            int sidetest = gelneigh.Side();
            if (sidetest==0) {
                gelneigh.SetSide(1);
            }
            else{
                gelneigh.SetSide(0);
            }
            gelneigh.AllNeighbours(allneigh);
            TPZStack<TPZGeoElSide> allneighdim;
            findElDim(allneigh, 1, allneighdim);
            int indexneig = gelneigh.Element()->Index();
            verificador[indexneig] =1;
            gelneigh.Element()->SetMaterialId(mat);
            els_cont1d[indexneig]=mat; 
            gelneigh =allneighdim[0];
        }
        nVug++; //TODO Based on the fact that each vug has its unique boundary
        mat++;
        int ok=0;
    }
    return nVug;
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
        if (gel->MaterialId() < EVugId || gel->MaterialId() > EVugBcId) continue;  //TODO melhorar isso

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
            for(int neigh = 0; neigh < nneighs; neigh++){
                TPZGeoElSide neighside = allneigh[neigh];
                TPZGeoEl *gelneigh = neighside.Element();
                // if (gelneigh->MaterialId() == EVugId) continue; // ignore neighbors that are part of the vug itself
                TPZCompEl *celneigh = gelneigh->Reference();
                celneigh->SetConnectIndex(neighside.Side(), connIndex);
            }
        }
    }
}

void SetMaterialIdVug(TPZGeoMesh *gmesh, TPZVec<int64_t> &els_cont1d){

    TPZVec<int64_t> AllVugsEls2d; 
    int nels = gmesh->NElements();
    // GROUP ALL 2D ELEMENTS OF VUGS (MATID=6) IN AllVugsEls2d
    for(int ind=0;ind<nels;ind++){ 
        TPZGeoEl *gEl=gmesh->Element(ind);
        if (!gEl) {
            continue;
        }
        if (gEl->Dimension() != 2) {
            continue;
        }
        if(gEl->MaterialId()!=EVugId){
            continue;
        }
        AllVugsEls2d.push_back(ind); // store the index of the element if it has matid=6 (Vugs)
    }
    std::cout<<"All vugs 2D size: "<<AllVugsEls2d.size()<<std::endl;
    int ElVugs_identified=0;

    // Identifies 2D vug elements (void/cavity elements) that lie on boundaries by detecting which 
    // ones have 1D neighbor elements, then assigns them material IDs (+500) that correspond to their associated boundary contour groups.
    for(int ind=0;ind<AllVugsEls2d.size();ind++){ // loop over all 2D vug elements
        TPZGeoEl *gEl=gmesh->Element(AllVugsEls2d[ind]);
        if (!gEl) {
            continue;
        }
        if (els_cont1d[AllVugsEls2d[ind]]!=-1) {  
            continue;
        }
        int sides=gEl->NSides();  // return the number of connectivities of the element, which is the number of sides
        for(int side = 3; side < sides; side++) {  // triangular elements, skips the corners side.
            TPZGeoElSide gelside(gEl, side);
            TPZStack<TPZGeoElSide> allneigh;
            gelside.AllNeighbours(allneigh);
            TPZStack<TPZGeoElSide> allneighdim, allneighdim2D;
            findElDim(allneigh, 1, allneighdim); // Find 1D neighbors of the current side and store them in allneighdim
            findElDim(allneigh, 2, allneighdim2D); // Find 2D neighbors and store them in allneighdim2D

            if(allneighdim.NElements()>0){ // it indicates this side (edge) touches a 1D boundary element.
                TPZGeoElSide gelneigh = allneighdim[0]; // Get the first 1D neighbor (should be only one since it's a boundary)
                auto matid1=gelneigh.Element()->MaterialId();
                std::cout << "  Lado " << side << " → Elemento " << AllVugsEls2d[ind]
                          << " (matID=" << matid1 << ")" << std::endl;
                gEl->SetMaterialId(EVugId+matid1); // Assign a new material ID to the 2D vug element based on the material ID of the neighboring 1D boundary element.
                els_cont1d[AllVugsEls2d[ind]] = EVugId+matid1;  // Fix: índice real
                ElVugs_identified++;
                break;  // Solo un contorno por elemento
            }
        }
    }

    // Propagates boundary contour assignments from already-identified 2D vug elements 
    // to their unassigned neighbors through an iterative flood-fill process, ensuring 
    // all connected vug elements within the same cavity receive matching material IDs.
    bool expanded = true;
    while(expanded) {
        expanded = false;
        for(int ind=0; ind<AllVugsEls2d.size(); ind++) { 
            int64_t elIdx = AllVugsEls2d[ind]; 
            if(els_cont1d[elIdx] != -1) continue; // only process unassigned 2D vug elements
            
            TPZGeoEl *gEl = gmesh->Element(elIdx); 
            int sides = gEl->NSides();
            for(int side=3; side<sides; side++) { 
                TPZGeoElSide gelside(gEl, side); 
                TPZStack<TPZGeoElSide> allneigh;
                gelside.AllNeighbours(allneigh); 
                TPZStack<TPZGeoElSide> neigh2D;
                findElDim(allneigh, 2, neigh2D); 
                
                for(int n=0; n<neigh2D.NElements(); n++) { 
                    int64_t nidx = neigh2D[n].Element()->Index(); 
                    if(els_cont1d[nidx] > EVugId) { 
                        els_cont1d[elIdx] = els_cont1d[nidx];  
                        gEl->SetMaterialId(els_cont1d[elIdx]); 
                        ElVugs_identified++;
                        expanded = true;
                        break; // assign the same contour ID as the neighbor and mark as expanded
                    }
                }
                if(expanded) break; // If the element was assigned a contour ID, do not process other sides
            if(expanded) break;
            }
        }
    }

    std::cout << "TOTAL ElVugs_identified: " << ElVugs_identified << std::endl;
}

void CreateInterfaceGeoEls(TPZGeoMesh *gmesh){
 
    int nEl = gmesh->NElements();
    for(int el = 0; el < nEl; el++){

        TPZGeoEl *gel = gmesh->Element(el);

        if(!gel || gel->MaterialId() < EVugBcId) continue;

        int nSides = gel->NSides();
        TPZGeoElSide gelSide(gel, nSides - 1);
        TPZGeoElSide neighSide = gelSide.HasNeighbour(EVugId); //TODO change this because which vug has its own matId

        if(!neighSide) DebugStop();

        TPZGeoElBC gelInterface(neighSide, EInterface);
    }
}

void InsertInterfaceEls(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh *gmesh){

    TPZLagrangeMultiplierCS<STATE> *matInterface = new TPZLagrangeMultiplierCS<STATE>(EInterface, gmesh->Dimension()-1, 1);
    cmesh->InsertMaterialObject(matInterface);

    int nEl = gmesh->NElements();
    for(int el = 0; el < nEl; el++){

        TPZGeoEl *gel = gmesh->Element(el);
        //TPZCompEl *cel = gel->Reference();
        if(!gel || gel->MaterialId() != EInterface) continue;

        int nSides = gel->NSides();
        TPZGeoElSide gelSide(gel, nSides - 1);

        TPZCompElSide neighVug = gelSide.HasNeighbour(EVugId).Reference();
        TPZCompElSide neighHdiv = gelSide.HasNeighbour(EVugBcId).Reference();

        if(!neighVug || !neighHdiv) DebugStop();

        TPZMultiphysicsInterfaceElement *interface = new TPZMultiphysicsInterfaceElement(*cmesh, gel, neighVug, neighHdiv);
    }
}


//-------------------------------------------------------------------------------------------------

void insertAtomicMaterials(TPZCompMesh *cmesh, std::set<int> matIdsVol, std::set<int> matIdsBcs, int typeMesh){
    int dim = cmesh->Dimension();
    if (typeMesh==0){//Se for malha de fluxo não insertar vugs
        for (auto iD:matIdsVol) {
            if(iD<99){
            TPZNullMaterial <STATE> *matDarcy = new TPZNullMaterial(iD, dim);
            cmesh->InsertMaterialObject(matDarcy);
            }
            else if(iD>99)continue;
            
        }
        for (auto iD:matIdsBcs) {
            if(iD<99){
            TPZNullMaterial<STATE> * face2 = new TPZNullMaterial(iD, dim-1);
            cmesh->InsertMaterialObject(face2);
            }
            else if(iD>99)continue;

        }
    }
    else if (typeMesh==1){//Se for malha de pressão
        for (auto iD:matIdsVol) {
            TPZNullMaterial <STATE> *matDarcy = new TPZNullMaterial(iD, dim);
            cmesh->InsertMaterialObject(matDarcy);
          
        }
        for (auto iD:matIdsBcs) {
            
            TPZNullMaterial<STATE> * face2 = new TPZNullMaterial(iD, dim-1);
            cmesh->InsertMaterialObject(face2);

        }
    
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

TPZCompMesh *CreateFluxMesh(TPZGeoMesh *gmesh, std::set<int> &volId, std::set<int> &bcId){
    int matId=1;
    int dim2d = 2;
    int typeMesh=0;
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
     
    //GetAtomicIds(gmesh, volId, bcId);
    insertAtomicMaterials(cmesh, volId, bcId,typeMesh);
        
    //int pOrder=1;
    int pOrder=1;
    cmesh->SetDefaultOrder(pOrder);
    int meshdim = gmesh->Dimension();

    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    //std::cout<<cmesh->NEquations() <<std::endl;
    return cmesh;
}

TPZCompMesh *CreatePressureMesh(TPZGeoMesh *gmesh, std::set<int> &volId, std::set<int> &bcId,int order){
    int TypeMesh=1;
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    //GetAtomicIds(gmesh, volId, bcId);
    insertAtomicMaterials(cmesh, volId, bcId,TypeMesh);
    SetUniqueVugConnect(gmesh, cmesh);

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
}

void Hdiv_MixedCT(){
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZGeoMesh *gmeshp = new TPZGeoMesh;

    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
    dim_name_and_physical_tagCoarse[2]["k11"] = EMatId;
    dim_name_and_physical_tagCoarse[2]["Vugs"] = EVugId;
    dim_name_and_physical_tagCoarse[1]["inlet"] = EbcInletId;
    dim_name_and_physical_tagCoarse[1]["outlet"] = EbcOutletId;
    dim_name_and_physical_tagCoarse[1]["noflux"] = EbcNoFlux;

    
    //std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/testskelSLICE77SP.msh";
    //std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/FewVugsMesh.msh";
    std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/FastMesh.msh";
    //std::string filename="/home/itopo/Stokes-Darcy_Research/Vugs/testskelSLICE77SP.msh";

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);

    
    std::ofstream file20("TestGeoMesh2D.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file20);
    
    int nels = gmesh->NElements();
    int nElVugBound = 0;
    TPZVec<int64_t> els_cont1d(nels,-1); // vector to store the mat ids of the 1d elements (contours), the rest will be -1
    int nVugs = CreateBoundaryElements(gmesh, els_cont1d); //MeshWithSegmentVugs(gmesh);

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
    cmesh_mult->InsertMaterialObject(matDarcy);
    matDarcy->SetConstantPermeability(0.01);
    
    int bc_typeN = 1;
    int bc_typeD = 0;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);
    int dim2d=2;
    int dim2d=2;

    TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,EbcNoFlux,bc_typeN,val1,val2);
    cmesh_mult->InsertMaterialObject(face2);

    val2[0]=100; 
    TPZBndCond * face = matDarcy->CreateBC(matDarcy,EbcInletId,bc_typeD,val1,val2);
    cmesh_mult->InsertMaterialObject(face);

    val2[0]=14; // Valor a ser impuesto como presión en la salida
    TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,EbcOutletId,bc_typeD,val1,val2);
    cmesh_mult->InsertMaterialObject(face1);

    cmesh_mult->ExpandSolution();
    cmesh_mult->ApproxSpace().Style()= TPZCreateApproximationSpace::EMultiphysics;
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

TPZCompMesh *CreateCompMeshFlux(TPZGeoMesh *gmesh, int pOrder, int nVugs){

    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(gmesh->Dimension());
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetAllCreateFunctionsHDiv();

    std::set<int> setMatID = {EMatId, EbcInletId, EbcOutletId, EbcNoFlux, EVugBcId};

    // Add materials (weak formulation)
    TPZNullMaterial<STATE> *mat = new TPZNullMaterial(EMatId, gmesh->Dimension());  
    cmesh->InsertMaterialObject(mat);

    //Create BCs
    int bc_typeN = 1;
    int bc_typeD = 0;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);
  
    TPZBndCond * face2 = mat->CreateBC(mat,EbcNoFlux,bc_typeN,val1,val2);
    cmesh->InsertMaterialObject(face2);
    
    val2[0]=100; 
    TPZBndCond * face = mat->CreateBC(mat,EbcInletId,bc_typeD,val1,val2);
    cmesh->InsertMaterialObject(face);
    
    val2[0]=10; 
    TPZBndCond * face1 = mat->CreateBC(mat,EbcOutletId,bc_typeD,val1,val2);
    cmesh->InsertMaterialObject(face1);

    //TODO Create comp elements of Vug Boundary
    for(int iel = 0; iel < nVugs; iel++) {
        int matid = EVugBcId + iel;
        TPZBndCond *faceVug = mat->CreateBC(mat,matid,bc_typeD,val1,val2);
        cmesh->InsertMaterialObject(faceVug);
    }

    // Set up the computational mesh
    cmesh->AutoBuild(setMatID);

    return cmesh;
}

TPZCompMesh *CreateCompMeshPressure(TPZGeoMesh *gmesh, int pOrder){
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(gmesh->Dimension());
    cmesh->SetDefaultOrder(pOrder);

    std::set<int> setMatID;
    
    if (pOrder < 1) {
    cmesh->SetAllCreateFunctionsDiscontinuous();
    } else {
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }

    // Add materials in atomic mesh
    TPZNullMaterial<STATE> *mat = new TPZNullMaterial(EMatId, gmesh->Dimension());
    setMatID.insert(EMatId);
    cmesh->InsertMaterialObject(mat);
    
    // Set up the computational mesh
    cmesh->AutoBuild(setMatID);
    gmesh->ResetReference();
    setMatID.clear();

    // // Add materials in atomic mesh
    // TPZNullMaterial<STATE> *mat = new TPZNullMaterial(EVugId, gmesh->Dimension());
    // setMatID.insert(EVugId);
    // cmesh->InsertMaterialObject(mat);
    
    // cmesh->AutoBuild(setMatID);
    // gmesh->ResetReference();

    // int ncon = cmesh->NConnects();
    // const int lagLevel = 1; // Lagrange multiplier level
    // for(int i=0; i<ncon; i++)
    // {
    //     TPZConnect &newnod = cmesh->ConnectVec()[i]; 
    //     newnod.SetLagrangeMultiplier(lagLevel);
    // }

    return cmesh;
}

void CreateCompMeshMP(TPZGeoMesh *gmesh, TPZManVector<TPZCompMesh *,2> &cmeshes){

    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDimModel(gmesh->Dimension());
    // cmesh->SetDefaultOrder(1);
    cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;

    TPZMixedDarcyFlow *mat = new TPZMixedDarcyFlow(EMatId, gmesh->Dimension());  
    cmesh->InsertMaterialObject(mat);


void CreateCompMeshMP(TPZGeoMesh *gmesh, TPZManVector<TPZCompMesh *,2> &cmeshes){

    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDimModel(gmesh->Dimension());
    // cmesh->SetDefaultOrder(1);
    cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;

    TPZMixedDarcyFlow *mat = new TPZMixedDarcyFlow(EMatId, gmesh->Dimension());  
    cmesh->InsertMaterialObject(mat);


}
