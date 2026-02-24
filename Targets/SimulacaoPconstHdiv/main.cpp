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

enum MatID{
    EMatId = 1,
    EVugId = 100,
    EVugBcId = 500,
    EbcInletId = 2,
    EbcOutletId = 3,
    EbcNoFlux = 4,
    ELagrange = 5,
    EInterfaceL = 6,
    EInterfaceR = 7
};

TPZCompMesh *CreateCompMeshFlux(TPZGeoMesh *gmesh);
TPZCompMesh *CreateCompMeshPressure(TPZGeoMesh *gmesh, int pOrder);
void CreateCompMeshMP(TPZMultiphysicsCompMesh *cmesh, TPZManVector<TPZCompMesh *,2> &cmeshes);

void CreateBoundaryElements(TPZGeoMesh *gmesh);
void CreateInterfaceElements(TPZGeoMesh *gmesh);
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

void CreateBoundaryElements(TPZGeoMesh *gmesh, TPZVec<int64_t> &els_cont1d){

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
    int nels2 = gmesh->NElements();
    //TPZVec<int64_t> els_cont1d(nels2,-1); // vector to store the mat ids of the 1d elements (contours), the rest will be -1
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

        if (gel->Dimension() != meshDim-1) continue; // only boundary elements (Dim-1 elems)
        if (gel->MaterialId() < EVugBcId) continue; 
        //int nVug = VugId(gel->MaterialId()); //TODO Melhorar o verificador

        int connIndex = 0;
        auto it = matId_connect.find(gel->MaterialId()); // check if the material ID of the vug boundary already has an associated connect index
        if(it != matId_connect.end()){ // if it has, use the same connect index for all boundaries of the same vug
            connIndex = matId_connect.at(gel->MaterialId());
        }
        else{ // if it doesn't, create a new connect index and associate it with the material ID of the vug boundary
            connIndex = cel->ConnectIndex(0);
            matId_connect.insert({gel->MaterialId(), connIndex});
        }

        //if(vugIndex[nVug] == -1) vugIndex[nVug] = connIndex; // updating vugIndex to tell that the n-th vug has updated its connect to coonIndex
        
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
                if (gelneigh->MaterialId() == EVugId) continue; // ignore neighbors that are part of the vug itself
                TPZCompEl *celneigh = gelneigh->Reference();
                celneigh->SetConnectIndex(neighside.Side(), connIndex);
            }
        }
    }

    //for(int vug = 0; vug < nVugs; vug++){
    //    int connIndex = vugIndex[vug];
    //   if(connIndex == -1) std::cout << "PROBLEM: ConnectID not set for Vug " << vug << std::endl;
    //    else std::cout << "ConnectID: " << connIndex << std::endl;
    //}
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

void CreateInterfaceElements(TPZGeoMesh *gmesh){
    //TODO Create Interface Geo Els

    //TODO Create Lagrange Geo Els

    
    int nel = gmesh->NElements();
    //Loop over elements
    TPZGeoEl *gel = gmesh->Element(iel);
    //Skip elements that are not domain elements

    //Loop over dim-1 sides (edges or faces)
    TPZGeoElSide gelside(gel, side);

    //If
    TPZGeoElBC gelsideWrap(gelside, EWrap);
    TPZGeoElBC gelsideIntR(gelsideWrap, EInterfaceR);
    if(!gelside.HasNeighbour(gBCIds)) {
        TPZGeoElBC gelsideLag(gelsideIntL, ELagrange);
    }   
    TPZGeoElBC gelsideWrap(gelside, EWrap);
    TPZGeoElBC gelsideIntR(gelsideWrap, EInterfaceR);

    //Verify !gelside.HasNeighbour(ELagrange) && !gelside.HasNeighbour(gBCIds)) DebugStop()

}

void Hdiv_MixedCT(){
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
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
