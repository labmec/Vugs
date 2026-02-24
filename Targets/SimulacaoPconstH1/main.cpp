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

enum MatID{
    EMatId = 1,
    EVugId = 100,
    EVugBcId = 500,
    EbcInletId = 2,
    EbcOutletId = 3,
    EbcNoFlux = 4
};

int CreateBoundaryElements(TPZGeoMesh *gmesh, TPZVec<int64_t> &els_cont1d);
void SetMaterialIdVug(TPZGeoMesh *gmesh, TPZVec<int64_t> &els_cont1d); //TODO 
void SetUniqueVugConnect(TPZGeoMesh *gmesh, TPZCompMesh *cmesh);

enum MatID{
    EMatId = 1,
    EVugId = 100,
    EVugBcId = 500,
    EbcInletId = 2,
    EbcOutletId = 3,
    EbcNoFlux = 4
};

int CreateBoundaryElements(TPZGeoMesh *gmesh, TPZVec<int64_t> &els_cont1d);
void SetMaterialIdVug(TPZGeoMesh *gmesh, TPZVec<int64_t> &els_cont1d); //TODO 
void SetUniqueVugConnect(TPZGeoMesh *gmesh, TPZCompMesh *cmesh);

//Function to generate a mesh using gmsh library
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void findElDim(TPZStack<TPZGeoElSide> &allneigh, int dim, TPZStack<TPZGeoElSide> &allneighdim);
void PrintCompMesh(TPZCompMesh *cmesh);

void H1Vugs();

//---------------------------MAIN-----------------------------------
int main (){

    H1Vugs();

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
            //cel->SetConnectIndex(side, connIndex);
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

void H1Vugs(){
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
    dim_name_and_physical_tagCoarse[2]["k11"] = EMatId;
    dim_name_and_physical_tagCoarse[2]["Vugs"] = EVugId;
    dim_name_and_physical_tagCoarse[1]["inlet"] = EbcInletId;
    dim_name_and_physical_tagCoarse[1]["outlet"] = EbcOutletId;
    dim_name_and_physical_tagCoarse[1]["noflux"] = EbcNoFlux;


    //std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/testskelSLICE77SP.msh";
    std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/FewVugsMesh.msh";
    //std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/FastMesh.msh";
    //std::string filename="/home/itopo/Stokes-Darcy_Research/Vugs/testskelSLICE77SP.msh";

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);
    int ncreated = 0;
    int nels = gmesh->NElements();
    int nElVugBound = 0;
    TPZVec<int64_t> els_cont1d(nels,-1); // vector to store the mat ids of the 1d elements (contours), the rest will be -1

    int nVugs = CreateBoundaryElements(gmesh, els_cont1d);
    
    std::ofstream file3("TestGeoMesh2Dskel.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file3);

    //Create CompMesh
    TPZCompMesh *cmesh =  new TPZCompMesh(gmesh);
  
    //Create Materials
    int dim2d = 2;
    int dim1d=1;
    TPZDarcyFlow *matDarcy = new TPZDarcyFlow(EMatId, dim2d);
    TPZDarcyFlow *matDarcySmallVug= new TPZDarcyFlow(EVugId,dim2d);

    matDarcy->SetConstantPermeability(0.01);
    matDarcySmallVug->SetConstantPermeability(1.0e6);
    
    cmesh->InsertMaterialObject(matDarcy);

    //Create BCs
    int bc_id=2;
    int bc_typeN = 1;
    int bc_typeD = 0;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);
  
    TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,EbcNoFlux,bc_typeN,val1,val2);
    cmesh->InsertMaterialObject(face2);
    
    val2[0]=100; // Valor a ser impuesto como presión en la entrada
    TPZBndCond * face = matDarcy->CreateBC(matDarcy,EbcInletId,bc_typeD,val1,val2);
    cmesh->InsertMaterialObject(face);
    
    val2[0]=10; // Valor a ser impuesto como presión en la salida
    TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,EbcOutletId,bc_typeD,val1,val2);
    cmesh->InsertMaterialObject(face1);

    val2[0]=1000;
    //TODO Create comp elements of Vug Boundary
    for(int iel = 0; iel < nVugs; iel++) {
        int matid = EVugBcId + iel;
        TPZBndCond *faceVug = matDarcySmallVug->CreateBC(matDarcySmallVug,matid,bc_typeD,val1,val2);
        cmesh->InsertMaterialObject(faceVug);
    }

    cmesh->AutoBuild();

    std::cout << "-----------------Assign unique connect to all vug elements-------------------\n";
    // A "connect" represents a degree of freedom (DOF) or interpolation point where solution values are computed. 
    // It contibutes in one place in the stiffness matrix

    SetUniqueVugConnect(gmesh, cmesh);

    //cmesh->ComputeNodElCon();  // Reconstruye conectividad
    cmesh->CleanUpUnconnectedNodes();

    //Esto hace que el espacio de aproxiación sea H1
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    //Esto hace que el espacio de aproxiación sea H1
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
      
    //Inicializa el tamaño del vector solución
    cmesh->ExpandSolution();
    
    PrintCompMesh(cmesh);
    
    //CreateAnalisys
    TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh);
    bool mustOptimizeBandwidth = false;
    //CreateAnalisys
    TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh);
    bool mustOptimizeBandwidth = false;
     
    //Carga la solución a la malla computacional
    Analisys->LoadSolution();
    //Carga la solución a la malla computacional
    Analisys->LoadSolution();
      
    // Selecciona el método numérico para resolver el problema algebraico
    TPZStepSolver<STATE> step;
    // Selecciona el método numérico para resolver el problema algebraico
    TPZStepSolver<STATE> step;
      
    //TPZSSpStructMatrix<STATE> matrix(cmesh);
    step.SetDirect(ELDLt);
    //TPZSSpStructMatrix<STATE> matrix(cmesh);
    step.SetDirect(ELDLt);
    
    //Analisys->SetStructuralMatrix(matrix);
    Analisys->SetSolver(step);
    //Analisys->SetStructuralMatrix(matrix);
    Analisys->SetSolver(step);
      
    //Ensamblaje de la matriz de rigidez y vector de carga
    Analisys->Assemble();
    //Ensamblaje de la matriz de rigidez y vector de carga
    Analisys->Assemble();

    //Resolución del sistema algebraico
    Analisys->Solve();
    //Resolución del sistema algebraico
    Analisys->Solve();

    //Definición de variables escalares y vectoriales a posprocesar
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    //Definición de variables escalares y vectoriales a posprocesar
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
      
    //Configuración del posprocesamiento
    int ref = 0; 
    std::string file_reservoir("Darcy_H1.vtk");
    std::string file_shape("Shape.vtk");
    //TPZVec<int64_t> equationindices (1,1); // indices of the equations to be postprocessed
    //Analisys->ShowShape(file_shape, equationindices); //TODO Verify if this is the correct way to show shape functions
    Analisys->DefineGraphMesh(dim2d,scalnames,vecnames,file_reservoir);
    
    Analisys->PostProcess(ref, dim2d);
}

void findElDim(TPZStack<TPZGeoElSide> &allneigh, int dim, TPZStack<TPZGeoElSide> &allneighdim){
    int nels = allneigh.size();
    for (int iel =0; iel<nels; iel++) {
        if (allneigh[iel].Element()->Dimension()==dim) {
            allneighdim.push_back(allneigh[iel]);
        }
    }
}

int VugId(int matID){
    return matID - 100; //TODO Calcular numero do vug a partir do matId dele (500 + algo)
}

void PrintCompMesh(TPZCompMesh *cmesh)
{
    std::cout << "\nPrinting multiphysics mesh in .txt and .vtk formats...\n";

    std::ofstream VTKCompMeshFile(cmesh->Name() + ".vtk");
    std::ofstream TextCompMeshFile(cmesh->Name() + ".txt");

    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, VTKCompMeshFile);
    cmesh->Print(TextCompMeshFile);
}