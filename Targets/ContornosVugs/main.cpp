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

using namespace std;
//Function to generate a mesh using gmsh library
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void findElDim(TPZStack<TPZGeoElSide> &allneigh, int dim, TPZStack<TPZGeoElSide> &allneighdim);


int main2DFracVug();



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
      //
      //Leitura de Malha com Vugs
      //
      TPZGeoMesh *gmesh = new TPZGeoMesh;
      TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
      dim_name_and_physical_tagCoarse[2]["k11"] = 1;
      dim_name_and_physical_tagCoarse[2]["Vugs"] = 6;
      dim_name_and_physical_tagCoarse[1]["inlet"] = 2;
      dim_name_and_physical_tagCoarse[1]["outlet"] = 3;
      dim_name_and_physical_tagCoarse[1]["noflux"] = 4;

      //
      //Malhas
      //
      std::string filename="/Users/victorvillegassalabarria/python-test/testskelSLICE77SP.msh";
      //std::string filename="/Users/victorvillegassalabarria/Documents/Github/Vugs/FastMesh.msh";
      //Criação de Malhas --- geometrica e computacional
      gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);
      TPZCompMesh *cmesh =  new TPZCompMesh(gmesh);
      void findElDim(TPZStack<TPZGeoElSide> &allneigh, int dim, TPZStack<TPZGeoElSide> &allneighdim);
      TPZVec<int64_t> elements_contorno;
      TPZVec<int64_t> elements_internos;
        int ncreated = 0;
        int nels = gmesh->NElements();

        // CREATION OF BOUNDARY ELEMENTS (MATID=100) IN THE INTERFACE BETWEEN DIFFERENT MATERIALS, Vugs and Porous Matrix
        for (int iel = 0; iel< nels; iel++) {
            TPZGeoEl *gel = gmesh->Element(iel);
            int dimension_el=gel->Dimension();
            int matid_el=gel->MaterialId();
            int nsides=gel->NSides();
            int ncorners=gel->NCornerNodes();
            if(dimension_el!=meshdim)continue;
            if(matid_el!=matid_vug)continue;
            int firstside= nsides-ncorners-1;
            for (int iside = firstside; iside<nsides; iside++) {
                  TPZGeoElSide gelside(gel, iside);
                  int matid = gelside.Element()->MaterialId();
                  auto hasDarcyNeigh=gelside.HasNeighbour(matDarcy);
                  auto hasVugBound=gelside.HasNeighbour(matid_Vugbound);
                  if(hasDarcyNeigh && !hasVugBound) {
                    gelside.Element()->CreateBCGeoEl(iside, 100);
                    ncreated++;
                  }
            }
        }
        std::cout<< "se crearon: " << ncreated << " elements"<<std::endl;

//    //gmesh->BuildConnectivity();
    int nels2 = gmesh->NElements(); // number of elements after creating BCs
    TPZVec<int64_t> els_cont1d(nels2,-1); // vector to store the mat ids of the 1d elements (contours), the rest will be -1 
    std::cout<<nels2<<std::endl;
    int nels1d=0;
    int nels2d=0;
    for (int i=0;i<nels2;i++){
        TPZGeoEl *gel=gmesh->Element(i);
        if (gel->Dimension()== 1) {
//            std::cout<<"Elemento 1d: "<<i<<std::endl;
            nels1d++;
            //break;
        }
        else{
//            std::cout<<"Elemento 2d: "<<i<<std::endl;
        nels2d++;
        }
    }
    std::cout<<"Elementos 1d: "<<nels1d<<std::endl;
    std::cout<<"Elementos 2d: "<<nels2d<<std::endl;

    //elements_contorno.Resize(const int64_t newsize)

    TPZVec<int> verificador(nels2, 0);
    int mat=100;

    // Iterates through newly created 1D boundary elements (those with material ID 100) 
    // and groups connected elements together, assigning each contiguous group a unique material ID
    for (int iel =nels-1; iel<nels2; iel++) { 

        TPZGeoEl * gel = gmesh->Element(iel);
        if (gel->MaterialId() ==100) {
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
        if (gel->MaterialId() != 100) {
            continue;
        }
        int side = 1;

        TPZGeoElSide gelside(gel, side);

        TPZStack<TPZGeoElSide> allneigh;
        gelside.AllNeighbours(allneigh); // Find all neighbors of the current side and store them in allneigh
        TPZStack<TPZGeoElSide> allneighdim;
        findElDim(allneigh, 1, allneighdim); // Only 1D neighbors, which are Boundary elements, are stored in allneighdim 
        int ntest = allneighdim.size();
        TPZGeoElSide gelneigh = allneighdim[0];
        gel->SetMaterialId(mat);
        els_cont1d[iel]=mat; 
        while (gel != gelneigh.Element()) { //loop over connected elements, until the loop returns to the original element, completing the contour group.
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
        mat++;
        int ok=0;
    }

    if(nels2==els_cont1d.size()){
        std::cout<<"Ok"<<std::endl;
    }
    TPZVec<int64_t> AllVugsEls2d; 

    // GROUP ALL 2D ELEMENTS OF VUGS (MATID=6) IN AllVugsEls2d
    for(int ind=0;ind<nels2;ind++){ 
        TPZGeoEl *gEl=gmesh->Element(ind);
        if (!gEl) {
            continue;
        }
        if (gEl->Dimension() != 2) {
            continue;
        }
        if(gEl->MaterialId()!=6){
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
        if (els_cont1d[ind]!=-1) {  // els_cont1d[AllVugsEls2d[ind]
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
                gEl->SetMaterialId(500+matid1); // Assign a new material ID to the 2D vug element based on the material ID of the neighboring 1D boundary element.
                els_cont1d[AllVugsEls2d[ind]] = 500+matid1;  // Fix: índice real
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
                    if(els_cont1d[nidx] > 500) { 
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

    std::cout << "TOTAL ElVugs_identified: " << ElVugs_identified << std::endl;

    

    std::ofstream file3("TestGeoMesh2Dskel.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file3);
//    std::cout<<"Vec_cont1D[1] : "<<els_cont1d[0]<<std::endl;
//    std::cout<<"Vec_cont1D[1] : "<<els_cont1d[1]<<std::endl;

}}


    //devuelve una malla L2

int main (){
    main2DFracVug();

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
