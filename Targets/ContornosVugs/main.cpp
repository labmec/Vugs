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
      TPZGeoMesh *gmesh = new TPZGeoMesh;
      TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
      dim_name_and_physical_tagCoarse[2]["k11"] = 1;
      dim_name_and_physical_tagCoarse[2]["Vugs"] = 6;
      dim_name_and_physical_tagCoarse[1]["inlet"] = 2;
      dim_name_and_physical_tagCoarse[1]["outlet"] = 3;
      dim_name_and_physical_tagCoarse[1]["noflux"] = 4;



      std::string filename="/Users/victorvillegassalabarria/python-test/testskelSLICE77SP.msh";

      gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);
      TPZCompMesh *cmesh =  new TPZCompMesh(gmesh);
      void findElDim(TPZStack<TPZGeoElSide> &allneigh, int dim, TPZStack<TPZGeoElSide> &allneighdim);
      TPZVec<int64_t> elements_contorno;
      TPZVec<int64_t> elements_internos;
        int ncreated = 0;
        int nels = gmesh->NElements();


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
    //            std::cout<<allneigh[0].Element()<<std::endl;
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
                       gelside.Element()->CreateBCGeoEl(iside, 100);
                       ncreated++;
                    }
                }
            }
        }

        std::cout<< "se crearon: " << ncreated << " elements"<<std::endl;

//    //gmesh->BuildConnectivity();
    int nels2 = gmesh->NElements();
    TPZVec<int64_t> els_cont1d(nels2,-1);
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
    // creador de contornos por ids
    int mat=100;
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
        mat++;
        int ok=0;
    }
//    std::ofstream file3("TestGeoMesh2Dskel.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file3);
//    std::cout<<"Vec_cont1D[1] : "<<els_cont1d[0]<<std::endl;
//    std::cout<<"Vec_cont1D[1] : "<<els_cont1d[1]<<std::endl;
    
    //
    //
    //Comprobacao 1D é consistente:
    //
    //
    ///
//    for(int in=0;in<nels2;in++){
//        TPZGeoEl *gEl=gmesh->Element(in);
//        if(gEl->Dimension()!=1){
//            continue;
//        }
//        if(gEl->MaterialId()<99){
//            continue;
//        }
//        if(gEl->MaterialId()==els_cont1d[in]);
//        std::cout<<"Correct Index : "<<in<<" Materialid = Value in Vec = "<<els_cont1d[in]<<std::endl;
//    }
    
    
    ///
    //////
    //////
    //////
//    if(nels2==els_cont1d.size()){
//        std::cout<<"Ok"<<std::endl;
//    }
//    TPZVec<int64_t> AllVugsEls2d;
//    for(int ind=0;ind<nels2;ind++){
//        TPZGeoEl *gEl=gmesh->Element(ind);
//        if (!gEl) {
//            continue;
//        }
//        if (gEl->Dimension() != 2) {
//            continue;
//        }
//        if(gEl->MaterialId()!=6){
//            continue;
//        }
//        AllVugsEls2d.push_back(ind);
//    }
//    std::cout<<"All vugs 2D size: "<<AllVugsEls2d.size()<<std::endl;
//    int ElVugs_identified=0;
//    for(int ind=0;ind<AllVugsEls2d.size();ind++){
//        TPZGeoEl *gEl=gmesh->Element(AllVugsEls2d[ind]);
//        if (!gEl) {
//            continue;
//        }
//        if (els_cont1d[ind]!=-1) {
//            continue;
//        }
//        //Vizinhança
//        int sides=gEl->NSides();
//        for(int side = 3; side < sides; side++) {
//            TPZGeoElSide gelside(gEl, side);
//            TPZStack<TPZGeoElSide> allneigh;
//            gelside.AllNeighbours(allneigh);
//            TPZStack<TPZGeoElSide> allneighdim;
//            TPZStack<TPZGeoElSide> allneighdim2D;
//            findElDim(allneigh, 1, allneighdim);
//            findElDim(allneigh, 2, allneighdim2D);
//
//            //TPZGeoElSide gelneigh = allneighdim[0];
//            if(allneighdim.NElements()>0){
//                //findElDim(allneigh, 1, allneighdim);
//                TPZGeoElSide gelneigh = allneighdim[0];
//                auto matid1=gelneigh.Element()->MaterialId();
//                std::cout << "  Lado " << side << " → Elemento " << AllVugsEls2d[ind]
//                            << " (matID=" << matid1 << ")" << std::endl;
//                gEl->SetMaterialId(500+matid1);
//                els_cont1d[ind]=500+matid1;
//                ElVugs_identified++;
//            }
//            if (allneighdim2D.NElements()>0) {
//
//
//                TPZGeoElSide gelneigh = allneighdim2D[0];
//                auto matid2=gelneigh.Element()->MaterialId();
//                if(matid2!=6){
//                    gEl->SetMaterialId(500+matid2);
//                    els_cont1d[ind]=500+matid2;
//                    ElVugs_identified++;
//                }
//                else{
//                    continue;
//                }
//
//                continue;
//            }////
//            ///
//        }
//
//        }

    
//
//
    // TERCER BUCLE: Agrupar vug por vug
    if(nels2==els_cont1d.size()){
        std::cout<<"Ok"<<std::endl;
    }
    TPZVec<int64_t> AllVugsEls2d;
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
        AllVugsEls2d.push_back(ind);
    }
    std::cout<<"All vugs 2D size: "<<AllVugsEls2d.size()<<std::endl;
    int ElVugs_identified=0;

    // TU PRIMER BUCLE (sin cambios)
    for(int ind=0;ind<AllVugsEls2d.size();ind++){
        TPZGeoEl *gEl=gmesh->Element(AllVugsEls2d[ind]);
        if (!gEl) {
            continue;
        }
        if (els_cont1d[ind]!=-1) {  // Fix: usa ind correcto
            continue;
        }
        int sides=gEl->NSides();
        //bool marked = false;
        for(int side = 3; side < sides; side++) {  // side=4 (no 3)
            TPZGeoElSide gelside(gEl, side);
            TPZStack<TPZGeoElSide> allneigh;
            gelside.AllNeighbours(allneigh);
            TPZStack<TPZGeoElSide> allneighdim, allneighdim2D;
            findElDim(allneigh, 1, allneighdim);
            findElDim(allneigh, 2, allneighdim2D);

            if(allneighdim.NElements()>0){
                TPZGeoElSide gelneigh = allneighdim[0];
                auto matid1=gelneigh.Element()->MaterialId();
                std::cout << "  Lado " << side << " → Elemento " << AllVugsEls2d[ind]
                          << " (matID=" << matid1 << ")" << std::endl;
                gEl->SetMaterialId(500+matid1);
                els_cont1d[AllVugsEls2d[ind]] = 500+matid1;  // Fix: índice real
                ElVugs_identified++;
                //marked = true;
                break;  // Solo un contorno por elemento
            }
        }
    }
//    int libres = 0;
//    for(int j=0; j<AllVugsEls2d.size(); j++) {
//        if(els_cont1d[AllVugsEls2d[j]] == -1) libres++;
//    }
//    std::cout << "2D LIBRES para flood-fill: " << libres << std::endl;
//
//    // *** FLOOD-FILL SOLO PARA LOS QUE FALTAN ***
//    int i = 0;
//    while(i < AllVugsEls2d.size() && libres > 0) {
//        int64_t seed = AllVugsEls2d[i];
//        if(els_cont1d[seed] != -1) { i++; continue; }
//
//        std::cout << "🔍 Buscando contorno para seed " << seed << std::endl;
//
//        // Busca contorno CON TU MÉTODO
//        TPZGeoEl *seedEl = gmesh->Element(seed);
//        int vugID = -1;
//        int sides = seedEl->NSides();
//        for(int side=3; side<sides; side++) {
//            TPZGeoElSide gelside(seedEl, side);
//            TPZStack<TPZGeoElSide> allneigh;
//            gelside.AllNeighbours(allneigh);
//            TPZStack<TPZGeoElSide> allneighdim;
//            findElDim(allneigh, 1, allneighdim);
////            if(allneighdim.NElements()>0) {
////                vugID = 500 + allneighdim[0].Element()->MaterialId();
////                std::cout << "  ✅ Seed " << seed << " → vug " << vugID-500 << std::endl;
////                break;
////            }
//        }
//
//        if(vugID == -1) {
////            std::cout << "  ❌ Seed " << seed << " sin contorno, saltando" << std::endl;
//            i++;
//            libres--;
//            continue;
//        }
//
//        // Flood-fill MÁS SIMPLE
//        TPZStack<int64_t> stack;
//        stack.Push(seed);
//        els_cont1d[seed] = vugID;
//        seedEl->SetMaterialId(vugID);
//        ElVugs_identified++;
//
//        int expanded = 0;
//        while(stack.NElements() != 0) {
//            int64_t curr = stack.Pop();
//            TPZGeoEl *curEl = gmesh->Element(curr);
//            for(int side=3; side<curEl->NSides(); side++) {
//                TPZGeoElSide gelside(curEl, side);
//                TPZStack<TPZGeoElSide> allneigh;
//                gelside.AllNeighbours(allneigh);
//                TPZStack<TPZGeoElSide> neigh2D;
//                findElDim(allneigh, 2, neigh2D);
//
//                for(int n=0; n<neigh2D.NElements(); n++) {
//                    int64_t nidx = neigh2D[n].Element()->Index();
//                    if(els_cont1d[nidx] == -1 && neigh2D[n].Element()->MaterialId()==6) {
//                        els_cont1d[nidx] = vugID;
//                        neigh2D[n].Element()->SetMaterialId(vugID);
//                        ElVugs_identified++;
//                        stack.Push(nidx);
//                        expanded++;
//                    }
//                }
//            }
//        }
//        std::cout << "✅ Vug " << vugID-500 << " expandido +" << expanded << " elementos" << std::endl;
//        i++;
//    }
//
//    std::cout << "TOTAL ElVugs_identified: " << ElVugs_identified << std::endl;
    // *** PROPAGAR DESDE LOS 1007 MARCADOS ***
    bool expanded = true;
    while(expanded) {
        expanded = false;
        for(int ind=0; ind<AllVugsEls2d.size(); ind++) {
            int64_t elIdx = AllVugsEls2d[ind];
            if(els_cont1d[elIdx] != -1) continue;  // Ya marcado
            
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
                    if(els_cont1d[nidx] > 500) {  // Tiene vecino MARCADO
                        els_cont1d[elIdx] = els_cont1d[nidx];  // Hereda
                        gEl->SetMaterialId(els_cont1d[elIdx]);
                        ElVugs_identified++;
                        expanded = true;
                        break;
                    }
                }
                if(expanded) break;
            }
            if(expanded) break;
        }
    }

    std::cout << "TOTAL ElVugs_identified: " << ElVugs_identified << std::endl;

    

    std::ofstream file3("TestGeoMesh2Dskel.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file3);
//    std::cout<<"Vec_cont1D[1] : "<<els_cont1d[0]<<std::endl;
//    std::cout<<"Vec_cont1D[1] : "<<els_cont1d[1]<<std::endl;

}


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
