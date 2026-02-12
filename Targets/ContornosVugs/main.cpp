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

//    int ncreated = 0;
//    int nels = gmesh->NElements();
//
//
//    for (int iel = 0; iel< nels; iel++) {
//        TPZGeoEl *gel = gmesh->Element(iel);
//        if (!gel){
//            continue;
//        }
//        if (gel->Dimension() != 2) {
//            continue;
//        }
//        int nsides= gel->NSides();
//        int ncorners= gel->NCornerNodes();
//        int firstside= nsides-ncorners-1;
//
//
//        for (int iside = firstside; iside<nsides; iside++) {
//            TPZGeoElSide gelside(gel, iside);
//            int matid = gelside.Element()->MaterialId();
//            TPZStack<TPZGeoElSide> allneigh;
//            gelside.AllNeighbours(allneigh);
////            std::cout<<allneigh[0].Element()<<std::endl;
//            int nneighs = allneigh.size();
//            //verify Dimension
//            int verify =0;
//
//            for (int ineigh=0; ineigh<nneighs; ineigh++) {
//                TPZGeoEl *gelneigh = allneigh[ineigh].Element();
//                int dimen = gelneigh->Dimension();
//
//                if (dimen == 1) {
//                    verify = 1;
//                }
//            }
//            if(verify == 1){
//                continue;
//            }
//
//            for (int ineigh=0; ineigh<nneighs; ineigh++) {
//                TPZGeoEl *gelneigh = allneigh[ineigh].Element();
//                int matNeigh = gelneigh->MaterialId();
//                if (matNeigh != matid && (gel->Dimension() == gelneigh->Dimension()) ) {
//                   gelside.Element()->CreateBCGeoEl(iside, 100);
//                   ncreated++;
//                }
//            }
//        }
//    }
//
//    std::cout<< "se crearon: " << ncreated << " elements"<<std::endl;
//
//    //gmesh->BuildConnectivity();
//    int nels2 = gmesh->NElements();
//    TPZVec<int> verificador(nels2, 0);
//    // creador de contornos por ids
//    int mat=100;
//    for (int iel =nels-1; iel<nels2; iel++) {
//
//        TPZGeoEl * gel = gmesh->Element(iel);
//        if (gel->MaterialId() ==100) {
//            std::cout<<"ok "<<std::endl;
//        }
//        if (!gel) {
//            continue;
//        }
//        if (gel->Dimension() != 1) {
//            continue;
//        }
//        if (verificador[iel]==1) {
//            continue;
//        }
//        if (gel->MaterialId() != 100) {
//            continue;
//        }
//        int side = 1;
//
//        TPZGeoElSide gelside(gel, side);
//
//        TPZStack<TPZGeoElSide> allneigh;
//        gelside.AllNeighbours(allneigh);
//        TPZStack<TPZGeoElSide> allneighdim;
//        findElDim(allneigh, 1, allneighdim);
//        int ntest = allneighdim.size();
//        TPZGeoElSide gelneigh = allneighdim[0];
//        gel->SetMaterialId(mat);
//        while (gel != gelneigh.Element()) {
//            TPZStack<TPZGeoElSide> allneigh;
//            int sidetest = gelneigh.Side();
//            if (sidetest==0) {
//                gelneigh.SetSide(1);
//            }
//            else{
//                gelneigh.SetSide(0);
//            }
//            gelneigh.AllNeighbours(allneigh);
//            TPZStack<TPZGeoElSide> allneighdim;
//            findElDim(allneigh, 1, allneighdim);
//            int indexneig = gelneigh.Element()->Index();
//            verificador[indexneig] =1;
//            gelneigh.Element()->SetMaterialId(mat);
//            gelneigh =allneighdim[0];
//        }
//
//
//        mat++;
//        int ok=0;
//    }
    //std::ofstream file3("TestGeoMesh2Dskel.vtk");
    //TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file3);

    
class VugManager {
private:
    struct VugData {
        int matid_contorno;
        TPZVec<int64_t> elements_contorno;
        TPZVec<int64_t> elements_internos;
        
        
        //VugData() : matid_contorno(-1), //master_connect_idx(-1) {}
    };
    
    std::map<int, VugData> vugs;
    TPZGeoMesh* gmesh;
    TPZCompMesh* cmesh;

public:
    VugManager(TPZGeoMesh* geomesh, TPZCompMesh* compmesh)
        : gmesh(geomesh), cmesh(compmesh) {}
    
    // FASE 1: Crear elementos de contorno
    int CreateBoundaryElements() {
        
        int nels = gmesh->NElements();
        int ncreated = 0;
        
        for (int iel = 0; iel < nels; iel++) {
            TPZGeoEl *gel = gmesh->Element(iel);
            if (!gel || gel->Dimension() != 2) continue;
            
            int nsides = gel->NSides();
            int ncorners = gel->NCornerNodes();
            int firstside = nsides - ncorners - 1;
            
            for (int iside = firstside; iside < nsides; iside++) {
                TPZGeoElSide gelside(gel, iside);
                int matid = gelside.Element()->MaterialId();
                
                TPZStack<TPZGeoElSide> allneigh;
                gelside.AllNeighbours(allneigh);
                
                // Verificar que NO hay vecinos 1D
                // Verificar que NO hay vecinos 1D (igual al código original)
                int verify = 0;
                for (int ineigh = 0; ineigh < allneigh.size(); ineigh++) {
                    TPZGeoEl *gelneigh = allneigh[ineigh].Element();
                    int dimen = gelneigh->Dimension();

                    if (dimen == 1) {
                        verify = 1;
                    }
                }
                if (verify == 1) continue;

                
                // Crear BC si hay vecinos con diferente matid
                for (int ineigh = 0; ineigh < allneigh.size(); ineigh++) {
                    TPZGeoEl *gelneigh = allneigh[ineigh].Element();
                    if (gelneigh->MaterialId() != matid &&
                        gel->Dimension() == gelneigh->Dimension()) {
                        gelside.Element()->CreateBCGeoEl(iside, 100);
                        ncreated++;
                        break;
                    }
                }
            }
        }
        std::cout << "Creados " << ncreated << " elementos de contorno" << std::endl;
        return ncreated;
    }
    
    // FASE 2: Identificar grupos de vugs
    void IdentifyVugs(int nels_before_bc) {

        int nels2 = gmesh->NElements();
        std::vector<int> verificador(nels2, 0);

        int mat = 100;

        for (int iel=nels_before_bc; iel < nels2; iel++) {

            TPZGeoEl *gel = gmesh->Element(iel);
            if (!gel) continue;
            if (gel->Dimension() != 1) continue;
            if (verificador[iel] == 1) continue;
            if (gel->MaterialId() != 100) continue;

            VugData vugdata;
            vugdata.matid_contorno = mat;

            //  guardar primer elemento
            vugdata.elements_contorno.push_back(iel);
            verificador[iel] = 1;

            TPZGeoElSide gelside(gel, 1);

            TPZStack<TPZGeoElSide> allneigh;
            gelside.AllNeighbours(allneigh);

            TPZStack<TPZGeoElSide> allneighdim;
            findElDim(allneigh, 1, allneighdim);

            if (allneighdim.size() == 0) continue;

            TPZGeoElSide gelneigh = allneighdim[0];

            gel->SetMaterialId(mat);

            while (gel != gelneigh.Element()) {

                int indexneig = gelneigh.Element()->Index();

                verificador[indexneig] = 1;
                vugdata.elements_contorno.push_back(indexneig);

                gelneigh.Element()->SetMaterialId(mat);

                int sidetest = gelneigh.Side();
                if (sidetest == 0)
                    gelneigh.SetSide(1);
                else
                    gelneigh.SetSide(0);

                allneigh.clear();
                allneighdim.clear();

                gelneigh.AllNeighbours(allneigh);
                findElDim(allneigh, 1, allneighdim);

                if (allneighdim.size() == 0) break;

                gelneigh = allneighdim[0];
            }

            vugs[mat] = vugdata;
            mat++;
        }
    }

    
    // FASE 3: Unificar connects usando SetConnectIndex
    void CreateInternalElements () {
        int total_internos=0;
        for (auto&[matid,vugdata]:vugs){
            vugdata.elements_internos.Resize(0);
            for(int64_t el_contorno: vugdata.elements_contorno){
                TPZGeoEl* gel_contorno=gmesh->Element(el_contorno);
                
            }
        };
    };
    void PrintBoundaryElements(int matid) {

        auto it = vugs.find(matid);

        if (it == vugs.end()) {
            std::cout << "No existe vug con matid " << matid << std::endl;
            return;
        }

        std::cout << "Vug " << matid << " - Elements contorno:\n";

        for (const auto& el : it->second.elements_contorno) {
            std::cout << el << " ";
        }

        std::cout << std::endl;
    }


    // Métodos de consulta
    void PrintVugSummary() const {
        std::cout << "\n=== RESUMEN DE VUGS ===" << std::endl;
        for (auto& [matid, vug] : vugs) {
            std::cout << "Vug " << matid << ":\n"
            << "  Contornos: " << vug.elements_contorno.size() << "\n"
            << std::endl;
            
     };
  };
};

    VugManager vug_manager(gmesh, cmesh);
    int nels=gmesh->NElements();
    // 1. Crear contornos
    vug_manager.CreateBoundaryElements();
    
    // 2. Identificar vugs
    vug_manager.IdentifyVugs(nels);
    vug_manager.PrintBoundaryElements(110);
 
    
    // 4. Ver resumen
    vug_manager.PrintVugSummary();
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
