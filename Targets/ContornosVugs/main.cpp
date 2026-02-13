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
    std::ofstream file3("TestGeoMesh2Dskel.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file3);
    std::cout<<"Vec_cont1D[1] : "<<els_cont1d[0]<<std::endl;
    std::cout<<"Vec_cont1D[1] : "<<els_cont1d[1]<<std::endl;
    for(int in=0;in<nels2;in++){
        TPZGeoEl *gEl=gmesh->Element(in);
        if(gEl->Dimension()!=1){
            continue;
        }
        if(gEl->MaterialId()<99){
            continue;
        }
        if(gEl->MaterialId()<99==els_cont1d[in]);
        std::cout<<"Correct Index : "<<in<<" Materialid = Value in Vec = "<<els_cont1d[in]<<std::endl;
    }
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
