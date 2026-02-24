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
      //
      //Parametros de malha
      //
      int ncreated = 0;
      int nels = gmesh->NElements();
      int matid_vug=6;
      int matDarcy=1;
      int meshdim=gmesh->Dimension();
      int matid_Vugbound=100;
      //
      //Criação de elementos 1D no contorno dos vugs
      //
      for(int iel=0;iel<nels;iel++){
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
    //
    //
    //
    //
    int nels_cont=gmesh->NElements();
    int mat=100;
    TPZVec<int64_t> verificador(nels_cont,-1);
    for(int el=0;el<nels;el++){
        TPZGeoEl *gel=gmesh->Element(el);
        if(gel->MaterialId()!=matid_vug)continue;
        if(verificador[el]!=-1)continue;;
        TPZStack<int> tocheck;
        tocheck.Push(el);
        while(tocheck.size()){
            int64_t elcheck=tocheck.Pop();
            if(verificador[elcheck]!=-1)continue;
            gel = gmesh->Element(elcheck); gel->SetMaterialId(mat+500);
            int nsides=gel->NSides();
            int ncorners=gel->NCornerNodes();
            int firstside= nsides-ncorners-1;
            verificador[elcheck]=mat;
            for (int iside = firstside; iside<nsides; iside++) {
                  TPZGeoElSide gelside(gel, iside);
                  auto hasDarcyNeigh=gelside.HasNeighbour(matDarcy);
                  auto hasVugBound=gelside.HasNeighbour(matid_Vugbound);
                  if(hasDarcyNeigh && !hasVugBound) DebugStop();
                  TPZGeoElSide neighbour=gelside.Neighbour();
                  if(hasVugBound) {
                    TPZGeoEl *boundel = neighbour.Element();
                    boundel->SetMaterialId(mat);
                  }
                  else if(!hasDarcyNeigh){
                      TPZGeoEl *boundel = neighbour.Element();

                      int64_t neighindex=boundel->Index();
                      tocheck.Push(neighindex);
                      }
                  
                
            }
            
            
            
        }
        mat++;
    }


    

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
