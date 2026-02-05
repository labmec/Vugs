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
int mainDarcy2d();
int mainDarcy3D();
TPZCompMesh* HdivMesh(TPZGeoMesh *);
TPZCompMesh* Pressuremesh(TPZGeoMesh *, int order);
void GetAtomicIds(TPZGeoMesh *geomesh, std::set<int> &volId, std::set<int> &bcId);
void insertAtomicMaterials(TPZCompMesh *cmesh, std::set<int> matIdsVol, std::set<int> matIdsBcs);

//using namespace cv;
//using namespace std;
//Function to generate a mesh using gmsh library
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void findElDim(TPZStack<TPZGeoElSide> &allneigh, int dim, TPZStack<TPZGeoElSide> &allneighdim);


int main3D();
int main2D();
int main2DFracVug();
int mainDarcy3D ();

//int main(){
//
//    return mainDarcy3D();
//}


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
     
      std::ofstream file3("TestGeoMesh2Dskel.vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file3);
      //Create CompMesh
      TPZCompMesh *cmesh =  new TPZCompMesh(gmesh);
  
      //Create Materials
      int matId=1;
      int dim2d = 2;
      int matIdsmallFract=2;
      int matIBigFract=2;
      int dim1d=1;
  
      TPZDarcyFlow *matDarcy = new TPZDarcyFlow(matId, dim2d);
      TPZDarcyFlow *matDarcySmallVug= new TPZDarcyFlow(6,dim2d);

  
      matDarcy->SetConstantPermeability(0.01);
      matDarcySmallVug->SetConstantPermeability(1e9);

  
      cmesh->InsertMaterialObject(matDarcy);
      int bc_id=2;
      int bc_typeN = 1;
      int bc_typeD = 0;
      TPZFMatrix<STATE> val1(1,1,0.0);
      TPZVec<STATE> val2(1,0.0);
      
      int bcinletId = 2;
      int bcOutletId = 3;
      int bcNoFlux = 4;
  
      TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,bcNoFlux,bc_typeN,val1,val2);
      cmesh->InsertMaterialObject(face2);
      
      val2[0]=100; // Valor a ser impuesto como presión en la entrada
      TPZBndCond * face = matDarcy->CreateBC(matDarcy,bcinletId,bc_typeD,val1,val2);
      cmesh->InsertMaterialObject(face);
      
      val2[0]=10; // Valor a ser impuesto como presión en la salida
      TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,bcOutletId,bc_typeD,val1,val2);
      cmesh->InsertMaterialObject(face1);

//      cmesh->InsertMaterialObject(matDarcySmallVug);
//      for (int p=100; p<190; p++) {
//          val2[0]=50;
//          TPZBndCond *contorno=matDarcy->CreateBC(matDarcy,p,bc_typeD,val1,val2);
//          cmesh->InsertMaterialObject(contorno);
//      }
      cmesh->AutoBuild();
      //Esto hace que el espacio de aproxiación sea H1
      cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
      
      //Inicializa el tamaño del vector solución
      cmesh->ExpandSolution();
      
    
      //CreateAnalisys
      TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh);
      bool mustOptimizeBandwidth = false;
     
      //Carga la solución a la malla computacional
      Analisys->LoadSolution();
      
     // Selecciona el método numérico para resolver el problema algebraico
      TPZStepSolver<STATE> step;
      
  //    TPZSSpStructMatrix<STATE> matrix(cmesh);
      step.SetDirect(ELDLt);
    
  //    Analisys->SetStructuralMatrix(matrix);
      

      Analisys->SetSolver(step);
      
      //Ensamblaje de la matriz de rigidez y vector de carga
      Analisys->Assemble();

      //Resolución del sistema algebraico
      Analisys->Solve();

      //Definición de variables escalares y vectoriales a posprocesar
      TPZStack<std::string,10> scalnames, vecnames;
      vecnames.Push("Flux");
      scalnames.Push("Pressure");
      
      //Configuración del posprocesamiento
      int ref =0; // Permite refinar la malla con la solucion obtenida
      std::string file_reservoir("Darcy_H1.vtk");
      Analisys->DefineGraphMesh(dim2d,scalnames,vecnames,file_reservoir);
      //Posprocesamiento
      Analisys->PostProcess(ref, dim2d);
   
      return 0;
}



int main (){
    main2DFracVug();
    return 0;
}
