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
void SetUniqueVugConnect(TPZGeoMesh *gmesh, TPZCompMesh *cmesh);
void MeshWithSegmentVugs(TPZGeoMesh *gmesh);

//using namespace cv;
//using namespace std;
//Function to generate a mesh using gmsh library
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void findElDim(TPZStack<TPZGeoElSide> &allneigh, int dim, TPZStack<TPZGeoElSide> &allneighdim);

int main3D();
int main2D();
int Hdiv_MixedCT();
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
int Hdiv_MixedCT(){
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZGeoMesh *gmeshp = new TPZGeoMesh;

    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
    dim_name_and_physical_tagCoarse[2]["k11"] = 1;
    dim_name_and_physical_tagCoarse[1]["inlet"] = 2;
    dim_name_and_physical_tagCoarse[1]["outlet"] = 3;
    dim_name_and_physical_tagCoarse[1]["noflux"] = 4;
    dim_name_and_physical_tagCoarse[2]["Vugs"] = 6;
    
    //std::string filename="/Users/victorvillegassalabarria/python-test/testskelSLICE77SP.msh";
    //std::string filename="/Users/victorvillegassalabarria/Documents/Github/Vugs/FastMesh.msh";
    std::string filename="/Users/victorvillegassalabarria/Documents/Github/Vugs/FewVugsMesh.msh";

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);
    gmeshp = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);

    
    void findElDim(TPZStack<TPZGeoElSide> &allneigh, int dim, TPZStack<TPZGeoElSide> &allneighdim);

    MeshWithSegmentVugs(gmesh);

 
    std::ofstream file20("TestGeoMesh2D.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file20);
    
    
    int order =1;
    TPZHDivApproxCreator hdivCreator(gmesh);
    hdivCreator.ProbType() = ProblemType::EDarcy;
    hdivCreator.SetDefaultOrder(order);
    hdivCreator.SetShouldCondense(false);
    
    
    // Add materials (weak formulation)
    TPZMixedDarcyFlow *matDarcy = new TPZMixedDarcyFlow(1,2);
    //TPZMixedDarcyFlow *matDarcyVugs= new TPZMixedDarcyFlow(100,2);

    matDarcy->SetConstantPermeability(0.01);
    hdivCreator.InsertMaterialObject(matDarcy);
    
    //matDarcyVugs->SetConstantPermeability(1e6);

//    for(int p=100;p<190; p++){
//        TPZMixedDarcyFlow *NewMat= new TPZMixedDarcyFlow(p,2);
//        hdivCreator.InsertMaterialObject(NewMat);
//    }
    
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
    TPZCompMesh *cmeshp =  new TPZCompMesh(gmeshp);
    //TPZCompMesh *cmesh1= new TPZCompMesh(gmesh);
    TPZManVector<TPZCompMesh *, 7> meshvec(2);
    //meshvec[0]= cmesh1;
    meshvec[1]= cmeshp;
    hdivCreator.CreateAtomicMeshes(meshvec, lagmultilevel); // This method increments the lagmultilevel
    TPZMultiphysicsCompMesh *cmesh = nullptr;
    hdivCreator.CreateMultiPhysicsMesh(meshvec, lagmultilevel, cmesh);

    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();

//    for (int iel = 0; iel < cmesh->NElements(); iel++) {
//        auto cel = cmesh->Element(iel);
//        if (!cel) continue;
//
//        if (cel->Material()->Id() == 3) {
//            std::cout << "BC element with ndof = "
//                      << cel->NConnects() << std::endl;
//        }
//    }

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
    return 0;
}



    //devuelve una malla L2

int main (){
    //main2DFracVug();
    //mainDarcy3D();
    //mainMixed();
    Hdiv_MixedCT();
    return 0;
}

void MeshWithSegmentVugs(TPZGeoMesh *gmesh ){
    int ncreated = 0;
    int nels = gmesh->NElements();
    int matid_vug=6;
    int matDarcytag=1;
    int meshdim=gmesh->Dimension();
    int matid_Vugbound=100;
  //
  //Os elementos 2D do primeiro vug vão ter id: mat+500=600
  //O contorno do primeiro vug será 100
    int mat = 100;
  //Vetor verificador onde sera salvada a informação
  //dos elemtentos de cada vug, é inicializado com -1.
    TPZVec<int64_t> verificador(nels,-1);
  //Percorrer todos os elementos da malha
    for(int el = 0; el < nels; el++){
      //Elemento geometrico gel
      TPZGeoEl *gel = gmesh->Element(el);
      if(!gel) continue;
      //Soamente considerar elementos com "matid_vug"
      if(gel->MaterialId() != matid_vug) continue;
      //Não considerar elementos que já foram analizados
      //valor deve ser igual a -1
      if(verificador[el] != -1) continue;
      //Pilha to check que sera utilizada para salvar elementos
      //De um mesmo vug.
      TPZStack<int64_t> tocheck;
      tocheck.Push(el);
      //Enquanto existirem elementos para verificar
      while(tocheck.size()){
          //Retira o último elemento inserido na pilha
          int64_t elcheck = tocheck.Pop();
          //Se já foi visitado, não precisamos processar novamente
          if(verificador[elcheck] != -1) continue;
          //Obtemos novamente o elemento geométrico correspondente
          TPZGeoEl *gelcheck = gmesh->Element(elcheck);
          verificador[elcheck] = mat;
          //Alteramos o MaterialId do elemento bidimensional
          //para mat+500, separando visualmente e numericamente
          //os elementos 2D dos contornos
          gelcheck->SetMaterialId(mat+500);

          int nsides   = gelcheck->NSides();
          int ncorners = gelcheck->NCornerNodes();
          int firstside = nsides - ncorners - 1;
          //Percorrer todas as faces relevantes do elemento
          for(int iside = firstside; iside < nsides; iside++){
              TPZGeoElSide gelside(gelcheck, iside);
              //Verificar a existencia de vizinhos do tipo Darcy
              bool hasDarcyNeigh = gelside.HasNeighbour(matDarcytag);
              //Verificar a existencia de um contorno criado
              bool hasVugBound   = gelside.HasNeighbour(matid_Vugbound);
              //  Crear BC aquí directamente
              if(hasDarcyNeigh && !hasVugBound){
                  gelside.Element()->CreateBCGeoEl(iside, mat);
                  hasVugBound = true;
              }
              //Obter vizinho da face atual
              TPZGeoElSide neighbour = gelside.Neighbour();
              TPZGeoEl *neighgel = neighbour.Element();
              //Se existe elemento de contorno na face, é
              //dado o mesmo material para esse elemento 2D.
              if(hasVugBound){
                  neighgel->SetMaterialId(mat);
              }
              //Os vizinhos internos, são agregados para a pila
              else if(!hasDarcyNeigh){
                  int64_t neighindex = neighgel->Index();
                  if(verificador[neighindex] == -1)
                      tocheck.Push(neighindex);
              }
          }
      }
      // incrementamos o valor de mat, para o seguinte vug.
      mat++;
    }
    
    
}




void SetUniqueVugConnect(TPZGeoMesh *gmesh, TPZCompMesh *cmesh){
    int EVugBcId=500;
    int EVugId=100;
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
