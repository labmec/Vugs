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
//TPZCompMesh* HdivMesh(TPZGeoMesh *);
//TPZCompMesh* Pressuremesh(TPZGeoMesh *, int order);
TPZCompMesh *CreateFluxMesh(TPZGeoMesh *, std::set<int> &volId, std::set<int> &bcId);
TPZCompMesh *CreatePressureMesh(TPZGeoMesh *,std::set<int> &volId, std::set<int> &bcId,int order);

void GetAtomicIds(TPZGeoMesh *geomesh, std::set<int> &volId, std::set<int> &bcId);
void insertAtomicMaterials(TPZCompMesh *cmesh, std::set<int> matIdsVol, std::set<int> matIdsBcs,int typeMesh);
void SetUniqueVugConnect(TPZGeoMesh *gmesh, TPZCompMesh *cmesh);
void MeshWithSegmentVugs(TPZGeoMesh *gmesh);
void PrintCompMesh(TPZCompMesh *cmesh);

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
    //gmeshp = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);

    
    //void findElDim(TPZStack<TPZGeoElSide> &allneigh, int dim, TPZStack<TPZGeoElSide> &allneighdim);
    MeshWithSegmentVugs(gmesh);
    std::set<int> volId, bcId;
    GetAtomicIds(gmesh, volId, bcId);
    TPZCompMesh *Flux_cmesh=CreateFluxMesh(gmesh,volId,bcId);
    
    TPZCompMesh *Pressure_cmesh=CreatePressureMesh(gmesh,volId,bcId,0);
    
    TPZMultiphysicsCompMesh *cmesh_mult= new TPZMultiphysicsCompMesh(gmesh);
    TPZVec<TPZCompMesh *> meshvec(2);
    
    meshvec[0]= Flux_cmesh;
    meshvec[1]= Pressure_cmesh;
    
    std::ofstream file20("TestGeoMesh2D.vtk");
    std::ofstream file21("Test_cmeshFlux.vtk");
    std::ofstream file22("Test_cmeshPressure.vtk");
    std::ofstream file23("Test_cmeshPressure.txt");

    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file20);
    TPZVTKGeoMesh::PrintCMeshVTK(Flux_cmesh,file21);
    TPZVTKGeoMesh::PrintCMeshVTK(Pressure_cmesh,file22);
    Pressure_cmesh->Print(file23);
//    int order =1;
//    TPZHDivApproxCreator hdivCreator(gmesh);
//    hdivCreator.ProbType() = ProblemType::EDarcy;
//    hdivCreator.SetDefaultOrder(order);
//    hdivCreator.SetShouldCondense(false);
    
    
    // Add materials (weak formulation)
    TPZMixedDarcyFlow *matDarcy = new TPZMixedDarcyFlow(1,2);
    //TPZMixedDarcyFlow *matDarcyVugs= new TPZMixedDarcyFlow(100,2);
    cmesh_mult->InsertMaterialObject(matDarcy);
    matDarcy->SetConstantPermeability(0.01);
//    hdivCreator.InsertMaterialObject(matDarcy);
    
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
    int dim2d=2;
    int bcinletId = 2;
    int bcOutletId = 3;
    int bcNoFlux = 4;

    TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,bcNoFlux,bc_typeN,val1,val2);
    cmesh_mult->InsertMaterialObject(face2);

//    hdivCreator.InsertMaterialObject(face2);

    
    val2[0]=100; // Valor a ser impuesto como presión en la entrada
    TPZBndCond * face = matDarcy->CreateBC(matDarcy,bcinletId,bc_typeD,val1,val2);
    cmesh_mult->InsertMaterialObject(face);

//    hdivCreator.InsertMaterialObject(face);

    val2[0]=14; // Valor a ser impuesto como presión en la salida
    TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,bcOutletId,bc_typeD,val1,val2);
    cmesh_mult->InsertMaterialObject(face1);

//    hdivCreator.InsertMaterialObject(face1);
    //cmesh_mult->ApproxSpace().SetAllCreateFunctionsContinuous();
    
          //Inicializa el tamaño del vector solución
    cmesh_mult->ExpandSolution();
    cmesh_mult->ApproxSpace().Style()= TPZCreateApproximationSpace::EMultiphysics;
    cmesh_mult->BuildMultiphysicsSpace(meshvec);
        //TPZManVector<int, 2> active_approx_spaces(2, 1);
        //cmesh_mult->BuildMultiphysicsSpace(active_approx_spaces, meshvec);
    cmesh_mult->InitializeBlock();
    std::cout<<cmesh_mult->Element(1)<<std::endl;
        
    bool mustOp = false;

    //CreateAnalisys
    TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh_mult);
        
    //new TPZLinearAnalysis(cmesh_mult);
        

    //TPZAnalysis *Analisys = new TPZAnalysis(cmesh_mult,true);
    //      bool mustOptimizeBandwidth = false;
         
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
          
          //Configuración del posprocesamiento
          int ref =0; // Permite refinar la malla con la solucion obtenida
          std::string file_reservoir("SolVictorCTmesh.vtk");
          Analisys->DefineGraphMesh(dim2d,scalnames,vecnames,file_reservoir);
          //Posprocesamiento
          Analisys->PostProcess(ref, dim2d);
//
//
//
//
//    int lagmultilevel = 1;
//    TPZCompMesh *cmeshp =  new TPZCompMesh(gmeshp);
//    //TPZCompMesh *cmesh1= new TPZCompMesh(gmesh);
//    TPZManVector<TPZCompMesh *, 7> meshvec(2);
//    //meshvec[0]= cmesh1;
////    meshvec[1]= cmeshp;
//    hdivCreator.CreateAtomicMeshes(meshvec, lagmultilevel); // This method increments the lagmultilevel
//    TPZMultiphysicsCompMesh *cmesh = nullptr;
//    hdivCreator.CreateMultiPhysicsMesh(meshvec, lagmultilevel, cmesh);
//
//    cmesh->Reference()->ResetReference();
//    cmesh->LoadReferences();
//
////    for (int iel = 0; iel < cmesh->NElements(); iel++) {
////        auto cel = cmesh->Element(iel);
////        if (!cel) continue;
////
////        if (cel->Material()->Id() == 3) {
////            std::cout << "BC element with ndof = "
////                      << cel->NConnects() << std::endl;
////        }
////    }
//
//    TPZLinearAnalysis anMixed(cmesh, RenumType::EMetis);
//  #ifdef PZ_USING_MKL
//    TPZSSpStructMatrix<STATE> matMixed(cmesh);
//  #else
//    TPZFStructMatrix<STATE> matMixed(cmesh);
//  #endif
//    matMixed.SetNumThreads(0);
//    anMixed.SetStructuralMatrix(matMixed);
//    TPZStepSolver<STATE> stepMixed;
//    stepMixed.SetDirect(ELDLt);
//    anMixed.SetSolver(stepMixed);
//    anMixed.Run();
//
//    // ---- Plotting ---
//
//    {
//      const std::string plotfile = "darcy_mixed";
//      constexpr int vtkRes{0};
//      TPZManVector<std::string, 2> fields = {"Flux", "Pressure"};
//      auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
//      vtk.Do();
//    }
//
//    // --- Clean up ---
//    delete cmesh;
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
    int EVugId = 500;
    int EVugBcId = 100;

    int nels = gmesh->NElements();
    //TPZVec<int64_t> vugIndex(nVugs, -1);
    TPZVec<int64_t> gelIndex(nels, -1);
    std::map<int, int> matId_connect;

    for (int64_t el = 0; el < cmesh->NElements(); el++){

        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int meshDim = gmesh->Dimension();

        //if (gel->Dimension() != meshDim-1) continue; // only boundary elements (Dim-1 elems)
        if (gel->MaterialId() < EVugId) continue;

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

void insertAtomicMaterials(TPZCompMesh *cmesh, std::set<int> matIdsVol, std::set<int> matIdsBcs, int typeMesh){
    int dim = cmesh->Dimension();
    if (typeMesh==0){//Se for malha de fluxo não insertar vugs
        for (auto iD:matIdsVol) {
            if(iD<499){
            TPZNullMaterial <STATE> *matDarcy = new TPZNullMaterial(iD, dim);
            cmesh->InsertMaterialObject(matDarcy);
            }
            else if(iD>499)continue;
            
        }
        for (auto iD:matIdsBcs) {
            
            TPZNullMaterial<STATE> * face2 = new TPZNullMaterial(iD, dim-1);
            cmesh->InsertMaterialObject(face2);
            


        }
    }
    else if (typeMesh==1){//Se for malha de pressão
        for (auto iD:matIdsVol) {
            
            TPZNullMaterial <STATE> *matDarcy = new TPZNullMaterial(iD, dim);
            cmesh->InsertMaterialObject(matDarcy);
          
        }
        for (auto iD:matIdsBcs) {
            if(iD<99){
            TPZNullMaterial<STATE> * face2 = new TPZNullMaterial(iD, dim-1);
            cmesh->InsertMaterialObject(face2);
            }
            else if(iD>99)continue;

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
    int dim2d = 2;
    int typeMesh=0;//Malha de fluxo
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
    int TypeMesh=1;//Malha de pressão
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    //GetAtomicIds(gmesh, volId, bcId);
    insertAtomicMaterials(cmesh, volId, bcId,TypeMesh);
    cmesh->AutoBuild();
    SetUniqueVugConnect(gmesh, cmesh);

    //cmesh->AutoBuild();
        
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
void PrintCompMesh(TPZCompMesh *cmesh)
{
    std::cout << "\nPrinting multiphysics mesh in .txt and .vtk formats...\n";

    std::ofstream VTKCompMeshFile(cmesh->Name() + ".vtk");
    std::ofstream TextCompMeshFile(cmesh->Name() + ".txt");

    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, VTKCompMeshFile);
    cmesh->Print(TextCompMeshFile);
}

