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
int VugId(int matid);

enum EMatIds {
    MatDarcy = 1,
    inlet = 2,
    outlet = 3,
    noflux = 4,
    MatVugs = 6,
    MatBoundVug = 100,
    MatVug = 500
};


int H1Vugs();



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

int H1Vugs(){
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
    dim_name_and_physical_tagCoarse[2]["k11"] = 1;
    dim_name_and_physical_tagCoarse[2]["Vugs"] = 6;
    dim_name_and_physical_tagCoarse[1]["inlet"] = 2;
    dim_name_and_physical_tagCoarse[1]["outlet"] = 3;
    dim_name_and_physical_tagCoarse[1]["noflux"] = 4;


    //std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/testskelSLICE77SP.msh";
    std::string filename="/home/itopo/Stokes-Darcy_Research/Vugs/testskelSLICE77SP.msh";

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);
    int ncreated = 0;
    int nels = gmesh->NElements();
    int nVugs = 0;

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
                    gelside.Element()->CreateBCGeoEl(iside, 100);
                    ncreated++;
                }
            }
        }
    }

    std::cout<< "se crearon: " << ncreated << " elements"<<std::endl;

    gmesh->BuildConnectivity();
    int nels2 = gmesh->NElements();
    TPZVec<int> verificador(nels2, 0);
    int mat = 100; //TODO
    // creador de contornos por ids
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
            gelneigh =allneighdim[0];
        }
        nVugs++; //TODO Based on the fact that each vug has its unique boundary
        mat++;
        int ok=0;
    }
    
    std::ofstream file3("TestGeoMesh2Dskel.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file3);

    //Create CompMesh
    TPZCompMesh *cmesh =  new TPZCompMesh(gmesh);
  
    //Create Materials
    int dim2d = 2;
    int dim1d=1;
    TPZDarcyFlow *matDarcy = new TPZDarcyFlow(MatDarcy, dim2d);
    TPZDarcyFlow *matDarcySmallVug= new TPZDarcyFlow(MatVugs,dim2d);

    matDarcy->SetConstantPermeability(0.01);
    matDarcySmallVug->SetConstantPermeability(1.0e6);
    
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


    //TODO Create comp elements of Vug Boundary
    for(int iel = 0; iel < nVugs; iel++) {
        int matid = MatBoundVug + iel;
        TPZBndCond *faceVug = matDarcySmallVug->CreateBC(matDarcySmallVug,matid,bc_typeD,val1,val2);
        cmesh->InsertMaterialObject(faceVug);
    }

    cmesh->AutoBuild();
    std::cout << "\n===== ASSIGNING FIRST CONNECT FOR ALL ELEMENTS IN EACH GROUP =====\n";
    // A "connect" represents a degree of freedom (DOF) or interpolation point where solution values are computed.
    std::map<int, int64_t> primer_connect_grupo;  // {{matid, first_connect}, ...}

    TPZVec<int64_t> vugIndex(nVugs, -1);

    for (int64_t el = 0; el < cmesh->NElements(); el++){

        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int meshDim = gmesh->Dimension();

        if (gel->Dimension() != meshDim-1) continue; // only boundary elements (Dim-1 elems)
        if (gel->MaterialId() < 100) continue; //TODO

        int nVug = VugId(gel->MaterialId()); //TODO Melhorar o verificador
        int connIndex = cel->ConnectIndex(0);

        if(vugIndex[nVug] == -1) vugIndex[nVug] = connIndex; // updating vugIndex to tell that the n-th vug has updated its connect to coonIndex
        
        int nsides = gel->NSides();
        int nVertex = gel->NCornerNodes();
        
        for(int side = 0; side < nVertex; side++){ // sides associated with vertices
            TPZGeoElSide gelside(gel, side);
            TPZStack<TPZGeoElSide> allneigh;
            gelside.AllNeighbours(allneigh);
            int nneighs = allneigh.size();

            for(int neigh = 0; neigh < nneighs; neigh++){
                TPZGeoElSide neighside = allneigh[neigh];
                TPZGeoEl *gelneigh = neighside.Element();
                if (gelneigh->MaterialId() == MatVugs) continue; //TODO Está pegando tambem o elemento do contorno vizinho, talvez esteja acontecendo retrabalho
                TPZCompEl *celneigh = gelneigh->Reference();
                celneigh->SetConnectIndex(neighside.Side(), connIndex);
            }
        }
    }

    for(int vug = 0; vug < nVugs; vug++){
        int connIndex = vugIndex[vug];
        if(connIndex == -1) std::cout << "PROBLEM: ConnectID not set for Vug " << vug << std::endl;
        else std::cout << "ConnectID: " << connIndex << std::endl;
    }

    //cmesh->ComputeNodElCon();  // Reconstruye conectividad
    cmesh->CleanUpUnconnectedNodes();

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
      
    //TPZSSpStructMatrix<STATE> matrix(cmesh);
    step.SetDirect(ELDLt);
    
    //Analisys->SetStructuralMatrix(matrix);
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
    H1Vugs();
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

int VugId(int matID){
    return matID - 100; //TODO Calcular numero do vug a partir do matId dele (500 + algo)
}