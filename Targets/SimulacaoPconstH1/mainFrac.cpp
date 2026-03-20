#include "sources.h"

void H1Vugs();

//---------------------------MAIN-----------------------------------
int main (){

    H1Vugs();

    return 0;
}
//-------------------------------------------------------------------

void SetUniqueVugConnect(TPZGeoMesh *gmesh, TPZCompMesh *cmesh){
    
    int nels = gmesh->NElements();
    //TPZVec<int64_t> vugIndex(nVugs, -1);
    TPZVec<int64_t> gelIndex(nels, -1);
    std::map<int, int> matId_connect;

    for (int64_t el = 0; el < cmesh->NElements(); el++){

        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int meshDim = gmesh->Dimension();

        if (gel->Dimension() != meshDim-1) continue; // only frac elements
        if (gel->MaterialId() < EFracId) continue; 

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
                if (gelneigh->MaterialId() == EFracBcId) continue; 
                TPZCompEl *celneigh = gelneigh->Reference();
                celneigh->SetConnectIndex(neighside.Side(), connIndex);
            }
        }
    }
}

void H1Vugs(){
   
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
    dim_name_and_physical_tagCoarse[2]["k11"] = 1;
    dim_name_and_physical_tagCoarse[1]["inlet"] = 2;
    dim_name_and_physical_tagCoarse[1]["outlet"] = 3;
    dim_name_and_physical_tagCoarse[1]["noflux"] = 4;
    dim_name_and_physical_tagCoarse[1]["SmallFract"] = 5;


    std::string filename="/home/marina/programming/Stokes-Darcy-Research/VUGS/Meshes/SingleFracture.msh";
    //std::string filename = "/home/marina/programming/Stokes-Darcy-Research/VUGS/Meshes/FewFractures.msh"

    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse, "SingleFractureMesh");
    
    int ncreated = 0;
    int nels = gmesh->NElements();
    int nElVugBound = 0;
    TPZVec<int64_t> els_cont1d(nels,-1); // vector to store the mat ids of the 1d elements (contours), the rest will be -1

    MeshWithSegmentFrac(gmesh);
    
    PrintGeoMesh(gmesh);

    //Create CompMesh
    TPZCompMesh *cmesh =  new TPZCompMesh(gmesh);
    cmesh->SetName("H1-Frac-cmesh");
  
    //Create Materials
    int dim2d = 2;
    int dim1d=1;

    TPZDarcyFlow *matDarcy = new TPZDarcyFlow(EMatId, dim2d);
    matDarcy->SetConstantPermeability(0.01);
    cmesh->InsertMaterialObject(matDarcy);

    TPZDarcyFlow *matDarcySmallVug = nullptr;
    for(auto Id: fracIds) {
        matDarcySmallVug = new TPZDarcyFlow(Id,dim1d);
        // matDarcySmallVug->SetConstantPermeability(1.0e6);
        cmesh->InsertMaterialObject(matDarcySmallVug);
    }

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

    // //TODO Create comp elements of Frac Boundary
    // val2[0] = 0;
    // for(auto bcId: fracBcIds) {
    //     TPZBndCond *faceFrac = matDarcySmallVug->CreateBC(matDarcySmallVug,bcId,bc_typeN,val1,val2);
    //     cmesh->InsertMaterialObject(faceFrac);
    // }

    cmesh->AutoBuild();

    std::cout << "-----------------Assign unique connect to all vug elements-------------------\n";
    // A "connect" represents a degree of freedom (DOF) or interpolation point where solution values are computed. 
    // It contibutes in one place in the stiffness matrix

    SetUniqueVugConnect(gmesh, cmesh);

    //cmesh->ComputeNodElCon();  // Reconstruye conectividad
    cmesh->CleanUpUnconnectedNodes();

    //Esto hace que el espacio de aproxiación sea H1
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
      
    //Inicializa el tamaño del vector solución
    cmesh->ExpandSolution();
    
    PrintCompMesh(cmesh);
    
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
    vecnames.Push("GradU");
    scalnames.Push("Pressure");
      
    //Configuración del posprocesamiento
    int ref = 0; 
    std::string file_reservoir("Darcy_H1_Frac.vtk");
    std::string file_shape("Shape.vtk");
    //TPZVec<int64_t> equationindices (1,1); // indices of the equations to be postprocessed
    //Analisys->ShowShape(file_shape, equationindices); 
    Analisys->DefineGraphMesh(dim2d,scalnames,vecnames,file_reservoir);
    
    Analisys->PostProcess(ref, dim2d);
}
