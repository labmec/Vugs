#include "sources.h"
#include "pzlog.h"

void Hdiv_MixedCT(){

    ReadJson inputData("/home/marina/programming/Stokes-Darcy-Research/VUGS/Inputs/FewFractures.json");
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;

    std::string meshName = inputData.MeshName();

    std::string filename = inputData.MeshFile();

    int problemDim = inputData.dim();

    int approxType = inputData.approxType();

    int problemType = inputData.problemType();

    int pressOrder = inputData.pressOrder();


    std::string approxName;

    if(approxType == 0){
        approxName = "_H1";
    }
    else if(approxType == 1){
        approxName = "_Mixed";
    }

    gmesh = generateGMeshWithPhysTagVec(inputData, filename, meshName);

    MeshWithSegment(inputData, gmesh);

    PrintGeoMesh(gmesh);

    std::set<int> elsId, bcId;
    GetAtomicIds(gmesh, elsId, bcId);


    if(approxType){
        
        TPZCompMesh *Flux_cmesh = CreateFluxMesh(gmesh, elsId, bcId, pressOrder, inputData);
        TPZCompMesh *Pressure_cmesh = CreatePressureMesh(gmesh, elsId, bcId, pressOrder, inputData);

        PrintCompMesh(Flux_cmesh);
        PrintCompMesh(Pressure_cmesh);
    
        TPZVec<TPZCompMesh *> meshvec(2);
        meshvec[0]= Flux_cmesh;
        meshvec[1]= Pressure_cmesh;
        TPZMultiphysicsCompMesh *cmesh_mult = CreateMultiMesh(gmesh, meshvec, inputData);

        PrintGeoMesh(gmesh);
        PrintCompMesh(cmesh_mult);

        // const std::string strShape = "Shape.vtk";
        // TPZVec<int64_t> eqIndices(4, 0);
        // Analisys->ShowShape(strShape, eqIndices);

        cmesh_mult->Reference()->ResetReference();
        cmesh_mult->LoadReferences();

        TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh_mult, RenumType::EMetis);

        Solve(Analisys, cmesh_mult, inputData);

        {
            const std::string plotfile = meshName + approxName;
            constexpr int vtkRes{0};
            TPZManVector<std::string, 2> fields = {"Flux", "Pressure", "GradFluxX"};
            auto vtk = TPZVTKGenerator(cmesh_mult, fields, plotfile, vtkRes);
            vtk.Do();
        }
    
        PrintCompMesh(cmesh_mult);
        PrintCompMesh(Flux_cmesh);
        PrintCompMesh(Pressure_cmesh);
        delete cmesh_mult;
    }
    else{

        TPZCompMesh *cmesh = CreateMesh(gmesh, inputData);
        
        PrintCompMesh(cmesh);
    
        //CreateAnalisys
        TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh);
        Analisys->LoadSolution();
        Solve(Analisys, cmesh, inputData);

        //Definición de variables escalares y vectoriales a posprocesar
        TPZStack<std::string,10> scalnames, vecnames;
        vecnames.Push("Flux");
        vecnames.Push("GradU");
        scalnames.Push("Pressure");
        
        //Configuración del posprocesamiento
        int ref = 0; 
        std::string plotfile = meshName + approxName + ".vtk";
  
        Analisys->DefineGraphMesh(problemDim, scalnames, vecnames, plotfile);
        
        Analisys->PostProcess(ref, problemDim);
    }
}


int main (){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    Hdiv_MixedCT();
    return 0;
}
