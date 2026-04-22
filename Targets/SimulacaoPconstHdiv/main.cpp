#include "sources.h"
#include "pzlog.h"

void Hdiv_MixedCT(){
 
    //ReadJson inputData("/home/marina/programming/Stokes-Darcy-Research/VUGS/Inputs/FewVugs.json");
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/FewVugs_V.json");//WORKS OK
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/SingleVug3D.json");// Breaking on side Sideoriented(flux_cmesh), after comment it. it works well
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/FewVugsFracts.json");// Testing both vugs&fract in a same mesh, for the first test it breaks on the TPZAnalyisis(Maybe I should make corrections at the .json file)
    //
    //
    //Hdiv test
    //
    //VUG - High Permebility
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/SingleVug3D_HdivHP.json");//RUN
    //VUG - Const Pressure
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/SingleVug3D_Hdiv.json");//RUN

    //FRATURA
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/FewVugsFracts.json");//
    //
    //
    //H1 test
    //
    //VUG - Const Pressure
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/SingleVug3D_H1.json");//RUN
    //VUG - High Permeability
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/SingleVug3D_H1HP.json");//RUN

    //
    //FRATURA
    ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/Fractures.json");//
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/SingleFracture.json");//
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/SingleFrature3D_H1.json");//


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

            // if (inputData.problemType() == 0){ // Frac boundary elements
            //     CondenseEndFrac(cmesh_mult);
            // }

        PrintCompMesh(Flux_cmesh);
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
//        std::set<int> matToProc;
//        matToProc.insert(700);
//        std::string file_reservoir2("Fract3D_Hdiv.vtk");
//
//        //Analisys->DefineGraphMesh(3,matToProc,scalnames, file_reservoir2, vtkRes);
//        constexpr int vtkRes{0};
//        //Definición de variables escalares y vectoriales a posprocesar
//        TPZStack<std::string,10> scalnames, vecnames;
//        vecnames.Push("Flux");
//        scalnames.Push("Pressure");
//        auto vtk2 = TPZVTKGenerator(cmesh_mult, matToProc,scalnames,file_reservoir2, vtkRes);
//        vtk2.Do();
        PrintCompMesh(cmesh_mult);
        PrintCompMesh(Flux_cmesh);
        PrintCompMesh(Pressure_cmesh);


        // --- Clean up ---
        
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
