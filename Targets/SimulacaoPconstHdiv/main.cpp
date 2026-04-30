#include "sources.h"
#include "pzlog.h"

void Hdiv_MixedCT(){

    //ReadJson inputData("/home/marina/programming/Stokes-Darcy-Research/VUGS/Inputs/SeveralFractures.json");
    //ReadJson inputData("/home/marina/programming/Stokes-Darcy-Research/VUGS/Inputs/SeveralVugs.json");
    //ReadJson inputData("/home/marina/programming/Stokes-Darcy-Research/VUGS/Inputs/FewFractures.json");
    ReadJson inputData("/home/marina/programming/Stokes-Darcy-Research/VUGS/Inputs/SingleVug.json");
    //ReadJson inputData("/home/marina/programming/Stokes-Darcy-Research/VUGS/Inputs/SingleFracture.json");
    //ReadJson inputData("/home/marina/programming/Stokes-Darcy-Research/VUGS/Inputs/FractureVug.json");
    
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
    if(problemType == 0){
        approxName += "_0";
    }
    else if(problemType == 1){
        approxName += "_1";
    }

    gmesh = generateGMeshWithPhysTagVec(inputData, filename, meshName);

    PrintGeoMesh(gmesh);

    MeshWithSegment(inputData, gmesh);

    PrintGeoMesh(gmesh);

    std::set<int> elsId, bcId;
    GetAtomicIds(gmesh, elsId, bcId);


    if(approxType == 1){
        
        TPZCompMesh *Flux_cmesh = CreateFluxMesh(gmesh, elsId, bcId, pressOrder, inputData);
        TPZCompMesh *Pressure_cmesh = CreatePressureMesh(gmesh, elsId, bcId, pressOrder, inputData);

        PrintCompMesh(Flux_cmesh);
        PrintCompMesh(Pressure_cmesh);
    
        TPZVec<TPZCompMesh *> meshvec(2);
        meshvec[0]= Flux_cmesh;
        meshvec[1]= Pressure_cmesh;
        TPZMultiphysicsCompMesh *cmesh_mult = CreateMultiMesh(gmesh, meshvec, inputData);

        PrintCompMesh(Flux_cmesh);
        cmesh_mult->ComputeNodElCon();
        PrintCompMesh(cmesh_mult);

        // const std::string strShape = "Shape.vtk";
        // TPZVec<int64_t> eqIndices(4, 0);
        // Analisys->ShowShape(strShape, eqIndices);

        cmesh_mult->Reference()->ResetReference();
        cmesh_mult->LoadReferences();

        TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh_mult, RenumType::ENone);

        SetAnalysis(Analisys, cmesh_mult, inputData);
        {
            TPZFMatrix<STATE> &rhs = Analisys->Rhs();
            rhs.Print(std::cout);
        }
        Analisys->Assemble();
    
        
        // BCInitialSolution(Analisys, cmesh_mult, bcId, inputData, 1);
        // ApplyEquationFilter(Analisys, cmesh_mult, bcId, inputData, 1);
        // NewtonMethod(cmesh_mult, 5, 1.0e-6, 1.0e-6, Analisys);
        TPZFMatrix<STATE> &rhs = Analisys->Rhs();
        Analisys->PrintVectorByElement(std::cout, rhs, 1.e-6);

        {
            const std::string plotfile = meshName + approxName;
            constexpr int vtkRes{0};
            TPZManVector<std::string, 3> fields = {"Flux", "Pressure", "GradFluxX"};
            auto vtk = TPZVTKGenerator(cmesh_mult, fields, plotfile, vtkRes);
            vtk.Do();
        }
    
        PrintCompMesh(cmesh_mult);
        PrintCompMesh(Flux_cmesh);
        PrintCompMesh(Pressure_cmesh);
        delete cmesh_mult;
    }
    else if(approxType == 2){
        
        TPZCompMesh *cmesh = CreateMesh(gmesh, inputData);
        TPZCompMesh *Flux_cmesh = CreateFluxMesh(gmesh, elsId, bcId, pressOrder, inputData);
        TPZCompMesh *Pressure_cmesh = CreatePressureMesh(gmesh, elsId, bcId, pressOrder, inputData);
    
        TPZVec<TPZCompMesh *> meshvec(2);
        meshvec[0]= Flux_cmesh;
        meshvec[1]= Pressure_cmesh;
        TPZMultiphysicsCompMesh *cmesh_mult = CreateMultiMesh(gmesh, meshvec, inputData);

        cmesh_mult->ComputeNodElCon();

        // const std::string strShape = "Shape.vtk";
        // TPZVec<int64_t> eqIndices(4, 0);
        // Analisys->ShowShape(strShape, eqIndices);

        cmesh_mult->Reference()->ResetReference();
        cmesh_mult->LoadReferences();

        TPZLinearAnalysis *AnH1 = new TPZLinearAnalysis(cmesh);

        TPZLinearAnalysis *AnHdiv = new TPZLinearAnalysis(cmesh_mult, RenumType::ENone);

        SetAnalysis(AnH1, cmesh, inputData);
        SetAnalysis(AnHdiv, cmesh_mult, inputData);
        // {
        //     TPZFMatrix<STATE> &rhs = Analisys->Rhs();
        //     rhs.Print(std::cout);
        // }

        {
            const std::string plotfile = meshName + "_H1";
            constexpr int vtkRes{0};
            TPZManVector<std::string, 3> fields = {"Flux", "Pressure", "GradU"};
            auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
            vtk.Do();
        }

        {
            const std::string plotfile = meshName + "_Mixed";
            constexpr int vtkRes{0};
            TPZManVector<std::string, 3> fields = {"Flux", "Pressure", "GradFluxX"};
            auto vtk = TPZVTKGenerator(cmesh_mult, fields, plotfile, vtkRes);
            vtk.Do();
        }

        PrintCompMesh(cmesh);
        PrintCompMesh(cmesh_mult);

        // Estimating Error
        int64_t nels = gmesh->NElements();
        int64_t nelsH1 = cmesh->NElements();
        TPZVec<TPZCompEl*> celH1(nels, nullptr); 
        TPZVec<TPZCompEl*> celHdiv(nels, nullptr);
        TPZSolutionMatrix &solMat = cmesh->ElementSolution();
        solMat.Redim(nelsH1,1);
        GetCompEls(gmesh, cmesh, Flux_cmesh, celH1, celHdiv);
        REAL error = ComputeErrorH1Hdiv(celH1, celHdiv, EMatId, solMat);
        std::cout << "Error: " << error << "\n";

        const char mychar = 'm';
        const char* name = &mychar;
        cmesh->ElementSolution().Print(name);
        //std::cout << CalcElementError(celH1[57], celHdiv[57], inputData) << "\n";

        {
            const std::string plotfile = meshName + "_H1";
            constexpr int vtkRes{0};
            TPZManVector<std::string, 3> fields = {"Flux", "Pressure", "GradU", "EstimatedError"};
            auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
            vtk.Do();
        }

        
        PrintCompMesh(cmesh);
        PrintCompMesh(cmesh_mult);
        PrintCompMesh(Flux_cmesh);
        PrintCompMesh(Pressure_cmesh);
        delete cmesh_mult;
    }
    else{

        TPZCompMesh *cmesh = CreateMesh(gmesh, inputData);
    
        TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh);

        SetAnalysis(Analisys, cmesh, inputData);

        {
            TPZFMatrix<STATE> &rhs = Analisys->Rhs();
            rhs.Print(std::cout);
        }
        Analisys->Assemble();
        
        
        // BCInitialSolution(Analisys, cmesh, bcId, inputData, 0);
        // ApplyEquationFilter(Analisys, cmesh, bcId, inputData, 0);
        // NewtonMethod(cmesh, 5, 1.0e-6, 1.0e-6, Analisys);

        TPZFMatrix<STATE> &rhs = Analisys->Rhs();
        Analisys->PrintVectorByElement(std::cout, rhs, 1.e-6);

        PrintCompMesh(cmesh);
        {
            const std::string plotfile = meshName + approxName;
            constexpr int vtkRes{0};
            TPZManVector<std::string, 3> fields = {"Flux", "Pressure", "GradU"};
            auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
            vtk.Do();
        }
    }
}


int main (){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    Hdiv_MixedCT();
    return 0;
}
