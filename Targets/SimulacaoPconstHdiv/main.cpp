#include "sources.h"
#include "pzlog.h"

void Hdiv_MixedCT(){
 
    //ReadJson inputData("/home/marina/programming/Stokes-Darcy-Research/VUGS/Inputs/FewVugs.json");
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/NewInputs/FewVugs_V.json");//WORKS OK
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/NewInputs/2Dsem.json");//WORKS OK

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
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/NewInputs/SingleVug3D.json");//RUN <---
    ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/NewInputs/CILAMCE/SingleVug3D_FluxZf6.json");//RUN <---

    //VUG - High Permeability
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/SingleVug3D_H1HP.json");//RUN

    //
    //FRATURA
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/FewFractures.json");//
//ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/SingleFrature3D_Hdiv.json");//
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/2Dsem.json");
    
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/NewInputs/3Dsem.json");

    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/FewVugs_V.json");//

    TPZGeoMesh *gmesh = new TPZGeoMesh;

    std::string meshName = inputData.MeshName();

    std::string filename = inputData.MeshFile();

    int problemDim = inputData.dim();

    int approxType = inputData.approxType();

    int problemType = inputData.problemType();

    int pressOrder = inputData.pressOrder();
    std::map<int, TPZVec<STATE>> IntegralFluxoNormalCmult(const std::set<int>& bcMatId, TPZMultiphysicsCompMesh *cmesh);
    std::map<int, TPZVec<STATE>> IntegralFluxoNormalCmesh(const std::set<int>& bcMatId, TPZCompMesh *cmesh);


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
        std::cout << "Flux mesh solution rows = "
                  << Flux_cmesh->Solution().Rows() << std::endl;

        std::cout << "Mult mesh solution rows = "
                  << cmesh_mult->Solution().Rows() << std::endl;

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
        //std::set<int> bcids={2,3,4};
        std::set<int> bcids={2,3,4,8,9,11};

        auto NormalfluxesCmul=IntegralFluxoNormalCmult(bcids, cmesh_mult);
        //auto NormalfluxesCmesh=IntegralFluxoNormalCmesh(bcids, Flux_cmesh);
        //Area
        //delete cmesh_mult;

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
//        const char mychar = "SolH1";
//        const char* name = &mychar;
        //const char *name = "SolH1.txt";

        //cmesh->ElementSolution().Print(name);
        //cmesh->ElementSolution().Print();
        //TPZSolutionMatrix &solmat=cmesh->ElementSolution();

    }
}


int main (){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    Hdiv_MixedCT();
    return 0;
}

std::map<int, TPZVec<STATE>> IntegralFluxoNormalCmult(const std::set<int>& bcMatId, TPZMultiphysicsCompMesh *cmesh){
    TPZGeoMesh* gmesh = cmesh->Reference();
    //const REAL zerotol = ZeroTolerance();
    std::map<int, TPZVec<STATE>> result;
    std::set<int> matidsInlet,matidsOutlet;
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    int nels=cmesh->NElements();
    //REAL area=0.0;
    REAL flux=0.0;
    std::cout<<"Integral Malha Multifisica"<<std::endl;
    std::cout<<" "<<std::endl;
    for(auto it: bcMatId){
        
        std::set<int> matid = {it};
        TPZVec<STATE> vecflux=cmesh->Integrate("BcNormalFlux",matid);
        //TPZVec<STATE> vecvol=cmesh->Integrate("Volume",matid);

        result[it] =vecflux;
        std::cout<<"Integral of normal flux for matid "<<it<<" = "<<vecflux<<std::endl;
        //std::cout<<"Volume for matid "<<it<<" = "<<vecvol<<std::endl;
        if (vecflux[0]>1e-5){
            //area=vecvol[0];
            flux=vecflux[0];
        }
        
       
        
    }
    REAL area=200.0;
    REAL deltap=100.0;
    REAL Lx=10.0;
    std::cout<<"Kequivalente "<<" = "<<(flux*Lx)/(area*deltap)<<std::endl;

    return result;
};

std::map<int, TPZVec<STATE>> IntegralFluxoNormalCmesh(const std::set<int>& bcMatId, TPZCompMesh *cmesh){
    TPZGeoMesh* gmesh = cmesh->Reference();
    //const REAL zerotol = ZeroTolerance();
    std::map<int, TPZVec<STATE>> result;
    std::set<int> matidsInlet,matidsOutlet;
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    int nels=cmesh->NElements();
    REAL area=0.0;
    REAL flux=0.0;
    std::cout<<"Integral Malha Computacional"<<std::endl;
    std::cout<<" "<<std::endl;

    for(auto it: bcMatId){
        
        std::set<int> matid = {it};
        TPZVec<STATE> vecflux=cmesh->Integrate("Pressure",matid);
        //TPZVec<STATE> vecvol=cmesh->Integrate("Volume",matid);

        result[it] =vecflux;
        std::cout<<"Integral of normal flux for matid "<<it<<" = "<<vecflux<<std::endl;
        //std::cout<<"Volume for matid "<<it<<" = "<<vecvol<<std::endl;
        if (it==4){
            //area=vecvol[0];
            flux=vecflux[0];
        }
       
        
    }
    std::cout<<"Kequivalente "<<" = "<<(flux)/(90)<<std::endl;

    return result;
};


//REAL IntegrateNormalFlux(TPZCompMesh *cmesh, int bc_matid){
//    REAL fluxo= 0.0;
//    int64_t nel=cmesh->NElements();
//    for(int64_t el=0;el<nel;el++){
//        TPZCompEl *cel =cmesh->Element(el);
//        if(!cel) continue;
//        TPZGeoEl *gel=cel->Reference();
//        if(!gel) continue;
//        if(gel->MaterialId()!=bc_matid)
//            continue;
//        TPZIntPoints *intrule=gel->CreateSideIntegrationRule(gel->NSides()-1, 4);
//        int npoints=intrule->NPoints();
//        TPZManVector<REAL,3> qsi(gel->Dimension());
//        TPZManVector<REAL,3> normal(3,0.);
//        TPZFNMatrix<9,REAL> jac,jacinv,axes;
//        REAL detjac;
//        for(int ip=0; ip<npoints; ip++)
//        {
//            REAL weight;
//            intrule->Point(ip,qsi,weight);
//            gel->Jacobian(qsi,jac,axes,detjac,jacinv);
            //
            // Obtener solución H(div)
                        //
            //TPZMaterialData data;
            //cel->InitMaterialData(data);
            //cel->ComputeRequiredData(data,qsi);
                        //
                        // data.sol[0] contiene el flujo
                        //
            //TPZVec<STATE> &sol = data.sol[0];
            //TPZGeoElSide gelside(gel,side);
            //gelside.Normal(qsi, normal);
//        }
//
//    }
//};
