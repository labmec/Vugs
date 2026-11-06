#include "sources.h"
#include "pzlog.h"

void CasePaper();
void Hdiv_MixedCT();

int main (){
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    CasePaper();
//    Hdiv_MixedCT();
    return 0;
}

void CasePaper(){
    
    
    ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/FewFractures.json");
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/NewInputs/2Dsem.json");
    TPZGeoMesh *gmesh = new TPZGeoMesh;

    std::string meshName = inputData.MeshName();
    std::string filename = inputData.MeshFile();

    int problemDim = inputData.dim();
    int approxType = inputData.approxType();
    int problemType = inputData.problemType();
    int pressOrder = inputData.pressOrder();
    std::map<int, TPZVec<STATE>> IntegralFluxoNormalCmult(const std::set<int>& bcMatId, TPZMultiphysicsCompMesh *cmesh);
    std::map<int, TPZVec<STATE>> IntegralFluxoNormalCmesh(const std::set<int>& bcMatId, TPZCompMesh *cmesh);
    void IntegrateElementVelocity(TPZCompEl* cel, TPZVec<STATE>& intFlux);
    REAL IntegrateElementEnergyCross(TPZCompEl* celH1_k, TPZCompEl* celH1_l);
    
    //REAL CalcularVelocidadeMedia(TPZCompMesh* cmesh_flux, const int direcao_alvo);


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
        TPZMultiphysicsCompMesh *cmesh_mult = CreateMultiMeshPaper(gmesh, meshvec, inputData);

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
            TPZManVector<std::string, 2> fields = {"Flux", "Pressure", "GradFluxX","GradFluxY"};
            auto vtk = TPZVTKGenerator(cmesh_mult, fields, plotfile, vtkRes);
            vtk.Do();
            
            cmesh_mult->SetDimModel(1);
            auto vtk1 = TPZVTKGenerator(cmesh_mult, fields, plotfile, vtkRes);
            vtk1.Do();
        }
    
        PrintCompMesh(cmesh_mult);
        PrintCompMesh(Flux_cmesh);
        PrintCompMesh(Pressure_cmesh);
        std::set<int> bcids={2,3,4};
        //std::set<int> bcids={2,3,4,8,9,11};

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
        
        
        std::string plotfile2 = meshName + approxName + "FRAC.vtk";
  
        Analisys->DefineGraphMesh(problemDim-1, scalnames, vecnames, plotfile2);
        
        Analisys->PostProcess(ref, problemDim-1);
//        const char mychar = "SolH1";
//        const char* name = &mychar;
        //const char *name = "SolH1.txt";

        //cmesh->ElementSolution().Print(name);
        //cmesh->ElementSolution().Print();
        //TPZSolutionMatrix &solmat=cmesh->ElementSolution();

    }
    
}
void Hdiv_MixedCT(){
 
    //ReadJson inputData("/home/marina/programming/Stokes-Darcy-Research/VUGS/Inputs/FewVugs.json");
    ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/NewInputs/ReproducedCase.json");

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
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/NewInputs/CILAMCE/SingleVug3D_FluxZf6.json");//RUN <---

    //VUG - High Permeability
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/SingleVug3D_H1HP.json");//RUN

    //
    //FRATURA
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/FewFractures.json");//
    //ReadJson inputData ("Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/SingleFrature3D_H1.json");
    //ReadJson inputData("/Users/victorvillegassalabarria/Documents/Github/Vugs/Inputs/SingleFrature3D_H1.json");
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
    void IntegrateElementVelocity(TPZCompEl* cel, TPZVec<STATE>& intFlux);
    REAL IntegrateElementEnergyCross(TPZCompEl* celH1_k, TPZCompEl* celH1_l);
    
    //REAL CalcularVelocidadeMedia(TPZCompMesh* cmesh_flux, const int direcao_alvo);


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
        std::set<int> bcids={2,3,4};
        //std::set<int> bcids={2,3,4,8,9,11};

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
        
        std::string plotfile2 = meshName + approxName + "Frac.vtk";
  
        Analisys->DefineGraphMesh(problemDim, scalnames, vecnames, plotfile2);
        
        Analisys->PostProcess(ref, problemDim-1);
//        const char mychar = "SolH1";
//        const char* name = &mychar;
        //const char *name = "SolH1.txt";

        //cmesh->ElementSolution().Print(name);
        //cmesh->ElementSolution().Print();
        //TPZSolutionMatrix &solmat=cmesh->ElementSolution();

    }
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
        
        std::set<int> matid = {700};
        TPZVec<STATE> vecflux=cmesh->Integrate("Flux",matid);
        
        //TPZVec<STATE> vecvol=cmesh->Integrate("Volume",matid);

        result[it] =vecflux;
        std::cout<<"Integral of normal flux x for matid "<<it<<" = "<<vecflux[0]<<std::endl;
        std::cout<<"Integral of normal flux y for matid "<<it<<" = "<<vecflux[1]<<std::endl;
        std::cout<<"Integral of normal flux z for matid "<<it<<" = "<<vecflux[2]<<std::endl;
        
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

void IntegrateElementVelocity(TPZCompEl* cel, TPZVec<STATE>& intFlux) {
    if (!cel) return;

    TPZGeoEl *gel = cel->Reference();
    if (!gel) return;

    int dim = cel->Dimension();
    intFlux.Fill(0.0);

    // Regla de integración sobre el elemento
    TPZIntPoints *intrule = gel->CreateSideIntegrationRule(gel->NSides() - 1, 4);
    int npoints = intrule->NPoints();

    TPZVec<REAL> xi(dim, 0.0);
    TPZFNMatrix<9, REAL> jac(dim, dim), jacinv(dim, dim), axes(dim, 3);
    REAL detjac;
    TPZVec<STATE> solFlux(dim, 0.0);

    for (int point = 0; point < npoints; point++) {
        REAL weight;
        intrule->Point(point, xi, weight);
        gel->Jacobian(xi, jac, axes, detjac, jacinv);

        // VarIndex 1 en TPZDarcyFlow / Hdiv suele ser el Vector de Flujo (Darcy Velocity)
        cel->Solution(xi, 1, solFlux);

        REAL weightDetJac = weight * std::fabs(detjac);
        for (int d = 0; d < dim; d++) {
            intFlux[d] += solFlux[d] * weightDetJac;
        }
    }
    delete intrule;
}
// Integra (Grad p_k)^T * Perm * (Grad p_l) en un elemento
REAL IntegrateElementEnergyCross(TPZCompEl* celH1_k, TPZCompEl* celH1_l) {
    if (!celH1_k || !celH1_l) return 0.0;

    TPZGeoEl *gel = celH1_k->Reference();
    if (!gel) return 0.0;

    int dim = celH1_k->Dimension();
    TPZMaterial *matH1 = celH1_k->Material();
    TPZDarcyFlow *matDarcy = dynamic_cast<TPZDarcyFlow*>(matH1);
    if (!matDarcy) DebugStop();

    TPZIntPoints *intrule = gel->CreateSideIntegrationRule(gel->NSides() - 1, 4);
    int npoints = intrule->NPoints();

    TPZVec<REAL> xi(dim, 0.0);
    TPZFNMatrix<9, REAL> jac(dim, dim), jacinv(dim, dim), axes(dim, 3);
    REAL detjac;

    TPZVec<STATE> gradP_k(dim, 0.0);
    TPZVec<STATE> gradP_l(dim, 0.0);

    REAL elementEnergy = 0.0;

    for (int point = 0; point < npoints; point++) {
        REAL weight;
        intrule->Point(point, xi, weight);
        gel->Jacobian(xi, jac, axes, detjac, jacinv);

        // VarIndex 2 corresponde al Gradiente de Presión Grad(p)
        celH1_k->Solution(xi, 2, gradP_k);
        celH1_l->Solution(xi, 2, gradP_l);

        const STATE perm = matDarcy->GetPermeability(xi);

        // Producto escalar con la permeabilidad (Grad p_k) * K * (Grad p_l)
        REAL grad_K_grad = 0.0;
        for (int d = 0; d < dim; d++) {
            grad_K_grad += gradP_k[d] * perm * gradP_l[d];
        }

        elementEnergy += grad_K_grad * weight * std::fabs(detjac);
    }
    delete intrule;
    return elementEnergy;
}
struct EffectivePermeabilityResult {
    TPZFNMatrix<4, REAL> K_eff;     // Tensor directo 2x2
    TPZFNMatrix<4, REAL> K_eff_star;// Tensor simétrico 2x2
    REAL alpha_star;                // Razón de anisotropía
    REAL beta_star;                 // Ángulo principal (radianes)
};

EffectivePermeabilityResult ComputeEffectivePermeability(
    TPZVec<TPZCompEl*>& celVecH1_sim1, TPZVec<TPZCompEl*>& celVecHdiv_sim1, // Simulación k=1 (Flujo X)
    TPZVec<TPZCompEl*>& celVecH1_sim2, TPZVec<TPZCompEl*>& celVecHdiv_sim2, // Simulación k=2 (Flujo Y)
    std::set<int> domainMatIds,
    REAL domainVolume)
{
    EffectivePermeabilityResult res;
    res.K_eff.Redim(2, 2);
    res.K_eff_star.Redim(2, 2);

    int nels = celVecH1_sim1.size();
    int dim = 2; // Para problemas 2D (Ejemplo 5.1)

    TPZVec<STATE> totalFlux_sim1(dim, 0.0);
    TPZVec<STATE> totalFlux_sim2(dim, 0.0);

    REAL energy_11 = 0.0;
    REAL energy_22 = 0.0;
    REAL energy_12 = 0.0;

    for (int el = 0; el < nels; el++) {
        TPZCompEl* celH1_1 = celVecH1_sim1[el];
        TPZCompEl* celHdiv_1 = celVecHdiv_sim1[el];
        TPZCompEl* celH1_2 = celVecH1_sim2[el];
        TPZCompEl* celHdiv_2 = celVecHdiv_sim2[el];

        if (!celH1_1 || !celH1_2) continue;

        int elMatId = celH1_1->Reference()->MaterialId();
        if (domainMatIds.find(elMatId) == domainMatIds.end()) continue;

        // -------------------------------------------------------------
        // A) Método Directo: Integración de Velocidades K_eff
        // -------------------------------------------------------------
        TPZVec<STATE> elFlux1(dim, 0.0), elFlux2(dim, 0.0);
        IntegrateElementVelocity(celHdiv_1 ? celHdiv_1 : celH1_1, elFlux1);
        IntegrateElementVelocity(celHdiv_2 ? celHdiv_2 : celH1_2, elFlux2);

        for (int d = 0; d < dim; d++) {
            totalFlux_sim1[d] += elFlux1[d];
            totalFlux_sim2[d] += elFlux2[d];
        }

        // -------------------------------------------------------------
        // B) Método Variacional Simétrico: Integración de Energía K_eff*
        // -------------------------------------------------------------
        energy_11 += IntegrateElementEnergyCross(celH1_1, celH1_1);
        energy_22 += IntegrateElementEnergyCross(celH1_2, celH1_2);
        energy_12 += IntegrateElementEnergyCross(celH1_1, celH1_2);
    }

    // 1. Matriz Directa K_eff = - (1/|Omega|) * Integral(u)
    res.K_eff(0, 0) = -totalFlux_sim1[0] / domainVolume; // K11
    res.K_eff(1, 0) = -totalFlux_sim1[1] / domainVolume; // K21
    res.K_eff(0, 1) = -totalFlux_sim2[0] / domainVolume; // K12
    res.K_eff(1, 1) = -totalFlux_sim2[1] / domainVolume; // K22

    // 2. Matriz Simétrica K_eff* = (1/|Omega|) * BilinearEnergy
    res.K_eff_star(0, 0) = energy_11 / domainVolume;
    res.K_eff_star(1, 1) = energy_22 / domainVolume;
    res.K_eff_star(0, 1) = energy_12 / domainVolume;
    res.K_eff_star(1, 0) = energy_12 / domainVolume;

    // 3. Autovalores, Anisotropía (alpha*) y Ángulo Principal (beta*)
    REAL a = res.K_eff_star(0, 0);
    REAL b = res.K_eff_star(0, 1);
    REAL c = res.K_eff_star(1, 1);

    REAL trace = a + c;
    REAL diff = a - c;
    REAL rad = std::sqrt(diff * diff + 4.0 * b * b);

    REAL lambda_max = (trace + rad) / 2.0;
    REAL lambda_min = (trace - rad) / 2.0;

    res.alpha_star = lambda_max / lambda_min;
    res.beta_star = 0.5 * std::atan2(2.0 * b, diff);

    return res;
}

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
