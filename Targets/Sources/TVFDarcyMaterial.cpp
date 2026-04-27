#include "TVFDarcyMaterial.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"
//#include "TPZLapack.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.material.darcy");
#endif

#define USEBLAS
#undef USEBLAS
TVFDarcyMaterial::TVFDarcyMaterial() : TPZRegisterClassId(&TVFDarcyMaterial::ClassId),
                                         TBase(), fDim(-1) {}

TVFDarcyMaterial::TVFDarcyMaterial(int id, int dim) : TPZRegisterClassId(&TVFDarcyMaterial::ClassId),
                                                        TBase(id), fDim(dim)
{
}

/**
         copy constructor
 */
TVFDarcyMaterial::TVFDarcyMaterial(const TVFDarcyMaterial &copy) : TBase(copy), fDim(copy.fDim)
{
    *this = copy;
}
/**
         copy constructor
 */
TVFDarcyMaterial &TVFDarcyMaterial::operator=(const TVFDarcyMaterial &copy)
{
    TBase::operator=(copy);
    fDim = copy.fDim;
    return *this;
}


void TVFDarcyMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef) {
    
    if (datavec[0].fShapeType == TPZMaterialData::EEmpty) {
        return;
    }

    // if(AssembleResidualOnly) {
    //     ContributeResidual(datavec, weight, ef);
    //     return;
    // }

    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].fDeformedDirections;
    TPZFMatrix<REAL> &phip = datavec[1].phi;
    TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
    TPZFMatrix<REAL> &divQ = datavec[0].divphi;
    TPZFNMatrix<9, REAL> dphiPXY(3, dphiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPXY, datavec[1].axes);

    int nphiQ, nphiP;
    nphiP = phip.Rows();
    nphiQ = datavec[0].fDeformedDirections.Cols();


    // Getting solution from the current step/iteration
    TPZVec<STATE> vQsol(3, 0.0);
    REAL vpsol = 0;
    REAL vdivQsol = 0;

    if(fNoLinearContext){
        vQsol = datavec[0].sol[0];
        vpsol = datavec[1].sol[0][0];
        vdivQsol = datavec[0].divsol[0][0];
    }

    TPZFMatrix<STATE> Qsol(vQsol.size(), 1);
    Qsol(0, 0) = vQsol[0];
    Qsol(1, 0) = vQsol[1];
    Qsol(2, 0) = vQsol[2];
    TPZFNMatrix<1, REAL> psol(1, 1, vpsol);
    TPZFNMatrix<1, REAL> divQsol(1, 1, vdivQsol);


    int nactive = 0;
    for (const auto &i : datavec) {
        if (i.fActiveApproxSpace) {
            nactive++;
        }
    }

    TPZFNMatrix<3,REAL> ForceTerm(fDim, 1, 0.0);
    TPZManVector<STATE> res(1);
    STATE force = 0;
    if (this->HasForcingFunction())
    {
        fForcingFunction(datavec[1].x, res);
        force = res[0];
        for (int64_t i = 0; i < fDim; i++){
            ForceTerm(i,0) = force;
        }
    }

    TPZFNMatrix<3,REAL> SourceTerm(3, 1, 0.0); //! Make as a function of fDim, in order to do so, only get fDim in phiQ vector funcs
    TPZVec<REAL> sourceAux(3);
    if (this->HasForcingFunction())
    {
        this->ForcingFunction()(datavec[0].x, sourceAux);
        for (int64_t i = 0; i < fDim; i++)
        {
            SourceTerm(i,0) = sourceAux[i];
        }
    }

    const STATE perm = GetPermeability(datavec[0].x);
    const STATE inv_perm = 1 / perm;

    // REAL axiFactor = 1.0;
    // if (fIsAxisymmetric) // Axisymmetric: assuming radius is aligned with the x axis
    // {
    // REAL r = datavec[0].x[0];
    // axiFactor = 1.0 / (2.0 * M_PI * r);
    // }
    

    //* Tangent Matrix
    //Calculate the matrix contribution for flux. Matrix A
    REAL factor = weight * inv_perm;
    ek.AddContribution(0, 0, phiQ, 1, phiQ, 0, factor); //A

    // Coupling terms between flux and pressure. Matrix B
    ek.AddContribution(0, nphiQ, divQ, 0, phip, 1, -weight); //B
    ek.AddContribution(nphiQ, 0, phip, 0, divQ, 1, -weight); //B^T

    //Source Term related to the flux equation
    ef.AddContribution(0, 0, phiQ, 1, SourceTerm, 0, weight);

    //Source Term related to the pressure equation //! See this
    //ef.AddContribution(nphiQ, 0, phip, 1, ForceTerm, 0, weight); // RHS


    //* Residual Vector
    ef.AddContribution(0, 0, phiQ, 1, Qsol, 0, -factor); // constitutive equation
    ef.AddContribution(0, 0, divQ, 0, psol, 0, weight); // constitutive equation 
    ef.AddContribution(nphiQ, 0, phip, 0, divQsol, 0, weight); // conservation equation
}


void TVFDarcyMaterial::ContributeResidual(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef) {
    return;
}


void TVFDarcyMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

    int dim = Dimension();

    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    int phrq = phiQ.Rows();

    REAL v2 = bc.Val2()[0];
    REAL v1 = bc.Val1()(0, 0);
    REAL u_D = 0;
    REAL normflux = 0.;

    TPZVec<STATE> Qsol(3, 0.0);
    REAL psol = 0;
    REAL divQsol = 0;

    if(fNoLinearContext){
        Qsol = datavec[0].sol[0];
        //psol = datavec[1].sol[0][0];
        //divQsol = datavec[0].divsol[0][0];
    }


    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(3, 1);
        bc.ForcingFunctionBC()(datavec[0].x, res, gradu);

        const STATE perm = GetPermeability(datavec[0].x);

        for (int i = 0; i < 3; i++) {
            normflux += datavec[0].normal[i] * perm * gradu(i, 0);
        }

        if (bc.Type() == 0 || bc.Type() == 4) {
            v2 = res[0];
            u_D = res[0];
            normflux *= (-1.);
        } else if (bc.Type() == 1 || bc.Type() == 2) {
            v2 = -normflux;
            if (bc.Type() == 2) {
                v2 = -res[0] + v2 / v1;
            }
        } else if (bc.Type() == 5) {
            v2 = res[0];
        } else {
            DebugStop();
        }
    } else {
        v2 = bc.Val2()[0];
    }

    switch (bc.Type()) {
        case 0 :    // Dirichlet condition
            for (int iq = 0; iq < phrq; iq++) { 
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq, 0) += (-1.) * v2 * phiQ(iq, 0) * weight;
            }
            break;

        case 1 :   // Neumann condition imposed in the initial solution
            for (int iq = 0; iq < phrq; iq++) {
                ef(iq, 0) += TPZMaterial::fBigNumber * (v2-Qsol[0]) * phiQ(iq, 0) * weight;
                for (int jq = 0; jq < phrq; jq++) {

                    ek(iq, jq) += TPZMaterial::fBigNumber * phiQ(iq, 0) * phiQ(jq, 0) * weight;
                }
            }
            break;

        case 2 :            // mixed condition
            for (int iq = 0; iq < phrq; iq++) {
                ef(iq, 0) += v2 * phiQ(iq, 0) * weight;
                for (int jq = 0; jq < phrq; jq++) {
                    ek(iq, jq) += weight / v1 * phiQ(iq, 0) * phiQ(jq, 0);
                }
            }
            break;

        case 4:
            //this case implemented the general Robin boundary condition
            // sigma.n = Km(u-u_D)+g
            //val1(0,0) = Km
            //val2(1,0) = g
            if (IsZero(bc.Val1()(0, 0))) {

                for (int iq = 0; iq < phrq; iq++) {
                    ef(iq, 0) += TPZMaterial::fBigNumber * normflux * phiQ(iq, 0) * weight;
                    for (int jq = 0; jq < phrq; jq++) {
                        ek(iq, jq) += TPZMaterial::fBigNumber * phiQ(iq, 0) * phiQ(jq, 0) * weight;
                    }
                }

            } else {

                REAL InvKm = 1. / bc.Val1()(0, 0);
                REAL g = normflux;
                for (int in = 0; in < phiQ.Rows(); in++) {
                    //<(InvKm g - u_D)*(v.n)
                    ef(in, 0) += (STATE) (InvKm * g - u_D) * phiQ(in, 0) * weight;
                    for (int jn = 0; jn < phiQ.Rows(); jn++) {
                        //InvKm(sigma.n)(v.n)
                        ek(in, jn) += (STATE) (InvKm * phiQ(in, 0) * phiQ(jn, 0) * weight);
                    }
                }
            }

            break;

        case 5:
            TPZFMatrix<REAL> &phi = datavec[0].fH1.fPhi;
            for (int in = 0; in < phi.Rows(); in++) {
                //<(InvKm g - u_D)*(v.n)
                ef(in, 0) += TPZMaterial::fBigNumber * v2 * phi(in, 0) * weight;
                for (int jn = 0; jn < phi.Rows(); jn++) {
                    //InvKm(sigma.n)(v.n)
                    ek(in, jn) += TPZMaterial::fBigNumber * phi(in, 0) * phi(jn, 0) * weight;
                }
            }

    }
}

void TVFDarcyMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &solOut) {
    solOut.Resize(this->NSolutionVariables(var));
    solOut.Fill(0.);
    TPZManVector<STATE, 10> SolP, SolQ;
    const STATE perm = GetPermeability(datavec[0].x);
    const STATE inv_perm = 1 / perm;

    // SolQ = datavec[0].sol[0];
    SolP = datavec[1].sol[0];
    if(SolP.size() == 0) SolP.Resize(1,0.);

    //TODO
    if(datavec[0].fShapeType == TPZMaterialData::EEmpty) {
        if (var == 2) solOut[0] = SolP[0];
        return;
    }

    if (var == 1) { //function (state variable Q)
        for (int i = 0; i < fDim; i++) {
            solOut[i] = datavec[0].sol[0][i];

        }
        return;
    }

    if (var == 2) {
        solOut[0] = SolP[0];//function (state variable p)
        return;
    }

    if (var == 3) {
        solOut[0] = datavec[0].dsol[0](0, 0);
        solOut[1] = datavec[0].dsol[0](1, 0);
        solOut[2] = datavec[0].dsol[0](2, 0);
        return;
    }

    if (var == 4) {
        solOut[0] = datavec[0].dsol[0](0, 1);
        solOut[1] = datavec[0].dsol[0](1, 1);
        solOut[2] = datavec[0].dsol[0](2, 1);
        return;
    }

    if (var == 5) {
        solOut[0] = datavec[0].divsol[0][0];
        return;
    }

    // Exact solution
    if (var == 6) {
        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> flux(3, 1);
        if (fExactSol) {
            fExactSol(datavec[0].x, exactSol, flux);
        }
        solOut[0] = exactSol[0];
        return;
    } // var6

    if (var == 7) {

        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> gradu(fDim, 1);

        if (fExactSol) {
            fExactSol(datavec[0].x, exactSol, gradu);
        }

        for (int i = 0; i < fDim; i++) {
            solOut[i] = -perm * gradu(i, 0);
        }

        return;
    } // var7

    if (var == 8) {
        solOut[0] = datavec[1].p;
        return;
    }

    if (var == 9) {

        if(datavec[1].fShapeType == TPZMaterialData::EEmpty) return;
        TPZFNMatrix<9, REAL> dsoldx(3, 1.,0.);
        TPZFNMatrix<9, REAL> dsoldaxes(fDim, 1,0.);

        dsoldaxes = datavec[1].dsol[0];
        TPZAxesTools<REAL>::Axes2XYZ(dsoldaxes, dsoldx, datavec[1].axes);

        for (int i = 0; i < fDim; i++) {
            solOut[i] = dsoldx(i, 0);
        }

        return;
    }

    if (var == 10) {
        solOut[0] = 0.;
        // solOut[0]=datavec[0].dsol[0](0,0)+datavec[0].dsol[0](1,1);
        solOut[0] = datavec[0].divsol[0][0];
        // for (int j = 0; j < fDim; j++) {
        //     solOut[0] += datavec[0].dsol[0](j, j);
        // }
        return;
    }

    if (var == 11) {
        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> flux(3, 1);
        fExactSol(datavec[0].x, exactSol, flux);
        solOut[0] = flux(2, 0);
        return;
    }
    if (var == 12) {
        for (int i = 0; i < fDim; i++) {
            solOut[i] = 0.;
        }
        for (int i = 0; i < fDim; i++) {
            solOut[i] -= inv_perm * datavec[0].sol[0][i];
        }
        return;
    }
    if (var == 13) {
        solOut[0] = perm;
        return;
    }

    if (datavec.size() == 4) {
        if (var == 14) {
            solOut[0] = datavec[2].sol[0][0];
            return;
        }
        if (var == 15) {
            solOut[0] = datavec[3].sol[0][0];
            return;
        }

    }

    if (var == 16) { //ExactFluxShiftedOrigin
        // Solution EArcTan returns NAN for (x,y) == (0,0). Replacing data.x by
        // inf solves this problem.
        STATE infinitesimal = 0.0000000001;
        TPZManVector<REAL, 3> inf = {infinitesimal, infinitesimal, infinitesimal};

        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> gradu(3, 1);

        if (fExactSol) {
            if (datavec[0].x[0] == 0. && datavec[0].x[1] == 0.) {
                fExactSol(inf, exactSol, gradu);
            } else {
                fExactSol(datavec[0].x, exactSol, gradu);
            }
        }
        for (int i = 0; i < 3; i++) {
            solOut[i] = -perm * gradu(i, 0);
        }

        return;
    }

    if (var == 17) {
        TPZVec<STATE> divsigma(1, 0.);
        if (fForcingFunction) {
            fForcingFunction(datavec[0].x, divsigma);
        }
        solOut[0] = divsigma[0];
        return;
    }
}

void TVFDarcyMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {

    /**
     * datavec[0]= Flux
     * datavec[1]= Pressure
     *
     * Errors:
     * [0] L2 for pressure
     * [1] L2 for flux
     * [2] L2 for div(flux)
     * [3] Grad pressure (Semi H1)
     * [4] Hdiv norm
    **/
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    TPZManVector<STATE, 3> fluxfem(3), pressurefem(1,0);
    fluxfem = data[0].sol[0];
    STATE divsigmafem = data[0].divsol[0][0];

    auto dsol = data[1].dsol;

    TPZManVector<STATE,1> divsigma(1,0.);

    TPZManVector<STATE,1> u_exact(1, 0);
    TPZFMatrix<STATE> du_exact(3, 1, 0);
    if (this->fExactSol) {
        this->fExactSol(data[0].x, u_exact, du_exact);
    }
    if (this->fForcingFunction) {
        this->fForcingFunction(data[0].x, divsigma);
    }

    REAL residual = (divsigma[0] - divsigmafem) * (divsigma[0] - divsigmafem);
    if(data[1].sol[0].size())
        pressurefem[0] = data[1].sol[0][0];

    const STATE perm = GetPermeability(data[0].x);
    const STATE inv_perm = 1 / perm;

    TPZManVector<STATE, 3> gradpressurefem(3, 0.);
    this->Solution(data, VariableIndex("GradPressure"), gradpressurefem);

    TPZManVector<STATE, 3> fluxexact(3, 0);
    TPZManVector<STATE, 3> gradpressure(3, 0);
    for (int i = 0; i < 3; i++) {
        gradpressure[i] = du_exact[i];
        fluxexact[i] = -perm * gradpressure[i];
    }

    REAL L2flux = 0., L2grad = 0.;
    for (int i = 0; i < 3; i++) {
        L2flux += (fluxfem[i] - fluxexact[i]) * inv_perm * (fluxfem[i] - fluxexact[i]);
        L2grad += (du_exact[i] - gradpressurefem[i]) * (du_exact[i] - gradpressurefem[i]);
    }
    errors[0] = (pressurefem[0] - u_exact[0]) * (pressurefem[0] - u_exact[0]);//L2 error for pressure
    errors[1] = L2flux;//L2 error for flux
    errors[2] = residual;//L2 for div
    errors[3] = L2grad;
    errors[4] = L2flux + residual;
#ifdef PZ_LOG
    if(logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "x " << data[0].x << " fluxfem " << fluxfem << " fluxexact " << fluxexact;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

int TVFDarcyMaterial::VariableIndex(const std::string &name) const {
    if (!strcmp("Flux", name.c_str())) return 1;
    if (!strcmp("Pressure", name.c_str())) return 2;
    if (!strcmp("GradFluxX", name.c_str())) return 3;
    if (!strcmp("GradFluxY", name.c_str())) return 4;
    if (!strcmp("DivFlux", name.c_str())) return 5;
    if (!strcmp("ExactPressure", name.c_str())) return 6;
    if (!strcmp("ExactFlux", name.c_str())) return 7;
    if (!strcmp("POrder", name.c_str())) return 8;
    if (!strcmp("GradPressure", name.c_str())) return 9;
    if (!strcmp("Divergence", name.c_str())) return 10;
    if (!strcmp("ExactDiv", name.c_str())) return 11;
    if (!strcmp("Derivative", name.c_str())) return 12;
    if (!strcmp("Permeability", name.c_str())) return 13;
    if (!strcmp("g_average", name.c_str())) return 14;
    if (!strcmp("u_average", name.c_str())) return 15;
    if (!strcmp("ExactFluxShiftedOrigin", name.c_str())) return 16;
    if (!strcmp("EstimatedError", name.c_str())) return 100;
    if (!strcmp("TrueError", name.c_str())) return 101;
    if (!strcmp("EffectivityIndex", name.c_str())) return 102;
    if (!strcmp("ExactDivSigma", name.c_str())) return 17;
    DebugStop();
    return -1;
}

int TVFDarcyMaterial::NSolutionVariables(int var) const {
    if (var == 1) return 3;
    if (var == 2) return 1;
    if (var == 3) return 3;
    if (var == 4) return 3;
    if (var == 5) return 1;
    if (var == 6) return 1;
    if (var == 7) return 3;
    if (var == 8) return 1;
    if (var == 9) return 3;
    if (var == 10 || var == 11) return 1;
    if (var == 12) return 3;
    if (var == 13) return 1;
    if (var == 14) return 1;
    if (var == 15) return 1;
    if (var == 16) return 3;
    if (var == 17) return 1;
    if (var == 100) return 1;
    if (var == 101) return 1;
    if (var == 102) return 1;

    DebugStop();
    return -1;
}

void TVFDarcyMaterial::SetDimension(int dim) {
    if (dim > 3 || dim < 1) DebugStop();
    fDim = dim;
}

int TVFDarcyMaterial::ClassId() const {
    return Hash("TVFDarcyMaterial") ^ TBase::ClassId() << 1;
}

void TVFDarcyMaterial::Print(std::ostream &out) const {
    out << "Material Name: " << this->Name() << "\n";
    out << "Material Id: " << this->Id() << "\n";
    out << "Dimension: " << this->Dimension() << "\n\n";
}

void TVFDarcyMaterial::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    int nref = datavec.size();
    for (int i = 0; i < nref; i++) {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsNeighborSol = false;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = false;
        datavec[i].fNeedsHSize = false;
        if(fNoLinearContext){
            datavec[i].fNeedsSol = true;
        }
    }
}

void TVFDarcyMaterial::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    // default is no specific data requirements
    int nref = datavec.size();
    for (int iref = 0; iref < nref; iref++) {
        datavec[iref].SetAllRequirements(false);
        // datavec[iref].fNeedsSol = false;
    }
    
    if(fNoLinearContext){
        datavec[0].fNeedsSol = true;
        datavec[1].fNeedsSol = true;
    }

    if (type == 50) {
        for (int iref = 0; iref < nref; iref++) {
            // datavec[iref].fNeedsSol = false;
        }
    }
}
#undef USEBLAS
