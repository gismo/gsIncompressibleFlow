/** @file uwbINS2SolverTestPruh.cpp

Author(s): J. Sourek
*/

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include <gismo.h>

#include "uwbINSSolverSteady.h"
#include "uwbINSSolverDecoupledInterface.h"
#include "uwbINSSolverUnsteady.h"

using namespace gismo;

std::string omegaXr_component(real_t omega, std::string variable);
gsBoundaryConditions<> defineBCs(gsFunctionExpr<> & Uin, gsFunctionExpr<> & Uwall, gsFunctionExpr<> & Ublade, char method);
gsBoundaryConditions<> defineBCs_nonPer(gsFunctionExpr<> & Uin, gsFunctionExpr<> & Uwall, gsFunctionExpr<> & Ublade, char method);

int main(int argc, char *argv[])
{

    // ========================================= Settings ========================================= 
    int numRefine = 0;
    real_t viscosity = 0.1;
    real_t omega = 0; // angular velocity
    real_t timeStep = 0.01;
    real_t alpha_u = 0.5; // relaxation parameter for projection
    real_t alpha_p = 1; // relaxation parameter for projection
    bool dg = false; // use DGFEM to connect patches

    gsAssemblerOptions opt;
    opt.dirStrategy = dirichlet::elimination;
    if (dg)
        opt.intStrategy = iFace::dg;
    

    gsInfo << "Solving the Kaplan turbine (stationary part) example.\n";

    // ========================================= Define problem ========================================= 
    gsBoundaryConditions<> bcInfo_C, bcInfo_P, bcInfo_nonPerC, bcInfo_nonPerP;
    gsFunctionExpr<> f("0", "0", "0", 3); // external force
    gsFunctionExpr<> *Uin, *Ublade, *Uwall; // relative velocity BCs for inlet/rotating part/stationary part
    
    // Boundary conditions
    Uin = new gsFunctionExpr<>("( -(100 / 9) * (sqrt(y^2+z^2)-0.7)^2 + 1)", "0", "0", 3);
    Uwall = new gsFunctionExpr<>("0", "0", "0", 3);

    if (omega)
        Ublade = new gsFunctionExpr<>("0", omegaXr_component(-omega, "z"), omegaXr_component(omega, "y"), 3);
    else
        Ublade = Uwall;

    bcInfo_C = defineBCs(*Uin, *Uwall, *Ublade, 'C');
    bcInfo_P = defineBCs(*Uin, *Uwall, *Ublade, 'P');
    bcInfo_nonPerC = defineBCs_nonPer(*Uin, *Uwall, *Ublade, 'C');
    bcInfo_nonPerP = defineBCs_nonPer(*Uin, *Uwall, *Ublade, 'P');

    // ========================================= Define geometry ========================================= 
    std::string input(MOTOR_DATA_DIR "/uwb-pilsen/pruh_withoutTrailingEdge.xml");

    gsFileData<> fileData(input);

    gsMultiPatch<>* patches = NULL;
    if (fileData.has< gsMultiPatch<> >())
    {
        gsInfo << "Loading geometry from " << input << "\n";
        patches = fileData.getFirst< gsMultiPatch<> >().release();
    }
    else
    {
        gsWarn << "Input file doesn't have a gsMultiPatch inside.\n";
        return -1;
    }
 
    gsInfo << *patches << "\n";

    // ========================================= Define basis ========================================= 

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> tbasis(*patches);
    for (int i = 0; i < numRefine; ++i)
        tbasis.uniformRefine();
    
    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure
    discreteBases[0].degreeElevate(1); //elevate the velocity space

    // ========================================= Define solver ========================================= 

    uwbINSPde<real_t> NSpde_C(*patches, bcInfo_C, f, viscosity);
    uwbINSPde<real_t> NSpde_P(*patches, bcInfo_P, f, viscosity);
    uwbINSPde<real_t> NSpde_nonPerC(*patches, bcInfo_nonPerC, f, viscosity);
    uwbINSPde<real_t> NSpde_nonPerP(*patches, bcInfo_nonPerP, f, viscosity);
    uwbINSSolverParams<real_t> params_C(NSpde_C, discreteBases, opt);
    uwbINSSolverParams<real_t> params_P(NSpde_P, discreteBases, opt);
    uwbINSSolverParams<real_t> params_nonPerC(NSpde_nonPerC, discreteBases, opt);
    uwbINSSolverParams<real_t> params_nonPerP(NSpde_nonPerP, discreteBases, opt);

    params_C.settings().set(constantsINS::timeStep, timeStep);
    params_nonPerC.settings().set(constantsINS::timeStep, timeStep);

    if (omega)
    {
        params_C.settings().set(constantsINS::omega, omega);
        params_P.settings().set(constantsINS::omega, omega);
        params_nonPerC.settings().set(constantsINS::omega, omega);
        params_nonPerP.settings().set(constantsINS::omega, omega);
    }

    // steady coupled solver
    uwbINSSolverSteady<real_t> navStokes_steady(params_C); 

    // decoupled1 solver
    params_nonPerP.settings().setProjVersion(decoupled::proj1);
    params_nonPerP.settings().set(constantsINS::alpha_u, alpha_u);
    params_nonPerP.settings().set(constantsINS::alpha_p, alpha_p);
    uwbINSSolverDecoupledInterface<real_t> navStokes_decoupled1(params_nonPerP);

    // decoupled2 solver
    params_nonPerP.settings().setProjVersion(decoupled::proj2);
    uwbINSSolverDecoupledInterface<real_t> navStokes_decoupled2(params_nonPerP);

    // decoupled1 solver, iterative for periodic conditions
    params_P.settings().setProjVersion(decoupled::proj1);
    params_P.settings().setDecoupledMethod(decoupled::iterative);
    params_P.settings().set(constantsINS::alpha_u, alpha_u);
    params_P.settings().set(constantsINS::alpha_p, alpha_p);
    uwbINSSolverDecoupledInterface<real_t> navStokes_decoupled1Iter(params_P);

    // decoupled1 solver, coupled for periodic conditions
    params_P.settings().setDecoupledMethod(decoupled::coupled);
    uwbINSSolverDecoupledInterface<real_t> navStokes_decoupled1Coup(params_P);

    // decoupled2 solver, iterative for periodic conditions
    params_P.settings().setProjVersion(decoupled::proj2);
    params_P.settings().setDecoupledMethod(decoupled::iterative);
    uwbINSSolverDecoupledInterface<real_t> navStokes_decoupled2Iter(params_P);

    // decoupled2 solver, coupled for periodic conditions
    params_P.settings().setDecoupledMethod(decoupled::coupled);
    uwbINSSolverDecoupledInterface<real_t> navStokes_decoupled2Coup(params_P);

    uwbINSSolverUnsteady<real_t> navStokes_unsteady(params_nonPerC);

    // ========================================= Solving ========================================= 
    
    gsInfo << "numDofs coupled periodic: " << navStokes_steady.numDofs() << "\n";
    gsInfo << "initialization...\n";

    navStokes_steady.initialize();
    navStokes_decoupled1.initialize();
    navStokes_decoupled2.initialize();
    navStokes_decoupled1Iter.initialize();
    navStokes_decoupled1Coup.initialize();
    navStokes_decoupled2Iter.initialize();
    navStokes_decoupled2Coup.initialize();
    navStokes_unsteady.initialize();
 
    navStokes_steady.solve(10,1e-5);
    navStokes_decoupled1.solve(10, 1e-5);
    navStokes_decoupled2.solve(10, 1e-5);
    navStokes_decoupled1Iter.solve(10, 1e-5);
    navStokes_decoupled1Coup.solve(10, 1e-5);
    navStokes_decoupled2Iter.solve(10, 1e-5);
    navStokes_decoupled2Coup.solve(10, 1e-5);
    navStokes_unsteady.solve(10, 1e-5);

    delete Uin;
    delete Uwall;
    delete patches;

    if (omega)
        delete Ublade;

    return 0; 
}

std::string omegaXr_component(real_t omega, std::string variable)
{
    std::ostringstream s;

    s << omega << " * " << variable;

    return s.str();
}

gsBoundaryConditions<> defineBCs(gsFunctionExpr<> & Uin, gsFunctionExpr<> & Uwall, gsFunctionExpr<> & Ublade, char method)
{
    gsBoundaryConditions<> bcInfo;

    
    bcInfo.addCondition(0, boundary::front, condition_type::dirichlet, Uin, 0);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uwall, 0);

    bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(1, boundary::east, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, Ublade, 0);

    bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, Ublade, 0);

    if (method == 'P')
    {
        gsFunctionExpr<real_t> * P = new gsFunctionExpr<real_t>("0", 3);
        bcInfo.addCondition(2, boundary::back, condition_type::dirichlet, P, 1);
    }

    bcInfo.addPeriodic(0, boundary::west, 0, boundary::east, Uin.domainDim());
    bcInfo.addPeriodic(2, boundary::west, 2, boundary::east, Uin.domainDim());
    
    // define transform matrix for periodic sides
    gsMatrix<real_t> transformMatrix(3, 3);

    const real_t phi = -(1. / 7)*EIGEN_PI;
    const real_t cos = math::cos(phi);
    const real_t sin = math::sin(phi);
    transformMatrix(0, 0) = 1;
    transformMatrix(0, 1) = 0;
    transformMatrix(0, 2) = 0;
    transformMatrix(1, 0) = 0;
    transformMatrix(1, 1) = cos;
    transformMatrix(1, 2) = sin;
    transformMatrix(2, 0) = 0;
    transformMatrix(2, 1) = -sin;
    transformMatrix(2, 2) = cos;

    bcInfo.setTransformMatrix(transformMatrix);

    return bcInfo;
}

gsBoundaryConditions<> defineBCs_nonPer(gsFunctionExpr<> & Uin, gsFunctionExpr<> & Uwall, gsFunctionExpr<> & Ublade, char method)
{
    gsBoundaryConditions<> bcInfo;

    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(0, boundary::front, condition_type::dirichlet, Uin, 0);
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, Ublade, 0);

    bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(1, boundary::east, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, Ublade, 0);

    bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(2, boundary::east, condition_type::dirichlet, Ublade, 0);

    if (method == 'P')
    {
        gsFunctionExpr<real_t> * P = new gsFunctionExpr<real_t>("0", "0", "0", 3);
        bcInfo.addCondition(2, boundary::back, condition_type::dirichlet, P, 1);
    }

    return bcInfo;
}
