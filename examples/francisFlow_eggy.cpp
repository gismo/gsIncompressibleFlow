#include <gismo.h>

//#include "uwbINSSolverSteady.h"
#include "gsINSSolver.h"

using namespace gismo;

int main(int argc, char *argv[])
{
    real_t viscosity = 0.01;
    int numRefine = 0;
    int uRefine = 0;
    int maxIter = 5;
    int nBlades = 15;
    real_t inVelMag = 1;
    real_t inVelAngle = 20; // in degrees

    int plotPts = 50000;
    bool plotMesh = false;
    std::string outPath = "";
    std::string inFile = FLOW_DATA_DIR "francisGeo.xml";

    gsCmdLine cmd(".");
    cmd.addInt("", "maxIt", "maximum Picard iterations", maxIter);
    cmd.addInt("s", "plotSamples", "", plotPts);
    cmd.addInt("r", "refine", "global uniform refine", numRefine);
    cmd.addInt("", "uRefine", "uniform refine in u-direction", uRefine);
    cmd.addInt("b", "nBlades", "number of blades", nBlades);
    cmd.addReal("v", "visc", "viscosity", viscosity);
    cmd.addReal("", "inVelMag", "magnitude of input velocity", inVelMag);
    cmd.addReal("", "inVelAng", "angle of input velocity", inVelAngle);
    cmd.addSwitch("", "mesh", "plot the computational mesh", plotMesh);
    cmd.addString("", "outPath", "path to the output directory", outPath);
    cmd.addString("", "inFile", "full path to the geometry file", inFile);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    real_t inVelAngleRad = - inVelAngle * (EIGEN_PI / 180); // in radians
    
    std::string cosIn = util::to_string(math::cos(inVelAngleRad));
    std::string sinIn = util::to_string(math::sin(inVelAngleRad));
    std::string magStr = "(" + util::to_string(inVelMag) + "/sqrt(x^2+y^2))*";
    std::string xIn = magStr + "(-" + cosIn + "*x + " + sinIn + "*y)";
    std::string yIn = magStr + "(-" + sinIn + "*x - " + cosIn + "*y)";

    gsFunctionExpr<> Uin(xIn, yIn, "0", 3);
    gsFunctionExpr<> Uwall("0", "0", "0", 3);
    gsFunctionExpr<> f("0", "0", "0", 3);

    real_t phi = -(2. / nBlades)*EIGEN_PI;
    real_t cos = math::cos(phi);
    real_t sin = math::sin(phi);
    gsMatrix<real_t> transformMatrix(3, 3);
    transformMatrix(0, 0) = cos;
    transformMatrix(0, 1) = -sin;
    transformMatrix(0, 2) = 0;
    transformMatrix(1, 0) = sin;
    transformMatrix(1, 1) = cos;
    transformMatrix(1, 2) = 0;
    transformMatrix(2, 0) = 0;
    transformMatrix(2, 1) = 0;
    transformMatrix(2, 2) = 1;

    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uin, 0);
    bcInfo.addCondition(0, boundary::front, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::back, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::front, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::back, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::front, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::back, condition_type::dirichlet, Uwall, 0);
    bcInfo.addPeriodic(0, boundary::west, 0, boundary::east, 3);
    bcInfo.addPeriodic(2, boundary::west, 2, boundary::east, 3);
    bcInfo.setTransformMatrix(transformMatrix);

    gsMultiPatch<> mp;
    gsReadFile<> mpFile(inFile, mp);

    gsMultiBasis<> tbasis(mp); 

    for (int i = 0; i < uRefine; ++i)
        tbasis.uniformRefine(1, 1, 0); // refine in u-direction

    for (int i = 0; i < numRefine; ++i)
        tbasis.uniformRefine();

    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(tbasis);
    discreteBases.push_back(tbasis);
    discreteBases[0].degreeElevate(1);

    gsAssemblerOptions opt;
    opt.dirStrategy = dirichlet::elimination;
    opt.dirValues = dirichlet::interpolation;
    opt.intStrategy = iFace::glue;

    gsNavStokesPde<real_t> NSpde(mp, bcInfo, &f, viscosity);
    gsFlowSolverParams<real_t> params(NSpde, discreteBases);
    params.options().setString("assemb.loop", "EbE");

    gsINSSolverSteady<real_t> solver(params);

    gsInfo << "numDofs: " << solver.numDofs() << "\n\n";
    gsInfo << "initialization...\n";

    solver.initialize();

    gsInfo << "done" << "\n";

    solver.solve(maxIter, 1e-4);

    gsInfo << "\nAssembly time:" << solver.getInitAssemblyTime() + solver.getAssemblyTime() << "\n";
    gsInfo << "Solve time:" << solver.getSolveTime() << "\n";
    gsInfo << "Solver setup time:" << solver.getSolverSetupTime() << "\n";


    // std::string viscStr = util::to_string(viscosity);
    // std::replace(viscStr.begin(), viscStr.end(), '.', '-');
    // viscStr.erase(viscStr.find_last_not_of('0') + 1, std::string::npos);

    // std::string velMagStr = util::to_string(inVelMag);
    // std::replace(velMagStr.begin(), velMagStr.end(), '.', '-');
    // velMagStr.erase(velMagStr.find_last_not_of('0') + 1, std::string::npos);

    // std::string idStr = "_ref" + util::to_string(numRefine) + "_uRef" + util::to_string(uRefine) + "_visc" + viscStr + "_velMag" + velMagStr + "_angle" + util::to_string(inVelAngle) + "_";

    gsField<> velocity = solver.constructSolution(0);
    gsField<> pressure = solver.constructSolution(1);
    gsWriteParaview<>(velocity, outPath + "francisVelocity", plotPts, plotMesh);
    gsWriteParaview<>(pressure, outPath + "francisPressure", plotPts, false);

    return 0;
}