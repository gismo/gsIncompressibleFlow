/** @file uwbPrecondTest.cpp

Solves a given linear system with a chosen iterative method and preconditioner(s).
Performs the given number of iterations and tracks the relative residual norm, relative error (compared to a LU solution) and computational time.
These values are written into an output file after every N iterations, where N is a given number.

Author(s): H. Hornikova
*/


#include <gismo.h>
#include "uwbLinSolvers.h"
#include "uwbPreconditioners.h"
#include "uwbPrecondTestUtils.h"

using namespace gismo;

void print(std::string str, std::ofstream& fd, std::vector<real_t> residuals, std::vector<real_t> errors, std::vector<real_t> times, std::vector<gsVector<> > errorsSep);
template <class SolverType> real_t solveIterative(const gsSparseMatrix<>& mat, uwbINSPreconditioner<real_t>::Ptr prec, gsMatrix<>& rhs, gsMatrix<>& x, gsMatrix<>& xExact, int itStep, int maxIter, real_t tol, int uSize, real_t precConstrT, std::ofstream& file, bool firstFiner = true);

int main(int argc, char** argv)
{
    std::string ifilePath; // path to the config file

    gsCmdLine cmd("Preconditioners test.");
    cmd.addString("f", "filePath", "Path to the configuration file (xml containing gsOptionList)", ifilePath);
    try { cmd.getValues(argc, argv); }
    catch (int rv) { return rv; }

    gsOptionList opt;
    gsReadFile<>(ifilePath, opt);

    // configuration params
    std::string ofilePath, matPath, matNS, matMv, matMp, matAp, matFp, solverName;
    int dim, maxIter, itStep;
    real_t viscosity, tol, gamma, gammaM;
    bool solveLU, firstFiner;
    std::vector<std::string> precVector;

    ofilePath = opt.getString("ofPath");
    matPath = opt.getString("matPath");
    matNS = opt.getString("matNS");
    matMv = opt.getString("matMv");
    matMp = opt.getString("matMp");
    matAp = opt.askString("matAp");
    matFp = opt.askString("matFp");
    solverName = opt.getString("solver");
    dim = opt.getInt("dim");
    maxIter = opt.getInt("maxIt");
    itStep = opt.getInt("itStep");
    viscosity = opt.getReal("visc");
    tol = opt.getReal("tol");
    gamma = opt.getReal("prec.gamma");
    gammaM = opt.askReal("prec.gammaM", gamma);
    solveLU = opt.getSwitch("lu");
    precVector = opt.getMultiString("precN");
    firstFiner = true;

    matNS = matPath + matNS;
    matMv = matPath + matMv;
    matMp = matPath + matMp;
    matAp = matPath + matAp;
    matFp = matPath + matFp;

    std::vector<real_t> residuals, errors, times;
    gsSparseMatrix<>  mat, velM, presM, presA, presF;
    std::map<std::string, gsSparseMatrix<> > matrices;
    gsMatrix<>        rhs, x, xExact;
    uwbINSPreconditioner<real_t>::Ptr prec;
    gsStopwatch clock;
    real_t time;
    std::ofstream file;
    file.open(ofilePath);

    bool isPCD = false;
    for (size_t i = 0; i < precVector.size(); i++)
    {
        if (precVector[i].substr(0, 3) == "PCD")
        {
            isPCD = true;
            break;
        }
    }

    if (isPCD)
    {
        if (readMatrices(matNS, matMv, matMp, matAp, matFp, mat, velM, presM, presA, presF, rhs) == -1)
        {
            gsInfo << "Reading matrices failed.\n";
            return -1;
        }
    }
    else
    {
        if (readMatrices(matNS, matMv, matMp, mat, velM, presM, rhs) == -1)
        {
            gsInfo << "Reading matrices failed.\n";
            return -1;
        }
    }

    matrices.insert(std::make_pair("matNS", mat));
    matrices.insert(std::make_pair("matMu", velM));
    matrices.insert(std::make_pair("matMp", presM));

    if (isPCD)
    {
        matrices.insert(std::make_pair("matAp", presA));
        matrices.insert(std::make_pair("matFp", presF));
    }

    file << "Input matrix:\n" << matNS << "\n\n";
    file << "Viscosity = " << viscosity << "\n";
    file << "Matrix size: " << mat.rows() << " x " << mat.cols() << "\n";
    file << "Non-zeros: " << mat.nonZeros() << "\n";
    file << "Max iterations: " << maxIter << "\n";
    file << "Iter step: " << itStep << "\n";
    file << "Tolerance: " << tol << "\n";
    file << "Gamma for AL: " << gamma << "\n";
    file << "Gamma for MAL: " << gammaM << "\n";
    file << "Alpha_p for SIMPLE-type: " << opt.getReal("prec.alphaP") << "\n";
    file << "Subsystems iteratively: " << opt.getSwitch("prec.iter") << "\n";
    file << "Subsystems max. iter.: " << opt.getInt("prec.maxIt") << "\n";
    file << "Subsystems tolerance: " << opt.getReal("prec.tol") << "\n";
    file << "Lumped diagonal mass matrices: " << opt.getSwitch("prec.lumpingM") << "\n";
    file << "Lumped diagonal blockA: " << opt.getSwitch("prec.lumpingA") << "\n";

    if(firstFiner)
        file << "First " << std::min(10, itStep - 1) << " iterations with step 1.\n\n";

    x.setZero(mat.rows(), 1);

    // multiply the part corresponding to continuity eq. by -1 to obtain block-symmetric matrix
    index_t pdofs = presM.rows();
    mat.bottomRows(pdofs) *= -1;
    rhs.bottomRows(pdofs) *= -1;

    if (solveLU)
    {
        gsInfo << "Solve Ax = b with Eigen's LU factorization (taken as exact solution).\n";

#ifdef GISMO_WITH_PARDISO
        gsSparseSolver<>::PardisoLU solverLU;
        uwbINSSolverBase<real_t>::pardisoSetup(solverLU);
#else
        gsSparseSolver<>::LU solverLU;
#endif

        clock.restart();
        solverLU.compute(mat);
        xExact = solverLU.solve(rhs);
        time = clock.stop();
        gsInfo << "Time: " << time << " s" << "\n\n";
        gsInfo << "xExact norm: " << xExact.norm() << "\n\n";
        file << "LU factorization time: " << time << "\nExact sol norm: " << xExact.norm() << "\n";
    }

    // ---------------------------------------- precond ----------------------------------------

    int udofs = velM.rows();

    if ((udofs + pdofs) == mat.rows())
        udofs /= dim;

    gsOptionList precOpt = opt.getGroup("prec");
    precOpt.addInt("dim", "Problem dimension", dim);
    precOpt.addInt("udofs", "Number of velocity dofs", udofs);
    precOpt.addInt("pdofs", "Number of pressure dofs", pdofs);
    precOpt.addReal("visc", "Viscosity", viscosity);

    gsInfo << "Using solver: " << solverName << "\n";
    file << "Using solver: " << solverName << "\n\n";

    // loop over preconditioners
    for (size_t i = 0; i < precVector.size(); i++)
    {
        std::string precType = precVector[i];

        gsInfo << "Using preconditioner: " << precType << "\n";

        matrices.at("matNS") = mat;
        
        gsMatrix<>* pRhs = &rhs;

        int numDofs = mat.rows();
        gsSparseMatrix<> mat1(numDofs, numDofs);
        gsMatrix<> rhs1(numDofs, 1);

        if (precType == "AL_Awhole" || precType == "AL_Amod")
        {
            if (precType == "AL_Awhole")
                precOpt.setReal("gamma", gamma);
            else
                precOpt.setReal("gamma", gammaM);

            gsInfo << "gamma = " << precOpt.getReal("gamma") << "\n";

            uwbINSBlockPrecondAL<real_t>::fillALmodifSystem_into(mat1, rhs1, matrices, rhs, precOpt);
            matrices.at("matNS") = mat1;
            pRhs = &rhs1;
        }

        clock.restart();
        prec = uwbINSPreconditioner<real_t>::make(precType, matrices, precOpt);
        real_t precConstrT = clock.stop();

        gsInfo << "Precond constructor time = " << precConstrT << "\n";
        file << "Precond constructor time = " << precConstrT << "\n";
            
        int uSize = dim * velM.rows();

        if (solverName == "GMRESright")
            solveIterative<uwbGMResRight<real_t> >(matrices.at("matNS"), prec, *pRhs, x, xExact, itStep, maxIter, tol, uSize, precConstrT, file);
        else if(solverName == "GMRES")
            solveIterative<gsGMRes<real_t> >(matrices.at("matNS"), prec, *pRhs, x, xExact, itStep, maxIter, tol, uSize, precConstrT, file);
        else if (solverName == "GCR")
            solveIterative<uwbGCR<real_t> >(matrices.at("matNS"), prec, *pRhs, x, xExact, itStep, maxIter, tol, uSize, precConstrT, file);
        else if (solverName == "BICG")
            solveIterative<uwbBiCGStab<real_t> >(matrices.at("matNS"), prec, *pRhs, x, xExact, itStep, maxIter, tol, uSize, precConstrT, file);
        else
        {
            gsInfo << "Invalid solver name, using GMRESright.\n";
            solveIterative<uwbGMResRight<real_t> >(matrices.at("matNS"), prec, *pRhs, x, xExact, itStep, maxIter, tol, uSize, precConstrT, file);
        }
    }

    file.close();

    return 0;
}


void print(std::string str, std::ofstream& fd, std::vector<real_t> residuals, std::vector<real_t> errors, std::vector<real_t> times, std::vector<gsVector<> > errorsSep)
{
    int size = residuals.size();

    fd << str + "Res = {";
    for (int i = 0; i < size - 1; i++)
        fd << residuals[i] << ", ";
    fd << residuals[size - 1] << "};\n";

    fd << str + "Err = {";
    for (int i = 0; i < size - 1; i++)
        fd << errors[i] << ", ";
    fd << errors[size - 1] << "};\n";

    /*fd << str + "ErrU = {";
    for (int i = 0; i < size - 1; i++)
        fd << errorsSep[i].at(0) << ", ";
    fd << errorsSep[size - 1].at(0) << "};\n";

    fd << str + "ErrP = {";
    for (int i = 0; i < size - 1; i++)
        fd << errorsSep[i].at(1) << ", ";
    fd << errorsSep[size - 1].at(1) << "};\n";*/

    fd << str + "Time = {";
    for (int i = 0; i < size - 1; i++)
        fd << times[i] << ", ";
    fd << times[size - 1] << "};\n";

    fd << "\n";
}

template <class SolverType> 
real_t solveIterative(const gsSparseMatrix<>& mat, uwbINSPreconditioner<real_t>::Ptr prec, gsMatrix<>& rhs, gsMatrix<>& x, gsMatrix<>& xExact, int itStep, int maxIter, real_t tol, int uSize, real_t precConstrT, std::ofstream& file, bool firstFiner)
{
    std::vector<real_t> residuals;
    std::vector<real_t> errors;
    std::vector<real_t> times;

    std::vector<gsVector<> > errorsSep;
    gsVector<> err0(2);
    err0 << 1, 1;

    // zero iteration
    residuals.push_back(1);
    errors.push_back(1);
    times.push_back(precConstrT);
    errorsSep.push_back(err0);

    gsStopwatch clock;
    real_t time;

    if (firstFiner)
    {
        for (int iter = 1; iter <= std::min(std::min(10, itStep - 1), maxIter); iter += 1)
        {
            x.setZero();
            SolverType solver(mat, prec);
            solver.setMaxIterations(iter);
            solver.setTolerance(tol);
            clock.restart();
            solver.solve(rhs, x);
            time = clock.stop();
            gsInfo << "Time: " << time << " s" << "\n";
            gsInfo << solver.detail();

            residuals.push_back(solver.error());
            errors.push_back(report(x, xExact));
            times.push_back(time + precConstrT);
            //errorsSep.push_back(reportSeparate(x, xExact,uSize));
        }
    }

    for (int iter = itStep; iter <= maxIter; iter += itStep)
    {
        x.setZero();
        SolverType solver(mat, prec);
        solver.setMaxIterations(iter);
        solver.setTolerance(tol);
        clock.restart();
        solver.solve(rhs, x);
        time = clock.stop();
        gsInfo << "Time: " << time << " s" << "\n";
        gsInfo << solver.detail();

        residuals.push_back(solver.error());
        errors.push_back(report(x, xExact));
        times.push_back(time + precConstrT);
        //errorsSep.push_back(reportSeparate(x, xExact, uSize));
    }

    print(prec->getName(), file, residuals, errors, times, errorsSep);

    return residuals.back();
}
