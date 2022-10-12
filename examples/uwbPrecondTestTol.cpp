/** @file uwbPrecondTestTol.cpp

Solves a given linear system with a chosen iterative method and preconditioner(s).
Stops when a given stopping tolerance or the maximum number of iterations is reached. Prints the number of iterations into an output file.

Author(s): H. Hornikova
*/

#include <gismo.h>
#include "uwbLinSolvers.h"
#include "uwbPreconditioners.h"
#include "uwbPrecondTestUtils.h"

using namespace gismo;

template <class SolverType> void solveIterative(const gsSparseMatrix<>& mat, uwbINSPreconditioner<real_t>::Ptr prec, gsMatrix<>& rhs, gsMatrix<>& x, int maxIter, real_t tol, int& numIter, real_t& solTime, bool resHist, gsMatrix<>& res);

int main(int argc, char** argv)
{
    std::string ifilePath; // path to the config file

    gsCmdLine cmd("Preconditioners test.");
    cmd.addString("f", "filePath", "Path to the configuration file", ifilePath);
    try { cmd.getValues(argc, argv); }
    catch (int rv) { return rv; }

    if (ifilePath == "")
    {
        gsInfo << "No configuration file given!\n";
        return -1;
    }

    gsOptionList opt;
    gsReadFile<>(ifilePath, opt);

    // configuration params
    std::string ofilePath, matPath, matNS, matMv, matMp, matAp, matFp, solverName;
    int dim, maxIter;
    real_t viscosity, tol, gamma, gammaM;
    bool resHist;
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
    viscosity = opt.getReal("visc");
    tol = opt.getReal("tol");
    gamma = opt.getReal("prec.gamma");
    gammaM = opt.askReal("prec.gammaM",gamma);
    resHist = opt.askSwitch("resHist", false);
    precVector = opt.getMultiString("precN");

    matNS = matPath + matNS;
    matMv = matPath + matMv;
    matMp = matPath + matMp;
    matAp = matPath + matAp;
    matFp = matPath + matFp;

    gsSparseMatrix<>  mat, velM, presM, presA, presF;
    std::map<std::string, gsSparseMatrix<> > matrices;
    gsMatrix<>        rhs, x, xExact;
    uwbINSPreconditioner<real_t>::Ptr prec;
    gsStopwatch clock;
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
    file << "Max iterarions: " << maxIter << "\n";
    file << "Tolerance: " << tol << "\n";
    file << "Gamma for AL: " << gamma << "\n";
    file << "Gamma for MAL: " << gammaM << "\n";
    file << "Alpha_p for SIMPLE-type: " << opt.getReal("prec.alphaP") << "\n";
    file << "Subsystems iteratively: " << opt.getSwitch("prec.iter") << "\n";
    file << "Subsystems max. iter.: " << opt.getInt("prec.maxIt") << "\n";
    file << "Subsystems tolerance: " << opt.getReal("prec.tol") << "\n";
    file << "Lumped diagonal mass matrices: " << opt.getSwitch("prec.lumpingM") << "\n";
    file << "Lumped diagonal blockA: " << opt.getSwitch("prec.lumpingA") << "\n\n";

    x.setZero(mat.rows(), 1);

    // multiply the part corresponding to continuity eq. by -1 to obtain block-symmetric matrix
    index_t pdofs = presM.rows();
    mat.bottomRows(pdofs) *= -1;
    rhs.bottomRows(pdofs) *= -1;

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
    file << "Using solver: " << solverName << "\n";

    std::vector<gsMatrix<> > resVec;

    // loop over preconditioners
    for (size_t i = 0; i < precVector.size(); i++)
    {
        std::string precType = precVector[i];

        gsInfo << "Using preconditioner: " << precType << "\n";
        file << "Using preconditioner: " << precType << "\n";

        matrices.at("matNS") = mat;

        gsMatrix<>* pRhs = &rhs;

        int numDofs = mat.rows();
        gsSparseMatrix<> mat1(numDofs, numDofs);
        gsMatrix<> rhs1(numDofs, 1);

        if (precType  == "AL_Awhole" || precType == "AL_Amod")
        {
            if (precType == "AL_Awhole")
                precOpt.setReal("gamma",gamma);
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

        int numIter = 0;
        real_t solTime = 0;
        gsMatrix<> res;

        if (solverName == "GMRESright")
            solveIterative<uwbGMResRight<real_t> >(matrices.at("matNS"), prec, *pRhs, x, maxIter, tol, numIter, solTime, resHist, res);
        else if(solverName == "GMRES")
            solveIterative<gsGMRes<real_t> >(matrices.at("matNS"), prec, *pRhs, x, maxIter, tol, numIter, solTime, resHist, res);
        else if (solverName == "GCR")
            solveIterative<uwbGCR<real_t> >(matrices.at("matNS"), prec, *pRhs, x, maxIter, tol, numIter, solTime, resHist, res);
        else if (solverName == "BICG")
            solveIterative<uwbBiCGStab<real_t> >(matrices.at("matNS"), prec, *pRhs, x, maxIter, tol, numIter, solTime, resHist, res);
        else
        {
            gsInfo << "Invalid solver name, using GMRESright.\n";
            solveIterative<uwbGMResRight<real_t> >(matrices.at("matNS"), prec, *pRhs, x, maxIter, tol, numIter, solTime, resHist, res);
        }

        resVec.push_back(res);

        gsInfo << "Number of iterations: " << numIter << "\n";
        gsInfo << "Precond construction time = " << precConstrT << "\n";
        gsInfo << "Solve time: " << solTime << "\n";
        gsInfo << "Total time: " << precConstrT+solTime << "\n\n";
        file << "Number of iterations: " << numIter << "\n";
        file << "Precond construction time = " << precConstrT << "\n";
        file << "Solve time: " << solTime << "\n";
        file << "Total time: " << precConstrT + solTime << "\n\n";
    }

    if (resHist)
    {
        file << "Residual records:\n";

        for (size_t i = 0; i < resVec.size(); i++)
        {
            file << "res" << i << " = {";

            for (int j = 0; j < resVec[i].rows() - 1; j++)
                file << resVec[i](j, 0) << ", ";

            file << resVec[i](resVec[i].rows()-1, 0) << "};\n";
        }
    }

    file.close();

    return 0;
}

template <class SolverType> 
void solveIterative(const gsSparseMatrix<>& mat, uwbINSPreconditioner<real_t>::Ptr prec, gsMatrix<>& rhs, gsMatrix<>& x, int maxIter, real_t tol, int& numIter, real_t& solTime, bool resHist, gsMatrix<>& res)
{
    gsStopwatch clock;
    x.setZero();
    SolverType solver(mat, prec);
    solver.setMaxIterations(maxIter);
    solver.setTolerance(tol);

    if (resHist)
    {
        clock.restart();
        solver.solveDetailed(rhs, x, res);
        solTime = clock.stop();
    }
    else
    {
        clock.restart();
        solver.solve(rhs, x);
        solTime = clock.stop();
    }

    numIter =  solver.iterations();
}