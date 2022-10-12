/** @file uwbPrecondOptGamma.cpp

Searches for an optimal parameter gamma for the AL preconditioner from a given interval with a given step.

Author(s): H. Hornikova
*/

#include <gismo.h>
#include "uwbLinSolvers.h"
#include "uwbPreconditioners.h"
#include "uwbPrecondTestUtils.h"

using namespace gismo;

template <class SolverType> int findGammaOpt(const gsSparseMatrix<>& mat, uwbINSPreconditioner<real_t>::Ptr prec, gsMatrix<>& rhs, gsMatrix<>& x, int maxIter, real_t gamma, real_t tol);
void fillALmodifSystem_into(gsSparseMatrix<>& mat1, gsMatrix<>& rhs1, const gsSparseMatrix<>& mat, const gsMatrix<>& rhs, const gsSparseMatrix<>& velM, const gsSparseMatrix<>& presM, int dim, real_t gamma);

int main(int argc, char** argv)
{
    std::string ifilePath; // path to the config file

    gsCmdLine cmd("Optimal gamma computation.");
    cmd.addString("f", "filePath", "Path to the configuration file", ifilePath);
    try { cmd.getValues(argc, argv); }
    catch (int rv) { return rv; }

    gsOptionList opt;
    gsReadFile<>(ifilePath, opt);

    // configuration params
    std::string ofilePath, matPath, matNS, matMv, matMp, solverName;
    int dim, maxIter;
    real_t viscosity, tol, gMin, gMax, gStep;
    std::vector<std::string> precVector;

    ofilePath = opt.getString("ofPath");
    matPath = opt.getString("matPath");
    matNS = opt.getString("matNS");
    matMv = opt.getString("matMv");
    matMp = opt.getString("matMp");
    solverName = opt.getString("solver");
    dim = opt.getInt("dim");
    maxIter = opt.getInt("maxIt");
    viscosity = opt.getReal("visc");
    tol = opt.getReal("tol");
    gMin = opt.getReal("gammaMin");
    gMax = opt.getReal("gammaMax");
    gStep = opt.getReal("gammaStep");
    precVector = opt.getMultiString("precN");

    matNS = matPath + matNS;
    matMv = matPath + matMv;
    matMp = matPath + matMp;

    std::vector<real_t> residuals, errors, times;
    gsSparseMatrix<>  mat, velM, presM;//, presA, presF;
    std::map<std::string, gsSparseMatrix<> > matrices;
    gsMatrix<>        rhs, x, xExact;
    uwbINSPreconditioner<real_t>::Ptr prec;
    gsStopwatch clock;
    std::ofstream file;
    file.open(ofilePath);

    if (readMatrices(matNS, matMv, matMp, mat, velM, presM, rhs) == -1)
    {
        gsInfo << "Reading matrices failed.\n";
        return -1;
    }

    matrices.insert(std::make_pair("matNS", mat));
    matrices.insert(std::make_pair("matMu", velM));
    matrices.insert(std::make_pair("matMp", presM));

    file << "Input matrix:\n" << matNS << "\n\n";
    file << "Viscosity = " << viscosity << "\n";
    file << "Matrix size: " << mat.rows() << " x " << mat.cols() << "\n";
    file << "Non-zeros: " << mat.nonZeros() << "\n";
    file << "Gamma min: " << gMin << "\n";
    file << "Gamma max: " << gMax << "\n";
    file << "Gamma step: " << gStep << "\n";
    file << "Max iterations: " << maxIter << "\n";
    file << "Tolerance: " << tol << "\n";
    file << "Subsystems iteratively: " << opt.getSwitch("prec.iter") << "\n";
    file << "Subsystems max. iter.: " << opt.getInt("prec.maxIt") << "\n";
    file << "Subsystems tolerance: " << opt.getReal("prec.tol") << "\n";
    file << "Lumped diagonal mass matrices: " << opt.getSwitch("prec.lumpingM") << "\n\n";

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
    precOpt.addReal("gamma", "Gamma", gMin);

    gsInfo << "Using solver: " << solverName << "\n";
    file << "Using solver: " << solverName << "\n";

    // loop over preconditioners
    for (size_t i = 0; i < precVector.size(); i++)
    {
        std::string precType = precVector[i];

        gsInfo << "Using preconditioner: " << precType << "\n";

        int numDofs = mat.rows();

        real_t gammaOpt = gMin;
        int minIter = maxIter;

        for (real_t gamma = gMin; gamma <= gMax; gamma += gStep)
        {
            gsInfo << "gamma = " << gamma << "\n";
            file << "gamma = " << gamma << "\n";

            matrices.at("matNS") = mat;
            precOpt.setReal("gamma", gamma);

            gsSparseMatrix<> mat1(numDofs, numDofs);
            gsMatrix<> rhs1(numDofs, 1);

            uwbINSBlockPrecondAL<real_t>::fillALmodifSystem_into(mat1, rhs1, matrices, rhs, precOpt);
            matrices.at("matNS") = mat1;
            
            prec = uwbINSPreconditioner<real_t>::make(precType, matrices, precOpt);

            int iter;
            if (solverName == "GMRESright")
                iter = findGammaOpt<uwbGMResRight<real_t> >(matrices.at("matNS"), prec, rhs1, x, maxIter, gamma, tol);
            else if (solverName == "GMRES")
                iter = findGammaOpt<gsGMRes<real_t> >(matrices.at("matNS"), prec, rhs1, x, maxIter, gamma, tol);
            else if (solverName == "GCR")
                iter = findGammaOpt<uwbGCR<real_t> >(matrices.at("matNS"), prec, rhs1, x, maxIter, gamma, tol);
            else if (solverName == "BICG")
                iter = findGammaOpt<uwbBiCGStab<real_t> >(matrices.at("matNS"), prec, rhs1, x, maxIter, gamma, tol);
            else
            {
                gsInfo << "Invalid solver name, using GMRESright.\n";
                iter = findGammaOpt<uwbGMResRight<real_t> >(matrices.at("matNS"), prec, rhs1, x, maxIter, gamma, tol);
            }

            file << "numIter = " << iter << "\n\n";

            if (iter < minIter)
            {
                gammaOpt = gamma;
                minIter = iter;
            }
        }

        file << "gamma opt: " << gammaOpt << "\n\n";

    }

    file.close();

    return 0;
}

template <class SolverType> 
int findGammaOpt(const gsSparseMatrix<>& mat, uwbINSPreconditioner<real_t>::Ptr prec, gsMatrix<>& rhs, gsMatrix<>& x, int maxIter, real_t gamma, real_t tol)
{
    x.setZero();
    SolverType solver(mat, prec);
    solver.setMaxIterations(maxIter);
    solver.setTolerance(tol);
    solver.solve(rhs, x);

    gsInfo << solver.detail();

    return solver.iterations();
}

void fillALmodifSystem_into(gsSparseMatrix<>& mat1, gsMatrix<>& rhs1, const gsSparseMatrix<>& mat, const gsMatrix<>& rhs, const gsSparseMatrix<>& velM, const gsSparseMatrix<>& presM, int dim, real_t gamma)
{
    gsInfo << "Filling the modified AL matrix... ";

    int udofs = velM.rows();
    int pdofs = presM.rows();
    int uSize = dim * udofs;
    int numDofs = uSize + pdofs;

    gsSparseMatrix<> presMinv(pdofs, pdofs); // approximation of pressure mass matrix inverse
    presMinv.setIdentity();
    for (index_t i = 0; i < pdofs; i++)
        presMinv.coeffRef(i, i) = 1 / presM.coeff(i, i);

    gsVector<int> nonZerosPerColumnVector;
    nonZerosPerColumnVector.setZero(numDofs);

    gsSparseMatrix<> Mgamma = gamma * mat.block(0, uSize, uSize, pdofs) * presMinv * mat.block(uSize, 0, pdofs, uSize);

    for (int i = 0; i < uSize; i++)
        nonZerosPerColumnVector(i) += Mgamma.col(i).nonZeros();

    mat1 = mat;
    mat1.reserve(nonZerosPerColumnVector);

    for (index_t col = 0; col < uSize; ++col)
        for (typename gsSparseMatrix<>::InnerIterator it(Mgamma, col); it; ++it)
            mat1.coeffRef(it.row(), col) += Mgamma(it.row(), col);

    mat1.makeCompressed();

    rhs1 = rhs;
    rhs1.topRows(uSize) += gamma * mat.block(0, uSize, uSize, pdofs) * presMinv * rhs.bottomRows(pdofs);

    gsInfo << "Done.\n";
}
