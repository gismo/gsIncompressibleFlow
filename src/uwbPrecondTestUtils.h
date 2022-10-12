/** @file uwbPrecondTestUtils.h

Contains useful methods for preconditioners testing.

Author(s): H. Hornikova
*/

#include <gismo.h>

using namespace gismo;

int readMatrices(std::string inputNS, std::string inputM, std::string inputMp, gsSparseMatrix<real_t>& matNS, gsSparseMatrix<real_t>& velM, gsSparseMatrix<real_t>& presM, gsMatrix<real_t>& rhs)
{
    matNS.clear();
    velM.clear();
    presM.clear();

    gsInfo << "Reading file " << inputNS << "\n";
    gsReadFile<>(inputNS, matNS);
    gsReadFile<>(inputNS, rhs);

    gsInfo << "Reading file " << inputM << "\n";
    gsReadFile<>(inputM, velM);

    gsInfo << "Reading file " << inputMp << "\n";
    gsReadFile<>(inputMp, presM);

    gsInfo << "Matrix size: " << matNS.rows() << " x " << matNS.cols() << "\n";
    gsInfo << "Non-zeros: " << matNS.nonZeros() << "\n";

    return 0;
}

int readMatrices(std::string inputNS, std::string inputM, std::string inputMp, std::string inputAp, std::string inputFp, gsSparseMatrix<real_t>& matNS, gsSparseMatrix<real_t>& velM, gsSparseMatrix<real_t>& presM, gsSparseMatrix<real_t>& matAp, gsSparseMatrix<real_t>& matFp, gsMatrix<real_t>& rhs)
{
    readMatrices(inputNS, inputM, inputMp, matNS, velM, presM, rhs);

    matAp.clear();
    matFp.clear();

    gsInfo << "Reading file " << inputAp << "\n";
    gsReadFile<>(inputAp, matAp);

    gsInfo << "Reading file " << inputFp << "\n";
    gsReadFile<>(inputFp, matFp);

    return 0;
}

int readMatrices(std::string inputNS, std::string inputM, std::string inputMp, std::vector<std::string> inputAp, std::vector<std::string> inputFp, gsSparseMatrix<real_t>& matNS, gsSparseMatrix<real_t>& velM, gsSparseMatrix<real_t>& presM, std::vector<gsSparseMatrix<real_t> >& matAp, std::vector<gsSparseMatrix<real_t> >& matFp, gsMatrix<real_t>& rhs)
{
    readMatrices(inputNS, inputM, inputMp, matNS, velM, presM, rhs);

    matAp.clear();
    matFp.clear();
    matAp.resize(inputAp.size());
    matFp.resize(inputFp.size());

    for (index_t i = 0; i < inputAp.size(); i++)
    {
        gsInfo << "Reading file " << inputAp[i] << "\n";
        gsReadFile<>(inputAp[i], matAp[i]);

        gsInfo << "Reading file " << inputFp[i] << "\n";
        gsReadFile<>(inputFp[i], matFp[i]);
    }

    return 0;
}

real_t report(const gsVector<>& computedSolution, const gsVector<>& exactSolution)
{
    real_t error = (computedSolution - exactSolution).norm() / exactSolution.norm();

    gsInfo << "Solution error: " << error << "\n\n";

    return error;
}

gsVector<> reportSeparate(const gsVector<>& computedSolution, const gsVector<>& exactSolution, int uSize)
{
    int pdofs = computedSolution.rows() - uSize;
    real_t errU = (computedSolution.topRows(uSize) - exactSolution.topRows(uSize)).norm() / exactSolution.topRows(uSize).norm();
    real_t errP = (computedSolution.bottomRows(pdofs) - exactSolution.bottomRows(pdofs)).norm() / exactSolution.bottomRows(pdofs).norm();

    gsInfo << "Velocity error: " << errU << "\n";
    gsInfo << "Pressure error: " << errP << "\n\n";

    gsVector<> errors(2);
    errors << errU, errP;

    return errors;
}

//void eraseSubStr(std::string & mainStr, const std::string & toErase)
//{
//    size_t pos = mainStr.find(toErase);
//
//    if (pos != std::string::npos)
//    {
//        mainStr.erase(pos, toErase.length());
//    }
//}

//void matFileNamesFromiFileName(std::string ifilePath, std::string& matNS, std::string& matMv, std::string& matMp)
//{
//    size_t pathEnd = ifilePath.find_last_of('/');
//    std::string ifileName = ifilePath.substr(pathEnd + 1);
//    ifileName = ifileName.substr(0, ifileName.length() - 4);
//
//    matNS = ifileName + "_matNS.xml";
//
//    matMv = ifileName + "_matVelMass.xml";
//    // erase what shouldn't be there (viscosity, omega, st/unst)
//
//    matMp = ifileName + "_matPresMass.xml";
//    // erase what shouldn't be there (viscosity, omega, st/unst)
//}