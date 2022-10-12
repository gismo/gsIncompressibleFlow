/*
    Author(s): J. Egermaier
*/

#ifndef UWBLINEARTESTING_H
#define UWBLINEARTESTING_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <gismo.h>
#include "uwbReadWriteOpt.h"
#include "uwbLinSolvers.h"
#include "uwbPreconditioners.h"
#include "uwbGeometryUtils.h"
#include "uwbGeometryProfile2D.h"
#include "uwbGeometryStep.h"
#include "uwbGenerateMatrices.h"

using namespace gismo;

template<class T>
class testMatrix
{
public:
    testMatrix(std::string inputFile)
    {
        m_inputFile = inputFile;
        paramList = readInputFile(m_inputFile);
        get_parameter(paramList,mainDIR,"mainDir");
        get_parameter(paramList,matrixDIR,"matrixDir"); matrixDIR = mainDIR + matrixDIR;
        get_parameter(paramList,domainDIR,"domainDir"); domainDIR = mainDIR + domainDIR;
    }

    testMatrix(std::map<std::string, gsVector<std::string>> pList)
    {
        //m_inputFile = inputFile;
        paramList = pList;
        get_parameter(paramList,mainDIR,"mainDir");
        get_parameter(paramList,matrixDIR,"matrixDir"); matrixDIR = mainDIR + matrixDIR;
        get_parameter(paramList,domainDIR,"domainDir"); domainDIR = mainDIR + domainDIR;
    }

public:
    gsMultiPatch<T> domain;
    gsMultiBasis<T> tbasis;
    std::string nameOfDomain, m_inputFile, mainDIR, matrixDIR, domainDIR, domainFile, tbasisFile, patchesAndBasisFile;
    std::map<std::string, gsVector<std::string>> paramList;


    void set_parameters()
    {
        get_parameter(paramList,nameOfDomain,"geometry_type");
    }

    void save_patches_and_basis()
    {
        if (nameOfDomain == "profile2D"){
            geometryProfile2D(paramList,domainDIR);
        }
        if (nameOfDomain == "runner3D"){
        //    geometryrunner3D(paramList);
        }
        if (nameOfDomain == "step2D"){
            geometryStep<T>(paramList,domainDIR);
        }
        if (nameOfDomain == "step3D"){
            geometryStep<T>(paramList,domainDIR);
        }
    }

    void load_patches_and_basis()
    {
        get_parameter(paramList,patchesAndBasisFile,"patchesAndBasisFile");
        domainFile = domainDIR + patchesAndBasisFile + "_patches.xml";
        tbasisFile = domainDIR + patchesAndBasisFile + "_tbasis.xml";
        if (nameOfDomain == "profile2D"){
            gsReadFile<T> (domainFile,domain);
            gsReadFile<T> (tbasisFile,tbasis);
        }
        if (nameOfDomain == "step2D"){
            gsReadFile<T> (domainFile,domain);
            gsReadFile<T> (tbasisFile,tbasis);
        }
        if (nameOfDomain == "step3D"){
            gsReadFile<T> (domainFile,domain);
            gsReadFile<T> (tbasisFile,tbasis);
        }
        if (nameOfDomain == "runner3D"){
            //domain = runnerBladeDomain(bool linearize, int numRefine = 1, bool plot = false, bool print_info = false);
        }
    }

    void generateNSmatrices()
    {
        if (nameOfDomain == "profile2D"){
            std::string paraviewDIR = matrixDIR + "paraview/";
            generateMatrices(paramList, domain, tbasis, paraviewDIR);
        }
        if (nameOfDomain == "step2D"){
            std::string paraviewDIR = matrixDIR + "paraview/";
            generateMatrices(paramList, domain, tbasis, paraviewDIR);
        }
        if (nameOfDomain == "step3D"){
            std::string paraviewDIR = matrixDIR + "paraview/";
            generateMatrices(paramList, domain, tbasis, paraviewDIR);
        }
    }

};

template<class T>
class linearTesting
{

public:
    linearTesting(std::string inputFile)
    {
        m_inputFile = inputFile;
        paramList = readInputFile(m_inputFile);
        get_parameter(paramList,mainDIR,"mainDir");
        get_parameter(paramList,matrixDIR,"matrixDir"); matrixDIR = mainDIR + matrixDIR;
        get_parameter(paramList,precondTestTolDIR,"precondTestTolDir"); precondTestTolDIR = mainDIR + precondTestTolDIR;
        //get_parameter(paramList,matrixFileNS,"matNS"); matrixFileNS = matrixDIR + matrixFileNS;
    }

public:
    gsSparseMatrix<T> m_NSmatrix, m_Vmatrix, m_Pmatrix, m_Apmatrix, m_Fpmatrix;
    gsMatrix<T> m_b;
    int dim, pdofs, udofs, maxIter;
    real_t viscosity, tol;
    std::map<std::string, gsSparseMatrix<T> > matrices;
    std::map<std::string, gsVector<std::string>> paramList;
    std::string m_inputFile, matrixFileNS, matrixFileV, matrixFileP, matrixFileAp, matrixFileFp, outputDIR, mainDIR, matrixDIR, precondTestTolDIR, solverName;
    gsVector<std::string> preconditioners, tests, precVector;
    gsOptionList precOpt;


public:

    void set_matrixFileNS(std::string fileName)
    {
        matrixFileNS = matrixDIR + fileName;
    }

    real_t readValue(std::string UK)
    {
        real_t out;
        std::string value = "";
        int pos = matrixFileNS.find(UK);
        if (pos!=std::string::npos){
            pos += UK.size();
            while (matrixFileNS[pos] != '_'){
                value += matrixFileNS[pos];
                pos++;
            }
            util::string_replace(value,"-",".");
            out = std::atof(value.c_str());
        } else {
            if (UK == "visc_"){
            get_parameter(paramList,out,"visc");
            }
        }
        return out;
    }

    void delString(std::string & name, std::string delText, int delNext)
    {
        int pos = name.find(delText);
        if (pos!=std::string::npos){
            name.erase(pos,delText.size());
            for (int i = 0; i < delNext; i++){
                while (name[pos] != '_'){
                    name.erase(pos,1);
                }
                name.erase(pos,1);
            }
        }
    }

    std::string getFileName(std::string endName)
    {
        std::string now = matrixFileNS;
        //delString(now,"visc_",1); // 1 - maze i dalsi text do dalsiho "_". Tj. zde hodnotu viskozity
        //delString(now,"om_",1);
        //delString(now,"st_",0);

        util::string_replace(now,"matNS",endName);
        return now;
    }

    void set_dim()
    {
        if (util::starts_with(gsFileManager::getBasename(matrixFileNS),"profile2D")){
            dim = 2;
        } else if (util::starts_with(gsFileManager::getBasename(matrixFileNS),"runner")){
            dim = 3;
        } else if (util::starts_with(gsFileManager::getBasename(matrixFileNS),"step2D")){
            dim = 2;
        } else if (util::starts_with(gsFileManager::getBasename(matrixFileNS),"step3D")){
            dim = 3;
        } else if (util::starts_with(gsFileManager::getBasename(matrixFileNS),"cavity")){
            dim = 2;
        } else {
            get_parameter(paramList,dim,"dim");
        }
    }

    void set_viscosity()
    {
        viscosity = readValue("visc_");
    }

    void set_matrixFromXML(std::string matFile)
    {
        if (util::ends_with(matFile,"matNS.xml")){
            matrixFileNS = matFile;
            gsReadFile<T> (matrixFileNS, m_NSmatrix);
            set_dim();
            set_viscosity();
        } else if (util::ends_with(matFile,"matVelMass.xml")){
            matrixFileV = matFile;
            gsReadFile<T> (matrixFileV, m_Vmatrix);
        } else if (util::ends_with(matFile,"matPresMass.xml")){
            matrixFileP = matFile;
            gsReadFile<T> (matrixFileP, m_Pmatrix);
        } else {
            gsInfo << "No matrix file was set, wrong file name! (possible: ...matNS.xml, ...matVelMass.xml, ...matPresMass.xml). Press any key..." << "\n";
            while (getchar() != '\n');
        }
    }

    void set_matricesFromXML()
    {
        gsReadFile<T> (matrixFileNS, m_NSmatrix);
        gsReadFile<T> (matrixFileNS, m_b);

        outputDIR = gsFileManager::getBasename(matrixFileNS); outputDIR = mainDIR + "inputs/" + outputDIR + "/";
        gsFileManager::mkdir(mainDIR + "inputs/");
        gsFileManager::mkdir(outputDIR);

        matrixFileV = getFileName("matVelMass"); gsInfo << "Velocity file: " << matrixFileV << "\n";
        matrixFileP = getFileName("matPresMass"); gsInfo << "Presure file: " << matrixFileP << "\n";
        matrixFileAp = getFileName("matAp"); gsInfo << "Ap file: " << matrixFileAp << "\n";
        matrixFileFp = getFileName("matFp"); gsInfo << "Fp file: " << matrixFileFp << "\n";
        gsReadFile<T> (matrixFileV, m_Vmatrix);
        gsReadFile<T> (matrixFileP, m_Pmatrix);
        gsReadFile<T> (matrixFileAp, m_Apmatrix);
        gsReadFile<T> (matrixFileFp, m_Fpmatrix);

        std::ifstream inFile(m_inputFile);
        std::ofstream cpFile(outputDIR + "/" + gsFileManager::getBasename(m_inputFile));
        cpFile << inFile.rdbuf();

        set_dim();
        set_viscosity();
    }

    void set_NSmatrix(gsSparseMatrix<> mat)
    {
        m_NSmatrix = &mat;
        matrixFileNS = "unsaved";
        outputDIR = mainDIR + "unsaved";
        set_dim();
        set_viscosity();
    }

    void set_Vmatrix(gsSparseMatrix<> mat)
    {
        m_Vmatrix = &mat;
    }

    void set_Pmatrix(gsSparseMatrix<> mat)
    {
        m_Pmatrix = &mat;
    }

    void set_aPmatrix(gsSparseMatrix<> mat)
    {
        m_Apmatrix = &mat;
    }

    void set_fPmatrix(gsSparseMatrix<> mat)
    {
        m_Fpmatrix = &mat;
    }

    void set_rhs(gsMatrix<> mat)
    {
        m_b = mat;
    }

    void set_parameters()
    {
        int PprecmaxIt, Pprecfill;
        real_t Pprectol, PprecdropTol, Pprecgamma, PprecgammaM, PprecalphaP;
        bool Ppreciter, PpreclumpingM, PpreclumpingA, saveSub;

        get_parameter(paramList,tests,"tests");
        get_parameter(paramList,solverName,"solver");
        get_parameter(paramList,maxIter,"PmaxIt");
        get_parameter(paramList,tol,"Ptol");
        get_parameter(paramList,precVector,"precN");

        pdofs = m_Pmatrix.rows();
        udofs = m_Vmatrix.rows();
        if ((udofs + pdofs) == m_NSmatrix.rows())
            udofs /= dim;

        std::string mP = matrixFileNS;
        util::string_replace(mP,"matNS.xml","");
        precOpt.addString("matrixPath","path and name of matrices to save", mP);
        precOpt.addInt("dim", "Problem dimension", dim);
        precOpt.addInt("udofs", "Number of velocity dofs", udofs);
        precOpt.addInt("pdofs", "Number of pressure dofs", pdofs);
        precOpt.addReal("visc", "Viscosity", viscosity);
        get_parameter(paramList,saveSub,"PprecSaveSubsystems"); precOpt.addSwitch("save_subsystems","",saveSub);
        get_parameter(paramList,PprecmaxIt,"PprecmaxIt");       precOpt.addInt("maxIt","",PprecmaxIt);
        get_parameter(paramList,Pprecfill,"Pprecfill");         precOpt.addInt("fill","",Pprecfill);
        get_parameter(paramList,Pprectol,"Pprectol");           precOpt.addReal("tol","",Pprectol);
        get_parameter(paramList,PprecdropTol,"PprecdropTol");   precOpt.addReal("dropTol","",PprecdropTol);
        get_parameter(paramList,Pprecgamma,"Pprecgamma");       precOpt.addReal("gamma","",Pprecgamma);
        get_parameter(paramList,PprecgammaM,"PprecgammaM");     precOpt.addReal("gammaM","",PprecgammaM);
        get_parameter(paramList,PprecalphaP,"PprecalphaP");     precOpt.addReal("alphaP","",PprecalphaP);
        get_parameter(paramList,Ppreciter,"Ppreciter");         precOpt.addSwitch("iter","",Ppreciter);
        get_parameter(paramList,PpreclumpingM,"PpreclumpingM"); precOpt.addSwitch("lumpingM","",PpreclumpingM);
        get_parameter(paramList,PpreclumpingA,"PpreclumpingA"); precOpt.addSwitch("lumpingA","",PpreclumpingA);

        bool isPCD = false;
        for (unsigned i = 0; i < precVector.size(); i++){
            if (precVector[i].substr(0, 3) == "PCD"){
                isPCD = true;
                break;
            }
        }

        matrices.insert(std::make_pair("matNS", m_NSmatrix));
        matrices.insert(std::make_pair("matMu", m_Vmatrix));
        matrices.insert(std::make_pair("matMp", m_Pmatrix));

        if (isPCD){
            matrices.insert(std::make_pair("matAp", m_Apmatrix));
            matrices.insert(std::make_pair("matFp", m_Fpmatrix));
        }

        /*get_parameter(paramList,solver,"solver");
        //============== Gamma =================================
        get_parameter(paramList,OmaxIt,"OmaxIt");
        get_parameter(paramList,gammaMax,"gammaMax");
        get_parameter(paramList,gammaMin,"gammaMin");
        get_parameter(paramList,gammaStep,"gammaStep");
        get_parameter(paramList,Otol,"Otol");
        get_parameter(paramList,OprecN,"OprecN");
        get_parameter(paramList,OprecmaxIt,"OprecmaxIt");
        get_parameter(paramList,Oprecfill,"Oprecfill");
        get_parameter(paramList,Oprectol,"Oprectol");
        get_parameter(paramList,OprecdropTol,"OprecdropTol");
        get_parameter(paramList,Opreciter,"Opreciter");
        get_parameter(paramList,OpreclumpingM,"OpreclumpingM");
        */
    }

    void make_symmetry()
    {
        // multiply the part corresponding to continuity eq. by -1 to obtain block-symmetric matrix
        m_NSmatrix.bottomRows(pdofs) *= -1;
        m_b.bottomRows(pdofs) *= -1;
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
        gsInfo << "Pocet iterací je: " << solver << "\n";
    }


    template <class SolverType>
    real_t solveIterative(const gsSparseMatrix<>& mat, uwbINSPreconditioner<real_t>::Ptr prec, gsMatrix<>& rhs, gsMatrix<>& x, gsMatrix<>& xExact, int itStep, int uSize, real_t precConstrT, std::ofstream& file, bool firstFiner = true)
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

        //for (int iter = itStep; iter <= maxIter; iter += itStep)
        real_t solverError = tol + 1.;
        int iter = itStep;
        while ((iter <= maxIter) && (solverError > tol))
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
            solverError = solver.error();
            iter += itStep;
        }

        print(prec->getName(), file, residuals, errors, times, errorsSep);

        return residuals.back();
    }

    T report(const gsVector<T>& computedSolution, const gsVector<T>& exactSolution)
    {
        T error = (computedSolution - exactSolution).norm() / exactSolution.norm();
        gsInfo << "Solution error: " << error << "\n\n";
        return error;
    }

    gsVector<T> reportSeparate(const gsVector<T>& computedSolution, const gsVector<T>& exactSolution, int uSize)
    {
        int pdofs = computedSolution.rows() - uSize;
        T errU = (computedSolution.topRows(uSize) - exactSolution.topRows(uSize)).norm() / exactSolution.topRows(uSize).norm();
        T errP = (computedSolution.bottomRows(pdofs) - exactSolution.bottomRows(pdofs)).norm() / exactSolution.bottomRows(pdofs).norm();

        gsInfo << "Velocity error: " << errU << "\n";
        gsInfo << "Pressure error: " << errP << "\n\n";

        gsVector<T> errors(2);
        errors << errU, errP;

        return errors;
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

        fd << str + "Time = {";
        for (int i = 0; i < size - 1; i++)
            fd << times[i] << ", ";
        fd << times[size - 1] << "};\n";

        fd << "\n";
    }

    void gammaOpt(){}

    void precondTest()
    {
        std::string ofilePath;
        int itStep;
        bool solveLU, firstFiner;
        gsMatrix<> x, xExact;
        uwbINSPreconditioner<real_t>::Ptr prec;
        gsStopwatch clock;
        real_t time;

        get_parameter(paramList,itStep,"itStep");
        real_t Gamma = precOpt.getReal("gamma");
        real_t GammaM = precOpt.getReal("gammaM");
        get_parameter(paramList,solveLU,"lu");
        get_parameter(paramList,firstFiner,"firstFiner");

        ofilePath = outputDIR + "/summaryPrecondTest.txt";
        std::ofstream file;
        file.open(ofilePath);
        file << "Input matrix:\n" << matrixFileNS << "\n\n";
        file << "Viscosity = " << viscosity << "\n";
        file << "Matrix size: " << m_NSmatrix.rows() << " x " << m_NSmatrix.cols() << "\n";
        file << "Non-zeros: " << m_NSmatrix.nonZeros() << "\n";
        file << "Max iterarions: " << maxIter << "\n";
        file << "Tolerance: " << tol << "\n";
        file << "Gamma for AL: " << Gamma << "\n";
        file << "Gamma for MAL: " << GammaM << "\n";
        file << "Alpha_p for SIMPLE-type: " << precOpt.getReal("alphaP") << "\n";
        file << "Subsystems iteratively: " << precOpt.getSwitch("iter") << "\n";
        file << "Subsystems max. iter.: " << precOpt.getInt("maxIt") << "\n";
        file << "Subsystems tolerance: " << precOpt.getReal("tol") << "\n";
        file << "Lumped diagonal mass matrices: " << precOpt.getSwitch("lumpingM") << "\n";
        file << "Lumped diagonal blockA: " << precOpt.getSwitch("lumpingA") << "\n\n";
        file << "Using solver: " << solverName << "\n";

        if(firstFiner)
            file << "First " << std::min(10, itStep - 1) << " iterations with step 1.\n\n";

        x.setZero(m_NSmatrix.rows(), 1);

        if (solveLU){
            gsInfo << "Solve Ax = b with Eigen's LU factorization (taken as exact solution).\n";
            gsSparseSolver<>::LU solverLU;
            clock.restart();
            solverLU.compute(m_NSmatrix);
            xExact = solverLU.solve(m_b);
            time = clock.stop();
            gsInfo << "Time: " << time << " s" << "\n\n";
            gsInfo << "xExact norm: " << xExact.norm() << "\n\n";
            file << "LU factorization time: " << time << "\nExact sol norm: " << xExact.norm() << "\n";
        }

        // ---------------------------------------- precond ----------------------------------------

        gsInfo << "Using solver: " << solverName << "\n";

        // loop over preconditioners
        for (unsigned i = 0; i < precVector.size(); i++){
            std::string precType = precVector[i];
            gsInfo << "Using preconditioner: " << precType << "\n";
            matrices.at("matNS") = m_NSmatrix;
            gsMatrix<>* pRhs = &m_b;

            int numDofs = m_NSmatrix.rows();
            gsSparseMatrix<> mat1(numDofs, numDofs);
            gsMatrix<> rhs1(numDofs, 1);

            if (precType == "AL_Awhole" || precType == "AL_Amod"){
                if (precType == "AL_Awhole")
                    precOpt.setReal("gamma", Gamma);
                else
                    precOpt.setReal("gamma", GammaM);

                gsInfo << "gamma = " << precOpt.getReal("gamma") << "\n";

                uwbINSBlockPrecondAL<real_t>::fillALmodifSystem_into(mat1, rhs1, matrices, m_b, precOpt);
                matrices.at("matNS") = mat1;
                pRhs = &rhs1;
            }

            clock.restart();
            prec = uwbINSPreconditioner<real_t>::make(precType, matrices, precOpt);
            real_t precConstrT = clock.stop();

            gsInfo << "Precond constructor time = " << precConstrT << "\n";
            file << "Precond constructor time = " << precConstrT << "\n";

            int uSize = dim * m_Vmatrix.rows();

            if (solverName == "GMRESright")
                solveIterative<uwbGMResRight<real_t> >(matrices.at("matNS"), prec, *pRhs, x, xExact, itStep, uSize, precConstrT, file, firstFiner);
            else if(solverName == "GMRES")
                solveIterative<gsGMRes<real_t> >(matrices.at("matNS"), prec, *pRhs, x, xExact, itStep, uSize, precConstrT, file, firstFiner);
            else if (solverName == "GCR")
                solveIterative<uwbGCR<real_t> >(matrices.at("matNS"), prec, *pRhs, x, xExact, itStep, uSize, precConstrT, file, firstFiner);
            else if (solverName == "BICG")
                    solveIterative<uwbBiCGStab<real_t> >(matrices.at("matNS"), prec, *pRhs, x, xExact, itStep, uSize, precConstrT, file, firstFiner);
            else {
                gsInfo << "Invalid solver name, using GMRESright.\n";
                solveIterative<uwbGMResRight<real_t> >(matrices.at("matNS"), prec, *pRhs, x, xExact, itStep, uSize, precConstrT, file, firstFiner);
            }
        }
        file.close();
    }

    void precondTestTol()
    {
        std::string ofilePath;
        gsMatrix<T> x;
        uwbINSPreconditioner<real_t>::Ptr prec;
        gsStopwatch clock;
        bool resHist;
        get_parameter(paramList,resHist,"resHist");

        x.setZero(m_NSmatrix.rows(), 1);

        // ---------------------------------------- precond ----------------------------------------

        real_t Gamma = precOpt.getReal("gamma");
        real_t GammaM = precOpt.getReal("gammaM");

        ofilePath = precondTestTolDIR + "TT_" + gsFileManager::getBasename(matrixFileNS) + ".txt";
        std::ofstream file;
        file.open(ofilePath);
        file << "Input matrix:\n" << matrixFileNS << "\n\n";
        file << "Viscosity = " << viscosity << "\n";
        file << "Matrix size: " << m_NSmatrix.rows() << " x " << m_NSmatrix.cols() << "\n";
        file << "Non-zeros: " << m_NSmatrix.nonZeros() << "\n";
        file << "Max iterarions: " << maxIter << "\n";
        file << "Tolerance: " << tol << "\n";
        file << "Gamma for AL: " << Gamma << "\n";
        file << "Gamma for MAL: " << GammaM << "\n";
        file << "Alpha_p for SIMPLE-type: " << precOpt.getReal("alphaP") << "\n";
        file << "Subsystems iteratively: " << precOpt.getSwitch("iter") << "\n";
        file << "Subsystems max. iter.: " << precOpt.getInt("maxIt") << "\n";
        file << "Subsystems tolerance: " << precOpt.getReal("tol") << "\n";
        file << "Lumped diagonal mass matrices: " << precOpt.getSwitch("lumpingM") << "\n";
        file << "Lumped diagonal blockA: " << precOpt.getSwitch("lumpingA") << "\n\n";
        file << "Using solver: " << solverName << "\n";

        gsInfo << "Using solver: " << solverName << "\n";
        std::vector<gsMatrix<> > resVec;

        for (unsigned i = 0; i < precVector.size(); i++)
        {
            std::string precType = precVector[i];

            gsInfo << "Using preconditioner: " << precType << "\n";
            file << "Using preconditioner: " << precType << "\n";

            matrices.at("matNS") = m_NSmatrix;

            gsMatrix<>* pRhs = &m_b;

            int numDofs = m_NSmatrix.rows();
            gsSparseMatrix<> mat1(numDofs, numDofs);
            gsMatrix<> rhs1(numDofs, 1);

            if (precType  == "AL_Awhole" || precType == "AL_Amod")
            {
                if (precType == "AL_Awhole")
                    precOpt.setReal("gamma",Gamma);
                else
                    precOpt.setReal("gamma",GammaM);

                gsInfo << "gamma = " << precOpt.getReal("gamma") << "\n";

                uwbINSBlockPrecondAL<real_t>::fillALmodifSystem_into(mat1, rhs1, matrices, m_b, precOpt);
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
            for (unsigned i = 0; i < resVec.size(); i++)
            {
                file << "res" << i << " = {";
                for (int j = 0; j < resVec[i].rows() - 1; j++)
                    file << resVec[i](j, 0) << ", ";
                file << resVec[i](resVec[i].rows()-1, 0) << "};\n";
            }
        }

        file.close();
    }

    void makeTests()
    {
        gsInfo << "Making tests" << "\n";
        for (int i = 0; i < tests.size(); i++){
            if (tests(i) == "gammaOpt"){gammaOpt();}
            if (tests(i) == "precondTest"){precondTest();}
            if (tests(i) == "precondTestTol"){precondTestTol();}
        }
    }

    void makeTests(gsVector<std::string> testy)
    {
        tests = testy;
        makeTests();
    }

};


#endif // UWBLINEARTESTING_H
