/*
    Author(s): J. Egermaier
*/

#pragma once

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <gismo.h>
#include "uwbINSSolverSteady.h"
#include "uwbINSSolverUnsteady.h"
#include "uwbINSSolverSteadyIterative.h"
#include "uwbINSSolverUnsteadyIterative.h"

using namespace gismo;

template<class T>
void printSTDvector(std::vector<T> v, std::string popis)
{
    gsInfo << popis;
    for (int i = 0; i < v.size(); i++){
        gsInfo << v[i];
        if (i < v.size()-1){
            gsInfo << ", ";
        }
    }
    gsInfo << "\n";
}

template<class T>
void printSTDvector(std::ofstream& file, std::vector<T> v, std::string popis)
{
    file << popis;
    for (int i = 0; i < v.size(); i++){
        file << v[i];
        if (i < v.size()-1){
            file << ", ";
        }
    }
    file << "\n";
}

template<class T>
void defineBCs_profile(gsBoundaryConditions<T>& bcInfo, bool periodic, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndOut, std::vector<std::pair<int, boxSide> >& bndWall)
{
    gsFunctionExpr<> Uin, Uwall;

        Uin = gsFunctionExpr<>("7.765004348","-0.28272",2);//Uin = gsFunctionExpr<>("(-4*(y-1.5)^2 + 1)", "0", 2);
        Uwall = gsFunctionExpr<>("0", "0", 2);

        if (periodic)
        {
            bndIn.push_back(std::make_pair(0, boundary::west));
            bndWall.push_back(std::make_pair(1, boundary::north));
            bndWall.push_back(std::make_pair(1, boundary::south));
            bndOut.push_back(std::make_pair(2, boundary::east));

            //bcInfo.addPeriodic(0, boundary::south, 0, boundary::north, 2);
            //bcInfo.addPeriodic(2, boundary::south, 2, boundary::north, 2);
            //bcInfo.setIdentityMatrix(2);
        }
        else
        {
            bndIn.push_back(std::make_pair(0, boundary::west));
            bndWall.push_back(std::make_pair(1, boundary::north));
            bndWall.push_back(std::make_pair(1, boundary::south));
            bndWall.push_back(std::make_pair(0, boundary::north));
            bndWall.push_back(std::make_pair(0, boundary::south));
            bndWall.push_back(std::make_pair(2, boundary::north));
            bndWall.push_back(std::make_pair(2, boundary::south));
            bndOut.push_back(std::make_pair(2, boundary::east));
        }

    for(size_t i = 0; i < bndWall.size(); i++)
        bcInfo.addCondition(bndWall[i].first, bndWall[i].second, condition_type::dirichlet, Uwall);

    for (size_t i = 0; i < bndIn.size(); i++)
        bcInfo.addCondition(bndIn[i].first, bndIn[i].second, condition_type::dirichlet, Uin);
}


template<class T>
void compute(uwbINSSolverBase<T>& solver, std::ofstream& file, std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, gsOptionList opt)
{
    solver.solve(opt.getInt("maxIter"), opt.getReal("solveTol"));
    std::string gType = opt.getString("geomType");

    std::vector<gsMultiBasis<real_t>> VB = solver.getAssembler()->getBlockAssembler().getBases();
    gsFileData<> dbp;
    dbp << VB[1];
    dbp.save(opt.getString("domainPath") + gType + opt.getString("idStr") + "_discreteBasis_P.xml");
    gsFileData<> dbu;
    dbu << VB[0];
    dbu.save(opt.getString("domainPath") + gType + opt.getString("idStr") + "_discreteBasis_U.xml");

    gsInfo << "Assembly time:" << solver.getAssemblyTime() << "\n";
    gsInfo << "Factorization time:" << solver.getSolverSetupTime() << "\n";
    gsInfo << "Solve time:" << solver.getSolveTime() << "\n";
    file << "Assembly time:" << solver.getAssemblyTime() << "\n";
    file << "Factorization time:" << solver.getSolverSetupTime() << "\n";
    file << "Solve time:" << solver.getSolveTime() << "\n\n";

    std::string matPath = opt.getString("matPath");
    std::string idStr = opt.getString("idStr");
    //std::string massIdStr = opt.getString("massIdStr");

    if (opt.getSwitch("saveMat"))
    {
        gsInfo << "Saving matrices... \n";

        gsFileData<> fd;
        fd << solver.getAssembler()->matrix();
        fd << solver.getAssembler()->rhs();
        fd.save(matPath + gType + idStr + "_matNS.xml");

        solver.getAssembler()->getBlockAssembler().assembleMassMatrix();
        solver.getAssembler()->getBlockAssembler().assemblePressureMassMatrix();

        fd.clear();
        gsSparseMatrix<> sm;
        gsMatrix<> mrhs;
        sm = solver.getAssembler()->getVelocityMassMatrix();
        int radku = sm.rows();
        mrhs = solver.getAssembler()->rhs();
        gsMatrix<> cut (radku, 1);
        cut.col(0) = mrhs.middleRows(0, radku);
        fd << sm;
        fd << cut;
        fd.save(matPath + gType + idStr + "_matVelMass.xml");

        fd.clear();
        fd << solver.getAssembler()->getPressureMassMatrix();
        fd.save(matPath + gType + idStr + "_matPresMass.xml");

        gsSparseMatrix<> Ap, Fp;
        solver.getAssembler()->getPCDblocks(Ap, Fp, bndIn, bndOut, bndWall);

        fd.clear();
        fd << Ap;
        fd.save(matPath + gType + idStr + "_matAp.xml");

        fd.clear();
        fd << Fp;
        fd.save(matPath + gType + idStr + "_matFp.xml");
    }

    if (opt.getSwitch("plot"))
    {
        gsField<> velocity = solver.constructSolution(0);
        gsField<> pressure = solver.constructSolution(1);

        // Write solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, opt.getString("outPath") + gType + idStr + "_vel", opt.getInt("plotPts"), true);
        gsWriteParaview<>(pressure, opt.getString("outPath") + gType + idStr + "_pres", opt.getInt("plotPts"));
    }
}

template <class T>
void computeIterSt(uwbINSSolverSteadyIterative<T, uwbGMResRight<T>>& solver, std::ofstream& file, std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, gsOptionList opt)
{
    solver.solve(opt.getInt("maxIter"), opt.getReal("solveTol"));
    gsInfo << "Assembly time:" << solver.getAssemblyTime() << "\n";
    gsInfo << "Factorization time:" << solver.getSolverSetupTime() << "\n";
    gsInfo << "Solve time:" << solver.getSolveTime() << "\n";
    std::vector<int> iCounts = solver.getLinIterVector();
    printSTDvector<int>(iCounts, "Iteration counts: ");

    file << "Assembly time:" << solver.getAssemblyTime() << "\n";
    file << "Factorization time:" << solver.getSolverSetupTime() << "\n";
    file << "Solve time:" << solver.getSolveTime() << "\n";
    printSTDvector<int>(file, iCounts, "Iteration counts: ");

    std::string gType = opt.getString("geomType");
    std::string matPath = opt.getString("matPath");
    std::string idStr = opt.getString("idStr");

    if (opt.getSwitch("saveMatIter"))
    {
        gsInfo << "Saving matrices... \n";

        gsFileData<> fd;
        fd << solver.getAssembler()->matrix();
        fd << solver.getAssembler()->rhs();
        fd.save(matPath + gType + idStr + "_matNS.xml");

        solver.getAssembler()->getBlockAssembler().assembleMassMatrix();
        solver.getAssembler()->getBlockAssembler().assemblePressureMassMatrix();

        fd.clear();
        gsSparseMatrix<> sm;
        gsMatrix<> mrhs;
        sm = solver.getAssembler()->getVelocityMassMatrix();
        int radku = sm.rows();
        mrhs = solver.getAssembler()->rhs();
        gsMatrix<> cut (radku, 1);
        cut.col(0) = mrhs.middleRows(0, radku);
        fd << sm;
        fd << cut;
        fd.save(matPath + gType + idStr + "_matVelMass.xml");

        fd.clear();
        fd << solver.getAssembler()->getPressureMassMatrix();
        fd.save(matPath + gType + idStr + "_matPresMass.xml");

        gsSparseMatrix<> Ap, Fp;
        solver.getAssembler()->getPCDblocks(Ap, Fp, bndIn, bndOut, bndWall);

        fd.clear();
        fd << Ap;
        fd.save(matPath + gType + idStr + "_matAp.xml");

        fd.clear();
        fd << Fp;
        fd.save(matPath + gType + idStr + "_matFp.xml");
    }

    if (opt.getSwitch("plot"))
    {
        gsField<> velocity = solver.constructSolution(0);
        gsField<> pressure = solver.constructSolution(1);

        // Write solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, opt.getString("outPath") + gType + idStr + "_vel", opt.getInt("plotPts"), true);
        gsWriteParaview<>(pressure, opt.getString("outPath") + gType + idStr + "_pres", opt.getInt("plotPts"));
    }
}

template <class T>
void computeIterUnst(uwbINSSolverUnsteadyIterative<T, uwbGMResRight<T>>& solver, std::ofstream& file, std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, gsOptionList opt)
{
    solver.solve(opt.getInt("maxIter"), opt.getReal("solveTol"));
    gsInfo << "Assembly time:" << solver.getAssemblyTime() << "\n";
    gsInfo << "Factorization time:" << solver.getSolverSetupTime() << "\n";
    gsInfo << "Solve time:" << solver.getSolveTime() << "\n";
    std::vector<int> iCounts = solver.getLinIterVector();
    printSTDvector<int>(iCounts, "Iteration counts: ");

    file << "Assembly time:" << solver.getAssemblyTime() << "\n";
    file << "Factorization time:" << solver.getSolverSetupTime() << "\n";
    file << "Solve time:" << solver.getSolveTime() << "\n";
    printSTDvector<int>(file, iCounts, "Iteration counts: ");

    std::string gType = opt.getString("geomType");
    std::string matPath = opt.getString("matPath");
    std::string idStr = opt.getString("idStr");

    if (opt.getSwitch("saveMatIter"))
    {
        gsInfo << "Saving matrices... \n";

        gsFileData<> fd;
        fd << solver.getAssembler()->matrix();
        fd << solver.getAssembler()->rhs();
        fd.save(matPath + gType + idStr + "_matNS.xml");

        solver.getAssembler()->getBlockAssembler().assembleMassMatrix();
        solver.getAssembler()->getBlockAssembler().assemblePressureMassMatrix();

        fd.clear();
        gsSparseMatrix<> sm;
        gsMatrix<> mrhs;
        sm = solver.getAssembler()->getVelocityMassMatrix();
        int radku = sm.rows();
        mrhs = solver.getAssembler()->rhs();
        gsMatrix<> cut (radku, 1);
        cut.col(0) = mrhs.middleRows(0, radku);
        fd << sm;
        fd << cut;
        fd.save(matPath + gType + idStr + "_matVelMass.xml");

        fd.clear();
        fd << solver.getAssembler()->getPressureMassMatrix();
        fd.save(matPath + gType + idStr + "_matPresMass.xml");

        gsSparseMatrix<> Ap, Fp;
        //solver.getAssembler()->getPCDblocks(Ap, Fp, bndIn, bndOut, bndWall);

        fd.clear();
        fd << Ap;
        fd.save(matPath + gType + idStr + "_matAp.xml");

        fd.clear();
        fd << Fp;
        fd.save(matPath + gType + idStr + "_matFp.xml");
    }

    if (opt.getSwitch("plot"))
    {
        gsField<> velocity = solver.constructSolution(0);
        gsField<> pressure = solver.constructSolution(1);

        // Write solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, opt.getString("outPath") + gType + idStr + "_vel", opt.getInt("plotPts"), true);
        gsWriteParaview<>(pressure, opt.getString("outPath") + gType + idStr + "_pres", opt.getInt("plotPts"));
    }
}

template<class T> void generateMatrices(std::map<std::string, gsVector<std::string>> paramList, gsMultiPatch<T> patches, gsMultiBasis<T> tbasis, std::string outPath)
{
    int numRefine, numRefineU, numRefineCorner, numRefineWalls, deg, numElevate, elemType;
    int maxIter, maxPicardIt, dim, linMaxIt;
    real_t viscosity, timeStep, solveTol, picardTol, gamma, linTol;
    bool plot, saveMat, steady, periodic, iterative, saveMatIter;
    std::string geomType, mainDir, matPath, domainPath, iterDir;
    gsVector<std::string> precType;

    get_parameter(paramList,numRefine,"GnumRef");;
    get_parameter(paramList,numRefineU,"GnumRefU");
    get_parameter(paramList,numRefineWalls,"GnumRefW");
    get_parameter(paramList,numRefineCorner,"GnumRefC");
    get_parameter(paramList,deg,"deg");
    get_parameter(paramList,numElevate,"degElev");
    get_parameter(paramList,elemType,"elemType");
    get_parameter(paramList,maxIter,"GmaxIt");
    get_parameter(paramList,maxPicardIt,"GmaxPicardIt");
    get_parameter(paramList,viscosity,"visc");
    get_parameter(paramList,timeStep,"timeStep");
    get_parameter(paramList,solveTol,"Gtol");
    get_parameter(paramList,linTol,"GlinTol");
    get_parameter(paramList,picardTol,"GtolPicard");
    get_parameter(paramList,linMaxIt,"GlinMaxIt");
    get_parameter(paramList,plot,"Gplot");
    get_parameter(paramList,saveMat,"saveMat");
    get_parameter(paramList,saveMatIter,"saveMatIter");
    get_parameter(paramList,steady,"steady");
    get_parameter(paramList,periodic,"periodic");
    get_parameter(paramList,iterative,"iterative");
    get_parameter(paramList,precType,"Gprec");
    get_parameter(paramList,geomType,"geometry_type");
    get_parameter(paramList,gamma,"Gprecgamma");
    get_parameter(paramList,iterDir,"precondTestTolDir");
    get_parameter(paramList,mainDir,"mainDir");
    get_parameter(paramList,matPath,"matrixDir"); matPath = mainDir + matPath;
    get_parameter(paramList,domainPath,"domainDir"); domainPath = mainDir + domainPath;

    int plotPts = 20000;

    std::string viscStr = "_visc_" + util::to_string(viscosity);
    std::replace(viscStr.begin(), viscStr.end(), '.', '-');
    viscStr.erase(viscStr.find_last_not_of('0') + 1, std::string::npos);

    std::string dtStr = "_dt_" + util::to_string(timeStep);
    std::replace(dtStr.begin(), dtStr.end(), '.', '-');
    dtStr.erase(dtStr.find_last_not_of('0') + 1, std::string::npos);

    std::string refMiddle;
    if (geomType == "profile2D"){
        refMiddle = util::to_string(numRefineU);
    } else {
        refMiddle = util::to_string(numRefineCorner);
    }
    std::string refStr = "_ref_" + util::to_string(numRefine) + "_" + refMiddle + "_" + util::to_string(numRefineWalls);
    std::string degStr = "_deg_" + util::to_string(deg);
    if (elemType == 0)
        degStr += "_elev_" + util::to_string(numElevate);
    else if (elemType == 1)
        degStr += "_SG";

    std::string idStr = viscStr + refStr + degStr;
    std::string massIdStr =  refStr + degStr;

    if (steady)
        idStr += "_st";
    else
        idStr += "_unst" + dtStr;

    if (numRefine < 0)
    {
        gsInfo << "Number of refinements must be non-negative, using 0.\n";
        numRefine = 0;
    }

    gsAssemblerOptions opt;
    opt.dirStrategy = dirichlet::elimination;
    opt.dirValues = dirichlet::l2Projection;
    opt.intStrategy = iFace::glue;

    std::ofstream file;
    if (iterative){
        file.open(mainDir+iterDir+"iter_"+geomType+idStr+".txt");
    } else {
        file.open(outPath+"iter_"+geomType+idStr+".txt");
    }
    file << "Generating matrices for the " + geomType + ".\n\n";
    file << "periodic domain: " << periodic << "\n";
    file << "uniform refine: " << numRefine << "\n";
    file << "refineU: " << numRefineU << "\n";
    file << "refineLocal: " << numRefineWalls << "\n";
    file << "refineCorner: " << numRefineCorner << "\n";
    file << "geometry basis degree: " << deg << "\n";
    file << "basis elevate: " << numElevate << "\n";
    file << "element type (0 - Taylor-Hood, 1 - SG): " << elemType << "\n";
    file << "viscosity: " << viscosity << "\n";
    file << "time step: " << timeStep << "\n";
    file << "max steps: " << maxIter << "\n";
    file << "max Picard iter: " << maxPicardIt << "\n";
    file << "tol: " << solveTol << "\n";
    file << "Picard tol: " << picardTol << "\n";
    file << "lin tol: " << linTol << "\n";

    // ========================================= Define geometry + BCs =========================================

    gsBoundaryConditions<> bcInfo, bcInfoP;
    std::vector<std::pair<int, boxSide> > bndIn;
    std::vector<std::pair<int, boxSide> > bndOut;
    std::vector<std::pair<int, boxSide> > bndWall;
    gsFunctionExpr<> f;

    if (geomType == "step3D"){
        dim = 3;
        f = gsFunctionExpr<>("0", "0", "0", 3);
        defineBCs_step<T>(bcInfo, bndIn, bndOut, bndWall, dim, false, "default");
    } else {
        dim = 2;
        f = gsFunctionExpr<>("0", "0", 2);
        if (geomType == "step2D"){
            defineBCs_step<T>(bcInfo, bndIn, bndOut, bndWall, dim, false, "default");
        } else {
            defineBCs_profile(bcInfo, periodic, bndIn, bndOut, bndWall);
        }
    }

    gsFileData<> bcU, bcP;
    bcU << bcInfo;
    bcU.save(domainPath + geomType + idStr + "_BCU.xml");
    bcP << bcInfoP;
    bcP.save(domainPath + geomType + idStr + "_BCP.xml");

    std::vector< gsMultiBasis<> >  discreteBases;

    switch (elemType)
    {
    case 0:
    default:
        tbasis.degreeElevate(numElevate);
        discreteBases.push_back(tbasis);//Basis for velocity
        discreteBases.push_back(tbasis);//Basis for pressure
        discreteBases[0].degreeElevate(1);//p-refine the velocity space (Taylor-Hood type)
        break;
    case 1:
        //refineBasis(tbasis, numRefine, numRefineU, numRefineLocal, addRefPart, a, b);
        discreteBases.push_back(tbasis);//Basis for velocity
        discreteBases.push_back(tbasis);//Basis for pressure
        discreteBases[0].degreeIncrease(1);
        discreteBases[0].uniformRefine();
        break;
    }

    gsInfo << "Velocity basis degree: " << discreteBases[0].degree() << "\n";
    gsInfo << "Pressure basis degree: " << discreteBases[1].degree() << "\n";
    file << "Velocity basis degree: " << discreteBases[0].degree() << "\n";

    // ========================================= Define solver ========================================= 
  
    uwbINSPde<real_t> NSpde(patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);

    params.settings().set(constantsINS::timeStep, timeStep);
    params.settings().set(constantsINS::unst_innerIt, maxPicardIt);
    params.settings().set(constantsINS::unst_innerTol, picardTol);

    gsOptionList computeOpt;
    computeOpt.addInt("maxIter", "", maxIter);
    computeOpt.addInt("plotPts", "", plotPts);
    computeOpt.addReal("solveTol", "", solveTol);
    computeOpt.addSwitch("plot", "", plot);
    computeOpt.addSwitch("saveMat", "", saveMat);
    computeOpt.addSwitch("saveMatIter", "", saveMatIter);
    computeOpt.addString("outPath", "", outPath);
    computeOpt.addString("matPath", "", matPath);
    computeOpt.addString("domainPath","",domainPath);
    computeOpt.addString("idStr", "", idStr);
    computeOpt.addString("massIdStr", "", massIdStr);
    computeOpt.addString("geomType","",geomType);


    // ========================================= Solving =========================================  

    if (steady){
        if (iterative){
            params.getPrecOptions().setReal("gamma",gamma);
            params.settings().set(constantsINS::iter_maxIt, linMaxIt);
            params.settings().set(constantsINS::iter_tol, linTol);
            for (int i = 0; i < precType.size(); i++){
                params.settings().setPrecondType(precType(i));
                uwbINSSolverSteadyIterative<real_t, uwbGMResRight<real_t>> navStokes(params);
                if (precType(i)=="PCD_AdiagEqual"){
                    navStokes.getAssembler()->preparePCDboundary(bndIn, bndOut, bndWall);
                }
                navStokes.initialize();
                gsInfo << "numDofs: " << navStokes.numDofs() << "\n";
                file << "numDofs: " << navStokes.numDofs() << "\n";
                gsInfo << "initialization...\n";
                file << "\n";
                file << "Using preconditioner: " << precType(i) << "\n";
                computeIterSt(navStokes, file, bndIn, bndOut, bndWall, computeOpt);
            }
      } else {
            uwbINSSolverSteady<real_t> navStokes(params);
            gsInfo << "numDofs: " << navStokes.numDofs() << "\n";
            file << "numDofs: " << navStokes.numDofs() << "\n\n";
            gsInfo << "initialization...\n";
            navStokes.initialize();
            compute(navStokes, file, bndIn, bndOut, bndWall, computeOpt);
        }
    } else {
        if (iterative){
            params.getPrecOptions().setReal("gamma",gamma);
            params.settings().set(constantsINS::iter_maxIt, linMaxIt);
            params.settings().set(constantsINS::iter_tol, linTol);
            for (int i = 0; i < precType.size(); i++){
                params.settings().setPrecondType(precType(i));
                uwbINSSolverUnsteadyIterative<real_t, uwbGMResRight<real_t> > navStokes(params);
                if (precType(i)=="PCD_AdiagEqual"){
                    navStokes.getAssembler()->preparePCDboundary(bndIn, bndOut, bndWall);
                }
                navStokes.initialize();
                navStokes.setStokesInitialCondition();
                gsInfo << "numDofs: " << navStokes.numDofs() << "\n";
                file << "numDofs: " << navStokes.numDofs() << "\n";
                gsInfo << "initialization...\n";
                file << "\n";
                file << "Using preconditioner: " << precType(i) << "\n";
                computeIterUnst(navStokes, file, bndIn, bndOut, bndWall, computeOpt);
            }
        } else {
        uwbINSSolverUnsteady<real_t> navStokes(params);
        gsInfo << "numDofs: " << navStokes.numDofs() << "\n";
        file << "numDofs: " << navStokes.numDofs() << "\n\n";
        gsInfo << "initialization...\n";
        navStokes.initialize();
        navStokes.setStokesInitialCondition();
        compute(navStokes, file, bndIn, bndOut, bndWall, computeOpt);
        }
    }

    file.close();
}

