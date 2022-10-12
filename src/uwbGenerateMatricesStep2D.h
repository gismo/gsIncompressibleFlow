#pragma once

#include <iostream>
#include <ctime>
#include <cstdlib>

#include <gismo.h>

#include "uwbINSSolverSteady.h"
#include "uwbINSSolverUnsteady.h"

using namespace gismo;

template<class T>
void compute(uwbINSSolverBase<T>& solver, std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, gsOptionList opt)
{
    solver.solve(opt.getInt("maxIter"), opt.getReal("solveTol"));

    std::vector<gsMultiBasis<real_t>> VB = solver.getAssembler()->getBlockAssembler().getBases();
    gsFileData<> dbp;
    dbp << VB[1];
    dbp.save(opt.getString("domainPath") + "step2D" + opt.getString("idStr") + "_discreteBasis_P.xml");
    gsFileData<> dbu;
    dbu << VB[0];
    dbu.save(opt.getString("domainPath") + "step2D" + opt.getString("idStr") + "_discreteBasis_U.xml");

    gsInfo << "Assembly time:" << solver.getAssemblyTime() << "\n";
    gsInfo << "Factorization time:" << solver.getSolverSetupTime() << "\n";
    gsInfo << "Solve time:" << solver.getSolveTime() << "\n";

    std::string matPath = opt.getString("matPath");
    std::string idStr = opt.getString("idStr");
    //std::string massIdStr = opt.getString("massIdStr");

    if (opt.getSwitch("saveMat"))
    {
        gsInfo << "Saving matrices... \n";

        gsFileData<> fd;
        fd << solver.getAssembler()->matrix();
        fd << solver.getAssembler()->rhs();
        fd.save(matPath + "step2D" + idStr + "_matNS.xml");

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
        fd.save(matPath + "step2D" + idStr + "_matVelMass.xml");

        fd.clear();
        fd << solver.getAssembler()->getPressureMassMatrix();
        fd.save(matPath + "step2D" + idStr + "_matPresMass.xml");

        gsSparseMatrix<> Ap, Fp;
        solver.getAssembler()->getPCDblocks(Ap, Fp, bndIn, bndOut, bndWall);

        fd.clear();
        fd << Ap;
        fd.save(matPath + "step2D" + idStr + "_matAp.xml");

        fd.clear();
        fd << Fp;
        fd.save(matPath + "step2D" + idStr + "_matFp.xml");
    }

    if (opt.getSwitch("plot"))
    {
        gsField<> velocity = solver.constructSolution(0);
        gsField<> pressure = solver.constructSolution(1);

        // Write solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, opt.getString("outPath") + "step2D" + idStr + "_vel", opt.getInt("plotPts"), true);
        gsWriteParaview<>(pressure, opt.getString("outPath") + "step2D" + idStr + "_pres", opt.getInt("plotPts"));
    }

}

template<class T> void defineBCs(gsBoundaryConditions<T>& bcInfo, bool periodic, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndOut, std::vector<std::pair<int, boxSide> >& bndWall)
{
    gsFunctionExpr<T> Uin, Uwall;

    Uin = gsFunctionExpr<T>("(-4*(y-1.5)^2 + 1)", "0", 2);
    Uwall = gsFunctionExpr<T>("0", "0", 2);

    if (periodic){
        bndIn.push_back(std::make_pair(2, boundary::west));
        bndWall.push_back(std::make_pair(0, boundary::west));
        bndWall.push_back(std::make_pair(2, boundary::south));
        bndWall.push_back(std::make_pair(2, boundary::north));
        bndOut.push_back(std::make_pair(0, boundary::east));
        bndOut.push_back(std::make_pair(1, boundary::east));

        bcInfo.addPeriodic(0, boundary::south, 1, boundary::north, 2);
        bcInfo.setIdentityMatrix(2);
    } else {
        bndIn.push_back(std::make_pair(2, boundary::west));
        bndWall.push_back(std::make_pair(0, boundary::west));
        bndWall.push_back(std::make_pair(0, boundary::south));
        bndWall.push_back(std::make_pair(1, boundary::north));
        bndWall.push_back(std::make_pair(2, boundary::south));
        bndWall.push_back(std::make_pair(2, boundary::north));
        bndOut.push_back(std::make_pair(0, boundary::east));
        bndOut.push_back(std::make_pair(1, boundary::east));
    }

    for(size_t i = 0; i < bndWall.size(); i++)
        bcInfo.addCondition(bndWall[i].first, bndWall[i].second, condition_type::dirichlet, Uwall);

    for (size_t i = 0; i < bndIn.size(); i++)
        bcInfo.addCondition(bndIn[i].first, bndIn[i].second, condition_type::dirichlet, Uin);
}


template<class T> void generateMatricesStep2D(std::map<std::string, gsVector<std::string>> paramList, gsMultiPatch<T> patches, gsMultiBasis<T> tbasis, std::string outPath)
{
    int numRefine, numRefineU, numRefineLocal, deg, numElevate, elemType;
    int maxIter, maxPicardIt;
    real_t viscosity, timeStep, solveTol, picardTol;
    bool plot, saveMat, steady, periodic;
    std::string geomType, mainDir, matPath, domainPath;

    get_parameter(paramList,numRefine,"GnumRef");;
    get_parameter(paramList,numRefineU,"GnumRefU");
    get_parameter(paramList,numRefineLocal,"GnumRefW");
    get_parameter(paramList,deg,"deg");
    get_parameter(paramList,numElevate,"degElev");
    get_parameter(paramList,elemType,"elemType");
    get_parameter(paramList,maxIter,"GmaxIt");
    get_parameter(paramList,maxPicardIt,"GmaxPicardIt");
    get_parameter(paramList,viscosity,"visc");
    get_parameter(paramList,timeStep,"timeStep");
    get_parameter(paramList,solveTol,"Gtol");
    get_parameter(paramList,picardTol,"GtolPicard");
    get_parameter(paramList,plot,"Gplot");
    get_parameter(paramList,saveMat,"saveMat");
    get_parameter(paramList,steady,"steady");
    get_parameter(paramList,periodic,"periodic");
    get_parameter(paramList,geomType,"geometry_type");
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

    std::string refStr = "_ref_" + util::to_string(numRefine) + "_" + util::to_string(numRefineU) + "_" + util::to_string(numRefineLocal);
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

    gsInfo << "Solving the step 2D.\n";

    std::ofstream file;
    file.open(outPath+geomType+idStr+".txt");
    file << "Generating matrices for the " + geomType + ".\n\n";
    file << "periodic domain: " << periodic << "\n";
    file << "uniform refine: " << numRefine << "\n";
    file << "refineU: " << numRefineU << "\n";
    file << "refineLocal: " << numRefineLocal << "\n";
    file << "geometry basis degree: " << deg << "\n";
    file << "basis elevate: " << numElevate << "\n";
    file << "element type (0 - Taylor-Hood, 1 - SG): " << elemType << "\n";
    file << "viscosity: " << viscosity << "\n";
    file << "time step: " << timeStep << "\n";
    file << "max steps: " << maxIter << "\n";
    file << "max Picard iter: " << maxPicardIt << "\n";
    file << "tol: " << solveTol << "\n";
    file << "Picard tol: " << picardTol << "\n\n";

    // ========================================= Define geometry + BCs ========================================= 

    gsBoundaryConditions<> bcInfo, bcInfoP;
    std::vector<std::pair<int, boxSide> > bndIn;
    std::vector<std::pair<int, boxSide> > bndOut;
    std::vector<std::pair<int, boxSide> > bndWall;

    defineBCs(bcInfo, periodic, bndIn, bndOut, bndWall);

    gsFileData<> bcU, bcP;
    bcU << bcInfo;
    bcU.save(domainPath + "step2D" + idStr + "_BCU.xml");
    bcP << bcInfoP;
    bcP.save(domainPath + "step2D" + idStr + "_BCP.xml");

    gsFunctionExpr<> f;
    f = gsFunctionExpr<>("0", "0", 2);

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
    computeOpt.addString("outPath", "", outPath);
    computeOpt.addString("matPath", "", matPath);
    computeOpt.addString("domainPath","",domainPath);
    computeOpt.addString("idStr", "", idStr);
    computeOpt.addString("massIdStr", "", massIdStr);
    computeOpt.addString("geomType","",geomType);


    // ========================================= Solving =========================================  

    if (steady)
    {
        uwbINSSolverSteady<real_t> navStokes(params);
        gsInfo << "numDofs: " << navStokes.numDofs() << "\n";
        file << "numDofs: " << navStokes.numDofs() << "\n\n";

        gsInfo << "initialization...\n";
        navStokes.initialize();
        compute(navStokes, bndIn, bndOut, bndWall, computeOpt);
    }
    else
    {
        uwbINSSolverUnsteady<real_t> navStokes(params);
        gsInfo << "numDofs: " << navStokes.numDofs() << "\n";
        file << "numDofs: " << navStokes.numDofs() << "\n\n";

        gsInfo << "initialization...\n";
        navStokes.initialize();
        navStokes.setStokesInitialCondition();
        compute(navStokes, bndIn, bndOut, bndWall, computeOpt);
    }

    file.close();
}

