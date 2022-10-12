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
#include "myGeometryUtils.h"

using namespace gismo;

template<class T> void compute_profile(uwbINSSolverBase<T>& solver, std::ofstream& file, std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, gsOptionList opt);
template<class T> void defineBCs_profile(gsBoundaryConditions<T>& bcInfo, bool periodic, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndOut, std::vector<std::pair<int, boxSide> >& bndWall);
//template<class T> void refineBasis(gsMultiBasis<T>& basis, int numRefine, int numRefineU, int numRefineLocal, real_t addRefPart, real_t a, real_t b);

void generateMatricesProfile2D(std::map<std::string, gsVector<std::string>> paramList, gsMultiPatch<> patches, gsMultiBasis<> tbasis, std::string outPath, std::string matPath)
{
    int numRefine, numRefineU, numRefineLocal, numElevate, elemType;
    int maxIter, maxPicardIt, deg;
    real_t viscosity, timeStep, solveTol, picardTol;
    bool plot, saveMat, steady, periodic;
    std::string geomType;

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

    gsBoundaryConditions<> bcInfo;
    std::vector<std::pair<int, boxSide> > bndIn;
    std::vector<std::pair<int, boxSide> > bndOut;
    std::vector<std::pair<int, boxSide> > bndWall;

    defineBCs_profile(bcInfo, periodic, bndIn, bndOut, bndWall);
    
    gsFunctionExpr<> f;
    //============================
    //gsReadFile<> ("/home/jirieggy/gismoQT/build/filePatches_beforeBasisRefine.xml", patches);
    f = gsFunctionExpr<>("0", "0", 2);
    //============================

    // Define discretization space by refining the basis of the geometry

    std::vector< gsMultiBasis<> >  discreteBases;

    switch (elemType)
    {
    case 0:
    default:
        //refineBasis(tbasis, numRefine, numRefineU, numRefineLocal, addRefPart, a, b);
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
        compute_profile(navStokes, file, bndIn, bndOut, bndWall, computeOpt);
    }
    else
    {
        uwbINSSolverUnsteady<real_t> navStokes(params);
        gsInfo << "numDofs: " << navStokes.numDofs() << "\n";
        file << "numDofs: " << navStokes.numDofs() << "\n\n";

        gsInfo << "initialization...\n";
        navStokes.initialize();
        navStokes.setStokesInitialCondition();
        compute_profile(navStokes, file, bndIn, bndOut, bndWall, computeOpt);
    }

    file.close();
}

template<class T>
void compute_profile(uwbINSSolverBase<T>& solver, std::ofstream& file, std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, gsOptionList opt)
{
    solver.solve(opt.getInt("maxIter"), opt.getReal("solveTol"));

    gsInfo << "Assembly time:" << solver.getAssemblyTime() << "\n";
    gsInfo << "Factorization time:" << solver.getSolverSetupTime() << "\n";
    gsInfo << "Solve time:" << solver.getSolveTime() << "\n";
    file << "Assembly time:" << solver.getAssemblyTime() << "\n";
    file << "Factorization time:" << solver.getSolverSetupTime() << "\n";
    file << "Solve time:" << solver.getSolveTime() << "\n\n";

    std::string matPath = opt.getString("matPath");
    std::string idStr = opt.getString("idStr");
    //std::string massIdStr = opt.getString("massIdStr");
    std::string geomType = opt.getString("geomType");

    if (opt.getSwitch("saveMat"))
    {
        gsInfo << "Saving matrices... \n";

        gsFileData<> fd;
        fd << solver.getAssembler()->matrix();
        fd << solver.getAssembler()->rhs();
        fd.save(matPath + geomType + idStr + "_matNS.xml");

        solver.getAssembler()->getBlockAssembler().assembleMassMatrix();
        solver.getAssembler()->getBlockAssembler().assemblePressureMassMatrix();

        fd.clear();
        fd << solver.getAssembler()->getVelocityMassMatrix();
        fd.save(matPath + geomType + idStr + "_matVelMass.xml");
//        fd.save(matPath + geomType + massIdStr + "_matVelMass.xml");

        fd.clear();
        fd << solver.getAssembler()->getPressureMassMatrix();
        fd.save(matPath + geomType + idStr + "_matPresMass.xml");
//        fd.save(matPath + geomType + massIdStr + "_matPresMass.xml");

        gsSparseMatrix<> Ap, Fp;
        solver.getAssembler()->getPCDblocks(Ap, Fp, bndIn, bndOut, bndWall);

        fd.clear();
        fd << Ap;
        fd.save(matPath + geomType + idStr + "_matAp.xml");

        fd.clear();
        fd << Fp;
        fd.save(matPath + geomType + idStr + "_matFp.xml");
    }

    if (opt.getSwitch("plot"))
    {
        gsField<> velocity = solver.constructSolution(0);
        gsField<> pressure = solver.constructSolution(1);

        // Write solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, opt.getString("outPath") + geomType + idStr + "_vel", opt.getInt("plotPts"), true);
        gsWriteParaview<>(pressure, opt.getString("outPath") + geomType + idStr + "_pres", opt.getInt("plotPts"));
    }

}

template<class T> void defineBCs_profile(gsBoundaryConditions<T>& bcInfo, bool periodic, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndOut, std::vector<std::pair<int, boxSide> >& bndWall)
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
/*
template<class T> void refineBasis(gsMultiBasis<T>& basis, int numRefine, int numRefineU, int numRefineLocal, real_t addRefPart, real_t a, real_t b)
{
    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();

    int bRefine = math::floor((a / b) - 1) + numRefineU;

        const gsTensorBSplineBasis<2, T>*  basis0 = dynamic_cast<const gsTensorBSplineBasis<2, T>*>(&(basis.basis(0))); //basis of the patch 0

        gsMatrix<> box(2, 2);
        box << 0, 1, 0, 0;
        for (int i = 0; i < bRefine; i++)
        {
            basis.refine(0, box);
            basis.refine(1, box);
        }

        box << 0, addRefPart, 0, 0;
        basis.refine(0, box);
        basis.refine(1, box);
        
        for (int i = 0; i < numRefineLocal; i++)
        {
            int sizeKnots_v = basis0->knots(1).size() - 1;
            real_t lastKnot_v = basis0->knot(1, sizeKnots_v);
            real_t vKnot_bottom = basis0->knot(1, basis0->degree(1) + 1);
            real_t vKnot_upper = basis0->knot(1, sizeKnots_v - (basis0->degree(1) + 1));

            box << 0, 0, 0, vKnot_bottom;
            basis.refine(0, box);
            basis.refine(1, box);
            basis.refine(2, box);

            box << 0, 0, vKnot_upper, lastKnot_v;
            basis.refine(0, box);
            basis.refine(1, box);
            basis.refine(2, box);
        }
        break;
    }
}
*/
