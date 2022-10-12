/** @file uwbRANSExamplesSetting.h

Author(s): E. Turnerova, K. Michalkova
*/

#pragma once

#include "uwbGeometryCreators.h"
#include "uwbTMSolverKOmega.h"
#include "uwbINSSolverSteady.h"
#include "uwbTMSolverKOmegaLinSteady.h"

namespace gismo
{

template<class T>
class uwbRANSExampleSetting
{
public:
    uwbRANSExampleSetting(T viscosity, int plot_pts) : m_viscosity(viscosity), m_plot_pts(plot_pts)
    {
        PI = 3.141592653589793238463;
        initMembers();
    }

    ~uwbRANSExampleSetting() { }

protected:
    virtual void initMembers()
    {
        m_sGeometry = "";
    }

public:
    virtual void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, gsFunctionExpr<T> Uin, gsFunctionExpr<T> UWall)
    {
        for (int i = 0; i < m_inletBoundaryPatches.rows(); i++)
            bcInfo.addCondition(m_inletBoundaryPatches(i), m_inletBoundarySides[i], condition_type::dirichlet, Uin, 0);
        for (int i = 0; i < m_wallBoundaryPatches.rows(); i++)
            bcInfo.addCondition(m_wallBoundaryPatches(i), m_wallBoundarySides[i], condition_type::dirichlet, UWall, 0);
    }

    virtual void defineBCs_TM(gsBoundaryConditions<T>& bcInfoTurb, T kIn, T kWall, T oIn, T oWall)
    {
        // Boundary conditions
        gsFunctionExpr<T> Kin(util::to_string(kIn), 2);
        gsFunctionExpr<T> Oin(util::to_string(oIn), 2);
        gsFunctionExpr<T> KWall(util::to_string(kWall), 2);
        gsFunctionExpr<T> OWall(util::to_string(oWall), 2);

        for (int i = 0; i < m_inletBoundaryPatches.rows(); i++)
        {
            bcInfoTurb.addCondition(m_inletBoundaryPatches(i), m_inletBoundarySides[i], condition_type::dirichlet, Kin, 0);
            bcInfoTurb.addCondition(m_inletBoundaryPatches(i), m_inletBoundarySides[i], condition_type::dirichlet, Oin, 1);
        }
        for (int i = 0; i < m_wallBoundaryPatches.rows(); i++)
        {
            bcInfoTurb.addCondition(m_wallBoundaryPatches(i), m_wallBoundarySides[i], condition_type::dirichlet, KWall, 0);
            bcInfoTurb.addCondition(m_wallBoundaryPatches(i), m_wallBoundarySides[i], condition_type::dirichlet, OWall, 1);
        }
    }

    void computeSteadyNS(uwbINSSolverSteady<T>& navStokes, int numIterSteadyNS)
    {
        gsInfo << "Solving Steady case: \n";
        gsInfo << "numDofs: " << navStokes.numDofs() << "\n";

        navStokes.initialize(); // steady solver
        navStokes.solve(numIterSteadyNS, 1e-5);

        gsField<T> velocitySteady = navStokes.constructSolution(0);
        gsField<T> pressureSteady = navStokes.constructSolution(1);

        plotSolutionField(velocitySteady, "steadyVelocity");
        plotSolutionField(pressureSteady, "steadyPressure");

        gsFileData<T> fd;
        fd << navStokes.getSolution();
        fd.save(m_sGeometry + "NSsteadySolution.xml");

        /*gsInfo << "Computing initial velocity by solving Stokes problem...\n";
        uwbINSPde<T> NSpde(*patches, bcInfo, f, m_viscosity);
        uwbINSSolverParams<T> params(NSpde, discreteBases, opt);
        uwbINSSolverSteady<T> stokesSolver(params);
        gsInfo << "numDofs: " << stokesSolver.numDofs() << "\n";
        stokesSolver.initialize();
        stokesSolver.setStokesSolution();*/
    }

    void computeSteadyKOmegaTM(uwbINSSolverSteady<T>& navStokes, uwbTMSolverKOmegaLinSteady<T>& turbSolver, int numIterKOmegaSteady, T numDofs)
    {
        gsStopwatch time;
        T Tassembly;
        T Tsolve;
        gsInfo << "\nComputing steady linearized k-omega...\n";
        gsInfo << "initialization...\n";
        gsField<T> velocity = navStokes.constructSolution(0);
        gsInfo << "TMnumDofs = " << numDofs << "\n";
        turbSolver.initialize(velocity);
        Tassembly = time.stop();
        gsInfo << "Assembly time:" << Tassembly << "\n";
        time.restart();
        turbSolver.solve(numIterKOmegaSteady, 1e-5); // solution change norm tol = 10^(-5)
        Tsolve = time.stop();
        gsInfo << "Solve time:" << Tsolve << "\n";
        gsField<T> kOmegaSteady = turbSolver.constructSolution();
        plotSolutionField(kOmegaSteady, "steadykOmega", true);

        gsFileData<T> fd;
        fd << turbSolver.getSolution();
        fd.save(m_sGeometry + "TMsteadySolution.xml");
    }

    virtual void plotSolutionField(gsField<T> solField, std::string const & fn, bool plotMesh = false)
    {
        gsWriteParaview<T>(solField, m_sGeometry + fn, m_plot_pts, plotMesh);
    }

    void refineBasisLocalZones(gsMultiBasis<T>& tbasis, int id_patch, int direction, int position, int knot, int num_of_refine)
    {
        gsMatrix<T> box(2, 2);
        for (int i = 0; i < num_of_refine; i++)
        {
            // in each refinement step take the first numKnot_v non_zero knots
            const gsTensorBSplineBasis<2, T>*  basis = dynamic_cast<const gsTensorBSplineBasis<2, T>*>(&(tbasis.basis(id_patch))); //basis of the first patch
            gsVector<T> RefineKnot((i + 1)*(knot)+1);
            RefineKnot.setZero((i + 1)*(knot)+1);
            int sizeKnots = basis->knots(direction).size() - 1;
            for (int k = 0; k < (i + 1)*(knot)+1; k++) //gsVector of first numKnots knots
            {
                if (position == 0)
                   RefineKnot(k) = basis->knot(direction, basis->degree(direction) + k); // first numKnot knots and 0 knot in direction v - all patches have same knots in v direction
                else
                   RefineKnot(k) = basis->knot(direction, sizeKnots - (basis->degree(direction) + k)); //last numKnot knots and 1 knot in direction v
            }
            for (int j = 0; j < knot; j++)
            {
               if (direction == 0)
               {
                   if (position == 0)
                       box << RefineKnot(j), RefineKnot(j+1), 0, 0;
                   else
                       box << RefineKnot(j+1), RefineKnot(j), 0, 0;
               }
               else
               {
                   if (position == 0)
                       box << 0,0, RefineKnot(j), RefineKnot(j+1);
                   else
                       box << 0,0, RefineKnot(j+1), RefineKnot(j);
               }

                tbasis.refine(id_patch,  box);
            }
        }
    }

public:
    virtual void setBoundaries() { GISMO_NO_IMPLEMENTATION }
    virtual gsMultiPatch<T> makeMultiPatch() { GISMO_NO_IMPLEMENTATION }
    virtual void refineBasis(gsMultiBasis<T>& tbasis, int numUniformRefine, gsMatrix<T> mRefZones, gsMatrix<T> mRefUniformZones) { GISMO_NO_IMPLEMENTATION }
    virtual T computeWallDistance(uwbINSSolverSteady<T>& navStokes) { GISMO_NO_IMPLEMENTATION }
    virtual void solvePoissonEquation(uwbTMSolverKOmega<T>& turbSolver) { GISMO_NO_IMPLEMENTATION }

    std::string getGeometryName() const { return m_sGeometry; }
    const gsVector<int>& getWallBoundaryPatches() const { return m_wallBoundaryPatches; }
    const std::vector<boxSide> getWallBoundarySides() const { return m_wallBoundarySides; }

protected:
    T m_viscosity;
    int m_plot_pts;
    T PI;

    gsMultiPatch<T> m_patches;
    std::vector<boxSide> m_wallBoundarySides;
    std::vector<boxSide> m_inletBoundarySides;
    gsVector<int> m_wallBoundaryPatches;
    gsVector<int> m_inletBoundaryPatches;

    std::string m_sGeometry;

}; //uwbRANSExampleSetting


template<class T>
class uwbRANSBackwardStep2DExample : public uwbRANSExampleSetting<T>
{
public:
    typedef uwbRANSExampleSetting<T> Base;

public:
    uwbRANSBackwardStep2DExample(T viscosity, T length, T height, T inletLength, T uMax, int plot_pts) :
        Base(viscosity, plot_pts), m_length(length), m_height(height), m_inletLength(inletLength), m_uMax(uMax)
    {
        initMembers();
        setBoundaries();
    }

    ~uwbRANSBackwardStep2DExample() { }

protected:
    void initMembers()
    {
        m_wallBoundaryPatches.setZero(5);
        m_inletBoundaryPatches.setZero(1);
        m_uFreeStream = m_uMax;
        m_sGeometry = "step_";
    }

public:
    void setBoundaries()
    {
        //walls at which we define Dirichlet boundary conditions
        m_wallBoundarySides.push_back(boundary::west);
        m_wallBoundarySides.push_back(boundary::south);
        m_wallBoundarySides.push_back(boundary::north);
        m_wallBoundarySides.push_back(boundary::north);
        m_wallBoundarySides.push_back(boundary::south);
        //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
        m_wallBoundaryPatches << 0, 0, 1, 2, 2;
        //inlet sides at which we define Dirichlet boundary conditions
        m_inletBoundarySides.push_back(boundary::west);
        //patches at which is the inlet boundary condition set (must have the same length and ordering as the vector of inletBoundarySides!!)
        m_inletBoundaryPatches << 2;
    }

    T computeWallDistance(uwbINSSolverSteady<T>& navStokes)
    {
        T inletWidth = m_height / 2;
        T Re = m_uFreeStream * inletWidth / m_viscosity;

        //vector of the sides of the patches from which the wall distance is computed
        std::vector<boxSide> distanceSides = m_wallBoundarySides;

        //vector of indexes of the patches corresponding to distanceSides
        //length of the vector distancePatches must be equal to the length of vector distanceSides
        gsVector<int> distancePatches = m_wallBoundaryPatches;

        int numSamplePts = 50; //number of sample points for which the distance to the boundary is computed
        T maxYplus = 2.5; //maximum dimensionless wall distance which is accepted

        //table with wall distance information for every chosen side is printed, if the second from the end input parameter is set as true
        //estimation of the wall distance will be computed, if the last input parameter is set as true
        return navStokes.computeDimensionlessWallDistance(distancePatches, distanceSides, m_viscosity, Re, m_uFreeStream, maxYplus, numSamplePts, true, true);
    }

    void defineBCs_NS(gsBoundaryConditions<T>& bcInfo)
    {
        /*        T coef = (-1) * m_uMax * 16. / math::pow(m_height, 2);
                T shift = 3. / 4. * m_height;
                std::string str_coef = std::to_string(coef);
                std::string str_shift = std::to_string(shift);
                std::string inletVelocity = str_coef + " * (y - " + str_shift + ")^2 + " + std::to_string(m_uMax);//"-4*(y-1.5)^2 + 1";//
        */

                T inletWidth = m_height / 2;
                T Re = m_uFreeStream * inletWidth / m_viscosity;
                T n = 1.03 * math::log(Re) - 3.6;
                std::string str_n = std::to_string(n);

                T shift = 3. / 4. * m_height;
                T R = m_height/2;
                std::string str_R = std::to_string(R);
                std::string str_shift = std::to_string(shift);
                std::string inletVelocity = std::to_string(m_uMax) + " * ( 1 - (abs(y - " + str_shift + ") /" + str_R + ") ^" + str_n + ")";

        gsFunctionExpr<T> Uin(inletVelocity, "0", 2); // inlet velocity
        gsFunctionExpr<T> Uwall("0", "0", 2); // wall velocity

        Base::defineBCs_NS(bcInfo, Uin, Uwall);
    }

    void solvePoissonEquation(uwbTMSolverKOmega<T>& turbSolver, bool plotMeshes = true)
    {
        int numRefinePoisson = 6;

        gsMultiBasis<T> tbasisPoisson(m_patches); // basis for RANS equations
        for (int i = 0; i < numRefinePoisson; ++i)
            tbasisPoisson.uniformRefine();
        gsMatrix<T> box_u0Poisson(2, 2);
        box_u0Poisson << 0, 1, 0, 0;
        int numBoxRef = 2;
        for (int i = 0; i < numBoxRef; i++)
        {
            tbasisPoisson.refine(0, box_u0Poisson);
            tbasisPoisson.refine(1, box_u0Poisson);
        }
        if (plotMeshes)
        {
            gsMesh<T> mesh;
            makeMesh(tbasisPoisson.at(0), mesh, 10);
            m_patches.patch(0).evaluateMesh(mesh);
            gsWriteParaview(mesh, m_sGeometry+"meshPoissonPatch0");
            gsMesh<T> mesh1;
            makeMesh(tbasisPoisson.at(1), mesh1, 10);
            m_patches.patch(1).evaluateMesh(mesh1);
            gsWriteParaview(mesh1, m_sGeometry+"meshPoissonPatch1");
            gsMesh<T> mesh2;
            makeMesh(tbasisPoisson.at(2), mesh2, 10);
            m_patches.patch(2).evaluateMesh(mesh2);
            gsWriteParaview(mesh2, m_sGeometry+"meshPoissonPatch2");
        }

        gsFunctionExpr<T> fw("1", 2);
        gsFunctionExpr<T> gw("0", 2);
        gsFunctionExpr<T> wallw("0.0", 2);
        gsBoundaryConditions<T> bcInfow;
        bcInfow.addCondition(0, boundary::east, condition_type::neumann, gw, 0);
        bcInfow.addCondition(1, boundary::east, condition_type::neumann, gw, 0);
        bcInfow.addCondition(2, boundary::west, condition_type::neumann, gw, 0);
        bcInfow.addCondition(0, boundary::west, condition_type::dirichlet, wallw, 0);
        bcInfow.addCondition(0, boundary::south, condition_type::dirichlet, wallw, 0);
        bcInfow.addCondition(1, boundary::north, condition_type::dirichlet, wallw, 0);
        bcInfow.addCondition(2, boundary::north, condition_type::dirichlet, wallw, 0);
        bcInfow.addCondition(2, boundary::south, condition_type::dirichlet, wallw, 0);

        gsInfo << "\nSolving Poisson equation.\n";

        turbSolver.setPoissonSolution(m_patches, tbasisPoisson, bcInfow, fw, true, m_plot_pts);

        gsInfo << "Poisson equation resolved.\n\n";
        turbSolver.plotWallDistance(m_sGeometry+"wall_distance", m_plot_pts);
    }

    gsMultiPatch<T> makeMultiPatch()
    {
        uwbGeometryCreators<T> geometryCreator;
        m_patches = geometryCreator.BSplineBackwardStep2D(m_length, m_height, m_inletLength);

        return m_patches;
    }

    void refineBasis(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal)
    {
        gsMatrix<> box(2, 2);
        box << 0, 1, 0, 0;

        for (int i = 0; i < numRefine; ++i)
            basis.uniformRefine();

        int numURef = 3;
        for (int i = 0; i < numURef; i++)
        {
            basis.refine(0, box);
            basis.refine(1, box);
            //if (i < 1)
            //    basis.refine(2, box);
        }

        // refinement near wal
        int numRefineLocalWal = 1;

        real_t parArea = 0.1;

        for (int j = 0; j < numRefineLocal; j++)
        {
            box << 0, 0, 0, parArea;
            for (int i = 0; i < numRefineLocalWal; i++)
            {
                basis.refine(0, box);
                basis.refine(1, box);
                basis.refine(2, box);
            }

            box << 0, 0, 1 - parArea, 1;
            for (int i = 0; i < numRefineLocalWal; i++)
            {
                basis.refine(0, box);
                basis.refine(1, box);
                basis.refine(2, box);
            }

            box << 0, parArea/math::pow(2, numURef), 0, 0;
            for (int i = 0; i < numRefineLocalWal; i++)
            {
                basis.refine(0, box);
                basis.refine(1, box);
            }

            box << 1 - parArea, 1, 0, 0;
            for (int i = 0; i < numRefineLocalWal; i++)
                basis.refine(2, box);

            parArea = parArea / 2;
        }
    }

public:
    const T getFreeStream() const { return m_uFreeStream; }

protected:
    T m_length;
    T m_height;
    T m_inletLength;
    T m_uMax;
    T m_uFreeStream;

    // members from uwbRANSExampleSetting
    using Base::m_viscosity;
    using Base::m_plot_pts;
    using Base::m_patches;
    using Base::m_wallBoundarySides;
    using Base::m_inletBoundarySides;
    using Base::m_wallBoundaryPatches;
    using Base::m_inletBoundaryPatches;
    using Base::m_sGeometry;

}; //uwbRANSBackwardStepExample

//========================================================================================================================================
//========================================================================================================================================
//========================================================================================================================================

template<class T>
class uwbRANSBackwardStep2D4PExample : public uwbRANSExampleSetting<T>
{
public:
    typedef uwbRANSExampleSetting<T> Base;

public:
    uwbRANSBackwardStep2D4PExample(T viscosity, T stepHeight, T preStepHeight, T preStepLength,
                                   T prePatchLength, T length, int deg, T uMax, int plot_pts, bool symmetry = false) :
        Base(viscosity, plot_pts), m_stepHeight(stepHeight), m_preStepHeight(preStepHeight),
        m_preStepLength(preStepLength), m_prePatchLength(prePatchLength),
        m_length(length), m_deg(deg), m_uMax(uMax), m_bSymmetry(symmetry)
    {
        initMembers();
        setBoundaries();
    }

    ~uwbRANSBackwardStep2D4PExample() { }

protected:
    void initMembers()
    {
        m_wallBoundaryPatches.setZero(5);
        m_inletBoundaryPatches.setZero(1);
        m_uFreeStream = m_uMax;
        m_sGeometry = "step4P_";
    }

public:
    void setBoundaries()
    {
        //walls at which we define Dirichlet boundary conditions
        m_wallBoundarySides.push_back(boundary::west);
        m_wallBoundarySides.push_back(boundary::south);
        m_wallBoundarySides.push_back(boundary::north);
        m_wallBoundarySides.push_back(boundary::north);
        m_wallBoundarySides.push_back(boundary::south);
        //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
        m_wallBoundaryPatches << 0, 0, 1, 2, 2;
        //inlet sides at which we define Dirichlet boundary conditions
        m_inletBoundarySides.push_back(boundary::west);
        //patches at which is the inlet boundary condition set (must have the same length and ordering as the vector of inletBoundarySides!!)
        if (m_bSymmetry)
            m_inletBoundaryPatches << 3;
        else
            m_inletBoundaryPatches << 2;
    }

    T computeWallDistance(uwbINSSolverSteady<T>& navStokes)
    {
        T inletWidth = m_preStepHeight;
        T Re = m_uFreeStream * inletWidth / m_viscosity;

        //vector of the sides of the patches from which the wall distance is computed
        std::vector<boxSide> distanceSides = m_wallBoundarySides;

        //vector of indexes of the patches corresponding to distanceSides
        //length of the vector distancePatches must be equal to the length of vector distanceSides
        gsVector<int> distancePatches = m_wallBoundaryPatches;

        int numSamplePts = 50; //number of sample points for which the distance to the boundary is computed
        T maxYplus = 2.5; //maximum dimensionless wall distance which is accepted

        //table with wall distance information for every chosen side is printed, if the second from the end input parameter is set as true
        //estimation of the wall distance will be computed, if the last input parameter is set as true
        return navStokes.computeDimensionlessWallDistance(distancePatches, distanceSides, m_viscosity, Re, m_uFreeStream, maxYplus, numSamplePts, true, true);
    }

    void defineBCs_NS(gsBoundaryConditions<T>& bcInfo)
    {
/*        T coef = (-1) * m_uMax * 16. / math::pow(m_height, 2);
        T shift = 3. / 4. * m_height;
        std::string str_coef = std::to_string(coef);
        std::string str_shift = std::to_string(shift);
        std::string inletVelocity = str_coef + " * (y - " + str_shift + ")^2 + " + std::to_string(m_uMax);//"-4*(y-1.5)^2 + 1";//
*/
        std::string inletVelocity;
        if (m_bSymmetry)
            inletVelocity = std::to_string(m_uMax);
        else
        {
            T inletWidth = m_preStepHeight;
            T Re = m_uFreeStream * inletWidth / m_viscosity;
            Re = Re / 100.;
            T n = 1.03 * math::log(Re) - 3.6;
            std::string str_n = std::to_string(n);

            T shift = (1. / 2. * m_preStepHeight) + m_stepHeight;
            T R = m_preStepHeight;
            std::string str_R = std::to_string(R);
            std::string str_shift = std::to_string(shift);
            inletVelocity = std::to_string(m_uMax) + " * ( 1 - (abs(y - " + str_shift + ") /" + str_R + ") ^" + str_n + ")";

/*            T shift = (1. / 2. * m_preStepHeight) + m_stepHeight;
            T coef = (-1) * m_uMax / math::pow(m_preStepHeight + m_stepHeight - shift, 2);
            std::string str_coef = std::to_string(coef);
            std::string str_shift = std::to_string(shift);
            inletVelocity = str_coef + " * (y - " + str_shift + ")^2 + " + std::to_string(m_uMax);//"-4*(y-1.5)^2 + 1";//
*/


        }

        gsFunctionExpr<T> Uin(inletVelocity, "0", 2); // inlet velocity
        gsFunctionExpr<T> Uwall("0", "0", 2); // wall velocity

        Base::defineBCs_NS(bcInfo, Uin, Uwall);
    }

    void defineBCs_TM(gsBoundaryConditions<T>& bcInfoTurb, T kIn, T kWall, T oIn, T oWall, bool bOFbc = false)
        {
        gsInfo << "OF BC = " << bOFbc << "\n";
        /*if (bOFbc)
            Base::defineBCs_TM(bcInfoTurb, kIn, kWall, oIn, oWall);
        else
        {*/
        std::string inletK, inletO;
        /*T shift = (1. / 2. * m_preStepHeight) + m_stepHeight;
        T coefK = (-1) * kIn / math::pow(m_preStepHeight + m_stepHeight - shift, 2);
        T coefO = (oWall - oIn) / math::pow(shift, 2);
        std::string str_coefK = std::to_string(coefK);
        std::string str_coefO = std::to_string(coefO);
        std::string str_shift = std::to_string(shift);
        inletK = str_coefK + " * (y - " + str_shift + ")^2 + " + std::to_string(kIn);//"-4*(y-1.5)^2 + 1";
        inletO = str_coefO + " * (y - " + str_shift + ")^2 + " + std::to_string(oIn);*/


        T inletWidth = m_preStepHeight;
        T Re = m_uFreeStream * inletWidth / m_viscosity;
        Re = Re / 100.;
        T n = 1.03 * math::log(Re) - 3.6;
        std::string str_n = std::to_string(n);

        T shift = (1. / 2. * m_preStepHeight) + m_stepHeight;
        T R = m_preStepHeight;
        std::string str_R = std::to_string(R);
        std::string str_shift = std::to_string(shift);
        inletK = std::to_string(kIn) + " * ( 1 - (abs(y - " + str_shift + ") /" + str_R + ") ^" + str_n + ")";
        inletO = "(-1) * " + std::to_string(oWall - oIn) + " * ( 1 - (abs(y - " + str_shift + ") /" + str_R + ") ^" + str_n + ") + " + std::to_string(oWall);

            // Boundary conditions
            //gsFunctionExpr<T> Kin(util::to_string(kIn), 2);
            //gsFunctionExpr<T> Oin(util::to_string(oIn), 2);

            gsFunctionExpr<T> Kin(inletK, 2);
            gsFunctionExpr<T> Oin(inletO, 2);
            gsFunctionExpr<T> KWall(util::to_string(kWall), 2);
            gsFunctionExpr<T> OWall(util::to_string(oWall), 2);

            for (int i = 0; i < m_inletBoundaryPatches.rows(); i++)
            {
                bcInfoTurb.addCondition(m_inletBoundaryPatches(i), m_inletBoundarySides[i], condition_type::dirichlet, Kin, 0);
                bcInfoTurb.addCondition(m_inletBoundaryPatches(i), m_inletBoundarySides[i], condition_type::dirichlet, Oin, 1);
            }
            for (int i = 0; i < m_wallBoundaryPatches.rows(); i++)
            {
                bcInfoTurb.addCondition(m_wallBoundaryPatches(i), m_wallBoundarySides[i], condition_type::dirichlet, KWall, 0);
                bcInfoTurb.addCondition(m_wallBoundaryPatches(i), m_wallBoundarySides[i], condition_type::dirichlet, OWall, 1);
            }
        //}
    }

    void solvePoissonEquation(uwbTMSolverKOmega<T>& turbSolver, bool plotMeshes = true)
    {
        int numRefinePoisson = 6;

        gsMultiBasis<T> tbasisPoisson(m_patches); // basis for RANS equations
        for (int i = 0; i < numRefinePoisson; ++i)
            tbasisPoisson.uniformRefine();
        gsMatrix<T> box_u0Poisson(2, 2);
        box_u0Poisson << 0, 1, 0, 0;
        int numBoxRef = 2;
        for (int i = 0; i < numBoxRef; i++)
        {
            tbasisPoisson.refine(0, box_u0Poisson);
            tbasisPoisson.refine(1, box_u0Poisson);
        }
        if (plotMeshes)
        {
            gsMesh<T> mesh;
            makeMesh(tbasisPoisson.at(0), mesh, 10);
            m_patches.patch(0).evaluateMesh(mesh);
            gsWriteParaview(mesh, m_sGeometry+"meshPoissonPatch0");
            gsMesh<T> mesh1;
            makeMesh(tbasisPoisson.at(1), mesh1, 10);
            m_patches.patch(1).evaluateMesh(mesh1);
            gsWriteParaview(mesh1, m_sGeometry+"meshPoissonPatch1");
            gsMesh<T> mesh2;
            makeMesh(tbasisPoisson.at(2), mesh2, 10);
            m_patches.patch(2).evaluateMesh(mesh2);
            gsWriteParaview(mesh2, m_sGeometry+"meshPoissonPatch2");
        }

        gsFunctionExpr<T> fw("1", 2);
        gsFunctionExpr<T> gw("0", 2);
        gsFunctionExpr<T> wallw("0.0", 2);
        gsBoundaryConditions<T> bcInfow;
        bcInfow.addCondition(0, boundary::east, condition_type::neumann, gw, 0);
        bcInfow.addCondition(1, boundary::east, condition_type::neumann, gw, 0);
        bcInfow.addCondition(2, boundary::west, condition_type::neumann, gw, 0);
        bcInfow.addCondition(0, boundary::west, condition_type::dirichlet, wallw, 0);
        bcInfow.addCondition(0, boundary::south, condition_type::dirichlet, wallw, 0);
        bcInfow.addCondition(1, boundary::north, condition_type::dirichlet, wallw, 0);
        bcInfow.addCondition(2, boundary::north, condition_type::dirichlet, wallw, 0);
        bcInfow.addCondition(2, boundary::south, condition_type::dirichlet, wallw, 0);

        gsInfo << "\nSolving Poisson equation.\n";

        turbSolver.setPoissonSolution(m_patches, tbasisPoisson, bcInfow, fw, true, m_plot_pts);

        gsInfo << "Poisson equation resolved.\n\n";
        turbSolver.plotWallDistance(m_sGeometry+"wall_distance", m_plot_pts);
    }

    gsMultiPatch<T> makeMultiPatch()
    {
        uwbGeometryCreators<T> geometryCreator;
        m_patches = geometryCreator.BSplineBackwardStep2D_4patches(m_deg, m_stepHeight, m_preStepHeight,
                                                                   m_preStepLength, m_prePatchLength, m_length);

        return m_patches;
    }

    void refineBasis(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal)
    {
        //gsInfo << "info1\n";

        //numRefine = 3;
        //numRefineLocal = 7;

        gsMatrix<> box(2, 2);
        box << 0, 1, 0, 0;
        gsMatrix<> box2(2, 2);
        box2 << 0, 0, 0, 1;

        //basis.uniformRefine();

        if (m_bSymmetry)
        {
            for (int p = 1; p < 4; p ++)
                basis.refine(p, box2);

            for (int i = 0; i < numRefine; ++i)
            {
                for (int p = 0; p < 3; p ++)
                {
                    basis.refine(p, box);
                    basis.refine(p, box2);
                }
                basis.refine(3, box2);
            }

            basis.refine(2, box);
            basis.refine(3, box);
            basis.refine(3, box);
            box2 << 0.5, 1, 0, 0;
            basis.refine(2, box2);
            box2 << 0.5, 1, 0, 0;
            basis.refine(3, box2);
        }
        else
        {
            box << 0, 1, 0, 0;
            box2 << 0, 0, 0, 1;

            for (int p = 1; p < 3; p ++)
                basis.refine(p, box2);

            for (int i = 0; i < numRefine; ++i)
            {
                for (int p = 0; p < 3; p ++)
                {
                    basis.refine(p, box);
                }
                basis.refine(0, box2);
            }
            for (int i = 0; i < numRefine-1; ++i)
                for (int p = 1; p < 3; p ++)
                    basis.refine(p, box2);

            box2 << 0, 0, 0, 0.25;
            for (int p = 1; p < 3; p ++)
                basis.refine(p, box2);
            box2 << 0, 0, 0.75, 1;
            for (int p = 1; p < 3; p ++)
                basis.refine(p, box2);


            basis.refine(2, box);
            box2 << 0.5, 1, 0, 0;
            basis.refine(2, box2);
        }

        //gsInfo << "info2\n";

            //basis.uniformRefine();

        int numURef = 1;//3;
        for (int i = 0; i < numURef; i++)
        {
            basis.refine(0, box);
            basis.refine(1, box);
            //if (i < 1)
            //    basis.refine(2, box);
        }
        //basis.refine(3, box);

        // refinement near wal
        int numRefineLocalWal = 1;

        int numPatches;
        if (m_bSymmetry)
            numPatches = 4;
        else
            numPatches = 3;

        real_t parArea = 0.1;//0.05;//0.1;//

        //for (int j = 0; j < numRefineLocal-3; j++)
        for (int j = 0; j < numRefineLocal-5; j++)
        {
            box << 0, 0, 0, parArea;
            for (int i = 0; i < numRefineLocalWal; i++)
                for (int p = 1; p < numPatches; p ++)
                    basis.refine(p, box);

            box << 0, 0, 1 - parArea, 1;
            for (int i = 0; i < numRefineLocalWal; i++)
                for (int p = 1; p < numPatches; p ++)
                    basis.refine(p, box);

            parArea = parArea / 2;
        }
        parArea = 0.1;//0.2;
        //for (int j = 0; j < numRefineLocal-3; j++)
        for (int j = 0; j < numRefineLocal-5; j++)
        {
            box << 0, 0, 0, parArea;
            for (int i = 0; i < numRefineLocalWal; i++)
                    basis.refine(0, box);

            parArea = parArea / 2;
        }
        parArea = 0.1;//0.2;
        for (int j = 0; j < numRefineLocal-5; j++)
        {
            box << 0, 0, 1 - parArea, 1;
            for (int i = 0; i < numRefineLocalWal; i++)
                    basis.refine(0, box);

            parArea = parArea / 2;
        }

        //refineBasisLocalZones(basis, patch, int direction, int position, int knot, int num_of_refine);
        for (int j = 0; j < 2; j++)
            Base::refineBasisLocalZones(basis, 0, 1, 0, 1, 1); //south
        Base::refineBasisLocalZones(basis, 0, 1, 1, 1, 1); //north
        for (int j = 0; j < 3; j++)
        {
            Base::refineBasisLocalZones(basis, 1, 1, 0, 1, 1);
            Base::refineBasisLocalZones(basis, 1, 1, 1, 1, 1);
            Base::refineBasisLocalZones(basis, 2, 1, 0, 1, 1);
            Base::refineBasisLocalZones(basis, 2, 1, 1, 1, 1);
        }

        //gsInfo << "info3\n";

        real_t parArea2 = 0.1;//0.2;
        for (int j = 0; j < numRefineLocal-3; j++)
        {
            box << 1 - parArea2, 1, 0, 0;
            for (int i = 0; i < numRefineLocalWal; i++)
                basis.refine(2, box);

            parArea2 = parArea2 / 2;
        }
        parArea2 = 0.1/8;//0.2;
        for (int j = 0; j < 2; j++)
        {
            box << 1 - parArea2, 1, 0, 0;
                basis.refine(2, box);

            parArea2 = parArea2 / 2;
        }

        box << 0.025, 0.5, 0, 0;
        basis.refine(2, box);

        parArea2 = 0.025;
        box << 0, parArea2, 0, 0;
        basis.refine(2, box);
        //refineBasisLocalZones(basis, patch, int direction, int position, int knot, int num_of_refine);
        for (int j = 0; j < 12; j++)
            Base::refineBasisLocalZones(basis, 2, 0, 0, 1, 1); //west

        for (int j = 1; j < numRefineLocal-2; j++)
        {

            box << 0, 1./math::pow(2, j), 0, 0;
            for (int i = 0; i < numRefineLocalWal; i++)
            {
                basis.refine(0, box);
                basis.refine(1, box);
            }
        }

        parArea2 = 0.03;
        for (int j = 0; j < 2; j++)
        {
            box << 0, parArea2, 0, 0;
            basis.refine(0, box);
            basis.refine(1, box);
            parArea2 = parArea2 / 2;
        }

        //gsInfo << "info4\n";

        //refineBasisLocalZones(basis, patch, int direction, int position, int knot, int num_of_refine);
        for (int j = 0; j < 2; j++)
            Base::refineBasisLocalZones(basis, 0, 1, 0, 1, 1); //south
        for (int j = 0; j < 3; j++)
        {
            Base::refineBasisLocalZones(basis, 1, 1, 0, 1, 1);
            Base::refineBasisLocalZones(basis, 1, 1, 1, 1, 1);
            Base::refineBasisLocalZones(basis, 2, 1, 0, 1, 1);
            Base::refineBasisLocalZones(basis, 2, 1, 1, 1, 1);
        }

        Base::refineBasisLocalZones(basis, 0, 0, 0, 1, 1);
        Base::refineBasisLocalZones(basis, 1, 0, 0, 1, 1);
        Base::refineBasisLocalZones(basis, 2, 0, 1, 1, 1);
    }

public:
    const T getFreeStream() const { return m_uFreeStream; }

protected:
    T m_stepHeight;
    T m_preStepHeight;
    T m_preStepLength;
    T m_prePatchLength;
    T m_length;
    int m_deg;
    T m_uMax;
    T m_uFreeStream;

    bool m_bSymmetry;

    // members from uwbRANSExampleSetting
    using Base::m_viscosity;
    using Base::m_plot_pts;
    using Base::m_patches;
    using Base::m_wallBoundarySides;
    using Base::m_inletBoundarySides;
    using Base::m_wallBoundaryPatches;
    using Base::m_inletBoundaryPatches;
    using Base::m_sGeometry;

}; //uwbRANSBackwardStep4PExample

//========================================================================================================================================
//========================================================================================================================================
//========================================================================================================================================

template<class T>
class uwbRANSBackwardStep2D1PatchExample : public uwbRANSExampleSetting<T>
{
public:
    typedef uwbRANSExampleSetting<T> Base;

public:
    uwbRANSBackwardStep2D1PatchExample(T viscosity, T length, T height, T inletLength, T uMax, int plot_pts) :
        Base(viscosity, plot_pts), m_length(length), m_height(height), m_inletLength(inletLength), m_uMax(uMax)
    {
        initMembers();
        setBoundaries();
    }

    ~uwbRANSBackwardStep2D1PatchExample() { }

protected:
    void initMembers()
    {
        m_wallBoundaryPatches.setZero(2);
        m_inletBoundaryPatches.setZero(1);
        m_uFreeStream = m_uMax;
        m_sGeometry = "step1Patch_";
    }

public:
    void setBoundaries()
    {
        //walls at which we define Dirichlet boundary conditions
        m_wallBoundarySides.push_back(boundary::south);
        m_wallBoundarySides.push_back(boundary::north);
        //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
        m_wallBoundaryPatches << 0, 0;
        //inlet sides at which we define Dirichlet boundary conditions
        m_inletBoundarySides.push_back(boundary::west);
        //patches at which is the inlet boundary condition set (must have the same length and ordering as the vector of inletBoundarySides!!)
        m_inletBoundaryPatches << 0;
    }

    T computeWallDistance(uwbINSSolverSteady<T>& navStokes)
    {
        T inletWidth = m_height / 2;
        T Re = m_uFreeStream * inletWidth / m_viscosity;

        //vector of the sides of the patches from which the wall distance is computed
        std::vector<boxSide> distanceSides = m_wallBoundarySides;

        //vector of indexes of the patches corresponding to distanceSides
        //length of the vector distancePatches must be equal to the length of vector distanceSides
        gsVector<int> distancePatches = m_wallBoundaryPatches;

        int numSamplePts = 50; //number of sample points for which the distance to the boundary is computed
        T maxYplus = 2.5; //maximum dimensionless wall distance which is accepted

        //table with wall distance information for every chosen side is printed, if the second from the end input parameter is set as true
        //estimation of the wall distance will be computed, if the last input parameter is set as true
        return navStokes.computeDimensionlessWallDistance(distancePatches, distanceSides, m_viscosity, Re, m_uFreeStream, maxYplus, numSamplePts, true, true);
    }

    void defineBCs_NS(gsBoundaryConditions<T>& bcInfo)
    {
        T coef = (-1) * m_uMax * 16. / math::pow(m_height, 2);
        T shift = 3. / 4. * m_height;
        std::string str_coef = std::to_string(coef);
        std::string str_shift = std::to_string(shift);
        std::string inletVelocity = str_coef + " * (y - " + str_shift + ")^2 + " + std::to_string(m_uMax);//"-4*(y-1.5)^2 + 1";//

        gsFunctionExpr<T> Uin(inletVelocity, "0", 2); // inlet velocity
        gsFunctionExpr<T> Uwall("0", "0", 2); // wall velocity

        Base::defineBCs_NS(bcInfo, Uin, Uwall);
    }

    void solvePoissonEquation(uwbTMSolverKOmega<T>& turbSolver, bool plotMeshes = true)
    {
        int numRefinePoisson = 6;

        gsMultiBasis<T> tbasisPoisson(m_patches); // basis for RANS equations
        for (int i = 0; i < numRefinePoisson; ++i)
            tbasisPoisson.uniformRefine();
        /*gsMatrix<T> box_u0Poisson(2, 2);
        box_u0Poisson << 0, 1, 0, 0;
        int numBoxRef = 2;
        for (int i = 0; i < numBoxRef; i++)
        {
            tbasisPoisson.refine(0, box_u0Poisson);
            tbasisPoisson.refine(1, box_u0Poisson);
        }*/
        if (plotMeshes)
        {
            gsMesh<T> mesh;
            makeMesh(tbasisPoisson.at(0), mesh, 10);
            m_patches.patch(0).evaluateMesh(mesh);
            gsWriteParaview(mesh, m_sGeometry+"meshPoissonPatch0");
        }

        gsFunctionExpr<T> fw("1", 2);
        gsFunctionExpr<T> gw("0", 2);
        gsFunctionExpr<T> wallw("0.0", 2);
        gsBoundaryConditions<T> bcInfow;
        bcInfow.addCondition(0, boundary::east, condition_type::neumann, gw, 0);
        bcInfow.addCondition(0, boundary::west, condition_type::neumann, gw, 0);
        bcInfow.addCondition(0, boundary::south, condition_type::dirichlet, wallw, 0);
        bcInfow.addCondition(0, boundary::north, condition_type::dirichlet, wallw, 0);

        gsInfo << "\nSolving Poisson equation.\n";

        turbSolver.setPoissonSolution(m_patches, tbasisPoisson, bcInfow, fw, true, m_plot_pts);

        gsInfo << "Poisson equation resolved.\n\n";
        turbSolver.plotWallDistance(m_sGeometry+"wall_distance", m_plot_pts);
    }

    gsMultiPatch<T> makeMultiPatch()
    {
        uwbGeometryCreators<T> geometryCreator;
        m_patches = geometryCreator.BSplineBackwardStep2D1Patch();//(m_length, m_height, m_inletLength);

        return m_patches;
    }

    void refineBasis(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal)
    {
        gsMatrix<> box(2, 2);
        box << 0, 1, 0, 0;

        for (int i = 0; i < numRefine; ++i)
            basis.uniformRefine();

        box << 0, 0, 0, 1;
        basis.refine(0, box);

        int numURef = 1;
        box << 0.35, 1, 0, 0;
        for (int i = 0; i < numURef; i++)
            basis.refine(0, box);
        numURef = 1;
        box << 0.66, 1, 0, 0;
        for (int i = 0; i < numURef; i++)
            basis.refine(0, box);

        // refinement near wal
        int numRefineLocalWal = 1;

        real_t parArea = 0.04;

        for (int j = 0; j < numRefineLocal; j++)
        {
            box << 0, 0, 0, parArea;
            for (int i = 0; i < numRefineLocalWal; i++)
            {
                basis.refine(0, box);
            }

            box << 0, 0, 1 - parArea, 1;
            for (int i = 0; i < numRefineLocalWal; i++)
            {
                basis.refine(0, box);
            }

            /*box << 0, parArea/math::pow(2, numURef), 0, 0;
            for (int i = 0; i < numRefineLocalWal; i++)
            {
                basis.refine(0, box);
            }

            box << 1 - parArea, 1, 0, 0;
            for (int i = 0; i < numRefineLocalWal; i++)
                basis.refine(2, box);*/

            parArea = parArea / 2;
        }
    }

public:
    const T getFreeStream() const { return m_uFreeStream; }

protected:
    T m_length;
    T m_height;
    T m_inletLength;
    T m_uMax;
    T m_uFreeStream;

    // members from uwbRANSExampleSetting
    using Base::m_viscosity;
    using Base::m_plot_pts;
    using Base::m_patches;
    using Base::m_wallBoundarySides;
    using Base::m_inletBoundarySides;
    using Base::m_wallBoundaryPatches;
    using Base::m_inletBoundaryPatches;
    using Base::m_sGeometry;

}; //uwbRANSBackwardStep2D1PatchExample

//========================================================================================================================================
//========================================================================================================================================
//========================================================================================================================================

template<class T>
class uwbRANSCircle2D4PatchesExample : public uwbRANSExampleSetting<T>
{
public:
    typedef uwbRANSExampleSetting<T> Base;

public:
    uwbRANSCircle2D4PatchesExample(T viscosity,
                                   const T length, T const width, T const widthExpand, T const radius,
                                   T const centreX, T const centreY, bool const prepatch, T const prepatchWidth,
                                   T uMax, int plot_pts, bool periodic = true) :
        Base(viscosity, plot_pts),
        m_length(length), m_width(width), m_widthExpand(widthExpand), m_radius(radius), m_centreX(centreX),
        m_centreY(centreY), m_prepatch(prepatch), m_prepatchWidth(prepatchWidth),
        m_uMax(uMax), m_periodic(periodic)
    {
        initMembers();
        setBoundaries();
    }

    ~uwbRANSCircle2D4PatchesExample() { }

protected:
    void initMembers()
    {
        if (m_periodic)
            m_wallBoundaryPatches.setZero(4);
        else
        {
            if (m_prepatch)
                m_wallBoundaryPatches.setZero(10);
            else
                m_wallBoundaryPatches.setZero(8);
        }
        m_circleBoundaryPatches.setZero(4);
        m_inletBoundaryPatches.setZero(1);
        m_uFreeStream = m_uMax;
        m_sGeometry = "circle4P_";
    }

public:
    void setBoundaries()
    {
        m_circleBoundarySides.push_back(boundary::east);
        m_circleBoundarySides.push_back(boundary::east);
        m_circleBoundarySides.push_back(boundary::east);
        m_circleBoundarySides.push_back(boundary::east);
        //walls at which we define Dirichlet boundary conditions
        if (m_periodic)
            m_wallBoundarySides = m_circleBoundarySides;
        else
        {
            m_wallBoundarySides.push_back(boundary::east);
            m_wallBoundarySides.push_back(boundary::west);
            m_wallBoundarySides.push_back(boundary::east);
            m_wallBoundarySides.push_back(boundary::west);
            m_wallBoundarySides.push_back(boundary::east);
            m_wallBoundarySides.push_back(boundary::east);
            m_wallBoundarySides.push_back(boundary::north);
            m_wallBoundarySides.push_back(boundary::south);
            if (m_prepatch)
            {
                m_wallBoundarySides.push_back(boundary::north);
                m_wallBoundarySides.push_back(boundary::south);
            }
        }

        m_circleBoundaryPatches << 0, 1, 2, 3;
        if (m_periodic)
            m_wallBoundaryPatches = m_circleBoundaryPatches;
        else
        {
            //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
            if (m_prepatch)
                m_wallBoundaryPatches << 0, 1, 1, 2, 2, 3, 4, 4, 5, 5;
            else
                m_wallBoundaryPatches << 0, 1, 1, 2, 2, 3, 4, 4;
        }

        //inlet sides at which we define Dirichlet boundary conditions
        m_inletBoundarySides.push_back(boundary::west);
        //patches at which is the inlet boundary condition set (must have the same length and ordering as the vector of inletBoundarySides!!)
        if (m_prepatch)
            m_inletBoundaryPatches << 5;
        else
            m_inletBoundaryPatches << 0;
    }

    T computeWallDistance(uwbINSSolverSteady<T>& navStokes)
    {
        T Re = m_uFreeStream * (2*m_radius) / m_viscosity;

        //vector of the sides of the patches from which the wall distance is computed
        std::vector<boxSide> distanceSides = m_wallBoundarySides;

        //vector of indexes of the patches corresponding to distanceSides
        //length of the vector distancePatches must be equal to the length of vector distanceSides
        gsVector<int> distancePatches = m_wallBoundaryPatches;

        int numSamplePts = 50; //number of sample points for which the distance to the boundary is computed
        T maxYplus = 2.5; //maximum dimensionless wall distance which is accepted

        //table with wall distance information for every chosen side is printed, if the second from the end input parameter is set as true
        //estimation of the wall distance will be computed, if the last input parameter is set as true
        return navStokes.computeDimensionlessWallDistance(distancePatches, distanceSides, m_viscosity, Re, m_uFreeStream, maxYplus, numSamplePts, true, true);
    }

    void defineBCs_NS(gsBoundaryConditions<T>& bcInfo)
    {
        std::string inletVelocity;

        if (m_periodic)
            inletVelocity = std::to_string(m_uMax);
        else
        {
            T coef = (-1) * m_uMax / math::pow(-m_width/2 - m_widthExpand/2, 2);
            T shift = m_widthExpand/2;
            std::string str_coef = std::to_string(coef);
            std::string str_shift = std::to_string(shift);
            inletVelocity = str_coef + " * (y - " + str_shift + ")^2 + " + std::to_string(m_uMax);



            /*T Re = m_uFreeStream * (2*m_radius) / m_viscosity;
            T n = 1.03 * math::log(Re) - 3.6;
            std::string str_n = std::to_string(n);

            T shift = m_widthExpand/2;
            T R = (m_widthExpand + m_width)/2;
            std::string str_R = std::to_string(R);
            std::string str_shift = std::to_string(shift);
            inletVelocity = std::to_string(m_uMax) + " * ( 1 - (abs(y - " + str_shift + ") /" + str_R + ") ^" + str_n + ")";*/

            //gsInfo << "Uinlet = " << inletVelocity << "\n";
            //gsInfo << m_uMax * ( 1 - math::pow((0.2 - shift) / R, n) );
            //getchar();
        }

        gsFunctionExpr<T> Uin(inletVelocity, "0", 2); // inlet velocity
        gsFunctionExpr<T> Uwall("0", "0", 2); // wall velocity

        Base::defineBCs_NS(bcInfo, Uin, Uwall);
    }

    void solvePoissonEquation(uwbTMSolverKOmega<T>& turbSolver, bool plotMeshes = true)
    {
        int numRefinePoisson = 6;

        gsMultiBasis<T> tbasisPoisson(m_patches); // basis for RANS equations
        for (int i = 0; i < numRefinePoisson; ++i)
            tbasisPoisson.uniformRefine();
        gsMatrix<T> box_u0Poisson(2, 2);
        box_u0Poisson << 0, 1, 0, 0;
        int numBoxRef = 2;
        for (int i = 0; i < numBoxRef; i++)
        {
            tbasisPoisson.refine(4, box_u0Poisson);
        }
        if (plotMeshes)
        {
            gsMesh<T> mesh;
            makeMesh(tbasisPoisson.at(0), mesh, 10);
            m_patches.patch(0).evaluateMesh(mesh);
            gsWriteParaview(mesh, m_sGeometry+"meshPoissonPatch0");
            gsMesh<T> mesh1;
            makeMesh(tbasisPoisson.at(1), mesh1, 10);
            m_patches.patch(1).evaluateMesh(mesh1);
            gsWriteParaview(mesh1, m_sGeometry+"meshPoissonPatch1");
            gsMesh<T> mesh2;
            makeMesh(tbasisPoisson.at(2), mesh2, 10);
            m_patches.patch(2).evaluateMesh(mesh2);
            gsWriteParaview(mesh2, m_sGeometry+"meshPoissonPatch2");
            gsMesh<T> mesh3;
            makeMesh(tbasisPoisson.at(3), mesh3, 10);
            m_patches.patch(3).evaluateMesh(mesh3);
            gsWriteParaview(mesh3, m_sGeometry+"meshPoissonPatch3");
            gsMesh<T> mesh4;
            makeMesh(tbasisPoisson.at(4), mesh4, 10);
            m_patches.patch(4).evaluateMesh(mesh4);
            gsWriteParaview(mesh4, m_sGeometry+"meshPoissonPatch4");
            if (m_prepatch)
            {
                gsMesh<T> mesh5;
                makeMesh(tbasisPoisson.at(5), mesh5, 10);
                m_patches.patch(5).evaluateMesh(mesh5);
                gsWriteParaview(mesh5, m_sGeometry+"meshPoissonPatch5");
            }
        }

        gsFunctionExpr<T> fw("1", 2);
        gsFunctionExpr<T> gw("0", 2);
        gsFunctionExpr<T> wallw("0.0", 2);
        gsBoundaryConditions<T> bcInfow;
        if (m_prepatch)
            bcInfow.addCondition(5, boundary::west, condition_type::neumann, gw, 0);
        else
            bcInfow.addCondition(0, boundary::west, condition_type::neumann, gw, 0);
        bcInfow.addCondition(4, boundary::east, condition_type::neumann, gw, 0);
        if (m_prepatch)
        {
            if (m_periodic)
            {
                bcInfow.addCondition(5, boundary::north, condition_type::neumann, gw, 0);
                bcInfow.addCondition(5, boundary::south, condition_type::neumann, gw, 0);
            }
            else
            {
                bcInfow.addCondition(5, boundary::north, condition_type::dirichlet, wallw, 0);
                bcInfow.addCondition(5, boundary::south, condition_type::dirichlet, wallw, 0);
            }
        }
        bcInfow.addCondition(0, boundary::east, condition_type::dirichlet, wallw, 0);
        bcInfow.addCondition(1, boundary::east, condition_type::dirichlet, wallw, 0);
        bcInfow.addCondition(2, boundary::east, condition_type::dirichlet, wallw, 0);
        bcInfow.addCondition(3, boundary::east, condition_type::dirichlet, wallw, 0);
        if (m_periodic)
        {
            bcInfow.addCondition(1, boundary::west, condition_type::neumann, gw, 0);
            bcInfow.addCondition(2, boundary::west, condition_type::neumann, gw, 0);
            bcInfow.addCondition(4, boundary::north, condition_type::neumann, gw, 0);
            bcInfow.addCondition(4, boundary::south, condition_type::neumann, gw, 0);
        }
        else
        {
            bcInfow.addCondition(1, boundary::west, condition_type::dirichlet, wallw, 0);
            bcInfow.addCondition(2, boundary::west, condition_type::dirichlet, wallw, 0);
            bcInfow.addCondition(4, boundary::north, condition_type::dirichlet, wallw, 0);
            bcInfow.addCondition(4, boundary::south, condition_type::dirichlet, wallw, 0);
        }

        gsInfo << "\nSolving Poisson equation.\n";

        turbSolver.setPoissonSolution(m_patches, tbasisPoisson, bcInfow, fw, true, m_plot_pts);

        gsInfo << "Poisson equation resolved.\n\n";
        turbSolver.plotWallDistance(m_sGeometry+"wall_distance", m_plot_pts);
    }

    gsMultiPatch<T> makeMultiPatch()
    {
        uwbGeometryCreators<T> geometryCreator;
        m_patches = geometryCreator.BSplineCircle2D4Patches(m_length, m_width, m_widthExpand, m_radius, m_centreX, m_centreY,
                                                            m_prepatch, m_prepatchWidth, m_periodic);

        return m_patches;
    }

    void refineBasis(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal, T wallRefineKnot)
    {
        for (int i = 0; i < numRefine; ++i)
            basis.uniformRefine();

        if (m_periodic)
        {
            gsMatrix<> box_v0(2, 2);
            box_v0 << 0, 0, 0, 1;
            for (size_t k = 0; k < m_patches.nPatches(); k++)
                basis.refine(k, box_v0);

            /*box_v0 << 0, 0, 0.25, 0.75;
            basis.refine(3, box_v0);
            basis.refine(4, box_v0);*/

            /*gsMatrix<> box_v1(2, 2);
            box_v1 << 0, 0, 0.4, 0.6;
            basis.refine(3, box_v1);
            basis.refine(4, box_v1);*/

            gsMatrix<> box_u0(2, 2);
            for (int k = 0; k < 4; k++){
                for (int i = 0; i < numRefineLocal; i++)
                {
                    box_u0 << 1 - wallRefineKnot / math::pow(2, i), 1, 0, 0;
                    basis.refine(k, box_u0);
                }
            }

            /*gsMatrix<> box_u1(2, 2);
            box_u1 << 0, 1, 0, 0;
            for (size_t k = 0; k < m_patches.nPatches(); k++)
                basis.refine(k, box_u1);*/

        }
        else
        {
            gsMatrix<> box_v0(2, 2);
            //gsMatrix<> box_v1(2, 2);
            box_v0 << 0, 0, 0, 1;
            for (size_t k = 0; k < m_patches.nPatches(); k++)
                basis.refine(k, box_v0);

            /*for (size_t k = 0; k < m_patches.nPatches(); k++){
                for (int i = 0; i < numRefineLocal-4; i++)
                {
                    box_v0 << 0, 0, 0, wallRefineKnot / math::pow(2, i);
                    box_v1 << 0, 0, 1 - wallRefineKnot / math::pow(2, i), 1;
                    basis.refine(k, box_v0);
                    basis.refine(k, box_v1);
                }
            }

            gsMatrix<> box_u0(2, 2);
            gsMatrix<> box_u1(2, 2);
            for (int k = 0; k < 4; k++){
                for (int i = 0; i < numRefineLocal; i++)
                {
                    box_u0 << 0, wallRefineKnot / math::pow(2, i), 0, 0;
                    box_u1 << 1 - wallRefineKnot / math::pow(2, i), 1, 0, 0;
                    basis.refine(k, box_u0);
                    basis.refine(k, box_u1);
                }
            }*/
        }
    }

public:
    const T getFreeStream() const { return m_uFreeStream; }
    const gsVector<int>& getCircleBoundaryPatches() const { return m_circleBoundaryPatches; }
    const std::vector<boxSide> getCirleBoundarySides() const { return m_circleBoundarySides; }

protected:
    T m_length;
    T m_width;
    T m_widthExpand;
    T m_radius;
    T m_centreX;
    T m_centreY;
    bool m_prepatch;
    T m_prepatchWidth;
    T m_uMax;
    T m_uFreeStream;

    std::vector<boxSide> m_circleBoundarySides;
    gsVector<int> m_circleBoundaryPatches;

    bool m_periodic;

    // members from uwbRANSExampleSetting
    using Base::m_viscosity;
    using Base::m_plot_pts;
    using Base::m_patches;
    using Base::m_wallBoundarySides;
    using Base::m_inletBoundarySides;
    using Base::m_wallBoundaryPatches;
    using Base::m_inletBoundaryPatches;
    using Base::m_sGeometry;

}; //uwbRANSCircle2D4PatchesExample

//========================================================================================================================================
//========================================================================================================================================
//========================================================================================================================================

template<class T>
class uwbRANSProfileExample : public uwbRANSExampleSetting<T>
{
public:
    typedef uwbRANSExampleSetting<T> Base;

public:
    uwbRANSProfileExample(T viscosity, std::string strProfile, unsigned indexOfProfile, unsigned numBladeprofiles,
                          unsigned numBlades, int geomChoice, std::vector<T> kvfit_knots, bool coarse,
                          gsVector<T> geomParams, int plot_pts) :
        Base(viscosity, plot_pts), m_strProfile(strProfile), m_index_of_profile(indexOfProfile),
        m_numBladeprofiles(numBladeprofiles), m_numBlades(numBlades), m_geomChoice(geomChoice),
        m_kvfit_knots(kvfit_knots), m_coarse(coarse), m_geomParams(geomParams)
    {
        initMembers();
        setBoundaries();
    }

    ~uwbRANSProfileExample() { }

protected:
    void initMembers()
    {
        switch(m_geomChoice){
            case 1:
            {
                m_wallBoundaryPatches.setZero(2);
                m_inletBoundaryPatches.setZero(1);
                m_sGeometry = "profileDomain1_";
            }
            break;

        case 2:
        {
            m_wallBoundaryPatches.setZero(5);
            m_inletBoundaryPatches.setZero(3);
            m_sGeometry = "profileDomain2_";
        }
        break;

        case 3:
        {
            m_wallBoundaryPatches.setZero(6);
            m_inletBoundaryPatches.setZero(1);
            m_sGeometry = "profileDomain3_";
        }
        break;

        case 32:
        {
            m_wallBoundaryPatches.setZero(7);
            m_inletBoundaryPatches.setZero(1);
            m_sGeometry = "profileDomain3b_";
        }
        break;

        case 33:
        {
            m_wallBoundaryPatches.setZero(6);
            m_inletBoundaryPatches.setZero(1);
            m_sGeometry = "profileDomain3c_";
        }
        break;


            case 4:
            {
                m_wallBoundaryPatches.setZero(3);
                m_inletBoundaryPatches.setZero(1);
                m_sGeometry = "profileDomain4_";
            }
            break;

            case 5:
            {
                m_wallBoundaryPatches.setZero(2);
                m_inletBoundaryPatches.setZero(1);
                m_sGeometry = "profileDomain5_";
            }
            break;

            default: {
                m_wallBoundaryPatches.setZero(2);
                m_inletBoundaryPatches.setZero(1);
                m_sGeometry = "profileDomain0_";
            }
        }

        m_velocity_absolute_x.setZero(m_numBladeprofiles);
        m_velocity_absolute_y.setZero(m_numBladeprofiles);
        m_rr.setZero(m_numBladeprofiles);
        m_leading_angle.setZero(m_numBladeprofiles);
        m_angle.setZero(m_numBladeprofiles);
        m_velocity_blade.setZero(m_numBladeprofiles);    //angular velocity
        m_flow_rate = 0.;
        m_uFreeStream = 0.;
        m_setting_of_blade_rotation = 0;
        m_setting_of_guide_vanes = 0;

        m_UyInZero = false;
    }

public:
    void setParameters(int setting_of_blade_rotation, int setting_of_guide_vanes)
    {
        m_setting_of_blade_rotation = setting_of_blade_rotation; //0-optimal, 1-maximal, 2-minimal
        m_setting_of_guide_vanes = setting_of_guide_vanes; //0-compute with the formula, 1- 56° GV for optimal and 60° GV for maximal, 2-68° GV for optimal a 72° GV for maximal
    }

    void setUyInZero(bool isZero) { m_UyInZero = isZero; }

    void setBoundaries()
    {
        switch(m_geomChoice){
            case 1:
            {
                //walls at which we define Dirichlet boundary conditions
                m_wallBoundarySides.push_back(boundary::north); // sunction side
                m_wallBoundarySides.push_back(boundary::south); // pressure side
                //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
                m_wallBoundaryPatches << 1, 1;
                //inlet sides at which we define Dirichlet boundary conditions
                m_inletBoundarySides.push_back(boundary::west);
                //patches at which is the inlet boundary condition set (must have the same length and ordering as the vector of inletBoundarySides!!)
                m_inletBoundaryPatches << 0;
            }
            break;

        case 2:
        {
            //walls at which we define Dirichlet boundary conditions
            m_wallBoundarySides.push_back(boundary::west); // sunction side
            m_wallBoundarySides.push_back(boundary::south); // leading side
            m_wallBoundarySides.push_back(boundary::south); // sunction side
            m_wallBoundarySides.push_back(boundary::south); // leading side
            m_wallBoundarySides.push_back(boundary::south); // sunction side
            //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
            m_wallBoundaryPatches << 0,1,2,9,10;
            //inlet sides at which we define Dirichlet boundary conditions
            m_inletBoundarySides.push_back(boundary::east);
            m_inletBoundarySides.push_back(boundary::north);
            m_inletBoundarySides.push_back(boundary::north);
            //patches at which is the inlet boundary condition set (must have the same length and ordering as the vector of inletBoundarySides!!)
            m_inletBoundaryPatches << 4,5,6;
        }
        break;

        case 3:
        {
            //walls at which we define Dirichlet boundary conditions
            m_wallBoundarySides.push_back(boundary::north); // sunction side
            m_wallBoundarySides.push_back(boundary::north); // pressure side
            m_wallBoundarySides.push_back(boundary::north); // sunction side
            m_wallBoundarySides.push_back(boundary::south); // pressure side
            m_wallBoundarySides.push_back(boundary::south); // pressure side
            m_wallBoundarySides.push_back(boundary::south); // pressure side
            //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
            m_wallBoundaryPatches << 0, 2,4,3,5,7;
            //inlet sides at which we define Dirichlet boundary conditions
            m_inletBoundarySides.push_back(boundary::west);
            //patches at which is the inlet boundary condition set (must have the same length and ordering as the vector of inletBoundarySides!!)
            m_inletBoundaryPatches << 1;
        }
        break;

        case 32:
        {
            //walls at which we define Dirichlet boundary conditions
            m_wallBoundarySides.push_back(boundary::east); // sunction side
            m_wallBoundarySides.push_back(boundary::east); // pressure side
            m_wallBoundarySides.push_back(boundary::east); // sunction side
            //m_wallBoundarySides.push_back(boundary::east); // pressure side
            m_wallBoundarySides.push_back(boundary::north); // pressure side
            m_wallBoundarySides.push_back(boundary::north); // pressure side
            m_wallBoundarySides.push_back(boundary::north); // pressure side
            //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
            m_wallBoundaryPatches << 0, 2,4,3,5,7;
            //inlet sides at which we define Dirichlet boundary conditions
            m_inletBoundarySides.push_back(boundary::west);
            //patches at which is the inlet boundary condition set (must have the same length and ordering as the vector of inletBoundarySides!!)
            m_inletBoundaryPatches << 1;
        }
        break;

        case 33:
        {
            //walls at which we define Dirichlet boundary conditions
            m_wallBoundarySides.push_back(boundary::east); // sunction side
            m_wallBoundarySides.push_back(boundary::east); // pressure side
            m_wallBoundarySides.push_back(boundary::east); // sunction side
            m_wallBoundarySides.push_back(boundary::west); // pressure side
            m_wallBoundarySides.push_back(boundary::west); // pressure side
            m_wallBoundarySides.push_back(boundary::west); // pressure side
            //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
            m_wallBoundaryPatches << 0, 2,4,3,5,7;
            //inlet sides at which we define Dirichlet boundary conditions
            m_inletBoundarySides.push_back(boundary::south);
            //patches at which is the inlet boundary condition set (must have the same length and ordering as the vector of inletBoundarySides!!)
            m_inletBoundaryPatches << 9;
        }
        break;

            case 4:
            {
            //walls at which we define Dirichlet boundary conditions
            m_wallBoundarySides.push_back(boundary::south); // sunction side
           // m_wallBoundarySides.push_back(boundary::south); // pressure side
            m_wallBoundarySides.push_back(boundary::south); // sunction side
            m_wallBoundarySides.push_back(boundary::west); // pressure side
            //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
            m_wallBoundaryPatches << 0, 1, 3;
            //inlet sides at which we define Dirichlet boundary conditions
            m_inletBoundarySides.push_back(boundary::north);
            //patches at which is the inlet boundary condition set (must have the same length and ordering as the vector of inletBoundarySides!!)
            m_inletBoundaryPatches << 0;
            }
            break;

            case 5:
            {
                //walls at which we define Dirichlet boundary conditions
                m_wallBoundarySides.push_back(boundary::north); // sunction side
                m_wallBoundarySides.push_back(boundary::east); // pressure side
                //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
                m_wallBoundaryPatches << 0, 1;
                //inlet sides at which we define Dirichlet boundary conditions
                m_inletBoundarySides.push_back(boundary::west);
                //patches at which is the inlet boundary condition set (must have the same length and ordering as the vector of inletBoundarySides!!)
                m_inletBoundaryPatches << 2;
            }
            break;

            default: {
                //walls at which we define Dirichlet boundary conditions
                m_wallBoundarySides.push_back(boundary::north); // sunction side
                m_wallBoundarySides.push_back(boundary::south); // pressure side
                //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
                m_wallBoundaryPatches << 1, 1;
                //inlet sides at which we define Dirichlet boundary conditions
                m_inletBoundarySides.push_back(boundary::west);
                //patches at which is the inlet boundary condition set (must have the same length and ordering as the vector of inletBoundarySides!!)
                m_inletBoundaryPatches << 0;
            }
        }
    }

    T computeWallDistance(uwbINSSolverSteady<T>& navStokes)
    {
        T inletWidth = (((2 * PI*m_rr[m_index_of_profile]) / m_numBlades) / 2.0);
        gsVector<T> Uin(2);
        Uin << m_velocity_absolute_x(m_index_of_profile), m_velocity_absolute_y(m_index_of_profile);
        m_uFreeStream = Uin.norm();
        T Re = m_uFreeStream * inletWidth / m_viscosity;
        std::vector<boxSide> distanceSides;
         gsVector<int> distancePatches;


        switch(m_geomChoice){
            case 1:
            {
                //vector of the sides of the patches from which the wall distance is computed
                distanceSides.push_back(boundary::north);
                distanceSides.push_back(boundary::south);
                //vector of indexes of the patches corresponding to distanceSides
                //length of the vector distancePatches must be equal to the length of vector distanceSides
                distancePatches.resize(2);
                distancePatches << 1, 1;

            }
            break;

        case 2:
        {
            //walls at which we define Dirichlet boundary conditions
            distanceSides.push_back(boundary::west); // sunction side
            distanceSides.push_back(boundary::south); // leading side
            distanceSides.push_back(boundary::south); // sunction side
            distanceSides.push_back(boundary::south); // leading side
            distanceSides.push_back(boundary::south); // sunction side
            //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
            distancePatches.resize(5);
            distancePatches << 0,1,2,9,10;

        }
        break;

        case 3:
        {

            //walls at which we define Dirichlet boundary conditions
            distanceSides.push_back(boundary::north); // sunction side
            distanceSides.push_back(boundary::north); // pressure side
            distanceSides.push_back(boundary::north); // sunction side
            distanceSides.push_back(boundary::south); // pressure side
            distanceSides.push_back(boundary::south); // pressure side
            distanceSides.push_back(boundary::south); // pressure side
            //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
            distancePatches.resize(6);
            distancePatches << 0, 2, 4, 3, 5, 7;

        }
        break;

        case 32:
        {

            //walls at which we define Dirichlet boundary conditions
            distanceSides.push_back(boundary::east); // sunction side
            distanceSides.push_back(boundary::east); // pressure side
            distanceSides.push_back(boundary::east); // sunction side
            //distanceSides.push_back(boundary::east); // pressure side
            distanceSides.push_back(boundary::north); // pressure side
            distanceSides.push_back(boundary::north); // pressure side
            distanceSides.push_back(boundary::north); // pressure side
            //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
            distancePatches.resize(6);
            distancePatches << 0, 2, 4, 3, 5, 7;

        }
        break;

        case 33:
        {

            //walls at which we define Dirichlet boundary conditions
            distanceSides.push_back(boundary::east); // sunction side
            distanceSides.push_back(boundary::east); // pressure side
            distanceSides.push_back(boundary::east); // sunction side
            distanceSides.push_back(boundary::west); // pressure side
            distanceSides.push_back(boundary::west); // pressure side
            distanceSides.push_back(boundary::west); // pressure side
            //patches at which is the wall boundary condition set (must have the same length and ordering as the vector of wallBoundarySides!!)
            distancePatches.resize(6);
            distancePatches << 0, 2, 4, 3, 5, 7;

        }
        break;

            case 4:
            {
                //vector of the sides of the patches from which the wall distance is computed

                distanceSides.push_back(boundary::south);
                distanceSides.push_back(boundary::south);
                distanceSides.push_back(boundary::west);

                //vector of indexes of the patches corresponding to distanceSides
                //length of the vector distancePatches must be equal to the length of vector distanceSides
                distancePatches.resize(3);
                distancePatches << 0, 1, 3;
            }
            break;

            case 5:
            {
                //vector of the sides of the patches from which the wall distance is computed

                distanceSides.push_back(boundary::north);
                distanceSides.push_back(boundary::east);

                //vector of indexes of the patches corresponding to distanceSides
                //length of the vector distancePatches must be equal to the length of vector distanceSides
                distancePatches.resize(2);
                distancePatches << 0, 1;
            }
            break;

            default: {
                //vector of the sides of the patches from which the wall distance is computed

                distanceSides.push_back(boundary::north);
                distanceSides.push_back(boundary::south);

                //vector of indexes of the patches corresponding to distanceSides
                //length of the vector distancePatches must be equal to the length of vector distanceSides
               distancePatches.resize(2);
                distancePatches << 1, 1;
            }
        }



        int numSamplePts = 50; //number of sample points for which the distance to the boundary is computed
        T maxYplus = 2.5; //maximum dimensionless wall distance which is accepted

        //table with wall distance information for every chosen side is printed, if the second from the end input parameter is set as true
        //estimation of the wall distance will be computed, if the last input parameter is set as true
        return navStokes.computeDimensionlessWallDistance(distancePatches, distanceSides, m_viscosity, Re, m_uFreeStream, maxYplus, numSamplePts, true, true);
    }

    void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, T velocity_x, T velocity_y, T velocity_blade)
        {

            /*std::ostringstream strs_vel_x;
            std::ostringstream strs_vel_y;
            std::ostringstream strs_velblade;

            strs_vel_x << velocity_x;
            std::string str_x = strs_vel_x.str();
            strs_vel_y << velocity_y;
            std::string str_y = strs_vel_y.str();
            strs_velblade << velocity_blade;
            std::string str2 = strs_velblade.str();*/

            std::string str_x = std::to_string(velocity_x);
            std::string str_y = std::to_string(velocity_y);
            std::string str2 = std::to_string(velocity_blade);

            gsFunctionExpr<T> Uin(str_x, str_y, 2);
            gsFunctionExpr<T> Ublade("0", str2, 2);

            Base::defineBCs_NS(bcInfo, Uin, Ublade);
        }

    void defineBCsADR(gsBoundaryConditions<T>& bcInfo, T bcInlet, T bcBlade)
    {
        std::string strInlet = std::to_string(bcInlet);
        std::string strBlade = std::to_string(bcBlade);

        gsFunctionExpr<T> Uin(strInlet, 2);
        gsFunctionExpr<T> Ublade(strBlade, 2);

        Base::defineBCs_NS(bcInfo, Uin, Ublade);
    }

    void solvePoissonEquation(uwbTMSolverKOmega<T>& turbSolver, bool plotMeshes = true)
    {
        int numRefinePoisson = 4;

        gsMultiBasis<T> tbasisPoisson(m_patches); // basis for RANS equations
        for (int i = 0; i < numRefinePoisson; ++i)
            tbasisPoisson.uniformRefine();
        gsMatrix<T> box_u0Poisson(2, 2);
        box_u0Poisson << 0, 1, 0, 0;
        switch(m_geomChoice)
        {
            case 1:
            {
              tbasisPoisson.refine(1, box_u0Poisson);
            }
            break;

            case 2:
            {
                //to do tbasisPoisson local refinement in case of oscillation
            }
            break;

            case 3:
            {
                //to do tbasisPoisson local refinement in case of oscillation
            }
            break;

        case 32:
        {
            //to do tbasisPoisson local refinement in case of oscillation
        }
        break;

        case 33:
        {
            //to do tbasisPoisson local refinement in case of oscillation
        }
        break;

            case 4:
            {
                //to do tbasisPoisson local refinement in case of oscillation
            }
            break;

            case 5:
            {
                //to do tbasisPoisson local refinement in case of oscillation
            }
            break;

            default: {
                tbasisPoisson.refine(1, box_u0Poisson);
            }

        }


        if (plotMeshes)
        {
            for (unsigned int index_of_mesh = 0; index_of_mesh < m_patches.nPatches(); index_of_mesh++)
            {
                 gsMesh<> mesh;
                std::ostringstream strs_mesh;
                strs_mesh << index_of_mesh;
                makeMesh(tbasisPoisson.at(index_of_mesh), mesh, 10);
                m_patches.patch(index_of_mesh).evaluateMesh(mesh);
                gsWriteParaview(mesh, m_sGeometry + "meshPoissonPatch" + strs_mesh.str() + m_strProfile);
            }

        }


        gsFunctionExpr<T> fw("1", 2);
        gsFunctionExpr<T> gw("0", 2);
        gsFunctionExpr<T> wallw("0.0", 2);
        gsBoundaryConditions<T> bcInfow;



        switch(m_geomChoice){
            case 1:
            {
                bcInfow.addCondition(0, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(0, boundary::south, condition_type::neumann, gw);
                bcInfow.addCondition(0, boundary::west, condition_type::neumann, gw);
                bcInfow.addCondition(2, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(2, boundary::south, condition_type::neumann, gw);
                bcInfow.addCondition(2, boundary::east, condition_type::neumann, gw);
                bcInfow.addCondition(1, boundary::north, condition_type::dirichlet, wallw);
                bcInfow.addCondition(1, boundary::south, condition_type::dirichlet, wallw);

            }
            break;

            case 2:
            {

                bcInfow.addCondition(4, boundary::east, condition_type::neumann, gw);
                bcInfow.addCondition(4, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(5, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(6, boundary::east, condition_type::neumann, gw);
                 bcInfow.addCondition(6, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(7, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(11, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(12, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(8, boundary::north, condition_type::neumann, gw);
                 bcInfow.addCondition(8, boundary::east, condition_type::neumann, gw);
                bcInfow.addCondition(3, boundary::east, condition_type::neumann, gw);
                bcInfow.addCondition(0, boundary::west, condition_type::dirichlet, wallw);
                bcInfow.addCondition(1, boundary::south, condition_type::dirichlet, wallw);
                bcInfow.addCondition(2, boundary::south, condition_type::dirichlet, wallw);
                bcInfow.addCondition(9, boundary::south, condition_type::dirichlet, wallw);
                bcInfow.addCondition(10, boundary::south, condition_type::dirichlet, wallw);

            }
            break;

            case 3:
            {
                bcInfow.addCondition(0, boundary::west, condition_type::neumann, gw);
                bcInfow.addCondition(1, boundary::south, condition_type::neumann, gw);
                bcInfow.addCondition(1, boundary::west, condition_type::neumann, gw);
                bcInfow.addCondition(6, boundary::north, condition_type::neumann, gw);
                 bcInfow.addCondition(8, boundary::east, condition_type::neumann, gw);
                bcInfow.addCondition(8, boundary::south, condition_type::neumann, gw);
                bcInfow.addCondition(0, boundary::north, condition_type::dirichlet, wallw);
                bcInfow.addCondition(2, boundary::north, condition_type::dirichlet, wallw);
                bcInfow.addCondition(4, boundary::north, condition_type::dirichlet, wallw);
                bcInfow.addCondition(3, boundary::south, condition_type::dirichlet, wallw);
                bcInfow.addCondition(5, boundary::south, condition_type::dirichlet, wallw);
                bcInfow.addCondition(7, boundary::south, condition_type::dirichlet, wallw);
            }
            break;

        case 32:
        {
            bcInfow.addCondition(0, boundary::south, condition_type::neumann, gw);
            bcInfow.addCondition(1, boundary::north, condition_type::neumann, gw);
            bcInfow.addCondition(1, boundary::west, condition_type::neumann, gw);
            bcInfow.addCondition(8, boundary::east, condition_type::neumann, gw);
            bcInfow.addCondition(8, boundary::north, condition_type::neumann, gw);
            bcInfow.addCondition(6, boundary::east, condition_type::neumann, gw);
            bcInfow.addCondition(0, boundary::east, condition_type::dirichlet, wallw);
            bcInfow.addCondition(2, boundary::east, condition_type::dirichlet, wallw);
            bcInfow.addCondition(4, boundary::east, condition_type::dirichlet, wallw);
            //bcInfow.addCondition(6, boundary::east, condition_type::dirichlet, wallw);
            bcInfow.addCondition(3, boundary::north, condition_type::dirichlet, wallw);
            bcInfow.addCondition(5, boundary::north, condition_type::dirichlet, wallw);
            bcInfow.addCondition(7, boundary::north, condition_type::dirichlet, wallw);
        }
        break;

        case 33:
        {
            bcInfow.addCondition(9, boundary::south, condition_type::neumann, gw);
            bcInfow.addCondition(9, boundary::east, condition_type::neumann, gw);
            bcInfow.addCondition(1, boundary::west, condition_type::neumann, gw);
            bcInfow.addCondition(6, boundary::east, condition_type::neumann, gw);
             bcInfow.addCondition(8, boundary::west, condition_type::neumann, gw);
            bcInfow.addCondition(8, boundary::north, condition_type::neumann, gw);
            bcInfow.addCondition(0, boundary::east, condition_type::dirichlet, wallw);
            bcInfow.addCondition(2, boundary::east, condition_type::dirichlet, wallw);
            bcInfow.addCondition(4, boundary::east, condition_type::dirichlet, wallw);
            bcInfow.addCondition(3, boundary::west, condition_type::dirichlet, wallw);
            bcInfow.addCondition(5, boundary::west, condition_type::dirichlet, wallw);
            bcInfow.addCondition(7, boundary::west, condition_type::dirichlet, wallw);
        }
        break;


            case 4:
            {
                bcInfow.addCondition(0, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(1, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(2, boundary::east, condition_type::neumann, gw);
                bcInfow.addCondition(2, boundary::north, condition_type::neumann, gw);
                //bcInfow.addCondition(3, boundary::south, condition_type::neumann, gw);
                bcInfow.addCondition(3, boundary::east, condition_type::neumann, gw);
                bcInfow.addCondition(4, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(4, boundary::east, condition_type::neumann, gw);
                bcInfow.addCondition(0, boundary::south, condition_type::dirichlet, wallw);
                bcInfow.addCondition(1, boundary::south, condition_type::dirichlet, wallw);
               // bcInfow.addCondition(3, boundary::north, condition_type::dirichlet, wallw);
                bcInfow.addCondition(3, boundary::west, condition_type::dirichlet, wallw);
            }
            break;

            case 5:
            {
                bcInfow.addCondition(2, boundary::west, condition_type::neumann, gw);
                bcInfow.addCondition(2, boundary::south, condition_type::neumann, gw);
                bcInfow.addCondition(2, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(0, boundary::west, condition_type::neumann, gw);
                bcInfow.addCondition(1, boundary::south, condition_type::neumann, gw);
                bcInfow.addCondition(4, boundary::east, condition_type::neumann, gw);

                bcInfow.addCondition(5, boundary::east, condition_type::neumann, gw);
                bcInfow.addCondition(5, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(6, boundary::east, condition_type::neumann, gw);
                bcInfow.addCondition(6, boundary::north, condition_type::neumann, gw);


                bcInfow.addCondition(0, boundary::north, condition_type::dirichlet, wallw);
                bcInfow.addCondition(1, boundary::east, condition_type::dirichlet, wallw);
            }
            break;

            default: {
                bcInfow.addCondition(0, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(0, boundary::south, condition_type::neumann, gw);
                bcInfow.addCondition(0, boundary::west, condition_type::neumann, gw);
                bcInfow.addCondition(2, boundary::north, condition_type::neumann, gw);
                bcInfow.addCondition(2, boundary::south, condition_type::neumann, gw);
                bcInfow.addCondition(2, boundary::east, condition_type::neumann, gw);
                bcInfow.addCondition(1, boundary::north, condition_type::dirichlet, wallw);
                bcInfow.addCondition(1, boundary::south, condition_type::dirichlet, wallw);
            }
        }

        gsInfo << "\nSolving Poisson equation.\n";

        turbSolver.setPoissonSolution(m_patches, tbasisPoisson, bcInfow, fw, true, m_plot_pts);

        gsInfo << "Poisson equation resolved.\n\n";
        turbSolver.plotWallDistance(m_sGeometry + "wall_distance" + m_strProfile, m_plot_pts);
    }

    gsMultiPatch<T> makeMultiPatch()
    {
        T rotation_of_blade;
        T uniformity_param;

        switch (m_setting_of_blade_rotation) {
          case 1: //maximal
             rotation_of_blade = 10; //in degrees
             uniformity_param = 0.005;   //parameter for mesh
          break;

          case 2: //minimal
            rotation_of_blade = -15; //in degrees
            uniformity_param = 0.01;
          break;

          default: //default is optimal
             rotation_of_blade = 0;
             uniformity_param = 0.01;
          break;
        };
        /*
        //------parameter definition------
        //parameters that can be changed
        T camberX;    //x-ova souradnice max polohy strednice
        T camberY;    //y-ova souradnice max polohy strednice
        T leadingAngle;      //vstupni uhel strednice
        T trailingAngle;        //vystupni uhel strednice
        T thicknessX;         //x-ova souradnice max prohnuti tloustky
        T thicknessY;            //y-ova souranice max prohnuti tloustky
        T outputAngle;         //vystupni uhel tloustky
        T radius;            //polomer osk. kruz v nabehu

        //parameters that can't be changed -> need to change file uwbProfileOptimization.h
        T endingOffset;     //odsazeni tloustky na konci
        T chordLength;       //delka lopatky -> v pripade zmeny potreba upravit ohranicujici kanal
        m_angle;   //otoceni do mrize
        T rotationCenterX;   //zde neni potreba, slouzi pro otoceni profilu do lop. mrize
        T rotationCenterY;
        */

        gsVector<T> camber_x(m_numBladeprofiles);
        gsVector<T> camber_y(m_numBladeprofiles);
        gsVector<T> trailing_angle(m_numBladeprofiles);
        gsVector<T> thickness_x(m_numBladeprofiles);
        gsVector<T> thickness_y(m_numBladeprofiles);
        gsVector<T> ending_offset(m_numBladeprofiles);
        gsVector<T> output_angle(m_numBladeprofiles);
        gsVector<T> radius(m_numBladeprofiles);
        gsVector<T> chord_length(m_numBladeprofiles);
        gsVector<T> rotation_center_x(m_numBladeprofiles);
        gsVector<T> rotation_center_y(m_numBladeprofiles);
        gsVector<T> right_angle((m_numBladeprofiles));
        right_angle.setConstant(m_numBladeprofiles, 90);
        gsVector<T> rotate_angle((m_numBladeprofiles));
        rotate_angle.setConstant(m_numBladeprofiles, rotation_of_blade);

        // SETTING PARAMETERS
        m_rr << 0.175, 0.229, 0.283, 0.338, 0.392, 0.446, 0.5;
        camber_x << 0.425908805, 0.419, 0.416039882, 0.416576484, 0.416716376, 0.419203131, 0.43280567; // 0.411886743 na originalne druhe pozici
        //gsInfo << camber_x << "\n";
        camber_y << 0.095896144, 0.064997436, 0.037083707, 0.02724709, 0.024356984, 0.023262639, 0.019802704;
        m_leading_angle << 38.47692, 21.399859, 12.641776, 9.28275, 8.38282, 8.338553, 8.446091;
        m_leading_angle = m_leading_angle*PI/180;
        trailing_angle << 20.717182, 10.721469, 6.371868, 4.702573, 4.217487, 4.164196, 4.233241;
        trailing_angle = trailing_angle*pi/180;
        thickness_x << 0.281408739, 0.275195139, 0.273161229, 0.272459059, 0.272310247, 0.272053329, 0.271637628;
        thickness_y << 0.064950942, 0.044591162, 0.030435219, 0.020799775, 0.01464564, 0.011033004, 0.009989223;
        //ending_offset << 0.001184, 0.000812, 0.000555, 0.000379, 0.000267, 0.000201, 0.000182;

        ending_offset << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
        output_angle << 8.381088, 5.891258, 4.042351, 2.765094, 1.950139, 1.466739, 1.329705;
        output_angle = output_angle*PI/180;
        radius << 0.020637, 0.014575, 0.007271, 0.003687, 0.001955, 0.001012, 0.000828;
        chord_length << 0.36686418, 0.42813715,	0.486298848, 0.545457321, 0.595371262, 0.615928719, 0.584588959;
        m_angle << 53.696487,	41.265848, 33.007703, 27.603276, 24.437586, 22.893162, 21.162381;
        m_angle = (right_angle - (m_angle + rotate_angle))*PI/180;
        rotation_center_x << 0.494758396, 0.469497406, 0.444542963, 0.417724545, 0.390108787, 0.361175154, 0.330805204;
        rotation_center_y << 0.060495569, 0.028225794, 0.00125711, -0.006884641, -0.010228889, -0.010435203, -0.00079539;

        gsVector<T> length_x1(m_numBladeprofiles);
        length_x1 << -0.186559, -0.177974, -0.169389, -0.160645, -0.152061, -0.143476, -0.134891;
        length_x1 *= 2.;
        T length_x2 = 0.433544;

        /*T length_x1_fixed = -0.8;
        T length_x2_fixed = 1;;*/
        //  T rotate_angle = 0.0 * PI/180;

        uwbGeometryCreators<T> geometryCreator;

        switch (m_geomChoice)
        {
            case 1:
            {
                    m_patches = geometryCreator.DomainBetweenBladeProfiles1(m_index_of_profile,
                    length_x1[m_index_of_profile],//length_x1_fixed,
                    length_x2, //length_x2_fixed,
                    ((2 * PI*m_rr[m_index_of_profile]) / m_numBlades),
                    camber_x[m_index_of_profile],
                    camber_y[m_index_of_profile],
                    m_leading_angle[m_index_of_profile],
                    trailing_angle[m_index_of_profile],
                    thickness_x[m_index_of_profile],
                    thickness_y[m_index_of_profile],
                    ending_offset[m_index_of_profile],
                    output_angle[m_index_of_profile],
                    radius[m_index_of_profile],
                    chord_length[m_index_of_profile],
                    m_angle[m_index_of_profile],
                    rotation_center_x[m_index_of_profile],
                    rotation_center_y[m_index_of_profile],
                    uniformity_param,
                    m_kvfit_knots, m_coarse, m_geomParams);
            }
                break;

            case 2:
            {
                    m_patches = geometryCreator.DomainAroundBladeProfile2(m_index_of_profile,
                    length_x1[m_index_of_profile],//length_x1_fixed,
                    length_x2, //length_x2_fixed,
                    ((2 * PI*m_rr[m_index_of_profile]) / m_numBlades),
                    camber_x[m_index_of_profile],
                    camber_y[m_index_of_profile],
                    m_leading_angle[m_index_of_profile],
                    trailing_angle[m_index_of_profile],
                    thickness_x[m_index_of_profile],
                    thickness_y[m_index_of_profile],
                    ending_offset[m_index_of_profile],
                    output_angle[m_index_of_profile],
                    radius[m_index_of_profile],
                    chord_length[m_index_of_profile],
                    m_angle[m_index_of_profile],
                    rotation_center_x[m_index_of_profile],
                    rotation_center_y[m_index_of_profile],
                    m_kvfit_knots, m_coarse, m_geomParams);
            }
                break;

            case 3:
            {
                    m_patches = geometryCreator.DomainBetweenBladeProfiles3(m_index_of_profile,
                    length_x1[m_index_of_profile],//length_x1_fixed,
                    length_x2, //length_x2_fixed,
                    ((2 * PI*m_rr[m_index_of_profile]) / m_numBlades),
                    camber_x[m_index_of_profile],
                    camber_y[m_index_of_profile],
                    m_leading_angle[m_index_of_profile],
                    trailing_angle[m_index_of_profile],
                    thickness_x[m_index_of_profile],
                    thickness_y[m_index_of_profile],
                    ending_offset[m_index_of_profile],
                    output_angle[m_index_of_profile],
                    radius[m_index_of_profile],
                    chord_length[m_index_of_profile],
                    m_angle[m_index_of_profile],
                    rotation_center_x[m_index_of_profile],
                    rotation_center_y[m_index_of_profile],
                    uniformity_param,
                    m_kvfit_knots, m_coarse, m_geomParams);

            }
                break;

        case 32:
        {
                m_patches = geometryCreator.DomainBetweenBladeProfiles3b(m_index_of_profile,
                length_x1[m_index_of_profile],//length_x1_fixed,
                length_x2, //length_x2_fixed,
                ((2 * PI*m_rr[m_index_of_profile]) / m_numBlades),
                camber_x[m_index_of_profile],
                camber_y[m_index_of_profile],
                m_leading_angle[m_index_of_profile],
                trailing_angle[m_index_of_profile],
                thickness_x[m_index_of_profile],
                thickness_y[m_index_of_profile],
                ending_offset[m_index_of_profile],
                output_angle[m_index_of_profile],
                radius[m_index_of_profile],
                chord_length[m_index_of_profile],
                m_angle[m_index_of_profile],
                rotation_center_x[m_index_of_profile],
                rotation_center_y[m_index_of_profile],
                uniformity_param,
                m_kvfit_knots, m_coarse, m_geomParams);

        }
            break;

        case 33:
        {
                m_patches = geometryCreator.DomainBetweenBladeProfiles3c(m_index_of_profile,
                length_x1[m_index_of_profile],//length_x1_fixed,
                length_x2, //length_x2_fixed,
                ((2 * PI*m_rr[m_index_of_profile]) / m_numBlades),
                camber_x[m_index_of_profile],
                camber_y[m_index_of_profile],
                m_leading_angle[m_index_of_profile],
                trailing_angle[m_index_of_profile],
                thickness_x[m_index_of_profile],
                thickness_y[m_index_of_profile],
                ending_offset[m_index_of_profile],
                output_angle[m_index_of_profile],
                radius[m_index_of_profile],
                chord_length[m_index_of_profile],
                m_angle[m_index_of_profile],
                rotation_center_x[m_index_of_profile],
                rotation_center_y[m_index_of_profile],
                uniformity_param,
                m_kvfit_knots, m_coarse, m_geomParams);

        }
            break;

            case 4:
            {
                    m_patches = geometryCreator.DomainAroundBladeProfile4(m_index_of_profile,
                    length_x1[m_index_of_profile],//length_x1_fixed,
                    length_x2, //length_x2_fixed,
                    ((2 * PI*m_rr[m_index_of_profile]) / m_numBlades),
                    camber_x[m_index_of_profile],
                    camber_y[m_index_of_profile],
                    m_leading_angle[m_index_of_profile],
                    trailing_angle[m_index_of_profile],
                    thickness_x[m_index_of_profile],
                    thickness_y[m_index_of_profile],
                    ending_offset[m_index_of_profile],
                    output_angle[m_index_of_profile],
                    radius[m_index_of_profile],
                    chord_length[m_index_of_profile],
                    m_angle[m_index_of_profile],
                    rotation_center_x[m_index_of_profile],
                    rotation_center_y[m_index_of_profile],
                    uniformity_param,
                    m_kvfit_knots, m_coarse, m_geomParams);

            }
                break;

            case 5:
            {
                    m_patches = geometryCreator.DomainBetweenBladeProfiles5(m_index_of_profile,
                    length_x1[m_index_of_profile],//length_x1_fixed,
                    length_x2, //length_x2_fixed,
                    ((2 * PI*m_rr[m_index_of_profile]) / m_numBlades),
                    camber_x[m_index_of_profile],
                    camber_y[m_index_of_profile],
                    m_leading_angle[m_index_of_profile],
                    trailing_angle[m_index_of_profile],
                    thickness_x[m_index_of_profile],
                    thickness_y[m_index_of_profile],
                    ending_offset[m_index_of_profile],
                    output_angle[m_index_of_profile],
                    radius[m_index_of_profile],
                    chord_length[m_index_of_profile],
                    m_angle[m_index_of_profile],
                    rotation_center_x[m_index_of_profile],
                    rotation_center_y[m_index_of_profile],
                    uniformity_param,
                    m_kvfit_knots, m_coarse, m_geomParams);
            }
                break;

            default:
            {
                    m_patches = geometryCreator.DomainBetweenBladeProfilesUniformKnotVector(m_index_of_profile,
                    length_x1[m_index_of_profile],//length_x1_fixed,
                    length_x2, //length_x2_fixed,
                    ((2 * PI*m_rr[m_index_of_profile]) / m_numBlades),
                    camber_x[m_index_of_profile],
                    camber_y[m_index_of_profile],
                    m_leading_angle[m_index_of_profile],
                    trailing_angle[m_index_of_profile],
                    thickness_x[m_index_of_profile],
                    thickness_y[m_index_of_profile],
                    ending_offset[m_index_of_profile],
                    output_angle[m_index_of_profile],
                    radius[m_index_of_profile],
                    chord_length[m_index_of_profile],
                    m_angle[m_index_of_profile],
                    rotation_center_x[m_index_of_profile],
                    rotation_center_y[m_index_of_profile],
                    uniformity_param);
            }
        }


        m_patches.addAutoBoundaries();
        return m_patches;
    }

    void setProfileInputVelocities(T omega)
    {
        switch (m_setting_of_blade_rotation) {
          case 1: //maximal
             m_flow_rate = 8.7;
          break;

          case 2: //minimal
            m_flow_rate = 1.58;
          break;

          default: //default is optimal
             m_flow_rate = 5.73;
          break;
        };

        T omega_rot = 2*PI*538.0/60.0;
        //T gravity = 9.81;
        //T etah = 0.945;
        //T H = 9.88929;
        //T Hn = H*etah;

        //T m_flow_rate = 5.73; //5.73 - optimal, 8.7 - maximal, 1.58 - minimal
        T D_out = 0.5 * 2;
        T d_out = 0.122926 * 2;
        T vel_mer_out = m_flow_rate / (PI*(D_out*D_out - d_out*d_out) / 4.0);
        T vel_mer_in = vel_mer_out;                //meridial velocity

        gsVector<T> velocity_absolute_y_data(m_numBladeprofiles);
        gsVector<T> velocity_relative_x(m_numBladeprofiles);
        gsVector<T> velocity_relative_y(m_numBladeprofiles);
        gsVector<T> v_u(m_numBladeprofiles);
        //gsVector<T> v_target(m_numBladeprofiles);
        //gsVector<T> v_t_input(m_numBladeprofiles);

        for (unsigned i = 0; i<m_numBladeprofiles; i++) {
            m_velocity_blade(i) = -m_rr(i)*omega;
        }

        for (unsigned i = 0; i<m_numBladeprofiles; i++) {
            v_u(i) = m_rr(i)*omega_rot;
        }

        switch(m_setting_of_guide_vanes)
        {
            case 1: //56° guide vanes for opt. runner blade, 68° guide vanes for max. runner blade
            {
                switch(m_setting_of_blade_rotation)
                {
                    case 1: // 68° guide vanes, maximal rotation of runner blade
                    {
                        velocity_relative_x << 8.04, 9.22, 10.14, 11.03, 12.01, 13.51, 14.29;
                        velocity_absolute_y_data << 5.99, 5.91, 5.29, 4.95, 4.8, 4.72, 4.27;
                        for (unsigned i = 0; i<m_numBladeprofiles; i++) {
                           velocity_relative_y(i) =  v_u(i) - velocity_absolute_y_data(i);
                        }
                    }
                    break;

                    default: //56° guide vanes, optimal rotation of runner blade
                    {
                        velocity_relative_x << 6.45, 6.75, 7.3, 7.77, 8.31, 8.31, 8.92;
                        velocity_absolute_y_data << 7.44, 7.21, 6.27, 5.72, 5.4, 5.17, 4.48;
                        for (unsigned i = 0; i<m_numBladeprofiles; i++) {
                           velocity_relative_y(i) =  v_u(i) - velocity_absolute_y_data(i);
                        }
                    }
                    break;

                }
            }
            break;

            case 2: //60° guide vanes for opt. runner blade, 72° guide vanes for max. runner blade
            {
                switch(m_setting_of_blade_rotation)
                {
                    case 1: // 72° guide vanes, maximal rotation of runner blade
                    {
                        velocity_relative_x << 8.21, 9.47, 10.48, 11.44, 12.52, 14.12, 15.01;
                        velocity_absolute_y_data << 4.88, 4.84, 4.38, 4.12, 4.03, 3.99, 3.64;
                        for (unsigned i = 0; i<m_numBladeprofiles; i++) {
                           velocity_relative_y(i) =  v_u(i) - velocity_absolute_y_data(i);
                        }
                    }
                    break;

                    default: // 60° guide vanes, maximal rotation of runner blade
                    {
                        velocity_relative_x << 6.54, 6.88, 7.5, 8.03, 8.64, 8.64, 9.33;
                        velocity_absolute_y_data << 6.39, 6.24, 5.47, 5.03, 4.8, 4.63, 4.05;
                        for (unsigned i = 0; i<m_numBladeprofiles; i++) {
                           velocity_relative_y(i) =  v_u(i) - velocity_absolute_y_data(i);
                        }
                    }
                    break;
                }
            }
            break;

            default: //compute with formula
            {
                gsVector<T> angle_input(m_numBladeprofiles);
                angle_input = m_angle - m_leading_angle;
                //      gsInfo << m_angle*180/PI;
                //      gsInfo<<m_leading_angle*180/PI;
                //    gsInfo<< "\n m_angle between input velocity and x-axis (degrees) \n";
                //    gsInfo << angle_input*180/PI;
                for (unsigned i = 0; i<m_numBladeprofiles; i++) {
                velocity_relative_x(i) = vel_mer_in;

                // //-------leading angle possibility
                //velocity_relative_y(i) = vel_mer_in*math::tan(angle_input(i));

                // //-------formula possibility
                //v_target(i) =(0.4*m_rr(i)+0.1)*vel_mer_out;
                //v_t_input(i) = v_target(i) + ((etah*gravity*Hn)/v_u(i));
                //velocity_relative_y(i) =  v_u(i) - v_t_input(i);

                // //-------data with deviation
                //velocity_absolute_y_data << -0.28569, 4.35998, 8.07131, 11.2644, 14.018, 16.5, 18.7736;  //94.98568% of water flow angle with formula
                //velocity_absolute_y_data << -0.27613, 4.18687, 7.64825, 10.5193, 12.9108, 15, 16.8598;  //91.80894% of water flow angle with formula
                //velocity_absolute_y_data << -0.27971, 4.25139, 7.80418, 10.7908, 13.3101, 15.5357, 17.5373; //93
                switch(m_setting_of_blade_rotation)
                {
                    case 1: //maximal
                    //velocity_absolute_y_data << -0.28273, 4.33645, 8.104, 11.4163, 14.3233, 16.9817, 19.4471; //94
                    velocity_absolute_y_data << -0.74342, 3.67571, 7.26132, 10.38162, 13.08805, 15.53386, 17.77616; //92
                    break;

                    default: //optimal
                    velocity_absolute_y_data << -0.28272, 4.30592, 7.93756, 11.0259, 13.6596, 16.0096, 18.1491; //94


                    break;
                }

                velocity_relative_y(i) =  velocity_absolute_y_data(i);
                }
            }
            break;
        }

        for (unsigned i = 0; i < m_numBladeprofiles; i++) {
            m_velocity_absolute_x(i) = velocity_relative_x(i);
            if (m_UyInZero)
                m_velocity_absolute_y(i) = 0.;
            else
                m_velocity_absolute_y(i) = velocity_relative_y(i) + m_velocity_blade(i);
        }

        /*
              //gsInfo << "vel_mer_in " << vel_mer_in;

                gsInfo<< "\n m_velocity_blade\n";
                gsInfo<< m_velocity_blade;
                gsInfo<< "\n";

                gsInfo<< "m_velocity_absolute_x\n";
                gsInfo<< m_velocity_absolute_x;
                gsInfo<< "\n";

                gsInfo<< "m_velocity_absolute_y\n";
                gsInfo<< m_velocity_absolute_y;
                gsInfo<< "\n";

                gsInfo<< "velocity_relative_x\n";
                gsInfo<< velocity_relative_x;
                gsInfo<< "\n";

                gsInfo<< "velocity_relative_y\n";
                gsInfo<< velocity_relative_y;
                gsInfo<< "\n";
        */
    }

     void refineBasisUniformZones(gsMultiBasis<T>& tbasis, int id_patch, int direction, int position, int knot, int num_of_refine)
     {
         gsMatrix<T> box(2, 2);
         for (int i = 0; i < num_of_refine; i++)
         {
             // in each refinement step take the first numKnot_v non_zero knots
             const gsTensorBSplineBasis<2, T>*  basis = dynamic_cast<const gsTensorBSplineBasis<2, T>*>(&(tbasis.basis(id_patch))); //basis of the first patch
             gsVector<T> RefineKnot((i + 1)*(knot)+1);
             RefineKnot.setZero((i + 1)*(knot)+1);
             int sizeKnots = basis->knots(direction).size() - 1;
             for (int k = 0; k < (i + 1)*(knot)+1; k++) //gsVector of first numKnots knots
             {
                 if (position == 0)
                    RefineKnot(k) = basis->knot(direction, basis->degree(direction) + k); // first numKnot knots and 0 knot in direction v - all patches have same knots in v direction
                 else
                    RefineKnot(k) = basis->knot(direction, sizeKnots - (basis->degree(direction) + k)); //last numKnot knots and 1 knot in direction v
             }
             for (int j = 0; j < (i + 1)*(knot); j++)
             {
                if (direction == 0)
                 {
                    if (position == 0)
                        box << RefineKnot(j), RefineKnot(j+1), 0, 0;
                    else
                        box << RefineKnot(j+1), RefineKnot(j), 0, 0;
                }
                else
                {
                    if (position == 0)
                        box << 0, 0, RefineKnot(j), RefineKnot(j+1);
                    else
                        box << 0, 0, RefineKnot(j+1), RefineKnot(j);
                }

                 tbasis.refine(id_patch,  box);
             }
         }
     }

     /*void refineBasisLocalZones(gsMultiBasis<T>& tbasis, int id_patch, int direction, int position, int knot, int num_of_refine)
     {
         gsMatrix<T> box(2, 2);
         for (int i = 0; i < num_of_refine; i++)
         {
             // in each refinement step take the first numKnot_v non_zero knots
             const gsTensorBSplineBasis<2, T>*  basis = dynamic_cast<const gsTensorBSplineBasis<2, T>*>(&(tbasis.basis(id_patch))); //basis of the first patch
             gsVector<T> RefineKnot((i + 1)*(knot)+1);
             RefineKnot.setZero((i + 1)*(knot)+1);
             int sizeKnots = basis->knots(direction).size() - 1;
             for (int k = 0; k < (i + 1)*(knot)+1; k++) //gsVector of first numKnots knots
             {
                 if (position == 0)
                    RefineKnot(k) = basis->knot(direction, basis->degree(direction) + k); // first numKnot knots and 0 knot in direction v - all patches have same knots in v direction
                 else
                    RefineKnot(k) = basis->knot(direction, sizeKnots - (basis->degree(direction) + k)); //last numKnot knots and 1 knot in direction v
             }
             for (int j = 0; j < knot; j++)
             {
                if (direction == 0)
                {
                    if (position == 0)
                        box << RefineKnot(j), RefineKnot(j+1), 0, 0;
                    else
                        box << RefineKnot(j+1), RefineKnot(j), 0, 0;
                }
                else
                {
                    if (position == 0)
                        box << 0,0, RefineKnot(j), RefineKnot(j+1);
                    else
                        box << 0,0, RefineKnot(j+1), RefineKnot(j);
                }

                 tbasis.refine(id_patch,  box);
             }
         }
     }*/

     void plotSolutionField(gsField<T> solField, std::string const & fn, bool plotMesh = false)
     {
         gsWriteParaview<T>(solField, m_sGeometry + fn + m_strProfile, m_plot_pts, plotMesh);
     }

public:
    const gsVector<T>& getVelocityAbsoluteX() const { return m_velocity_absolute_x; }
    const gsVector<T>& getVelocityAbsoluteY() const { return m_velocity_absolute_y; }
    const gsVector<T>& getVelocityBlade() const { return m_velocity_blade; }
    const gsVector<T>& getRr() const { return m_rr; }
    const T getRr(int i) const { return m_rr(i); }
    const T getFlowRate() const { return m_flow_rate; }
    const T getFreeStream() const { return m_uFreeStream; }

protected:
    std::string m_strProfile;
    unsigned m_index_of_profile;
    unsigned m_numBladeprofiles;
    unsigned m_numBlades;
    int m_geomChoice;
    std::vector<T> m_kvfit_knots;
    bool m_coarse;
    gsVector<T> m_geomParams;
    T m_flow_rate;
    T m_uFreeStream;

    bool m_UyInZero;

    gsVector<T> m_velocity_absolute_x;
    gsVector<T> m_velocity_absolute_y;
    gsVector<T> m_rr;
    gsVector<T> m_leading_angle;
    gsVector<T> m_angle;
    gsVector<T> m_velocity_blade;

    //parameters Examples
    int m_setting_of_blade_rotation; //0-optimal, 1-maximal, 2-minimal
    int m_setting_of_guide_vanes;

    // members from uwbRANSExampleSetting
    using Base::m_viscosity;
    using Base::m_plot_pts;
    using Base::m_patches;
    using Base::m_wallBoundarySides;
    using Base::m_inletBoundarySides;
    using Base::m_wallBoundaryPatches;
    using Base::m_inletBoundaryPatches;
    using Base::m_sGeometry;
    using Base::PI;

}; //uwbRANSProfileExample

} //namespace gismo
