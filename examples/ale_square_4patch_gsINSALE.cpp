// 4-patch rotating-square ALE example using the new gsINSSolverUnsteadyALE
// - Rebuilds parameterization each time step and splits into 4 patches
// - Uses gsBarrierPatch for parametrization quality when requested

#include <gismo.h>
#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSSolverALE.h>

using namespace gismo;

// Build a single patch annulus-like surface bounded by an outer square [-1,2]^2
// and an inner unit square [0,1]^2 rotated around (0.5,0.5) by 'angle' (radians).
static gsMultiPatch<> makeConnectedBoundarySinglePatch(real_t angleRad)
{
    // Outer boundary: closed polyline square
    gsMatrix<> outerCoefs(5, 2);
    outerCoefs << -1, -1,
                   2, -1,
                   2,  2,
                  -1,  2,
                  -1, -1;
    gsKnotVector<> kv(0, 1, 3, 2);
    gsBSpline<> outerCurve(kv, outerCoefs);

    // Inner boundary: unit square, rotated by 'angleRad' around (0.5,0.5)
    gsMatrix<> innerCoefs(5, 2);
    innerCoefs << 0, 0,
                  1, 0,
                  1, 1,
                  0, 1,
                  0, 0;
    gsBSpline<> innerCurve(kv, innerCoefs);

    if (std::abs(angleRad) > 1e-14)
    {
        innerCurve.translate(gsEigen::Vector2d(-0.5, -0.5));
        innerCurve.rotate(angleRad);
        innerCurve.translate(-gsEigen::Vector2d(-0.5, -0.5));
    }

    // Fixed seam: do not rotate parameter start; keep u seam fixed across time
    // Ensure compatible knot vectors between outer and inner curves
    std::vector<real_t> diff;
    outerCurve.knots(0).difference(innerCurve.knots(0), diff); innerCurve.insertKnots(diff.begin(), diff.end());
    innerCurve.knots(0).difference(outerCurve.knots(0), diff); outerCurve.insertKnots(diff.begin(), diff.end());

    // Build surface with OUTER at v=1 (north), INNER at v=0 (south)
    gsMatrix<> surfCoefs(2 * innerCurve.numCoefs(), 2);
    surfCoefs.topRows(innerCurve.numCoefs())    = outerCurve.coefs();  // v=1: outer
    surfCoefs.bottomRows(innerCurve.numCoefs()) = innerCurve.coefs();  // v=0: inner

    gsKnotVector<> uKnot = innerCurve.knots(0);
    gsKnotVector<> vKnot(0, 1, 0, 2);
    uKnot.affineTransformTo(0, 4); // map u in [0,4]

    gsTensorBSpline<2> surf(uKnot, vKnot, surfCoefs);
    surf.degreeElevate(1);

    gsMultiPatch<> mp;
    mp.addPatch(surf.clone());
    return mp;
}

// Split the single patch (u in [0,4]) into 4 patches at u = 1,2,3
static gsMultiPatch<> makeFourPatches(real_t angleRad)
{
    gsMultiPatch<> single = makeConnectedBoundarySinglePatch(angleRad);
    gsTensorBSpline<2, real_t> surf = static_cast<gsTensorBSpline<2, real_t>&>(single.patch(0));

    gsTensorBSpline<2, real_t> p0, pRest;
    surf.splitAt(0, 1.0, p0, pRest);
    gsTensorBSpline<2, real_t> p1, pRest2;
    pRest.splitAt(0, 2.0, p1, pRest2);
    gsTensorBSpline<2, real_t> p2, p3;
    pRest2.splitAt(0, 3.0, p2, p3);

    gsMultiPatch<> mp;
    mp.addPatch(p0.clone());
    mp.addPatch(p1.clone());
    mp.addPatch(p2.clone());
    mp.addPatch(p3.clone());
    mp.computeTopology();

    // Mark both sides as boundaries (outer: north, inner: south)
    for (index_t i = 0; i < mp.nPatches(); ++i) {
        mp.addBoundary(patchSide(i, boundary::north));
        mp.addBoundary(patchSide(i, boundary::south));
    }

    return mp;
}

// Insert extra knots only near boundary::south (v ≈ 0) to increase DOFs along the inner interface
static void applySouthRefinement(gsMultiPatch<>& mp)
{
    const std::vector<real_t> kvSouthWest = {0.02, 0.05, 0.1};

    for (index_t p = 0; p < mp.nPatches(); ++p)
    {
        for (real_t kv : kvSouthWest)
            mp.patch(p).insertKnot(kv, /*dir=*/1, /*mult=*/1); // dir=1 -> v-direction
    }
}

// Classify which outer-side (north) of a patch corresponds to physical left/right/top/bottom
// Returns one of: 0=left, 1=right, 2=bottom, 3=top
static int classifyOuterSide(const gsGeometry<>& geo)
{
    // Sample midpoints along v=1 boundary
    const int N = 5;
    real_t xsum = 0, ysum = 0;
    for (int i = 0; i < N; ++i)
    {
        real_t u = (i + 0.5) / N * 1.0; // u in [0,1]
        gsMatrix<> uv(2,1); uv << u, 1.0;
        gsMatrix<> pt = geo.eval(uv);
        xsum += pt(0,0);
        ysum += pt(1,0);
    }
    real_t xavg = xsum / N;
    real_t yavg = ysum / N;
    // Domain is approximately [-1,2] x [-1,2]
    if (xavg < -0.7) return 0; // left
    if (xavg >  1.7) return 1; // right
    if (yavg < -0.7) return 2; // bottom
    return 3;                  // top
}

int main(int argc, char* argv[])
{
    gsInfo << "4-patch rotating-square with gsINSALE (re-param every step)\n";

    // Physical and numerical parameters
    real_t meanVelocity = 2.0; // inlet mean speed
    real_t Re = 200.0;         // Reynolds number
    real_t L = 1.0;            // length scale
    real_t nu = meanVelocity * L / Re; // kinematic viscosity

    // Time control
    real_t timeSpan = 1.0;
    real_t timeStep = 0.02;
    real_t maxAngleDeg = 20.0; // rotation amplitude

    // Discretization
    int    numRefine = 3;
    bool   useMeshOpt = true;   // apply BarrierPatch each step on the 4-patch domain
    int    outEvery   = 10;

    gsCmdLine cmd("4-patch ALE with gsINSAle");
    cmd.addReal("t","time","Simulation time span", timeSpan);
    cmd.addReal("s","step","Time step size", timeStep);
    cmd.addReal("a","angle","Max rotation angle (deg)", maxAngleDeg);
    cmd.addInt ("r","refine","Uniform refinements", numRefine);
    cmd.addInt ("o","output","Output every N steps", outEvery);
    cmd.addSwitch("m","meshopt","Use BarrierPatch optimization each step", useMeshOpt);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Geometry at t=0: 4 patches
    gsMultiPatch<> fluidDomain = makeFourPatches(0.0);
    for (int i = 0; i < numRefine; ++i) fluidDomain.uniformRefine();
    // Add local DOFs only near south boundary (inner interface)
    applySouthRefinement(fluidDomain);

    // Build Taylor–Hood spaces: elevate for velocity only
    gsMultiPatch<> pressureDomain = fluidDomain; // copy before elevating
    fluidDomain.degreeElevate(1);
    gsMultiBasis<> basisU(fluidDomain);
    gsMultiBasis<> basisP(pressureDomain);

    // Zero vector (0,0) for body force and wall BCs
    gsConstantFunction<> zeroVec(0.0, 0.0, 2);

    // Prepare inlet (U,0) function (lifetime: whole main)
    class InletVec : public gsFunction<>
    {
        real_t U;
    public:
        explicit InletVec(real_t u):U(u){}
        short_t domainDim() const { return 2; }
        short_t targetDim() const { return 2; }
        void eval_into(const gsMatrix<>& x, gsMatrix<>& y) const {
            y.setZero(2, x.cols());
            for (index_t i=0;i<x.cols();++i){ y(0,i)=U; y(1,i)=0.0; }
        }
        gsFunction<>::uPtr clone() const { return memory::make_unique(new InletVec(*this)); }
    } inlet(meanVelocity);
    // Boundary mesh velocity computed by difference of current and previous geometry
    class BoundaryMeshVel : public gsFunction<>
    {
        const gsGeometry<>* m_cur;
        const gsGeometry<>* m_prev;
        real_t m_dt;
    public:
        BoundaryMeshVel(const gsGeometry<>* cur, const gsGeometry<>* prev, real_t dt)
        : m_cur(cur), m_prev(prev), m_dt(dt) {}
        void setCurrent(const gsGeometry<>* cur){ m_cur = cur; }
        void setPrevious(const gsGeometry<>* prev){ m_prev = prev; }
        void setDt(real_t dt){ m_dt = dt; }
        short_t domainDim() const { return 2; }
        short_t targetDim() const { return 2; }
        void eval_into(const gsMatrix<>& x, gsMatrix<>& y) const {
            const index_t N = x.cols();
            y.resize(2,N);
            for (index_t i=0;i<N;++i)
            {
                gsMatrix<> phys(2,1); phys.col(0) = x.col(i);
                gsMatrix<> uv; m_cur->invertPoints(phys, uv);
                gsMatrix<> xold = m_prev->eval(uv);
                y.col(i) = (phys.col(0) - xold) / m_dt;
            }
        }
        gsFunction<>::uPtr clone() const { return memory::make_unique(new BoundaryMeshVel(*this)); }
    };

    // prev geometry for boundary velocity at t=0 equals current (zero mesh velocity)
    gsMultiPatch<> prevPatches = fluidDomain;

    // Build fixed-numbering BCs once
    std::vector<std::unique_ptr<BoundaryMeshVel>> northFuncs;
    northFuncs.reserve(fluidDomain.nPatches());
    gsBoundaryConditions<> bcInfo;
    
    // North boundaries: ALL patches have ALE no-slip boundary conditions
    for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
    {
        northFuncs.emplace_back(new BoundaryMeshVel(&fluidDomain.patch(p), &prevPatches.patch(p), timeStep));
        bcInfo.addCondition(p, boundary::north, condition_type::dirichlet, northFuncs.back().get(), 0);
    }
    
    // South boundaries: patch 3 is inflow, patch 1 is outflow, others are walls
    for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
    {
        if (p == 3)
            bcInfo.addCondition(p, boundary::south, condition_type::dirichlet, &inlet, 0); // inflow
        else if (p == 1)
        {
            // outflow: do-nothing (no BC added)
        }
        else
            bcInfo.addCondition(p, boundary::south, condition_type::dirichlet, &zeroVec, 0); // walls
    }
    bcInfo.addCornerValue(boundary::corner::southwest, 0.0, 0, 1);

    // PDE and solver parameters (initial)
    std::vector<gsMultiBasis<>> bases; bases.push_back(basisU); bases.push_back(basisP);
    gsNavStokesPde<real_t> nsPde0(fluidDomain, bcInfo, &zeroVec, nu);
    gsFlowSolverParams<real_t> params0(nsPde0, bases);
    params0.options().setReal("timeStep", timeStep);
    params0.options().setInt ("nonlin.maxIt", 50);
    params0.options().setReal("nonlin.tol", 1e-6);
    params0.options().setInt ("lin.maxIt", 200);
    params0.options().setReal("lin.tol", 1e-8);
    params0.options().setSwitch("quiet", true);

    gsINSSolverUnsteadyALE<> solver(memory::make_shared_not_owned(&params0));
    solver.initialize();
    solver.setALEActive(true);

    // Output collections
    gsParaviewCollection collV("ale4p_velocity");
    gsParaviewCollection collP("ale4p_pressure");
    gsParaviewCollection collM("ale4p_mesh");

    index_t nSteps = math::ceil(timeSpan / timeStep);
    for (index_t step = 0; step <= nSteps; ++step)
    {
        real_t t = step * timeStep;
        real_t angleDeg = maxAngleDeg * std::sin(2 * M_PI * t / timeSpan);
        real_t angleRad = angleDeg * M_PI / 180.0;

        // Rebuild 4-patch geometry with current rotation
        gsMultiPatch<> newGeo = makeFourPatches(angleRad);
        for (int i = 0; i < numRefine; ++i) newGeo.uniformRefine();
        // Apply the same south-only refinement pattern to keep basis/geometry consistent
        applySouthRefinement(newGeo);
        newGeo.degreeElevate(2); // keep degree consistent with initial elevate(2)

        if (useMeshOpt)
        {
            gsBarrierPatch<2, real_t> opt(newGeo, false);
            opt.options().setInt("Verbose", 0);
            opt.options().setInt("ParamMethod", 0);
            opt.options().setInt("AAPreconditionType", 0);
            opt.compute();
            newGeo = opt.result();
        }

        // Update solver geometry in-place
        gsMultiPatch<>& patches = const_cast<gsMultiPatch<>&>(solver.getAssembler()->getPatches());
        patches = newGeo;

        // Update north boundary velocity functors with current and previous geometries
        for (index_t p = 0; p < patches.nPatches(); ++p)
        {
            northFuncs[p]->setCurrent(&patches.patch(p));
            northFuncs[p]->setPrevious(&prevPatches.patch(p));
            northFuncs[p]->setDt(timeStep);
        }

        // Solve
        if (step > 0)
            solver.nextIteration();
        else
            solver.solveStokes();

        // Update prev geometry
        prevPatches = patches;

        // Output
        if (step % outEvery == 0)
        {
            gsField<> vField = solver.constructSolution(0);
            gsField<> pField = solver.constructSolution(1);
            std::string vs = "ale4p_velocity_" + std::to_string(step);
            std::string ps = "ale4p_pressure_" + std::to_string(step);
            std::string ms = "ale4p_mesh_" + std::to_string(step);
            gsWriteParaview(vField, vs, 150, false);
            gsWriteParaview(pField, ps, 150, false);
            gsWriteParaview(patches, ms, 150, false);
            for (size_t p = 0; p < patches.nPatches(); ++p)
            {
                collV.addPart(vs + std::to_string(p) + ".vts", t, std::to_string(p));
                collP.addPart(ps + std::to_string(p) + ".vts", t, std::to_string(p));
                collM.addPart(ms + std::to_string(p) + ".vts", t, std::to_string(p));
            }
        }
    }

    collV.save();
    collP.save();
    collM.save();
    gsInfo << "Done.\n";
    return 0;
}
