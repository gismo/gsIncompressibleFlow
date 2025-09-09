/** @file perpendicular_flap_fluid_precice.cpp
    @brief Fluid participant for a preCICE FSI example (perpendicular flap style)

    - Uses ALE Navier–Stokes solver in G+Smo
    - Coupling via preCICE: Fluid writes traction (vector "Stress"), reads structure displacement (vector "Displacement")
    - Geometry: reuses flapping-beam flow domain (Turek–Hron style) from optional/gsElasticity/filedata
    - Interface detection: picks boundary sides not on the outer box as the fluid–structure interface

    This file is part of the G+Smo library under MPL-2.0.
*/

#include <gismo.h>

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSSolverALE.h>

#include <gsPreCICE/gsPreCICE.h>
#include <gsNurbs/gsNurbsCreator.h>

using namespace gismo;

namespace
{
// Construct a 2D fluid domain that matches the preCICE "perpendicular flap" setup:
// - Global channel: x in [-3, 3], y in [0, 4]
// - Flap (solid) occupies x in [-0.05, 0.05], y in [0, 1]; this is a "hole" in the fluid
// The fluid domain is built as a 5-patch MultiPatch:
//   0: Lower-left     [-3,   -0.05] x [0, 1]
//   1: Lower-right    [ 0.05,  3   ] x [0, 1]
//   2: Upper-left     [-3,   -0.05] x [1, 4]
//   3: Upper-right    [ 0.05,  3   ] x [1, 4]
//   4: Middle-top     [-0.05, 0.05] x [1, 4]   (its bottom y=1 is the flap top interface)
static gsMultiPatch<> makePerpendicularFlapFluidDomain()
{
    const real_t xl = -3.0, xr = 3.0;
    const real_t yb =  0.0, yt = 4.0;
    const real_t fxL = -0.05, fxR = 0.05; // flap thickness 0.1
    const real_t fyT =  1.0;              // flap height 1.0

    gsMultiPatch<> mp;
    // degree-2 rectangles (3x3 ctrl pts)
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(xl, yb, fxL, fyT)); // 0
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(fxR, yb, xr, fyT)); // 1
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(xl, fyT, fxL, yt)); // 2
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(fxR, fyT, xr, yt)); // 3
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(fxL, fyT, fxR, yt)); // 4 (middle top)

    mp.computeTopology();
    return mp;
}
// Helper: classify boundary sides into external box vs internal (FSI interface)
static void detectInterfaceSides(
    const gsMultiPatch<>& patches,
    std::vector<patchSide>& interfaceSides,
    std::vector<patchSide>& externalSides)
{
    // Compute overall bounding box of the domain (from control points)
    real_t xmin =  std::numeric_limits<real_t>::infinity();
    real_t xmax = -std::numeric_limits<real_t>::infinity();
    real_t ymin =  std::numeric_limits<real_t>::infinity();
    real_t ymax = -std::numeric_limits<real_t>::infinity();
    for (index_t p = 0; p < patches.nPatches(); ++p)
    {
        const gsMatrix<>& C = patches.patch(p).coefs();
        for (index_t i = 0; i < C.rows(); ++i)
        {
            xmin = std::min(xmin, C(i,0)); xmax = std::max(xmax, C(i,0));
            ymin = std::min(ymin, C(i,1)); ymax = std::max(ymax, C(i,1));
        }
    }

    const real_t tol = 1e-8;

    // Iterate over declared boundary sides
    for (gsMultiPatch<>::const_biterator it = patches.bBegin(); it != patches.bEnd(); ++it)
    {
        const index_t p = it->patch;
        const boxSide  s = it->side();

        // Sample a few physical points on this boundary side
        const int ns = 9;
        gsMatrix<> uv(2, ns);
        for (int k = 0; k < ns; ++k)
        {
            const real_t t = static_cast<real_t>(k) / static_cast<real_t>(ns - 1);
            if (s == boundary::west)  { uv(0,k) = 0.0; uv(1,k) = t; }
            if (s == boundary::east)  { uv(0,k) = 1.0; uv(1,k) = t; }
            if (s == boundary::south) { uv(0,k) = t;  uv(1,k) = 0.0; }
            if (s == boundary::north) { uv(0,k) = t;  uv(1,k) = 1.0; }
        }
        gsMatrix<> XY; patches.patch(p).eval_into(uv, XY);

        auto nearVal = [&](int dim, real_t val)
        {
            real_t a = 0.0; for (int k = 0; k < ns; ++k) a += std::abs(XY(dim,k) - val);
            a /= static_cast<real_t>(ns);
            return a < tol;
        };

        const bool onLeft   = nearVal(0, xmin);
        const bool onRight  = nearVal(0, xmax);
        const bool onBottom = nearVal(1, ymin);
        const bool onTop    = nearVal(1, ymax);

        if (onLeft || onRight || onBottom || onTop)
            externalSides.emplace_back(p, s);
        else
            interfaceSides.emplace_back(p, s);
    }
}

// Build the list of interface points using boundary anchors on the detected sides.
// Also build a mapping from these points to velocity DOF indices for ALE displacement application.
static void buildInterfaceAnchors(
    const gsMultiPatch<>& patches,
    const gsMultiBasis<>& velBasis,
    const std::vector<patchSide>& interfaceSides,
    gsMatrix<real_t>& X_if,                 // 2 x N: interface point coordinates (physical)
    std::vector<std::pair<index_t,index_t>>& dofLoci, // (patchID, coefIdx) for each interface point
    std::vector<std::pair<index_t, boxSide>>& lociSides, // (patchID, side) parallel to X_if columns
    gsMatrix<real_t>& uv_if                 // 2 x N: parametric coordinates for each interface point
)
{
    std::vector<gsVector<real_t>> coords;
    std::vector<std::pair<index_t,index_t>> dofs;
    std::vector<std::pair<index_t, boxSide>> sides;
    std::vector<gsVector<real_t>> uvs;

    for (const auto& ps : interfaceSides)
    {
        const index_t p = ps.patch;
        const boxSide  s = ps.side();

        // Boundary DOF indices on this side (in patch coefficient numbering)
        gsMatrix<index_t> bIdx = velBasis.basis(p).boundary(s);

        // Boundary 1D basis anchors on this side (param values in 1D)
        typename gsBasis<>::uPtr bndBasis = velBasis.basis(p).boundaryBasis(s);
        gsMatrix<> anchors = bndBasis->anchors(); // 1 x n

        GISMO_ASSERT(bIdx.size() == anchors.cols(), "Mismatch between DOF count and anchor count on a side");

        for (index_t i = 0; i < bIdx.size(); ++i)
        {
            // Parametric 2D coordinates along this side
            gsVector<> uv(2);
            const real_t t = anchors(0, i);
            if (s == boundary::west)  { uv[0] = 0.0; uv[1] = t; }
            if (s == boundary::east)  { uv[0] = 1.0; uv[1] = t; }
            if (s == boundary::south) { uv[0] = t;  uv[1] = 0.0; }
            if (s == boundary::north) { uv[0] = t;  uv[1] = 1.0; }

            // Physical coordinate at this param location
            gsMatrix<> uvM(2,1); uvM.col(0) = uv;
            gsMatrix<> P; patches.patch(p).eval_into(uvM, P); // 2 x 1

            coords.push_back(P.col(0));
            dofs.emplace_back(p, bIdx(i,0));
            sides.emplace_back(p, s);
            uvs.push_back(uv);
        }
    }

    // Deduplicate corner points (and any accidental duplicates)
    const real_t tol = 1e-12;
    std::vector<gsVector<real_t>> uniqueCoords;
    std::vector<gsVector<real_t>> uniqueUvs;
    std::vector<std::pair<index_t,index_t>> uniqueDofs;
    std::vector<std::pair<index_t, boxSide>> uniqueSides;

    auto isClose = [&](const gsVector<real_t>& a, const gsVector<real_t>& b)->bool
    { return (a - b).template lpNorm<gsEigen::Infinity>() < tol; };

    for (index_t k = 0; k < static_cast<index_t>(coords.size()); ++k)
    {
        bool dup = false;
        for (index_t m = 0; m < static_cast<index_t>(uniqueCoords.size()); ++m)
        {
            if (isClose(coords[k], uniqueCoords[m])) { dup = true; break; }
        }
        if (!dup)
        {
            uniqueCoords.push_back(coords[k]);
            uniqueUvs.push_back(uvs[k]);
            uniqueDofs.push_back(dofs[k]);
            uniqueSides.push_back(sides[k]);
        }
    }

    const index_t N = static_cast<index_t>(uniqueCoords.size());
    X_if.resize(2, N);
    uv_if.resize(2, N);
    dofLoci.resize(N);
    lociSides.resize(N);
    for (index_t k = 0; k < N; ++k)
    {
        X_if.col(k) = uniqueCoords[k];
        uv_if.col(k) = uniqueUvs[k];
        dofLoci[k] = uniqueDofs[k];
        lociSides[k] = uniqueSides[k];
    }
}

// Compute unit outer normal of the fluid domain at parametric point uv on side s, for patch p
static gsVector<real_t> boundaryUnitNormal(
    const gsMultiPatch<>& patches,
    index_t p, boxSide s, const gsVector<real_t>& uv)
{
    gsMapData<> md; md.flags = NEED_OUTER_NORMAL; md.patchId = p; { gsMatrix<> uvM(2,1); uvM.col(0) = uv; md.points = uvM; } md.side = s;
    patches.patch(p).computeMap(md);
    gsVector<> n = md.outNormal(0);
    real_t ln = n.norm();
    if (ln > 0) n /= ln;
    return n;
}

// Assemble traction = -p * n (pressure-only by default) at interface anchor points
static void computeFluidTraction(
    const gsField<>& velField,
    const gsField<>& presField,
    const gsMultiPatch<>& patches,
    const gsMatrix<real_t>& uv_if,
    const std::vector<std::pair<index_t, boxSide>>& lociSides,
    gsMatrix<real_t>& traction,
    real_t mu,
    bool includeViscous)
{
    const index_t N = uv_if.cols();
    traction.resize(2, N);

    // Evaluate pressure values at the parametric points (patch-wise)
    // We need to group by patch to call value(..., patchID)
    // For simplicity, loop point-wise (N typically modest).
    for (index_t k = 0; k < N; ++k)
    {
        const index_t p = lociSides[k].first;
        const boxSide  s = lociSides[k].second;
        const gsVector<> uvk = uv_if.col(k);

        // pressure at (uv) on patch p
        gsMatrix<> uvM(2,1); uvM.col(0) = uvk;
        const gsMatrix<> pk = presField.value(uvM, p); // 1 x 1

        // unit outer normal of fluid boundary
        gsVector<> n = boundaryUnitNormal(patches, p, s, uvk);

        // pressure traction
        gsVector<> t = -pk(0,0) * n;

        if (includeViscous)
        {
            // Optional: add viscous contribution tau*n ≈ mu * (∇u + ∇u^T) n
            // NOTE: evaluating ∇u in physical coordinates requires more setup (expr evaluator). Omitted for now.
        }

        traction.col(k) = t;
    }
}

} // namespace

int main(int argc, char* argv[])
{
    // ---------------------- Options ----------------------
    std::string preciceConfig;
    std::string geomPath = "auto"; // default: build perpendicular-flap fluid domain
    real_t viscosity = 1.0;   // kinematic viscosity (m^2/s)
    real_t totalTime = 5.0;
    real_t dt = 0.01;
    int numRefine = 0;
    bool plot = false;
    bool pressureOnly = true; // write traction = -p n

    gsCmdLine cmd("Fluid participant (ALE) for perpendicular-flap FSI via preCICE");
    cmd.addString("c","config","preCICE configuration file", preciceConfig);
    cmd.addString("g","geom","Fluid geometry XML (flappingBeam_flow.xml)", geomPath);
    cmd.addReal  ("v","nu","Kinematic viscosity", viscosity);
    cmd.addReal  ("t","time","Total simulation time", totalTime);
    cmd.addReal  ("s","step","Time step", dt);
    cmd.addInt   ("r","refine","Uniform h-refinements", numRefine);
    cmd.addSwitch("plot","Write Paraview outputs", plot);
    cmd.addSwitch("P","pressureOnly","Only pressure traction (no viscous)", pressureOnly);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Try to locate a default preCICE config if not provided
    if (preciceConfig.empty())
    {
        auto fileExists = [](const std::string &p)->bool{
            std::ifstream ifs(p.c_str());
            return static_cast<bool>(ifs);
        };
        const char* candidates[] = {
            "../optional/gsPreCICE/examples/perpendicular-flap/precice-config.xml",
            "../../optional/gsPreCICE/examples/perpendicular-flap/precice-config.xml",
            "../../../optional/gsPreCICE/examples/perpendicular-flap/precice-config.xml",
            "optional/gsPreCICE/examples/perpendicular-flap/precice-config.xml",
            "precice-config.xml"
        };
        for (const char* c : candidates)
        {
            if (fileExists(c)) { preciceConfig = c; break; }
        }
    }

    if (preciceConfig.empty())
    {
        gsWarn << "No preCICE config passed. Use -c <precice-config.xml>.\n";
        gsWarn << "Tried common locations but none found. Aborting.\n";
        return -1;
    }

    gsInfo << "[FSI-Fluid] Geometry: " << geomPath << "\n";
    gsInfo << "[FSI-Fluid] preCICE config: " << preciceConfig << "\n";

    // ---------------------- Geometry & spaces ----------------------
    gsMultiPatch<> fluidDomain;
    if (geomPath == "auto")
    {
        gsInfo << "[FSI-Fluid] Building perpendicular-flap fluid domain (5 patches)...\n";
        fluidDomain = makePerpendicularFlapFluidDomain();
    }
    else
    {
        GISMO_ASSERT(gsFileManager::fileExists(geomPath), "Geometry file not found: " << geomPath);
        gsReadFile<>(geomPath, fluidDomain);
    }

    gsMultiBasis<> basisVelocity(fluidDomain);
    gsMultiBasis<> basisPressure(fluidDomain);

    // Taylor–Hood: elevate deg of velocity basis
    basisVelocity.setDegree(basisVelocity.maxCwiseDegree() + 1);
    for (int r = 0; r < numRefine; ++r)
    {
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
    }

    // Detect interface sides (internal boundaries) vs external box boundaries
    std::vector<patchSide> interfaceSides, externalSides;
    detectInterfaceSides(fluidDomain, interfaceSides, externalSides);
    // Explicitly restrict FSI interface to flap contour: (0,east), (4,south), (1,west)
    interfaceSides.clear();
    interfaceSides.emplace_back(0, boundary::east);
    interfaceSides.emplace_back(4, boundary::south);
    interfaceSides.emplace_back(1, boundary::west);
    gsInfo << "[FSI-Fluid] Using interface sides: (0,east), (4,south), (1,west); externals: " << externalSides.size() << "\n";

    // ---------------------- Boundary conditions ----------------------
    gsBoundaryConditions<> bcInfo;
    // Inlet parabolic at left boundary (x = xmin). Detect using externalSides + sampling
    // We will set: inlet at minimum x, no-slip at top/bottom, natural at outlet (x = xmax).
    // Interface sides: no Dirichlet; handled via ALE moving mesh.

    // Determine bbox again for classification
    real_t xmin =  std::numeric_limits<real_t>::infinity();
    real_t xmax = -std::numeric_limits<real_t>::infinity();
    real_t ymin =  std::numeric_limits<real_t>::infinity();
    real_t ymax = -std::numeric_limits<real_t>::infinity();
    for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
    {
        const gsMatrix<>& C = fluidDomain.patch(p).coefs();
        for (index_t i = 0; i < C.rows(); ++i)
        {
            xmin = std::min(xmin, C(i,0)); xmax = std::max(xmax, C(i,0));
            ymin = std::min(ymin, C(i,1)); ymax = std::max(ymax, C(i,1));
        }
    }

    // Left boundary: constant inflow 10 m/s in x-direction
    gsConstantFunction<> inletVec(10.0, 0.0, 2);
    gsConstantFunction<> zeroVel(0.0, 0.0, 2);

    // Tag external sides
    const real_t tol = 1e-8;
    for (const auto& ps : externalSides)
    {
        const index_t p = ps.patch; const boxSide s = ps.side();
        // Sample a few points to classify
        const int ns = 7;
        gsMatrix<> uv(2, ns);
        for (int k = 0; k < ns; ++k)
        {
            const real_t t = static_cast<real_t>(k) / static_cast<real_t>(ns - 1);
            if (s == boundary::west)  { uv(0,k) = 0.0; uv(1,k) = t; }
            if (s == boundary::east)  { uv(0,k) = 1.0; uv(1,k) = t; }
            if (s == boundary::south) { uv(0,k) = t;  uv(1,k) = 0.0; }
            if (s == boundary::north) { uv(0,k) = t;  uv(1,k) = 1.0; }
        }
        gsMatrix<> XY; fluidDomain.patch(p).eval_into(uv, XY);

        auto nearVal = [&](int dim, real_t val)
        {
            real_t a = 0.0; for (int k = 0; k < ns; ++k) a += std::abs(XY(dim,k) - val);
            a /= static_cast<real_t>(ns);
            return a < tol;
        };

        const bool onLeft   = nearVal(0, xmin);
        const bool onRight  = nearVal(0, xmax);
        const bool onBottom = nearVal(1, ymin);
        const bool onTop    = nearVal(1, ymax);

        if (onLeft)
        {
            // Inlet Dirichlet on velocity (both components)
            bcInfo.addCondition(p, s, condition_type::dirichlet, &inletVec, 0);
        }
        else if (onBottom || onTop)
        {
            // No-slip walls
            bcInfo.addCondition(p, s, condition_type::dirichlet, &zeroVel, 0);
        }
        // onRight: natural outlet (do nothing)
    }

    // ---------------------- PDE & solver ----------------------
    gsConstantFunction<> fZero(0.0, 0.0, 2);
    gsNavStokesPde<real_t> nsPde(fluidDomain, bcInfo, &fZero, viscosity);

    std::vector<gsMultiBasis<>> spaces;
    spaces.push_back(basisVelocity);
    spaces.push_back(basisPressure);

    gsFlowSolverParams<real_t> params(nsPde, spaces);
    params.options().setReal("timeStep", dt);
    params.options().setInt ("nonlin.maxIt", 50);
    params.options().setReal("nonlin.tol", 1e-6);
    params.options().setInt ("lin.maxIt", 250);
    params.options().setReal("lin.tol", 1e-8);
    params.options().setSwitch("quiet", false);

    gsINSSolverUnsteadyALE<> solver(memory::make_shared_not_owned(&params));
    solver.initialize();
    solver.setALEActive(true);
    // Use re-parameterization to smoothly extend boundary motion
    solver.setMeshOptimization(true);
    solver.getMeshOptOptions().setInt("Verbose", 0);
    solver.getMeshOptOptions().setInt("ParamMethod", 1);
    solver.getMeshOptOptions().setInt("AAPreconditionType", 0);

    // ---------------------- Coupling mesh (anchors on interface sides) ----------------------
    gsMatrix<> X_if, uv_if;
    std::vector<std::pair<index_t,index_t>> dofLoci; // (patch, coefIdx)
    std::vector<std::pair<index_t, boxSide>> lociSides;
    buildInterfaceAnchors(fluidDomain, basisVelocity, interfaceSides, X_if, dofLoci, lociSides, uv_if);

    gsPreCICE<real_t> participant("Fluid", preciceConfig);
    const std::string FluidMesh = "Fluid-Mesh";
    const std::string StressData = "Stress";
    const std::string DisplData  = "Displacement";

    gsVector<index_t> vertexIDs;
    participant.addMesh(FluidMesh, X_if, vertexIDs);

    real_t precice_dt = participant.initialize();
    gsInfo << "[FSI-Fluid] preCICE dt_max = " << precice_dt << "\n";

    // ---------------------- Time loop ----------------------
    const index_t udofs = solver.getAssembler()->getUdofs();
    const index_t dim   = fluidDomain.geoDim();
    GISMO_ASSERT(dim==2, "This example is 2D");

    gsMatrix<> Uchk, Pchk; // checkpoints
    real_t time = 0.0;
    index_t step = 0;

    // Output collections
    gsParaviewCollection pvU("fsi_flap_fluid_U");
    gsParaviewCollection pvP("fsi_flap_fluid_P");

    while (participant.isCouplingOngoing())
    {
        if (participant.requiresWritingCheckpoint())
        {
            Uchk = solver.solutionCoefs(0);
            Pchk = solver.solutionCoefs(1);
        }

        // 1) Read structure displacement at interface points
        gsMatrix<> disp_if; // 2 x N
        participant.readData(FluidMesh, DisplData, vertexIDs, disp_if);

        // 2) Build full mesh displacement vector for ALE (size = dim*udofs)
        gsMatrix<> meshDisp(dim * udofs, 1); meshDisp.setZero();

        const gsDofMapper& uMap = solver.getAssembler()->getMappers()[0];
        for (index_t k = 0; k < static_cast<index_t>(dofLoci.size()); ++k)
        {
            const index_t p = dofLoci[k].first;
            const index_t i = dofLoci[k].second;
            if (!uMap.is_free(i, p))
                continue; // skip fixed dof if any

            const index_t gi = uMap.index(i, p);
            meshDisp(gi, 0)          = disp_if(0, k);
            meshDisp(gi + udofs, 0)  = disp_if(1, k);
        }

        // 3) Update mesh for current time
        solver.updateMesh(meshDisp);

        // 4) Solve one fluid iteration (ALE)
        // Ensure fluid dt does not exceed preCICE dt
        const real_t dt_step = std::min(dt, precice_dt);
        params.options().setReal("timeStep", dt_step);
        solver.nextIteration();

        // 5) Evaluate traction at interface points and write to preCICE
        gsField<> velField = solver.constructSolution(0);
        gsField<> prsField = solver.constructSolution(1);

        gsMatrix<> traction;
        computeFluidTraction(velField, prsField, fluidDomain, uv_if, lociSides, traction, viscosity, !pressureOnly);

        participant.writeData(FluidMesh, StressData, vertexIDs, traction);

        // 6) Advance coupling
        precice_dt = participant.advance(dt_step);

        // 7) Handle checkpoints
        if (participant.requiresReadingCheckpoint())
        {
            solver.setSolutionCoefs(Uchk, 0);
            solver.setSolutionCoefs(Pchk, 1);
        }
        else
        {
            time += dt_step;
            ++step;

            if (plot && (step % 5 == 0))
            {
                gsField<> vel = velField;
                gsField<> pre = prsField;
                std::string suf = util::to_string(step);
                gsWriteParaview(vel, "fsi_flap_U_" + suf, 400, false);
                gsWriteParaview(pre, "fsi_flap_P_" + suf, 400, false);
                for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
                {
                    pvU.addPart("fsi_flap_U_" + suf + std::to_string(p) + ".vts", time, std::to_string(p));
                    pvP.addPart("fsi_flap_P_" + suf + std::to_string(p) + ".vts", time, std::to_string(p));
                }
            }
        }
    }

    if (plot)
    {
        pvU.save();
        pvP.save();
    }

    participant.finalize();
    gsInfo << "[FSI-Fluid] Completed.\n";
    return 0;
}
