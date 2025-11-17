/** @file perpendicular_flap_fluid_precice_remesh.cpp
    @brief Fluid participant for preCICE FSI coupling with ALE formulation

    Features:
    - ALE Navier-Stokes solver with moving mesh
    - preCICE coupling: fluid writes traction, reads structure displacement
    - 3-patch perpendicular flap geometry
    - Optional local mesh refinement near FSI interface
    - Mesh motion handled via boundary velocity conditions u = w

    This file is part of the G+Smo library under MPL-2.0.
*/

#include <gismo.h>

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSSolverALE.h>

#include <gsPreCICE/gsPreCICE.h>
#include <gsNurbs/gsNurbsCreator.h>
#include <gsAssembler/gsQuadrature.h>
#include <gsModeling/gsBarrierPatch.h>

using namespace gismo;

namespace
{
/// @brief Local h-refinement near FSI interface by inserting knots
/// @param mp Multi-patch geometry to refine
/// @param sides FSI interface sides
/// @param levels Number of refinement levels (>=1)
/// @param dist0 Initial parametric distance from interface (e.g., 0.10)
/// @param ratio Geometric ratio for subsequent levels (e.g., 0.5)
static void refineNearInterface(
    gsMultiPatch<>& mp,
    const std::vector<patchSide>& sides,
    int levels,
    real_t dist0,
    real_t ratio)
{
    if (levels <= 0)
        return;
    // Clamp parameters to sane ranges
    ratio = std::max<real_t>(std::min<real_t>(ratio, 0.99), 0.05);
    dist0 = std::max<real_t>(std::min<real_t>(dist0, 0.49), 1e-6);

    for (const auto& ps : sides)
    {
        const index_t p = ps.patch;
        const boxSide s = ps.side();

        // u-direction for west/east; v-direction for south/north
        const int dir = (s == boundary::west || s == boundary::east) ? 0 : 1;
        const bool towardsOne = (s == boundary::east || s == boundary::north);

        real_t d = dist0;
        for (int lv = 0; lv < levels; ++lv)
        {
            real_t pos = towardsOne ? (1.0 - d) : d;
            pos = std::max<real_t>(std::min<real_t>(pos, 1.0 - 1e-10), 1e-10);
            mp.patch(p).insertKnot(pos, dir, /*mult=*/1);
            d *= ratio;
        }
    }

    mp.computeTopology();
}

/// @brief Construct perpendicular flap fluid domain (3-patch)
/// Channel: x ∈ [-3,3], y ∈ [0,4]; Flap: x ∈ [-0.05,0.05], y ∈ [0,1]
static gsMultiPatch<> makePerpendicularFlapFluidDomain()
{
    const real_t xl = -3.0, xr = 3.0;
    const real_t yb =  0.0, yt = 4.0;
    const real_t fxL = -0.05, fxR = 0.05; // flap thickness 0.1
    const real_t fyT =  1.0;              // flap height 1.0

    gsMultiPatch<> mp;
    gsKnotVector<> uKnots(0, 1, 0, 2), vKnots(0, 1, 0, 2); // Bilinear patches

    // Patch 0: Left trapezoid
    gsMatrix<> coefs0(4, 2);
    coefs0 << xl,  yb, fxL, yb, xl,  yt, fxL, fyT;
    mp.addPatch(gsTensorBSpline<2>(uKnots, vKnots, coefs0));

    // Patch 1: Right trapezoid
    gsMatrix<> coefs1(4, 2);
    coefs1 << fxR, yb, xr,  yb, fxR, fyT, xr,  yt;
    mp.addPatch(gsTensorBSpline<2>(uKnots, vKnots, coefs1));

    // Patch 2: Top trapezoid
    gsMatrix<> coefs2(4, 2);
    coefs2 << fxL, fyT, fxR, fyT, xl,  yt, xr,  yt;
    mp.addPatch(gsTensorBSpline<2>(uKnots, vKnots, coefs2));

    mp.computeTopology();
    return mp;
}

//==============================================================================
// FSI Interface Coupling Functions
//==============================================================================

/// @brief Build FSI interface anchor points from velocity basis
static void buildInterfaceAnchors(
    const gsMultiPatch<>& patches,
    const gsMultiBasis<>& velBasis,
    const std::vector<patchSide>& fsiInterfaceSides,
    gsMatrix<real_t>& interfacePoints,
    std::vector<std::pair<index_t,index_t>>& dofLocations,
    std::vector<std::pair<index_t, boxSide>>& sideLocations,
    gsMatrix<real_t>& parametricCoords)
{
    std::vector<gsVector<real_t>> coords;
    std::vector<std::pair<index_t,index_t>> dofs;
    std::vector<std::pair<index_t, boxSide>> sides;
    std::vector<gsVector<real_t>> uvs;

    for (const auto& ps : fsiInterfaceSides)
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

    // Deduplicate corner points
    const real_t tol = 1e-12;
    auto isClose = [tol](const gsVector<real_t>& a, const gsVector<real_t>& b) {
        return (a - b).template lpNorm<gsEigen::Infinity>() < tol;
    };

    std::vector<gsVector<real_t>> uniqueCoords, uniqueUvs;
    std::vector<std::pair<index_t,index_t>> uniqueDofs;
    std::vector<std::pair<index_t, boxSide>> uniqueSides;

    for (size_t k = 0; k < coords.size(); ++k)
    {
        bool dup = false;
        for (size_t m = 0; m < uniqueCoords.size() && !dup; ++m)
            dup = isClose(coords[k], uniqueCoords[m]);

        if (!dup) {
            uniqueCoords.push_back(coords[k]);
            uniqueUvs.push_back(uvs[k]);
            uniqueDofs.push_back(dofs[k]);
            uniqueSides.push_back(sides[k]);
        }
    }

    const index_t N = static_cast<index_t>(uniqueCoords.size());
    interfacePoints.resize(2, N);
    parametricCoords.resize(2, N);
    dofLocations.resize(N);
    sideLocations.resize(N);
    for (index_t k = 0; k < N; ++k)
    {
        interfacePoints.col(k) = uniqueCoords[k];
        parametricCoords.col(k) = uniqueUvs[k];
        dofLocations[k] = uniqueDofs[k];
        sideLocations[k] = uniqueSides[k];
    }
}

/// @brief Build geometry control points on FSI interface (for IGA coupling)
static void buildGeomInterfaceCPs(
    const gsMultiPatch<>& patches,
    const std::vector<patchSide>& fsiInterfaceSides,
    gsMatrix<real_t>& cpPoints,
    std::vector<std::pair<index_t,index_t>>& cpDofLocations,
    std::vector<std::pair<index_t, boxSide>>& cpSideLocations,
    gsMatrix<real_t>& cpParametric)
{
    std::vector<gsVector<real_t>> coords;
    std::vector<std::pair<index_t,index_t>> dofs;
    std::vector<std::pair<index_t, boxSide>> sides;
    std::vector<gsVector<real_t>> uvs;

    for (const auto& ps : fsiInterfaceSides)
    {
        const index_t p = ps.patch;
        const boxSide  s = ps.side();

        // Boundary control point indices on this side (geometry basis)
        gsMatrix<index_t> bIdx = patches.patch(p).basis().boundary(s);

        // Boundary 1D anchors for geometry boundary basis
        typename gsBasis<>::uPtr bndBasis = patches.patch(p).basis().boundaryBasis(s);
        gsMatrix<> anchors = bndBasis->anchors(); // (1 x n)

        GISMO_ASSERT(bIdx.size() == anchors.cols(), "Mismatch between CP count and anchor count on a side");

        for (index_t i = 0; i < bIdx.size(); ++i)
        {
            gsVector<> uv(2);
            const real_t t = anchors(0, i);
            if (s == boundary::west)  { uv[0] = 0.0; uv[1] = t; }
            if (s == boundary::east)  { uv[0] = 1.0; uv[1] = t; }
            if (s == boundary::south) { uv[0] = t;  uv[1] = 0.0; }
            if (s == boundary::north) { uv[0] = t;  uv[1] = 1.0; }

            // Physical coordinate on boundary at this anchor
            gsMatrix<> uvM(2,1); uvM.col(0) = uv;
            gsMatrix<> P; patches.patch(p).eval_into(uvM, P); // 2x1

            coords.push_back(P.col(0));
            dofs.emplace_back(p, bIdx(i,0));
            sides.emplace_back(p, s);
            uvs.push_back(uv);
        }
    }

    // Deduplicate corner points
    const real_t tol = 1e-12;
    auto isClose = [tol](const gsVector<real_t>& a, const gsVector<real_t>& b) {
        return (a - b).template lpNorm<gsEigen::Infinity>() < tol;
    };

    std::vector<gsVector<real_t>> uniqueCoords, uniqueUvs;
    std::vector<std::pair<index_t,index_t>> uniqueDofs;
    std::vector<std::pair<index_t, boxSide>> uniqueSides;

    for (size_t k = 0; k < coords.size(); ++k) {
        bool dup = false;
        for (size_t m = 0; m < uniqueCoords.size() && !dup; ++m)
            dup = isClose(coords[k], uniqueCoords[m]);
        if (!dup) {
            uniqueCoords.push_back(coords[k]);
            uniqueDofs.push_back(dofs[k]);
            uniqueSides.push_back(sides[k]);
            uniqueUvs.push_back(uvs[k]);
        }
    }

    const index_t N = static_cast<index_t>(uniqueCoords.size());
    cpPoints.resize(2, N);
    cpParametric.resize(2, N);
    cpDofLocations.resize(N);
    cpSideLocations.resize(N);
    for (index_t k = 0; k < N; ++k)
    {
        cpPoints.col(k) = uniqueCoords[k];
        cpParametric.col(k) = uniqueUvs[k];
        cpDofLocations[k] = uniqueDofs[k];
        cpSideLocations[k] = uniqueSides[k];
    }
}

//==============================================================================
// Traction Computation Functions
//==============================================================================

/// @brief Compute ALE-consistent fluid traction on FSI interface
/// Piola stress: P = J * σ * F^{-T}, Traction: t = P * n_ref
static void computeFluidTractionALE(
    const gsField<>& velField,
    const gsField<>& presField,
    const gsField<>& meshDispField,           // ALE displacement field
    const gsMultiPatch<>& refPatches,         // reference (undeformed) geometry
    const gsMatrix<real_t>& parametricCoords, // boundary param coords on reference geometry
    const std::vector<std::pair<index_t, boxSide>>& sideLocations,
    gsMatrix<real_t>& traction,
    real_t rho,
    real_t nu,
    bool includeViscous)
{
    const index_t N = parametricCoords.cols();
    const index_t dim = 2;
    traction.resize(dim, N);

    // Identity
    const gsMatrix<real_t> I = gsMatrix<real_t>::Identity(dim, dim);

    for (index_t k = 0; k < N; ++k)
    {
        const index_t p = sideLocations[k].first;
        const boxSide  s = sideLocations[k].second;

        // Param location on boundary (reference geometry param domain)
        const gsVector<> uvk = parametricCoords.col(k);
        gsMatrix<> uvM(2,1); uvM.col(0) = uvk;

        // Reference geometry map data (need grad transform and normal on side)
        gsMapData<> mdGeo(NEED_GRAD_TRANSFORM | NEED_MEASURE);
        mdGeo.patchId = p; mdGeo.points = uvM; mdGeo.side = s;
        refPatches.patch(p).computeMap(mdGeo);

        // Unit normal in reference configuration
        gsVector<> nref = mdGeo.outNormal(0);
        real_t nlen = nref.norm();
        if (nlen > 0) nref /= nlen;

        // Pressure at this point (from fluid field)
        const gsMatrix<> pk = presField.value(uvM, p); // 1x1

        // Velocity gradient wrt reference physical coordinates: grad(u) = (du/dxi) * (dxi/dX)
        gsMapData<> mdVel(NEED_DERIV);
        mdVel.patchId = p; mdVel.points = uvM;
        velField.patches().patch(p).computeMap(mdVel);
        gsMatrix<> gradU = mdVel.jacobian(0) * mdGeo.jacobian(0).cramerInverse();

        // Displacement gradient wrt reference: Grad(d) = (dd/dxi) * (dxi/dX)
        gsMapData<> mdDisp(NEED_DERIV);
        mdDisp.patchId = p; mdDisp.points = uvM;
        meshDispField.patches().patch(p).computeMap(mdDisp);
        gsMatrix<> GradD = mdDisp.jacobian(0) * mdGeo.jacobian(0).cramerInverse();

        // Deformation gradient F and its inverse
        gsMatrix<> F = I + GradD;
        gsMatrix<> Finv = F.cramerInverse();
        real_t J = F.determinant();

        // Cauchy stress with ALE-consistent viscous term
        // Following gsFsiLoad convention: sigma = p*I - mu*(grad(u)*F^{-1} + F^{-T}*grad(u)^T)
        gsMatrix<> sigma = pk(0,0) * I;  // Pressure contribution (positive)
        if (includeViscous && rho * nu > 0)
        {
            const real_t muDyn = rho * nu; // dynamic viscosity
            sigma -= muDyn * (gradU * Finv + Finv.transpose() * gradU.transpose());
        }

        // Piola (1st) stress: P = J * sigma * F^{-T}
        gsMatrix<> P = J * sigma * Finv.transpose();

        // Traction transmitted to solid (per reference unit normal)
        // Following gsFsiLoad: traction = P * n_ref (no negative sign)
        traction.col(k) = P * nref;
    }
}

/// @brief Compute Eulerian fluid traction: t = -σ * n
static void computeFluidTractionEulerian(
    const gsField<>& velField,
    const gsField<>& presField,
    const gsMultiPatch<>& currentGeometry,
    const gsMatrix<real_t>& parametricCoords,
    const std::vector<std::pair<index_t, boxSide>>& sideLocations,
    gsMatrix<real_t>& traction,
    real_t rho,
    real_t nu,
    bool includeViscous)
{
    const index_t N = parametricCoords.cols();
    const index_t dim = 2;
    traction.resize(dim, N);

    const gsMatrix<real_t> I = gsMatrix<real_t>::Identity(dim, dim);

    for (index_t k = 0; k < N; ++k)
    {
        const index_t p = sideLocations[k].first;
        const boxSide s = sideLocations[k].second;

        const gsVector<> uvk = parametricCoords.col(k);
        gsMatrix<> uvM(2,1); uvM.col(0) = uvk;

        // Current geometry map data (for normal)
        gsMapData<> mdGeo(NEED_GRAD_TRANSFORM | NEED_MEASURE);
        mdGeo.patchId = p; mdGeo.points = uvM; mdGeo.side = s;
        currentGeometry.patch(p).computeMap(mdGeo);

        // Unit outer normal
        gsVector<> n = mdGeo.outNormal(0);
        real_t nlen = n.norm();
        if (nlen > 0) n /= nlen;

        // Pressure
        const gsMatrix<> pk = presField.value(uvM, p);

        // Velocity gradient in current configuration
        gsMapData<> mdVel(NEED_DERIV);
        mdVel.patchId = p; mdVel.points = uvM;
        velField.patches().patch(p).computeMap(mdVel);
        gsMatrix<> gradU = mdVel.jacobian(0) * mdGeo.jacobian(0).cramerInverse();

        // Cauchy stress: sigma = -p I + mu * (grad(u) + grad(u)^T)
        gsMatrix<> sigma = -pk(0,0) * I;
        if (includeViscous && rho * nu > 0)
        {
            const real_t muDyn = rho * nu;
            sigma += muDyn * (gradU + gradU.transpose());
        }

        // Traction transmitted TO solid: t_solid = -sigma_fluid * n_fluid
        // (Action-reaction: fluid stress acts on solid with opposite sign)
        traction.col(k) = -sigma * n;
    }
}

} // namespace

int main(int argc, char* argv[])
{
    //==========================================================================
    // Command-line Options
    //==========================================================================

    // Simulation parameters
    std::string preciceConfig;
    std::string geomPath = "auto";
    real_t viscosity = 1.0;      // kinematic viscosity nu (m^2/s)
    real_t density   = 1.0;      // fluid density rho (kg/m^3)
    real_t totalTime = 5.0;      // total simulation time (s)
    real_t dt = 0.01;            // time step (s)

    // Discretization parameters
    int numRefine = 0;           // uniform h-refinements
    int degreeElevate = 0;       // degree elevation
    int ifaceRefine = 0;         // local refinement levels near FSI interface
    real_t ifaceDist = 0.10;     // initial parametric distance from interface
    real_t ifaceRatio = 0.5;     // geometric ratio for local refinement

    // Output and coupling parameters
    bool plot = false;           // write Paraview output
    real_t plotEvery = 1.0;      // plot interval (seconds)
    bool pressureOnly = false;   // compute only pressure traction (no viscous)
    bool splineCoupling = false; // use IGA-based coupling on control points

    // Parse command line
    gsCmdLine cmd("Fluid participant (ALE) for perpendicular-flap FSI via preCICE");
    cmd.addString("c", "config",  "preCICE configuration file", preciceConfig);
    cmd.addString("g", "geom",    "Fluid geometry XML", geomPath);
    cmd.addReal  ("v", "nu",      "Kinematic viscosity", viscosity);
    cmd.addReal  ("R", "density", "Fluid density", density);
    cmd.addReal  ("t", "time",    "Total simulation time", totalTime);
    cmd.addReal  ("s", "step",    "Time step", dt);
    cmd.addInt   ("r", "refine",  "Uniform h-refinements", numRefine);
    cmd.addInt   ("p", "degree",  "Degree elevation", degreeElevate);
    cmd.addInt   ("i", "ifaceRefine", "Local refinement levels near FSI interface", ifaceRefine);
    cmd.addReal  ("d", "ifaceDist",   "Parametric distance for local refinement", ifaceDist);
    cmd.addReal  ("q", "ifaceRatio",  "Geometric ratio for local refinement", ifaceRatio);
    cmd.addSwitch("plot",       "Write Paraview outputs", plot);
    cmd.addReal  ("E", "plot-every", "Plot interval in seconds", plotEvery);
    cmd.addSwitch("P", "pressureOnly", "Only pressure traction (no viscous)", pressureOnly);
    cmd.addSwitch("igaCoupling", "Use IGA-based coupling on control points", splineCoupling);
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

    gsInfo << "[FSI-Fluid-ALE] Geometry: " << geomPath << "\n";
    gsInfo << "[FSI-Fluid-ALE] preCICE config: " << preciceConfig << "\n";

    //==========================================================================
    // Geometry Setup
    //==========================================================================

    gsMultiPatch<> fluidDomain;
    if (geomPath == "auto")
    {
        gsInfo << "[FSI-Fluid|NoRemesh] Building perpendicular-flap fluid domain (3 patches)...\n";
        fluidDomain = makePerpendicularFlapFluidDomain();
    }
    else
    {
        GISMO_ASSERT(gsFileManager::fileExists(geomPath), "Geometry file not found: " << geomPath);
        gsReadFile<>(geomPath, fluidDomain);
    }

    // Apply degree elevation if requested
    for (int d = 0; d < degreeElevate; ++d)
        fluidDomain.degreeElevate(1);

    // Optional: refine near the known FSI interface sides before creating bases
    if (ifaceRefine > 0)
    {
        std::vector<patchSide> refineSides;
        refineSides.emplace_back(0, boundary::east);   // left side of flap
        refineSides.emplace_back(2, boundary::south);  // top of flap
        refineSides.emplace_back(1, boundary::west);   // right side of flap

        gsInfo << "[FSI-Fluid|NoRemesh] Local interface refinement: levels="
               << ifaceRefine << ", dist0=" << ifaceDist << ", ratio=" << ifaceRatio << "\n";
        refineNearInterface(fluidDomain, refineSides, ifaceRefine, ifaceDist, ifaceRatio);
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

    // Directly specify FSI interface and outer boundary sides for 3-patch domain
    // FSI interface: left side, top, and right side of flap
    std::vector<patchSide> fsiInterfaceSides;
    fsiInterfaceSides.emplace_back(0, boundary::east);   // left side of flap
    fsiInterfaceSides.emplace_back(2, boundary::south);  // top of flap
    fsiInterfaceSides.emplace_back(1, boundary::west);   // right side of flap

    // Outer boundaries: inlet (left), outlet (right), walls (top/bottom)
    std::vector<patchSide> outerBoundarySides;
    outerBoundarySides.emplace_back(0, boundary::west);   // inlet (left wall)
    outerBoundarySides.emplace_back(0, boundary::south);  // bottom wall (left patch)
    outerBoundarySides.emplace_back(1, boundary::east);   // outlet (right wall)
    outerBoundarySides.emplace_back(1, boundary::south);  // bottom wall (right patch)

    outerBoundarySides.emplace_back(2, boundary::north);  // top wall


    gsInfo << "[FSI-Fluid] FSI interface: (0,east), (2,south), (1,west)\n";
    gsInfo << "[FSI-Fluid] Outer boundaries: " << outerBoundarySides.size() << " sides\n";

    // ---------------------- Boundary conditions ----------------------
    gsBoundaryConditions<> bcInfo;

    // Inlet: constant velocity 10 m/s in x-direction at left boundary
    gsConstantFunction<> inletVec(10.0, 0.0, 2);
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, &inletVec, 0);

    // No-slip walls: zero velocity at top and bottom boundaries
    gsConstantFunction<> zeroVel(0.0, 0.0, 2);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &zeroVel, 0);  // bottom left
    bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, &zeroVel, 0);  // bottom right

    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, &zeroVel, 0);  // top wall


    // Outlet: natural boundary condition at right boundary (patch 1, east) - no condition needed
    // FSI interface: will be handled by ALE boundary conditions

    // FSI Interface Kinematic Condition in ALE Formulation:
    // Attach u = w (mesh velocity) on FSI sides BEFORE building the PDE, so
    // the Dirichlet mapper includes these boundary DOFs.
    std::vector<std::shared_ptr< gsFSIMeshVelocityFunction<> >> ifaceMeshVelFuncs;
    ifaceMeshVelFuncs.reserve(fsiInterfaceSides.size());
    for (const auto& ps : fsiInterfaceSides)
    {
        auto f = std::make_shared< gsFSIMeshVelocityFunction<> >(
            &fluidDomain.patch(ps.patch), &fluidDomain.patch(ps.patch), dt);
        ifaceMeshVelFuncs.push_back(f);
        bcInfo.addCondition(ps.patch, ps.side(), condition_type::dirichlet, f.get(), 0);
    }
    gsInfo << "[FSI-Fluid] Applied Dirichlet BC u=w on "
           << fsiInterfaceSides.size() << " FSI interface sides\n";

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

    gsINSSolverUnsteadyALE<> solver(memory::make_shared_not_owned(&params));

    // Inform ALE/mesh optimizer about FSI interface sides (keep them as moving interface)
    solver.setFSIInterfaceSides(fsiInterfaceSides);

    solver.initialize();

    // NOTE: We do NOT use setMeshUpdateFunction() / updateMesh() approach.
    // Instead, we directly modify geometry control points (following ale_square_4patch_gsINSALE.cpp).
    // This is simpler and avoids issues with Dirichlet BC on velocity DOFs.

    solver.setALEActive(true);

    // Enable mesh optimization to smooth interior after boundary displacement
    solver.setMeshOptimization(true);
    solver.getMeshOptOptions().setInt("Verbose", 0);
    // ParamMethod: 0=PDEPatch, 1=PenaltyPatch, 2=PenaltyPatch2, 3=BarrierPatch, 4=VariationalHarmonicPatch
    solver.getMeshOptOptions().setInt("ParamMethod", 3);  // Use BarrierPatch method with HLBFGS
    solver.getMeshOptOptions().setInt("AAPreconditionType", 0);
    // Note: Using GLOBAL optimization mode (all patches together, interfaces FREE to move)
    // Modified in gsINSSolverALE.h: gsBarrierPatch<2, T> opt(patches, false) for global HLBFGS
    gsInfo << "[FSI-Fluid|NoRemesh] Using ALE with GLOBAL mesh optimization (HLBFGS, interfaces free)\n";

    // Store reference (undeformed) geometry for ALE-consistent traction
    gsMultiPatch<> refFluidGeometry = solver.getAssembler()->getPatches();
    // Previous step geometry (for boundary velocity computation)
    gsMultiPatch<> prevPatches = solver.getAssembler()->getPatches();

    // ---------------------- Build FSI coupling mesh ----------------------
    // Option 1: Use velocity basis anchors (current approach - for vertex-based coupling)
    // Option 2: Use geometry control points on boundary (for IGA-IGA coupling)

    gsMatrix<> interfacePoints, parametricCoords;
    std::vector<std::pair<index_t,index_t>> dofLocations; // (patch, coefIdx)
    std::vector<std::pair<index_t, boxSide>> sideLocations;

    // Build interface anchor points from velocity basis
    buildInterfaceAnchors(fluidDomain, basisVelocity, fsiInterfaceSides,
                          interfacePoints, dofLocations, sideLocations, parametricCoords);

    // Enhanced deduplication to avoid RBF matrix singularity
    const real_t tol_dedup = 1e-10;
    std::vector<index_t> uniqueIndices;
    uniqueIndices.reserve(interfacePoints.cols());

    for (index_t i = 0; i < interfacePoints.cols(); ++i)
    {
        bool isDuplicate = false;
        for (index_t j : uniqueIndices)
        {
            if ((interfacePoints.col(i) - interfacePoints.col(j)).template lpNorm<gsEigen::Infinity>() < tol_dedup)
            {
                isDuplicate = true;
                break;
            }
        }
        if (!isDuplicate)
        {
            uniqueIndices.push_back(i);
        }
    }

    // Create deduplicated data
    gsMatrix<> uniquePoints(2, uniqueIndices.size());
    gsMatrix<> uniqueParams(2, uniqueIndices.size());
    std::vector<std::pair<index_t,index_t>> uniqueDofLocs;
    std::vector<std::pair<index_t, boxSide>> uniqueSideLocs;

    for (index_t k = 0; k < static_cast<index_t>(uniqueIndices.size()); ++k)
    {
        index_t origIdx = uniqueIndices[k];
        uniquePoints.col(k) = interfacePoints.col(origIdx);
        uniqueParams.col(k) = parametricCoords.col(origIdx);
        uniqueDofLocs.push_back(dofLocations[origIdx]);
        uniqueSideLocs.push_back(sideLocations[origIdx]);
    }

    // Replace original data with deduplicated versions
    interfacePoints = uniquePoints;
    parametricCoords = uniqueParams;
    dofLocations = uniqueDofLocs;
    sideLocations = uniqueSideLocs;

    gsInfo << "[FSI-Fluid|NoRemesh] FSI interface points: " << interfacePoints.cols()
           << " (after deduplication)\n";

    // (already set above) Inform ALE/mesh optimizer about FSI interface sides

    gsPreCICE<real_t> participant("Fluid", preciceConfig);
    const std::string FluidMesh = "Fluid-Mesh";         // vertex/anchor-based
    const std::string FluidCPMesh = "Fluid-CP-Mesh";    // geometry control point-based (IGA)
    const std::string StressData = "Stress";
    const std::string DisplData  = "Displacement";

    gsVector<index_t> vertexIDs;
    gsVector<index_t> cpIDs;

    // If using spline/IGA coupling, build CP-based interface mesh
    gsMatrix<> cpPoints, cpParams;
    std::vector<std::pair<index_t,index_t>> cpDofLocs;
    std::vector<std::pair<index_t, boxSide>> cpSideLocs;
    if (splineCoupling)
    {
        buildGeomInterfaceCPs(fluidDomain, fsiInterfaceSides,
                              cpPoints, cpDofLocs, cpSideLocs, cpParams);
        participant.addMesh(FluidCPMesh, cpPoints, cpIDs);
        gsInfo << "[FSI-Fluid|NoRemesh] Using spline/IGA coupling: CP points = "
               << cpPoints.cols() << "\n";
    }
    else
    {
        participant.addMesh(FluidMesh, interfacePoints, vertexIDs);
    }

    real_t precice_dt = participant.initialize();
    gsInfo << "[FSI-Fluid|NoRemesh] preCICE dt_max = " << precice_dt << "\n";

    // ---------------------- Time loop ----------------------
    GISMO_ASSERT(fluidDomain.geoDim()==2, "This example is 2D");

    gsMatrix<> Uchk, Pchk; // checkpoints
    real_t time = 0.0;
    real_t nextPlotTime = plotEvery; // first plot at 1*plotEvery seconds
    index_t step = 0;

    // Output collections
    gsParaviewCollection pvU("fsi_flap_fluid_U_norem");
    gsParaviewCollection pvP("fsi_flap_fluid_P_norem");

    while (participant.isCouplingOngoing())
    {
        if (participant.requiresWritingCheckpoint())
        {
            Uchk = solver.solutionCoefs(0);
            Pchk = solver.solutionCoefs(1);
        }

        // 1) Read structure displacement at interface points
        gsMatrix<> disp_if; // 2 x N (vertex-based) or 2 x Nc (CP-based)
        if (splineCoupling)
            participant.readData(FluidCPMesh, DisplData, cpIDs, disp_if);
        else
            participant.readData(FluidMesh, DisplData, vertexIDs, disp_if);

        // 2) Apply interface displacement to mesh
        real_t maxDispMag = disp_if.template lpNorm<gsEigen::Infinity>();
        if (step % 20 == 0)
            gsInfo << "[FSI|Step " << step << "] Displacement: ||d||∞ = " << maxDispMag << "\n";

        // Create a temporary displacement field on reference geometry
        // by directly setting interface control point displacements
        gsMultiPatch<> dispField;
        for (index_t p = 0; p < refFluidGeometry.nPatches(); ++p)
        {
            gsMatrix<> dispCoefs = refFluidGeometry.patch(p).coefs(); // copy
            dispCoefs.setZero(); // will store displacement

            // Clone geometry structure
            dispField.addPatch(refFluidGeometry.patch(p).clone());
            dispField.patch(p).coefs() = dispCoefs;
        }

        int totalSet = 0;
        if (splineCoupling)
        {
            // Directly set boundary control point displacements from CP-coupling data
            GISMO_ASSERT(static_cast<index_t>(cpDofLocs.size()) == disp_if.cols(),
                         "Mismatch between CP locations and received displacement size");
            for (index_t k = 0; k < disp_if.cols(); ++k)
            {
                const index_t p = cpDofLocs[k].first;
                const index_t coefIdx = cpDofLocs[k].second;
                dispField.patch(p).coefs()(coefIdx, 0) = disp_if(0, k);
                dispField.patch(p).coefs()(coefIdx, 1) = disp_if(1, k);
                totalSet++;
            }
        }
        else
        {
            // Vertex-based: fit interface displacement to boundary control points (LS)
            // Group by (patch, side) to handle each boundary separately
            std::map<std::pair<index_t,boxSide>, std::vector<index_t>> sideToCols;
            for (index_t k = 0; k < parametricCoords.cols(); ++k)
            {
                sideToCols[sideLocations[k]].push_back(k);
            }

            for (const auto& entry : sideToCols)
            {
                const index_t p = entry.first.first;
                const boxSide s = entry.first.second;
                const std::vector<index_t>& sidePoints = entry.second;

                if (sidePoints.empty()) continue;

                // Get boundary control points
                gsMatrix<index_t> bIdx = refFluidGeometry.patch(p).basis().boundary(s);

                // Build collocation
                const index_t m = bIdx.size();
                const index_t n = sidePoints.size();

                gsMatrix<> quadUV(1, n);
                gsMatrix<> targetX(n, 1), targetY(n, 1);

                for (index_t j = 0; j < n; ++j)
                {
                    const index_t kcol = sidePoints[j];
                    targetX(j, 0) = disp_if(0, kcol);
                    targetY(j, 0) = disp_if(1, kcol);

                    const gsVector<> uvk = parametricCoords.col(kcol);
                    if (s == boundary::west || s == boundary::east) {
                        quadUV(0, j) = uvk[1];
                    } else {
                        quadUV(0, j) = uvk[0];
                    }
                }

                typename gsBasis<>::uPtr bndBasis = refFluidGeometry.patch(p).basis().boundaryBasis(s);
                gsMatrix<> N; bndBasis->eval_into(quadUV, N);

                gsMatrix<> NtN = N * N.transpose();
                gsMatrix<> Nttarget_x = N * targetX;
                gsMatrix<> Nttarget_y = N * targetY;

                gsVector<> cx = NtN.colPivHouseholderQr().solve(Nttarget_x);
                gsVector<> cy = NtN.colPivHouseholderQr().solve(Nttarget_y);

                for (index_t j = 0; j < m; ++j)
                {
                    const index_t coefIdx = bIdx(j, 0);
                    dispField.patch(p).coefs()(coefIdx, 0) = cx[j];
                    dispField.patch(p).coefs()(coefIdx, 1) = cy[j];
                    totalSet++;
                }
            }
        }

        // Apply displacement to geometry DIRECTLY (similar to ale_square_4patch_gsINSALE.cpp)
        // Key insight: preCICE sends CUMULATIVE displacement from reference,
        // so we need to reset to reference geometry first, then apply displacement
        gsMultiPatch<>& patches = const_cast<gsMultiPatch<>&>(solver.getAssembler()->getPatches());

        // Store geometry BEFORE deformation for comparison
        gsMultiPatch<> geoBeforeUpdate = patches;

        // Reset to reference geometry, then apply cumulative displacement
        for (index_t p = 0; p < patches.nPatches(); ++p)
        {
            const gsMatrix<>& dispCoefs = dispField.patch(p).coefs();
            gsMatrix<>& geoCoefs = patches.patch(p).coefs();
            const gsMatrix<>& refCoefs = refFluidGeometry.patch(p).coefs();

            GISMO_ASSERT(dispCoefs.rows() == geoCoefs.rows(),
                        "Displacement and geometry control point count mismatch");
            GISMO_ASSERT(refCoefs.rows() == geoCoefs.rows(),
                        "Reference and current geometry control point count mismatch");

            // Set geometry = reference + displacement
            for (index_t i = 0; i < geoCoefs.rows(); ++i)
            {
                geoCoefs(i, 0) = refCoefs(i, 0) + dispCoefs(i, 0);
                geoCoefs(i, 1) = refCoefs(i, 1) + dispCoefs(i, 1);
            }
        }

        const gsMultiPatch<>& curGeo = solver.getAssembler()->getPatches();

        // 3) Update FSI boundary velocity (u = w) functions with current/previous geometry
        const real_t dt_step = std::min(dt, precice_dt);
        for (index_t i = 0; i < static_cast<index_t>(fsiInterfaceSides.size()); ++i)
        {
            const index_t p = fsiInterfaceSides[i].patch;
            ifaceMeshVelFuncs[i]->setCurrent(&curGeo.patch(p));
            ifaceMeshVelFuncs[i]->setPrevious(&prevPatches.patch(p));
            ifaceMeshVelFuncs[i]->setDt(dt_step);
        }

        // 4) Solve fluid step with ALE
        params.options().setReal("timeStep", dt_step);
        solver.nextIteration();

        // 5) Evaluate traction at interface points (ALE-consistent) and write to preCICE
        gsField<> velField = solver.constructSolution(0);
        gsField<> prsField = solver.constructSolution(1);
        gsField<> meshDispField = solver.getMeshDisplacementField();

        // 5) Compute and write fluid traction
        gsMatrix<> traction;
        if (splineCoupling) {
            computeFluidTractionALE(velField, prsField, meshDispField, refFluidGeometry,
                                    cpParams, cpSideLocs, traction, density, viscosity, !pressureOnly);
            participant.writeData(FluidCPMesh, StressData, cpIDs, traction);
        } else {
            computeFluidTractionALE(velField, prsField, meshDispField, refFluidGeometry,
                                    parametricCoords, sideLocations, traction, density, viscosity, !pressureOnly);
            participant.writeData(FluidMesh, StressData, vertexIDs, traction);
        }

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

            // Update previous geometry for next step boundary velocity
            prevPatches = solver.getAssembler()->getPatches();

            // Plot based on simulation time interval, not step modulo
            if (plot && time + 1e-14 >= nextPlotTime)
            {
                gsField<> vel = velField;
                gsField<> pre = prsField;
                std::string suf = util::to_string(step);
                gsWriteParaview(vel, "fsi_flap_U_norem_" + suf, 400, false);
                gsWriteParaview(pre, "fsi_flap_P_norem_" + suf, 400, false);

                // Export deformed mesh geometry for visualization
                gsMultiPatch<> deformedMesh = solver.getAssembler()->getPatches();
                gsWriteParaview(deformedMesh, "fsi_mesh_deformed_" + suf, 400, false);
                // Advance next plot threshold
                nextPlotTime += plotEvery;

                for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
                {
                    pvU.addPart("fsi_flap_U_norem_" + suf + std::to_string(p) + ".vts", time, std::to_string(p));
                    pvP.addPart("fsi_flap_P_norem_" + suf + std::to_string(p) + ".vts", time, std::to_string(p));
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
    gsInfo << "[FSI-Fluid|NoRemesh] Completed.\n";
    return 0;
}
