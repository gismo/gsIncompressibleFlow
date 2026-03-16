import pygismo as gismo
#import pygismo.iflow as iflow
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.animation import FuncAnimation

# Build a 2D backward-facing step (deg=1, a=8, b=2, a_in=1, h=1)
mp = gismo.core.gsMultiPatch()
patch0 = gismo.nurbs.gsNurbsCreator.BSplineSquare()
patch0.coefs()[:,0] *= 8.0

patch1 = gismo.nurbs.gsNurbsCreator.BSplineSquare()
patch1.coefs()[:,0] *= 8.0
patch1.coefs()[:,1] += 1.0

patch2 = gismo.nurbs.gsNurbsCreator.BSplineSquare()
patch2.coefs()[:,0] -= 1.0
patch2.coefs()[:,1] += 1.0

mp.addPatch(patch0)
mp.addPatch(patch1)
mp.addPatch(patch2)
mp.computeTopology(1e-6, True, True)

# Apply knot refinement like BSplineStep2D does
# aNumElem = floor(2*a/b) = floor(2*8/2) = 8
# ainNumElem = floor(2*a_in/b) = floor(2*1/2) = 1
import math

aNumElem = 8
ainNumElem = 1
aStep = 1.0 / aNumElem
ainStep = 1.0 / ainNumElem

# Insert knots in patches 0 and 1 (x-direction, direction=0)
for p in range(2):
    for i in range(1, aNumElem):
        mp.patch(p).insertKnot(i * aStep, 0)

# Insert knots in patch 2 (x-direction, direction=0)
for i in range(1, ainNumElem):
    mp.patch(2).insertKnot(i * ainStep, 0)

basis_velocity = gismo.core.gsMultiBasis(mp)
basis_velocity.degreeElevate()
basis_velocity.uniformRefine()
basis_velocity.uniformRefine()
basis_velocity.uniformRefine()
basis_pressure = gismo.core.gsMultiBasis(mp)
basis_pressure.uniformRefine()
basis_pressure.uniformRefine()
basis_pressure.uniformRefine()

discreteBases = []
discreteBases.append(basis_velocity)
discreteBases.append(basis_pressure)


# Boundary conditions following defineBCs_step (2D, non-periodic)
bcs = gismo.pde.gsBoundaryConditions()
Uin = gismo.core.gsFunctionExpr("(-4*(y-1.5)^2 + 1)", "0", 2)
Uwall = gismo.core.gsFunctionExpr("0", "0", 2)
Pout = gismo.core.gsFunctionExpr("0", 2)

# inflow: patch 2, west
bcs.addCondition(2, gismo.core.side.west, gismo.pde.bctype.dirichlet, Uin, 0, False, -1)
# walls
bcs.addCondition(0, gismo.core.side.west, gismo.pde.bctype.dirichlet, Uwall, 0, False, -1)
bcs.addCondition(0, gismo.core.side.south, gismo.pde.bctype.dirichlet, Uwall, 0, False, -1)
bcs.addCondition(1, gismo.core.side.north, gismo.pde.bctype.dirichlet, Uwall, 0, False, -1)
bcs.addCondition(2, gismo.core.side.south, gismo.pde.bctype.dirichlet, Uwall, 0, False, -1)
bcs.addCondition(2, gismo.core.side.north, gismo.pde.bctype.dirichlet, Uwall, 0, False, -1)
# outflow (pressure): patches 0 and 1, east, unknown = 1
bcs.addCondition(0, gismo.core.side.east, gismo.pde.bctype.dirichlet, Pout, 1, False, -1)
bcs.addCondition(1, gismo.core.side.east, gismo.pde.bctype.dirichlet, Pout, 1, False, -1)

# attach geometry map
bcs.setGeoMap(mp)

# zero force
force = gismo.core.gsFunctionExpr("0", "0", 2)

# Create Navier-Stokes PDE and solver params (omit logger)
time = 0.0
dt = 0.5
maxTime = 10.0

viscosity = 0.01
NSpde = gismo.iflow.gsNavStokesPde(mp, bcs, force, viscosity)
params = gismo.iflow.gsFlowSolverParams(NSpde, discreteBases)
params.options().setReal("timeStep", dt)
print(params.options())

print("Created step geometry, BCs, bases, gsNavStokesPde and gsFlowSolverParams")

assembler = gismo.iflow.gsINSAssemblerUnsteady(params)
print("Created assembler")

assembler.initialize()
print("Initialized assembler")
# SOLVER STAGE
if not assembler.isInitialized():
    assembler.initialize()

numdofs = assembler.numDofs()


minIterations = 0
maxIterations = 10
epsilon = 1e-4

u = np.zeros((numdofs, 1))
frames = []
time += dt
while time < maxTime:
    print(f"Time step {time}")
    relNorm = 1
    it = 0
    assembler.update(u,True)  # Update assembler for new time step (with nonlinearity)
    A = assembler.matrix()
    b = assembler.rhs()
    utmp = sp.linalg.spsolve(A,b)
    while ((it < minIterations) or (relNorm > epsilon)) and (it < maxIterations):
        uold = utmp.copy()
        # Update assembler
        assembler.update(utmp,False)
        # Solve system
        A = assembler.matrix()
        b = assembler.rhs()
        utmp = sp.linalg.spsolve(A,b)
        print(
                # "Matrix norm:", np.linalg.norm(A)
            #   , 
              "RHS norm:", np.linalg.norm(b)
              , "Solution norm:", np.linalg.norm(utmp))
        # Compute relative norm
        relNorm = np.linalg.norm(utmp-uold)/np.linalg.norm(utmp)
        # Print info
        print(f"Iteration {it}, relative norm = {relNorm}")

        it += 1
    
    u = utmp.copy()
    time += dt

    # Construct solution
    velocity = assembler.constructSolution(u,0)
    pressure = assembler.constructSolution(u,1)

    # Sample velocity field at a grid of points (per patch)
    N = 20
    [xi, eta] = np.meshgrid(np.linspace(0, 1, N), np.linspace(0, 1, N))

    frame_data = []
    for p in range(mp.nPatches()):
        X = np.zeros((N,N))
        Y = np.zeros((N,N))
        u_vals = np.zeros((N,N))
        v_vals = np.zeros((N,N))
        p_vals = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                point = [xi[i, j], eta[i, j]]
                ppoint = mp.patch(p).eval(point)
                X[i, j] = ppoint[0,0]
                Y[i, j] = ppoint[1,0]
                uv = velocity.value(point,p)
                u_vals[i, j] = uv[0,0]
                v_vals[i, j] = uv[1,0]
                P = pressure.value(point,p)
                p_vals[i, j] = P[0,0]   
        frame_data.append({'X':X, 'Y':Y, 'u':u_vals, 'v':v_vals, 'p':p_vals})
    frames.append(frame_data)

# Animation
fig, ax = plt.subplots()
def animate(frame_num):
    ax.clear()
    ax.set_title(f"Velocity field and pressure at time {(frame_num+1)*dt:.2f}")
    for p_data in frames[frame_num]:
        ax.contourf(p_data['X'], p_data['Y'], p_data['p'], levels=np.linspace(-1e-1,1e-1,100))
        ax.streamplot(p_data['X'], p_data['Y'], p_data['u'], p_data['v'])
    ax.axis("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

ani = FuncAnimation(fig, animate, frames=len(frames), interval=500, repeat=True)
# Add colorbar
cs = ax.contourf(frames[-1][0]['X'], frames[-1][0]['Y'], frames[-1][0]['p'], levels=np.linspace(-1e-1,1e-1,100))
plt.colorbar(cs, label="Pressure")
plt.show()

# X = np.zeros((N*npatches,N*npatches))
# Y = np.zeros((N*npatches,N*npatches))
# u_vals = np.zeros((N*npatches,N*npatches))
# v_vals = np.zeros((N*npatches,N*npatches))
# p_vals = np.zeros((N*npatches,N*npatches))
# for p in range(mp.nPatches()):
#     for i in range(N):
#         for j in range(N):
#             point = [xi[i, j], eta[i, j]]
#             ppoint = mp.patch(p).eval(point)
#             X[p*N + i, p*N + j] = ppoint[0,0]
#             Y[p*N + i, p*N + j] = ppoint[1,0]
#             uv = velocity.value(point,p)
#             u_vals[p*N + i, p*N + j] = uv[0,0]
#             v_vals[p*N + i, p*N + j] = uv[1,0]
#             P = pressure.value(point,p)
#             p_vals[p*N + i, p*N + j] = P[0,0]

# # Plot in the physical domain
# import matplotlib.pyplot as plt
# plt.figure()
# plt.quiver(X, Y, u_vals, v_vals)
# plt.title("Velocity field")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.axis("equal")
# plt.contourf(X, Y, p_vals, levels=20)
# plt.colorbar(label="Pressure")
# plt.show()


    # # Plot velocity field using quiver

    # plt.figure()
    # plt.quiver(X, Y, u_vals, v_vals)
    # plt.title(f"Velocity field on patch {p}")
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.axis("equal")
    # plt.show()

# # Assemble Stokes system
# assembler.fillStokesSystem(A, b)
# print("Assembled Stokes system")

# # Use gsINSSolverSteady instead of gsINSAssembler directly
# # gsINSSolver manages object lifetimes internally, avoiding the GC issue
# solver = gismo.iflow.gsINSSolverSteady(params)

# # Initialize and solve Stokes to get the matrix/rhs
# solver.initialize()
# solver.solveStokes()

# # Get the assembler (which has the assembled system)
# assembler = solver.getAssembler()

# # Get the assembled matrix and rhs
# matrix = assembler.matrix()
# rhs = assembler.rhs()

# print(f"Assembled matrix size: {matrix.rows()} x {matrix.cols()}")
# print(f"Assembled rhs size: {rhs.rows()} x {rhs.cols()}")
# print("Matrix and rhs assembled successfully!")

# sys.stdout.flush()

# # Create gsINSAssemblerUnsteady (just create, don't initialize yet)
# print("Creating assembler...")
# sys.stdout.flush()
# assembler = gismo.iflow.gsINSAssemblerUnsteady(params)
# print("Assembler created successfully!")
# sys.stdout.flush()

# # Try to initialize
# print("Initializing assembler...")
# assembler.initialize()
# print("Assembler initialized successfully!")

# # Get initial solution (zeros)
# num_dofs = assembler.getUdofs() + assembler.getPdofs()
# print(f"Number of DOFs: {num_dofs}")

# # Create zero solution vector as numpy array
# import numpy as np
