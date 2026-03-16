# Required for finding pygismo
import os, sys

# Obtain pygismo
gismo_path = os.path.join(os.getcwd(), "../../../")
print("G+Smo path:", gismo_path, "(change if needed).")
sys.path.append(gismo_path + "build/lib")

import pygismo as gismo
import pygismo.iflow as iflow
import numpy as np
import scipy.sparse as sp

fd = gismo.io.gsFileData();

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

fd.add(mp)


basis_velocity = gismo.core.gsMultiBasis(mp)
basis_velocity.degreeElevate()
basis_pressure = gismo.core.gsMultiBasis(mp)

discreteBases = []
discreteBases.append(basis_velocity)
discreteBases.append(basis_pressure)

fd.add(discreteBases[0])
fd.add(discreteBases[1])

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

fd.add(bcs)

# zero force
force = gismo.core.gsFunctionExpr("0", "0", 2)

# Create Navier-Stokes PDE and solver params (omit logger)
viscosity = 0.01
NSpde = gismo.iflow.gsNavStokesPde(mp, bcs, force, viscosity)
print(NSpde)
params = gismo.iflow.gsFlowSolverParams(NSpde, discreteBases)

print("Created step geometry, BCs, bases, gsNavStokesPde and gsFlowSolverParams")

assembler = gismo.iflow.gsINSAssemblerSteady(params)
print("Created assembler")

assembler.initialize()
print("Initialized assembler")

# Get DOF counts
numdofs = assembler.numDofs()
# Prepare output
A = sp.csr_matrix((numdofs, numdofs))
b = np.zeros((numdofs, 1))
print("Prepared matrix and rhs")


# u = np.zeros((numdofs, 1))
# assembler.update(u,True)
# A = assembler.matrix()
# b = assembler.rhs()
# fd.add(A)
# fd.add(b)
# SOLVER STAGE
if not assembler.isInitialized():
    assembler.initialize()

# relNorm = np.norm(solNew-solOld)/np.norm(solNew)

minIterations = 0
maxIterations = 10
epsilon = 1e-4

u = np.zeros((numdofs, 1))
relNorm = 1.0
it = 0
while ((it < minIterations) or (relNorm > epsilon)) and (it < maxIterations):
    uold = u.copy()
    print(f"Iteration {it}")
    # Update assembler
    assembler.update(u,True)
    # Solve system
    A = assembler.matrix()
    b = assembler.rhs()
    # fd.add(A)
    # fd.add(b)

    # np.savetxt("matrix.txt", A.todense())

    u = sp.linalg.spsolve(A,b)
    # Compute relative norm
    relNorm = np.linalg.norm(u-uold)/np.linalg.norm(u)
    print(f"Iteration {it}: relNorm = {relNorm}")
    it += 1

# fd.dump()


# Construct solution
velocity = assembler.constructSolution(u,0)
pressure = assembler.constructSolution(u,1)

# gismo.io.gsWriteParaview(sol,"solution")

# Plot results using matplotlib
# Sample velocity field at a grid of points (per patch)
N = 20
npatches = mp.nPatches()
[xi, eta] = np.meshgrid(np.linspace(0, 1, N), np.linspace(0, 1, N))
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

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

plt.figure()
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
    # plt.contourf(X,Y,p_vals,levels=np.logspace(-9,-1,100),norm=LogNorm())
    plt.contourf(X,Y,p_vals,levels=np.linspace(-1e-1,1e-1,100))
    plt.streamplot(X,Y,u_vals,v_vals)
    # plt.quiver(X,Y,u_vals,v_vals,scale=2,scale_units='x')

plt.axis("equal")
plt.xlabel("x")
plt.ylabel("y")
plt.colorbar(label="Pressure")  
plt.show()

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

# # Create gsINSAssemblerSteady (just create, don't initialize yet)
# print("Creating assembler...")
# sys.stdout.flush()
# assembler = gismo.iflow.gsINSAssemblerSteady(params)
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
