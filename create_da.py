from petsc4py import PETSc
da = PETSc.DA().create([16,16,16], stencil_width=1)

#da.view()
print da.comm.getRank(), da.getGhostCorners()
#print da.getGhostCorners()