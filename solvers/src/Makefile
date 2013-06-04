all: FDTDSolvers2D.so PMLBuild.so

FDTDSolvers2D.so: FDTDSolvers2D.f90
	f2py  -m $(basename $(notdir $@)) -c $^

PMLBuild.so: FDTDPMLBuild.f90
	f2py  -m $(basename $(notdir $@)) -c $^
