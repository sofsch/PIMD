OBJECTS = Global.o Staging.o Init.o Force.o Langevin.o Main.o 
#F90COMP = mpif90
F90COMP = gfortran
#OPT = -O3
#OPT = -fopenmp -O3
OPT = -fbounds-check -Wall -Wno-tabs

PI.x: $(OBJECTS)
	$(F90COMP) $(OPT) $(OBJECTS) -o PI.x
global.mod: Global.o Global.f90
	$(F90COMP) -c $(OPT) Global.f90
Global.o: Global.f90
	$(F90COMP) -c $(OPT) Global.f90
langevin.mod: Langevin.o Langevin.f90
	$(F90COMP) -c $(OPT) Langevin.f90
Langevin.o: Langevin.f90
	$(F90COMP) -c $(OPT) Langevin.f90
init.mod: Init.o Init.f90
	$(F90COMP) -c $(OPT) Init.f90
staging.mod: Staging.o Staging.f90
	$(F90COMP) -c $(OPT) Staging.f90
Staging.o: Staging.f90
	$(F90COMP) -c $(OPT) Staging.f90
Init.o: Init.f90
	$(F90COMP) -c $(OPT) Init.f90
%.o: %.f90
	$(F90COMP) -c $(OPT) $<

clean:
	rm *.mod
	rm $(OBJECTS)
	rm PI.x
