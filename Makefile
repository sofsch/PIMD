OBJECTS = Constants.o Global.o Staging.o Init.o Force.o Langevin.o Estimator.o Dynamics.o Distributions.o Main.o 
#F90COMP = mpif90
F90COMP = gfortran
OPT = -O3
#OPT = -fopenmp -O3
#OPT = -fbounds-check -Wall -Wno-tabs

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
estimator.mod: Estimator.o Estimator.f90
	$(F90COMP) -c $(OPT) Estimator.f90
Estimator.o: Estimator.f90
	$(F90COMP) -c $(OPT) Estimator.f90
distributions.mod: Distributions.o Distributions.f90
	$(F90COMP) -c $(OPT) Distributions.f90
Distributions.o: Distributions.f90
	$(F90COMP) -c $(OPT) Distributions.f90
dynamics.mod: Dynamics.o Dynamics.f90
	$(F90COMP) -c $(OPT) Dynamics.f90
Dynamics.o: Distributions.f90
	$(F90COMP) -c $(OPT) Dynamics.f90	
constants.mod: Constants.o Constants.f90
	$(F90COMP) -c $(OPT) Constants.f90
Constants.o: Constants.f90
	$(F90COMP) -c $(OPT) Constants.f90
Init.o: Init.f90
	$(F90COMP) -c $(OPT) Init.f90
%.o: %.f90
	$(F90COMP) -c $(OPT) $<

clean:
	rm *.mod
	rm $(OBJECTS)
	rm PI.x
