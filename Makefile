OBJECTS = Global.o Init.o Read_namelist.o Staging.o Force.o Langevin.o Main.o 
F90COMP = gfortran
OPT = -fbounds-check

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
staging.mod: Staging.o Staging.f90
	$(F90COMP) -c $(OPT) Staging.f90
Staging.o: Staging.f90
	$(F90COMP) -c $(OPT) Staging.f90
%.o: %.f90
	$(F90COMP) -c $(OPT) $<

clean:
	rm *.mod
	rm $(OBJECTS)
