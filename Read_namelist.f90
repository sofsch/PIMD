SUBROUTINE Read_namelist
USE Global
IMPLICIT NONE

TYPE list_ions
	CHARACTER(3)	:: Atom_name
	REAL(8)		:: Atom_position
END TYPE list_ions

TYPE (list_ions) ions_name_pos

NAMELIST /simulation/ nbeads, nstep, dt
NAMELIST /dynamics/ Temperature, init_gamma_lang
NAMELIST /ions/ Mass_ref,ions_name_pos

OPEN(10,FILE="PI.in")
READ(10,NML=simulation)
READ(10,NML=dynamics)
READ(10,NML=ions)
CLOSE(10)
init_pos = ions_name_pos%Atom_position

END SUBROUTINE Read_namelist
