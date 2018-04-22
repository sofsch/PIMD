!---------------------------------------------------------------------!
! Module Init_close :	read namelist input file		      !
!			allocate/deallocate                           !
!			Initialize				      !
!			INPUT FILE : PI.in			      !
!---------------------------------------------------------------------!
MODULE Init_close
USE Constants, ONLY : DP
CONTAINS

SUBROUTINE INITIALIZE()
USE constants, ONLY : AMU_AU, BOHR_RADIUS_SI, K_BOLTZMANN_AU, T_AU, AU_SEC
USE Global
USE Staging

IMPLICIT NONE
INTEGER							:: i,j
REAL(DP)						:: init_gamma_lang
NAMELIST /simulation/nat,ntyp,nbeads,nstep,dt,output
NAMELIST /dynamics/Temperature, init_gamma_lang
NAMELIST /species/species_label_mass
NAMELIST /ions/ions_name_pos
TYPE list_species
	CHARACTER(3)					:: Atom_label
	REAL(DP)					:: Mass_ref
END TYPE list_species	
TYPE list_ions
	CHARACTER(3)					:: Atom_name
	REAL(DP), DIMENSION(3)				:: Atom_position
	INTEGER, DIMENSION(3)				:: Force_c
END TYPE list_ions
TYPE (list_species), 	DIMENSION(:), ALLOCATABLE	:: species_label_mass
TYPE (list_ions), 	DIMENSION(:), ALLOCATABLE 	:: ions_name_pos


OPEN(10,FILE="PI.in")
READ(10,NML=simulation)
ALLOCATE(species_label_mass(ntyp))
ALLOCATE(ions_name_pos(nat))
READ(10,NML=dynamics)
READ(10,NML=species)
READ(10,NML=ions)
CLOSE(10)
dt=dt*2._DP

ALLOCATE(U(nat,nbeads,3))
ALLOCATE(P(nat,nbeads,3))
ALLOCATE(tau(nat,nbeads,3))
ALLOCATE(Mass(nat))
ALLOCATE(Mp(nat,nbeads),Ma(nat,nbeads))
ALLOCATE(gamma_lang(nbeads))
ALLOCATE(force_constraint(nat,3))


DO i=1,nat
	DO j=1,ntyp
		IF ( ions_name_pos(i)%Atom_name == species_label_mass(j)%Atom_label ) then
			Mass(i)=species_label_mass(j)%Mass_ref*AMU_AU
		ENDIF
	ENDDO
ENDDO

DO i=1,nat
	DO j=1,nbeads
		tau(i,j,:) = ions_name_pos(i)%Atom_position(:)/(BOHR_RADIUS_SI*1E+10_DP)
	ENDDO
ENDDO

DO i=1,nat
	force_constraint(i,:)= ions_name_pos(i)%Force_c(:)
ENDDO


CALL TRANSFORM()

P(:,:,:)=0._DP

Mp(:,1)=Mass(:)
Ma(:,1)=0._DP
output=trim(output)
DO i=1,nat
	DO j=2,nbeads
		Mp(i,j)=(j/(j-1._DP))*Mass(i)
		Ma(i,j)=Mp(i,j)
	ENDDO
ENDDO
!Temperature=Temperature*T_AU
Beta=1._DP/(K_BOLTZMANN_AU*Temperature)
wp=sqrt(1._DP*nbeads)/(Beta)!*hbar)
IF (init_gamma_lang .ne. 0) then
	gamma_lang(:)=init_gamma_lang*AU_SEC
ELSE
	gamma_lang(:)=wp
ENDIF

OPEN(15,FILE=TRIM(output)//".dist")
OPEN(16,FILE=TRIM(output)//".pos")
OPEN(17,FILE=TRIM(output)//".E")



END SUBROUTINE INITIALIZE

SUBROUTINE FINALIZE()
USE Global
IMPLICIT NONE

DEALLOCATE(U)
DEALLOCATE(P)
DEALLOCATE(tau)
DEALLOCATE(Mp,Ma,Mass)
DEALLOCATE(gamma_lang)
DEALLOCATE(force_constraint)
CLOSE(15)
CLOSE(16)
CLOSE(17)

END SUBROUTINE FINALIZE


END MODULE Init_close
