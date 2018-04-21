!---------------------------------------------------------------------!
! Module Init_close :	read namelist input file		      !
!			allocate/deallocate                           !
!			Initialize				      !
!			INPUT FILE : PI.in			      !
!---------------------------------------------------------------------!
MODULE Init_close

CONTAINS

SUBROUTINE INITIALIZE()

USE Global
USE Staging

IMPLICIT NONE
INTEGER					:: i,j
NAMELIST /simulation/nat,ntyp,nbeads,nstep,dt,output
NAMELIST /dynamics/Temperature, init_gamma_lang
NAMELIST /species/species_label_mass
NAMELIST /ions/ions_name_pos
TYPE list_species
	CHARACTER(3)			:: Atom_label
	REAL(8)				:: Mass_ref
END TYPE list_species
TYPE list_ions
	CHARACTER(3)			:: Atom_name
	REAL(8), DIMENSION(3)		:: Atom_position
	INTEGER, DIMENSION(3)		:: Force_c
END TYPE list_ions
TYPE (list_species), DIMENSION(:), ALLOCATABLE :: species_label_mass
TYPE (list_ions), DIMENSION(:), ALLOCATABLE :: ions_name_pos

OPEN(10,FILE="PI.in")
READ(10,NML=simulation)
ALLOCATE(species_label_mass(ntyp))
ALLOCATE(ions_name_pos(nat))
READ(10,NML=dynamics)
READ(10,NML=species)
READ(10,NML=ions)
CLOSE(10)


ALLOCATE(U(nat,nbeads,3))
ALLOCATE(P(nat,nbeads,3))
ALLOCATE(tau(nat,nbeads,3))
ALLOCATE(Mass(nat))
ALLOCATE(Mp(nat,nbeads),Ma(nat,nbeads))
ALLOCATE(gamma_lang(nbeads))
ALLOCATE(force_constraint(nat,3))


do i=1,nat
	do j=1,ntyp
		if ( ions_name_pos(i)%Atom_name == species_label_mass(j)%Atom_label ) then
			Mass(i)=species_label_mass(j)%Mass_ref
		endif
	enddo
enddo

do i=1,nat
	do j=1,nbeads
		tau(i,j,:) = ions_name_pos(i)%Atom_position(:)
	enddo
enddo

do i=1,nat
	force_constraint(i,:)= ions_name_pos(i)%Force_c(:)
enddo


CALL TRANSFORM()

P(:,:,:)=0._8

Mp(:,1)=Mass(:)
Ma(:,1)=0._8
output=trim(output)
do i=1,nat
	do j=2,nbeads
		Mp(i,j)=(j/(j-1._8))*Mass(i)
		Ma(i,j)=Mp(i,j)
	enddo
enddo

beta=1._8/(Kb*Temperature)
wp=sqrt(1._8*nbeads)/(beta*hbar)
if (init_gamma_lang .ne. 0) then
	gamma_lang(:)=init_gamma_lang
else
	gamma_lang(:)=wp
endif

OPEN(15,FILE=TRIM(output)//".dist")
OPEN(16,FILE=TRIM(output)//".pos")



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

END SUBROUTINE FINALIZE


END MODULE Init_close
