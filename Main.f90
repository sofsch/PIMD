!---------------------------------------------------------------------!
! Main program							      !
!---------------------------------------------------------------------!
PROGRAM Principal

USE Langevin
USE Global!, ONLY : tau,dt,nstep,nbeads
USE Staging
USE Init_close
USE Distributions

IMPLICIT NONE

INTEGER						:: i,j,k,l
REAL(8), DIMENSION(:,:), ALLOCATABLE 		:: proba
REAL(8), DIMENSION(:,:,:), ALLOCATABLE		:: pos_tot


CALL INITIALIZE()
ALLOCATE(pos_tot(nstep*nbeads,nat,2))
l=0
DO i=1,nstep
	CALL B(dt/2._8)
	CALL A(dt/2._8)
	CALL O(dt)
	CALL A(dt/2._8)
	CALL B(dt/2._8)
	
	CALL ITRANSFORM()
	write(16,*) i*dt,tau(1,:,1)
	
	DO j=1,nat
		DO k=1,nbeads
			l=l+1
			pos_tot(l,j,1)=tau(j,k,1)
			pos_tot(l,j,2)=tau(j,k,2)
		ENDDO
	ENDDO

	IF (MODULO((100.*(1.*i/nstep)),1.)==0.) then
		write(*,*) "test",100.*i/nstep
	ENDIF
ENDDO



CALL distribution2d(pos_tot(:,1,1),pos_tot(:,1,2),proba,1.e-12_8,1.e-12_8,-2.e-10_8,2.e-10_8,-2.e-10_8,2.e-10_8,.TRUE.)


DO i=1,size(proba,1)-1
	WRITE(15,*) proba(i,:)
	if (proba(i+1,1) .NE. proba(i,1)) THEN
		WRITE(15,*) ""
	endif
ENDDO

 
CALL FINALIZE()
DEALLOCATE(pos_tot)
DEALLOCATE(proba)

END PROGRAM Principal
