MODULE Dynamics

CONTAINS

SUBROUTINE BAOAB()
USE Langevin
USE Estimator
USE Staging
USE Global, ONLY : dt,nstep,nat,nbeads,KineticEnergy,pos_tot,tau

IMPLICIT NONE

INTEGER		:: i,j,k,l

l=0
DO i=1,nstep
	CALL B(dt/2._DP)
	CALL A(dt/2._DP)
	CALL O(dt)
	CALL A(dt/2._DP)
	CALL B(dt/2._DP)
	
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
	CALL KineticV()
	WRITE(17,*) i*dt, KineticEnergy
ENDDO

END SUBROUTINE BAOAB


SUBROUTINE VERLET()
USE Langevin
USE Estimator
USE Staging
USE Global, ONLY : dt,nstep,nat,nbeads,KineticEnergy,pos_tot,tau

IMPLICIT NONE

INTEGER		:: i,j,k,l

l=0
DO i=1,nstep
	CALL B(dt/2._DP)
	CALL A(dt)
	CALL B(dt/2._DP)
	
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
	CALL KineticV()
	WRITE(17,*) i*dt, KineticEnergy
ENDDO
END SUBROUTINE VERLET

END MODULE Dynamics




