!---------------------------------------------------------------------!
! Main program							      !
!---------------------------------------------------------------------!
PROGRAM Principal

USE Langevin
USE Global!, ONLY : tau,dt,nstep,nbeads
USE Staging
USE Init_close

IMPLICIT NONE

INTEGER			:: i,j,k,nproba
REAL(8), DIMENSION(:,:), ALLOCATABLE :: proba
REAL(8)					:: xmin,xmax,dx,norm

xmin=0._8
xmax=2.e-10_8
dx=1.e-12_8
nproba=nint((xmax-xmin)/dx)+1

CALL INITIALIZE()

ALLOCATE(proba(0:nproba,2))
proba(:,:)=0
DO i=1,nstep
	CALL B(dt/2._8)
	CALL A(dt/2._8)
	CALL O(dt)
	CALL A(dt/2._8)
	CALL B(dt/2._8)
	
	CALL ITRANSFORM()
	write(16,*) i*dt,tau(1,:,1)
	
	DO j=1,nbeads
		DO k=0,nproba
			IF (abs(tau(1,j,1)) .GT. (xmin + k*dx) .AND. abs(tau(1,j,1)) .LT. (xmin + (k+1)*dx) ) then
				proba(k,2)=proba(k,2)+(1._8)
			ENDIF
		ENDDO
	ENDDO
	IF (MODULO((100.*(1.*i/nstep)),1.)==0.) then
		write(*,*) "test",100.*i/nstep
	ENDIF
ENDDO
proba(:,:)=proba(:,:)/(nbeads)
norm=0._8
DO i=0,nproba
	norm=norm + (proba(i,2)*dx)
ENDDO

proba(:,2)=(proba(:,2))/norm
DO i=0,nproba
	proba(i,1)= (xmin + i*dx)*1.e10_8
	write(15,*) proba(i,:)
ENDDO
 
CALL FINALIZE()
DEALLOCATE(proba)

END PROGRAM Principal
