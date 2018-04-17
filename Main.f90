PROGRAM Principal
USE Langevin
USE Global
USE Staging
USE Init_close
IMPLICIT NONE
INTEGER			:: i,j,nproba
REAL(8), DIMENSION(:,:), ALLOCATABLE :: proba
REAL(8)					:: xmin,xmax,dx,norm

xmin=0.
xmax=3.e-9
dx=1e-12
nproba=nint((xmax-xmin)/dx)+1



CALL Read_namelist()
Call INITIALIZE


ALLOCATE(proba(0:nproba,2))

proba(:,:)=0

do i=1,nstep
	CALL B(dt/2)
	CALL A(dt/2)
	CALL O(dt)
	CALL A(dt/2)
	CALL B(dt/2)
	
	CALL ITRANSFORM(X,U)
	!WRITE(15,*) i*dt, X(:), P(:)
	do j=1,nbeads
		!write(16,*) X(j)
		!write(17,*) sum(X(:))/nbeads
		
		proba(nint(abs(X(j)/dx)),2)=proba(nint(abs(X(j)/dx)),2)+1.
	enddo
	write(*,*) i*dt
enddo

norm=0
do i=0,nproba
	norm=norm + (proba(i,2)*dx)
enddo

proba(:,2)=proba(:,2)/norm
do i=0,nproba
	proba(i,1)= xmin + i*dx
	write(18,*) proba(i,:)
enddo
 
CALL FINALIZE()
DEALLOCATE(proba)

END PROGRAM Principal
