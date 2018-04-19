PROGRAM Principal
USE Langevin
USE Global, ONLY : X,dt,nstep,nbeads
USE Staging
USE Init_close
IMPLICIT NONE
INTEGER			:: i,j,k,nproba
REAL(8), DIMENSION(:,:), ALLOCATABLE :: proba
REAL(8)					:: xmin,xmax,dx,norm

xmin=-3.e-10_8
xmax=3.e-10_8
dx=1.e-12_8
nproba=nint((xmax-xmin)/dx)+1



CALL Read_namelist()
Call INITIALIZE


ALLOCATE(proba(0:nproba,2))
proba(:,:)=0
do i=1,nstep
	CALL B(dt/2._8)
	CALL A(dt)
	!CALL O(dt)
	!CALL A(dt/2._8)
	CALL B(dt/2._8)
	
	CALL ITRANSFORM()
	do j=1,nbeads
		!write(16,*) X(j)
		do k=0,nproba
			if (X(j) .GT. (xmin + k*dx) .AND. X(j) .LT. (xmin + (k+1)*dx) ) then
				proba(k,2)=proba(k,2)+(1._8)
			endif
		enddo
	enddo
	if (MODULO((100.*(1.*i/nstep)),1.)==0.) then
		write(*,*) "test",100.*i/nstep
	endif
enddo
proba(:,:)=proba(:,:)/(nbeads)
norm=0._8
do i=0,nproba
	norm=norm + (proba(i,2)*dx)
enddo

proba(:,2)=(proba(:,2)*1.e-10_8)/norm
do i=0,nproba
	proba(i,1)= (xmin + i*dx)*1.e10_8
	write(15,*) proba(i,:)
enddo
 
CALL FINALIZE()
DEALLOCATE(proba)

END PROGRAM Principal
