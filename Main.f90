PROGRAM Principal
USE Langevin
USE Global
USE Staging
IMPLICIT NONE
INTEGER			:: i,j

X(:)=1.e-12
U(:)=X(:)
P(:)=0.
Mp(1)=1.e-27
T=300.
dt=10.*4.8e-17
beta=1./(Kb*T)
wp=sqrt(1.*nbeads)/(beta*hbar)
gamma_lang(:)=1.e12
write(*,*) wp
write(*,*) beta
do i=2,nbeads
	Mp(i)=(i/(i-1))*1.e-27
enddo


OPEN(15,FILE="Beads.pos")

do i=1,nstep
	CALL B(dt/2)
	CALL A(dt/2)
	CALL O(dt)
	CALL A(dt/2)
	CALL B(dt/2)
	
	CALL ITRANSFORM(X,U)
	WRITE(15,*) i*dt, X(:)
	do j=1,nbeads
		write(16,*) X(j)
	enddo
enddo
 CLOSE(15)

END PROGRAM Principal
