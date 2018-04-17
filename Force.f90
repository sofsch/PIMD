SUBROUTINE FORCES(X,F)
USE Global, ONLY : nbeads,beta
IMPLICIT NONE
REAL(8), DIMENSION(nbeads), INTENT(in)	:: X
REAL(8), DIMENSION(nbeads), INTENT(out)	:: F
REAL(8)					:: E,b
INTEGER					:: i

E=0.1/beta
b=1.e-10
!write(*,*) E
do i=1,nbeads
	F(i)=-(E/(b**4))*(4.*(X(i)**3) - 4.*(b**2)*X(i)) !-Ma(i)*(freq**2)*X(i)
enddo
END SUBROUTINE FORCES
