SUBROUTINE FORCES(X,F)
USE Global, ONLY : nbeads,Kb
IMPLICIT NONE
REAL(8), DIMENSION(nbeads), INTENT(in)	:: X
REAL(8), DIMENSION(nbeads), INTENT(out)	:: F
REAL(8)					:: h,d
INTEGER					:: i

h=1000.*Kb
d=0.6e-10
do i=1,nbeads
	F(i)=-8.*h*(8.*((X(i)**3)/(d**4)) - 2.*(X(i)/(d**2)) ) !-Ma(i)*(freq**2)*X(i)
enddo
END SUBROUTINE FORCES
