SUBROUTINE FORCES(X,F)
USE Global, ONLY : nbeads,Ma,wp
IMPLICIT NONE
REAL(8), DIMENSION(nbeads), INTENT(in)	:: X
REAL(8), DIMENSION(nbeads), INTENT(out)	:: F
REAL(8)					:: freq
INTEGER					:: i

freq=100.e12
do i=1,nbeads
	F(i)=-Ma(i)*(freq**2)*X(i)
enddo
END SUBROUTINE FORCES
