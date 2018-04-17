SUBROUTINE FORCES(X,F)
USE Global, ONLY : nbeads,wp,Mp
IMPLICIT NONE
REAL(8), DIMENSION(nbeads), INTENT(in)	:: X
REAL(8), DIMENSION(nbeads), INTENT(out)	:: F
REAL(8)					:: k
INTEGER					:: i


do i=1,nbeads
	F(i)=-Mp(i)*(wp**2)*X(i)
enddo
END SUBROUTINE FORCES
