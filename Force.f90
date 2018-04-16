MODULE Force
IMPLICIT NONE
CONTAINS

SUBROUTINE FORCES(X,F)
USE Global, ONLY : nbeads
IMPLICIT NONE
REAL(8), DIMENSION(nbeads), INTENT(in)	:: X
REAL(8), DIMENSION(nbeads), INTENT(out)	:: F
REAL(8)					:: k
INTEGER					:: i

k=0.5e-2

do i=1,nbeads
	F(i)=-k*X(i)
enddo
END SUBROUTINE FORCES

END MODULE Force

