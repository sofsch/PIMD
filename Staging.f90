MODULE Staging

CONTAINS

SUBROUTINE TRANSFORM(X,U)
USE Global, ONLY : nbeads
IMPLICIT NONE
REAL(8), DIMENSION(nbeads), INTENT(in)	:: X
REAL(8), DIMENSION(nbeads), INTENT(out)	:: U
INTEGER					:: i

U=0
U(1)=X(1)

do i=2,nbeads
	U(i)=X(i) - ( ( (i-1)*X(i+1) + X(1) )/( i ) )
enddo
END SUBROUTINE TRANSFORM


SUBROUTINE ITRANSFORM(X,U)
USE Global, ONLY : nbeads
IMPLICIT NONE
REAL(8), DIMENSION(nbeads), INTENT(out)	:: X
REAL(8), DIMENSION(nbeads), INTENT(in)	:: U
INTEGER					:: i,j

X=0
X(1)=U(1)

do i=2,nbeads
	do j=i,nbeads
		X(i)=X(i) + ( ( (i-1)/(j-1) )*U(j) )
	enddo
	X(i)=X(i) + U(1)
enddo
END SUBROUTINE ITRANSFORM

END MODULE Staging
