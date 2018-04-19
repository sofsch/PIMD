MODULE Staging

CONTAINS

SUBROUTINE TRANSFORM()
USE Global, ONLY : nbeads,X,U
IMPLICIT NONE
INTEGER					:: i

U(:)=0._8
U(1)=X(1)

do i=2,nbeads-1
	U(i)=X(i) - ( ( (i-1._8)*X(i+1) + X(1) )/( DBLE(i) ) )
enddo
U(nbeads)=X(nbeads) - ( ( (nbeads-1._8)*X(1) + X(1) )/( DBLE(nbeads) ) ) 
END SUBROUTINE TRANSFORM


SUBROUTINE ITRANSFORM()
USE Global, ONLY : nbeads,U,X
IMPLICIT NONE
INTEGER					:: i,j

X(:)=0._8
X(1)=U(1)

do i=2,nbeads
	do j=i,nbeads
		X(i)=X(i) + ( ( (i-1._8)/(j-1._8) )*U(j) )
	enddo
	X(i)=X(i) + U(1)
enddo
END SUBROUTINE ITRANSFORM

END MODULE Staging
