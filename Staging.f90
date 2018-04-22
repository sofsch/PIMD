!---------------------------------------------------------------------!
! Module staging : Define transformation from beads to normal modes   !
!---------------------------------------------------------------------!
MODULE Staging
USE Constants, ONLY : DP
CONTAINS

SUBROUTINE TRANSFORM()

USE Global, ONLY : nbeads,nat,tau,U

IMPLICIT NONE

INTEGER					:: i,j,k

U(:,:,:)=0._DP
DO i=1,nat
	DO j=1,3
		U(i,1,j)=tau(i,1,j)
	ENDDO
ENDDO

DO i=1,nat
	DO j=2,nbeads-1
		DO k=1,3
			U(i,j,k)=tau(i,j,k) - ( ( (i-1._DP)*tau(i,j+1,k) + tau(i,1,k) )/( DBLE(j) ) )
		ENDDO
	ENDDO
ENDDO

DO i=1,nat
	DO j=1,3
		U(i,nbeads,j)=tau(i,nbeads,j) - ( ( (nbeads-1._DP)*tau(i,1,j) + tau(i,1,j) )/( DBLE(nbeads) ) ) 
	ENDDO
ENDDO
END SUBROUTINE TRANSFORM


SUBROUTINE ITRANSFORM()

USE Global, ONLY : nbeads,nat,U,tau

IMPLICIT NONE

INTEGER					:: i,j,k,l

tau(:,:,:)=0._DP
DO i=1,nat
	DO j=1,3
		tau(i,1,j)=U(i,1,j)
	ENDDO
ENDDO

DO i=1,nat
	DO j=2,nbeads
		DO k=j,nbeads
			DO l=1,3
				tau(i,j,l)=tau(i,j,l) + ( ( (j-1._DP)/(k-1._DP) )*U(i,k,l) )
			ENDDO
		ENDDO
		DO l=1,3
			tau(i,j,l)=tau(i,j,l) + U(i,1,l)
		ENDDO
	ENDDO
ENDDO

END SUBROUTINE ITRANSFORM

END MODULE Staging
