!---------------------------------------------------------------------!
! Module staging : Define transformation from beads to normal modes   !
!---------------------------------------------------------------------!
MODULE Staging

CONTAINS

SUBROUTINE TRANSFORM()

USE Global, ONLY : nbeads,nat,tau,U

IMPLICIT NONE

INTEGER					:: i,j,k

U(:,:,:)=0._8
do i=1,nat
	do j=1,3
		U(i,1,j)=tau(i,1,j)
	enddo
enddo

do i=1,nat
	do j=2,nbeads-1
		do k=1,3
			U(i,j,k)=tau(i,j,k) - ( ( (i-1._8)*tau(i,j+1,k) + tau(i,1,k) )/( DBLE(j) ) )
		enddo
	enddo
enddo

do i=1,nat
	do j=1,3
		U(i,nbeads,j)=tau(i,nbeads,j) - ( ( (nbeads-1._8)*tau(i,1,j) + tau(i,1,j) )/( DBLE(nbeads) ) ) 
	enddo
enddo
END SUBROUTINE TRANSFORM


SUBROUTINE ITRANSFORM()

USE Global, ONLY : nbeads,nat,U,tau

IMPLICIT NONE

INTEGER					:: i,j,k,l

tau(:,:,:)=0._8
do i=1,nat
	do j=1,3
		tau(i,1,j)=U(i,1,j)
	enddo
enddo

do i=1,nat
	do j=2,nbeads
		do k=j,nbeads
			do l=1,3
				tau(i,j,l)=tau(i,j,l) + ( ( (j-1._8)/(k-1._8) )*U(i,k,l) )
			enddo
		enddo
		do l=1,3
			tau(i,j,l)=tau(i,j,l) + U(i,1,l)
		enddo
	enddo
enddo

END SUBROUTINE ITRANSFORM

END MODULE Staging
