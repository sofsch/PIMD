!---------------------------------------------------------------------!
! Return the inter-atomic force component F                           !
!---------------------------------------------------------------------!
SUBROUTINE FORCES(F)

USE Global, ONLY : nbeads,nat,Kb,tau,force_constraint,Ma

IMPLICIT NONE
REAL(8), DIMENSION(nat,nbeads,3), INTENT(out)	:: F
REAL(8)						:: h,d,freq
INTEGER						:: i,j,k

!Initialize
F(:,:,:)=0._8
h=1000._8*Kb
d=0.6*1e-10_8
freq=100.e12_8

DO i=1,nat
	DO j=1,nbeads
		DO k=1,3
			F(i,j,k)=-Ma(i,j)*(freq**2)*tau(i,j,k)
			F(i,j,k)=F(i,j,k)*force_constraint(i,k)
		ENDDO
		
	ENDDO
ENDDO

END SUBROUTINE FORCES

!Double well from Ceriotti
!-16._8*h*((4._8*(tau(i,j,k)**3_8)/(d**4_8)) - (tau(i,j,k)/(d**2_8)) ) 
