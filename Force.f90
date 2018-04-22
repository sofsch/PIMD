!---------------------------------------------------------------------!
! Return the inter-atomic force component F                           !
!---------------------------------------------------------------------!
SUBROUTINE FORCES(F)
USE Constants, ONLY : DP,K_BOLTZMANN_AU, BOHR_RADIUS_SI, AU_SEC, hbar
USE Global, ONLY 	: nbeads,nat,tau,force_constraint,Ma

IMPLICIT NONE
REAL(DP), DIMENSION(nat,nbeads,3), INTENT(out)		:: F
REAL(DP)						:: h,d,freq,s,a
INTEGER							:: i,j,k

!Initialize
F(:,:,:)=0._DP
h=1000._DP*K_BOLTZMANN_AU
d=0.6*1e-10_DP/BOHR_RADIUS_SI
freq=100.e12_DP*AU_SEC
s=0.5e-10_DP/BOHR_RADIUS_SI
a=0.02_DP

DO i=1,nat
	DO j=1,nbeads
		DO k=1,3
			F(i,j,k)=-Ma(i,j)*(freq**2)*tau(i,j,k) + &
			&(a*tau(i,j,k)/(s**2))*exp( - (tau(i,j,1)**2 + tau(i,j,2)**2)/(2._DP*(s**2)))
			F(i,j,k)=F(i,j,k)*force_constraint(i,k)
		ENDDO
		
	ENDDO
ENDDO

END SUBROUTINE FORCES

!2D Harmonic
!F(i,j,k)=-Ma(i,j)*(freq**2)*tau(i,j,k)
!F(i,j,k)=F(i,j,k)*force_constraint(i,k)


!Double well from Ceriotti
!-16._DP*h*((4._DP*(tau(i,j,k)**3_DP)/(d**4_DP)) - (tau(i,j,k)/(d**2_DP)) ) 
