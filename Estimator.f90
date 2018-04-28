MODULE Estimator
USE Constants, ONLY : DP, hbar
CONTAINS

SUBROUTINE KineticP()
USE Global

IMPLICIT NONE
INTEGER						:: i,j,k

KineticEnergy=0._DP
DO i=1,nat
	DO j=1,nbeads-1
		DO k=1,3
			KineticEnergy = KineticEnergy - 0.5_DP*Mass(i)*(wp**2)*((tau(i,j,k)-tau(i,j+1,k))**2)
		ENDDO
	ENDDO
	DO k=1,3
		KineticEnergy = KineticEnergy - 0.5_DP*Mass(i)*(wp**2)*((tau(i,nbeads,k)-tau(i,1,k))**2)
		KineticEnergy = KineticEnergy + ((Force_constraint(i,k)*nbeads)/(2._DP*beta))
	ENDDO
	
ENDDO


END SUBROUTINE KineticP


SUBROUTINE KineticV()
USE Global

IMPLICIT NONE

REAL(8), DIMENSION(nat,nbeads,3)		:: F
INTEGER						:: i,j,k

CALL FORCES(F)

KineticEnergy=0._DP
DO i=1,nat
	DO j=1,nbeads
		do k=1,3
			KineticEnergy = KineticEnergy - (1._DP/(2._DP*nbeads))*tau(i,j,k)*F(i,j,k)
		ENDDO
	ENDDO
ENDDO


END SUBROUTINE KineticV

SUBROUTINE KineticCV()
USE Global,	ONLY : tau, tau_c, nbeads, nat, KineticEnergy, beta, Force_constraint
IMPLICIT NONE

REAL(8), DIMENSION(nat,nbeads,3)		:: F
INTEGER						:: i,j,k

CALL Centroid()
CALL FORCES(F)

KineticEnergy=0._DP
DO i=1,nat
	DO j=1,nbeads
		do k=1,3
			KineticEnergy = KineticEnergy - (1._DP/(2._DP*nbeads))*(tau(i,j,k)-tau_c(i,k))*F(i,j,k)
		ENDDO
	ENDDO
ENDDO
KineticEnergy = KineticEnergy + ((SUM(Force_constraint(:,:)))/(2._DP*beta))


END SUBROUTINE KineticCV

SUBROUTINE Centroid()
USE Global, ONLY : tau_c, tau, nbeads, nat

IMPLICIT NONE

INTEGER		:: i,j

tau_c(:,:)=0._DP
DO i=1,nat
	DO j=1,3
		tau_c(i,j)=SUM(tau(i,:,j))/nbeads
	ENDDO
ENDDO

END SUBROUTINE Centroid

END MODULE Estimator
