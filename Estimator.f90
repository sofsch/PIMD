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
		do k=1,3
			KineticEnergy = KineticEnergy - ( (Mass(i)*nbeads)/( 2._DP*(beta**2)*(hbar**2) ))*((tau(i,j+1,k)-tau(i,j,k))**2)
		ENDDO
	ENDDO
	KineticEnergy = KineticEnergy - ( (Mass(i)*nbeads)/( 2._DP*(beta**2)*(hbar**2) ))*((tau(i,1,k)-tau(i,nbeads,k))**2)
	KineticEnergy = KineticEnergy + sum(Force_constraint(i,:))*nat*(nbeads/(2._DP*Beta))
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
	DO j=1,nbeads-1
		do k=1,3
			KineticEnergy = KineticEnergy - (1._DP/(2._DP*nbeads))*tau(i,j,k)*F(i,j,k)
		ENDDO
	ENDDO
ENDDO


END SUBROUTINE KineticV

END MODULE Estimator
