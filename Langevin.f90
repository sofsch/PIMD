!---------------------------------------------------------------------!
! Module Langevin : Langevin dynamics (BAOAB integrator)	      !
!---------------------------------------------------------------------!
MODULE Langevin
USE Constants, ONLY : DP
CONTAINS

SUBROUTINE B(deltaT)

USE Staging
USE Global, ONLY : nbeads,nat,P

IMPLICIT NONE

REAL(DP),	DIMENSION(nat,nbeads,3)		:: F,F_staging
REAL(DP),	INTENT(in)			:: deltaT
INTEGER						:: i,j,k

CALL ITRANSFORM()
CALL FORCES(F)

F_staging(:,:,:)=0._DP

DO i=1,nat
	DO j=1,nbeads
		DO k=1,3
			F_staging(i,1,k)=F_staging(i,1,k)+(F(i,j,k)/DBLE(nbeads))
		ENDDO
	ENDDO
ENDDO

DO i=1,nat
	DO j=2,nbeads
		DO k=1,3
			F_staging(i,j,k)= (F(i,j,k)/DBLE(nbeads)) + ((( (j-2._DP)/(j-1._DP) )*F_staging(i,j-1,k)))
		ENDDO
	ENDDO
ENDDO


DO i=1,nat
	DO j=1,nbeads
		DO k=1,3
			P(i,j,k)=P(i,j,k) + (F_staging(i,j,k)*deltaT)
		ENDDO
	ENDDO
ENDDO
END SUBROUTINE B


SUBROUTINE A(deltaT)

USE Global, ONLY : nbeads,nat,U,P,Mp,wp,Mass

IMPLICIT NONE

REAL(DP),	DIMENSION(nat,nbeads,3)		:: Utmp, Ptmp
REAL(DP),	INTENT(in)			:: deltaT
INTEGER						:: i,j,k


Utmp(:,:,:)=0._DP
Ptmp(:,:,:)=0._DP

DO i=1,nat
	DO j=1,3
		Utmp(i,1,j)=U(i,1,j) + ( ( P(i,1,j)/Mp(i,j) )*deltaT )
		Ptmp(i,1,j)=P(i,1,j)
	ENDDO
ENDDO

DO i=1,nat
	DO j=2,nbeads
		DO k=1,3
			Utmp(i,j,k) = cos(wp*deltaT)*U(i,j,k) + ( sin(wp*deltaT)/(wp*Mp(i,j)) )*P(i,j,k)
			Ptmp(i,j,k) = -wp*sin(wp*deltaT)*Mp(i,j)*U(i,j,k) + cos(wp*deltaT)*P(i,j,k)
		ENDDO
	ENDDO
ENDDO

U(:,:,:)=Utmp(:,:,:)
P(:,:,:)=Ptmp(:,:,:)
END SUBROUTINE A


SUBROUTINE O(deltaT)

USE Global, ONLY : nbeads,nat,P,gamma_lang,Beta,Mp,force_constraint

IMPLICIT NONE

REAL(DP),	DIMENSION(nat,nbeads,3)		:: R
REAL(DP),	INTENT(in)			:: deltaT
INTEGER						:: i,j,k

DO i=1,nat
	DO j=1,nbeads
		DO k=1,3
			R(i,j,k)=GAUSS()*force_constraint(i,k)
		ENDDO
	ENDDO
ENDDO

DO i=1,nat
	DO j=1,nbeads
		DO k=1,3
			P(i,j,k)=exp(-gamma_lang(j)*deltaT)*P(i,j,k) + sqrt( 1._DP- exp(-2._DP*gamma_lang(j)*deltaT) )*sqrt(Mp(i,j)/Beta)*R(i,j,k)
		ENDDO
	ENDDO
ENDDO
END SUBROUTINE O

REAL(DP) FUNCTION GAUSS()

IMPLICIT NONE

REAL(DP),	PARAMETER	:: a1 = 3.949846138_DP, a3 = 0.252408784_DP, a5 = 0.076542912_DP, a7 = 0.008355968_DP, a9 = 0.029899776_DP
REAL(DP)			:: r, r2
REAL(DP),	DIMENSION(12)	:: s

CALL RANDOM_NUMBER(s) ; r = (sum(s)-6.0_DP)/4.0_DP ; r2 = r*r
GAUSS =  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 )*r
END FUNCTION GAUSS

END MODULE Langevin
