MODULE Langevin

CONTAINS

SUBROUTINE B(deltaT)
USE Force
USE Global, ONLY : nbeads,U,P
IMPLICIT NONE
REAL(8), DIMENSION(nbeads)			:: X,F,F_staging
REAL(8), INTENT(in)				:: deltaT
INTEGER						:: i

CALL ITRANSFORM(X,U)
CALL FORCES(X,F)

F_staging=0

do i=1,nbeads
	F_staging(1)=F_staging(1)+F(i)
enddo

do i=2,nbeads
	F_staging(i)=F(i) + ( ( (i-2)/(i-1) )*F_staging(i-1) )
enddo

do i=1,nbeads
	P(i)=P(i) + F_staging(i)*deltaT
enddo
END SUBROUTINE B


SUBROUTINE A(deltaT)
USE Global, ONLY : nbeads,U,P,Mp,wp
IMPLICIT NONE
REAL(8), DIMENSION(nbeads)			:: Utmp, Ptmp
REAL(8), DIMENSION(nbeads)			:: Mp
REAL(8), INTENT(in)				:: deltaT
INTEGER						:: i


Utmp=0
Ptmp=0

Utmp(1)=U(1) + ( ( P(1)/Mp(1) )*deltaT
Ptmp(1)=P(1)

do i=2,nbeads
	Utmp(i) = cos(wp*deltaT)*U(i) + ( sin(wp*deltaT)/(wp*Mp(i)) )*P(i)
	Ptmp(i) = -wp*sin(wp*deltaT)*Mp(i)*U(i) + cos(wp*deltaT)*P(i)
enddo

U=Utmp
P=Ptmp

END SUBROUTINE A


SUBROUTINE O(deltaT)
USE Global, ONLY : nbeads,P,gamma_lang,beta,Mp
IMPLICIT NONE
REAL(8), DIMENSION(nbeads)			:: Mp,R
REAL(8), INTENT(in)				:: deltaT
INTEGER						:: i

do i=1,n_beads
	R(i)=GAUSS()
enddo

R=R-(sum(R)/(n_beads))

do i=1,nbeads
	P(i)=exp(-gamma_lang(i)*deltaT)*P(i) + sqrt( 1.- exp(-2*gamma_lang(i)*deltaT) )*sqrt(Mp(i))*R(i)
enddo

END SUBROUTINE O

REAL(8) FUNCTION GAUSS()
implicit none
REAL(8), PARAMETER	 :: a1 = 3.949846138, a3 = 0.252408784, a5 = 0.076542912, a7 = 0.008355968, a9 = 0.029899776
REAL			 :: r, r2
REAL, DIMENSION(12) 	 :: s

CALL RANDOM_NUMBER(s) ; r = (sum(s)-6.0)/4.0 ; r2 = r*r
GAUSS =  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 )*r
END FUNCTION GAUSS

END MODULE Langevin
