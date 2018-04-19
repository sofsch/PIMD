MODULE Langevin

CONTAINS

SUBROUTINE B(deltaT)
USE Staging
USE Global, ONLY : nbeads,P
IMPLICIT NONE
REAL(8), DIMENSION(nbeads)			:: F,F_staging
REAL(8), INTENT(in)				:: deltaT
INTEGER						:: i

CALL ITRANSFORM()
CALL FORCES(F)

F_staging(:)=0._8

do i=1,nbeads
	F_staging(1)=F_staging(1)+(F(i)/DBLE(nbeads))
enddo

do i=2,nbeads
	F_staging(i)=(F(i)/DBLE(nbeads)) + ((( (i-2._8)/(i-1._8) )*F_staging(i-1)))
enddo

do i=1,nbeads
	P(i)=P(i) + (F_staging(i)*deltaT)
enddo
END SUBROUTINE B


SUBROUTINE A(deltaT)
USE Global, ONLY : nbeads,U,P,Mp,wp
IMPLICIT NONE
REAL(8), DIMENSION(nbeads)			:: Utmp, Ptmp
REAL(8), INTENT(in)				:: deltaT
INTEGER						:: i


Utmp=0._8
Ptmp=0._8

Utmp(1)=U(1) + ( ( P(1)/Mp(1) )*deltaT )
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
REAL(8), DIMENSION(nbeads)			:: R
REAL(8), INTENT(in)				:: deltaT
INTEGER						:: i

do i=1,nbeads
	R(i)=GAUSS()
enddo

!R=R-(sum(R)/(nbeads))

do i=1,nbeads
	P(i)=exp(-gamma_lang(i)*deltaT)*P(i) + sqrt( 1._8- exp(-2._8*gamma_lang(i)*deltaT) )*sqrt(Mp(i)/beta)*R(i)
enddo

END SUBROUTINE O

REAL(8) FUNCTION GAUSS()
IMPLICIT NONE
REAL(8), PARAMETER	 :: a1 = 3.949846138_8, a3 = 0.252408784_8, a5 = 0.076542912_8, a7 = 0.008355968_8, a9 = 0.029899776_8
REAL			 :: r, r2
REAL, DIMENSION(12) 	 :: s

CALL RANDOM_NUMBER(s) ; r = (sum(s)-6.0_8)/4.0_8 ; r2 = r*r
GAUSS =  (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 )*r
END FUNCTION GAUSS

END MODULE Langevin
