SUBROUTINE FORCES(F)
USE Global, ONLY : nbeads,Kb,X,Ma
IMPLICIT NONE
REAL(8), DIMENSION(nbeads), INTENT(out)	:: F
REAL(8)					:: h,d,freq
INTEGER					:: i
F(:)=0
h=1000._8*Kb
d=0.6*1e-10_8
freq=100.e12_8
do i=1,nbeads
	F(i)=-Ma(i)*(freq**2)*X(i)!-16._8*h*((4._8*(X(i)**3_8)/(d**4_8)) - (X(i)/(d**2_8)) ) !-Ma(i)*(freq**2)*X(i)
enddo
END SUBROUTINE FORCES
