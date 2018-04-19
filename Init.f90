MODULE Init_close

CONTAINS

SUBROUTINE INITIALIZE()
USE Global
USE Staging
IMPLICIT NONE
INTEGER			:: i

ALLOCATE(U(nbeads))
ALLOCATE(P(nbeads))
ALLOCATE(X(nbeads))
ALLOCATE(Mp(nbeads),Ma(nbeads))
ALLOCATE(gamma_lang(nbeads))

X(:)=init_pos
CALL TRANSFORM()
P(:)=0._8
Mp(1)=Mass_ref
Ma(1)=0._8
output=trim(output)
do i=2,nbeads
	Mp(i)=(i/(i-1._8))*Mp(1)
	Ma(i)=Mp(i)
enddo

beta=1._8/(Kb*Temperature)
wp=sqrt(1._8*nbeads)/(beta*hbar)
if (init_gamma_lang .ne. 0) then
	gamma_lang(:)=init_gamma_lang
else
	gamma_lang(:)=wp
endif

OPEN(15,FILE=TRIM(output)//".dist")
!OPEN(16,FILE=TRIM(output)//".pos")



END SUBROUTINE INITIALIZE

SUBROUTINE FINALIZE()
USE Global
IMPLICIT NONE

DEALLOCATE(U)
DEALLOCATE(P)
DEALLOCATE(X)
DEALLOCATE(Mp,Ma)
DEALLOCATE(gamma_lang)

CLOSE(15)
!CLOSE(16)

END SUBROUTINE FINALIZE


END MODULE Init_close
