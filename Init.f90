MODULE Init_close

CONTAINS

SUBROUTINE INITIALIZE()
USE Global
IMPLICIT NONE
INTEGER			:: i

ALLOCATE(U(nbeads))
ALLOCATE(P(nbeads))
ALLOCATE(X(nbeads))
ALLOCATE(Mp(nbeads),Ma(nbeads))
ALLOCATE(gamma_lang(nbeads))

X(:)=init_pos
U(:)=X(:)
P(:)=0.
Mp(1)=Mass_ref
Ma(1)=0.
do i=2,nbeads
	Mp(i)=(i/(i-1))*Mp(1)
	Ma(i)=Mp(i)
enddo

beta=1./(Kb*Temperature)
wp=sqrt(1.*nbeads)/(beta*hbar)
gamma_lang(:)=wp!init_gamma_lang

OPEN(15,FILE="Beads.pos")
OPEN(16,FILE="16Beads.pos")

!write(15,*) 0, X(:), P(:)


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
CLOSE(16)

END SUBROUTINE FINALIZE


END MODULE Init_close
