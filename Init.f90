MODULE Init_close

CONTAINS

SUBROUTINE INITIALIZE()
USE Global
IMPLICIT NONE
INTEGER			:: i

ALLOCATE(U(nbeads))
ALLOCATE(P(nbeads))
ALLOCATE(X(nbeads))
ALLOCATE(Mp(nbeads))
ALLOCATE(gamma_lang(nbeads))

X(:)=init_pos
U(:)=X(:)
P(:)=0.
Mp(1)=Mass_ref

do i=2,nbeads
	Mp(i)=(i/(i-1))*Mp(1)
enddo

beta=1./(Kb*Temperature)
wp=sqrt(1.*nbeads)/(beta*hbar)
gamma_lang(:)=init_gamma_lang!wp

OPEN(15,FILE="Beads.pos")
OPEN(16,FILE="16Beads.pos")

write(15,*) 0, X(:)


END SUBROUTINE INITIALIZE

SUBROUTINE FINALIZE()
USE Global
IMPLICIT NONE

DEALLOCATE(U)
DEALLOCATE(P)
DEALLOCATE(X)
DEALLOCATE(Mp)
DEALLOCATE(gamma_lang)

CLOSE(15)
CLOSE(16)

END SUBROUTINE FINALIZE


END MODULE Init_close
