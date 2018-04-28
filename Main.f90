!---------------------------------------------------------------------!
! Main program							      !
!---------------------------------------------------------------------!
PROGRAM Principal
USE Constants, ONLY : DP
USE Global, ONLY : pos_tot, ions_dynamics
USE Init_close
USE Distributions
USE Dynamics

IMPLICIT NONE

INTEGER							:: i,j,k,l
REAL(DP),	DIMENSION(:,:),		ALLOCATABLE 	:: proba


CALL INITIALIZE()

SELECT CASE(ions_dynamics)
	CASE("BAOAB")
		CALL BAOAB()
	CASE("Verlet")
		CALL VERLET()
END SELECT

!CALL distribution2d(pos_tot(:,1,1),pos_tot(:,1,2),proba,0.05_DP,0.05_DP,-5._DP,5._DP,-5._DP,5._DP,.TRUE.)


!DO i=1,size(proba,1)-1
!	WRITE(15,*) proba(i,:)
!	if (proba(i+1,1) .NE. proba(i,1)) THEN
!		WRITE(15,*) ""
!	endif
!ENDDO

 
CALL FINALIZE()
DEALLOCATE(proba)

END PROGRAM Principal
