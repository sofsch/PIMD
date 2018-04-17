PROGRAM Principal
USE Langevin
USE Global
USE Staging
USE Init_close
IMPLICIT NONE
INTEGER			:: i,j

CALL Read_namelist()
Call INITIALIZE

do i=1,nstep
	CALL B(dt/2)
	CALL A(dt/2)
	CALL O(dt)
	CALL A(dt/2)
	CALL B(dt/2)
	
	CALL ITRANSFORM(X,U)
	WRITE(15,*) i*dt, X(:), P(:)
	do j=1,nbeads
		write(16,*) X(j)
	enddo
enddo
 
CALL FINALIZE()

END PROGRAM Principal
