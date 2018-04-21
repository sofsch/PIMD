!---------------------------------------------------------------------!
! Module Distributions						      !
!---------------------------------------------------------------------!
MODULE Distributions

CONTAINS

SUBROUTINE DISTRIBUTIONS1D(array,proba,dp,xmin,xmax)

	IMPLICIT NONE
	
	REAL(8), DIMENSION(:,:), INTENT(in) 			:: array
	REAL(8), INTENT(in) 					:: xmin, xmax
	REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(out) 	:: proba
	REAL(8) 						:: dp,norm
	INTEGER 						:: nb_step,i,j
	
	nb_step=nint((xmax-xmin)/dp)
	ALLOCATE(proba(nb_step+1,2))
	proba(:,:)=0.
	norm=0.
		
	
	DO i=1,nb_step+1
		DO j=1,size(array,1)
			IF ( array(j,1) <= (xmin + (i-1)*dp + dp)  .AND. array(j,1) > (xmin + (i-1)*dp) ) THEN
				proba(i,2)=proba(i,2)+(1./array(j,1))
			ENDIF
		ENDDO
		proba(i,1)=xmin + (i-1)*dp
	ENDDO
	norm=sum(proba(:,2)*dp)
	proba(:,2)=proba(:,2)/norm
	
END SUBROUTINE DISTRIBUTIONS1D


SUBROUTINE distribution2d(array1,array2,proba,dp1,dp2,min_1,max_1,min_2,max_2,gnu)

	IMPLICIT NONE
	
	REAL(8), DIMENSION(:), INTENT(in) 			:: array1,array2
	REAL(8), INTENT(in) 					:: min_1,min_2,max_1,max_2
	REAL(8) 						:: xmin1, xmax1, xmin2, xmax2
	REAL(8), INTENT(in)	 				:: dp1,dp2
	LOGICAL, INTENT(in) 					:: gnu
	REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(out) 	:: proba
	REAL(8), DIMENSION(:,:), ALLOCATABLE 			:: proba_prov
	REAL(8) 						:: norme
	INTEGER 						:: nb_step1,nb_step2,i,j,x,y,k
	INTEGER 						:: start1,stop1,start2,stop2
	
	norme=0
	IF(min_1 == 0. .AND. min_2 == 0. .AND. max_1 == 0. .AND. max_2 == 0.) THEN
		xmin1=minval(array1(:))
		xmax1=maxval(array1(:))
		xmin2=minval(array2(:))
		xmax2=maxval(array2(:))
	else
		xmin1=min_1
		xmax1=max_1
		xmin2=min_2
		xmax2=max_2
	ENDIF
	

	nb_step1=nint(((xmax1-xmin1)/dp1)+1)
	nb_step2=nint(((xmax2-xmin2)/dp2)+1)
	start1= nint(xmin1/dp1)
	start2= nint(xmin2/dp2)
	stop1= nint(xmax1/dp1)
	stop2= nint(xmax2/dp2)
	ALLOCATE(proba_prov(start1:stop1,start2:stop2))
	proba_prov(:,:)=0

	k=0
	DO i=1,size(array1,1)
			x = nint((array1(i))/dp1)
			y = nint((array2(i))/dp2)
			proba_prov(x,y)=proba_prov(x,y)+1
	ENDDO
	norme = sum(proba_prov(:,:))*dp1*dp2
	proba_prov(:,:)=proba_prov(:,:)/norme
	
	IF (gnu .EQV. .TRUE.) THEN
	ALLOCATE(proba((nb_step1)*(nb_step2),3))
	proba(:,:)=0
	k=0
	DO i=0, size(proba_prov,1)-1
		DO j=0, size(proba_prov,2)-1
			k=k+1
			proba( k,1)= xmin1 + ((i)*dp1)
			proba( k,2)= xmin2 + ((j)*dp2)
			proba( k,3)= proba_prov(i+start1,j+start2)		
		ENDDO
	ENDDO

		
	else IF (gnu .NEQV. .TRUE.) THEN		
		ALLOCATE(proba(nb_step1,nb_step2))
		proba=proba_prov/maxval(proba_prov)
	ENDIF
	
	DEALLOCATE(proba_prov)
	

END SUBROUTINE distribution2d

END MODULE Distributions
