!---------------------------------------------------------------------!
! Module Global :	contains variables			      !
!---------------------------------------------------------------------!
MODULE Global
USE Constants, ONLY : DP
IMPLICIT NONE

INTEGER						:: nbeads,nat,nstep,ntyp
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE		:: U,P,tau
REAL(DP), DIMENSION(:,:), ALLOCATABLE		:: Mp,Ma,force_constraint
REAL(DP), DIMENSION(:), ALLOCATABLE		:: gamma_lang,Mass
REAL(DP)					:: dt,wp,Temperature, beta
REAL(DP)					:: KineticEnergy
CHARACTER(50)					:: output

END MODULE Global
