!---------------------------------------------------------------------!
! Module Global :	contains variables			      !
!---------------------------------------------------------------------!
MODULE Global

IMPLICIT NONE

INTEGER					:: nbeads,nat,nstep,ntyp
REAL(8), DIMENSION(:,:,:), ALLOCATABLE	:: U,P,tau
REAL(8), DIMENSION(:,:), ALLOCATABLE	:: Mp,Ma,force_constraint
REAL(8), DIMENSION(:), ALLOCATABLE	:: gamma_lang,Mass
REAL(8)					:: dt,wp,Temperature, beta, init_gamma_lang
REAL(8), dimension(:,:), ALLOCATABLE	:: init_pos
REAL(8), PARAMETER			:: Kb=1.38064852e-23_8, hbar=1.054571800e-34_8
CHARACTER(50)				:: output

END MODULE Global
