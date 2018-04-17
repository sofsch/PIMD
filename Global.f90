MODULE Global
IMPLICIT NONE

INTEGER					:: nbeads,nstep
REAL(8), DIMENSION(:), ALLOCATABLE	:: U,P,X,Mp,gamma_lang,Ma
REAL(8)					:: dt,wp,Temperature, beta, Mass_ref, init_pos
REAL(8)					:: init_gamma_lang
REAL(8), PARAMETER			:: Kb=1.38064852e-23, hbar=1.054571800e-34

END MODULE Global
