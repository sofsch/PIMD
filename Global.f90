MODULE Global
IMPLICIT NONE

INTEGER, PARAMETER		:: nbeads=16,nstep=1000000
REAL(8), DIMENSION(nbeads)	:: U,P,X,Mp,gamma_lang
REAL(8)				:: dt,wp,T, beta
REAL(8), PARAMETER		:: Kb=1.38064852e-23, hbar=1.054571800e-34

END MODULE Global
