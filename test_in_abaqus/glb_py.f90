module global

! aba_param.inc inclusion is commented out in uexternaldb and getoutdir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the control parameters to run the material-related routines
INTEGER NELEM, NSDV, NTERM, FACTOR, NDIR, NGP
PARAMETER (NELEM=1)
PARAMETER (NSDV=4)
PARAMETER (NTERM=60) ! 60
PARAMETER (FACTOR=1)
PARAMETER (NDIR=20 * FACTOR**2)
PARAMETER (NGP=8)

DOUBLE PRECISION  ONE, TWO, THREE, FOUR, SIX, ZERO
PARAMETER (ZERO=0.D0, ONE=1.0D0,TWO=2.0D0)
PARAMETER (THREE=3.0D0,FOUR=4.0D0,SIX=6.0D0)
DOUBLE PRECISION  HALF,THIRD
PARAMETER (HALF=0.5d0,THIRD=1.d0/3.d0)

CHARACTER(256) DIR2, DIR3
PARAMETER (DIR2='prefdir.inp')
PARAMETER (DIR3='etadir.inp')

END module global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE global_read_py(glob_vars)

    !use,intrinsic :: ISO_Fortran_env
    use global

    INTEGER, INTENT(OUT) :: glob_vars(4)

    glob_vars(1) = nsdv
    glob_vars(2) = nelem
    glob_vars(3) = ndir
    glob_vars(4) = ngp

END SUBROUTINE
