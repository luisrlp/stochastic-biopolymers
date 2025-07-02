SUBROUTINE run_umat_py(NPROPS, NDI, PROPS, DFGRD1, STRESS, TIME, DTIME, KSTEP)

use,intrinsic :: ISO_Fortran_env
use global
!INCLUDE 'aba_param.inc'

!C     ADD COMMON BLOCKS HERE IF NEEDED ()
!C      COMMON /KBLOCK/KBLOCK
      COMMON /KFIL/MF0
      COMMON /KFILR/RW
      COMMON /KFILP/PREFDIR
      !DOUBLE PRECISION PREFDIR(NELEM,4)
      DOUBLE PRECISION PREFDIR(1,4)

PARAMETER(NTENS = 6, NSTATEV = NSDV, NSHR=3)
PARAMETER(NOEL = 1, NPT = 8)

INTEGER, INTENT(IN) :: NPROPS, NDI, KSTEP
DOUBLE PRECISION, INTENT(IN) :: PROPS(NPROPS), DFGRD1(NDI,NDI), TIME(2), DTIME
DOUBLE PRECISION, INTENT(OUT) :: STRESS(NTENS)
!
integer ii
CHARACTER*8 CMNAME, stri
CHARACTER*40 FILENAME
DIMENSION STATEV(NSTATEV),DDSDDE(NTENS,NTENS),DDSDDT(NTENS),      &
DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),PREDEF(1),DPRED(1),            &
COORDS(3),DROT(3,3),DFGRD0(3,3)


i=1.0d0
j=1.0d0
DO i=1,NTENS
    DO j=1,NTENS
        DDSDDE(i,j)=0.0D0
    ENDDO
    STRESS(i)=0.0D0
ENDDO
!
STATEV=0.D0
erf=0.d0
RHO=0.D0
!
!################################################################################################!


!time(1)=0.d0
!time(2)=0.d0
!dtime = 0.1d0
!kstep = 1
!
!!!
call UEXTERNALDB(0,0,time,0.D0,0,0)
!
CALL UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT, DRPLDE,DRPLDT,STRAN,     &
DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,NSTATEV,PROPS,  &
NPROPS,COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

 write(*,*) DFGRD1
 write(*,*)
 write(*,*) STRESS
 write(*,*)
 write(*,*) DDSDDE
close(150)
!################################################################################################!

END SUBROUTINE run_umat_py