SUBROUTINE sdvwrite(det,etac_sdv,statev)
!>    VISCOUS DISSIPATION: WRITE STATE VARS
use global
implicit none

INTEGER :: pos1, i
!
DOUBLE PRECISION, INTENT(IN)             :: det
DOUBLE PRECISION, INTENT(IN)             :: etac_sdv(nsdv-1)
DOUBLE PRECISION, INTENT(OUT)            :: statev(nsdv)
!
pos1=0
statev(pos1+1)=det
IF (nsdv .GE. 2) THEN
    DO i = pos1+2, nsdv
        statev(i)=etac_sdv(i-1)
    END DO
END IF

RETURN

END SUBROUTINE sdvwrite
