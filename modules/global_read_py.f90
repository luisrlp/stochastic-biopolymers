SUBROUTINE global_read_py(glob_vars)

    !use,intrinsic :: ISO_Fortran_env
    use global

    INTEGER, INTENT(OUT) :: glob_vars(4)

    glob_vars(1) = nsdv
    glob_vars(2) = nelem
    glob_vars(3) = ndir
    glob_vars(4) = ngp

END SUBROUTINE
