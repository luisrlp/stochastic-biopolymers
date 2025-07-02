subroutine filpce(q0, lambda0, r0, force, dw, ddw)

    implicit none
    real(8), intent(in) :: q0, lambda0, r0
    real(8), intent(out) :: force, dw, ddw

    ! Compute the force based on the given polynomial
    force = 3628.4111205007407*q0**6 - 23035.252924661596*q0**5 + &
                    61011.76210946901*q0**4 - 86271.28882401937*q0**3 + &
                    68673.50070298258*q0**2 - 29173.565910204914*q0 + &
                    5166.433769616657
    write(*,*) "Force: ", force
    ! Compute the first derivative of the strain energy
    dw = lambda0 * r0 * force
    ddw = 1.0

end subroutine filpce

