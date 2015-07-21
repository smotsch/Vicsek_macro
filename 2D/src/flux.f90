module flux

  !- computation of the value at the interfaces (+ upwind)
  !- Contains :
  !-    Upwind_x, FluxSemi_x, FluxCons_x, FluxSplit_x

  use elementary                ! to have AngleModulo, ThetaMean
  use matrix                    ! to have all the matrix necessary

  implicit none

  
contains
  
  Function FluxSplit_x(rho_l,rho_r,u_l,u_r,v_l,v_r,c1,c2,ld,degree)
    ! value of the flux (ρ,ρu,ρv) for the splitting method
    implicit none
    Double precision, intent(in)                 :: rho_l, rho_r, u_l, u_r, v_l, v_r
    Double precision, intent(in)                 :: c1,c2, ld
    Integer, Intent(in)                          :: degree
    Double precision, Dimension(3,3)             :: mat_abs_A_roe
    Double precision, Dimension(3,1)             :: Temp
    Double precision, Dimension(3)               :: FluxSplit_x
    !- the function
    !- the matrix abs_A_roe
    mat_abs_A_roe = Abs_A_roe_Poly(rho_l,rho_r,u_l,u_r,v_l,v_r,c1,c2,ld,degree)
    !- the value
    Temp = 1/2d0*reshape( (/ &
         c1*rho_l*u_l + c1*rho_r*u_r , &
         c2*( rho_l*u_l**2  + rho_r*u_r**2 ) + ld*(rho_l+rho_r) , &
         c2*( rho_l*u_l*v_l + rho_r*u_r*v_r) /), (/3,1/) ) &
         - &
         1/2d0*matmul(mat_abs_A_roe, &
         reshape( (/ &
         rho_r-rho_l , &
         rho_r*u_r - rho_l*u_l , &
         rho_r*v_r - rho_l*v_l /), (/3,1/) ))

    FluxSplit_x = (/ Temp(1,1) , Temp(2,1), Temp(3,1) /)

  End Function FluxSplit_x

end module flux
