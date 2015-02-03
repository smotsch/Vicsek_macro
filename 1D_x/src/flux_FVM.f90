module flux_FVM

  !- computation of the value at the interfaces (+ upwind)
  !- Contains :
  !-    Upwind_x, FluxSemi_x, FluxCons_x, FluxSplit_x

  use elementary             ! to have AngleModulo, ThetaMean
  use matrix_FVM             ! to have all the matrix necessary

  implicit none

  
contains

  function Upwind_x( rho_im,rho_i,rho_ip, theta_im,theta_i,theta_ip, c1,c2,ld, dt,dx)
    !- upwind in x
    implicit none
    Double Precision, intent(in)               :: c1,c2,ld, dt,dx
    Double Precision, intent(in)               :: rho_im, rho_i, rho_ip
    Double Precision, intent(in)               :: theta_im, theta_i, theta_ip
    Double Precision, Dimension(2,1)           :: Diff_im, Diff_ip, U_i, U_i_p
    Double Precision, Dimension(2)             :: Upwind_x
    !-- init state
    U_i    = reshape( (/ rho_i, AngleModulo(theta_i) /), (/2,1/) )
    Diff_im = reshape( (/ rho_i-rho_im, AngleModulo(theta_i-theta_im) /), (/2,1/) )
    Diff_ip = reshape( (/ rho_ip-rho_i, AngleModulo(theta_ip-theta_i) /), (/2,1/) )
    !- update
    U_i_p = U_i - dt/dx*( matmul(A_plus(rho_i,theta_i,c1,c2,ld) , Diff_im) &
         + matmul(A_minus(rho_i,theta_i,c1,c2,ld) , Diff_ip) )
    Upwind_x = (/ U_i_p(1,1) , U_i_p(2,1) /)
  end function Upwind_x

  
  function FluxSemi_x(rho_l,rho_r,theta_l,theta_r,c1,c2,ld)
    !--- Flux in ρ for the semi_conservative method
    implicit none
    Double precision, intent(in)                 :: c1,c2, ld
    Double precision, intent(in)                 :: rho_l, rho_r, theta_l, theta_r
    Double precision                             :: rho_m, theta_m
    Double precision, Dimension(2,2)             :: abs_A_m
    Double precision                             :: FluxSemi_x
    !- init
    rho_m   = 5d-1*(rho_l+rho_r)
    theta_m = ThetaMean(theta_l,theta_r)
    abs_A_m = abs_A(rho_m,theta_m,c1,c2,ld)
    !-- Finally
    !FluxSemi_x = rho_m*cos(theta_m) - abs_A_m(1,1)*(rho_r-rho_l)/2d0 &
    !     - abs_A_m(1,2)*AngleModulo(theta_r-theta_l)/2d0
    FluxSemi_x = 1/2d0*(c1*rho_l*cos(theta_l) + c1*rho_r*cos(theta_r)) - abs_A_m(1,1)*(rho_r-rho_l)/2d0 &
         - abs_A_m(1,2)*AngleModulo(theta_r-theta_l)/2d0
    
    ! Old technic (not good...)
    !!  !-- matrice signe
    !!  sgn_A_m = sgn_A(rho_m,theta_m,c,ld)
    !!  ! valeur de U alpha_2*r_2 + beta_1*r_1
    !!  U_x = reshape( (/ (rho_l+rho_r)/2d0 , ThetaMean(theta_l,theta_r) /) , (/2,1/)) - &
    !!       matmul(sgn_A_m , &
    !!       reshape( (/ (rho_r-rho_l)/2d0, AngleModulo(theta_r-theta_l)/2d0 /) , (/2,1/)))
    !!  !- Finally
    !!  FluxSemi_x = U_x(1,1)*cos(U_x(2,1))
  end function FluxSemi_x
  
  
  function FluxCons_x(rho_l,rho_r,theta_l,theta_r,c1,c2,ld,degree)
    !---  flux at the interface in (ρ,θ) for the conservative method
    implicit none
    Double precision, intent(in)                 :: rho_l, rho_r, theta_l, theta_r
    Double precision, intent(in)                 :: c1,c2, ld
    Integer, Intent(in)                          :: degree
    Double precision                             :: rho_m, theta_m
    Double precision, Dimension(2,2)             :: mat_abs_Acons
    Double precision, Dimension(2,1)             :: Diff_lr
    Double precision, Dimension(2)               :: FluxCons_x
    !- init
    rho_m   = 5d-1*(rho_l+rho_r)
    theta_m = ThetaMean(theta_l,theta_r)
    !--- l'approximation de |A1c_m| ---!
    mat_abs_Acons = abs_Acons_Poly(rho_m,theta_m,c1,c2,ld,degree)
    !--- value at the interface ---!
    Diff_lr = matmul(mat_abs_Acons, &
         reshape( (/rho_r-rho_l, f1(theta_r)-f1(theta_l) /),(/2,1/) ))

    FluxCons_x = 1d0/2*( (/ &
         c1*rho_l*cos(theta_l) + c1*rho_r*cos(theta_r) , &
         c2*f2(theta_l)-ld*log(rho_l) + c2*f2(theta_r)-ld*log(rho_r) /) &
         -  (/ &
         Diff_lr(1,1), &
         Diff_lr(2,1)/) &
         )
    !-ou bien
    !FluxCons_x = (/ c1*rho_m*cos(theta_m) , c2*f2(theta_m)-ld*log(rho_m) /) &
    !     - (/Diff_lr(1,1),Diff_lr(2,1)/)
  end function FluxCons_x

  
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

end module flux_FVM
