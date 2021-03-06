module matrix_FVM

  !- contains
  !-     Det
  !-     ValeurPropre_x
  !-     Matrix_A, Abs_A, A_plus, A_minus
  !-     Matrix_Acons, Abs_Acons, Abs_Acons_Poly
  !-     Matrix_A_roe, Abs_A_roe, Abs_A_roe_Poly
  
  use elementary              ! for PolyAbs
  
  implicit none
  
contains

  !-------  tools
  !--------------
  function Det(A)
    !-- Determinant matrix 2x2
    implicit none
    Double Precision, Dimension(2,2), intent(in) :: A
    Double Precision                             :: Det
    !- the function
    Det = A(1,1)*A(2,2) - A(2,1)*A(1,2)
  End function Det
  
  subroutine ValeurPropre_x(theta,c1,c2,ld, vp1_x,vp2_x)
    !-- vp of A
    implicit none
    Double Precision, intent(in)        :: theta, c1,c2,ld
    Double Precision                    :: discri
    Double Precision, intent(out)       :: vp1_x,vp2_x
    discri = (c1-c2)**2*cos(theta)**2 + 4*c1*ld*sin(theta)**2
    vp1_x = 5d-1*( (c1+c2)*cos(theta) - sqrt(discri) )
    Vp2_x = 5d-1*( (c1+c2)*cos(theta) + sqrt(discri) )
  end subroutine ValeurPropre_x

  function ValeurPropreSplit_x(u,c1,c2,ld)
    implicit none
    Double precision, intent(in)                 :: u,c1,c2,ld
    Double precision, dimension(3)               :: ValeurPropreSplit_x
    Double precision                             :: sq_Delta
    sq_Delta = sqrt( (c2**2-c1*c2)*u**2 + c1*ld )
    ValeurPropreSplit_x(1) = c2*u - sq_Delta
    ValeurPropreSplit_x(2) = c2*u
    ValeurPropreSplit_x(3) = c2*u + sq_Delta
  endfunction ValeurPropreSplit_x

  
  !------------------------------------------------------------------------!
  !--------------------    1) A splitting method     ----------------------!
  !------------------------------------------------------------------------!

  function Abs_A_roe(rho_l,rho_r,u_l,u_r,v_l,v_r, c1,c2,ld)
    !- absolute value of the Roe matrix
    implicit none
    Double Precision, intent(in)               :: rho_l,rho_r,u_l,u_r,v_l,v_r
    Double Precision, intent(in)               :: c1,c2,ld
    Double precision                           :: u_m,v_m
    Double precision                           :: z1, z3, det_P_inv
    Double Precision, Dimension(3,3)           :: P, P_inv, abs_D_roe
    Double Precision, Dimension(3,3)           :: Abs_A_roe
    !- entropy fix
    Double Precision, Dimension(3,1)           :: alpha,beta,q_r_1,q_r_2
    Double Precision, Dimension(3)             :: vpL,vpM,vpR
    Double precision                           :: vp1,vp2,vp3
    Double precision                           :: vp1_fix,vp2_fix,vp3_fix
    !--- (u,v) at the middle
    u_m = (sqrt(rho_l)*u_l + sqrt(rho_r)*u_r)/(sqrt(rho_l) + sqrt(rho_r))
    v_m = (sqrt(rho_l)*v_l + sqrt(rho_r)*v_r)/(sqrt(rho_l) + sqrt(rho_r))
    !--- eigenvalues (vp)
    vpM = ValeurPropreSplit_x(u_m,c1,c2,ld)
    vp1 = vpM(1); vp2 = vpM(2); vp3 = vpM(3);
    !-- eigenvectors [V_1,V_2,V_3]
    z1 = (c1*c2*u_m*v_m-c2*v_m*vp1)/(c2*u_m-vp1)
    z3 = (c1*c2*u_m*v_m-c2*v_m*vp3)/(c2*u_m-vp3)
    P_inv = reshape( (/ &
         c1  , 0d0 , c1  , &
         vp1 , 0d0 , vp3 , &
         z1  , 1d0 , z3 /) , (/3,3/) )
    det_P_inv = c1*(-vp3 + vp1)
    P = 1/det_P_inv*reshape( (/ &
         -vp3  , c1 , 0d0 , &
         -(vp1*z3-vp3*z1) , c1*(z3-z1) , -c1*(vp3-vp1), &
         vp1   , -c1 , 0d0 /) , (/3,3/) )
    !- α,β coordinates of U_l,U_r in V_1,V_2,V_3
    alpha = matmul(P,reshape( (/rho_l, rho_l*u_l, rho_l*v_l/),(/3,1/)))
    beta  = matmul(P,reshape( (/rho_r, rho_r*u_r, rho_r*v_r/),(/3,1/)))
    !- q_r_1 (=q_l_2), q_r_2 (=q_l_3)
    q_r_1 = matmul(P_inv,reshape( (/beta(1,1),alpha(2,1),alpha(3,1) /),(/3,1/)))
    q_r_2 = matmul(P_inv,reshape( (/beta(1,1),beta(2,1) ,alpha(3,1) /),(/3,1/)))
    !-------------    
    !- entropy-fix
    !-------------
    !  We increase the eigenvalues to avoid sonic rarefaction wave (LLF method)
    !      λ_m = max(|λ_l|,|λ_r|)
    ! for each eigenvalue.
    vpL = ValeurPropreSplit_x(u_l,c1,c2,ld)
    vpR = ValeurPropreSplit_x(u_r,c1,c2,ld)
    vp1_fix = max(abs(vpL(1)),abs(vpR(1)))
    vp2_fix = max(abs(vpL(2)),abs(vpR(2)))
    vp3_fix = max(abs(vpL(3)),abs(vpR(3)))
    ! the absolute value of D
    abs_D_roe = reshape( (/ &
         vp1_fix, 0d0    , 0d0, &
         & 0d0  , vp2_fix, 0d0, &
         & 0d0  , 0d0    , vp3_fix /) , (/3,3/) )
    !-- the mix matrix A^+
    Abs_A_roe = matmul(transpose(P_inv),matmul(abs_D_roe,transpose(P)))       

  end function Abs_A_roe

  function Abs_A_roe_Poly(rho_l,rho_r,u_l,u_r,v_l,v_r,c1,c2,ld,degree)
    !---  compute the approximate (or not) of |A_roe| using (or not) the polynome scheme
    implicit none
    Double precision, intent(in)                 :: rho_l,rho_r,u_l,u_r,v_l,v_r
    Double precision, intent(in)                 :: c1,c2, ld
    Double precision                             :: u,v
    Integer, Intent(in)                          :: degree
    Double precision                             :: vp1, vp2, vp3, sq_Delta
    Double precision, Dimension(3,3)             :: A_roe
    Double precision                             :: a_m, a_p, a_max
    Double precision, Dimension(3,3)             :: Abs_A_roe_Poly
    !- (u,v) at the middle
    u = (sqrt(rho_l)*u_l + sqrt(rho_r)*u_r)/(sqrt(rho_l) + sqrt(rho_r))
    v = (sqrt(rho_l)*v_l + sqrt(rho_r)*v_r)/(sqrt(rho_l) + sqrt(rho_r))
    !--- l'approximation de |A| ---!
    select case(degree)
    case(0,1,2)
       ! approximation avec methode P
       !--  init  --!
       ! valeurs propres
       sq_Delta = sqrt( (c2**2-c1*c2)*u**2 + c1*ld )
       vp1  = c2*u - sq_Delta
       vp2  = c2*u
       vp3  = c2*u + sq_Delta
       ! la jacobienne
       A_roe = reshape( (/ 0d0 , c1, 0d0, &
            ld - c2*u**2 , 2*c2*u , 0d0, &
            -c2*u*v, c2*v, c2*u /) , (/3,3/) )
       A_roe = transpose(A_roe)
       a_m = min(vp1,vp2,vp3)
       a_p = max(vp1,vp2,vp3)
       if (abs(a_p)>abs(a_m)) then
          a_max = a_p
       else
          a_max = a_m
       endif
       !-- The matrix
       Call PolyAbs(A_roe,a_m,a_p,a_max,degree,Abs_A_roe_poly)
    case default
       ! the exact absolute matrix
       Abs_A_roe_Poly = Abs_A_roe(rho_l,rho_r,u_l,u_r,v_l,v_r,c1,c2,ld)
    end select
    
  end function Abs_A_roe_Poly


  !-------------------------------------------------------------------------------!
  !--------------------    2) A the conservative method     ----------------------!
  !-------------------------------------------------------------------------------!

  function Matrix_Acons(rho,theta,c1,c2,ld)
    !- matrix Ac (the matrix involved in the conservative scheme)
    implicit none
    Double Precision, intent(in)               :: c1,c2,ld
    Double Precision, intent(in)               :: rho,theta
    Double Precision, Dimension(2,2)           :: Matrix_Acons
    !-the function
    Matrix_Acons = reshape( (/ c1*cos(theta) , -c1*rho*sin(theta)**2 , &
         -ld/rho , c2*cos(theta) /), (/2,2/) )
    Matrix_Acons = transpose(Matrix_Acons)
  end function Matrix_Acons
  
  function Abs_Acons(rho,theta,c1,c2,ld)
    !- matrix abs(A): |A| = 1/2 ( A^+ - A^- )
    !                  A  = 1/2 ( A^+ + A^- )
    implicit none
    Double Precision, intent(in)               :: c1,c2,ld
    Double Precision, intent(in)               :: rho,theta
    Double precision                           :: vp1, vp2
    Double Precision, Dimension(2,2)           :: P, P_inv, abs_Dcons
    Double Precision, Dimension(2,2)           :: Abs_Acons
    ! init
    Call ValeurPropre_x(theta,c1,c2,ld, vp1,vp2)
    if (vp1*vp2>0) Then
       !- same sign: |A| = +/- A
       Abs_Acons = sign(1d0,vp1)*Matrix_Acons(rho,theta,c1,c2,ld)
    else
       !- decomposition
       P_inv = reshape( (/ c1*rho*sin(theta)**2 , c2*cos(theta)-vp2, &
            c1*cos(theta)-vp1 , ld/rho /), (/2,2/))
       P     = 1/Det(P_inv)*reshape( (/ ld/rho , -c2*cos(theta)+vp2, &
            -c1*cos(theta)+vp1 , c1*rho*sin(theta)**2  /), (/2,2/))
       abs_Dcons = reshape( (/ -vp1, 0d0, 0d0, vp2 /) , (/2,2/) )
       !-- |A|
       Abs_Acons = matmul(transpose(P_inv),matmul(abs_Dcons,transpose(P)))
    end if
  end function Abs_Acons

  function Abs_Acons_Poly(rho,theta,c1,c2,ld,degree)
    !---  compute the approximate (or not) of |Acons| using (or not) polynome scheme
    implicit none
    Double precision, intent(in)                 :: c1,c2, ld
    Integer, Intent(in)                          :: degree
    Double precision, intent(in)                 :: rho, theta
    Double precision                             :: vp1, vp2
    Double precision, Dimension(2,2)             :: Acons
    Double precision                             :: a_m, a_p, a_max
    Double precision, Dimension(2,2)             :: Abs_Acons_Poly
    !--- l'approximation de |A| ---!
    select case(degree)
    case(0,1,2)
       ! approximation with polynomial method
       Call ValeurPropre_x(theta,c1,c2,ld, vp1,vp2)
       Acons = Matrix_Acons(rho,theta,c1,c2,ld)
       a_m = vp1; a_p = vp2
       if (abs(a_p)>abs(a_m)) Then
          a_max = a_p
       else
          a_max = a_m
       endif
       Call PolyAbs(Acons,a_m,a_p,a_max,degree,Abs_Acons_Poly)
    case default
       ! the exact absolute matrix
       Abs_Acons_Poly = Abs_Acons(rho,theta,c1,c2,ld)
    end select
    
  end function Abs_Acons_Poly


  !-------------------------------------------------------------------------------!
  !------------------    3-4) A upwind or semiconservative     -------------------!
  !-------------------------------------------------------------------------------!

  function Matrix_A(rho,theta,c1,c2,ld)
    implicit none
    Double Precision, intent(in)               :: c1,c2,ld
    Double Precision, intent(in)               :: rho,theta
    Double Precision, Dimension(2,2)           :: Matrix_A
    Matrix_A = reshape( (/ c1*cos(theta) , -c1*rho*sin(theta) , &
         -ld*sin(theta)/rho , c2*cos(theta) /), (/2,2/) )
    Matrix_A = transpose(Matrix_A)
  end function Matrix_A
  
  function Abs_A(rho,theta,c1,c2,ld)
    !- matrix abs(A): |A| = 1/2 ( A^+ - A^- )
    !                  A  = 1/2 ( A^+ + A^- )
    implicit none
    Double Precision, intent(in)               :: rho,theta,c1,c2,ld
    Double precision                           :: vp1, vp2
    Double Precision, Dimension(2,2)           :: P, P_inv, D
    Double Precision, Dimension(2,2)           :: Abs_A
    ! init
    Call ValeurPropre_x(theta,c1,c2,ld, vp1,vp2)
    if (vp1*vp2>0) Then
       !- same sign: |A| = +/- A
       Abs_A = sign(1d0,vp1)*Matrix_A(rho,theta,c1,c2,ld)
    else
       !- decomposition
       P_inv = reshape( (/ c1*rho*sin(theta) , c2*cos(theta)-vp2, &
            c1*cos(theta)-vp1 , ld*sin(theta)/rho /) , (/2,2/))
       P  = 1/Det(P_inv)*reshape( (/ ld*sin(theta)/rho , -c2*cos(theta)+vp2, &
            -c1*cos(theta)+vp1 , c1*rho*sin(theta) /) , (/2,2/) )
       D = reshape( (/ -vp1, 0d0, 0d0, vp2 /) , (/2,2/) )
       !-- |A|
       Abs_A = matmul(transpose(P_inv),matmul(D,transpose(P)))
    end if
  end function Abs_A

  function A_plus(rho,theta,c1,c2,ld)
    implicit none
    Double Precision, intent(in)               :: rho,theta,c1,c2,ld
    Double Precision, Dimension(2,2)           :: A_plus
    A_plus = 1d0/2*( Abs_A(rho,theta,c1,c2,ld) + Matrix_A(rho,theta,c1,c2,ld) )
  end function A_plus

  function A_minus(rho,theta,c1,c2,ld)
    implicit none
    Double Precision, intent(in)               :: rho,theta,c1,c2,ld
    Double Precision, Dimension(2,2)           :: A_minus
    A_minus = 1d0/2*( Matrix_A(rho,theta,c1,c2,ld) - Abs_A(rho,theta,c1,c2,ld) )
  end function A_minus





  function Abs_A_roe_bidon(rho_l,rho_r,u_l,u_r,v_l,v_r, c1,c2,ld)
    !- absolute value of the Roe matrix
    ! bidon: l'entropie fix est trop petit
    implicit none
    Double Precision, intent(in)               :: rho_l,rho_r,u_l,u_r,v_l,v_r
    Double Precision, intent(in)               :: c1,c2,ld
    Double precision                           :: u_m,v_m
    Double precision                           :: z1, z3, det_P_inv
    Double Precision, Dimension(3,3)           :: P, P_inv, abs_D_roe
    Double Precision, Dimension(3,3)           :: Abs_A_roe_bidon
    !- entropy fix
    Double Precision, Dimension(3,1)           :: alpha,beta,q_r_1,q_r_2
    Double Precision, Dimension(3)             :: vp_m,vp_temp
    Double precision                           :: vp1,vp2,vp3
    Double precision                           :: vp1_l,vp2_l,vp3_l
    Double precision                           :: vp1_r,vp2_r,vp3_r
    Double precision                           :: vp1_fix,vp2_fix,vp3_fix
    !--- (u,v) at the middle
    u_m   = (sqrt(rho_l)*u_l + sqrt(rho_r)*u_r)/(sqrt(rho_l) + sqrt(rho_r))
    v_m   = (sqrt(rho_l)*v_l + sqrt(rho_r)*v_r)/(sqrt(rho_l) + sqrt(rho_r))
    !--- eigenvalues (vp)
    vp_m = ValeurPropreSplit_x(u_m,c1,c2,ld)
    vp1  = vp_m(1); vp2 = vp_m(2); vp3 = vp_m(3);
    !-- eigenvectors [V_1,V_2,V_3]
    z1 = (c1*c2*u_m*v_m-c2*v_m*vp1)/(c2*u_m-vp1)
    z3 = (c1*c2*u_m*v_m-c2*v_m*vp3)/(c2*u_m-vp3)
    P_inv = reshape( (/ &
         c1  , 0d0 , c1  , &
         vp1 , 0d0 , vp3 , &
         z1  , 1d0 , z3 /) , (/3,3/) )
    det_P_inv = c1*(-vp3 + vp1)
    P = 1/det_P_inv*reshape( (/ &
         -vp3  , c1 , 0d0 , &
         -(vp1*z3-vp3*z1) , c1*(z3-z1) , -c1*(vp3-vp1), &
         vp1   , -c1 , 0d0 /) , (/3,3/) )
    !- α,β coordinates of U_l,U_r in V_1,V_2,V_3
    alpha = matmul(P,reshape( (/rho_l, rho_l*u_l, rho_l*v_l/),(/3,1/)))
    beta  = matmul(P,reshape( (/rho_r, rho_r*u_r, rho_r*v_r/),(/3,1/)))
    !- q_r_1 (=q_l_2), q_r_2 (=q_l_3)
    q_r_1 = matmul(P_inv,reshape( (/beta(1,1),alpha(2,1),alpha(3,1) /),(/3,1/)))
    q_r_2 = matmul(P_inv,reshape( (/beta(1,1),beta(2,1) ,alpha(3,1) /),(/3,1/)))
    !-------------    
    !- entropy-fix
    !-------------
    ! Change the value of the eigenvalues of |A| in the case of
    ! a sonic rarefaction
    !  We increase the eigenvalue (when it is needed, i.e. λ_l<0<λ_r)
    ! according to the following rule (Harten–Hyman Entropy Fix):
    !   The value of
    !      λ_m = (1-β)λ_r  +  βλ_l
    !   with β = (λ_r-λ_m)/(λ_r-λ_l) is changed to:
    !     λ_m~ = (1-β)λ_r -  βλ_l
    ! 
    ! vp1
    vp_temp = ValeurPropreSplit_x(u_l,c1,c2,ld)
    vp1_l   = vp_temp(1)
    vp_temp = ValeurPropreSplit_x(q_r_1(2,1)/q_r_1(1,1),c1,c2,ld)
    vp1_r   = vp_temp(1)
    ! vp2 (q_l_2=q_r_1)
    vp2_l   = vp_temp(2)
    vp_temp = ValeurPropreSplit_x(q_r_2(2,1)/q_r_2(1,1),c1,c2,ld)
    vp2_r   = vp_temp(2)
    ! vp3 (q_l_3=q_r_2)
    vp3_l   = vp_temp(3)
    vp_temp = ValeurPropreSplit_x(u_r,c1,c2,ld)
    vp3_r   = vp_temp(3)
    !- change the eigenvalues if needed for the entropy fix
    !- eigenvalues middle
    vp1_fix = vp_m(1); vp2_fix = vp_m(2); vp3_fix = vp_m(3)
    if (vp1_l*vp1_r<0) then
       vp1_fix = EntropyFixHartenHyman(vp1_l,vp_m(1),vp1_r)
       !-vp1_fix = max(-vp1_l,vp1_r) ! Méthode bourine
    endif
    if (vp2_l*vp2_r<0) then
       vp2_fix = EntropyFixHartenHyman(vp2_l,vp_m(2),vp2_r)
       !-vp2_fix = max(-vp2_l,vp2_r) ! Méthode bourine
    endif
    if (vp3_l*vp3_r<0) then
       vp3_fix = EntropyFixHartenHyman(vp3_l,vp_m(3),vp3_r)
       !-vp3_fix = max(-vp3_l,vp3_r) ! Méthode bourine
    endif
    ! the absolute value of D
    abs_D_roe = reshape( (/ &
         abs(vp1_fix)  , 0d0  , 0d0, &
         & 0d0  , abs(vp2_fix) , 0d0, &
         & 0d0  , 0d0  , abs(vp3_fix) /) , (/3,3/) )
    !-- the mix matrix A^+
    Abs_A_roe_bidon = matmul(transpose(P_inv),matmul(abs_D_roe,transpose(P)))       

  contains
    function EntropyFixHartenHyman(ld_l,ld_m,ld_r)
      implicit none
      Double precision, intent(in)     :: ld_l,ld_m,ld_r
      Double precision                 :: EntropyFixHartenHyman
      Double precision                 :: beta
      !- init
      beta = (ld_r - ld_m)/(ld_r-ld_l)
      !- the value
      EntropyFixHartenHyman = (1-beta)*ld_r - beta*ld_l
    end function EntropyFixHartenHyman

  end function Abs_A_roe_bidon


  
end module matrix_FVM
