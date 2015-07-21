module elementary

  !-- contains :
  ! Heaviside, Norm_2, AngleModulo, ThetaMean, AngleVec*
  ! init_rand_seed
  ! poly_abs
  
  implicit none

contains
  

  !------------------------------!
  !--   elementary function    --!
  !------------------------------!

  function Heaviside(s)
    !-- la fonction d'heaviside
    implicit none
    Double Precision, intent(in)  :: s
    Double Precision              :: Heaviside
    !-- the function
    If (s>0) Then
       Heaviside = 1d0;
    else
       Heaviside = 0d0;
    End If
  end Function Heaviside

  function Norm_2(V)
    !-- the square of the euclidien norm of a 2D vector
    implicit none
    Double Precision, dimension(2), intent(in)  :: V
    Double Precision                            :: Norm_2
    !-- the function
    Norm_2 = V(1)**2 + V(2)**2
  end Function Norm_2


  !------------------------------!
  !--           angle          --!
  !------------------------------!
  
  function AngleModulo(x)
    !- renvoie un angle compris entre -pi et pi
    implicit none
    Double Precision, intent(in)   :: x
    Double Precision, PARAMETER    :: PI = 3.14159265358979323846
    Double Precision               :: AngleModulo
    !-- the function
    AngleModulo = x - floor( (x+PI)/(2d0*PI) )*(2d0*PI)
  end function AngleModulo

  function ThetaMean(theta_l,theta_r)
    !- the mean angle
    implicit none
    Double Precision, intent(in)   :: theta_l, theta_r
    Double Precision               :: ThetaMean
    !-- the function
    ThetaMean = 5d-1*theta_l + 5d-1*(theta_l+AngleModulo(theta_r-theta_l))
    ThetaMean = AngleModulo(ThetaMean)
  end function ThetaMean

  function AngleVec(Vec)
    !- compute the angle of the vector Vec
    implicit none
    Double Precision, Dimension(2), intent(in) :: Vec
    Double Precision                           :: AngleVec
    Double Precision, PARAMETER                :: PI = 3.14159265358979323846
    !- the function
    !-- arctan
    if (Vec(1)==0) Then
       AngleVec = PI/2*sign(1d0,Vec(2));
    else
       AngleVec = atan( Vec(2)/Vec(1) );
    endif
    !-- on remet entre -pi et pi
    if (Vec(1)<0d0) Then
       AngleVec = AngleVec + PI*sign(1d0,Vec(2));
    endif
  endfunction AngleVec

  function AngleVec_N2(u,v)
    !- compute the angle of a matrix N*M  (u,v)
    implicit none
    Double Precision, Dimension(:,:), intent(in)  :: u,v
    Integer                                       :: i,j, m,n
    Double Precision, Dimension(:,:), Allocatable :: AngleVec_N2
    !- the function
    !- init
    m = size(u(:,1))
    n = size(u(1,:))
    Allocate(AngleVec_N2(m,n))
    !- loop
    Do i=1,m
       Do j=1,n
          AngleVec_N2(i,j) = AngleVec( (/ u(i,j), v(i,j) /) )
       end Do
    end Do
  endfunction AngleVec_N2


  !--------------------------------!
  !--          Random            --!
  !--------------------------------!

  Subroutine InitRandomSeed()
    !- initialise le random (??)
    implicit none
    Integer                             :: i, n, clock
    Integer, Dimension(:), Allocatable  :: seed

    Call Random_seed(size = n)
    Allocate(seed(n))
    Call System_clock(COUNT=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    Call Random_seed(PUT = seed)
    Deallocate(seed)
  End Subroutine InitRandomSeed

  
  !---------------------------------!
  !--   For the polynom scheme    --!
  !---------------------------------!

  Subroutine PolyAbs(A,a_m,a_p,a_max,degree,abs_A_poly)
    !- compute the approximate absolute value of A
    !   . A     : the matrix we have to find the absolute value
    !   . a_m   : the minimum eigenvalues (min(vp_1,...,vp_n))
    !   . a_p   : the maximum eigenvalues (max(vp_1,...,vp_n))
    !   . a_max : the maximum of absolute eigenvalues (max(|a_p|,|a_m|))
    !   . degree : the degree of the polynom we use to approach |x|
    !   . A_abs_poly : the approximation of |A|
    implicit none
    Double precision, Dimension(:,:), intent(in)  :: A
    Double precision, intent(in)                  :: a_m, a_p, a_max
    Integer, intent(in)                           :: degree
    Double precision, Dimension(:,:), intent(out) :: abs_A_poly
    !-- method poly
    Integer                                       :: i,N
    Double precision, Dimension(:,:), Allocatable :: Id
    Double precision                              :: alpha, beta, gamma
    Double precision                              :: det_M
    Double precision, Dimension(3,3)              :: M_inv
    Double precision, Dimension(3,1)              :: B, U

    !- the subroutine
    !- Init
    N = size(A(:,1))
    Allocate(Id(N,N))
    Id = 0d0
    Do i=1,N
       Id(i,i) = 1d0
    end Do
    
    !--- l'approximation de |A| ---!
    select case(degree)
    case(0)
       ! approximation d'ordre 0
       !   P(x) = abs(a_max)
       abs_A_poly = abs(a_max)*Id
    case(1)
       ! approximation d'ordre 1
       !   P(a_) = |a_|  ,  P(a+) = |a+|
       alpha = (abs(a_p) - abs(a_m))/(a_p - a_m)
       beta  = abs(a_p) - alpha*a_p
       abs_A_poly = alpha*A + beta*Id
    case(2)
       ! approximation d'ordre 2
       ! On résoud un système linéaire à 3 inconnues
       !   P(a_) = |a_|  ,  P(a+) = |a+|  ,  P'(a_max) = sign(a_max)
       !
       ! Avec Maple
       !   M = reshape( (/ a_p^2,a_p,1, a_m^2,a_m,1, 2*a_max,1,0/), (/3,3/))
       !   M := matrix( [[a_p^2,a_p,1], [a_m^2,a_m,1], [2*a_max,1,0]]);
       ! puis det(M) et adj(M)
       det_M = -a_p**2 + a_m**2 + 2*a_max*(a_p-a_m)

       M_inv = 1/det_M*reshape( (/ &
            -1d0, 1d0, a_p-a_m, &
            2*a_max, -2*a_max, -a_p**2+a_m**2, &
            a_m**2-2*a_m*a_max, -a_p**2+2*a_p*a_max, a_p**2*a_m-a_p*a_m**2 &
            /),(/3,3/))
       M_inv = transpose(M_inv)
       B = reshape( (/ &
            abs(a_p), &
            abs(a_m), &
            sign(1d0,a_max) &
              /),(/3,1/))
       U = matmul(M_inv,B)
       alpha = U(1,1); beta = U(2,1); gamma = U(3,1)
       ! l'approximation
       abs_A_poly = alpha*matmul(A,A) + beta*A + gamma*Id
    case default
       !- no better approximation...
       print *,"Pbm : in PolyAbs routine, the choice of the degree is 0, 1 or 2."
       print *,"      His value is",degree
       stop
    end select

    !- to finish
    Deallocate(Id)
    
  end Subroutine PolyAbs

  
  
end module elementary
