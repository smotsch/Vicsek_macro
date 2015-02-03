module initial_condition_1D

  use type_PARAM                !- for PARAM_MacroVic, PARAM_InitMV
  use elementary                !- to have random
  
contains

  subroutine InitCond_1D(rho,theta,P,Pinit)
    !- give the initial condition for the 1D problem
    implicit none
    Double Precision, Dimension(:), intent(inout) :: rho, theta
    TYPE(PARAM_MacroVic), intent(in)              :: P
    TYPE(PARAM_InitMV), intent(in)                :: Pinit
    Integer                                       :: nX1, nRadius1, nX2, nRadius2, nX3, nRadius3

    
    select case (Pinit%choiceInit)
    case(1)
       !---------------------------------------!
       !------------  1) random  --------------!
       !---------------------------------------!
       Call InitRandomSeed()
       Call Random_number(rho)
       Call Random_number(theta)
       !- init
       rho   = Pinit%meanRho + Pinit%rangeRho*(rho - .5d0)
       theta = Pinit%meanTheta + Pinit%rangeTheta*(theta - .5d0)
       
    case(2)
       !--------------------------------------------!
       !----------   2) Riemann problem   ----------!
       !--------------------------------------------!
       rho(   1:(P%nSpace/2)            )  = Pinit%rhoL
       rho(   ((P%nSpace/2)+1):P%nSpace )  = Pinit%rhoR
       theta( 1:(P%nSpace/2)            )  = Pinit%thetaL
       theta( ((P%nSpace/2)+1):P%nSpace )  = Pinit%thetaR

    case(31,32,33)
       !------------------------------------------!
       !--------        3) Blocks        ---------!
       !------------------------------------------!
       rho   = Pinit%rhoBG
       theta = Pinit%thetaBG
       
       ! The first block
       nX1      = floor(Pinit%x1Center/P%dx)
       nRadius1 = floor(Pinit%radius1/P%dx)
       rho(   (nX1-nRadius1):(nX1+nRadius1) )   = Pinit%rho1
       theta( (nX1-nRadius1):(nX1+nRadius1) )   = Pinit%theta1

       if (Pinit%choiceInit/=31) then
          ! The second block
          nX2      = floor(Pinit%x2Center/P%dx)
          nRadius2 = floor(Pinit%radius2/P%dx)
          rho(   (nX2-nRadius2):(nX2+nRadius2) )   = Pinit%rho2
          theta( (nX2-nRadius2):(nX2+nRadius2) )   = Pinit%theta2
       endif

       if (Pinit%choiceInit==33) then
          ! The third block
          nX3      = floor(Pinit%x3Center/P%dx)
          nRadius3 = floor(Pinit%radius3/P%dx)
          rho(   (nX3-nRadius3):(nX3+nRadius3) )   = Pinit%rho3
          theta( (nX3-nRadius3):(nX3+nRadius3) )   = Pinit%theta3
       endif

    case(4)
       !-------------------------------------------!
       !-------  4) Read from extern files  -------!
       !-------------------------------------------!
       !-
       !- Warning! The initial condition should
       !-            have the right format (dx,dy,Lx,Ly)
       Open(Unit=10,file=Pinit%cheminRho) !kop: marche pas...
       Open(Unit=11,file=Pinit%cheminTheta)
       read(10,*) rho
       read(11,*) theta
       close(10)
       close(11)

    case default
       !----------------------------------------------!
       !-------------   Bad choice...   --------------!
       !----------------------------------------------!
       print *,"Wrong initial condition...\n"
       stop                     ! we stop the program
       
    end select


       !if (I_init_cond==901) Then
       !   !- ghost shock with smooth transition
       !   !- get ep
       !   write(*,"(a)",advance='no') "  Espacement (epsilon) ? "
       !   read *,ep
       !   !- modify initial condition
       !   n_ep = floor(ep/P%dx + 5d-1) + 1
       !   Do i = -n_ep,n_ep
       !      theta( nSpace/2 + i ) = -i*1d0/n_ep + 1d-9
       !   end Do
       !endif

    
  end subroutine InitCond_1D


end module initial_condition_1D
