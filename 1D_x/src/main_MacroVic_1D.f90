program MacroVic_1D

  use elementary
  use flux_FVM
  use input_output_1D
  use boundary_1D
  use initial_condition_1D
  
  
  !------------- Declare variables and constants ---------------!
  implicit none
  
  !- Model
  TYPE(PARAM_MacroVic)               :: P
  TYPE(PARAM_InitMV)                 :: Pinit
  !- Data
  Double Precision, Dimension(:), Allocatable :: rho, theta, thetaCons
  Double Precision, Dimension(:), Allocatable :: u, v
  Double Precision, Dimension(:), Allocatable :: rhoUpwind, thetaUpwind
  !- Flux
  Double Precision, Dimension(:), Allocatable :: FluxSemiRho, FluxSplitRho
  Double Precision, Dimension(:), Allocatable :: FluxConsRho, FluxConsTheta
  Double Precision, Dimension(:), Allocatable :: FluxSplitU, FluxSplitV
  !- Temp
  Integer                            :: nSpace, i, iTime, nTime
  Double Precision, Dimension(2)     :: Utemp2, Ftemp2, Omega
  Double Precision, Dimension(3)     :: Utemp3, Ftemp3
  Double Precision                   :: thetaAbs, dtDivDx, C_0, f1temp
  !- time to execut the code
  real                               :: start, finish
  Double Precision, PARAMETER        :: PI = 3.14159265358979323846


  
  !- let's count the time it takes -!
  Call Cpu_time(start)
  
  !-------------------  Read parameters  ---------------------!
  Call Lecture_1D(P,Pinit)
  
  !------------------- initialisation ------------------------!
  nSpace  = P%nSpace
  dtDivDx = P%dt/P%dx
  nTime   = floor(P%Time/P%dt + 51d-2)
  if (nTime>=1d6) Then
     !- warning
     print *,"Warning : the number of iteration in time is too big."
     stop
  end if
  !-------- Allocate --------!
  !- Allocate mass and speed
  Allocate(rho(nSpace),theta(nSpace))
  select case (P%methodNum)
  case(1)
     Allocate(u(nSpace),v(nSpace))
  case(2)
     Allocate(thetaCons(nSpace))
  end select
  !- Allocate for upwind method
  if (P%methodNum/=1) Then
     Allocate(rhoUpwind(nSpace),thetaUpwind(nSpace))
  endif
  !- Allocate flux
  select case (P%methodNum)
     ! explanation :   !   Flux*(i)  |   rho(i)   |   Flux*(i+1)
  case(1)
     !- flux for the splitting method
     Allocate(FluxSplitRho(nSpace+1))
     Allocate(FluxSplitU(nSpace+1),FluxSplitV(nSpace+1))
  case(2)
     !- conservative
     Allocate(FluxConsRho(nSpace+1),FluxConsTheta(nSpace+1))
  case(3)
     !- semi-conservative
     Allocate(FluxSemiRho((nSpace+1)))
  end select
  !-------- initial condition --------!
  Call InitCond_1D(rho,theta,P,Pinit)

  !- apply boundary conditions
  Call BoundaryCondition_1D(rho,P%boundCond)
  Call BoundaryCondition_1D(theta,P%boundCond)
  !--- write initial condition ---!
  Call FilePrint_1D(rho,theta,0) ! au temps 0
  
  !- for others methods
  select case (P%methodNum)
  case(1)
     u = cos(theta)
     v = sin(theta)
     if (P%epsilon/=0d0) Then
        Call FilePrint_1D_relax(rho,u,v,0) ! au temps 0
     end if
  case(2)
     thetaCons = theta   ! Dos theta at the same time
  end select

  
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
  !------------------------- The big loop ------------------------!
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
  
  Do iTime=1,nTime

     !--      1) Flux (at the interface i+1/2             --!
     !------------------------------------------------------!
     select case (P%methodNum)
     case(1)
        !- FluxSplit* (splitting method)
        Do i = 2,nSpace
           !- at i-1/2
           Ftemp3 = FluxSplit_x( &
                rho(i-1), rho(i), u(i-1), u(i), v(i-1), v(i), &
                P%c1,P%c2,P%ld,P%degree)
           FluxSplitRho(i) = Ftemp3(1)
           FluxSplitU(i)   = Ftemp3(2)
           FluxSplitV(i)   = Ftemp3(3)
        end Do
     case(2)
        !- FluxCons* (conservative method)
        Do i = 2,nSpace
           !- at i-1/2
           Ftemp2 = FluxCons_x( &
                rho(i-1), rho(i), thetaCons(i-1), thetaCons(i), &
                P%c1,P%c2,P%ld,P%degree)
           FluxConsRho(i)   = Ftemp2(1)
           FluxConsTheta(i) = Ftemp2(2)
        end Do
     case(3)
        !- FluxSemiRho (semi-conservative method)
        Do i = 2,nSpace
           !- at i-1/2
           FluxSemiRho(i) = FluxSemi_x(rho(i-1),rho(i),theta(i-1),theta(i),P%c1,P%c2,P%ld)
        end Do        
     end select

     
     !--    1bis) Upwind     --!
     !-------------------------!
     If (P%methodNum/=1) Then
        Do i = 2,(nSpace-1)
           !- upwind
           Utemp2 = Upwind_x(rho(i-1),rho(i),rho(i+1), &
                theta(i-1),theta(i),theta(i+1), &
                P%c1,P%c2,P%ld,P%dt,P%dx)
           rhoUpwind(i)   = Utemp2(1)
           thetaUpwind(i) = Utemp2(2)
        end Do
     End If
     

     !--       2) New value (update in time)               -!
     !------------------------------------------------------!
     Do i = 2,(nSpace-1)

        select case (P%methodNum)
        case(1)
           !-----  splitting  -----!
           !- conservative variable (masse et quantité de mvt)
           Utemp3 = (/ rho(i) , rho(i)*u(i), rho(i)*v(i) /)
           !- Roe method (U^{n+1/2})
           Utemp3 = Utemp3 - dtDivDx*(/ &
                FluxSplitRho(i+1) - FluxSplitRho(i), &
                FluxSplitU(i+1)   - FluxSplitU(i)  , &
                FluxSplitV(i+1)   - FluxSplitV(i)   /)
           !- relaxation term (Omega^{n+1})
           Omega = (/ Utemp3(2)/Utemp3(1) , Utemp3(3)/Utemp3(1) /)
           if (P%epsilon==0d0) Then
              !- Omega of unit length
              Omega = Omega/Norm2(Omega)
           else
              !- Omega relax
              C_0 = 1/Norm2(Omega)**2 - 1d0
              Omega =  (1d0/sqrt(1d0+C_0*exp(-2d0/P%epsilon/P%dt)))/Norm2(Omega)*Omega
           endif
           !- update
           rho(i) = Utemp3(1)
           u(i)   = Omega(1)
           v(i)   = Omega(2)
        case(2)
           !-----  conservative  -----!
           ! rho
           rho(i) = rho(i) - dtDivDx*(FluxConsRho(i+1) - FluxConsRho(i));
           ! theta
           f1temp = f1(thetaCons(i)) - dtDivDx*(FluxConsTheta(i+1) - FluxConsTheta(i));
           thetaAbs = f1_inv(f1temp);
           !- au final (idée Pierre pour détermination de l'angle)
           theta(i)     = thetaUpwind(i) 
           thetaCons(i) = thetaAbs*sign(1d0,AngleModulo(theta(i)))
        case(3)
           !-----  semi-conservative  -----!
           rho(i)   = rho(i) - dtDivDx*(FluxSemiRho(i+1) - FluxSemiRho(i))
           !- upwind pour theta...
           theta(i) = thetaUpwind(i)
        case(4)
           !-----  upwind  -----!
           rho(i)   = rhoUpwind(i)
           theta(i) = thetaUpwind(i)

        end select
        
     end Do
     
     !--          3) Boundary condition                              --!
     !-----------------------------------------------------------------!

     Call BoundaryCondition_1D(rho,P%boundCond)
     select case (P%methodNum)
     case(1)
        Call BoundaryCondition_1D(u,P%boundCond)
        Call BoundaryCondition_1D(v,P%boundCond)
     case(2)
        Call BoundaryCondition_1D(theta,P%boundCond)
        Call BoundaryCondition_1D(thetaCons,P%boundCond)
     case(3,4)
        Call BoundaryCondition_1D(theta,P%boundCond)
     end select

     
     !--         4) Output (print)                                   --!
     !-----------------------------------------------------------------!

     if (P%shouldSaveAll .or. iTime==nTime) Then
        !-- pour shouldSaveAll=false, on garde que le 1er et dernier plot
        select case (P%methodNum)
        case(1)
           Call FilePrint_1D(rho,AngleVec_N(u,v),iTime)
           if (P%epsilon /= 0d0) Then
              Call FilePrint_1D_relax(rho,u,v,iTime)
           end if
        case(2)
           Call FilePrint_1D(rho,thetaCons,iTime)
        case(3,4)
           Call FilePrint_1D(rho,theta,iTime)
        end select
     end if
     
     !- progress...
     call BarProgress(iTime,nTime)
     
  end Do

  !-----------------------------------------------------------------!
  !------------------------ end loop in time -----------------------!
  !-----------------------------------------------------------------!

  
  !- finishing line...
  Call Cpu_time(finish)
  print "(A24,f11.3)"," Time to compute (s)   = ",finish-start
  print *,"******************************************************"
  print *,""

  !- Deallocate
  deallocate(rho)
  if (P%methodNum/=1) Then
     deallocate(theta,rhoUpwind,thetaUpwind)
  end if
  select case (P%methodNum)
  case(1)
     deallocate(u,v,FluxSplitRho,FluxSplitU,FluxSplitV)
  case(2)
     deallocate(thetaCons,FluxConsRho,FluxConsTheta)
  case(3)
     deallocate(FluxSemiRho)
  end select
  
end program MacroVic_1D
