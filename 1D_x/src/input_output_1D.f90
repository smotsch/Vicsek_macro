module input_output_1D

  use type_PARAM                !- for PARAM_MacroVic, PARAM_InitMV
  use matrix_FVM                !- to have ValeurPropre_1/2

contains

  Subroutine Lecture_1D(P,Pinit)
    !- take the initial  parameters
    implicit none
    
    !- Model
    TYPE(PARAM_MacroVic), intent(out)       :: P
    TYPE(PARAM_InitMV), intent(out)         :: Pinit
    !- temp
    character(8)                            :: temp
    integer                                 :: i
    !- CFL
    Double Precision                        :: tanThetaPw2, theta1, theta2
    Double Precision                        :: vp1_x,vp2_x, discri, cfl
    Double Precision, PARAMETER             :: PI = 3.14159265358979323846
    

    !-  Parameters de la simulation  !
    !--------------------------------!
    open(unit=15,file='PARAMETER_1D.txt',status='old')

    Do i=1,6
       read(15,*)temp
    End Do
    read(15,*)P%c1
    read(15,*)P%c2
    read(15,*)P%ld
    read(15,*)temp
    read(15,*)P%Lx
    read(15,*)P%dx
    read(15,*)temp
    read(15,*)P%Time
    read(15,*)P%dt
    read(15,*)temp
    read(15,*)P%boundCond
    read(15,*)temp
    read(15,*)temp
    read(15,*)P%methodNum
    read(15,*)temp
    read(15,*)P%degree
    read(15,*)P%epsilon
    read(15,*)temp
    read(15,*)P%shouldSaveAll
    
    close(unit=15)

    P%nSpace = floor(P%Lx/P%dx+.5d0) + 1

    !-- mise en garde...
    !-- condition CFL : P%dt max(|vp|)<dx
    if (P%methodNum/=1) then
       tanThetaPw2 = 1/(4*P%c1*P%ld)*( ((P%c1-P%c2)**2-4*P%c1*P%ld)**2/(P%c1+P%c2)**2 - (P%c1-P%c2)**2 )
       theta2 = atan( sqrt(tanThetaPw2) )
       theta1 = PI-theta2
       call ValeurPropre_x(theta1,P%c1,P%c2,P%ld, vp1_x,vp2_x)
       cfl = P%dt/P%dx*max( abs(vp1_x), abs(vp2_x))
    else
       !- different cfl for the splitting method
       discri = P%c1*P%ld-(P%c1-P%c2)*P%c2
       cfl = P%dt/P%dx*(P%c2 + sqrt(discri))
    endif


    !-  Parameters pour la condition initiale  !
    !------------------------------------------!
    open(unit=16,file='PARAMETER_init.txt',status='old')

    read(16,*)temp
    read(16,*)temp
    read(16,*)Pinit%choiceInit
    read(16,*)temp
    !- 1) Random
    read(16,*)temp
    read(16,*)temp
    read(16,*)Pinit%meanRho
    read(16,*)Pinit%meanTheta
    read(16,*)Pinit%rangeRho
    read(16,*)Pinit%rangeTheta
    !- 2) Riemann pbm
    read(16,*)temp
    read(16,*)temp
    read(16,*)Pinit%rhoL
    read(16,*)Pinit%thetaL
    read(16,*)Pinit%rhoR
    read(16,*)Pinit%thetaR
    !- 3) Blocks
    read(16,*)temp
    read(16,*)temp
    read(16,*)temp
    read(16,*)Pinit%rhoBG
    read(16,*)Pinit%thetaBG
    read(16,*)temp
    read(16,*)Pinit%rho1
    read(16,*)Pinit%theta1
    read(16,*)Pinit%x1Center
    read(16,*)Pinit%radius1
    read(16,*)temp
    read(16,*)Pinit%rho2
    read(16,*)Pinit%theta2
    read(16,*)Pinit%x2Center
    read(16,*)Pinit%radius2
    read(16,*)temp
    read(16,*)Pinit%rho3
    read(16,*)Pinit%theta3
    read(16,*)Pinit%x3Center
    read(16,*)Pinit%radius3
    !- Predefined
    read(16,*)temp
    read(16,*)temp
    read(16,*)Pinit%cheminRho
    read(16,*)Pinit%cheminTheta
    
    close(unit=16)



    !-----  Information on the terminal  -----!
    !-----------------------------------------!
    print *,"******************************************************"
    print *,"***********  Parameters of the simulation  ***********"
    print *,"******************************************************"

    print *,"|--------  Parameters model  --------|"
    print "(A25,f8.2)","      c1        =  ",P%c1
    print "(A25,f8.2)","      c2        =  ",P%c2
    print "(A25,f8.2)","      lambda    =  ",P%ld
    !-- boundary condition
    select case (P%boundCond)
    case(1)
       print *,"      Boundary cond  : periodic"
    case(2)
       print *,"      Boundary cond  : Neumann homogenous"
    case default
       print *,"  Boundary condition : uncorrect choice !!"
       stop
    end select
    
    print *,"|----  Parameters computation  ------|"
    !-- Numerical method
    select case (P%methodNum)
    case(1)
       print *,"   Num scheme  : splitting"
       print "(A17,f10.5)","    epsilon     = ",P%epsilon
    case(2)
       print *,"   Num scheme  : conservative"
    case(3)
       print *,"   Num scheme  : semi-conservative"
    case(4)
       print *,"   Num scheme  : upwind scheme"
    end select
    print "(A17,f7.2)", "   Total time   = ",P%time
    print "(A17,f10.5)","        dt      = ",P%dt
    print "(A17,f10.5)","        dx      = ",P%dx
    print "(A17,f7.2)", "        Lx      = ",P%Lx
    print "(A17,f7.2)", "       CFL      = ",cfl
    !-- shouldSaveAll
    if (P%shouldSaveAll) Then
       print *,"    Print      : all the data"
    else
       print *,"    Print      : only first and final step"
    end if
    print *,""
    
    
  end subroutine Lecture_1D


  
  !###################  print  ###################!
  !###############################################!

  Subroutine FilePrint_1D(rho,theta,numberIter)
    !- to write data
    implicit none
    Double Precision, Dimension(:), intent(in) :: rho, theta
    Integer, intent(in)                        :: numberIter
    Character(9)                               :: Extension
    !- the subroutine
99  format(I9.9)                  !--- format for the variable "Extension"
    !- test
    write(Extension,99) numberIter
    !- open, write and close
    Open(Unit=10,file='../data/rho_'//trim(Extension)//trim('.dat'),form='UNFORMATTED')
    Open(Unit=11,file='../data/theta_'//trim(Extension)//trim('.dat'),form='UNFORMATTED')
    write(10) rho
    write(11) theta
    close(10)
    close(11)

  end Subroutine FilePrint_1D
 
  
  Subroutine FilePrint_1D_relax(rho,u,v,numberIter)
    !- to write data for the relax splitting method
    implicit none
    Double Precision, Dimension(:), intent(in) :: rho, u,v
    Integer, intent(in)                        :: numberIter
    Character(9)                               :: Extension
    !- the subroutine
99  format(I9.9)                  !--- format for the variable "Extension"
    !- test
    write(Extension,99) numberIter
    !- open, write and close
    Open(Unit=10,file='../data/rho_'//trim(Extension)//trim('.dat'),form='UNFORMATTED')
    Open(Unit=11,file='../data/u_'//trim(Extension)//trim('.dat'),form='UNFORMATTED')
    Open(Unit=12,file='../data/v_'//trim(Extension)//trim('.dat'),form='UNFORMATTED')
    write(10) rho
    write(11) u
    write(12) v
    close(10)
    close(11)
    close(12)

  end Subroutine FilePrint_1D_relax

  

  
  !###################  progress bar  ###################!
  !######################################################!


  Subroutine BarProgress(i,imax)
    implicit none
    integer          :: i,imax,k
    character(len=1) :: bar, back
    ! the subroutine
    !-init
    back = char(8)
    bar  = '='
    ! print the percentage and the bar
    write(*,'(256a1)', advance='no') (back, k =1,(30*i/imax)+9)
    write(*,'(2x,1i3,1a1,2x,1a1,256a1)', advance='no') &
         100*i/imax,'%','|', (bar, k =1,30*i/imax)
    if (i==imax) Then
       write(*,'(a)') '| done.\n'
    end if
  End subroutine BarProgress




end module input_output_1D
