module initial_condition_2D

  use elementary                !- to have random
  use input_output_2D           !- for PARAM_Vic, PARAM_InitMV
  
contains

  subroutine InitCond_2D(rho,u,v,P,Pinit)
    !- prescribe the initial condition for the 2D model
    implicit none
    Double Precision, Dimension(:,:), intent(inout) :: rho,u,v
    TYPE(PARAM_MacroVic), intent(in)                :: P
    TYPE(PARAM_InitMV), intent(in)                  :: Pinit
    Double Precision, Dimension(:,:), Allocatable   :: theta
    Double Precision, PARAMETER                     :: PI = 3.14159265358979323846
    Integer                                         :: i,j
    Integer                                         :: nX,nY
    ! mill formation
    Double Precision, dimension(2)                  :: millCenter
    Double Precision                                :: minRho, maxRho, maxR, r, xx, yy

    ! init
    nX = P%nCellX+2
    nY = P%nCellY+2
    Allocate( theta(nX,nY) )
    
    !- the subroutine
    select case (Pinit%choiceInit)
    case(1)
       !---------------------------------------!
       !------------  1) random  --------------!
       !---------------------------------------!
       Call InitRandomSeed()
       Call Random_number(rho)
       Call Random_number(theta)
       !- init
       rho   = Pinit%meanRho   + Pinit%rangeRho*(rho - .5d0)
       theta = Pinit%meanTheta + Pinit%rangeTheta*(theta - .5d0)

    case(2)
       !--------------------------------------------!
       !----------   2) Riemann problem   ----------!
       !--------------------------------------------!
       rho(   1:(nX/2)    ,:)  = Pinit%rhoL
       rho(   (nX/2+1):nX ,:)  = Pinit%rhoR
       theta( 1:(nX/2)    ,:)  = Pinit%thetaL
       theta( (nX/2+1):nX ,:)  = Pinit%thetaR

    case(31,32,33)
       !------------------------------------------!
       !--------        3) Blocks        ---------!
       !------------------------------------------!
       rho   = Pinit%rhoBG
       theta = Pinit%thetaBG
       
       !--- The first block ---!
       Call AddBlockInit(rho,theta,Pinit%rho1,Pinit%theta1,&
            Pinit%x1Center,Pinit%y1Center,Pinit%radius1,Pinit%shape1,P%dx,P%dy)

       if (Pinit%choiceInit/=31) then
          ! The second block
          Call AddBlockInit(rho,theta,Pinit%rho2,Pinit%theta2,&
               Pinit%x2Center,Pinit%y2Center,Pinit%radius2,Pinit%shape2,P%dx,P%dy)
       endif
       If ( Pinit%choiceInit==33 ) Then
          ! The second block
          Call AddBlockInit(rho,theta,Pinit%rho3,Pinit%theta3,&
               Pinit%x3Center,Pinit%y3Center,Pinit%radius3,Pinit%shape3,P%dx,P%dy)
       End If
              
    case(4)
       !------------------------------------------!
       !--------------   4) Mill   ---------------!
       !------------------------------------------!
       !- center
       millCenter = (/ P%Lx/2d0 , P%Ly/2d0 /);
       minRho = 1d-9                  ! pour r=0
       maxR   = 5d-1*sqrt(P%Lx**2+P%Ly**2);
       maxRho = minRho + maxR**(P%c2/P%ld)         ! pour r=maxR , alpha = c/ld
       do i = 1,nX
          do j = 1,nX
             xx = (i-1)*P%dx-millCenter(1)
             yy = (j-1)*P%dy-millCenter(2)
             r = sqrt( xx**2 + yy**2 );
             rho(i,j) = minRho + (maxRho-minRho)*(r/maxR)**(P%c2/P%ld)
             ! ou encore par exemple 2-1.8*(maxR-r)/maxR;
             if (xx > 0) Then
                theta(i,j) = atan(yy/xx) + PI/2d0;
             elseif (xx<0) Then
                theta(i,j) = -atan(-yy/xx) - PI/2d0;
             else
                if (yy>0) Then
                   theta(i,j) = PI
                else
                   theta(i,j) = 0
                end if
             end if
          end do
       end do

    case(5)
       !-------------------------------------------!
       !-------  5) Read from extern files  -------!
       !-------------------------------------------!
       !-
       !- Warning! The initial condition should
       !-            have the right format (dx,dy,Lx,Ly)
       Open(Unit=10,file=Pinit%cheminRho) !kop: marche pas...
       Open(Unit=11,file=Pinit%cheminTheta)
       Do j=1,nY
          read(10,*) rho(j,:)
          read(11,*) theta(j,:)
       end Do
       close(10)
       close(11)

    case default
       !----------------------------------------------!
       !-------------   Bad choice...   --------------!
       !----------------------------------------------!
       print *,"Wrong initial condition...\n"
       stop                     ! we stop the program

    end select

    ! get back with u,v
    u = cos(theta);
    v = sin(theta);

    deallocate(theta)
    
  end subroutine InitCond_2D


  Subroutine AddBlockInit(rho,theta,valueRho,valueTheta,xCenter,yCenter,radius,shape,dx,dy)
    !- Rajoute un "block" de densité valueRho avec direction valueTheta à la donnée rho, theta
    implicit none
    Double Precision, Dimension(:,:), intent(inout) :: rho, theta
    Double Precision                                :: valueRho,valueTheta
    Double Precision                                :: xCenter,yCenter,radius,dx,dy
    Integer                                         :: shape
    Integer                                         :: nXcenter,nYcenter,nXradius,nYradius
    Double Precision                                :: pythaY
    Integer                                         :: i,j,nPythaY
    
    !--- The first block ---!
    nXcenter      = floor(xCenter/dx)
    nYcenter      = floor(yCenter/dy)
    nXradius      = floor(radius/dx)
    nYradius      = floor(radius/dy)
    Select Case (shape)
    case(0)
       ! circle 0
       Do i=-nXradius,nXradius
          pythaY = sqrt(radius**2-(i*dx)**2)
          nPythaY = floor(pythaY/dy)
          Do j=-nPythaY,nPythaY
             rho(   nXcenter+i , nYcenter+j ) = valueRho
             theta( nXcenter+i , nYcenter+j ) = valueTheta
          End Do
       End Do
    case(1)
       ! rectangle 1
       rho( (nXcenter-nXradius):(nXcenter+nXradius), &
            (nYcenter-nYradius):(nYcenter+nYradius) ) = valueRho
       theta( (nXcenter-nXradius):(nXcenter+nXradius), &
            (nYcenter-nYradius):(nYcenter+nYradius) ) = valueTheta
    case(2)
       ! une bande 2
       rho( (nXcenter-nXradius):(nXcenter+nXradius)   ,:) = valueRho
       theta( (nXcenter-nXradius):(nXcenter+nXradius) ,:) = valueTheta
    End Select
    
  End Subroutine AddBlockInit

  
  
end module initial_condition_2D
