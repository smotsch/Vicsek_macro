module input_output_2D

  use Matrix                    ! to have ValeurPropre*

  TYPE PARAM_MacroVic
     Double Precision                   :: c1,c2,ld
     Double Precision                   :: Lx, Ly, dx, dy
     Double Precision                   :: Time, dt
     Integer                            :: boundCondX, boundCondY
     Integer                            :: degreePolyScheme
     Double Precision                   :: etaRelax
     Logical                            :: shouldSaveAll
     Logical                            :: isFormatVtk
     
     Integer                            :: nCellX, nCellY
  end TYPE PARAM_MacroVic

  TYPE PARAM_InitMV
     Integer                            :: choiceInit
     Double Precision                   :: meanRho, meanTheta, rangeRho, rangeTheta
     Double Precision                   :: rhoL, rhoR, thetaL, thetaR
     Double Precision                   :: rhoBG, thetaBG
     Double Precision                   :: rho1, theta1, x1Center, y1Center, radius1
     Double Precision                   :: rho2, theta2, x2Center, y2Center, radius2
     Double Precision                   :: rho3, theta3, x3Center, y3Center, radius3
     Integer                            :: shape1, shape2, shape3
     Character*80                       :: cheminRho, cheminTheta
  end TYPE PARAM_InitMV


contains

  Subroutine Lecture_2D(P,Pinit)
    !- take the initial  parameters
    implicit none

    !- Model
    TYPE(PARAM_MacroVic), intent(out)      :: P
    TYPE(PARAM_InitMV), intent(out)        :: Pinit
    !- temp
    character(8)                           :: temp
    integer                                :: i
    !- CFL
    Double Precision                       :: discri,cfl
    Double Precision, PARAMETER            :: PI = 3.14159265358979323846

    
    !-  Parameters de la simulation  !
    !--------------------------------!
    open(unit=15,file='PARAMETER_2D.txt',status='old')

    Do i=1,6
       read(15,*)temp
    End Do
    read(15,*)P%c1
    read(15,*)P%c2
    read(15,*)P%ld
    read(15,*)temp
    read(15,*)P%Lx
    read(15,*)P%Ly
    read(15,*)P%dx
    read(15,*)P%dy
    read(15,*)temp
    read(15,*)P%Time
    read(15,*)P%dt
    read(15,*)temp
    read(15,*)P%boundCondX
    read(15,*)P%boundCondY
    read(15,*)temp
    read(15,*)temp
    read(15,*)P%degreePolyScheme
    read(15,*)temp
    read(15,*)P%etaRelax
    read(15,*)temp
    read(15,*)P%shouldSaveAll
    read(15,*)P%isFormatVtk

    close(unit=15)
    
    P%nCellX = floor(P%Lx/P%dx+.5d0)
    P%nCellY = floor(P%Ly/P%dy+.5d0)

    !-- mise en garde...
    !-- condition CFL : dt max(|vp|)<dx, en particulier dt max(c,sqrt(ld))<dx
    !- different cfl for the splitting method
    discri = P%c1*P%ld-(P%c1-P%c2)*P%c2
    cfl = P%dt/P%dx*(P%c2 + sqrt(discri))

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
    read(16,*)Pinit%rhoR
    read(16,*)Pinit%thetaL
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
    read(16,*)Pinit%y1Center
    read(16,*)Pinit%radius1
    read(16,*)Pinit%shape1
    read(16,*)temp
    read(16,*)Pinit%rho2
    read(16,*)Pinit%theta2
    read(16,*)Pinit%x2Center
    read(16,*)Pinit%y2Center
    read(16,*)Pinit%radius2
    read(16,*)Pinit%shape2
    read(16,*)temp
    read(16,*)Pinit%rho3
    read(16,*)Pinit%theta3
    read(16,*)Pinit%x3Center
    read(16,*)Pinit%y3Center
    read(16,*)Pinit%radius3
    read(16,*)Pinit%shape3
    !- Predefined
    read(16,*)temp
    read(16,*)temp
    read(16,*)Pinit%cheminRho
    read(16,*)Pinit%cheminTheta

    close(unit=16)


    !--  Information on the terminal
    !-------------------------------
    print *,"******************************************************"
    print *,"***********  Parameters of the simulation  ***********"
    print *,"******************************************************"

    print *,"|--------  Parameters model  --------|"
    print "(A20,f7.2)","      c1        =  ",P%c1
    print "(A20,f7.2)","      c2        =  ",P%c2
    print "(A20,f7.2)","      lambda    =  ",P%ld
    !-- boundary condition
    select case (P%boundCondX)
    case(1)
       print *,"  Boundary in x :  periodic"
    case(2,3)
       print *,"  Boundary in x :  Neumann"
    case default
       print *,"  Boundary condition : uncorrect choice !!"
       stop
    end select
    select case (P%boundCondY)
    case(1)
       print *,"  Boundary in y :  periodic"
    case(2,3)
       print *,"  Boundary in y :  Neumann"
    case default
       print *,"  Boundary condition : uncorrect choice !!"
       stop
    end select

    print *,"|----  Parameters computation  ------|"
    !-- Numerical method
    if (P%etaRelax>0) then
       print "(A18,f7.4)","    eta         = ",P%etaRelax
    endif
    print "(A20,f7.1)","  Total time   =  ",P%time
    print "(A20,f7.4)","      dt       =  ",P%dt
    print "(A20,f7.4)","      dx       =  ",P%dx
    print "(A20,f7.4)","      dy       =  ",P%dy
    print "(A20,f7.1)","      Lx       =  ",P%Lx
    print "(A20,f7.1)","      Ly       =  ",P%Ly
    print "(A20,f7.2)","     CFL       =  ",CFL
    !-- shouldSaveAll
    if (P%shouldSaveAll) Then
       print *,"   Print        :   all data"
    else
       print *,"   Print        :   only first and final"
    endif
    !-- Octave or VTK
    if (P%isFormatVtk) Then
       print *,"   Format save  :   VTK"
    else
       print *,"   Format save  :   Octave/Matlab"
    endif
    print *,""

  end subroutine Lecture_2D



  Subroutine FilePrint_2D(rho,theta,numberIter)
    !- to write data
    implicit none
    Double Precision, Dimension(:,:), intent(in) :: rho, theta
    Integer, intent(in)                          :: numberIter
    Integer                                      :: nX,nY
    Character(9)                                 :: Extension
    !- the subroutine
99  format(I9.9)                  !--- format for the variable "Extension"
    nX = size(rho(:,1))
    nY = size(rho(1,:))
    write(Extension,99) numberIter
    !- open, write and close
    Open(Unit=10,file='../data/rho_'//trim(Extension)//trim('.dat'),form='UNFORMATTED')
    Open(Unit=11,file='../data/theta_'//trim(Extension)//trim('.dat'),form='UNFORMATTED')
    write(10) rho(2:(nX-1),2:(nY-1))
    write(11) theta(2:(nX-1),2:(nY-1))
    close(10)
    close(11)

  end Subroutine FilePrint_2D


  Subroutine FilePrintArray2D_vtk(rho,u,v,dx,dy,fileName,isUnFormatted,nbIteration)
    !- Same as FilePrint for array
    implicit none
    Double Precision, Dimension(:,:), intent(in) :: rho,u,v
    Double Precision, intent(in)                 :: dx,dy
    Character (len=*), intent(in)                :: fileName
    Logical, intent(in)                          :: isUnFormatted
    Integer, intent(in), optional                :: nbIteration
    real, Dimension(:,:), allocatable            :: uv_combined
    Integer                                      :: i,j, nCellX, nCellY
    Character(9)                                 :: cNbIteration
    Character(9)                                 :: cFormat
    Character (len=80)                           :: fileName_ext,cbuffer
    character(1), parameter                      :: newline=char(10)

    !- Init: change the fileName
    nCellX = size(rho(:,1))-2
    nCellY = size(rho(1,:))-2
    if (present(nbIteration)) Then
       write(cNbIteration,'(I9.9)') nbIteration
    else
       cNbIteration = ''
    endif
    if (isUnFormatted) then
       cFormat='_bin.vtk'
    else
       cFormat='.vtk'
    endif
    fileName_ext = trim(fileName)//trim(cNbIteration)//cFormat

    allocate(uv_combined(3,nCellX*nCellY))
    uv_combined(1,:) = pack(u(2:(nCellX+1),2:(nCellY+1)),.true.) !real( reshape( u  )
    uv_combined(2,:) = pack(v(2:(nCellX+1),2:(nCellY+1)),.true.)
    uv_combined(3,:) = 0

    !- Open and write
    if (isUnFormatted) then
       ! binary format
       open(unit       = 10,       &
         file       = trim(fileName_ext), &
         form       = 'UNFORMATTED', &
         access     = 'STREAM', &
         action     = 'WRITE',        &
         convert    = 'BIG_ENDIAN')
    else    
       Open(Unit=10,file=fileName_ext)
20     format(100000(I5))
40     format(100000(F6.2))
50     format(100000(F12.6))
    endif

    !-------------------------------------!
    !-----------  header file  -----------!
    !-------------------------------------!
    if (isUnFormatted) then
       write(10) "# vtk DataFile Version 1.0"//newline
       write(10) "vtk output"//newline
       write(10) "BINARY"//newline
       write(10) "DATASET RECTILINEAR_GRID"//newline
       write(cbuffer,fmt='(A,3I7)') 'DIMENSIONS ',nCellX,nCellY,1
       write(unit=10)trim(cbuffer)//newline
       write(cbuffer,fmt='(A,I7,A)') 'X_COORDINATES ',nCellX,' float'
       write(10) trim(cbuffer)//newline
       write(10) real(dx*(/ (i, i=1,nCellX) /) -dx/2)
       write(10) newline
       write(cbuffer,fmt='(A,I7,A)') 'Y_COORDINATES ',nCellY,' float'
       write(10) trim(cbuffer)//newline
       write(10) real(dy*(/ (j, j=1,nCellY) /) -dy/2)
       write(10) newline
       write(cbuffer,fmt='(A)') 'Z_COORDINATES 1 float'
       write(10) trim(cbuffer)//newline
       write(10) real(0)
    else
       write(10,"(A)") "# vtk DataFile Version 1.0"
       write(10,"(A)") "vtk output"
       write(10,"(A)") "ASCII"
       write(10,"(A)") "DATASET RECTILINEAR_GRID"
       write(10,"(A)",advance='no') "DIMENSIONS "
       write(10,20) nCellX, nCellY, 1
       write(10,"(A)",advance='no') "X_COORDINATES "
       write(10,'(I0)',advance='no') nCellX
       write(10,"(A)") " float"
       write(10,40) dx*(/ (i, i=1,nCellX) /) -dx/2
       write(10,"(A)",advance='no') "Y_COORDINATES "
       write(10,'(I0)',advance='no') nCellY
       write(10,"(A)") " float"
       write(10,40) dy*(/ (j, j=1,nCellY) /) -dy/2
       write(10,"(A)") "Z_COORDINATES 1 float"
       write(10,40) 0
    endif

    !-------------------------------------!
    !-----------      data     -----------!
    !-------------------------------------!

    if (isUnFormatted) then
       write(10) newline//newline
       write(cbuffer,fmt='(A,I14)')'POINT_DATA ', nCellX*nCellY
       write(10) trim(cbuffer)//newline
       write(10) "SCALARS Density float"//newline
       write(10) "LOOKUP_TABLE default"//newline
       write(10) real(rho(2:(nCellX+1),2:(nCellY+1)))
       write(10) newline//newline
       write(10) "VECTORS vectors float"//newline
       write(10) uv_combined
    else
       write(10,"(A)")
       write(10,"(A)",advance='no') "POINT_DATA "
       write(10,"(I0)") nCellX*nCellY
       write(10,"(A)") "SCALARS Density float"
       write(10,"(A)") "LOOKUP_TABLE default"
       write(10,50) rho(2:(nCellX+1),2:(nCellY+1))
       write(10,"(A)")
       write(10,"(A)") "VECTORS vectors float"
       write(10,50) uv_combined
    endif

    close(10)

    Deallocate(uv_combined)
    
  end Subroutine FilePrintArray2D_vtk

  
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



end module input_output_2D
