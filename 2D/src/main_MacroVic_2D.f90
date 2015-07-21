program MacroVic_2D

  use elementary
  use flux
  use input_output_2D
  use boundary_2D
  use initial_condition_2D


  !------------- Declare variables and constants ---------------!
  implicit none

  !- Parameters model
  TYPE(PARAM_MacroVic)               :: P
  TYPE(PARAM_InitMV)                 :: Pinit
  !- Data
  Double Precision, Dimension(:,:), Allocatable :: rho,u,v
  !- Flux
  Double Precision, Dimension(:,:), Allocatable :: FluxSplitRhoX
  Double Precision, Dimension(:,:), Allocatable :: FluxSplitUx, FluxSplitVx
  Double Precision, Dimension(:,:), Allocatable :: FluxSplitRhoY
  Double Precision, Dimension(:,:), Allocatable :: FluxSplitUy, FluxSplitVy
  !- Temp
  Integer                            :: nCellX, nCellY, i,j, iTime, nTime
  Double Precision, Dimension(2)     :: Omega
  Double Precision, Dimension(3)     :: Ftemp3, Utemp3
  Double Precision                   :: dtDivDx,dtDivDy, C_0
  !- time to execut the code
  real                               :: start, finish
  Double Precision, PARAMETER        :: PI = 3.14159265358979323846



  !- let's count the time it takes -!
  Call Cpu_time(start)

  !---------------- Lecture des paramètres -------------------!
  Call Lecture_2D(P,Pinit)

  !------------------- initialisation ------------------------!
  nCellX = P%nCellX; dtDivDx = P%dt/P%dx
  nCellY = P%nCellY; dtDivDy = P%dt/P%dy
  nTime = floor(P%Time/P%dt + 51d-2)
  if (nTime>=1d6) Then
     !- warning
     print *,"Warning : the number of iteration in time is too large."
     stop
  end if
  !- Allocate ρ,u,v: size (nCellX,nCellY) + 2 extra cells for the BC
  Allocate( rho(nCellX+2,nCellY+2))
  Allocate( u(nCellX+2,nCellY+2), v(nCellX+2,nCellY+2) )
  !- Allocate flux
  Allocate( FluxSplitRhoX(nCellX+1,nCellY+1), FluxSplitRhoY(nCellX+1,nCellY+1) )
  Allocate( FluxSplitUx(nCellX+1,nCellY+1),   FluxSplitUy(nCellX+1,nCellY+1) )
  Allocate( FluxSplitVx(nCellX+1,nCellY+1),   FluxSplitVy(nCellX+1,nCellY+1) )

  !- explanation: ρ,u,v
  !--------------
  !
  !         ⋮     extra     ⋮       ...     ⋮       extra        ⋮            
  !         -------------------------------------------------------
  !  ⋮ extra |  (2,CellY+1) |      ...     | (CellX+1,CellY+1) | extra ⋮
  !  ⋮ extra |  (2,CellY)   |      ...     | (CellX+1,CellY)   | extra ⋮
  !    
  !          ⋮              ⋮               ⋮                    ⋮
  !    
  !  ⋮ extra |    (2,3)     |      ...     |   (CellX+1,3)     | extra ⋮
  !  ⋮ extra |    (2,2)     |      ...     |   (CellX+1,2)     | extra ⋮
  !         -----------------------------------------------------
  !         ⋮     extra     ⋮       ...     ⋮    extra           ⋮
  ! 
  ! Flux:
  !------
  !  ⋮ extra -F(1)-  Cell 2  -F(2)-  Cell 3 -  ...  -  Cell+1  -F(Cell+1)-  extra ⋮
  
  !- initial condition
  Call InitCond_2D(rho,u,v,P,Pinit)
  !--- apply the boundary conditions ---!
  Call BoundaryCondition_2D(rho,P%boundCondX,P%boundCondY)
  if (P%boundCondX==3 .and. P%boundCondY==3) then
     Call BC_Vortex_2D(u,v)
  endif
  Call BoundaryCondition_2D(u,P%boundCondX,P%boundCondY)
  Call BoundaryCondition_2D(v,P%boundCondX,P%boundCondY)
  !--- writing ---!
  if (P%isFormatVtk) then
     Call FilePrintArray2D_vtk(rho,u,v,P%dx,P%dy,"../data/rho_uv_2D_",.true.,0)
  else
     Call FilePrint_2D(rho,AngleVec_N2(u,v),0)
  endif

  
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
  !------------------------- The big loop ------------------------!
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  Do iTime=1,nTime
     
     !----------------------------------------------------!
     !---------------------   A) xXx   -------------------!
     !----------------------------------------------------!

     !--    A.1) Flux (calcul des valeurs aux interfaces)   --!
     !--------------------------------------------------------!
     Do i = 1,nCellX+1
        Do j = 2,nCellY+1
           !- en i-1/2
           Ftemp3 = FluxSplit_x( &
                rho(i,j),rho(i+1,j),u(i,j),u(i+1,j),v(i,j),v(i+1,j), &
                P%c1,P%c2,P%ld,P%degreePolyScheme)
           FluxSplitRhoX(i,j) = Ftemp3(1)
           FluxSplitUx(i,j)   = Ftemp3(2)
           FluxSplitVx(i,j)   = Ftemp3(3)
        End do
     End do
     !--     A.2) New value (on avance d'un 1/2 pas de temps)      --!
     !---------------------------------------------------------------!
     Do i = 2,(nCellX+1)
        Do j = 2,(nCellY+1)
           !- conservative variable (masse et quantité de mvt)
           Utemp3 = (/ rho(i,j) , rho(i,j)*u(i,j), rho(i,j)*v(i,j) /)
           !- Roe method (U^{n+1/2})
           Utemp3 = Utemp3 - dtDivDx*(/ &
                FluxSplitRhoX(i,j) - FluxSplitRhoX(i-1,j), &
                FluxSplitUx(i,j)   - FluxSplitUx(i-1,j)  , &
                FluxSplitVx(i,j)   - FluxSplitVx(i-1,j)   /)
           !- relaxation term (Omega^{n+1})
           Omega = (/ Utemp3(2)/Utemp3(1) , Utemp3(3)/Utemp3(1) /)
           if (P%etaRelax==0d0) Then
              !- Omega of unit length
              Omega = Omega/sqrt(Omega(1)**2 + Omega(2)**2)
           else
              !- Omega relax
              C_0 = 1/Norm_2(Omega) - 1d0
              Omega =  (1d0/sqrt(1d0+C_0*exp(-2d0/P%etaRelax/P%dt)))/sqrt(Norm_2(Omega))*Omega
           endif
           !- update
           rho(i,j) = Utemp3(1)
           u(i,j)   = Omega(1)
           v(i,j)   = Omega(2)
        end Do
     end Do
     !--      A.3) Boundary condition      --!
     !---------------------------------------!
     Call BoundaryCondition_2D(rho,P%boundCondX,P%boundCondY)
     if (P%boundCondX==3 .and. P%boundCondY==3) then
        Call BC_Vortex_2D(u,v)
     endif     
     Call BoundaryCondition_2D(u,P%boundCondX,P%boundCondY)
     Call BoundaryCondition_2D(v,P%boundCondX,P%boundCondY)

    
     !----------------------------------------------------!
     !---------------------   B) yYy   -------------------!
     !----------------------------------------------------!

     !--    B.1) Flux (calcul des valeurs aux interfaces)   --!
     !--------------------------------------------------------!
     Do i = 2,nCellX+1
        Do j = 1,nCellY+1
           !- en j-1/2: in the basis {e_y,-e_x}
           Ftemp3 = FluxSplit_x( &
                rho(i,j),rho(i,j+1),v(i,j),v(i,j+1),-u(i,j),-u(i,j+1), &
                P%c1,P%c2,P%ld,P%degreePolyScheme)
           ! we get back in the basis {e_x,e_y}
           FluxSplitRhoY(i,j) = Ftemp3(1)
           FluxSplitVy(i,j)   = Ftemp3(2)
           FluxSplitUy(i,j)   = -Ftemp3(3)
        End do
     End do
     !--     B.2) New value (on avance d'un 1/2 pas de temps)      --!
     !---------------------------------------------------------------!
     Do i = 2,(nCellX+1)
        Do j = 2,(nCellY+1)
           !- conservative variable (masse et quantité de mvt)
              Utemp3 = (/ rho(i,j) , rho(i,j)*u(i,j), rho(i,j)*v(i,j) /)
              !- Roe method (U^{n+1/2})
              Utemp3 = Utemp3 - dtDivDy*(/ &
                   FluxSplitRhoY(i,j) - FluxSplitRhoY(i,j-1), &
                   FluxSplitUy(i,j)   - FluxSplitUy(i,j-1)  , &
                   FluxSplitVy(i,j)   - FluxSplitVy(i,j-1)   /)
              !- relaxation term (Omega^{n+1})
              Omega = (/ Utemp3(2)/Utemp3(1) , Utemp3(3)/Utemp3(1) /)
              if (P%etaRelax==0d0) Then
                 !- Omega of unit length
                 Omega = Omega/sqrt(Omega(1)**2 + Omega(2)**2)
              else
                 !- Omega relax
                 C_0 = 1/Norm_2(Omega) - 1d0
                 Omega =  (1d0/sqrt(1d0+C_0*exp(-2d0/P%etaRelax/P%dt)))/sqrt(Norm_2(Omega))*Omega
              endif
              !- update
              rho(i,j) = Utemp3(1)
              u(i,j)   = Omega(1)
              v(i,j)   = Omega(2)
        end Do
     end Do
     !--      B.3) Boundary condition      --!
     !---------------------------------------!
     Call BoundaryCondition_2D(rho,P%boundCondX,P%boundCondY)
     if (P%boundCondX==3 .and. P%boundCondY==3) then
        Call BC_Vortex_2D(u,v)
     endif
     Call BoundaryCondition_2D(u,P%boundCondX,P%boundCondY)
     Call BoundaryCondition_2D(v,P%boundCondX,P%boundCondY)


     !---------------------------------!
     !--        Output (print)       --!
     !---------------------------------!
     if (P%shouldSaveAll .or. iTime==nTime) Then
        !-- pour shouldSaveAll=false, on garde que le 1er et dernier plot
        if (P%isFormatVtk) then
           Call FilePrintArray2D_vtk(rho,u,v,P%dx,P%dy,"../data/rho_uv_2D_",.true.,iTime)
        else
           Call FilePrint_2D(rho,AngleVec_N2(u,v),iTime)
        endif
     end if
     !- progress...
     call BarProgress(iTime,nTime)

  end Do

  !---------------------------------------------------------------!
  !----------------------- end loop in time ----------------------!
  !---------------------------------------------------------------!

  !- finishing line...
  Call Cpu_time(finish)
  print "(A24,f11.3)"," Time to compute (s)   = ",finish-start
  print *,"******************************************************"
  print *,""

  !- Deallocate
  deallocate(rho,u,v)
  deallocate(FluxSplitRhoX,FluxSplitUx,FluxSplitVx)
  deallocate(FluxSplitRhoY,FluxSplitUy,FluxSplitVy)

end program MacroVic_2D
