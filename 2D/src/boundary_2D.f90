module boundary_2D
  
contains
  
  Subroutine BoundaryCondition_2D(density,boundCondX,boundCondY)
    !- change the values of the 2 extra cells to have the correct BC
    implicit none
    Double Precision, Dimension(:,:), intent(inout) :: density
    Integer, intent(in)                             :: boundCondX, boundCondY
    Integer                                         :: nCellX, nCellY
    !- the subroutine
    nCellX = size(density(:,1)) - 2 ! 2 extra values for the BC
    nCellY = size(density(1,:)) - 2

    Select case (boundCondX)
    case (1)
       !- periodic
       density(1,:)        = density(nCellX+1,:)
       density(nCellX+2,:) = density(2,:)
    case (2,3)
       !- Neumann homogeneous
       density(1,:)        = density(2,:)
       density(nCellX+2,:) = density(nCellX+1,:)
    end Select

    Select case (boundCondY)
    case (1)
       !- periodic
       density(:,1)        = density(:,nCellY+1)
       density(:,nCellY+2) = density(:,2)
    case (2,3)
       !- Neumann homogeneous
       density(:,1)        = density(:,2)
       density(:,nCellY+2) = density(:,nCellY+1)
    end Select
    
  end Subroutine BoundaryCondition_2D

  Subroutine BC_Vortex_2D(u,v)
    !- modifie u,v pour que ça tourne !
    implicit none
    Double Precision, Dimension(:,:), intent(inout) :: u,v
    Integer                                         :: nCellX, nCellY
    !- the subroutine
    nCellX = size(u(:,1)) - 2 ! 2 extra values for the BC
    nCellY = size(u(1,:)) - 2

    ! the order is important !!
    u(2:(nCellX+1),2)        =  1d0 !   bottom  (x   , y=0)
    v(2:(nCellX+1),2)        =  0d0
    u(nCellX+1,2:(nCellY+1)) =  0d0 !   right   (x=l , y  )
    v(nCellX+1,2:(nCellY+1)) =  1d0    
    u(2:(nCellX+1),nCellY+1) = -1d0 !   top     (x   , y=l)
    v(2:(nCellX+1),nCellY+1) =  0d0    
    u(2,2:(nCellY+1))        =  0d0 !   left    (x=0 , y  )
    v(2,2:(nCellY+1))        = -1d0
    ! The corners
    u(2,2)                   =  0d0
    v(2,2)                   =  0d0
    u(nCellX+1,2)            =  0d0
    v(nCellX+1,2)            =  0d0
    u(nCellX+1,nCellY+1)     =  0d0
    v(nCellX+1,nCellY+1)     =  0d0
    u(2,nCellY+1)            =  0d0
    v(2,nCellY+1)            =  0d0
    
    
  end Subroutine BC_Vortex_2D

  
End module boundary_2D
