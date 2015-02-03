module boundary_1D

  !- contains :
  !    boundary condition
  
contains
  
  Subroutine BoundaryCondition_1D(density,boundCond)
    !- modifie les densités pour avoir les conditions aux bords
    implicit none
    Double Precision, Dimension(:), intent(inout) :: density
    Integer, intent(in)                           :: boundCond
    Integer                                       :: n_x
    !- the subroutine
    n_x = size(density)

    Select case (boundCond)
    case(1)
       !- periodic
       density(1)   = density(n_x-1)
       density(n_x) = density(2)
    case(2)
       !- Neumann homogeneous
       density(1)   = density(2)
       density(n_x) = density(n_x-1)
    end Select

  end Subroutine BoundaryCondition_1D

  
End module boundary_1D
