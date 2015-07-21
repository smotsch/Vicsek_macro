module type_PARAM

    TYPE PARAM_MacroVic
     Double Precision                   :: c1, c2, ld
     Double Precision                   :: Lx, dx
     Double Precision                   :: Time, dt
     Integer                            :: boundCond
     Integer                            :: methodNum, degree
     Double Precision                   :: epsilon
     Logical                            :: shouldSaveAll
     Logical                            :: isFormatVtk

     Integer                            :: nSpace
  end TYPE PARAM_MacroVic

  TYPE PARAM_InitMV
     Integer                            :: choiceInit
     Double Precision                   :: meanRho, meanTheta, rangeRho, rangeTheta
     Double Precision                   :: rhoL, rhoR, thetaL, thetaR
     Double Precision                   :: rhoBG, thetaBG
     Double Precision                   :: rho1, theta1, x1Center, radius1
     Double Precision                   :: rho2, theta2, x2Center, radius2
     Double Precision                   :: rho3, theta3, x3Center, radius3
     Character*80                       :: cheminRho, cheminTheta
  end TYPE PARAM_InitMV

end module type_PARAM
