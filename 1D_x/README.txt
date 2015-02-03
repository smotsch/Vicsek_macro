       **************************************************
       *                 Documentation                  *
       *                                                *
       *	        for the program                 *
       *                                                *
       *	            "MacroVic_1D"               *
       *                                                *
       *                Seb, 02/02/2015                 *
       **************************************************

       
1) Installation :
-----------------
You just have to execute the makefile:
 > make
The compiler creates the executable file"./MacroVic_1D" in the folder 'bin'.

Remark: by default, the makefile uses the gfortran compiler. If you want to use
ifort, open the 'makefile' and use the instruction written at the bottom. 


2) Program execution:
---------------------
To execute the program, you just have to run MacroVic_1D:
 > ./MacroVic_1D
The program is going to read the parameters of the model in the files:
   . PARAMETER_1D: for most of the parameters (numerical method, Δt...)
   . PARAMETER_init: for the initial condition (uniform, Gaussian...)
The parameters are written in external files because we do not need to recompile when we
change one parameter by doing so. The problem with this method is that we have to write
the value of the parameters at a given line. We cannot add a comment line in the file
PARAMETER_1D, otherwise the numbering is ruined.
 During the execution, the program  "./MacroVic_1D" is going to compute the density (ρ) and
the average velocity (theta) of the solutions at each time step using the algorithm describing in
the article [1]. In the terminal, the program gives some information about the parameters
used for the simulation. At the end, it displays the computation time.

[1] S. Motsch, L. Navoret, "Numerical simulations of a non-conservative hyperbolic system
with geometric constraints describing swarming behavior" (2011).


3) Output:
----------
The program MacroVic_1D can store two types of data. First, it can save the trajectories and
the velocities of each particle over time. At each time step, it creates in the directory
"data" the files:
   . rho_******     : x coordinate of the particles
   . theta_******   : y coordinate of the particles
with ****** a counter of time step. The program writes only if the parameter
"shouldSaveAll" is "True" (line 25 in PARAMETER_1D).


4) Visual output:
-----------------
To display the results of the computations, run the script:
 > Display_MacroVic_1D.m
in the folder 'visualization' with Octave. For Matlab, the script is called
'Display_MacroVic_1D_Matlab.m'.

5) The parameters:
------------------
 PARAMETER_1D
   + c1, c2, ld    : coefficients of the macroscopic model
   + Lx, dx	   : size of the domain in x and meshgrid
   + Time, dt	   : total time and Δt
   + boundCond	   : boundary condition (periodic or Neumann)
   + methodNum	   : 4 choices of numerical methods, degree is
   + degree        : order for the polynomial upwind (0,1,2) in estimating |DF|
                     -> degree==3 means no approximation
   + epsilon       : in the splitting method, source term is in 1/ε
                     -> epsilon==0 leads to a normalization operator
   + shouldSaveAll : save all the data in time
 PARAMETER_init
   + choiceInit     : choice for the initial condition for X
   + meanRho,meanTheta,rangeRho,rangeTheta: random values in x
   + rhoL,rhoR,thetaL,thetaR : values for the Riemann problem
   + ...

6) More details about the program:
----------------------------------
The main program is the file:
  MacroVic_1D.f90  
It uses different modules (defined in separated files) :
  . elementary           : contains the usual functions (AngleVec, RandNorm...)
  . input_output_1D      : for the input/output (Lecture, FilePrint...)
  . boundary_1D          : the effect of the wall (Wall)
  . initial_condition_1D : to initialize with the proper initial conditions (InitCond)
  . matrix_FVM/flux_FVM  : compute the flux at the interface between cells
The architecture of the program is the following:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%    main_MacroVic_1D                                                        %%
  %%      -> declaration of variables                                           %%
  %%      -> lecture of parameters ("Lecture")                                  %%
  %%      -> initialization of variables ("InitCond")                           %%
  %%                                                                            %%
  %%      -> loop in time                                                       %%
  %%        1) computation of flux at the interface of the cells                %%
  %%        2) new values of ρ,θ                                                %%
  %%        3) boundary conditions                                              %%
  %%        4) write data ("FilePrint")                                         %%
  %%      -<                                                                    %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
To make the program easier to read, two structures are used to store the parameters ('PARAM_MacroVic' and 'PARAM_InitMV'). This avoids to write 36 arguments each time a subroutine is called.

7) Legal notice information:
----------------------------
 This (modest) program is distributed under the GNU GPL license version 2. More
information are available in the file COPYING.txt. For any information or
bugs, please contact me at: <smotsch at asu.edu>.
