Program *MacroVic_2D*
=======================


## Installation

Execute the `makefile`:
```bash
	> make
```
The compiler creates the executable file `./MacroVic_2D` and move it to the
directory `bin`. By default, the makefile uses the *gfortran* compiler. If you
want to use *ifort*, open the `makefile` and use the instruction written at the
bottom.

## Program execution

To execute the program, run the binary file:
```bash
	> ./MacroVic_2D
```
The program is going to read the parameters of the model in the files:
* `PARAMETER_MacroVic_2D.txt`: for most of the parameters
* `PARAMETER_init.txt`: for the initial condition (uniform, Gaussian...)

The parameters are written in external files because we do not need to recompile
when we change one parameter by doing so. The problem with this method is that
we have to write the value of the parameters at a given line. We cannot add a
comment line in the file `PARAMETER_MacroVic_2D.txt`, otherwise the numbering is
ruined.

During the execution, the program `MacroVic_2D` is going to compute the density
`ρ` and average velocity (`θ`) of the particles at each time step using the
algorithm describing in the article [1]. In the terminal, the program gives some
information about the parameters used for the simulation. At the end, it
displays the computation time.

[1] S. Motsch, L. Navoret, *"Numerical simulations of a non-conservative
hyperbolic system with geometric constraints describing swarming behavior"* (2011).


Output
------
The program `MacroVic_2D` can store two types of data. First, it can save the
trajectories and the velocities of each particle over time. At each time step,
it creates in the directory `data` the files:
* `rho_******`     : density 
* `theta_******`   : direction of mean velocity

with `******` a counter of time step. The program writes only if the parameter
`shouldSaveAll` is `True` (line 27 in `PARAMETER_2D.txt`).


Visual output
-------------
To display the results of the computations with **Octave**
```bash
	> Display_HydroVic_2D.m
```
in the folder `visualization`.

The parameters
--------------
* `PARAMETER_2D.txt`
 * `c1`, `c2`, `ld`    : coefficients of the macroscopic model
 * `Lx`, `Ly`	   : size of the domain in x and y
 * `dx`, `dy`	   : meshgrid in x an dy
 * `Time`, `dt`	   : total time and Δt
 * `boundCond*`	   : boundary condition (periodic or Neumann)
 * `degree`        : order for the polynomial upwind (0,1,2) in estimating |DF| (`degree=3` means no approximation)
 * `epsilon`       : in the splitting method, source term is in 1/ε (`ε=0` leads to a normalization operator)
 * `shouldSaveAll` : save all the data in time
 * `isFormatVtk`   : save in VTK format or Octave
* `PARAMETER_init.txt`
 * `choiceInit`     : choice for the initial condition for `X`
 * `meanRho`,`meanTheta`,`rangeRho`,`rangeTheta`: random values in `x`
 * `rhoL`,`rhoR`,`thetaL`,`thetaR` : values for the Riemann problem
 
More details about the program
------------------------------
The main program is the file `MacroVic_1D.f90`. It uses different modules
(defined in separated files):

| File                         | Description   |
| -----------------------------|:-------------:|
| `elementary`           | contains the usual functions (AngleVec, RandNorm...) |
| `input_output_2D`      | for the input/output (Lecture, FilePrint...)         |
| `boundary_2D`          | the effect of the wall (Wall)                        |
| `initial_condition_2D` | to initialize with the proper initial conditions (InitCond)|
| `matrix_FVM/flux_FVM`  | compute the flux at the interface between cells      |

The architecture of the program is the following:
```bash
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%    main_MacroVic_2D                                                 %%
  %%      -> declaration of variables                                    %%
  %%      -> lecture of parameters ("Lecture")                           %%
  %%      -> initialization of variables ("InitCond")                    %%
  %%                                                                     %%
  %%      -> loop in time                                                %%
  %%        1) computation of flux at the interface of the cells         %%
  %%        2) new values of ρ,θ                                         %%
  %%        3) boundary conditions                                       %%
  %%        4) write data ("FilePrint")                                  %%
  %%      -<                                                             %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```
To make the program easier to read, two structures are used to store the
parameters (`PARAM_MacroVic` and `PARAM_InitMV`). This avoids to write 36
arguments each time a subroutine is called.

Legal notice information
------------------------
 This (modest) program is distributed under the GNU GPL license version 2. More
information are available in the file COPYING.txt. For any information or bugs,
please contact me at: <smotsch at asu.edu>.
