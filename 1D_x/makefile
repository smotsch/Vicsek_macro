################################
##    makefile MacroVic_1D    ##
################################

##----------------------------------------##
##         main files to compile          ##
##----------------------------------------##
## objects
OBJECTS = elementary.o        \
       matrix_FVM.o           \
       flux_FVM.o             \
       type_PARAM.o           \
       input_output_1D.o      \
       boundary_1D.o          \
       initial_condition_1D.o \
       main_MacroVic_1D.o
## compiler
CC = gfortran
## options for the compiler (-g for debugger, parano: -check all -warn all),
FLAGS = -Wall -fbackslash -O2
## name of the executable
EXEC = MacroVic_1D


##----------------------------------------##
##            jobs for 'make'             ##
##----------------------------------------##

## Tricks to use make with several folders
##----------------------------------------
# Source tree:
#    Makefile
#    -- src/   *.f90
#    -- build/ *.o *.mod
#    -- bin/   $(EXE)
# 
# Paths
OBJDIR = build
vpath %.f90 src
vpath %.o $(OBJDIR)
# to help make find the object in the folder $(OBJDIR)
OBJECTS2 = $(addprefix $(OBJDIR)/, $(OBJECTS))

## command
all: init_folder $(EXEC) move init_parameters

$(EXEC): $(OBJECTS)
	$(CC) $(FLAGS) -o $(EXEC) $(OBJECTS2)

%.o : %.f90
	$(CC) $(FLAGS) -c $^ -o $(OBJDIR)/$@ -J$(OBJDIR)

.PHONY: init_folder init_parameters move clean
init_folder:
	mkdir -p $(OBJDIR) data images videos
init_parameters:
	@if [ ! -f 'bin/PARAMETER_init.txt' ]; then \
	  cp bin/parameters/PARAMETER_init_bak.txt bin/PARAMETER_init.txt; \
	fi
	@if [ ! -f 'bin/PARAMETER_1D.txt' ]; then \
	  cp bin/parameters/PARAMETER_1D_bak.txt bin/PARAMETER_1D.txt; \
	fi
move:
	mv $(EXEC) bin
clean:
	rm -rf $(OBJDIR)/*.mod
	rm -rf $(OBJDIR)/*.o


##----------------------------------------##
##      remark: changes for ifort         ##
##----------------------------------------##
# 1. CC = ifort
# 2. FLAGS = -assume bscc -msse2 -O2
# 3. $(OBJECTS): %.f90
#	$(CC) -c $< -o $(OBJDIR)/$@ -module $(OBJDIR)
