FC=gfortran

HOMEDIR=../../../

CFLAGS += -O3
CFLAGS += -fbounds-check -fcheck-array-temporaries -Wall
#CFLAGS += -DDEBUG

PROGRAM=$(HOMEDIR)bin/solver

.PHONY: all cmp run gridgen clean cp

	
all: cmp run

cmp:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

run: git.bin
	@echo 'Running Grid-Adaper on simple Chamber'
	@$(PROGRAM) #-debug
	@$(HOMEDIR)bin/show_sol

git.bin: gridgen
	@./$<
	
gridgen: gridgen.F90
	@$(FC) $(CFLAGS) $< -o $@
	
clean:
	@rm -vf git_out.bin *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin gridgen
	@rm -vrf *~

##### CLEAN PROJECT
cp: 
	@echo "CLEANING SOLVER DIR"
	@make -C $(HOMEDIR) clean --no-print-directory
