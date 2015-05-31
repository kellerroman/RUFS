FC=gfortran

HOMEDIR=../../../

CFLAGS += -O3
CFLAGS += -fbounds-check -fcheck-array-temporaries -Wall
#CFLAGS += -DDEBUG
PRE_CONFIG=$(HOMEDIR)tools/pre_config.py
PROGRAM=$(HOMEDIR)bin/solver

.PHONY: all cmp run gridgen clean cp

	
all: cmp run

cmp:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

str.dx: ../stubex
	@./$<
	
../stubex: ../stubex.f
	@$(FC) $(CFLAGS) $< -o $@


run: git.bin str.dx
	@echo 'Running Grid-Adaper on simple Chamber'
	@$(PRE_CONFIG)
	@$(PROGRAM) #-debug
	@$(HOMEDIR)bin/extract_1D x
	@$(HOMEDIR)bin/sol2visit

git.bin: gridgen
	@./$<
	

	
gridgen: gridgen.F90
	@$(FC) $(CFLAGS) $< -o $@
	
clean:
	@rm -vf git_out.bin *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.silo *.dx stubex
	@rm -vrf *~

##### CLEAN PROJECT
cp: 
	@echo "CLEANING SOLVER DIR"
	@make -C $(HOMEDIR) clean --no-print-directory
