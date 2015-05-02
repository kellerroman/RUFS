FC=gfortran

HOMEDIR=../../

CFLAGS += -O3
CFLAGS += -fbounds-check -fcheck-array-temporaries -Wall
CFLAGS += -DDEBUG

PROGRAM=$(HOMEDIR)bin/solver

.PHONY: all cmp run gridgen clean cp

	
all: cmp run

cmp:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

run: gridgen
	@echo 'Running Grid-Adaper on simple Chamber'
	@ln -sf gridgen/git.bin
	@ln -sf gridgen/bc.bin
	@$(PROGRAM) #-debug
	@$(HOMEDIR)bin/sol2visit

gridgen:
	@$(MAKE) -C gridgen run --no-print-directory
	
clean:
	@rm -vf git_out.bin *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin *.silo
	@rm -vrf *~
	@$(MAKE) -C gridgen clean --no-print-directory

##### CLEAN PROJECT
cp: 
	@echo "CLEANING SOLVER DIR"
	@make -C $(HOMEDIR) clean --no-print-directory
