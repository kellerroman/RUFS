FC=gfortran

HOMEDIR=../../

CFLAGS += -O3
CFLAGS += -fbounds-check -fcheck-array-temporaries -Wall
CFLAGS += -DDEBUG
PRE_CONFIG=$(HOMEDIR)tools/pre_config.py
PROGRAM=$(HOMEDIR)bin/solver

.PHONY: all cmp run gridgen clean cp

	
all: cmp run

cmp:
	@$(MAKE) -C $(HOMEDIR) --no-print-directory

preconfig: config.bin

config.bin: config.cfg
	@$(PRE_CONFIG)

run: gridgen preconfig
	@echo 'Running Grid-Adaper on simple Chamber'
	@$(PROGRAM) #-debug
	@$(HOMEDIR)bin/sol2visit

gridgen:
	@$(MAKE) -C gridgen run --no-print-directory
	@ln -sf gridgen/git.bin
	@ln -sf gridgen/bc.bin
	
clean:
	@rm -vf git_out.bin *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin *.silo
	@rm -vrf *~
	@$(MAKE) -C gridgen clean --no-print-directory

##### CLEAN PROJECT
cp: 
	@echo "CLEANING SOLVER DIR"
	@make -C $(HOMEDIR) clean --no-print-directory
