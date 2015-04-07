OBJECTS_DIR= obj
REALMAKEFILE=../src/Makefile.in
TOOLMAKEFILE=../tools/Makefile.tools

all: solver

solver: FORCE
	@(cd $(OBJECTS_DIR) && $(MAKE) -f $(REALMAKEFILE) --no-print-directory)

clean: FORCE
	@(cd $(OBJECTS_DIR) && $(MAKE) -f $(REALMAKEFILE) clean --no-print-directory)
	@rm -rf obj bin *~
	@$(MAKE) -C test clean --no-print-directory
	
FORCE:
	@mkdir -p obj bin

test:
	@$(MAKE) -C test --no-print-directorys

dep:
	@(cd src && ../tools/make_dependencies.py .) 
