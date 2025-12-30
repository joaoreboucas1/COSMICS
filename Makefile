# Makefile for initial conditions codes.

SHELL   = /bin/sh

# Subdirectories to make recursively
SUBDIR  = linger_con linger_syn deltat grafic

all help:
	@echo
	@echo "Usage: make sys-type"
	@echo
	@echo "  where sys-type is one of the following:"
	@echo
	@echo "dec alpha axp osf convex cray unicos hp hpux"
	@echo "linux ibm rs6k rs6000 aix sgi irix sun sunos generic"
	@echo


dec alpha axp osf:
	@make ARCH=DEC doall

convex:
	@make ARCH=CONVEX doall

cray unicos:
	@make ARCH=CRAY doall

hp hpux:
	@make ARCH=HP doall

ibm rs6k rs6000 aix:
	@make ARCH=RS6K doall

pc linux:
	@make ARCH=LINUX doall

sgi irix:
	@make ARCH=SGI doall

sun sunos:
	@make ARCH=SUN doall

generic:
	@make ARCH=GENERIC doall

doall:
	@if [ ! -d ./bin ]; then ( mkdir -p ./bin ) ; fi;
	@echo "****************** MAKING SUBDIRECTORIES ****************";
	@for i in ${SUBDIR}; do (echo "*** COMPILING $$i DIRECTORY"; cd $$i; \
           $(MAKE) "ARCH=$(ARCH)"; cp $$i ../bin ); done;
	@echo "************************* DONE **************************";

test:	
	@if [ ! -d ./test_results ]; then ( mkdir -p ./test_results ) ; fi;
	@csh test.csh
	@echo "test results are in ./test_results"

realclean:	clean cleangood

superclean:	clean cleangood cleanclean

clean:
	@for i in ${SUBDIR}; do (cd $$i; $(MAKE) -i clean); done

cleangood:
	@if [ ! -d bin ] ; then true; \
	      else /bin/rm -r bin ; fi
	@if [ ! -d test_results ] ; then true; \
	      else /bin/rm -r test_results ; fi
	@echo "Deleting *~ #* core *.o"
	@find . \( -name \*~ -o -name \#\* -o -name core \) \
	   -exec /bin/rm {} \; -print

cleanclean:
	@echo "Deleting *.dat fort.* *.0 *.1 *.log cosmics_temp.f"
	@find . \( -name \*.dat -o -name fort.\* -o -name \*.0 -o -name \*.1  -o -name \*.log -o -name cosmics_temp.f \) \
	   -exec /bin/rm {} \; -print

