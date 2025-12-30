# type 'make ARCH=xxx'
# where xxx is one of: CONVEX,CRAY,DEC,HP,LINUX,RS6K,SGI,SUN,GENERIC
#

include ../Make_files/Make.$(ARCH)

#No need to touch anything below this line.
#-------------------------------------------------------------------------------

ARCH_FLAG = -DARCH_IS_$(ARCH)

F77_OBJS = linger_syn.o dverk.o

all:	fix_generic linger_syn

fix_generic:
	@if [ -f cosmics_temp.f ] ; then \
	  rm cosmics_temp.f; fi
	@if [ "$(ARCH)" = "GENERIC" ] ; then \
	  mv linger_syn.f.generic linger_syn.f; fi

linger_syn:	$(F77_OBJS)
		$(F77) $(F77FLAGS) $(F77_OBJS) -o linger_syn

clean:
	@if [ -f linger_syn.f ] ; then \
	  mv linger_syn.f linger_syn.f.generic; fi
	rm linger_syn *.o cosmics_temp.f

.SUFFIXES: .f .F

.f.o:
	$(F77) $(F77FLAGS) -c $<

# the ibm compiler won't pass $(ARCH_FLAG) to the preprocessor.
# the cray compiler won't allow both -c and -o together.
# thus we apply the preprocessor seperately and mv the object file.
.F.o:
	@if [ "$(ARCH)" = "LINUX" ]; then \
	  $(CPP) -E $(ARCH_FLAG) $< -o cosmics_temp.f; \
	else \
	  $(CPP) $(ARCH_FLAG) $< cosmics_temp.f; fi
	$(F77) $(F77FLAGS) -c cosmics_temp.f
	mv cosmics_temp.o $*.o
	\rm -rf cosmics_temp.f
