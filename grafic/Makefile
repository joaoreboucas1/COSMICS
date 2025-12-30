# type 'make ARCH=xxx'
# where xxx is one of: CONVEX,CRAY,DEC,HP,RS6K,SGI,SUN,GENERIC
#

include ../Make_files/Make.$(ARCH)

#No need to touch anything below this line.
#-------------------------------------------------------------------------------

ARCH_FLAG = -DARCH_IS_$(ARCH)
F77_OBJS = grafic.o ic3.o constr.o power.o time.o util.o random8.o $(FFT_OBJ)

all:	grafic

grafic:	$(F77_OBJS)
	$(F77) $(F77FLAGS) $(F77_OBJS) $(LIBS) -o grafic

clean:
	rm grafic *.o

.SUFFIXES: .f

.f.o:
	$(F77) $(F77FLAGS) -c $<

grafic.o:	grafic.inc
ic3.o:		grafic.inc
constr.o:	grafic.inc
