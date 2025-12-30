c  Parameters for grafic.  User must set np1,np2,np3,ncon.
c  ** Important: if no constraints are desired, set ncon=1 **
c  ** (not 0) and input ido=1. **
c
	parameter (np1=32,np2=32,np3=32,ncon=1)
c
	parameter (npart=np1*np2*np3)
c
c
c 
c np[1-3]:      particle lattice size for computing initial conditions.
c                 All must be divisible by 2.
c ncon:         number of constraints applied to the initial density field.
c                 Used only if ido>1; if so, user must supply correct constr.f.
