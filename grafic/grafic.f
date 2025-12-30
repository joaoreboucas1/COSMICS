cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	program grafic
c  Generate initial conditions for cosmological N-body integration
c  as a Gaussian random field.
c
	include 'grafic.inc'
c
	real x(3,npart),v(3,npart),delta(npart)
	double precision sigma,sigd,ekin0
	common /cosmoparms/ omegam,omegav,h0
c
c  npart is the number of particles; initial conditions give displacements
c  from a grid of size np1*np2*np3=npart.  dx is the mean interparticle
c  spacing in comoving Mpc; omegam, omegav, and h0 are at a=1.
c  Initial positions and velocities are provided at a=astart<1.
c
c  Initialize power spectrum.
	call pini
c
	write(*,*) 'Particle lattice size: np1,np2,np3=',np1,np2,np3
c
	write(*,*) 'Enter dx (Mpc, =1 for scale-free), epsilon (Mpc,',
     &             ' =0 for PM), etat (=0.05)'
	read(*,*) dx,epsilon,etat
	write(*,1001) 'npart, L_x, L_y, L_z=',npart,np1*dx,np2*dx,
     2             np3*dx,' Mpc'
1001	format(1x,a21,i9,3(2x,f7.2),a4)
	amass0=2.776e7*omegam*h0**2*dx**3
	write(*,1002) 'Particle mass= ',amass0,' solar masses'
1002	format(1x,a15,e9.4,a13)
c  Choose your own random number seed.
	write(*,*) 'Enter random number seed (9-digit integer)'
	read(*,*) iseed
	call randini(iseed)
	write(*,*) 'Enter ido (=1,2,3) for action:'
	write(*,*) '  ido=1 Compute unconstrained realization'
	write(*,*) '  ido=2 Compute mean field of constraints'
	write(*,*) '  ido=3 Compute constrained realization'
	read(*,*) ido
	if (ido.lt.1.or.ido.gt.3) stop
	if (ido.gt.1) then
	  write(*,*)
     2 'Make sure you have provided constraints in subroutine constr'
	end if
c  Compute realizations of density perturbation and displacement fields.
	call ic3(ido,dx,delta,v,delm,dispm)
c  Scale density and displacement back to starting redshift.
c  Require that max |delta_rho/rho| <= delmax at a=astart.
	delmax=1.0
	dpl=delmax/delm
cc  Require that max displacement <= deltax at a=astart.
c	dpl=dx/dispm
c  No, actually scale it so that maximum delta=1 at a=astart.
c  To change it to displacement scaling, uncomment the previous line.
	dpls=dpl*dplus(1.0,omegam,omegav)
	astart=adp(dpls,omegam,omegav)
	write(*,*)
	write(*,*) 'Scaling density and displacements to a=',astart
	xfact=dpl/dx
	vfact=xfact*fomega(astart,omegam,omegav)*
     2        dladt(astart,omegam,omegav)/astart
c
	sigma=0.0
	delm=0.0
	sigd=0.0
	dispm=0.0
	ekin0=0.0
	j=0
	  do 30 j3=1,np3
	  x3=(j3-1)
	    do 20 j2=1,np2
	    x2=(j2-1)
	      do 10 j1=1,np1
	      x1=(j1-1)
	      j=j+1
c  Compute statistics.
	      sigma=sigma+delta(j)**2
	      delm=max(delm,abs(delta(j)))
	      dsigd=v(1,j)**2+v(2,j)**2+v(3,j)**2
	      sigd=sigd+dsigd
	      dispm=max(dispm,dsigd)
c  Rescale density perturbation.
	      delta(j)=delta(j)*dpl
c  Apply displacements.
	      x(1,j)=x1+v(1,j)*xfact
	      x(2,j)=x2+v(2,j)*xfact
	      x(3,j)=x3+v(3,j)*xfact
c  Periodic boundary conditions.
	      if (x(1,j).lt.0.0) x(1,j)=x(1,j)+np1
	      if (x(1,j).ge.np1) x(1,j)=x(1,j)-np1
	      if (x(2,j).lt.0.0) x(2,j)=x(2,j)+np2
	      if (x(2,j).ge.np2) x(2,j)=x(2,j)-np2
	      if (x(3,j).lt.0.0) x(3,j)=x(3,j)+np3
	      if (x(3,j).ge.np3) x(3,j)=x(3,j)-np3
c  Compute velocity.
	      v(1,j)=v(1,j)*vfact
	      v(2,j)=v(2,j)*vfact
	      v(3,j)=v(3,j)*vfact
c  Compute kinetic energy.
	      ekin0=ekin0+v(1,j)**2+v(2,j)**2+v(3,j)**2
10	    continue
20	  continue
30	continue
	sigma=sqrt(sigma/npart)*dpl
	delm=delm*dpl
	sigd=sqrt(sigd/npart)*dpl
	dispm=sqrt(dispm)*dpl
	write(*,*)
	write(*,*) 'For a=astart: linear sigma, delmax=',real(sigma),delm
	write(*,*) 'RMS, max. 3-D displacement=',real(sigd),
     2             dispm,' Mpc'
c  Set parameters for output p3m file.
	nstep=0
c  First timestep: max. displacement=epsilon (or 0.5*dx for PM).
c  Note that dt := H0*dtau/a where dtau is conformal time -- this t is
c  not proper time!
	dterm=fomega(astart,omegam,omegav)*dladt(astart,omegam,omegav)
	if (epsilon.gt.0.0) then
	  dt=(epsilon/dispm)/dterm
	else
	  dt=(0.5*dx/dispm)/dterm
	end if
c  Cosmological timestep criterion: limit timestep to 0.1*Hubble time.
	dterm=0.1/dladt(astart,omegam,omegav)
	dt=min(dt,dterm)
c  Proper specific energies in units of (km/s)**2.
	ekin=0.5*ekin0/npart*(h0*dx)**2
	egrav=-1.5*ekin
	egint=egrav/3.0
c
c  Output delta.dat file.
	open(10,file='delta.dat',status='unknown',form='unformatted')
	rewind 10
	write(10) np1,np2,np3,dx,astart,omegam,omegav,h0
	write(10) (delta(j),j=1,npart)
	close(10)
c  Output p3m.dat file.
	open(10,file='p3m.dat',status='unknown',form='unformatted')
	rewind 10
	write(10) npart,np1,np2,np3,dx,epsilon,astart,omegam,omegav,h0,
     2                  dt,etat,nstep,ekin,egrav,egint
c  Note that do 30 has changed the units of x to dx, while v has units
c  of (H0*dx).  Correct them back to dimensional units.
	write(10) (x(1,j)*dx,j=1,npart)
	write(10) (x(2,j)*dx,j=1,npart)
	write(10) (x(3,j)*dx,j=1,npart)
	write(10) (v(1,j)*h0*dx,j=1,npart)
	write(10) (v(2,j)*h0*dx,j=1,npart)
	write(10) (v(3,j)*h0*dx,j=1,npart)
	close(10)
c
	stop
	end
