C
C Developed by Edmund Bertschinger and Paul Bode
C              Deparment of Physics
C              MIT, Room 6-207
C              Cambridge, MA 02139
C              edbert@mit.edu
C
C  For conditions of distribution and use, see the accompanying 
C  cosmics/LICENSE file.
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	program deltot
	implicit double precision (a-h,o-z)
	parameter (lmax0=5000,lone=20,nkmax=10000)
	dimension deltat(nkmax,lone),dt(0:lmax0)
	dimension y(0:nkmax),dy(0:nkmax)
	integer kdone(nkmax)
	character*80 filename
	common /deltal/ y,dy,dk,an,l,lprint
	external ddelta

	  do i=1,nkmax
	  kdone(i)=-1
	end do
	nkdone=0

	write(*,*)
     2 'Enter lmax for C_l; spectral index n (scale-invariant is n=1)'
	read(*,*) lmax1,an
	if (lmax1.lt.2.or.lmax1.gt.lmax0) then
	  write(*,*) 'Illegal value of lmax: should lie between 2 and',
     2               lmax0
	  stop
	end if
	write(*,*)
     2   'As a diagnostic, deltat can write out one of the photon'
	write(*,*)
     2   '    transfer functions Delta_l(k) (i.e., for one value of l)'
	write(*,*) 'Enter l to write delta_l(k) to file (-1 for none)'
	read(*,*) lprint

	write(*,*) 'Enter lingerg.dat filename'
	read(*,'(a80)') filename

	open(10,file=filename,status='old',form='unformatted')
	rewind 10
	read(10) omegab,omegac,omegav,omegan,h0,tcmb,yhe,nnur,nnunr,
     2           initfl
	read(10) ifulbol
	read(10) akmin,akmax,nk,zend,tau
	if (ifulbol.eq.1) then
	  dk=akmax/nk
	else
	  dlk=log(akmax/akmin)/(nk-1)
	end if
c  Find which waves are done.
	  do i=1,nkmax
	  read(10,end=5) ik,ak,tau,lmax
	  lmax=min(lmax,lmax1)
	  read(10,end=5) (dt(l),l=0,lmax)
	  read(10,end=5) (dtp,l=0,lmax)
	  kdone(ik)=1
	end do
5	  do ik=1,nkmax
	  if (nkdone.eq.ik-1.and.kdone(ik).eq.1) nkdone=ik
	end do
	nk=nkdone
	write(*,*) 'nkdone=',nkdone
	fourpi=4.0*3.1415926535898d0
	tol=1.0d-7
	call splini
	open(11,file='deltat.dat',status='unknown')
	rewind 11
	write(11,'(4(1pe12.6,1x))') omegab,omegac,omegav,omegan
	write(11,'(3(1pe12.6,1x),3(i2,2x))') h0,tcmb,yhe,nnur,nnunr,
     2           initfl
	write(11,'(1pe12.6)') an

	if (lprint.ge.0) then
	  open(12,file='deltal.dat',status='unknown')
	  rewind 12
	  write(12,*) lprint
	end if

c  Outer loop over l, doing lone at a time.
	do lstart=1,lmax1,lone

c  Initialize array.
	  do l=1,lone
	    do k=1,nk
	    deltat(k,l)=0.0
	  end do
	end do
	lend=min(lstart+lone-1,lmax1)
	rewind 10
	read(10) omegab,omegac,omegav,omegan,h0,tcmb,yhe,nnur,nnunr,
     2           initfl
	read(10) ifulbol
	read(10) akmin,akmax,nk1,zend,tau

c  Read delta_l(k).
	  do ik=1,nk
	  read(10) k,ak,tau,lmax
	  lmax=min(lmax,lmax1)
	  read(10) (dt(l),l=0,lmax)
	  read(10) (dtp,l=0,lmax)
	    do l=lstart,lend
	    if (l.le.lmax) deltat(k,l-lstart+1)=dt(l)
	  end do
	end do

c  Inner loop over l.
	  do l=lstart,lend
	      do k=1,nk
	      y(k)=deltat(k,l-lstart+1)
	    end do
	if (l.eq.lprint) then
	  do k=1,nk
	  write(12,*) real(k*dk),real(y(k))
	  end do
	end if
	    y(0)=0.0d0
	    call splder(y,dy,nk+1)
	    cl=fourpi*fourpi*rombint(ddelta,0.0d0,nk*dk,tol)
	    write(11,*) l,l*(l+1)*cl
#ifndef ARCH_IS_CRAY 
#ifndef ARCH_IS_HP 
	    call flush(11)
#endif
#endif
c  End inner loop on l.
	  end do
  
c  End outer loop on l.
	end do
	close(11)
	if (lprint.ge.0) close(12)
c
	stop
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function ddelta(ak)
	implicit double precision (a-h,o-z)
	parameter (nkmax=10000)
	dimension y(0:nkmax),dy(0:nkmax)
	common /deltal/ y,dy,dk,an,l,lprint

	if (ak.eq.0.0d0) then
	  ddelta=0.0d0
	  return
	end if

	d=ak/dk
	i=d
	d=d-i

	deltl=y(i)+d*(dy(i)+d*(3.0*(y(i+1)-y(i))-2.0*dy(i)
     2   -dy(i+1)+d*(dy(i)+dy(i+1)+2.0*(y(i)-y(i+1)))))

	ddelta=deltl*deltl*ak**(an-2.0)

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine splder(y,dy,n)
c  Splder fits a cubic spline to y and returns the first derivatives at
c  the grid points in dy.  Dy is equivalent to a 4th-order Pade
c  difference formula for dy/di.
c
	implicit double precision (a-h,o-z)
	dimension f(10001),g(10001),y(n),dy(n)
	common /spline/ g
c
	n1=n-1
	if (n1.gt.10000)
     2       write(*,*) 'Spline array overflow!!! n1=',n1,'>10000'
c  Quartic fit to dy/di at boundaries, assuming d3y/di3=0.
	f(1)=(-10.0d0*y(1)+15.0d0*y(2)-6.0d0*y(3)+y(4))/6.0d0
	f(n)=(10.0d0*y(n)-15.0d0*y(n1)+6.0d0*y(n-2)-y(n-3))/6.0d0
c  Solve the tridiagonal system
c  dy(i-1)+4*dy(i)+dy(i+1)=3*(y(i+1)-y(i-1)), i=2,3,...,n1,
c  with dy(1)=f(1), dy(n)=f(n).
	  do 10 i=2,n1
	  f(i)=g(i)*(3.0d0*(y(i+1)-y(i-1))-f(i-1))
10	continue
	dy(n)=f(n)
	  do 20 i=n1,1,-1
	  dy(i)=f(i)-g(i)*dy(i+1)
20	continue
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine splini
c  Splini must be called before splder to initialize array g in common.
c
	implicit double precision (a-h,o-z)
	dimension g(10001)
	common /spline/ g
	save /spline/
c
	g(1)=0.0d0
	  do 10 i=2,10001
	  g(i)=1.0d0/(4.0d0-g(i-1))
10	continue
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function rombint(f,a,b,tol)
c  Rombint returns the integral from a to b of f(x)dx using Romberg integration.
c  The method converges provided that f(x) is continuous in (a,b).  The function
c  f must be double precision and must be declared external in the calling
c  routine.  tol indicates the desired relative accuracy in the integral.
c
	parameter (MAXITER=15,MAXJ=5)
	implicit double precision (a-h,o-z)
	dimension g(MAXJ+1)
	external f
c
	h=0.5d0*(b-a)
	gmax=h*(f(a)+f(b))
	g(1)=gmax
	nint=1
	error=1.0d20
	i=0
10	  i=i+1
	  if (i.gt.MAXITER.or.(i.gt.10.and.abs(error).lt.tol))
     2      go to 40
c  Calculate next trapezoidal rule approximation to integral.
	  g0=0.0d0
	    do 20 k=1,nint
	    g0=g0+f(a+(k+k-1)*h)
20	  continue
	  g0=0.5d0*g(1)+h*g0
	  h=0.5d0*h
	  nint=nint+nint
	  jmax=min(i,MAXJ)
	  fourj=1.0d0
	    do 30 j=1,jmax
c  Use Richardson extrapolation.
	    fourj=4.0d0*fourj
	    g1=g0+(g0-g(j))/(fourj-1.0d0)
	    g(j)=g0
	    g0=g1
30	  continue
C	  if (abs(g0).gt.tol) then
	    error=1.0d0-gmax/g0
C	  else
C	    error=gmax
C	  end if
	  gmax=g0
	  g(jmax+1)=g0
	go to 10
40	rombint=g0
	if (i.gt.MAXITER.and.abs(error).gt.1.e-4)
     2   write(*,*) 'Rombint relative error > 1.e-4; integral, error=',
     3   rombint,error
	return
	end
