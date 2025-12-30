cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ic3(ido,dx,rho,psi,delm,dispm)
c  Generate a constrained (ido=3) or unconstrained (ido=1) gaussian
c  random field rho (or the mean constrained field if ido=2) using
c  the power spectrum p(ak) at a=1, and the corresponding displacement
c  field psi in Mpc.
c
c  Note: i,i1,i2 label constraints; j labels spatial grid points; k,k1,k2,k3
c  label grid points in Fourier space.  The rectangular grid of size
c  (np1,np2,np3) has even spacings dx=dy=dz.  np1 should be divisible by 4
c  while np2 and np3 should be divisible by 2.  In Fourier space there are
c  (np1/2+1)*np2*np3 grid points (k1 < 0 is redundant because rho is real).
c
	include 'grafic.inc'
c
	parameter (n12=np1/2,n22=np2/2,n32=np3/2,n23=np2*np3)
	parameter (ngr2=n23*n12,nmax=ngr2+n23)
	real rho(npart),psi(3,npart)
	double precision twopi,sigma,sigmav,chisq
	parameter (twopi=6.283185307179586d0,sqr2=1.41421356)
	real power(nmax),con(2*nmax,ncon)
	real g(ncon),g0(ncon)
	double precision q(ncon,ncon),qinv(ncon,ncon)
	complex rhoc(nmax),conc(nmax,ncon),z
	complex psi1(nmax),psi2(nmax),psi3(nmax)
	real p
	equivalence (con,conc)
c
	dk1=twopi/(np1*dx)
	dk2=twopi/(np2*dx)
	dk3=twopi/(np3*dx)
	d3k=dk1*dk2*dk3
	akmax=twopi/dx
	sigma=0.0
	sigmav=0.0
	chisq=0.0
	nharm=0
c  Generate unconstrained sample of rho in Fourier transform space.
	  do 30 k3=1,np3
	  ak3=(k3-1)*dk3
	  if (k3.gt.n32) ak3=ak3-akmax
	  ak33=ak3*ak3
	    do 20 k2=1,np2
	    ak2=(k2-1)*dk2
	    if (k2.gt.n22) ak2=ak2-akmax
	    ak23=ak2*ak2+ak33
	    k23=(k2-1+(k3-1)*np2)*n12
	      do 10 k1=2,n12
c  Do k1=1 and k1=n12+1 separately below.
	      ak1=(k1-1)*dk1
	      k=k1+k23
	      akk=ak1*ak1+ak23
	      ak=sqrt(akk)
c  Mean power per mode (k1,k2,k3).
	      power(k)=p(ak)*d3k
	      sigft=sqrt(power(k))
	      call randa(xr)
	      arg=twopi*xr
	      z=cmplx(cos(arg),sin(arg))
5	      call randa(xr)
	      if (xr.eq.0.0) go to 5
c  Reduce z by sqrt(2) for real f since need sum of |z|^2 over all
c  wavenumbers (positive and negative) to be Chi-square for N degrees
c  of freedom, not 2N as in the complex case.  Equivalently, the power
c  is shared with the conjugate harmonics with k1 > n12+1 (k1 < 0).
	      z=sqrt(-log(xr))*z
	      rhoc(k)=z*sigft
c  Double the contribution to account for modes with k1 > n12+1 (k1 < 0).
	      sigma=sigma+2.0*power(k)
	      sigmav=sigmav+2.0*power(k)/akk
	      chisq=chisq+2.0*conjg(z)*z
	      nharm=nharm+2
10	    continue
	    k23=k2-1+(k3-1)*np2
	    k23c=np2+1-k2+(np3+1-k3)*np2
	    if (k2.eq.1) k23c=k23c-np2
	    if (k3.eq.1) k23c=k23c-n23
	    if (k23.gt.k23c) go to 20
	      do 15 i=1,2
c  Do k1=1 and k1=n12+1.
	      if (i.eq.1) then
		k=1+n12*k23
		kc=1+n12*k23c
		ak1=0.0
	      else
		k=ngr2+1+k23
		kc=ngr2+1+k23c
		ak1=-0.5*akmax
	      end if
	      akk=ak1*ak1+ak23
	      ak=sqrt(akk)
	      power(k)=p(ak)*d3k
	      power(kc)=power(k)
	      sigft=sqrt(power(k))
	      call randa(xr)
	      arg=twopi*xr
	      if (k23.ne.k23c) then
	        z=cmplx(cos(arg),sin(arg))
	      else
	        z=sqr2*cmplx(cos(arg),0.0)
	      end if
12	      call randa(xr)
	      if (xr.eq.0.0) go to 12
	      z=sqrt(-log(xr))*z
	      rhoc(k)=z*sigft
c  Determine conjugate harmonic using symmetry.
	      rhoc(kc)=conjg(rhoc(k))
	      if (k23.ne.k23c) then
	        sigma=sigma+2.0*power(k)
	        sigmav=sigmav+2.0*power(k)/akk
	        chisq=chisq+2.0*conjg(z)*z
	        nharm=nharm+2
	      else if (akk.ne.0.0) then
	        sigma=sigma+power(k)
	        chisq=chisq+conjg(z)*z
	        nharm=nharm+1
	      end if
15	    continue
20	  continue
30	continue
	rhoc(1)=cmplx(0.0,0.0)
	sigma=sqrt(sigma)
	sigmav=sqrt(sigmav)
	write(*,*)
	write(*,*) 'For unconstrained field at a=1:'
	write(*,*) '     Mean sigma_delta, sigma_psi=',real(sigma),
     2    real(sigmav),' Mpc'
	anu=(chisq-nharm)/sqrt(2.0*nharm)
	write(*,*) '     Chisq, dof, nu=',real(chisq),nharm,anu
	if (ido.eq.1) go to 165
c
c  Get constraints.
c  Constraints are of the form g(i)=sum from j=1 to 2*ngr2 of con(j,i)*rho(j).
	call constr(dx,con,g)
c  Fourier transform the constraints so that the convolutions with the
c  covariance matrix required to compute the mean field and the covariance
c  of the constraints becomes multiplication by the power spectrum.
	write(*,*)
	  do 50 i=1,ncon
	  call fft3r(conc(1,i),np1,np2,np3)
c  Compute the values of the constraints for the unconstrained realization.
	  sigma=0.0
	    do 40 k=1,nmax
c  Take twice the real part for harmonics with k1 > 0 to account for the
c  conjugate harmonics with k1 < 0.
	    sigma=sigma+2.0*conjg(conc(k,i))*rhoc(k)
40	  continue
c  Correct for k1=1 and k1=n12+1, which are self-conjugate.
	    do 46 k3=1,np3
	      do 43 k2=1,np2
	      k23=k2-1+(k3-1)*np2
	      k=1+n12*k23
	      sigma=sigma-conjg(conc(k,i))*rhoc(k)
	      k=ngr2+1+k23
	      sigma=sigma-conjg(conc(k,i))*rhoc(k)
43	    continue
46	  continue
	  g0(i)=sigma
	  write(*,1001) 'Constraint ',i,':  Sampled, desired=',g0(i),g(i)
50	continue
1001	format(1x,a,i3,a,2g15.8)
c  Compute covariance matrix of constraints.
	  do 80 i2=1,ncon
c	  write(*,*) 'i2=',i2
	    do 70 i1=1,ncon
	    q(i1,i2)=0.0d0
	      do 60 k=1,nmax
	      q(i1,i2)=q(i1,i2)+2.0*conjg(conc(k,i1))*conc(k,i2)*power(k)
60	    continue
c  Correct q for k1=1 and k1=n12+1.
	    do 66 k3=1,np3
	      do 63 k2=1,np2
	      k23=k2-1+(k3-1)*np2
	      k=1+n12*k23
	      q(i1,i2)=q(i1,i2)-conjg(conc(k,i1))*conc(k,i2)*power(k)
	      k=ngr2+1+k23
	      q(i1,i2)=q(i1,i2)-conjg(conc(k,i1))*conc(k,i2)*power(k)
63	    continue
66	  continue
70	  continue
80	continue
c  Invert double precision matrix q to give qinv.  Warning: q is destroyed!
	call matinv(q,ncon,ncon,qinv)
c  Form chi-squared for the constraints.
	chisq=0.0
	chisq0=0.0
	  do 100 i2=1,ncon
	    do 90 i1=1,ncon
	    chisq=chisq+g(i1)*qinv(i1,i2)*g(i2)
	    chisq0=chisq0+g0(i1)*qinv(i1,i2)*g0(i2)
90	  continue
100	continue
	write(*,*)
	write(*,1001) 'Chi-square for the ',ncon,' constraints:'
	write(*,*) '     Sampled, desired=',chisq0,real(chisq)
	chisq0=chisq
	write(*,*)
	if (ido.eq.2) then
	  write(*,*) 'Computing mean field only'
	    do 105 k=1,nmax
105	  rhoc(k)=cmplx(0.0,0.0)
	    do 106 i=1,ncon
106	  g0(i)=0.0
	else
	  write(*,*) 'Computing a constrained realization'
	end if
c  Compute corrected density field.
	  do 130 i2=1,ncon
	    do 120 i1=1,ncon
	    fact=(g(i2)-g0(i2))*qinv(i1,i2)
	      do 110 k=1,nmax
	      rhoc(k)=rhoc(k)+fact*conc(k,i1)*power(k)
110	    continue
120	  continue
130	continue
	  do 160 i=1,ncon
c  Check the values of the constraints for the constrained realization.
	  sigma=0.0
	    do 150 k=1,nmax
c  Take twice the real part for harmonics with k1 > 0 to account for the
c  conjugate harmonics with k1 < 0.
	    sigma=sigma+2.0*conjg(conc(k,i))*rhoc(k)
150	  continue
c  Correct for k1=1 and k1=n12+1, which are self-conjugate.
	    do 156 k3=1,np3
	      do 153 k2=1,np2
	      k23=k2-1+(k3-1)*np2
	      k=1+n12*k23
	      sigma=sigma-conjg(conc(k,i))*rhoc(k)
	      k=ngr2+1+k23
	      sigma=sigma-conjg(conc(k,i))*rhoc(k)
153	    continue
156	  continue
	  write(*,1001) 'Constraint ',i,':  Final=',real(sigma)
160	continue
c  Compute displacement field.
165	continue
	sigma=0.0
	sigmav=0.0
	chisq=0.0
	nharm=0
	psi1(1)=cmplx(0.0,0.0)
	psi2(1)=cmplx(0.0,0.0)
	psi3(1)=cmplx(0.0,0.0)
	  do 190 k3=1,np3
	  ak3=(k3-1)*dk3
	  if (k3.gt.n32) ak3=ak3-akmax
	  ak33=ak3*ak3
	    do 180 k2=1,np2
	    ak2=(k2-1)*dk2
	    if (k2.gt.n22) ak2=ak2-akmax
	    ak23=ak2*ak2+ak33
	    k23=(k2-1+(k3-1)*np2)*n12
	      do 170 k1=2,n12
c  Do k1=1 and k1=n12+1 separately below.
	      ak1=(k1-1)*dk1
	      k=k1+k23
	      akk=ak1*ak1+ak23
	      z=cmplx(0.0,1.0)*rhoc(k)/akk
	      psi1(k)=ak1*z
	      psi2(k)=ak2*z
	      psi3(k)=ak3*z
	      if (k2.eq.n22+1) psi2(k)=cmplx(0.0,0.0)
	      if (k3.eq.n32+1) psi3(k)=cmplx(0.0,0.0)
	      dsig=2.0*conjg(rhoc(k))*rhoc(k)
	      sigma=sigma+dsig
	      dsigv=conjg(psi1(k))*psi1(k)+conjg(psi2(k))*psi2(k)
     2             +conjg(psi3(k))*psi3(k)
	      sigmav=sigmav+2.0*dsigv
	      chisq=chisq+dsig/power(k)
	      nharm=nharm+2
170	    continue
	    k23=k2-1+(k3-1)*np2
	    k23c=np2+1-k2+(np3+1-k3)*np2
	    if (k2.eq.1) k23c=k23c-np2
	    if (k3.eq.1) k23c=k23c-n23
	    if (k23.gt.k23c) go to 180
	      do 175 i=1,2
c  Do k1=1 and k1=n12+1.
	      if (i.eq.1) then
		k=1+n12*k23
		kc=1+n12*k23c
		ak1=0.0
	      else
		k=ngr2+1+k23
		kc=ngr2+1+k23c
		ak1=-0.5*akmax
	      end if
	      akk=ak1*ak1+ak23
	      if (k23.ne.k23c) then
	        z=cmplx(0.0,1.0)*rhoc(k)/akk
	        dsig=2.0*conjg(rhoc(k))*rhoc(k)
	        sigma=sigma+dsig
		chisq=chisq+dsig/power(k)
		nharm=nharm+2
	      else if (akk.ne.0.0) then
	        z=cmplx(0.0,0.0)
		dsig=conjg(rhoc(k))*rhoc(k)
	        sigma=sigma+dsig
		chisq=chisq+dsig/power(k)
		nharm=nharm+1
	      end if
	      psi1(k)=cmplx(0.0,0.0)
	      psi2(k)=ak2*z
	      psi3(k)=ak3*z
	      if (k2.eq.n22+1) psi2(k)=cmplx(0.0,0.0)
	      if (k3.eq.n32+1) psi3(k)=cmplx(0.0,0.0)
c  Determine conjugate harmonic using symmetry.
	      psi1(kc)=conjg(psi1(k))
	      psi2(kc)=conjg(psi2(k))
	      psi3(kc)=conjg(psi3(k))
	      dsigv=psi1(kc)*psi1(k)+psi2(kc)*psi2(k)
     2             +psi3(kc)*psi3(k)
	      sigmav=sigmav+2.0*dsigv
175	    continue
180	  continue
190	continue
	sigma=sqrt(sigma)
	sigmav=sqrt(sigmav)
	write(*,*)
	write(*,*) 'For random realization at a=1:'
	write(*,*) '     sigma_delta, sigma_psi=',real(sigma),
     2    real(sigmav),' Mpc'
	if (ido.eq.1) then
	  ndof=nharm
	else
	  ndof=nharm-ncon
	end if
	write(*,*) '     Chisq, dof=',real(chisq-chisq0),ndof
c  Transform to position space.
	call fft3rinv(rhoc,np1,np2,np3)
	call fft3rinv(psi1,np1,np2,np3)
	call fft3rinv(psi2,np1,np2,np3)
	call fft3rinv(psi3,np1,np2,np3)
	delm=0.0
	dispm=0.0
	ii=0
	  do 200 j=1,ngr2
	  del=abs(real(rhoc(j)))
	  delm=max(delm,del)
	  del=abs(aimag(rhoc(j)))
	  delm=max(delm,del)
	  ddisp=real(psi1(j))**2+real(psi2(j))**2+real(psi3(j))**2
	  dispm=max(dispm,ddisp)
	  ddisp=aimag(psi1(j))**2+aimag(psi2(j))**2+aimag(psi3(j))**2
	  dispm=max(dispm,ddisp)
	  ii=ii+1
	  rho(ii)=real(rhoc(j))
	  psi(1,ii)=real(psi1(j))
	  psi(2,ii)=real(psi2(j))
	  psi(3,ii)=real(psi3(j))
	  ii=ii+1
	  rho(ii)=aimag(rhoc(j))
	  psi(1,ii)=aimag(psi1(j))
	  psi(2,ii)=aimag(psi2(j))
	  psi(3,ii)=aimag(psi3(j))
200	continue
	dispm=sqrt(dispm)
	write(*,*) '     Maximum delta, displacement=',delm,dispm,' Mpc'
c
	return
	end
