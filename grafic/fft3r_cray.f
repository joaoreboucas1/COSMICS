cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fft3r(a,n1,n2,n3)
c  fft3r performs a 3-D FFT from the spatial domain to the spectral domain.
c  Array a is real of length (n1+2)*n2*n3.  n1 must be divisible by 4;
c  otherwise there are no restrictions on n1,n2,n3.
c  On input, the first n1*n2*n3 elements are filled with values of a(i,j,k)
c  in the spatial domain.  On output, the first n1*n2*n3 elements are filled
c  with values a(i,j,k), i=1,...,n1/2 of the complex Fourier transform.
c  The last 2*n2*n3 values are filled with complex transform values for
c  i=n/2+1.
c  N.B. The transform from spatial to spectral domain is defined with a
c  minus sign in the complex exponential.  To change this, change the sign
c  of twopi.  After transforming and then performing the inverse transform,
c  a(i,j,k) must be divided by n1*n2*n3.
c
	parameter (twopi=-6.283185307179586d0,ndim=3,n0=1024)
	implicit real (a-h,o-z)
	real a(*),wra(n0),wia(n0),wff1(n0),wff2(2*n0),wff3(2*n0)
	complex work(n0*n0)
	integer ifax1(19),ifax2(19),ifax3(19)
	double precision theta,wr,wi
	equivalence (wra(1),wff2(1)),(wia(1),wff2(n0+1))
c
	if (n1.gt.n0.or.n2.gt.n0.or.n3.gt.n0) then
	  write(*,*) 'Need to increase dimension n0 in fft3r!'
	  write(*,*) 'n0,n1,n2,n3=',n0,n1,n2,n3
	  stop
	end if
c
c  Assume contiguous storage, i.e., no strides allowed.
	n1d2=n1/2
	inc2y=n1/2
	inc3y=inc2y*n2
	n3eff=1+(n1-1)+2*((n2-1)*inc2y+(n3-1)*inc3y)
c
c  3-D complex FFT as a series of 1-D complex FFTs.
c  Initialize trignometric arrays.
	call cftfax(n1d2,ifax1,wff1)
	if (ifax1(1).eq.-99) then
	  write(*,*) 'n12=',n1d2,' is not factorizable by cfftmlt!'
	  stop
	end if
	call cftfax(n2,ifax2,wff2)
	if (ifax2(1).eq.-99) then
	  write(*,*) 'n2=',n2,' is not factorizable by cfftmlt!'
	  stop
	end if
	call cftfax(n3,ifax3,wff3)
	if (ifax3(1).eq.-99) then
	  write(*,*) 'n3=',n3,' is not factorizable by cfftmlt!'
	  stop
	end if
c  FFT x-y planes.
	  do 10 i3=1,n3
	  ind=1+2*(i3-1)*inc3y
c  FFT along x-direction.
	  call cfftmlt(a(ind),a(ind+1),work,wff1,ifax1,2,
     2                 2*inc2y,n1d2,n2,-1)
c  FFT along y-direction.
	  call cfftmlt(a(ind),a(ind+1),work,wff2,ifax2,2*inc2y,
     2                 2,n2,n1d2,-1)
10	continue
c  FFT in x-z plane.
	  do 20 i2=1,n2
	  ind=1+2*(i2-1)*inc2y
c  FFT along z-direction.
	  call cfftmlt(a(ind),a(ind+1),work,wff3,ifax3,2*inc3y,
     2                 2,n3,n1d2,-1)
20	continue
c
c  Use symmetries to obtain the true transform.  This part is specially
c  written for 3-D.
	theta=twopi/n1
	  do 24 k1=1,n1
	  wra(k1)=cos(theta*(k1-1))
	  wia(k1)=sin(theta*(k1-1))
24	continue
c  Save k1=1,n14+1 for later.
	n14=n1d2/2
	  do 30 k3=1,n3
	  k3c=n3+2-k3
	  if (k3.eq.1) k3c=1
c  Swap k1 and k2 loops so that the inner loop is longer (for vectorization).
	    do 28 k1=2,n14
	    k1c=n1d2+2-k1
	    wr=wra(k1)
	    wi=wia(k1)
CDIR$ IVDEP
	      do 26 k2=1,n2
	      k2c=n2+2-k2
	      if (k2.eq.1) k2c=1
	      ind=1+2*((k1-1)+(k2-1)*inc2y+(k3-1)*inc3y)
	      indc=1+2*((k1c-1)+(k2c-1)*inc2y+(k3c-1)*inc3y)
	      ffer=0.5*(a(ind)+a(indc))
	      ffei=0.5*(a(ind+1)-a(indc+1))
	      ffor=0.5*(a(ind+1)+a(indc+1))
	      ffoi=-0.5*(a(ind)-a(indc))
	      tempr=wr*ffor-wi*ffoi
	      tempi=wi*ffor+wr*ffoi
	      a(ind)=ffer+tempr
	      a(ind+1)=ffei+tempi
	      a(indc)=ffer-tempr
	      a(indc+1)=-ffei+tempi
26	    continue
28	  continue
30	continue
	  do 34 k3=1,n3
	  k3c=n3+2-k3
	  if (k3.eq.1) k3c=1
	    do 32 k2=1,n2
	    k2c=n2+2-k2
	    if (k2.eq.1) k2c=1
c  k1=n14.
	    k1=n14+1
	    k1c=n1d2+2-k1
	    ind=1+2*((k1-1)+(k2-1)*inc2y+(k3-1)*inc3y)
	    indc=1+2*((k1c-1)+(k2c-1)*inc2y+(k3c-1)*inc3y)
	    if (ind.gt.indc) go to 31
	    ffer=0.5*(a(ind)+a(indc))
	    ffei=0.5*(a(ind+1)-a(indc+1))
	    ffor=0.5*(a(ind+1)+a(indc+1))
	    ffoi=-0.5*(a(ind)-a(indc))
	    wr=wra(k1)
	    wi=wia(k1)
	    tempr=wr*ffor-wi*ffoi
	    tempi=wi*ffor+wr*ffoi
	    a(ind)=ffer+tempr
	    a(ind+1)=ffei+tempi
	    a(indc)=ffer-tempr
	    a(indc+1)=-ffei+tempi
c  k1=1.
31	    k1=1
	    ind=1+2*((k2-1)*inc2y+(k3-1)*inc3y)
	    indc=1+2*((k2c-1)*inc2y+(k3c-1)*inc3y)
	    if (ind.gt.indc) go to 32
	    k23=(k2-1)+(k3-1)*n2
	    k23c=(k2c-1)+(k3c-1)*n2
	    i3=n3eff+1+2*k23
	    i4=n3eff+1+2*k23c
	    ffer=0.5*(a(ind)+a(indc))
	    ffei=0.5*(a(ind+1)-a(indc+1))
	    ffor=0.5*(a(ind+1)+a(indc+1))
	    ffoi=-0.5*(a(ind)-a(indc))
	    a(ind)=ffer+ffor
	    a(ind+1)=ffei+ffoi
	    a(indc)=a(ind)
	    a(indc+1)=-a(ind+1)
	    a(i4)=ffer-ffor
	    a(i4+1)=-ffei+ffoi
	    a(i3)=a(i4)
	    a(i3+1)=-a(i4+1)
32	 continue
34	continue
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fft3rinv(a,n1,n2,n3)
c  fft3rinv performs a 3-D inverse FFT from the spectral domain to the spatial
c  domain.  Array a is real of length (n1+2)*n2*n3.  n1 must be divisible by 4;
c  otherwise there are no restrictions on n1,n2,n3.
c  On input, the first n1*n2*n3 elements are filled with values a(i,j,k),
c  i=1,...,n1/2 of the complex Fourier transform.  The last 2*n2*n3 values
c  are filled with complex transform values for i=n1/2+1.  On output, the
c  first n1*n2*n3 elements are filled with values of a(i,j,k) in the spatial
c  domain.
c  N.B. The transform from spectral to spatial domain is defined with a
c  plus sign in the complex exponential.  To change this, change the sign
c  of twopi.  After transforming and then performing the inverse transform,
c  a(i,j,k) must be divided by n1*n2*n3.
c
	parameter (twopi=6.283185307179586d0,ndim=3,n0=1024)
	implicit real (a-h,o-z)
	real a(*),wra(n0),wia(n0),wff1(n0),wff2(2*n0),wff3(2*n0)
	complex work(n0*n0)
	integer ifax1(19),ifax2(19),ifax3(19)
	double precision theta,wr,wi
	equivalence (wra(1),wff2(1)),(wia(1),wff2(n0+1))
c
	if (n1.gt.n0.or.n2.gt.n0.or.n3.gt.n0) then
	  write(*,*) 'Need to increase dimension n0 in fft3rinv!'
	  write(*,*) 'n0,n1,n2,n3=',n0,n1,n2,n3
	  stop
	end if
c
c  Assume contiguous storage, i.e., no strides allowed.
	n1d2=n1/2
	inc2y=n1/2
	inc3y=inc2y*n2
	n3eff=1+(n1-1)+2*((n2-1)*inc2y+(n3-1)*inc3y)
c
c  Use symmetries to obtain the true transform.  This part is specially
c  written for 3-D.
	theta=twopi/n1
	  do 24 k1=1,n1
	  wra(k1)=cos(theta*(k1-1))
	  wia(k1)=sin(theta*(k1-1))
24	continue
c  Save k1=1,n14+1 for later.
	n14=n1d2/2
	  do 30 k3=1,n3
	  k3c=n3+2-k3
	  if (k3.eq.1) k3c=1
c  Swap k1 and k2 loops so that the inner loop is longer (for vectorization).
	    do 28 k1=2,n14
	    k1c=n1d2+2-k1
	    wr=wra(k1)
	    wi=wia(k1)
CDIR$ IVDEP
	      do 26 k2=1,n2
	      k2c=n2+2-k2
	      if (k2.eq.1) k2c=1
	      ind=1+2*((k1-1)+(k2-1)*inc2y+(k3-1)*inc3y)
	      indc=1+2*((k1c-1)+(k2c-1)*inc2y+(k3c-1)*inc3y)
c  Multiply by 2 since this is a real FFT summing over only half
c  of the harmonics.
	      ffer=a(ind)+a(indc)
	      ffei=a(ind+1)-a(indc+1)
	      ffor=-(a(ind+1)+a(indc+1))
	      ffoi=a(ind)-a(indc)
	      tempr=wr*ffor-wi*ffoi
	      tempi=wi*ffor+wr*ffoi
	      a(ind)=ffer+tempr
	      a(ind+1)=ffei+tempi
	      a(indc)=ffer-tempr
	      a(indc+1)=-ffei+tempi
26	    continue
28	  continue
30	continue
	  do 34 k3=1,n3
	  k3c=n3+2-k3
	  if (k3.eq.1) k3c=1
	    do 32 k2=1,n2
	    k2c=n2+2-k2
	    if (k2.eq.1) k2c=1
c  k1=n14.
	    k1=n14+1
	    k1c=n1d2+2-k1
	    ind=1+2*((k1-1)+(k2-1)*inc2y+(k3-1)*inc3y)
	    indc=1+2*((k1c-1)+(k2c-1)*inc2y+(k3c-1)*inc3y)
	    if (ind.gt.indc) go to 31
	    ffer=a(ind)+a(indc)
	    ffei=a(ind+1)-a(indc+1)
	    ffor=-(a(ind+1)+a(indc+1))
	    ffoi=a(ind)-a(indc)
	    wr=wra(k1)
	    wi=wia(k1)
	    tempr=wr*ffor-wi*ffoi
	    tempi=wi*ffor+wr*ffoi
	    a(ind)=ffer+tempr
	    a(ind+1)=ffei+tempi
	    a(indc)=ffer-tempr
	    a(indc+1)=-ffei+tempi
c  k1=1.
31	    k1=1
	    ind=1+2*((k2-1)*inc2y+(k3-1)*inc3y)
	    indc=1+2*((k2c-1)*inc2y+(k3c-1)*inc3y)
	    if (ind.gt.indc) go to 32
	    k23=(k2-1)+(k3-1)*n2
	    k23c=(k2c-1)+(k3c-1)*n2
	    i4=n3eff+1+2*k23c
	    ffer=a(ind)+a(i4)
	    ffei=a(ind+1)-a(i4+1)
	    ffor=-(a(ind+1)+a(i4+1))
	    ffoi=a(ind)-a(i4)
	    a(ind)=ffer+ffor
	    a(ind+1)=ffei+ffoi
	    a(indc)=ffer-ffor
	    a(indc+1)=-ffei+ffoi
32	 continue
34	continue
c
c  3-D complex FFT as a series of 1-D complex FFTs.
c  Initialize trignometric arrays.
	call cftfax(n1d2,ifax1,wff1)
	if (ifax1(1).eq.-99) then
	  write(*,*) 'n12=',n1d2,' is not factorizable by cfftmlt!'
	  stop
	end if
	call cftfax(n2,ifax2,wff2)
	if (ifax2(1).eq.-99) then
	  write(*,*) 'n2=',n2,' is not factorizable by cfftmlt!'
	  stop
	end if
	call cftfax(n3,ifax3,wff3)
	if (ifax3(1).eq.-99) then
	  write(*,*) 'n3=',n3,' is not factorizable by cfftmlt!'
	  stop
	end if
c  FFT x-y planes.
	  do 10 i3=1,n3
	  ind=1+2*(i3-1)*inc3y
c  FFT along x-direction.
	  call cfftmlt(a(ind),a(ind+1),work,wff1,ifax1,2,
     2                 2*inc2y,n1d2,n2,+1)
c  FFT along y-direction.
	  call cfftmlt(a(ind),a(ind+1),work,wff2,ifax2,2*inc2y,
     2                 2,n2,n1d2,+1)
10	continue
c  FFT in x-z plane.
	  do 20 i2=1,n2
	  ind=1+2*(i2-1)*inc2y
c  FFT along z-direction.
	  call cfftmlt(a(ind),a(ind+1),work,wff3,ifax3,2*inc3y,
     2                 2,n3,n1d2,+1)
20	continue
c
	return
	end
