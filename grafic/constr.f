cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine constr(dx,con,g)
c  constr provides the constraint matrix for IC3.f.  con(j1,j2,j3,i)=
c  con((i-1)*n1*n2*n3+j), j=j1+(j2-1)*n1+(j3-1)*n1*n2, is defined so that the
c  constraints are
c		g(i)=sum from j=1 to n1*n2*n3 of con(j,i)*f(j),
c  for i=1 to ncon.  E.g., if n=64 and there are two constraints
c  f(5)=2, f(7)=3,  then con(5,1)=con(5)=2, con(7,2)=con(71)=3, and all
c  other values of con vanish.
c
	include 'grafic.inc'
c
	parameter (n12=np1/2,n22=np2/2,n32=np3/2,n23=np2*np3)
	parameter (ngr2=n23*n12,nmax=ngr2+n23)
	real con(2*nmax,ncon),g(ncon)
	double precision rnum
c
c  Here follows an example: two constraints: a cluster and a void.
c  The constraints are set by the value of delta averaged over a
c  gaussian ball of smoothing radius 20 Mpc.
	grad=20.0/dx
	  do 5 i=1,ncon
	  do 5 j=1,2*nmax
	  con(j,i)=0.0
5	continue
	  do 40 i=1,ncon
	  rnum=0.0d0
	  if (i.eq.1) then
c  Constraint: void centered at X=Lx/4, Y=Ly/2, Z=Lz/4 (Lx,Ly,Lz are the
c  box lengths in the three coordinate directions and X=Y=Z=0 is a corner
c  of the box).
	    x1c=0.25*np1
	    y1c=0.5*np2
	    z1c=0.25*np3
c  N.B. this is the linear density contrast, which may be less than -1.0.
	    delta=-1.0
	  else
c  Constraint: cluster centered in the middle of the box.
	    x1c=0.5*np1
	    y1c=0.5*np2
	    z1c=0.5*np3
	    delta=0.5
	  end if
	  g(i)=0.0
	  j=0
	    do 30 j3=1,np3
	    x3=j3-1
	      do 20 j2=1,np2
	      x2=j2-1
	        do 10 j1=1,np1
	        x1=j1-1
		j=j+1
c  Apply a constraint on a Gaussian sphere of radius grad centered
c  at x1c,y1c,z1c.
	        x1r=x1-x1c
	        x2r=x2-y1c
	        x3r=x3-z1c
	        rad=sqrt(x1r*x1r+x2r*x2r+x3r*x3r)/grad
	        fact=exp(-rad*rad/2.0)
	        con(j,i)=fact
		rnum=rnum+fact
10	      continue
20	    continue
30	  continue
	    do 35 j=1,2*nmax
	    con(j,i)=con(j,i)/rnum
35	  continue
	g(i)=delta
40	continue
	return
	end
