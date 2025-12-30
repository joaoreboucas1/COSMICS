c  Compare deltal.dat, deltat.dat, linger.con, linger.syn,
c  power.dat, p3m.0a, and p3m.0b with standard case.
	parameter (nparmax=32*32*32)
	real x1(3,nparmax),v1(3,nparmax),x2(3,nparmax)

	open(10,file='deltal.dat',status='old')
	rewind 10
	open(11,file='../test_results/deltal.dat',status='old')
	rewind 11
	write(*,*)
	write(*,*) '---------------------------------------------------'
	write(*,*) 'Checking deltal.dat:'
	write(*,*)
	read(10,*,err=10) l1
	read(11,*,err=10) l2
	if (l1.ne.l2) then
	  write(*,*) 'deltal.dat differs from the standard case:'
	  write(*,*) 'Expected, actual l=',l2,l1
	  go to 10
	end if
	ekmax=0.0
	edmax=0.0
	do k=1,10000
	  read(10,*,end=5) ak1,del1
	  read(11,*,end=5) ak2,del2
	  ekmax=max(ekmax,abs(ak1/ak2)-1.0)
	  edmax=max(edmax,abs(del1-del2))
	end do
5	if (ekmax.gt.1.0e-4) then
	  write(*,*) 'deltal.dat differs from the standard case:'
	  write(*,*) 'Max. rel. error in k=',ekmax
	  go to 10
	end if
	write(*,*) 'Maximum absolute error in Delta_l(k)=',edmax

10	open(10,file='deltat.dat',status='old')
	rewind 10
	open(11,file='../test_results/deltat.dat',status='old')
	rewind 11
	write(*,*)
	write(*,*) '---------------------------------------------------'
	write(*,*) 'Checking deltat.dat:'
	write(*,*)
	read(10,*) omb1,omc1,omv1,omn1
	read(11,*) omb2,omc2,omv2,omn2
	emax=abs(omb1-omb2)+abs(omc1-omc2)+abs(omv1-omv2)
     &       +abs(omn1-omn2)
	if (emax.gt.1.0e-4) then
	  write(*,*) 'deltat.dat header differs from the standard case:'
	  write(*,*) 'Error in omega=',emax
	  go to 20
	end if
	read(10,*) h01,tcmb1,yhe1,nnur1,nnunr1,inf1
	read(11,*) h02,tcmb2,yhe2,nnur2,nnunr2,inf2
	read(10,*) an1
	read(11,*) an2
	emax=abs(h01-h02)+abs(tcmb1-tcmb2)+abs(yhe1-yhe2)
     &      +abs(nnur1-nnur2)+abs(nnunr1-nnunr2)+abs(inf1-inf2)
     &      +abs(an1-an2)
	if (emax.gt.1.0e-3) then
	  write(*,*) 'deltat.dat header differs from the standard case:'
	  write(*,*) 'Error in H0, Tcmb, etc.=',emax
	  go to 20
	end if
	lmax=0.0
	ecmax=0.0
	do k=1,10000
	  read(10,*,end=15) l1,cl1
	  read(11,*,end=15) l2,cl2
	  lmax=max(lmax,abs(l1-l2))
	  ecmax=max(ecmax,abs(cl1-cl2))
	end do
15	if (lmax.ne.0) then
	  write(*,*) 'l-values not in order in deltat.dat!'
	  go to 20
	end if
	write(*,*) 'Maximum absolute error in l(l+1)C_l=',ecmax

20	open(10,file='linger.con',status='old')
	rewind 10
	open(11,file='../test_results/linger.con',status='old')
	rewind 11
	write(*,*)
	write(*,*) '---------------------------------------------------'
	write(*,*) 'Checking linger.con:'
	write(*,*)
	read(10,*) omb1,omc1,omv1,omn1
	read(11,*) omb2,omc2,omv2,omn2
	emax=abs(omb1-omb2)+abs(omc1-omc2)+abs(omv1-omv2)
     &       +abs(omn1-omn2)
	if (emax.gt.1.0e-4) then
	  write(*,*) 'linger.con header differs from the standard case:'
	  write(*,*) 'Error in omega=',emax
	  go to 30
	end if
	read(10,*) h01,tcmb1,yhe1,nnur1,nnunr1,inf1
	read(11,*) h02,tcmb2,yhe2,nnur2,nnunr2,inf2
	emax=abs(h01-h02)+abs(tcmb1-tcmb2)+abs(yhe1-yhe2)
     &      +abs(nnur1-nnur2)+abs(nnunr1-nnunr2)+abs(inf1-inf2)
	if (emax.gt.1.0e-3) then
	  write(*,*) 'linger.con header differs from the standard case:'
	  write(*,*) 'Error in H0, Tcmb, etc.=',emax
	  go to 30
	end if
	err=0
	do k=1,10000
	  read(10,*,end=25) lk1,ak1,a1,t1,ps1,ph1,dc1,db1,dg1,dr1,dn1,
     &     thc1,thb1,thg1,thr1,thn1,shg1,shr1,shn1,ec1
	  read(11,*,end=25) lk2,ak2,a2,t2,ps2,ph2,dc2,db2,dg2,dr2,dn1,
     &     thc2,thb2,thg2,thr2,thn2,shg2,shr2,shn2,ec2
	  if (ik1.ne.ik2) then
	    write(*,*) 'k-values not in order in linger.con!'
	    go to 30
	  end if
	  err=err+abs(ps1-ps2)+abs(ph1-ph2)+abs(dc1-dc2)+abs(db1-db2)
     &     +abs(dg1-dg2)+abs(dr1-dr2)+abs(dn1-dn2)+abs(thc1-thc2)
     &     +abs(thb1-thb2)+abs(thg1-thg2)+abs(thr1-thr2)+abs(thn1-thn2)
     &     +abs(shg1-shg2)+abs(shr1-shr2)+abs(shn1-shn2)+abs(ec1-ec2)
	end do
25	err=err/(15.*k)
	write(*,*) 'Mean absolute error in linger.con=',err

30	open(10,file='linger.syn',status='old')
	rewind 10
	open(11,file='../test_results/linger.syn',status='old')
	rewind 11
	write(*,*)
	write(*,*) '---------------------------------------------------'
	write(*,*) 'Checking linger.syn:'
	write(*,*)
	read(10,*) omb1,omc1,omv1,omn1
	read(11,*) omb2,omc2,omv2,omn2
	emax=abs(omb1-omb2)+abs(omc1-omc2)+abs(omv1-omv2)
     &       +abs(omn1-omn2)
	if (emax.gt.1.0e-4) then
	  write(*,*) 'linger.syn header differs from the standard case:'
	  write(*,*) 'Error in omega=',emax
	  go to 40
	end if
	read(10,*) h01,tcmb1,yhe1,nnur1,nnunr1,inf1
	read(11,*) h02,tcmb2,yhe2,nnur2,nnunr2,inf2
	emax=abs(h01-h02)+abs(tcmb1-tcmb2)+abs(yhe1-yhe2)
     &      +abs(nnur1-nnur2)+abs(nnunr1-nnunr2)+abs(inf1-inf2)
	if (emax.gt.1.0e-3) then
	  write(*,*) 'linger.syn header differs from the standard case:'
	  write(*,*) 'Error in H0, Tcmb, etc.=',emax
	  go to 40
	end if
	err=0
	do k=1,10000
	  read(10,*,end=35) ik1,ak1,a1,t1,ps1,ph1,dc1,db1,dg1,dr1,dn1,
     &     thc1,thb1,thg1,thr1,thn1,shg1,shr1,shn1,ec1
	  read(11,*,end=35) ik2,ak2,a2,t2,ps2,ph2,dc2,db2,dg2,dr2,dn2,
     &     thc2,thb2,thg2,thr2,thn2,shg2,shr2,shn2,ec2
	  if (ik1.ne.ik2) then
	    write(*,*) 'k-values not in order in linger.syn!'
	    go to 40
	  end if
	  err=err+abs(ps1-ps2)+abs(ph1-ph2)+abs(dc1-dc2)+abs(db1-db2)
     &     +abs(dg1-dg2)+abs(dr1-dr2)+abs(dn1-dn2)+abs(thc1-thc2)
     &     +abs(thb1-thb2)+abs(thg1-thg2)+abs(thr1-thr2)+abs(thn1-thn2)
     &     +abs(shg1-shg2)+abs(shr1-shr2)+abs(shn1-shn2)+abs(ec1-ec2)
	end do
35	err=err/(15.*k)
	write(*,*) 'Mean absolute error in linger.syn=',err

40	open(10,file='power.dat',status='old')
	rewind 10
	open(11,file='../test_results/power.dat',status='old')
	rewind 11
	write(*,*)
	write(*,*) '---------------------------------------------------'
	write(*,*) 'Checking power.dat:'
	write(*,*)
	err=0
	do k=1,10000
	  read(10,*,end=45) ak1,p1
	  read(11,*,end=45) ak2,p2
	  err=max(err,abs(p1/p2-1.0))
	end do
45	write(*,*) 'Maximum relative error in P(k)=',err

	stop
	end
