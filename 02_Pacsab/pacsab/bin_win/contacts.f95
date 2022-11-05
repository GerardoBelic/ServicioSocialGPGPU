!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C PROGRAMA DMDOOCKINGHB
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C RESOLUCIO ATOMICA A L'INTERFICIE I NOMES A NIVELL RESIDU A LA RESTA.
!C
!C FA SERVIR PERFILS DE POTENCIAL AMB NSTEP ESGLAONS.
!C
!C POT CREAR I DESTRUIR PONTS D'HIDROGEN ENTRE RECEPTOR I LLIGAND. LA RESTA FIXES.
!C
!C DISTANCIES I VELOCITATS EN ANGSTROMS
!C TEMPERATURA EN K
!C ENERGIES EN kcal/mol
!C CONSIDERA TOTS ELS ATOMS EXCEPTE ELS H
!C LLEGEIX LA TOPOLOGIA DE L'ARXIU file7, QUE HA SE SER UNA TAULA DE DUES COLUMNES
!C AMB TOTES LES PARELLES D'ATOMS UNITS PER BONDS O PSEUDOBONDS
!C CONSIDERA INTERACCIONS HIDROFOBIQUES A TRAVES DELS 
!C POTENCIALS DE SOLVATACIO DE LAZARIDIS-KARPLUS
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	program dmdockinghb
	implicit real*8(a-h,o-z)
	real*4 rstep,estep,xdir
	integer*2 nstep,ixoc,icut,istruct
	parameter(natmax=4000,facte=4167.)
	parameter(nmax=10000)
	character*5 atom,res
	common/xoc/r(natmax,3),v(natmax,3),xm(natmax),ierr
	common/pous/rstep(natmax,natmax,5),estep(natmax,natmax,5)
	common/inter/nstep(natmax,natmax),istruct(natmax,natmax)
	common/cov/rb(natmax,natmax),drb(natmax,natmax),rhc(natmax)
	common/pdb/atom(nmax),res(nmax),ind2(nmax)
	common/fisic/evdw(natmax),rvdw(natmax),qq(natmax),gfree(natmax),gfrpy(natmax),vol(natmax),xlamb(natmax)
	common/param/fvdw,fsolv,eps,hps,dpsext,dpsint,dhf,dcut,disrad,isolv
	common/moment/deltak(natmax,natmax),timp(natmax,natmax)
	dimension rexp(natmax,3),rant(natmax,3),r0(nmax,3),igo(natmax),irs(natmax),ipart(natmax),tpart(natmax),itg(natmax),isec(natmax)
	dimension ica(natmax),ico(natmax),in(natmax),io(natmax),ind1(natmax),idock(natmax),iid(nmax),ihb(natmax)
	dimension rcm(3),vcm(3),vd(3),vm1(3),vm2(3),argk(2)
	dimension ethf(13,13)
	dimension rb0(nmax,nmax)
	dimension inb(natmax),nblist(natmax,natmax),ncnat(natmax,natmax),nct(natmax,natmax)
	character*100 file7,file8,file9,file10,file11,file12,file13,file16,file17,file18,file19,file20
	character*5 c1,c2,c3,c4,c5,tipus(10),atp(natmax),cad(natmax)
	namelist /input/ tsnap,tcorr,tact,sigma,temp,kkk,file7,file8,file9,file10,file12,file16,file17,file18,file19,file20,nbloc,rcutgo,sigmago,tmin,ego,ehb,coshb,dhf,isolv,fsolv,fvdw,eps,hps,dpsint,dpsext,dcut,irig,irigd,dr,rfixa,rdock,natom1,rcutoff
	namelist /distancies/rnomax,rnomin,rncmax,rncmin,rcomax,rcomin
 2000	FORMAT (A4,2X,I5,2X,A3,1X,A3,1X,A1,1X,I3,5X,F7.2,1X,F7.2,1X,F7.2)
 1000	FORMAT (A4,2X,I5,2X,A3,1X,A3,3X,I3,5X,F7.2,1X,F7.2,1X,F7.2)
 333	FORMAT('MODEL',8X,I4)
 334	FORMAT('ENDMDL')
 1015	FORMAT(E12.6,20(1X,F7.2))
 111	FORMAT('TER')


!C default values
	file9='input.pdb'
	file11='distancia.dat'
	file12='energia.dat'
	file19='patro.pdb'
	file20='trajectory.pdb'
	file7='topologia.dat'
	file16='atomtypes.dat'
	file17='potentials.dat'
	file18='potpydock.dat'
	file8='res'
	file10='complex.pdb'
	pi=atan(1.d0)*4.0d0
	a=1.d-10
	tmin=1.d-22
	kkk=2381
	temp=300.
	nbloc=2
	tact=5.d-14
	tcorr=1.d-13
	tsnap=1.d-11
	rcutgo=10.
	rcutoff=5.
	rco=1.4
	rca=1.6
	rcb=1.6
	rs=1.5
	rn=1.2
	ro=1.2
	rvdwc=2.
	rvdwo=1.4
	rvdwn=1.4
	rvdws=1.7
	xmc=.012
	xmo=.016
	xmn=.014
	xms=.032
	xmres=.1
	sigma=0.01
	sigmago=0.1
!C interacciones fisiques
	eps=1.
	hps=0.4
	dpsext=7.0
	dpsint=5.0
	dhf=5.0
	dcut=10.
	fsolv=.5
	fvdw=.6
	disrad=0.
	ego=1.
	ehb=1.
	coshb=0.8
!C docking
	rfixa=12.
	rdock=8.
	dr=1.
	irig=0
	irigd=0

	isolv=1

	rnomax=4.1
	rnomin=2.7
	rncmax=5.
	rncmin=3.6
	rcomax=5.
	rcomin=3.6

	tipus(1)='CA'
	natom1=0

	read(5,INPUT)
	read(5,DISTANCIES)


	if(tcorr.gt.tsnap)tcorr=tsnap
	if(tact.gt.tcorr)tact=tcorr

!C	file16='/home/agusti/dat/atomtypes.dat'
!C	file17='/home/agusti/dat/potpydock.dat'

	dene=0.

!C converteixo la temperatura en energia per particula
!C utilitzo la constant de Boltzmann per mol
	temp0=1.5*8.32*temp
	tempcal=temp/facte

!C llegeix les coordenades (fitxer pdb)
	open(unit=9,file=file9)
	n=1
 50	read(9,*,end=51)c1,j,c2,c3,c4,k,x,y,z
	if(c2.eq.'OXT'.and.natom1.eq.0)natom1=n
	atom(n)=c2
	res(n)=c3
	cad(n)=c4
	ind2(n)=k
	r0(n,1)=x
	r0(n,2)=y
	r0(n,3)=z
	if(atom(n).eq.'N')in(k)=n
	if(atom(n).eq.'CA')ica(k)=n
	if(atom(n).eq.'C')ico(k)=n
	if(atom(n).eq.'O')io(k)=n
	n=n+1
	goto 50
 51	natom=n-1

	close(9)
	nres=ind2(natom)
	nres1=ind2(natom1)

	open(unit=10,file=file10)
	do n=1,natom
	read(10,*)c1,j,c2,c3,c4,k,x0,y0,z0
	rexp(n,1)=x0
	rexp(n,2)=y0
	rexp(n,3)=z0
	enddo
	close(10)

	do i=1,nres1
	do j=nres1+1,nres
	ncnat(i,j)=0
	enddo
	enddo

	do i=1,natom1
	ii=ind2(i)
	do j=natom1+1,natom
	jj=ind2(j)
	rij1=rexp(i,1)-rexp(j,1)
	rij2=rexp(i,2)-rexp(j,2)
	rij3=rexp(i,3)-rexp(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(rij.lt.rcutoff)ncnat(ii,jj)=1
	enddo
	enddo

	ncontnat=0
	do i=1,nres1
	do j=nres1+1,nres
	if(ncnat(i,j).eq.1)ncontnat=ncontnat+1
	enddo
	enddo





	do i=1,nres
	idock(i)=0
	enddo

!C reconeix els residus aprop de l'interficie per establir els lligams que estabilitzen el complex
	do i=1,natom1
	ii=ind2(i)
	do j=natom1+1,natom
	jj=ind2(j)
	rij1=r0(i,1)-r0(j,1)
	rij2=r0(i,2)-r0(j,2)
	rij3=r0(i,3)-r0(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	 if(rij.lt.rfixa)then
	 idock(ii)=1
	 idock(jj)=1
	 endif
	enddo
	enddo

!C reconeix els residus de l'interficie
	do i=1,natom1
	ii=ind2(i)
	do j=natom1+1,natom
	jj=ind2(j)
	rij1=r0(i,1)-r0(j,1)
	rij2=r0(i,2)-r0(j,2)
	rij3=r0(i,3)-r0(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	 if(rij.lt.rdock)then
	 idock(ii)=2
	 idock(jj)=2
	 endif
	enddo
	enddo

!C es reescriu l'estructura
	i=0
	do n=1,nres
	ii=in(n)
	if(idock(n).lt.1)then
!C	i=i+1
!C	iid(i)=in(n)
	i=i+1
	iid(i)=ica(n)
!C	i=i+1
!C	iid(i)=ico(n)
!C	i=i+1
!C	iid(i)=io(n)
	else
	k=0
 234	i=i+1
	iid(i)=ii+k
	k=k+1
	if(ind2(ii+k).eq.n)goto 234
	endif
	enddo

	natred=i


	do i=1,nres
	in(i)=0
	ica(i)=0
	ico(i)=0
	io(i)=0
	enddo

	iin=0
	iica=0
	iico=0
	iio=0

	do i=1,natred
	ii=iid(i)
	r(i,1)=r0(ii,1)
	r(i,2)=r0(ii,2)
	r(i,3)=r0(ii,3)
!C ii mes gran o igual que i sempre
	c2=atom(ii)
	atom(i)=c2
	c2=res(ii)
	res(i)=c2
	k=ind2(ii)
	ind2(i)=k
	k=ind2(i)
	if(atom(i).eq.'N')then
	in(k)=i
	iin=iin+1
	endif
	if(atom(i).eq.'CA')then
	ica(k)=i
	iica=iica+1
	endif
	if(atom(i).eq.'C')then
	ico(k)=i
	iico=iico+1
	endif
	if(atom(i).eq.'O')then
	io(k)=i
	iio=iio+1
	endif
	enddo

	if(idock(nres1+1).lt.1)then
	natred1=ica(nres1+1)-1
	elseif(idock(nres1).lt.1)then
	natred1=ica(nres1)
	else
	  n1=ica(nres1+1)
	  if(res(n1).eq.'PRO')then
	  natred1=n1-5
	  else
	  natred1=n1-2
	  endif
	endif



!Ccccccccccccccccccccccccccccccccccccccc
!C llegeix la matriu de topologia
!Ccccccccccccccccccccccccccccccccccccccc
	do i=1,natom-1
	do j=i+1,natom
	rb0(i,j)=0.
	enddo
	enddo

	open(unit=7,file=file7)
 700	read(7,*,end=701)i,j
	  rij1=r0(i,1)-r0(j,1)
	  rij2=r0(i,2)-r0(j,2)
	  rij3=r0(i,3)-r0(j,3)
	  rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	rb0(i,j)=rij
	goto 700
 701	close(7)



!C tradueix la matriu de topologia
	do i=1,natred-1
	ii=iid(i)
	do j=i+1,natred
	jj=iid(j)
	rb(i,j)=rb0(ii,jj)
	drb(i,j)=sigma*rb(i,j)
	nstep(i,j)=0
	istruct(i,j)=0
	enddo
	enddo

	natom=natred
	natom1=natred1

	open(unit=19,file=file19)
        do i=1,natom
	write(19,1000)c1,i,atom(i),res(i),ind2(i),(r(i,j),j=1,3)
	enddo
	close(19)


!C potencials de Go entre els CA
	do i=1,natom1-1
	ii=ind2(i)
	do j=i+1,natom1
	jj=ind2(j)
	if(idock(ii).lt.2.and.idock(jj).lt.2)then
!C	if(atom(i).eq.'CA'.and.atom(j).eq.'CA')then
        rij1=r(i,1)-r(j,1)
        rij2=r(i,2)-r(j,2)
        rij3=r(i,3)-r(j,3)
        rij=sqrt(rij1**2+rij2**2+rij3**2)
	if(rij.lt.rcutgo.and.rb(i,j).lt.1.d-10)then
	rb(i,j)=rij
	drb(i,j)=rij*sigmago
	endif
!C	endif
	endif
	enddo
	enddo

!C potencials de Go entre els CA
	do i=natom1,natom-1
	ii=ind2(i)
	do j=i+1,natom
	jj=ind2(j)
	if(idock(ii).lt.2.and.idock(jj).lt.2)then
!C	if(atom(i).eq.'CA'.and.atom(j).eq.'CA')then
        rij1=r(i,1)-r(j,1)
        rij2=r(i,2)-r(j,2)
        rij3=r(i,3)-r(j,3)
        rij=sqrt(rij1**2+rij2**2+rij3**2)
	if(rij.lt.rcutgo.and.rb(i,j).lt.1.d-10)then
	rb(i,j)=rij
	drb(i,j)=rij*sigmago
	endif
!C	endif
	endif
	enddo
	enddo


!C opcio que fixa la posicio relativa de les proteines 
	if(irigd.eq.1)then
	do i=1,natom1
	ii=ind2(i)
	do j=natom1+1,natom
	jj=ind2(j)
        rij1=r(i,1)-r(j,1)
        rij2=r(i,2)-r(j,2)
        rij3=r(i,3)-r(j,3)
        rij=sqrt(rij1**2+rij2**2+rij3**2)
	if(idock(ii).lt.2.and.idock(jj).lt.2)then
	if(atom(i).eq.'CA'.and.atom(j).eq.'CA')then
	rb(i,j)=rij
	drb(i,j)=dr
	if((rij-dr).lt.4.)drb(i,j)=rij-4.
	endif
	endif
	enddo
	enddo
	endif

!Cccccccccccccccccccccccccccccccccccccccccccccccc
!C assigna un tipus a cada atom
!Cccccccccccccccccccccccccccccccccccccccccccccccc

	open(unit=7,file=file16)

	do i=1,natom
	if(atom(i).eq.'OXT')then
	atp(i)='o-'
	atp(i-1)='o-'
	atp(i-2)='cr2'
	else
 600	read(7,*,end=501)c3,c2,c4
	  if(res(i).eq.c3.and.atom(i).eq.c2)then
	  atp(i)=c4
!C		write(6,*)atp(i),atom(i),res(i),ind2(i)
	  rewind(7)
	  else
	  goto 600
	  endif
	endif
	enddo
	close(7)

!C nitrogen terminal
	if(atom(1).eq.'N')atp(1)='nk'
	if(atom(natom1+1).eq.'N')atp(natom1+1)='nk'


!Cccccccccccccccccccccccccccccccccccccccccccccccc
!C carrega els parametres de cada tipus d'atom
!Cccccccccccccccccccccccccccccccccccccccccccccccc

	open(unit=7,file=file17)
	open(unit=8,file=file18)
	do i=1,natom
	xlamb(i)=3.5
 500	read(7,*,end=501)c1,xq,xfree,xvol,xevdw,xrvdw,xrhc,xmassa
	read(8,*)c1,xq,xfrpy,xvol,xevdw,xrvdw,xrhc,xmassa
	if(atp(i).eq.c1)then
	qq(i)=xq
	gfree(i)=xfree
	gfrpy(i)=xfrpy
	vol(i)=xvol
	evdw(i)=xevdw
	rvdw(i)=xrvdw
	rhc(i)=xrhc
	xm(i)=xmassa
	rewind(7)
	rewind(8)
	else
	goto 500
	endif
	enddo
	close(7)
	close(8)

	goto 502

 501	write(6,*)'Atom sense tipus definit',i,atom(i),res(i),ind2(i)
	stop

 502	continue

	write(6,*)nres1,nres
	write(6,*)


!C bloc que reconeix els ponts disulfur 
	do i=1,natom
	if(atom(i).eq.'SG')then
	do j=i+1,natom
	if(atom(j).eq.'SG')then
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(rij.lt.2.5)then
	rb(i,j)=rij
	drb(i,j)=sigma*rij
	rij1=r(i-1,1)-r(j,1)
	rij2=r(i-1,2)-r(j,2)
	rij3=r(i-1,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	rb(i-1,j)=rij
	drb(i-1,j)=sigma*rij
	rij1=r(i,1)-r(j-1,1)
	rij2=r(i,2)-r(j-1,2)
	rij3=r(i,3)-r(j-1,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	rb(i,j-1)=rij
	drb(i,j-1)=sigma*rij
	endif
	endif
	enddo
	endif
	enddo


!C assigna radis de hardcore en funcio de l'element
	xmassa=0.
	do i=1,natom
	xmassa=xmassa+xm(i)
	igo(i)=0
	xlamb(i)=3.5
	v(i,1)=0.
	v(i,2)=0.
	v(i,3)=0.
	enddo


	c1='ATOM'




	temps=0.

	dhf2=dhf*dhf
	dpsint2=dpsint*dpsint
	dpsext2=dpsext*dpsext

	open(unit=20,file=file20)

!C	read(20,*)
	do i=1,natom
	read(20,2000)c1,j,c2,c3,c4,k,x,y,z
	r(i,1)=x
	r(i,2)=y
	r(i,3)=z
	enddo
!C	read(20,*)
	close(20)

	do i=1,nres1
	do j=nres1+1,nres
	nct(i,j)=0
	enddo
	enddo

	do i=1,natom1
	ii=ind2(i)
	do j=natom1+1,natom
	jj=ind2(j)
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(rij.lt.rcutoff)nct(ii,jj)=1
	enddo
	enddo

	ncont=0
	do i=1,nres1
	do j=nres1+1,nres
	if(ncnat(i,j).eq.1.and.nct(i,j).eq.1)then
	write(6,*)i,j
	ncont=ncont+1
	endif
	enddo
	enddo

	write(6,*)'Contactes complex:',ncontnat,ncont

	stop
	end

