!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C PROGRAMA DMDHB
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C DMD PER UNA PROTEINA. 
!C
!C FA SERVIR PERFILS DE POTENCIAL AMB NSTEP ESGLAONS
!C CONSIDERA ELS HIDROGENS AMINO
!C CREA I DESTRUEIX PONTS D'HIDROGEN
!C
!C DISTANCIES I VELOCITATS EN ANGSTROMS
!C TEMPERATURA EN K
!C ENERGIES EN kcal/mol
!C LLEGEIX LA TOPOLOGIA DE L'ARXIU file7, QUE HA DE SER UNA TAULA DE DUES COLUMNES
!C AMB TOTES LES PARELLES D'ATOMS UNITS PER BONDS O PSEUDOBONDS
!C CONSIDERA INTERACCIONS HIDROFOBIQUES A TRAVES DELS 
!C POTENCIALS DE SOLVATACIO DE LAZARIDIS-KARPLUS
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	program tdmdmstep
	implicit real*8(a-h)
	implicit real*8(o-z)
	integer*2 inter,istruct,ireg,icov
	parameter(natmax=5000,facte=4186.0d0)
	parameter(npairmax=40000)
	parameter(natp=20)
	character*5 atom,res
	common/xoc/r(natmax,3),v(natmax,3),xm(natmax),rbox,ierr
	common/pous/rstep(natmax,natmax,2),estep(natmax,natmax,2)
	common/intr/nstep(natmax,natmax),istruct(natmax,natmax),inter(natmax,natmax)
	common/cov/icov(natmax,natmax),rbound(npairmax),ibound(npairmax,2),rhc(natmax),sigma
	common/pdb/atom(natmax),res(natmax),ind2(natmax),nat(natmax),imol(natmax)
	common/atpres/ ihb(natmax),ica(natmax),io(natmax),ih(natmax),ico(natmax),in(natmax),icb(natmax)
	common/shake/ishk(natmax),nshk(natmax,natmax)
	common/npt/ipot(natmax),npot(natmax,natmax)
	common/fisic/evdw(natmax),rvdw(natmax),qq(natmax),gfree(natmax),vol(natmax)
	common/param/fvdw,fsolv,eps,xlamb
	common/parmsolv/rsolv,asolv,bsolv,dwat,icont(natmax),fcont(natmax)
	common/random_c/pseudo

	dimension qa(natmax,natp),gfreea(natmax,natp),va(natmax,natp)
	dimension evdwa(natmax,natp),rvdwa(natmax,natp),rhca(natmax,natp),xma(natmax,natp)
	dimension rant(natmax,3),ipart(natmax),tpart(natmax),ibeta(natmax),ind1(natmax)
	dimension rcm(3),vcm(3),vd(3),vm1(3),vm2(3)
	dimension inb1(natmax),inb2(natmax),nblist1(natmax,natmax),nblist2(natmax,natmax)
	dimension ireg(natmax,natmax),timp(natmax,natmax)
	dimension v1(20),v2(20)
	character*100 file7,file8,file9,file10,file11,file12,file13,file15,file16,file17,file19,file20,file21
	character*5 c1,c2,c3,c4,c5,cad(natmax),atp(natmax,natp)
	namelist /input/ tsnap,tene,tterm,tact,sigma,temp,kkk,file7,file9,file10,file12,file15,file16,file17,file19,file20,file21,nbloc,rcutgo,tmin,dcut,ego,isolv,fpot,ehb,ehbc,idab,igoab,irig,rshake,rpot,beta,ebond,dstep,dijmin,isec,tpush,iterm,factm,ihbr,fvdw,fsolv,eps,rbox,iwr,fterm,facthc,factr,icons,iprint,rsolv,asolv,bsolv,dwat
	namelist /distancies/rohmin,rohmax,rnomin,rnomax,rchmin,rchmax
 1000	FORMAT (A4,2X,I5,2X,A3,1X,A3,1X,A1,1X,I3,4X,F8.3,F8.3,F8.3)
 2000	FORMAT (A4,2X,I5,2X,A3,1X,A3,1X,A1,1X,I3,4X,F8.3,F8.3,F8.3,F8.3,1X,I2)
 333	FORMAT('MODEL',7X,I5)
 334	FORMAT('ENDMDL')
 1015	FORMAT(E12.6,20(1X,F7.2))


!C default values
	file9='nativain.pdb'
	file10='res'
	file11='distancia.dat'
	file12='energia.dat'
	file15='res'
	file19='output.pdb'
	file20='snapcg.pdb'
	file21='snapca.pdb'
	file7='topologia.dat'
	file16='atomtypes.dat'
	file17='potentials.dat'
	pi=atan(1.d0)*4.0d0
	a=1.d-10
	tmin=1.d-30
	dijmin=1.d-4
	kkk=2381
	temp=300.0d0
	nbloc=10000
	iwr=0
	iprint=1

	tact=2.d-14
	tene=1.d-12
	tsnap=1.d-11

	rcutgo=8.d0
	sigma=.05d0
	sigmago=0.1d0
	idab=0
	irig=0
	fterm=4.d0
	ebond=1000.d0
	dstep=1.d-4
	tpush=5.d-4
	isolv=1
	iterm=1
	factm=1.d0
	icm=0
	rbox=0.d0
	rshake=50.d0
	rpot=50.d0
	dcut=10.d0
	icons=1



!C Interaccions
	eps=16.5d0
	fvdw=8.d0
	fsolv=15.d0
	xlamb=3.5d0
	facthc=0.8d0
	factr=0.9d0


!C Ponts d'hidrogen
	isec=0
	ehbc=4.d0
	ehb=3.d0
	rohmin  = 1.75d0
	roha=2.15d0
	rohb=2.34d0
	rohmax  = 2.50d0
	rnomin  = 2.75d0
	rnoa=3.1d0
	rnob=3.2d0
	rnomax  = 3.50d0
	rchmin  = 2.90d0
	rcha=3.27d0
	rchb=3.4d0
	rchmax  = 3.75d0


	rsolv=3.5d0
	asolv=10.d0
	bsolv=0.5d0
	dwat=6.d0


	read(5,INPUT)
	read(5,DISTANCIES)

!C	write(6,INPUT)

	if(tene.gt.tsnap)tene=tsnap
	if(tact.gt.tene)tact=tene
	if(rbox.lt.1.d-10)then
	icm=1
	rbox=300.d0
	endif

	rbox2=0.5d0*rbox
	rshake2=rshake*rshake
	rpot2=rpot*rpot

	call random_seed()


!C llegeix les coordenades (fitxer pdb)
	open(unit=9,file=file9)
	kk=0
	n=1
	im=0
 50	read(9,1000,end=51)c1,j,c2,c3,c4,k,x,y,z
	atom(n)=c2
	res(n)=c3
	cad(n)=c4
	ind1(n)=k
!C	if(k.lt.ind1(n-1))then
!C	kk=kk+ind1(n-1)
!C	im=im+1
!C	endif
	if(k.lt.ind1(n-1))kk=kk+ind1(n-1)
	if(cad(n).ne.cad(n-1))im=im+1
	ind2(n)=k+kk
	imol(n)=im
	k1=ind2(n)
	r(n,1)=x
	r(n,2)=y
	r(n,3)=z
	if(atom(n).eq.'N')in(k1)=n
	if(atom(n).eq.'H')ih(k1)=n
	if(atom(n).eq.'CA')ica(k1)=n
	if(atom(n).eq.'C')ico(k1)=n
	if(atom(n).eq.'O')io(k1)=n
	n=n+1
	goto 50
 51	natom=n-1
	close(9)
	nres=ind2(natom)

	if(natom.gt.natmax)then
	write(6,*)'El numero de particules supera el limit maxim de ',natmax
	stop
	endif

!Cccccccccccccccccccccccccccccccccccccccccccccccc
!C assigna un tipus a cada atom
!Cccccccccccccccccccccccccccccccccccccccccccccccc

	open(unit=7,file=file16)
	do i=1,natom
	nat(i)=0
	if(atom(i).eq.'N'.or.atom(i).eq.'H'.or.atom(i).eq.'C'.or.atom(i).eq.'O'.or.atom(i).eq.'OXT')then
	nat(i)=1
	else
 70	read(7,*,end=71)c2,c3,c4
	if(atom(i).eq.c2.and.c3.eq.res(i))then
	  nat(i)=nat(i)+1
	  j=nat(i)
	  atp(i,j)=c4
	endif
	goto 70
 71	rewind(7)
	endif
	enddo
	close(7)

	do i=1,natom
	if(atom(i).eq.'N')atp(i,1)='nh'
	if(atom(i).eq.'H')atp(i,1)='h'
	if(atom(i).eq.'C')atp(i,1)='co'
	if(atom(i).eq.'O'.or.atom(i).eq.'OXT')atp(i,1)='oc'
	enddo

!Cccccccccccccccccccccccccccccccccccccccccccccccc
!C carrega els parametres de cada tipus d'atom
!Cccccccccccccccccccccccccccccccccccccccccccccccc

	open(unit=7,file=file17)
	do i=1,natom
	do j=1,nat(i)
 500	read(7,*)c1,xq,xfree,xvol,xevdw,xrvdw,xrhc,xmassa
	if(atp(i,j).eq.c1)then
	qa(i,j)=xq
	gfreea(i,j)=xfree
	va(i,j)=xvol
	evdwa(i,j)=xevdw
	rvdwa(i,j)=xrvdw
	rhca(i,j)=0.8d0*xrvdw
	xma(i,j)=xmassa
	rewind(7)
	else
	goto 500
	endif
	enddo
	enddo
	close(7)

!Ccccccccccccccccccccccccccccccccccccccc
!C propietats de les boles
!Ccccccccccccccccccccccccccccccccccccccc
	xmassa=0.d0
	do i=1,natom
	xm(i)=0.d0
	qq(i)=0.d0
	vol(i)=0.d0
	gfree(i)=0.d0
	evdw(i)=0.d0
	sumrhc=0.d0
	sumrvdw=0.d0
	  if(nat(i).eq.1)then
	  xm(i)=xma(i,1)
	  qq(i)=qa(i,1)
	  vol(i)=va(i,1)
	  gfree(i)=gfreea(i,1)
	  evdw(i)=evdwa(i,1)
	  rvdw(i)=rvdwa(i,1)
	  rhc(i)=0.8d0*rvdw(i)
	  else
	    do j=1,nat(i)
	    xm(i)=xm(i)+xma(i,j)
	    qq(i)=qq(i)+qa(i,j)
	    vol(i)=vol(i)+va(i,j)
	    gfree(i)=gfree(i)+gfreea(i,j)
	    evdw(i)=evdw(i)+evdwa(i,j)
	    sumrhc=sumrhc+rhca(i,j)**3.d0
	    sumrvdw=sumrvdw+rvdwa(i,j)**3.d0
	    enddo
!C	    xm(i)=0.012
	  rvdw(i)=factr*sumrvdw**0.3333d0
	  rhc(i)=0.8d0*rvdw(i)
	  endif
	xmassa=xmassa+xm(i)
	v(i,1)=0.d0
	v(i,2)=0.d0
	v(i,3)=0.d0
	enddo

!Ccccccccccccccccccccccccccccccccccccccc
!C llegeix la matriu de topologia
!Ccccccccccccccccccccccccccccccccccccccc
	do i=1,natom-1
	do j=i+1,natom
	icov(i,j)=0
!C	rb(i,j)=0.
	nstep(i,j)=0
	inter(i,j)=0
	istruct(i,j)=0
	enddo
	enddo

	k=0
	open(unit=7,file=file7)
 700	read(7,*,end=701)i,j,rij
	icov(i,j)=1
!C	rb(i,j)=rij
	k=k+1
	ibound(k,1)=i
	ibound(k,2)=j
	rbound(k)=rij
	goto 700
 701	close(7)

	npair=k

	!write(6,*)npair
	!read(5,*)debug

	open(unit=8,file='dmd.out')
	write(8,INPUT)


	c1='ATOM'

!C reconeix estructura secundaria i estableix ponts d'hidrogen

	nhb=0

!C	do i=1,nres
!C	ii=io(i)
!C	natp(ii)=9
!C	enddo

	do i=1,natom
	ihb(i)=0
	enddo


!C	write(6,*)'ponts hidrogen O --> H'
	do i=1,nres-4
	ii=io(i)
	do j=i+4,nres
	if(res(ica(j)).ne.'PRO')then
	jj=ih(j)
	if(ihb(ii).eq.0.and.ihb(jj).eq.0)then
	n1=io(i)
	n2=ih(j)
	rij1=dbox(n2,n1,1)
	rij2=dbox(n2,n1,2)
	rij3=dbox(n2,n1,3)
	  roh=(rij1*rij1+rij2*rij2+rij3*rij3)**0.5d0
	if(roh.lt.rohmax.and.roh.gt.rohmin)then
	n1=io(i)
	n2=in(j)
	rij1=dbox(n2,n1,1)
	rij2=dbox(n2,n1,2)
	rij3=dbox(n2,n1,3)
	  rno=(rij1*rij1+rij2*rij2+rij3*rij3)**0.5d0
	if(rno.lt.rnomax.and.rno.gt.rnomin)then
	n1=ico(i)
	n2=ih(j)
	rij1=dbox(n2,n1,1)
	rij2=dbox(n2,n1,2)
	rij3=dbox(n2,n1,3)
	  rch=(rij1*rij1+rij2*rij2+rij3*rij3)**0.5d0
	if(rch.lt.rchmax.and.rch.gt.rchmin)then
	nhb=nhb+1
	n1=io(i)
	n2=ih(j)
	ihb(n1)=1
	ihb(n2)=1
	write(6,*)'HBOND ',atom(ii),res(ii),ind2(ii),atom(jj),res(jj),ind2(jj)
	write(8,*)'HBOND ',atom(ii),res(ii),ind2(ii),atom(jj),res(jj),ind2(jj)
	endif
	endif
	endif
	endif
	endif
	enddo
	enddo

!C	write(6,*)'ponts hidrogen H --> O'
	do i=1,nres-4
	if(res(ica(i)).ne.'PRO')then
	ii=ih(i)
	do j=i+4,nres
	jj=io(j)
	if(ihb(ii).eq.0.and.ihb(jj).eq.0)then
	n1=ih(i)
	n2=io(j)
	rij1=dbox(n2,n1,1)
	rij2=dbox(n2,n1,2)
	rij3=dbox(n2,n1,3)
	  roh=(rij1*rij1+rij2*rij2+rij3*rij3)**0.5d0
	if(roh.lt.rohmax.and.roh.gt.rohmin)then
	n1=in(i)
	n2=io(j)
	rij1=dbox(n2,n1,1)
	rij2=dbox(n2,n1,2)
	rij3=dbox(n2,n1,3)
	  rno=(rij1*rij1+rij2*rij2+rij3*rij3)**0.5d0
	if(rno.lt.rnomax.and.rno.gt.rnomin)then
	n1=ih(i)
	n2=ico(j)
	rij1=dbox(n2,n1,1)
	rij2=dbox(n2,n1,2)
	rij3=dbox(n2,n1,3)
	  rch=(rij1*rij1+rij2*rij2+rij3*rij3)**0.5d0
	if(rch.lt.rchmax.and.rch.gt.rchmin)then
	nhb=nhb+1
	n1=ih(i)
	n2=io(j)
	ihb(n1)=1
	ihb(n2)=1
	write(6,*)'HBOND ',atom(ii),res(ii),ind2(ii),atom(jj),res(jj),ind2(jj)
	write(8,*)'HBOND ',atom(ii),res(ii),ind2(ii),atom(jj),res(jj),ind2(jj)
	endif
	endif
	endif
	endif
	enddo
	endif
	enddo
	
	!write(6,*)nhb
	!read(5,*)debug

!C	do i=1,natom-1
!C	do j=i+1,natom
!C	if(inter(i,j).eq.1)write(6,*)i,j,inter(i,j)
!C	enddo
!C	enddo
!C	stop


	call potencial(natom)

!C llista de solapaments plausibles
	do i=1,natom-1
	ishk(i)=0
	ipot(i)=0
	do j=i+1,natom
	if(icov(i,j).eq.0)then
	rij1=dbox(i,j,1)
	rij2=dbox(i,j,2)
	rij3=dbox(i,j,3)
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	  if(rmod2.lt.rshake2)then
	  ishk(i)=ishk(i)+1
	  k=ishk(i)
	  nshk(i,k)=j
	  endif
!C
	  if(istruct(i,j).ne.1)then
	  if(rmod2.lt.rpot2)then
!write(6,fmt="(a)",advance="no")'x'
	  ipot(i)=ipot(i)+1
	  k=ipot(i)
	  npot(i,k)=j
	  endif
	  endif
	endif
	enddo
!write(6,*)''
	enddo
!read(5,*)debug


	if(isolv.ne.0)call enchufa(natom,dcut)

!C assigna la regio on es troba la interaccio entre dues particules
	do i=1,natom-1
!write(6,fmt="(i4,a)",advance="no")i,": "
	do j=i+1,natom
	if(icov(i,j).eq.0)then
	rij1=dbox(j,i,1)
	rij2=dbox(j,i,2)
	rij3=dbox(j,i,3)
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	rij=(rmod2**0.5d0)
	k=1
	  do while (rij.gt.rstep(i,j,k).and.k.le.nstep(i,j))
	  k=k+1
	  enddo
	ireg(i,j)=k
!write(6,fmt="(i4,1x)",advance="no")k
!C	if(istruct(i,j).eq.1)write(6,*)i,j,(estep(i,j,k),k=1,nstep(i,j))
	endif
	enddo
!write(6,*)''
	enddo
!read(5,*)debug



!C suma l'energia potencial de la conformacio inicial
	epot0=0.d0
	epotmol0=0.d0
	epothb0=0.d0
	epothbmol0=0.d0
	do i=1,natom-1
	do j=i+1,natom
	rij1=dbox(i,j,1)
	rij2=dbox(i,j,2)
	rij3=dbox(i,j,3)
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	dist=(rmod2**0.5d0)
	if(inter(i,j).eq.1)then
	  k=nstep(i,j)
	  do while(dist.lt.rstep(i,j,k).and.k.gt.0)
	  epot0=epot0-estep(i,j,k)
	  if(imol(i).ne.imol(j))epotmol0=epotmol0-estep(i,j,k)
	  if(istruct(i,j).eq.1)then
	  epothb0=epothb0-estep(i,j,k)
	  if(imol(i).ne.imol(j))epothbmol0=epothbmol0-estep(i,j,k)
	  endif
	  k=k-1
	  enddo
	endif
	enddo
	enddo


!C assigna velocitats aleatories
	do j=1,3
	vcm(j)=0.d0
	do i=1,natom
	call random_number_new(fi)
	v(i,j)=fi
	vcm(j)=vcm(j)+xm(i)*v(i,j)
!write(6,*)v(i,j),vcm(j)
	enddo
	vcm(j)=vcm(j)/xmassa
	enddo
!read(5,*)debug
!write(6,fmt="(f20.16,1x,f20.16,1x,f20.16)")vcm(1),vcm(2),vcm(3)
!call imprimir_v(natom,1)
!C ajusta l'energia cinetica a la temperatura requerida
	ekin=0.0d0
	do j=1,3
	do i=1,natom
	v(i,j)=v(i,j)-vcm(j)
	ekin=ekin+0.5d0*xm(i)*(v(i,j)*a)**2.d0
!write(6,fmt="(f20.16,1x,f20.16,1x,f20.16,1x,f20.16)")(0.5d0*xm(i)*(v(i,j)*a)**2.d0)*1.d24,((v(i,j)*a))*1.d11,((v(i,j)*a)**2)*1.d21,v(i,j)
!write(6,*)v(i,j),ekin
	enddo
	enddo

!write(6,fmt="(f20.16)")ekin*1.d20
	sto=1.5d0*8.314d0*natom*temp/ekin
	ekin0=0.0d0
!write(6,fmt="(f20.16,1x,f20.16)")sto*1.d-26,sqrt(sto)*1.d-13
!write(6,fmt="(f40.35,1x,f40.35,1x,f40.35,1x,f40.35)")1.e-11,1.d-11,dble(1.e-11),real(1.d-11)
!read(5,*)debug
	do j=1,3
	do i=1,natom
	v(i,j)=v(i,j)*(sto**0.5d0)
	ekin0=ekin0+0.5d0*xm(i)*(v(i,j)*a)**2.d0
!write(6,*)v(i,j),ekin0
	enddo
	enddo
	ekin0=ekin0/facte

	etot0=epot0+ekin0
!write(6,*)ekin0,etot0
!read(5,*)debug
!call imprimir_v(natom,1)
!C ara busca el CM
	do j=1,3
	rcm(j)=0.0d0
	do i=1,natom
	rcm(j)=rcm(j)+xm(i)*r(i,j)
	enddo
        rcm(j)=rcm(j)/xmassa
	enddo

	ibloc=0
!call imprimir_r(natom,1)
!do i=1,3
!	write(6,fmt="(f50.45,1x,f50.45)")rcm(i),rbox2
!enddo
!read(5,*)debug
!C escriu les coordenades en el SRCM
	open(unit=12,file='input.pdb')
	open(unit=20,file=file20)
	open(unit=21,file=file21)
	write(20,333)ibloc
	write(21,333)ibloc
	do i=1,natom
	do j=1,3
	r(i,j)=r(i,j)-rcm(j)+rbox2
	if(r(i,j).gt.rbox)r(i,j)=r(i,j)-rbox
	if(r(i,j).lt.0.d0)r(i,j)=r(i,j)+rbox
	enddo
	c2=atom(i)
	write(12,2000)c1,i,c2,res(i),cad(i),ind1(i),(r(i,j),j=1,3),fcont(i),icont(i)
	write(20,1000)c1,i,c2,res(i),cad(i),ind1(i),(r(i,j),j=1,3)
	if(atom(i).eq.'CA')write(21,1000)c1,i,c2,res(i),cad(i),ind1(i),(r(i,j),j=1,3)
	enddo
!call imprimir_r(natom,1)
	close(12)
	write(20,334)
	write(21,334)

	  do i=1,natom
	  do j=1,3
	  rant(i,j)=r(i,j)
	  enddo
	  enddo



	do i=1,natom-1
	do j=i+1,natom
	timp(i,j)=1.d0
	enddo
	enddo



	write(6,*)'natom,nres,tsnap,nbloc',natom,nres,tsnap,nbloc
	write(6,*)'epot,ekin,nhb',epot0,ekin0,etot0,nhb,epothb0

	open(unit=8,file='dmd.out')
	write(8,INPUT)
	write(8,DISTANCIES)
	write(8,*)'natom,nres,tsnap,nbloc',natom,nres,tsnap,nbloc
	write(8,*)'epot,ekin,nhb',epot0,ekin0,etot0,nhb


!C comença a iterar per buscar l'event mes proper
!C la variable ixoc indica quin tipus d'event ocorre:
!C 0 -> enllaç
!C 1 -> xoc entre atoms no enllaçats

!C	open(unit=10,file='events.dat')
	open(unit=11,file=file11)
	open(unit=12,file=file12)

	temps=0.d0
	temps0=0.d0
	tevent=0.d0

	write(12,*)'# Energia inicial ',epot0,epotmol0,epothb0,epothbmol0,ekin0,etot0,nhb


	icrash=0

	do 200 ibloc=1,nbloc

	tacum=0.d0

	icg=0

icounter_1 = 0
icounter_2 = 0

	do 40 while(tacum.lt.tsnap)
!write(6,*)'counter_1: ',icounter_1
	do i=1,natom-1
	do j=i+1,natom
	if(istruct(i,j).eq.0)inter(i,j)=0
	enddo
	enddo

	if(isec.eq.0)then


	do i=1,natom-1
	do j=i+1,natom
	inter(i,j)=0
	enddo
	enddo

	call potencial(natom)

	nhb=0


	do i=1,natom
	ihb(i)=0
	enddo

!C	write(6,*)'ponts hidrogen O --> H'
	do i=1,nres-4
	ii=io(i)
	do j=i+4,nres
	if(res(ica(j)).ne.'PRO')then
	jj=ih(j)
	if(ihb(ii)+ihb(jj).eq.0)then
	ic=0
	n1=io(i)
	n2=ih(j)
	rij1=dbox(n2,n1,1)
	rij2=dbox(n2,n1,2)
	rij3=dbox(n2,n1,3)
	  roh=(rij1*rij1+rij2*rij2+rij3*rij3)**0.5d0
	if(roh.lt.rohmax.and.roh.gt.rohmin)then
	  ic=ic+1
	else
	  if(istruct(n1,n2).eq.1)then
	    if(ireg(n1,n2).le.nstep(n1,n2).and.ireg(n1,n2).gt.1) then
		ic=ic+1
		endif
	  endif
	endif
	n1=io(i)
	n2=in(j)
	rij1=dbox(n2,n1,1)
	rij2=dbox(n2,n1,2)
	rij3=dbox(n2,n1,3)
	  rno=(rij1*rij1+rij2*rij2+rij3*rij3)**0.5d0
	if(rno.lt.rnomax.and.rno.gt.rnomin)then
	  ic=ic+1
	else
	  if(istruct(n1,n2).eq.1)then
	    if(ireg(n1,n2).le.nstep(n1,n2).and.ireg(n1,n2).gt.1) then
		ic=ic+1
		endif
	  endif
	endif
	n1=ico(i)
	n2=ih(j)
	rij1=dbox(n2,n1,1)
	rij2=dbox(n2,n1,2)
	rij3=dbox(n2,n1,3)
	  rch=(rij1*rij1+rij2*rij2+rij3*rij3)**0.5d0
	if(rch.lt.rchmax.and.rch.gt.rchmin)then
	  ic=ic+1
	else
	  if(istruct(n1,n2).eq.1)then
	    if(ireg(n1,n2).le.nstep(n1,n2).and.ireg(n1,n2).gt.1) then
		ic=ic+1
		endif
	  endif
	endif
	if(ic.eq.3)then
!write(6,*)'ic==3 in ',i-1,j-1
	nhb=nhb+1
	n1=io(i)
	n2=ih(j)
	ihb(n1)=1
	ihb(n2)=1
	sto=fcont(n1)*fcont(n2)
	ene=ehb*sto+ehbc*(1.-sto)
	call creapouhb(n1,n2,rohmin,roha,rohb,rohmax,ene)
	ireg(n1,n2)=2
	n1=io(i)
	n2=in(j)
	sto=fcont(n1)*fcont(n2)
	ene=ehb*sto+ehbc*(1.-sto)
	call creapouhb(n1,n2,rnomin,rnoa,rnob,rnomax,ene)
	ireg(n1,n2)=2
	n1=ico(i)
	n2=ih(j)
	sto=fcont(n1)*fcont(n2)
	ene=ehb*sto+ehbc*(1.-sto)
	call creapouhb(n1,n2,rchmin,rcha,rchb,rchmax,ene)
	ireg(n1,n2)=2
	else
	n1=io(i)
	n2=ih(j)
!write(6,*)'a: io[',i-1,']=',n1,' ih[',j-1,']=',n2
	istruct(n1,n2)=0
	ireg(n1,n2)=0
	n1=io(i)
	n2=in(j)
!write(6,*)'b: io[',i-1,']=',n1,' in[',j-1,']=',n2
	istruct(n1,n2)=0
	ireg(n1,n2)=0
	n1=ico(i)
	n2=ih(j)
!write(6,*)'c: ico[',i-1,']=',n1,' ih[',j-1,']=',n2
	istruct(n1,n2)=0
	ireg(n1,n2)=0
	endif
	endif
	endif
	enddo
!if(i.eq.20)read(5,*)debug
	enddo
!read(5,*)debug

!C	write(6,*)'ponts hidrogen H --> O'
	do i=1,nres-4
	if(res(ica(i)).ne.'PRO')then
	ii=ih(i)
	do j=i+4,nres
	jj=io(j)
	if(ihb(ii)+ihb(jj).eq.0)then
	ic=0
	n1=ih(i)
	n2=io(j)
	rij1=dbox(n2,n1,1)
	rij2=dbox(n2,n1,2)
	rij3=dbox(n2,n1,3)
	  roh=(rij1*rij1+rij2*rij2+rij3*rij3)**0.5d0
	if(roh.lt.rohmax.and.roh.gt.rohmin)then
	  ic=ic+1
	else
	  if(istruct(n1,n2).eq.1)then
	    if(ireg(n1,n2).le.nstep(n1,n2).and.ireg(n1,n2).gt.1)ic=ic+1
	  endif
	endif
	n1=in(i)
	n2=io(j)
	rij1=dbox(n2,n1,1)
	rij2=dbox(n2,n1,2)
	rij3=dbox(n2,n1,3)
	  rno=(rij1*rij1+rij2*rij2+rij3*rij3)**0.5d0
	if(rno.lt.rnomax.and.rno.gt.rnomin)then
	  ic=ic+1
	else
	  if(istruct(n1,n2).eq.1)then
	    if(ireg(n1,n2).le.nstep(n1,n2).and.ireg(n1,n2).gt.1)ic=ic+1
	  endif
	endif
	n1=ih(i)
	n2=ico(j)
	rij1=dbox(n2,n1,1)
	rij2=dbox(n2,n1,2)
	rij3=dbox(n2,n1,3)
	  rch=(rij1*rij1+rij2*rij2+rij3*rij3)**0.5d0
	if(rch.lt.rchmax.and.rch.gt.rchmin)then
	  ic=ic+1
	else
	  if(istruct(n1,n2).eq.1)then
	    if(ireg(n1,n2).le.nstep(n1,n2).and.ireg(n1,n2).gt.1)ic=ic+1
	  endif
	endif
	if(ic.eq.3)then
!write(6,*)'ic==3 in ',i-1,j-1
	nhb=nhb+1
	n1=ih(i)
	n2=io(j)
	ihb(n1)=1
	ihb(n2)=1
	sto=fcont(n1)*fcont(n2)
	ene=ehb*sto+ehbc*(1.-sto)
	call creapouhb(n1,n2,rohmin,roha,rohb,rohmax,ene)
	ireg(n1,n2)=2
	n1=in(i)
	n2=io(j)
	sto=fcont(n1)*fcont(n2)
	ene=ehb*sto+ehbc*(1.d0-sto)
	call creapouhb(n1,n2,rnomin,rnoa,rnob,rnomax,ene)
	ireg(n1,n2)=2
	n1=ih(i)
	n2=ico(j)
	sto=fcont(n1)*fcont(n2)
	ene=ehb*sto+ehbc*(1.d0-sto)
	call creapouhb(n1,n2,rchmin,rcha,rchb,rchmax,ene)
	ireg(n1,n2)=2
	else
	n1=ih(i)
	n2=io(j)
!write(6,*)'a: ih[',i-1,']=',n1,' io[',j-1,']=',n2
	istruct(n1,n2)=0
	ireg(n1,n2)=0
	n1=in(i)
	n2=io(j)
!write(6,*)'b: in[',i-1,']=',n1,' io[',j-1,']=',n2
	istruct(n1,n2)=0
	ireg(n1,n2)=0
	n1=ih(i)
	n2=ico(j)
!write(6,*)'c: ih[',i-1,']=',n1,' ico[',j-1,']=',n2
	istruct(n1,n2)=0
	ireg(n1,n2)=0
	endif
	endif
	enddo
	endif
!if(i.eq.20)read(5,*)debug
	enddo
	endif
!read(5,*)debug
	call potencial(natom)

	do i=1,natom-1
!write(6,fmt="(i4,a)", advance="no")i-1,': '
	do j=i+1,natom
!write(6,fmt="(i4,1x)", advance="no")ireg(i,j)
	if(ireg(i,j).eq.0.and.icov(i,j).eq.0)then
	rmin2=rstep(i,j,1)*rstep(i,j,1)
	rmax2=rstep(i,j,2)*rstep(i,j,2)
	rij1=dbox(i,j,1)
	rij2=dbox(i,j,2)
	rij3=dbox(i,j,3)
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
!write(6,fmt="(f16.6,1x,f16.6,1x,f16.6,1x)", advance="no")rmod2,rmin2,rmax2
	if(rmod2.lt.rmin2)then
	ireg(i,j)=1
	elseif(rmod2.gt.rmax2)then
	ireg(i,j)=3
	else
	ireg(i,j)=2
	endif
	endif
!write(6,fmt="(i4,1x)", advance="no")ireg(i,j)
	enddo
!write(6,*)''
	enddo
!read(5,*)debug
	
	
!C llista de solapaments plausibles
	do i=1,natom-1
	ishk(i)=0
	ipot(i)=0
	do j=i+1,natom
	if(icov(i,j).eq.0)then

	rij1=dbox(i,j,1)
	rij2=dbox(i,j,2)
	rij3=dbox(i,j,3)
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	  if(rmod2.lt.rshake2)then
!write(6,fmt="(a)",advance="no")'x'
	  ishk(i)=ishk(i)+1
	  k=ishk(i)
	  nshk(i,k)=j
	  endif
!C
	  if(istruct(i,j).ne.1)then
	  if(rmod2.lt.rpot2)then
!write(6,fmt="(a)",advance="no")'y'
	  ipot(i)=ipot(i)+1
	  k=ipot(i)
	  npot(i,k)=j
	  endif
	  endif
	endif
	enddo
!write(6,*)''
	enddo
!read(5,*)

	mem1=0
	mem2=0


	iev=0
	ierr=0
	ierr2=0

 20	temps0=temps

	  do i=1,natom
	  do j=1,3
	  r(i,j)=rant(i,j)
	  enddo
	  enddo

	tacene=0.

	do 30 while(tacene.lt.tene)
!write(6,*)'counter_2: ',icounter_2
!if((ibloc.eq.2).and.(icounter_2.eq.37))call imprimir_v(natom,1)

!write(6,*)'counter_2: ',icounter_2
!if (icounter_2.eq.1)call imprimir_inter(natom)
	tacact=0.d0

	icint=0
!write(6,fmt="(a,f12.4)")'<-->',v(38,1)*1E-10
!write(6,fmt="(a,f12.4)")'<-->',v(569,1)*1E-10
	call dmdshake(natom,npair,mem1,mem2)
!write(6,fmt="(a,f12.4)")'>--<',v(38,1)*1E-10
!write(6,fmt="(a,f12.4)")'>--<',v(569,1)*1E-10
	do i=1,natom-1
	do j=i+1,natom
	if(istruct(i,j).eq.0)inter(i,j)=0
	enddo
	enddo
!call imprimir_v(natom,1)
	if(isolv.ne.0)call enchufa(natom,dcut)
!if (icounter_2.eq.1)call imprimir_inter(natom)
!if (icounter_2.eq.1)read(5,*)debug
!call imprimir_v(natom,1)
!C fa tornar a la seva regio els parells que han canviat de regio sense que el programa se n'adoni
	icont=0
	icont_2=0

	do i=1,natom-1
!if(icounter_2.eq.14)write(6,fmt="(i4,1x)", advance="no")i-1
!if(i.eq.38)write(6,fmt="(i4,a,f12.4)")i-1,'->',v(i,1)*1E-10
	do j=i+1,natom
!if(j.eq.38)write(6,fmt="(i4,a,f12.4)")j-1,'|>',v(j,1)*1E-10
	if(inter(i,j).eq.1)then
	rij1=dbox(j,i,1)
	rij2=dbox(j,i,2)
	rij3=dbox(j,i,3)
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	vij1=v(i,1)-v(j,1)
	vij2=v(i,2)-v(j,2)
	vij3=v(i,3)-v(j,3)
	vmod2=vij1*vij1+vij2*vij2+vij3*vij3
	prod=rij1*vij1+rij2*vij2+rij3*vij3
	rij=(rmod2**0.5d0)
!if((ibloc.eq.2).and.(icounter_2.eq.37).and.(i.eq.1052))write(6,fmt="(i4,1x,i4,1x,f40.35,1x,f40.35,1x,f40.35,1x,f40.35,1x,f40.35,1x,f40.35)")i,j,v(i,1)*1.0d-11,v(i,2)*1.0d-11,v(i,3)*1.0d-11,v(j,1)*1.0d-11,v(j,2)*1.0d-11,v(j,3)*1.0d-11
	k=1
	  do while (rij.gt.rstep(i,j,k).and.k.le.nstep(i,j))
	  k=k+1
	  enddo
!C	write(6,*)i,j,k,rstep(i,j,k-1),rij,rstep(i,j,k),timp(i,j)
	k0=ireg(i,j)
!write(6,fmt="(i4,1x)", advance="no")k0
	ich=0
!if((i.eq.12).and.(j.eq.15))write(6,fmt="(i4,1x,i4,a,i1,1x,i1,2x,f12.4)")i,j,"==>",k,k0,prod*1E-10
	if(k0.gt.k.and.prod.lt.0)then
!if (icounter_2.eq.14)write(6,fmt="(i4,a,f12.4,1x,f12.4,1x,f12.4)")j-1,'=',prod*1E-10,v(i,1)*1E-10,v(j,1)*1E-10
	icont=icont+1
	icont_2 = icont_2 + 1
	dpot=-estep(i,j,k0-1)*facte
!if((ibloc.eq.2).and.(icounter_2.eq.37).and.(i.eq.1052).and.(j.eq.1077))write(6,fmt="(i4,1x,i4,1x,f45.35)")k,k0,dpot
!if((ibloc.eq.2).and.(icounter_2.eq.37).and.(i.eq.1052).and.(j.eq.1077))write(6,fmt="(f40.35,1x,f40.35,1x,f40.35,1x,f40.35,1x,f40.35,1x,f40.35)")v(i,1)*1.0d-11,v(i,2)*1.0d-11,v(i,3)*1.0d-11,v(j,1)*1.0d-11,v(j,2)*1.0d-11,v(j,3)*1.0d-11
	call chgmomene(i,j,rij1,rij2,rij3,dpot,ich)
!if((ibloc.eq.2).and.(icounter_2.eq.37).and.(i.eq.1052).and.(j.eq.1077))write(6,fmt="(f40.35,1x,f40.35,1x,f40.35,1x,f40.35,1x,f40.35,1x,f40.35)")v(i,1)*1.0d-11,v(i,2)*1.0d-11,v(i,3)*1.0d-11,v(j,1)*1.0d-11,v(j,2)*1.0d-11,v(j,3)*1.0d-11
	if(ich.eq.1)ireg(i,j)=ireg(i,j)-1
	ierr2=ierr2+1
!if(j.eq.38)write(6,fmt="(i4,a,f12.4)")j-1,')>',v(j,1)*1E-10
	elseif(k0.lt.k.and.prod.gt.0)then
	icont=icont+1
	icont_2 = icont_2 + 1
	dpot=estep(i,j,k0)*facte
	call chgmomene(i,j,rij1,rij2,rij3,dpot,ich)
	if(ich.eq.1)ireg(i,j)=ireg(i,j)+1
	ierr2=ierr2+1
!if(j.eq.38)write(6,fmt="(i4,a,f12.4)")j-1,'(>',v(j,1)*1E-10
	endif
!if((ibloc.eq.2).and.(icounter_2.eq.37).and.(i.eq.1052))write(6,fmt="(i4,1x,i4,1x,f40.35,1x,f40.35,1x,f40.35,1x,f40.35,1x,f40.35,1x,f40.35)")i,j,v(i,1)*1.0d-11,v(i,2)*1.0d-11,v(i,3)*1.0d-11,v(j,1)*1.0d-11,v(j,2)*1.0d-11,v(j,3)*1.0d-11
	endif
!if(j.eq.38)write(6,fmt="(i4,a,f12.4)")j-1,'/>',v(j,1)*1E-10
	enddo
!if(i.eq.38)write(6,fmt="(i4,a,f12.4)")i-1,'=>',v(i,1)*1E-10
!if (icounter_2.eq.14)write(6,*)''
	enddo
!if((ibloc.eq.21).and.(icounter_2.eq.0))call imprimir_v(natom,1)
!read(5,*)debug

!write(6,*)icont(1),' ',ierr2

!if (icounter_2.eq.14) then
!do i=1,natom
!write(6,fmt="(f12.2,1x,f12.2,1x,f12.2)")v(i,1)*1E-10,v(i,2)*1E-10,v(i,3)*1E-10
!enddo
!read(5,*)debug
!endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C translacio i variacio dels temps
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!write(6,fmt="(f12.4,1x,f12.4,1x,f12.4)")dbox(12,471,1),dbox(12,471,2),dbox(12,471,3)
!write(6,fmt="(f12.4,1x,f12.4,1x,f12.4)")v(12,1)*1E-10,v(12,2)*1E-10,v(12,3)*1E-10
!write(6,fmt="(f12.4,1x,f12.4,1x,f12.4)")v(471,1)*1E-10,v(471,2)*1E-10,v(471,3)*1E-10
!write(6,fmt="(f12.4,1x,f12.4,1x,f12.4)")v(15,1)*1E-10,v(15,2)*1E-10,v(15,3)*1E-10
	do j=1,3
	do i=1,natom
	r(i,j)=r(i,j)+tact*v(i,j)
	enddo
	enddo
!write(6,fmt="(f12.4,1x,f12.4,1x,f12.4)")dbox(12,471,1),dbox(12,471,2),dbox(12,471,3)
!write(6,fmt="(f12.4,1x,f12.4,1x,f12.4)")v(12,1)*1E-10,v(12,2)*1E-10,v(12,3)*1E-10
!write(6,fmt="(f12.4,1x,f12.4,1x,f12.4)")v(471,1)*1E-10,v(471,2)*1E-10,v(471,3)*1E-10
!write(6,fmt="(f12.4,1x,f12.4,1x,f12.4)")v(15,1)*1E-10,v(15,2)*1E-10,v(15,3)*1E-10
	if(icm.eq.0)then
	do i=1,natom
	do j=1,3
	if(r(i,j).gt.rbox)r(i,j)=r(i,j)-rbox
	if(r(i,j).lt.0.d0)r(i,j)=r(i,j)+rbox
	enddo
	enddo
	endif

!write(6,fmt="(f12.4,1x,f12.4,1x,f12.4)")dbox(12,471,1),dbox(12,471,2),dbox(12,471,3)
!write(6,fmt="(f12.4,1x,f12.4,1x,f12.4)")v(12,1)*1E-10,v(12,2)*1E-10,v(12,3)*1E-10
!write(6,fmt="(f12.4,1x,f12.4,1x,f12.4)")v(471,1)*1E-10,v(471,2)*1E-10,v(471,3)*1E-10
!write(6,fmt="(f12.4,1x,f12.4,1x,f12.4)")v(15,1)*1E-10,v(15,2)*1E-10,v(15,3)*1E-10


 110	tacum=tacum+tact
	tacene=tacene+tact
	tacterm=tacterm+tact
	temps=temps+tact

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
iterm=0
	if(iterm.eq.1)then

!C termostat Andersen
!C selecciona una particula que termalitzar
	call random_number_new(fi)
	i=int(natom*fi)+1
!!write(6,*)'random: ',i,' ',fi
	do j=1,3
	call rnd_gauss(fi,xm(i),temp)
	v(i,j)=fi/a
	enddo

	endif
	
	icounter_2 = icounter_2 + 1
!if(icounter_2.eq.16)read(5,*)debug

 30	continue
 
!do i=1,natom
!write(6,*)v(i,1),v(i,2),v(i,3)
!enddo
!read(5,*)debug

!C	stop

	ekin=0.d0
	do j=1,3
	do i=1,natom
	ekin=ekin+0.5d0*xm(i)*(v(i,j)*a)**2.d0
!write(6,fmt="(f14.6,1x,f14.6)")xm(i),v(i,j)*1E-10
	enddo
	enddo
!read(5,*)debug
	ekin2=ekin/facte


	  if(icm.eq.1)then
	  do j=1,3
	  rcm(j)=0.d0
	  do i=1,natom
	  rcm(j)=rcm(j)+xm(i)*r(i,j)
	  enddo
	  rcm(j)=rcm(j)/xmassa
	  enddo
	  do i=1,natom
	  do j=1,3
	  rant(i,j)=r(i,j)-rcm(j)+rbox2
	  enddo
	  enddo
	  else
	  do i=1,natom
	  do j=1,3
	  rant(i,j)=r(i,j)
	  enddo
	  enddo
	  endif



!C suma l'energia potencial de la conformacio inicial
	epothb=0.d0
	epot=0.d0
	do i=1,natom-1
!write(6,fmt="(i4,a)", advance="no")i,": "
	do j=i+1,natom
	rij1=dbox(j,i,1)
	rij2=dbox(j,i,2)
	rij3=dbox(j,i,3)
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	dist=(rmod2**0.5d0)
!if ((i.eq.12).and.(j.eq.471))write(6,fmt="(i4,1x,f12.4,1x)", advance="no")inter(i,j),dist
	if(inter(i,j).eq.1)then
!write(6,fmt="(i4,1x)", advance="no")j
!write(6,fmt="(f12.6,1x)", advance="no")dist
	  k=nstep(i,j)
!write(6,*)k
!write(6,fmt="(a,i4,1x)",advance="no")"j=",j-1
!if (dist.lt.rstep(i,j,k)) then
!write(6,fmt="(a)",advance="no")"1"
!else
!write(6,fmt="(a)",advance="no")"0"
!endif
!write(6,fmt="(a)",advance="no")"&"
!if (k.gt.0) then
!write(6,fmt="(a)",advance="no")"1 "
!else
!write(6,fmt="(a)",advance="no")"0 "
!endif

	  do while(dist.lt.rstep(i,j,k).and.k.gt.0)
!write(6,fmt="(f14.3,1x)", advance="no")epot
	  epot=epot-estep(i,j,k)
	  if(imol(i).ne.imol(j))epotmol=epotmol-estep(i,j,k)
	  if(istruct(i,j).eq.1)then
	  epothb=epothb-estep(i,j,k)
	  if(imol(i).ne.imol(j))epothbmol=epothbmol-estep(i,j,k)
	  endif
	  k=k-1
	  enddo
	endif
	enddo
!write(6,*)''
	enddo
!if (ibloc.eq.1)read(5,*)debug	
!write(6,*)epot,ekin2

!do i=1,natom-1
!write(6,fmt="(i4,a)", advance="no")i-1,"="
!do j=i+1,natom
!if (inter(i,j).eq.1)write(6,fmt="(a)", advance="no")'x'
!enddo
!write(6,*)''
!enddo
!read(5,*)debug

!write(6,*)epotmol


	etot=epot+ekin2
!write(6,*)epot,etot
!if(ibloc.eq.2)call imprimir_v(natom,1)
	if(iprint.eq.1)write(12,*)temps,epot,epotmol,epothb,epothbmol,ekin2,etot,nhb

	icounter_1 = icounter_1 + 1

  40	continue

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


	write(20,333)ibloc
	write(21,333)ibloc
	do i=1,natom
	c2=atom(i)
	write(20,1000)c1,i,c2,res(i),cad(i),ind1(i),(r(i,j),j=1,3)
	if(atom(i).eq.'CA')write(21,1000)c1,i,c2,res(i),cad(i),ind1(i),(rant(i,j),j=1,3)
	enddo
	write(20,334)
	write(21,334)

	
	write(6,*)'Temps',temps,' hbonds',nhb,' epot',epot,etot
	write(8,*)'Temps',temps,' hbonds',nhb,' epot',epot,etot

 200	continue

!C	close(10)
	close(11)
	close(12)
	close(20)
	close(8)

	ekin=ekin/facte


        open(unit=19,file=file19)
	do i=1,natom
	write(19,1000)c1,i,atom(i),res(i),cad(i),ind1(i),(r(i,j),j=1,3)
	enddo
	close(19)
!C

       stop
       end



	subroutine potencial(natom)
	implicit real*8(a-h)
	implicit real*8(o-z)
	integer*2 inter,istruct,icov
	parameter(natmax=5000)
	parameter(npairmax=40000)
	character*5 atom,res
	common/xoc/r(natmax,3),v(natmax,3),xm(natmax),rbox,ierr
	common/pous/rstep(natmax,natmax,2),estep(natmax,natmax,2)
	common/intr/nstep(natmax,natmax),istruct(natmax,natmax),inter(natmax,natmax)
	common/cov/icov(natmax,natmax),rbound(npairmax),ibound(npairmax,2),rhc(natmax),sigma
	common/pdb/atom(natmax),res(natmax),ind2(natmax),nat(natmax),imol(natmax)
	common/fisic/evdw(natmax),rvdw(natmax),qq(natmax),gfree(natmax),vol(natmax)
	common/param/fvdw,fsolv,eps,xlamb
	common/atpres/ ihb(natmax),ica(natmax),io(natmax),ih(natmax),ico(natmax),in(natmax),icb(natmax)
	common/parmsolv/rsolv,asolv,bsolv,dwat,icont(natmax),fcont(natmax)
	
	rsolv2=rsolv*rsolv
	rbox2=0.5d0*rbox
	nres=ind2(natom)
	dwatd=dwat/(3.d0**0.5d0)
	

	do i=1,natom
	r1x=r(i,1)+dwat
	r2x=r(i,1)-dwat
	r3x=r(i,1)
	r4x=r(i,1)
	r5x=r(i,1)
	r6x=r(i,1)
	r1y=r(i,2)
	r2y=r(i,2)
	r3y=r(i,2)+dwat
	r4y=r(i,2)-dwat
	r5y=r(i,2)
	r6y=r(i,2)
	r1z=r(i,3)
	r2z=r(i,3)
	r3z=r(i,3)
	r4z=r(i,3)
	r5z=r(i,3)+dwat
	r6z=r(i,3)-dwat
	rv1x=r(i,1)+dwatd
	rv2x=r(i,1)+dwatd
	rv3x=r(i,1)+dwatd
	rv4x=r(i,1)+dwatd
	rv5x=r(i,1)-dwatd
	rv6x=r(i,1)-dwatd
	rv7x=r(i,1)-dwatd
	rv8x=r(i,1)-dwatd
	rv1y=r(i,2)+dwatd
	rv2y=r(i,2)+dwatd
	rv3y=r(i,2)-dwatd
	rv4y=r(i,2)-dwatd
	rv5y=r(i,2)+dwatd
	rv6y=r(i,2)+dwatd
	rv7y=r(i,2)-dwatd
	rv8y=r(i,2)-dwatd
	rv1z=r(i,3)+dwatd
	rv2z=r(i,3)-dwatd
	rv3z=r(i,3)+dwatd
	rv4z=r(i,3)-dwatd
	rv5z=r(i,3)+dwatd
	rv6z=r(i,3)-dwatd
	rv7z=r(i,3)+dwatd
	rv8z=r(i,3)-dwatd
	icont(i)=0
	do j=1,natom
	if(imol(j).ne.imol(i))cycle
	if(j.eq.i)cycle
	rij1=r1x-r(j,1)
	  if(rij1.gt.rbox2)then
	  rij1=rij1-rbox
	  elseif(rij1.lt.-rbox2)then
	  rij1=rij1+rbox
	  endif
	rij2=r1y-r(j,2)
	  if(rij2.gt.rbox2)then
	  rij2=rij2-rbox
	  elseif(rij2.lt.-rbox2)then
	  rij2=rij2+rbox
	  endif
	rij3=r1z-r(j,3)
	  if(rij3.gt.rbox2)then
	  rij3=rij3-rbox
	  elseif(rij3.lt.-rbox2)then
	  rij3=rij3+rbox
	  endif
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	if(rmod2.lt.rsolv2)then
	icont(i)=icont(i)+1
	exit
	endif
	enddo
	do j=1,natom
	if(imol(j).ne.imol(i))cycle
	if(j.eq.i)cycle
	rij1=r2x-r(j,1)
	  if(rij1.gt.rbox2)then
	  rij1=rij1-rbox
	  elseif(rij1.lt.-rbox2)then
	  rij1=rij1+rbox
	  endif
	rij2=r2y-r(j,2)
	  if(rij2.gt.rbox2)then
	  rij2=rij2-rbox
	  elseif(rij2.lt.-rbox2)then
	  rij2=rij2+rbox
	  endif
	rij3=r2z-r(j,3)
	  if(rij3.gt.rbox2)then
	  rij3=rij3-rbox
	  elseif(rij3.lt.-rbox2)then
	  rij3=rij3+rbox
	  endif
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	if(rmod2.lt.rsolv2)then
	icont(i)=icont(i)+1
	exit
	endif
	enddo
	do j=1,natom
	if(imol(j).ne.imol(i))cycle
	if(j.eq.i)cycle
	rij1=r3x-r(j,1)
	  if(rij1.gt.rbox2)then
	  rij1=rij1-rbox
	  elseif(rij1.lt.-rbox2)then
	  rij1=rij1+rbox
	  endif
	rij2=r3y-r(j,2)
	  if(rij2.gt.rbox2)then
	  rij2=rij2-rbox
	  elseif(rij2.lt.-rbox2)then
	  rij2=rij2+rbox
	  endif
	rij3=r3z-r(j,3)
	  if(rij3.gt.rbox2)then
	  rij3=rij3-rbox
	  elseif(rij3.lt.-rbox2)then
	  rij3=rij3+rbox
	  endif
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	if(rmod2.lt.rsolv2)then
	icont(i)=icont(i)+1
	exit
	endif
	enddo
	do j=1,natom
	if(imol(j).ne.imol(i))cycle
	if(j.eq.i)cycle
	rij1=r4x-r(j,1)
	  if(rij1.gt.rbox2)then
	  rij1=rij1-rbox
	  elseif(rij1.lt.-rbox2)then
	  rij1=rij1+rbox
	  endif
	rij2=r4y-r(j,2)
	  if(rij2.gt.rbox2)then
	  rij2=rij2-rbox
	  elseif(rij2.lt.-rbox2)then
	  rij2=rij2+rbox
	  endif
	rij3=r4z-r(j,3)
	  if(rij3.gt.rbox2)then
	  rij3=rij3-rbox
	  elseif(rij3.lt.-rbox2)then
	  rij3=rij3+rbox
	  endif
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	if(rmod2.lt.rsolv2)then
	icont(i)=icont(i)+1
	exit
	endif
	enddo
	do j=1,natom
	if(imol(j).ne.imol(i))cycle
	if(j.eq.i)cycle
	rij1=r5x-r(j,1)
	  if(rij1.gt.rbox2)then
	  rij1=rij1-rbox
	  elseif(rij1.lt.-rbox2)then
	  rij1=rij1+rbox
	  endif
	rij2=r5y-r(j,2)
	  if(rij2.gt.rbox2)then
	  rij2=rij2-rbox
	  elseif(rij2.lt.-rbox2)then
	  rij2=rij2+rbox
	  endif
	rij3=r5z-r(j,3)
	  if(rij3.gt.rbox2)then
	  rij3=rij3-rbox
	  elseif(rij3.lt.-rbox2)then
	  rij3=rij3+rbox
	  endif
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	if(rmod2.lt.rsolv2)then
	icont(i)=icont(i)+1
	exit
	endif
	enddo
	do j=1,natom
	if(imol(j).ne.imol(i))cycle
	if(j.eq.i)cycle
	rij1=r6x-r(j,1)
	  if(rij1.gt.rbox2)then
	  rij1=rij1-rbox
	  elseif(rij1.lt.-rbox2)then
	  rij1=rij1+rbox
	  endif
	rij2=r6y-r(j,2)
	  if(rij2.gt.rbox2)then
	  rij2=rij2-rbox
	  elseif(rij2.lt.-rbox2)then
	  rij2=rij2+rbox
	  endif
	rij3=r6z-r(j,3)
	  if(rij3.gt.rbox2)then
	  rij3=rij3-rbox
	  elseif(rij3.lt.-rbox2)then
	  rij3=rij3+rbox
	  endif
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	if(rmod2.lt.rsolv2)then
	icont(i)=icont(i)+1
	exit
	endif
	enddo
!C vertexs
	do j=1,natom
	if(imol(j).ne.imol(i))cycle
	if(j.eq.i)cycle
	rij1=rv1x-r(j,1)
	  if(rij1.gt.rbox2)then
	  rij1=rij1-rbox
	  elseif(rij1.lt.-rbox2)then
	  rij1=rij1+rbox
	  endif
	rij2=rv1y-r(j,2)
	  if(rij2.gt.rbox2)then
	  rij2=rij2-rbox
	  elseif(rij2.lt.-rbox2)then
	  rij2=rij2+rbox
	  endif
	rij3=rv1z-r(j,3)
	  if(rij3.gt.rbox2)then
	  rij3=rij3-rbox
	  elseif(rij3.lt.-rbox2)then
	  rij3=rij3+rbox
	  endif
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	if(rmod2.lt.rsolv2)then
	icont(i)=icont(i)+1
	exit
	endif
	enddo
	do j=1,natom
	if(imol(j).ne.imol(i))cycle
	if(j.eq.i)cycle
	rij1=rv2x-r(j,1)
	  if(rij1.gt.rbox2)then
	  rij1=rij1-rbox
	  elseif(rij1.lt.-rbox2)then
	  rij1=rij1+rbox
	  endif
	rij2=rv2y-r(j,2)
	  if(rij2.gt.rbox2)then
	  rij2=rij2-rbox
	  elseif(rij2.lt.-rbox2)then
	  rij2=rij2+rbox
	  endif
	rij3=rv2z-r(j,3)
	  if(rij3.gt.rbox2)then
	  rij3=rij3-rbox
	  elseif(rij3.lt.-rbox2)then
	  rij3=rij3+rbox
	  endif
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	if(rmod2.lt.rsolv2)then
	icont(i)=icont(i)+1
	exit
	endif
	enddo
	do j=1,natom
	if(imol(j).ne.imol(i))cycle
	if(j.eq.i)cycle
	rij1=rv3x-r(j,1)
	  if(rij1.gt.rbox2)then
	  rij1=rij1-rbox
	  elseif(rij1.lt.-rbox2)then
	  rij1=rij1+rbox
	  endif
	rij2=rv3y-r(j,2)
	  if(rij2.gt.rbox2)then
	  rij2=rij2-rbox
	  elseif(rij2.lt.-rbox2)then
	  rij2=rij2+rbox
	  endif
	rij3=rv3z-r(j,3)
	  if(rij3.gt.rbox2)then
	  rij3=rij3-rbox
	  elseif(rij3.lt.-rbox2)then
	  rij3=rij3+rbox
	  endif
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	if(rmod2.lt.rsolv2)then
	icont(i)=icont(i)+1
	exit
	endif
	enddo
	do j=1,natom
	if(imol(j).ne.imol(i))cycle
	if(j.eq.i)cycle
	rij1=rv4x-r(j,1)
	  if(rij1.gt.rbox2)then
	  rij1=rij1-rbox
	  elseif(rij1.lt.-rbox2)then
	  rij1=rij1+rbox
	  endif
	rij2=rv4y-r(j,2)
	  if(rij2.gt.rbox2)then
	  rij2=rij2-rbox
	  elseif(rij2.lt.-rbox2)then
	  rij2=rij2+rbox
	  endif
	rij3=rv4z-r(j,3)
	  if(rij3.gt.rbox2)then
	  rij3=rij3-rbox
	  elseif(rij3.lt.-rbox2)then
	  rij3=rij3+rbox
	  endif
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	if(rmod2.lt.rsolv2)then
	icont(i)=icont(i)+1
	exit
	endif
	enddo
	do j=1,natom
	if(imol(j).ne.imol(i))cycle
	if(j.eq.i)cycle
	rij1=rv5x-r(j,1)
	  if(rij1.gt.rbox2)then
	  rij1=rij1-rbox
	  elseif(rij1.lt.-rbox2)then
	  rij1=rij1+rbox
	  endif
	rij2=rv5y-r(j,2)
	  if(rij2.gt.rbox2)then
	  rij2=rij2-rbox
	  elseif(rij2.lt.-rbox2)then
	  rij2=rij2+rbox
	  endif
	rij3=rv5z-r(j,3)
	  if(rij3.gt.rbox2)then
	  rij3=rij3-rbox
	  elseif(rij3.lt.-rbox2)then
	  rij3=rij3+rbox
	  endif
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	if(rmod2.lt.rsolv2)then
	icont(i)=icont(i)+1
	exit
	endif
	enddo
	do j=1,natom
	if(imol(j).ne.imol(i))cycle
	if(j.eq.i)cycle
	rij1=rv6x-r(j,1)
	  if(rij1.gt.rbox2)then
	  rij1=rij1-rbox
	  elseif(rij1.lt.-rbox2)then
	  rij1=rij1+rbox
	  endif
	rij2=rv6y-r(j,2)
	  if(rij2.gt.rbox2)then
	  rij2=rij2-rbox
	  elseif(rij2.lt.-rbox2)then
	  rij2=rij2+rbox
	  endif
	rij3=rv6z-r(j,3)
	  if(rij3.gt.rbox2)then
	  rij3=rij3-rbox
	  elseif(rij3.lt.-rbox2)then
	  rij3=rij3+rbox
	  endif
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	if(rmod2.lt.rsolv2)then
	icont(i)=icont(i)+1
	exit
	endif
	enddo
	do j=1,natom
	if(imol(j).ne.imol(i))cycle
	if(j.eq.i)cycle
	rij1=rv7x-r(j,1)
	  if(rij1.gt.rbox2)then
	  rij1=rij1-rbox
	  elseif(rij1.lt.-rbox2)then
	  rij1=rij1+rbox
	  endif
	rij2=rv7y-r(j,2)
	  if(rij2.gt.rbox2)then
	  rij2=rij2-rbox
	  elseif(rij2.lt.-rbox2)then
	  rij2=rij2+rbox
	  endif
	rij3=rv7z-r(j,3)
	  if(rij3.gt.rbox2)then
	  rij3=rij3-rbox
	  elseif(rij3.lt.-rbox2)then
	  rij3=rij3+rbox
	  endif
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	if(rmod2.lt.rsolv2)then
	icont(i)=icont(i)+1
	exit
	endif
	enddo
	do j=1,natom
	if(imol(j).ne.imol(i))cycle
	if(j.eq.i)cycle
	rij1=rv8x-r(j,1)
	  if(rij1.gt.rbox2)then
	  rij1=rij1-rbox
	  elseif(rij1.lt.-rbox2)then
	  rij1=rij1+rbox
	  endif
	rij2=rv8y-r(j,2)
	  if(rij2.gt.rbox2)then
	  rij2=rij2-rbox
	  elseif(rij2.lt.-rbox2)then
	  rij2=rij2+rbox
	  endif
	rij3=rv8z-r(j,3)
	  if(rij3.gt.rbox2)then
	  rij3=rij3-rbox
	  elseif(rij3.lt.-rbox2)then
	  rij3=rij3+rbox
	  endif
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	if(rmod2.lt.rsolv2)then
	icont(i)=icont(i)+1
	exit
	endif
	enddo

	enddo
	
	do i=1,natom
	fcont(i)=1.d0/(1.d0+exp((icont(i)-asolv)/bsolv))
!write(6,fmt="(f12.6)")fcont(i)
	enddo
!read(5,*)debug
	
	
	do i=1,natom-1
!write(6,fmt="(i4,a)", advance="no")i,': '
	ii=ind2(i)
	do j=i+1,natom
	jj=ind2(j)
	if(icov(i,j).eq.0)then
	icount_1 = icount_1 + 1
	if(istruct(i,j).eq.0)then
!write(6,fmt="(i4,1x)", advance="no")j
	icount_2 = icount_2 + 1
	rvdwij=rvdw(i)+rvdw(j)
!if (i.eq.1000)write(6,fmt="(f12.6)")rvdwij
	sto=(2.d0/(nat(i)**0.33d0+nat(j)**0.33d0))**6.d0
!if (i.eq.1000)write(6,fmt="(f12.6)")sto
	potvdw=((evdw(i)*evdw(j))**0.5d0)*sto*(sto-2.d0)
!if (i.eq.1000)write(6,fmt="(f12.6)")potvdw
	potlk=-0.09d0/xlamb*(gfree(i)*vol(j)+gfree(j)*vol(i))/(rvdwij**2.d0*exp((rvdwij/xlamb)**2.d0))
!if (i.eq.1000)write(6,fmt="(f12.6)")potlk
	eij=fvdw*potvdw+fsolv*potlk*fcont(i)*fcont(j)+eps*qq(i)*qq(j)/rvdwij
!if (i.eq.1134)write(6,fmt="(f14.8)")eij
	nstep(i,j)=2
	   rstep(i,j,1)=0.9d0*rvdwij
	   rstep(i,j,2)=1.1d0*rvdwij
	 if(eij.lt.0.d0)then
	   estep(i,j,1)=3.d0*eij
	   estep(i,j,2)=-eij
	  else
	   estep(i,j,1)=-eij
	   estep(i,j,2)=-eij
	 endif
!write(6,fmt="(f14.8,a,f14.8,1x)", advance="no")rstep(i,j,1),'-',estep(i,j,1)
	endif
	endif
	enddo
!write(6,*)''
	enddo
!read(5,*)debug
	return
	end






	subroutine enchufa(natom,dcut)
	implicit real*8(a-h)
	implicit real*8(o-z)
	integer*2 inter,istruct,icov
	parameter(natmax=5000)
	parameter(npairmax=40000)
	character*5 atom,res
	common/xoc/r(natmax,3),v(natmax,3),xm(natmax),rbox,ierr
	common/intr/nstep(natmax,natmax),istruct(natmax,natmax),inter(natmax,natmax)
	common/cov/icov(natmax,natmax),rbound(npairmax),ibound(npairmax,2),rhc(natmax),sigma
	common/atpres/ ihb(natmax),ica(natmax),io(natmax),ih(natmax),ico(natmax),in(natmax),icb(natmax)
	common/pdb/atom(natmax),res(natmax),ind2(natmax),nat(natmax),imol(natmax)
	common/npt/ipot(natmax),npot(natmax,natmax)
	dcut2=dcut*dcut
	do i=1,natom-1
!write(6,fmt="(i4,1x)",advance="no")i
!write(6,*)i,' ',ipot(i)
!C	do j=i+1,natom
	do l=1,ipot(i)
	j=npot(i,l)
!if ((i.eq.12).and.(j.eq.471))write(6,fmt="(i4,1x,f12.4,1x,f12.4,1x,f12.4)")inter(i,j),dbox(i,j,1),dbox(i,j,2),dbox(i,j,3)
!C	if(istruct(i,j).eq.0)then
!C	if(icov(i,j).eq.0)then
!if ((i.eq.12).and.(j.eq.471))write(6,*)nat(i),' ',nat(j)
	if(nat(i).gt.1.and.nat(j).gt.1)then
!write(6,fmt="(i4,1x)",advance="no")j
!C	  inter(i,j)=0
!C mira si hi ha definits potencials d'estructura
	  rij1=dbox(i,j,1)
	  rij2=dbox(i,j,2)
	  rij3=dbox(i,j,3)
	  rmod2=rij1*rij1+rij2*rij2+rij3*rij3
!write(6,fmt="(f10.5,1x)",advance="no")rmod2
	  if(rmod2.lt.dcut2) then
	  inter(i,j)=1
!if ((i.eq.12).and.(j.eq.471))write(6,fmt="(i4,1x,f12.4,1x,f12.4,1x,f12.4,2x,f12.4,1x,f12.4)")inter(i,j),dbox(i,j,1),dbox(i,j,2),dbox(i,j,3),rmod2,dcut2
!write(6,fmt="(i4,1x)",advance="no")j
	  endif
	endif
!C	endif
!C	endif
	enddo
!write(6,*)''
	enddo
	
	nres=ind2(natom)
!write(6,*)'inicio'
	do i=1,nres-1
	n1=in(i)
	n2=in(i+1)
!write(6,*)n1-1,' ',n2-1
	inter(n1,n2)=0
	enddo
	do i=2,nres
	n1=ico(i-1)
	n2=ico(i)
!write(6,*)n1-1,' ',n2-1
	inter(n1,n2)=0
	enddo

	return
	end



	subroutine dmdshake(natom,nbound,mem1,mem2)
	implicit real*8(a-h)
	implicit real*8(o-z)
	integer*2 inter,istruct,icov
	parameter(natmax=5000,npairmax=40000,facte=4186.d0)
	character*5 atom,res
	common/xoc/r(natmax,3),v(natmax,3),xm(natmax),rbox,ierr
	common/cov/icov(natmax,natmax),rbound(npairmax),ibound(npairmax,2),rhc(natmax),sigma
	common/shake/ishk(natmax),nshk(natmax,natmax)

	ierr=0
!write(6,*)'-->',v(38,1)*1E-10
!C particules NO enllaçades
	do i=1,natom-1
!write(6,*)ishk(i)
!C	if(i.eq.mem1.or.i.eq.mem2)cycle
	do l=1,ishk(i)
	j=nshk(i,l)
!C	if(j.eq.mem1.or.j.eq.mem2)cycle
	rij1=dbox(j,i,1)
	rij2=dbox(j,i,2)
	rij3=dbox(j,i,3)
        rij=(rij1*rij1+rij2*rij2+rij3*rij3)**0.5d0
        vij1=v(i,1)-v(j,1)
        vij2=v(i,2)-v(j,2)
        vij3=v(i,3)-v(j,3)
        prod=rij1*vij1+rij2*vij2+rij3*vij3
	dmin=rhc(i)+rhc(j)
!write(6,fmt="(f40.35)")rij
!if (i.eq.569)write(6,*)'=>',v(569,1)*1E-10
!if (j.eq.38)write(6,*)'->',v(38,1)*1E-10
!C xoc frontal entre particules no enllaçades
	if(rij.lt.dmin.and.prod.lt.0.0d0)then
!write(6,fmt="(f40.35,1x,f40.35,1x,f40.35,1x,f40.35,1x,f40.35,1x,f40.35)")v(i,1)*1d-11,v(i,2)*1d-11,v(i,3)*1d-11,v(j,1)*1d-11,v(j,2)*1d-11,v(j,3)*1d-11
	call chgmom(i,j,rij1,rij2,rij3)
	ierr=ierr+1
!write(6,fmt="(f40.35,1x,f40.35,1x,f40.35,1x,f40.35,1x,f40.35,1x,f40.35)")v(i,1)*1d-11,v(i,2)*1d-11,v(i,3)*1d-11,v(j,1)*1d-11,v(j,2)*1d-11,v(j,3)*1d-11
	endif
!if (i.eq.2) read(5,*)debug
!if (i.eq.51) call imprimir_v(natom,1)
!if (i.eq.569)write(6,*)'<=',v(569,1)*1E-10
!if (j.eq.38)write(6,*)'<-',v(38,1)*1E-10
	enddo
	enddo

!write(6,*)'==>',v(38,1)*1E-10
!write(6,*)ierr
!read(5,*)debug

!C particules enllaçades
	do k=1,nbound
	i=ibound(k,1)
	j=ibound(k,2)
!C	if(i.eq.mem1.or.i.eq.mem2)cycle
!C	if(j.eq.mem1.or.j.eq.mem2)cycle
	rbmin=rbound(k)*(1.d0-sigma)
	rbmax=rbound(k)*(1.d0+sigma)
	rbmin2=rbmin*rbmin
	rbmax2=rbmax*rbmax
	rij1=dbox(j,i,1)
	rij2=dbox(j,i,2)
	rij3=dbox(j,i,3)
        rmod2=rij1*rij1+rij2*rij2+rij3*rij3
        vij1=v(i,1)-v(j,1)
        vij2=v(i,2)-v(j,2)
        vij3=v(i,3)-v(j,3)
        prod=rij1*vij1+rij2*vij2+rij3*vij3
!if(((i.gt.36).and.(i.lt.47)).or.(i.eq.587).or.(j.eq.587))write(6,fmt="(i4,1x,i4,1x,f12.4,1x,f12.4,1x,f12.4,1x,f12.4)")i-1,j-1,rmod2,prod*1E-10,rbmin2,rbmax2
!if(((i.gt.36).and.(i.lt.47)).or.(i.eq.587).or.(j.eq.587))write(6,fmt="(f12.4,1x,f12.4,1x,f12.4,2x,f12.4,1x,f12.4,1x,f12.4)")v(i,1)*1E-10,v(i,2)*1E-10,v(i,3)*1E-10,v(j,1)*1E-10,v(j,2)*1E-10,v(j,3)*1E-10
	if(rmod2.gt.rbmax2.and.prod.gt.0.d0)then
	ierr=ierr+1
!if((i.eq.38).or.(j.eq.38))write(6,*)i-1,' ',j-1,' ||> ',v(38,1)*1E-10
	call chgmom(i,j,rij1,rij2,rij3)
!if((i.eq.38).or.(j.eq.38))write(6,*)i-1,' ',j-1,' ++> ',v(38,1)*1E-10
	else
	  if(rmod2.lt.rbmin2.and.prod.lt.0.d0)then
	  ierr=ierr+1
!write(6,*)'))>',i-1,' ',j-1,' ',v(38,1)*1E-10
	  call chgmom(i,j,rij1,rij2,rij3)
!write(6,*)'((>',i-1,' ',j-1,' ',v(38,1)*1E-10
	  endif
	endif
	enddo
	
!write(6,*)nbound
!read(5,*)debug

	return
	end





	subroutine chgmom(mem1,mem2,rij1,rij2,rij3)
	implicit real*8(a-h)
	implicit real*8(o-z)
	parameter(natmax=5000)
	common/xoc/r(natmax,3),v(natmax,3),xm(natmax),rbox,ierr
!if((mem1.eq.38).or.(mem2.eq.38))write(6,*)mem1-1,' ',mem2-1,' => ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.587).or.(mem2.eq.587))write(6,*)mem1-1,' ',mem2-1,' => ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.569).or.(mem2.eq.569))write(6,*)mem1-1,' ',mem2-1,' => ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.12).or.(mem2.eq.12))write(6,*)mem1-1,' ',mem2-1,' => ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.471).or.(mem2.eq.471))write(6,*)mem1-1,' ',mem2-1,' => ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.15).or.(mem2.eq.15))write(6,*)mem1-1,' ',mem2-1,' => ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
	vdmod=0.d0
	vdmod=vdmod+(v(mem2,1)-v(mem1,1))*rij1
	vdmod=vdmod+(v(mem2,2)-v(mem1,2))*rij2
	vdmod=vdmod+(v(mem2,3)-v(mem1,3))*rij3
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
	vdmod=vdmod/rmod2
	xsum=0.5d0*(1.d0/xm(mem1)+1.d0/xm(mem2))
!C modul del moment transferit en la colisio
	dp=vdmod/xsum
        v(mem1,1)=v(mem1,1)+dp/xm(mem1)*rij1
        v(mem1,2)=v(mem1,2)+dp/xm(mem1)*rij2
        v(mem1,3)=v(mem1,3)+dp/xm(mem1)*rij3
        v(mem2,1)=v(mem2,1)-dp/xm(mem2)*rij1
        v(mem2,2)=v(mem2,2)-dp/xm(mem2)*rij2
        v(mem2,3)=v(mem2,3)-dp/xm(mem2)*rij3
!write(6,fmt="(f45.35,1x,f45.35,1x,f45.35,1x,f45.35)")vdmod*1d-10,rmod2,xsum,dp*1d-10
!if((mem1.eq.38).or.(mem2.eq.38))write(6,*)mem1-1,' ',mem2-1,' -> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.587).or.(mem2.eq.587))write(6,*)mem1-1,' ',mem2-1,' -> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.569).or.(mem2.eq.569))write(6,*)mem1-1,' ',mem2-1,' -> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.12).or.(mem2.eq.12))write(6,*)mem1-1,' ',mem2-1,' -> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.471).or.(mem2.eq.471))write(6,*)mem1-1,' ',mem2-1,' -> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.15).or.(mem2.eq.15))write(6,*)mem1-1,' ',mem2-1,' -> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!write(6,*)vdmod,rmod2,xsum,dp
!write(6,*)v(mem1,1), v(mem1,2), v(mem1,3), v(mem2,1),v(mem2,2),v(mem2,3)
!read(5,*)debug
		
	return
	end

	subroutine chgmomene(mem1,mem2,rij1,rij2,rij3,dpot,ich)
	implicit real*8(a-h)
	implicit real*8(o-z)
	parameter(natmax=5000)
	common/xoc/r(natmax,3),v(natmax,3),xm(natmax),rbox,ierr
!if((mem1.eq.38).or.(mem2.eq.38))write(6,*)mem1-1,' ',mem2-1,' /\> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.587).or.(mem2.eq.587))write(6,*)mem1-1,' ',mem2-1,' /\> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.569).or.(mem2.eq.569))write(6,*)mem1-1,' ',mem2-1,' /\> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.12).or.(mem2.eq.12))write(6,*)mem1-1,' ',mem2-1,' /\> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.471).or.(mem2.eq.471))write(6,*)mem1-1,' ',mem2-1,' /\> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.15).or.(mem2.eq.15))write(6,*)mem1-1,' ',mem2-1,' /\> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
	a=1.d-10
	vdmod=0
	vdmod=vdmod+(v(mem2,1)-v(mem1,1))*rij1
	vdmod=vdmod+(v(mem2,2)-v(mem1,2))*rij2
	vdmod=vdmod+(v(mem2,3)-v(mem1,3))*rij3
	rmod2=rij1*rij1+rij2*rij2+rij3*rij3
!C projeccio del moment en l'eix que uneix les dues particules
	vdmod=vdmod/rmod2
	xsum=0.5d0*(1.d0/xm(mem1)+1.d0/xm(mem2))
!C modul del moment transferit/distancia en un xoc elastic
	dp=vdmod/xsum
	sto=dp*dp/4.d0-dpot/(rmod2*xsum*a*a)
!if((mem1.eq.1052).and.(mem2.eq.1077))write(6,fmt="(f45.40,1x,f45.40)")dp*1.0d-10,(sto**0.5d0)*1.0d-10
        if(sto.gt.0.d0)then
!C sempre es la resta dels dos valors absoluts
        if(vdmod.gt.0.d0)then
        dp=dp/2.d0-(sto**0.5d0)
        else
        dp=dp/2.d0+(sto**0.5d0)
        endif
        ich=1
        else
!C no traspassa la barrera
        ich=0
        endif
!if((mem1.eq.1052).and.(mem2.eq.1077))write(6,fmt="(f45.40)")dp*1.0d-10
!if((mem1.eq.1052).and.(mem2.eq.1077))write(6,fmt="(f45.40,1x,f45.40,1x,f45.40,2x,f45.40,1x,f45.40)")rmod2,vdmod*1.0d-10,sto*1.0d-20,xm(mem1),xm(mem2)
		v(mem1,1)=v(mem1,1)+(dp/xm(mem1))*rij1
        v(mem1,2)=v(mem1,2)+(dp/xm(mem1))*rij2
        v(mem1,3)=v(mem1,3)+(dp/xm(mem1))*rij3
        v(mem2,1)=v(mem2,1)-(dp/xm(mem2))*rij1
        v(mem2,2)=v(mem2,2)-(dp/xm(mem2))*rij2
        v(mem2,3)=v(mem2,3)-(dp/xm(mem2))*rij3
!if((mem1.eq.38).or.(mem2.eq.38))write(6,*)mem1-1,' ',mem2-1,' \/-> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.587).or.(mem2.eq.587))write(6,*)mem1-1,' ',mem2-1,' \/-> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.569).or.(mem2.eq.569))write(6,*)mem1-1,' ',mem2-1,' \/> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.12).or.(mem2.eq.12))write(6,*)mem1-1,' ',mem2-1,' \/> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.471).or.(mem2.eq.471))write(6,*)mem1-1,' ',mem2-1,' \/> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
!if((mem1.eq.15).or.(mem2.eq.15))write(6,*)mem1-1,' ',mem2-1,' \/> ',v(mem1,1)*1E-10,v(mem2,1)*1E-10
	return
	end


	
	
	subroutine creapouhb(n1,n2,rmin,r0,r1,rmax,ehb)
	implicit real*8(a-h)
	implicit real*8(o-z)
	integer*2 inter,istruct
	parameter(natmax=5000)
	common/pous/rstep(natmax,natmax,2),estep(natmax,natmax,2)
	common/intr/nstep(natmax,natmax),istruct(natmax,natmax),inter(natmax,natmax)
	inter(n1,n2)=1
	istruct(n1,n2)=1
	nstep(n1,n2)=2
	rstep(n1,n2,1)=rmin
	rstep(n1,n2,2)=rmax
	estep(n1,n2,1)=-1.5d0*ehb
	estep(n1,n2,2)=1.5d0*ehb
	return
	end


	function dbox(n1,n2,k)
	implicit real*8(a-h)
	implicit real*8(o-z)
	parameter(natmax=5000)
	common/xoc/r(natmax,3),v(natmax,3),xm(natmax),rbox,ierr
	rbox2=0.5d0*rbox
	r12=r(n2,k)-r(n1,k)
	if(r12.gt.rbox2)then
	r12=r12-rbox
	elseif(r12.lt.-rbox2)then
	r12=r12+rbox
	endif
	dbox=r12
	end
	
	subroutine random_number_new(ra)
	implicit real*8(a-h)
	implicit real*8(o-z)
	common/random_c/pseudo

		pseudo = pseudo + 0.0625d0
		
		if(pseudo > 0.99d0)pseudo = 0.0625d0
		
		ra = pseudo

	end subroutine random_number_new

      subroutine rnd_gauss ( fi, xm, T)
	implicit real*8(a-h)
	implicit real*8(o-z)
	R=8.314d0
        std_dev = ((T*R) / xm)**0.5d0
        pi = 4.0d0 * atan ( 1.0d0 )
        call random_number_new ( rnd1 )
        call random_number_new ( rnd2 )
        fi = std_dev * (( -2.0d0 * log ( rnd1 ) )**0.5d0) * cos ( 2.0d0 * pi * rnd2 )
      end subroutine rnd_gauss

	subroutine imprimir_inter(natom)
	implicit real*8(a-h)
	implicit real*8(o-z)
	integer*2 inter,istruct
	parameter(natmax=5000)
	common/intr/nstep(natmax,natmax),istruct(natmax,natmax),inter(natmax,natmax)
	
	do i=1,natom-1
		write(6,fmt="(i4,a)", advance="no")i-1,"="
		
		do j=i+1,natom
			if (inter(i,j).eq.1)write(6,fmt="(a)", advance="no")'x'
		enddo
		
		write(6,*)''
	enddo

	read(5,*)debug
	
	return
	end subroutine imprimir_inter
	
	subroutine imprimir_v(natom,nterminar)
	implicit real*8(a-h)
	implicit real*8(o-z)
	integer*2 inter,istruct
	parameter(natmax=5000)
	common/xoc/r(natmax,3),v(natmax,3),xm(natmax),rbox,ierr
	
	do i=1,natom
		write(6,fmt="(i4,1x,f22.16,1x,f22.16,1x,f22.16)")i,v(i,1)*1.d-10,v(i,2)*1.d-10,v(i,3)*1.d-10
	enddo

	if (nterminar.eq.1)read(5,*)debug
	
	return
	end subroutine imprimir_v
	
	subroutine imprimir_r(natom,nterminar)
	implicit real*8(a-h)
	implicit real*8(o-z)
	integer*2 inter,istruct
	parameter(natmax=5000)
	common/xoc/r(natmax,3),v(natmax,3),xm(natmax),rbox,ierr
	
	do i=1,natom
		write(6,fmt="(i4,1x,f45.40,1x,f45.40,1x,f45.40)")i,r(i,1)*1.d-10,r(i,2)*1.d-10,r(i,3)*1.d-10
	enddo

	if (nterminar.eq.1)read(5,*)debug
	
	return
	end subroutine imprimir_r