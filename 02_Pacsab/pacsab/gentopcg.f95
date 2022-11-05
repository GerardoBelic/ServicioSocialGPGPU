!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C PROGRAMA  TOPOLOGIACG
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	program topologiacg

!c	VARIABLES

	implicit real(a-h)
	implicit real*8(o-z)
	parameter(natmax=100000,nmax=10000)
       	dimension r(natmax,3),rnat(natmax,3),rant(natmax,3),v(natmax,3),igo(natmax),ipot(natmax),iq(natmax),ihf(natmax),irs(natmax),ipart(natmax),tpart(natmax),b(natmax,3)
	dimension ica(natmax),ico(natmax),in(natmax),icb(natmax),io(natmax),ih(natmax),ibead(natmax)
	dimension rhc(natmax),xm(natmax),ind1(natmax),ind2(natmax),inb1(natmax),qq(natmax)
	dimension gfree(natmax),lamb(natmax),evdw(natmax),ibn(natmax),ibca(natmax),ibco(natmax),ibh(natmax),ibo(natmax),xmb(natmax)
	dimension ipsb(nmax,nmax)
	dimension rcm(3),vcm(3),vd(3),vm1(3),vm2(3)
	dimension ethf(13,13)
	character*20 file7,file8,file9,file10,file11,file12,file13,file15,file16,file17,file19,file20
	character*5 c1,c2,c3,c4,part(natmax),bead(natmax)
	character*3 atom(natmax),res(natmax),cad(natmax),cbead(natmax),cadind(26)
	namelist /input/ tsnap,rco,rca,rn,rcb,ro,rs,sigma,sigmahb,temp,kkk,file7,file8,file9,file10,file13,file19,file20,nbloc,id,rcut,ixerra,rcutgo,sigmago,file12,tmin,irand,trnd,timpr,itots,dps,dhf,tref,ehb,eps,ego,ipotgo,factd,isolv,irig
	namelist /distancies/rnomax,rnomin,rncmax,rncmin,rcomax,rcomin
 1000	FORMAT (A4,3X,I4,2X,A3,1X,A3,3X,I3,4X,F8.3,F8.3,F8.3)
 2000	FORMAT (A4,3X,I4,2X,A3,1X,A3,1X,A1,I4,4X,F8.3,F8.3,F8.3)
 333	FORMAT('MODEL',8X,I4)
 334	FORMAT('ENDMDL')
 1015	FORMAT(E12.6,20(1X,F7.2))
	DATA (cadind(i),i=1,26)/'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/


!C default values
!C 	Archivos 
	file11='distancia.dat'
	file12='energia.dat'
	file20='snapshots.pdb'
	file13='velocitats.out'
	file7='topologia.dat'
	file10='res'
	file15='pontshidrogen.dat'
	file16='pontssalins.dat'
	file17='hidrofobiques.dat'
	file19='structurecg.pdb'

!C	Constantes

	pi=atan(1.d0)*4.0d0
	a=1.d-10
	tmin=1.d-22
	kkk=2381
	temp=1.
	nbloc=2
	rcut=10.
	rcutgo=8.
	rco=1.4
	rca=2.
	rcb=2.
	rs=1.5
	rn=1.2
	ro=1.2
	rh=0.7
	xmc=.012
	xmo=.016
	xmn=.014
	xms=.032
	xmh=.01
	xmres=.1
	sigma=.05
	sigmago=0.1
	ehb=10.
	dhf=6.
!C per defecte no considera interaccions electrostatiques
	eps=0.
	dps=6.
	ego=1.
	t0=0.
	tref=1.d-10
	trnd=1.d-10
	timpr=1.d-10
	id=0
	ires1=1
	ires2=2
	iread=0
	ixerra=0
	irand=1
	ipotgo=0
	factd=0.2
	isolv=1
	irig=0
	rnomax=3.4
	rnomin=2.7
	rncmax=4.6
	rncmin=3.7
	rcomax=4.3
	rcomin=3.4


	ehb=-ehb
	ia0=ia
	ib0=ib
	if(timpr.gt.tsnap)timpr=tsnap
	if(file10.eq.'res')file10=file9
!C converteixo la temperatura en energia per particula
!C utilitzo la constant de Boltzmann per mol
	temp=1.5*8.32*temp
	facte=4167.
	tempcal=temp/facte
	natom=0
	nres=0

!C llegeix les coordenades (fitxer pdb)

	nbead=0
	imol=0
	
	k0=0
	kk=0
	n=1
	nres0=0
	
 50	read(5,2000,end=51)c1,j,c2,c3,c4,k,x,y,z
	atom(n)=c2
	res(n)=c3
	if(c4.ne.cad(n-1))nres0=ind2(n-1)
	if(k.ne.k0)kk=kk+1
	cad(n)=c4
	ind2(n)=kk
	k0=k
	r(n,1)=x
	r(n,2)=y
	r(n,3)=z
	rnat(n,1)=x
	rnat(n,2)=y
	rnat(n,3)=z
	qq(n)=0
	irs(n)=0
	igo(n)=0
	lamb(n)=3.5
	gfree(n)=0.
	evdw(n)=0.
	ind1(n)=ind2(n)-nres0
	k=ind2(n)
	if(atom(n).eq.'N')in(k)=n
	if(atom(n).eq.'CA')ica(k)=n
	if(atom(n).eq.'C')ico(k)=n
	if(atom(n).eq.'O')io(k)=n
	if(atom(n).eq.'H')ih(k)=n
!C	write(6,2000)c1,n,atom(n),res(n),cad(n),ind1(n),x,y,z
	n=n+1
	goto 50
 51	natom=n-1
	imol=imol+1
!C	stop

	nres=ind2(natom)

	write(6,*)natom,nres

!C assigna radis de hardcore en funcio de l'element
	xmassa=0.
	do i=1,natom
	xm(i)=xmc
	rhc(i)=rco
	if(atom(i).eq.'CA')then
	rhc(i)=rca
	endif
	if(atom(i).eq.'CB')then
	rhc(i)=rcb
	endif
	if(atom(i)(1:1).eq.'N')then
	xm(i)=xmn
	rhc(i)=rn
	endif
	if(atom(i)(1:1).eq.'O')then
	xm(i)=xmo
	rhc(i)=ro
	endif
	if(atom(i)(1:1).eq.'S')then
	xm(i)=xms
	rhc(i)=rs
	endif
	if(atom(i).eq.'H')then
	xm(i)=xmh
	rhc(i)=rh
	endif
	xmassa=xmassa+xm(i)
	v(i,1)=0.
	v(i,2)=0.
	v(i,3)=0.
	enddo


	m=nbead+1

	do i=1,nres
!C primer bead corresponent al backbone
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
!C
	ii=in(i)
	b(m,1)=r(ii,1)
	b(m,2)=r(ii,2)
	b(m,3)=r(ii,3)
	part(m)='N'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	inb1(m)=ind1(ii)
	ibn(i)=m
	m=m+1
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
!C
	if(res(ii).ne.'PRO')then
	ii=ih(i)
	b(m,1)=r(ii,1)
	b(m,2)=r(ii,2)
	b(m,3)=r(ii,3)
	part(m)='H'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	inb1(m)=ind1(ii)
	ibh(i)=m
	m=m+1
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	endif
!C
	ii=ica(i)
	b(m,1)=r(ii,1)
	b(m,2)=r(ii,2)
	b(m,3)=r(ii,3)
	part(m)='CA'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	inb1(m)=ind1(ii)
	ibca(i)=m
	
	m=m+1
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
!C
	ii=ico(i)
	b(m,1)=r(ii,1)
	b(m,2)=r(ii,2)
	b(m,3)=r(ii,3)
	part(m)='C'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	inb1(m)=ind1(ii)
	ibco(i)=m
	m=m+1
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)

	ii=io(i)
	b(m,1)=r(ii,1)
	b(m,2)=r(ii,2)
	b(m,3)=r(ii,3)
	part(m)='O'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	inb1(m)=ind1(ii)
	ibo(i)=m
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	xmb(m-2)=xmbead
	m=m+1
!C
!C beads corresponents a la cadena lateral	
!C
	if(res(ii).eq.'ARG')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,4
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=5,7
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S2'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'ASN')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,4
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	xmb(m)=xmbead
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'ASP')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,4
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii)(1:2).eq.'CY')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,2
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'GLN')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,5
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'GLU')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,5
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii)(1:2).eq.'HI')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,2
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=3,4
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S2'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=5,6
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S3'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'ILE')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,4
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'LEU')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,4
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'LYS')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,3
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=4,5
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S2'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'MET')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,4
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'PHE')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,3
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=4,5
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S2'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=6,7
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S3'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'PRO')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,3
	ii=ica(i)-j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'SER')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,2
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'THR')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,3
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'TRP')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,2
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=3,5
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S2'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=8,10
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S3'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=6,7
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S4'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'TYR')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,3
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=4,5
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S2'
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	xmb(m)=xmbead
	inb1(m)=ind1(ii)
	m=m+1
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=6,8
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S3'
	xmb(m)=xmbead
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	inb1(m)=ind1(ii)
	m=m+1
	endif
!C
	if(res(ii).eq.'VAL')then
	rcm(1)=0.
	rcm(2)=0.
	rcm(3)=0.
	xmbead=0.
	do j=1,3
	ii=ica(i)+j
	rcm(1)=rcm(1)+xm(ii)*r(ii,1)
	rcm(2)=rcm(2)+xm(ii)*r(ii,2)
	rcm(3)=rcm(3)+xm(ii)*r(ii,3)
	xmbead=xmbead+xm(ii)
	enddo
	b(m,1)=rcm(1)/xmbead
	b(m,2)=rcm(2)/xmbead
	b(m,3)=rcm(3)/xmbead
	part(m)='S1'
	xmb(m)=xmbead
	bead(m)=res(ii)
	cbead(m)=cad(ii)
	inb1(m)=ind1(ii)
	m=m+1
	endif

	enddo


	nbead=m-1



!C topologia dels enllaços covalents
	do k=1,nres
!C	write(6,*)bead(k)
	i0=ibca(k-1)
	j0=ibn(k-1)
	i1=ibca(k)
	j1=ibn(k)
	i2=ibca(k+1)
	j2=ibn(k+1)
	
	write(6,*)i0,i1,i2," j ",j0,j1,j2
!C bonds backbone
	if(bead(i1).ne.'PRO')then
	ipsb(i1-2,i1-1)=1
	ipsb(i1-2,i1)=1
	ipsb(i1-2,i1+1)=1
	ipsb(i1-1,i1)=1
	else
!C si es prolina, no hi ha H despres del N
	ipsb(i1-1,i1)=1
	ipsb(i1-1,i1+1)=1
!C fixa el diedre phi
	ipsb(i0+1,i1+1)=1
	endif
	ipsb(i1,i1+1)=1
	ipsb(i1,i1+2)=1
	ipsb(i1,i2-1)=1
	ipsb(i1,i2-2)=1
	ipsb(i1,i2)=1
	ipsb(i1+1,i1+2)=1
	ipsb(i1+1,i2-2)=1
	ipsb(i1+1,i2-1)=1
	ipsb(i1+1,i2)=1
	ipsb(i1+2,i2-2)=1
	ipsb(i1+2,i2)=1
!C bonds sidechain
	imax=j2-1
	if(imax.gt.i1+2)then
	ipsb(j1,i1+3)=1
	ipsb(i1+1,i1+3)=1
	endif
	if(k.eq.nres)imax=nbead
	do i=i1+3,imax
	ipsb(i1,i)=1
	do j=i+1,imax
	ipsb(i,j)=1
	enddo
	enddo
	if(part(i1+3).eq.'S1')then
	if(bead(k).ne.'PRO')then
	ipsb(i1-2,i1+3)=1
	else
	ipsb(i1-1,i1+3)=1
	endif
	ipsb(i1+1,i1+3)=1
	endif
	enddo

!C bloc que reconeix els ponts disulfur 
	do i=1,natom-1
	if(atom(i).eq.'SG')then
	do j=i+1,natom
	if(atom(j).eq.'SG')then
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(rij.lt.2.5)then
	ii=ind2(i)
	jj=ind2(j)
	i1=ibca(ii)
	i2=ibca(jj)
	ipsb(i1+3,i2+3)=1
	ipsb(i1+3,i2)=1
	ipsb(i1,i2+3)=1
	endif
	endif
	enddo
	endif
	enddo


	c1='ATOM'



	write(6,*)'natom0,natom',natom0,natom,nbead

	if(nbead.gt.nmax)then
	write(6,*)'Massa particules!',nbead
	write(6,*)'El maxim possible es',nmax
	stop
	endif


	write(6,*)natom0,natom
	open(unit=10,file='topcg.dat')
	do i=1,nbead-1
	do j=i+1,nbead
	rij1=b(i,1)-b(j,1)
	rij2=b(i,2)-b(j,2)
	rij3=b(i,3)-b(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(cbead(i).eq.cbead(j))then
	  if(ipsb(i,j).eq.1)write(10,*)i,j,rij
	endif
	enddo
	enddo
	close(10)

	c1='ATOM'
	l=1
	open(unit=19,file=file19)
	do i=1,nbead
	write(19,2000)c1,i,part(i),bead(i),cbead(i),inb1(i),b(i,1),b(i,2),b(i,3)
	enddo
	close(19)


	stop

	end
