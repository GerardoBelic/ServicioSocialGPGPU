	program topol
	parameter (natmax=10000)
	character*5 c1,c2,c3,atom(natmax),res(natmax),tipus(100),cadind(5)
	dimension imol(natmax)
	dimension ind2(natmax),r(natmax,3),rb(natmax,natmax),iread(2)
 1000	FORMAT (A4,3X,I4,2X,A3,1X,A3,3X,I3,4X,F8.3,F8.3,F8.3)
 2000	FORMAT (A4,3X,I4,2X,A3,1X,A3,1X,A1,1X,I3,4X,F8.3,F8.3,F8.3)
	DATA (tipus(i),i=1,36)/'CA','C','N','O','OXT','CB','CG','CG1','CG2','NG','NG1','NG2','SG','OG','OG1','OG2','CD','CD1','CD2','ND','ND1','ND2','SD','OD','OD1','OD2','CE','CE1','CE2','OE','OE1','OE2','NE','NE1','NE2','NZ'/
	DATA (tipus(i),i=37,46)/'CZ','NH1','NH2','OH','CE3','CZ2','CZ3','CH2','H','H1'/
	DATA (cadind(i),i=1,5)/'A','B','C','D','E'/

	nmoc=1
	natom0=0
	n=1
	im=1
c	read(5,*)
 1	iseg=0 
 50	read(5,1000,end=51)c1,j,c2,c3,k,x,y,z
	if(c1.eq.'END')goto 2
	iflag=0
	do i=1,100
	if(c2.eq.tipus(i))iflag=1
	enddo
	if(iflag.eq.1)then
	atom(n)=c2
	if(c2.eq.'H1')atom(n)='H'
	res(n)=c3
	ind2(n)=k
	r(n,1)=x
	r(n,2)=y
	r(n,3)=z
	imol(n)=im
	n=n+1
	endif
	if(c1.eq.'TER')then
	iseg=1
	im=im+1
c	write(6,*)n
c	stop
	goto 51
	endif
	goto 50
 52  	iseg=0
 51	natom=n-1

ccccccccccccccccccccccccccccccccccccccc
c defineix la matriu de topologia
ccccccccccccccccccccccccccccccccccccccc
	do i=natom0+1,natom-1
	do j=i+1,natom
	rb(i,j)=0.
	enddo
	enddo

	do i=natom0+1,natom
	if(atom(i).eq.'CA')then
	j=i
 150	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'C')rb(i,j)=rij
	if(atom(j).eq.'O')rb(i,j)=rij
	if(atom(j).eq.'OXT')goto 151
	if(atom(j).eq.'N')rb(i,j)=rij
	if(atom(j).eq.'H')rb(i,j)=rij
	if(atom(j).eq.'CB')rb(i,j)=rij
	if(atom(j).eq.'CG')rb(i,j)=rij
	if(atom(j).eq.'CG1')rb(i,j)=rij
	if(atom(j).eq.'CG2')rb(i,j)=rij
	if(atom(j).eq.'OG')rb(i,j)=rij
	if(atom(j).eq.'OG1')rb(i,j)=rij
	if(atom(j).eq.'OG2')rb(i,j)=rij
	if(atom(j).eq.'SG')rb(i,j)=rij
	if(atom(j).eq.'CA')goto 151
	goto 150
 151	rb(i,j)=rij
 155	continue
	elseif(atom(i).eq.'N')then
	j=i
 250	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'CD'.and.res(j).eq.'PRO')rb(i,j)=rij
	if(atom(j).eq.'CG'.and.res(j).eq.'PRO')rb(i,j)=rij
	if(atom(j).eq.'CA')rb(i,j)=rij
	if(atom(j).eq.'CB')rb(i,j)=rij
	if(atom(j).eq.'H')rb(i,j)=rij
	if(atom(j).eq.'C')goto 251
	goto 250
 251	rb(i,j)=rij
	elseif(atom(i).eq.'CB')then
	j=i
 350	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'CG')rb(i,j)=rij
	if(atom(j).eq.'CG1')rb(i,j)=rij
	if(atom(j).eq.'CG2')rb(i,j)=rij
	if(atom(j).eq.'OG')rb(i,j)=rij
	if(atom(j).eq.'OG1')rb(i,j)=rij
	if(atom(j).eq.'OG2')rb(i,j)=rij
	if(atom(j).eq.'CD')rb(i,j)=rij
	if(atom(j).eq.'CD1')rb(i,j)=rij
	if(atom(j).eq.'CD2')rb(i,j)=rij
	if(atom(j).eq.'OD')rb(i,j)=rij
	if(atom(j).eq.'OD1')rb(i,j)=rij
	if(atom(j).eq.'OD2')rb(i,j)=rij
	if(atom(j).eq.'ND')rb(i,j)=rij
	if(atom(j).eq.'ND1')rb(i,j)=rij
	if(atom(j).eq.'ND2')rb(i,j)=rij
	if(atom(j).eq.'SG')rb(i,j)=rij
	if(atom(j).eq.'SD')rb(i,j)=rij
	if(atom(j).eq.'OH')rb(i,j)=rij
	if(res(j)(1:2).eq.'HI'.and.(atom(j).eq.'NE2'.or.atom(j).eq.'CE1'))rb(i,j)=rij
	if(atom(j).eq.'CE2')rb(i,j)=rij
	if(atom(j).eq.'CE3')rb(i,j)=rij
	if(atom(j)(1:2).eq.'CZ'.and.res(j).ne.'ARG')rb(i,j)=rij
	if(atom(j).eq.'CH2')rb(i,j)=rij
c cas de les prolines
	if(atom(j).eq.'CA')rb(i,j)=rij
	if(atom(j).eq.'C')goto 351
	goto 350
 351	rb(i,j)=rij
	elseif(atom(i).eq.'CG')then
	j=i
 450	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'CA')rb(i,j)=rij
	if(atom(j).eq.'CB')rb(i,j)=rij
	if(atom(j).eq.'CD')rb(i,j)=rij
	if(atom(j).eq.'CD1')rb(i,j)=rij
	if(atom(j).eq.'CD2')rb(i,j)=rij
	if(atom(j).eq.'OD')rb(i,j)=rij
	if(atom(j).eq.'OD1')rb(i,j)=rij
	if(atom(j).eq.'OD2')rb(i,j)=rij
	if(atom(j).eq.'ND')rb(i,j)=rij
	if(atom(j).eq.'ND1')rb(i,j)=rij
	if(atom(j).eq.'ND2')rb(i,j)=rij
	if(atom(j).eq.'CE')rb(i,j)=rij
	if(atom(j).eq.'CE1')rb(i,j)=rij
	if(atom(j).eq.'CE2')rb(i,j)=rij
	if(atom(j).eq.'OE')rb(i,j)=rij
	if(atom(j).eq.'OE1')rb(i,j)=rij
	if(atom(j).eq.'OE2')rb(i,j)=rij
	if(atom(j).eq.'NE')rb(i,j)=rij
	if(atom(j).eq.'NE1')rb(i,j)=rij
	if(atom(j).eq.'NE2')rb(i,j)=rij
	if(atom(j).eq.'SD')rb(i,j)=rij
	if(atom(j).eq.'CE2')rb(i,j)=rij
	if(atom(j).eq.'CE3')rb(i,j)=rij
	if(atom(j)(1:2).eq.'CZ'.and.res(j).ne.'ARG')rb(i,j)=rij
	if(atom(j).eq.'CH2')rb(i,j)=rij
	if(atom(j).eq.'C')goto 451
	goto 450
 451	continue
	elseif(atom(i).eq.'CG1'.or.atom(i).eq.'OG1')then
	j=i
 550	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'OG2')rb(i,j)=rij
	if(atom(j).eq.'CG2')rb(i,j)=rij
	if(atom(j).eq.'NG2')rb(i,j)=rij
	if(atom(j).eq.'CD1')rb(i,j)=rij
	if(atom(j).eq.'OD1')rb(i,j)=rij
	if(atom(j).eq.'ND1')rb(i,j)=rij
	if(atom(j).eq.'CE1')rb(i,j)=rij
	if(atom(j).eq.'OE1')rb(i,j)=rij
	if(atom(j).eq.'NE1')rb(i,j)=rij
	if(atom(j).eq.'C')goto 551
	goto 550
 551	continue
	elseif(atom(i).eq.'CG2'.or.atom(i).eq.'OG2')then
	j=i
 650	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'OG1')rb(i,j)=rij
	if(atom(j).eq.'CG1')rb(i,j)=rij
	if(atom(j).eq.'NG1')rb(i,j)=rij
	if(atom(j).eq.'CD2')rb(i,j)=rij
	if(atom(j).eq.'OD2')rb(i,j)=rij
	if(atom(j).eq.'ND2')rb(i,j)=rij
	if(atom(j).eq.'CE2')rb(i,j)=rij
	if(atom(j).eq.'OE2')rb(i,j)=rij
	if(atom(j).eq.'NE2')rb(i,j)=rij
	if(atom(j).eq.'C')goto 651
	goto 650
 651	continue
	elseif(atom(i).eq.'CD'.or.atom(i).eq.'SD')then
	j=i
 750	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'CE')rb(i,j)=rij
	if(atom(j).eq.'CE1')rb(i,j)=rij
	if(atom(j).eq.'OE1')rb(i,j)=rij
	if(atom(j).eq.'NE1')rb(i,j)=rij
	if(atom(j).eq.'CE2')rb(i,j)=rij
	if(atom(j).eq.'OE2')rb(i,j)=rij
	if(atom(j).eq.'NE2')rb(i,j)=rij
	if(atom(j).eq.'NE')rb(i,j)=rij
	if(atom(j).eq.'NZ')rb(i,j)=rij
	if(atom(j).eq.'CZ')rb(i,j)=rij
c cas de les prolines
	if(atom(j).eq.'CG')rb(i,j)=rij
	if(atom(j).eq.'CB')rb(i,j)=rij
	if(atom(j).eq.'CA')rb(i,j)=rij
	if(atom(j).eq.'C')goto 751
	goto 750
 751	continue
	elseif(atom(i).eq.'CD1'.or.atom(i).eq.'OD1'.or.atom(i).eq.'ND1'.or.atom(i).eq.'CD2'.or.atom(i).eq.'OD2'.or.atom(i).eq.'ND2')then
	j=i
 850	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'CD1'.or.atom(j).eq.'OD1'.or.atom(j).eq.'ND1')rb(i,j)=rij
	if(atom(j).eq.'CD2'.or.atom(j).eq.'OD2'.or.atom(j).eq.'ND2')rb(i,j)=rij
	if(atom(j).eq.'CE1'.or.atom(j).eq.'OE1'.or.atom(j).eq.'NE1')rb(i,j)=rij
	if(atom(j).eq.'CE2'.or.atom(j).eq.'OE2'.or.atom(j).eq.'NE2')rb(i,j)=rij
	if(atom(j)(1:2).eq.'CZ')rb(i,j)=rij
	if(atom(j).eq.'CE3')rb(i,j)=rij
	if(atom(j).eq.'CH2')rb(i,j)=rij
	if(atom(j).eq.'C')goto 851
	goto 850
 851	continue
	elseif(atom(i).eq.'NH1'.or.atom(i).eq.'NH2')then
	j=i
 950	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'NH1')rb(i,j)=rij
	if(atom(j).eq.'NH2')rb(i,j)=rij
	if(atom(j).eq.'C')goto 951
	goto 950
 951	continue
	elseif(atom(i).eq.'CE'.or.atom(i).eq.'CE1'.or.atom(i).eq.'CE2'.or.atom(i).eq.'OE1'.or.atom(i).eq.'OE2'.or.atom(i).eq.'NE1'.or.atom(i).eq.'NE2')then
	j=i
 1050	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'NZ')rb(i,j)=rij
	if(atom(j).eq.'CE1')rb(i,j)=rij
	if(atom(j).eq.'CE2')rb(i,j)=rij
	if(atom(j).eq.'OE1')rb(i,j)=rij
	if(atom(j).eq.'OE2')rb(i,j)=rij
	if(atom(j).eq.'NE1')rb(i,j)=rij
	if(atom(j).eq.'NE2')rb(i,j)=rij
	if(atom(j).eq.'CD1')rb(i,j)=rij
	if(atom(j).eq.'CD2')rb(i,j)=rij
	if(atom(j).eq.'OD1')rb(i,j)=rij
	if(atom(j).eq.'OD2')rb(i,j)=rij
	if(atom(j).eq.'ND1')rb(i,j)=rij
	if(atom(j).eq.'ND2')rb(i,j)=rij
	if(atom(j)(1:2).eq.'CZ')rb(i,j)=rij
	if(atom(j).eq.'CH2')rb(i,j)=rij
	if(atom(j).eq.'CE3')rb(i,j)=rij
	if(atom(j).eq.'C')goto 1051
	goto 1050
 1051	continue
	elseif(atom(i).eq.'NE'.or.atom(i).eq.'CZ')then
	j=i
 1150	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'NH1')rb(i,j)=rij
	if(atom(j).eq.'NH2')rb(i,j)=rij
	if(atom(j).eq.'OH')rb(i,j)=rij
	if(atom(j).eq.'CE2')rb(i,j)=rij
	if(atom(j).eq.'CD2')rb(i,j)=rij
	if(atom(j).eq.'CZ')rb(i,j)=rij
	if(atom(j).eq.'C')goto 1151
	goto 1150
 1151	continue
c topologia del TRP
	elseif(atom(i).eq.'CZ2'.or.atom(i).eq.'CZ3'.or.atom(i).eq.'CH2'.or.atom(i).eq.'CE3')then
	j=i
 1250	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'CZ2'.or.atom(j).eq.'CZ3'.or.atom(j).eq.'CH2'.or.atom(j).eq.'CE3'.or.atom(j).eq.'CD2')rb(i,j)=rij	
	if(atom(j).eq.'C')goto 1251
	goto 1250
 1251	continue
	elseif(atom(i).eq.'O')then
	j=i
 1350	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'N')rb(i,j)=rij
	if(atom(j).eq.'CA')goto 1351
	if(atom(j).eq.'OXT')goto 1351
	goto 1350
 1351	rb(i,j)=rij
 1352	continue
	elseif(atom(i).eq.'C')then
	j=i
 1450	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'O')rb(i,j)=rij
	if(atom(j).eq.'OXT')goto 1451
	if(atom(j).eq.'N')rb(i,j)=rij
	if(atom(j).eq.'H')rb(i,j)=rij
	if(atom(j).eq.'CA')goto 1451
	goto 1450
 1451	rb(i,j)=rij
 1452	continue
	elseif(atom(i).eq.'H')then
	j=i
 1550	j=j+1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	if(atom(j).eq.'OXT')goto 1552
	if(atom(j).eq.'CA')goto 1551
	goto 1550
 1551	rb(i,j)=rij
	endif
	enddo
c s'ha acabat la proteina
 1552	rb(i,j)=rij
	natom=j
c connexio O-OXT
	i=j-1
	rij1=r(i,1)-r(j,1)
	rij2=r(i,2)-r(j,2)
	rij3=r(i,3)-r(j,3)
	rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
	rb(i,j)=rij
c	write(6,*)natom0,natom
	natom0=natom
	if(iseg.eq.1)then
	nmoc=nmoc+1
	goto 1
	endif
c	natom=natom-1
 2	c1='ATOM'
	open(unit=9,file='structure.pdb')
	do i=1,natom
	k=ind2(i)-ind2(1)+1
	kk=imol(i)
	write(9,2000)c1,i,atom(i),res(i),cadind(kk),k,(r(i,j),j=1,3)
	enddo
	close(9)
c	write(6,*)nmoc,' molecules'
	open(unit=9,file='topology.dat')
	do i=1,natom-1
	do j=i+1,natom
	if(imol(i).eq.imol(j))then
	if(rb(i,j).gt.1.d-10)write(9,*)i,j,rb(i,j)
	endif
	enddo
	enddo
	close(9)
	stop
	end
