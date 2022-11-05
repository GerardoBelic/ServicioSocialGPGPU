	


	character*5 c1
	dimension i(10000),j(10000),rij(10000)
	read(5,*)c1,nmol
	write(6,*)nmol
	l=1
 100  	read(5,*,end=101)i(l),j(l),rij(l)
c	write(6,*)i(l),j(l),rij(l)
	l=l+1
	goto 100
 101	ncov=l-1
  	natom=j(ncov)
!c	stop
!c	write(6,*)l,natom
	do k=0,nmol
	do l=1,ncov
	write(6,*)i(l)+k*natom,j(l)+k*natom,rij(l)
	enddo
	enddo
	stop
	end
