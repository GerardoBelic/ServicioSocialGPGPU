


program cycles



a=10
b=1000

do j=1,20

write(6,*)j

if(j.ge.10)then
write(6,*)"j >= 10"
   if(j.eq.10)then
   b=j+b   !!si j==10 b=1010 si no b=1000

	write(6,*)"Neverland j = 10"
	write(6,*)b,a,j
   endif
endif

enddo

write(6,*)b,a

end
