

!!cycle omite todo lo posterior en un ciclo do

!!


program cyclesa



a=10
b=1000

do j=1,20

write(6,*)j

if(j.lt.10)cycle
write(6,*)"j >= 10"
   if(j.gt.10)cycle
   b=j+b   !!si j==10 b=1010 si no b=1000

write(6,*)"Neverland j = 10"
write(6,*)b,a,j

enddo

write(6,*)b,a

end


