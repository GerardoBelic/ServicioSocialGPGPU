
using DelimitedFiles
using DataStructures

using BenchmarkTools

input0=readdlm(stdin,String)



function caso(atom::String)
	
	m=UInt16(0)	
	
	if 	atom=="CA"
		m=UInt16(1)
	elseif  atom=="N" 
		m=UInt16(2)
	elseif  atom=="CB"
		m=UInt16(3)
	elseif  atom=="CG"
		m=UInt16(4)
	elseif  atom=="CG1"||atom=="OG1"
		m=UInt16(5)
	elseif  atom=="CG2"||atom=="OG2"
		m=UInt16(6)
	elseif  atom=="CD"||atom=="SD"
		m=UInt16(7)
	elseif  atom=="CD1"||atom=="OD1"||atom=="ND1"||atom=="CD2"||atom=="OD2"||atom=="ND2"
		m=UInt16(8)
	elseif  atom=="NH1"||atom=="NH2"
		m=UInt16(9)
	elseif  atom=="CE"||atom=="CE1"||atom=="CE2"||atom=="OE1"||atom=="OE2"||atom=="NE1"||atom=="NE2"
		m=UInt16(10)
	elseif  atom=="NE"||atom=="CZ"
		m=UInt16(11)
	elseif  atom=="CZ2"||atom=="CZ3"||atom=="CH2"||atom=="CE3"
		m=UInt16(12)
	elseif  atom=="O"
		m=UInt16(13)
	elseif  atom=="C"
		m=UInt16(14)
	elseif  atom=="H"	
		m=UInt16(15)
	end

	
	return m

end

function dtipus(a::String) 

	f =Bool(false)
	 tipus =["CA","C","N","O","OXT","CB","CG","CG1","CG2",
              "NG","NG1","NG2","SG","OG","OG1","OG2","CD",
              "CD1","CD2","ND","ND1","ND2","SD","OD","OD1",
              "OD2","CE","CE1","CE2","OE","OE1","OE2","NE",
              "NE1","NE2","NZ","CZ","NH1","NH2","OH","CE3",
              "CZ2","CZ3","CH2","H","H1"]
       for i=1:46
           if tipus[i] == a
               f=true
           end
       end
	return f
end


function cadind(i::UInt16)

	f = String("")

	 c =["A","B","C","D","E"]

	if i<5 && i>0
               f=c[i]
           end
	
	return f

end



function topol(input::Array{String})

      natmax=UInt16(6000)              ##Número máximo de elementos
      atom=Array{String}(undef,natmax)  ##Almacena átomos o moléculas
      res=Array{String}(undef,natmax)   ##Almacena alias de los residuos
      

      imol=Array{UInt16}(zeros(natmax))   ##Índice de la molécula
      ind2=Array{UInt16}(zeros(natmax))   ##Índice del residuo
      r=Array{Float64}(zeros(natmax,3))   ##Magnitudes
      rb=Array{Float64}(zeros(natmax,natmax))  ##Distancias en el espacio

    

	rij=Int64(0)


	


	##nmoc=1   inecesario hace lo mismo que ims
	natom0=Int64(0)   ##Índice de referencia para leer la proteína
	natom=Int64(0)
	n=Int64(1)
	i=Int64(0)
	j=Int64(0)
	ims=Int64(1) #=im le cambio su nombre a ims ya que im
		se le suele utilizar para tipos con números complejos
		
		+++ims cuenta moléculas++++++	
	      =# 

	iseg=Int64(0)   ##Bandera para indicar que se encontró TER en la STDIN

##Se realiza la lectura de la Entrada Estándar
	
	

	limite=natmax ##Agregado para dar robuztez
	c1=String("0")
	c2=String("0")
	c3=String("0")
	k=UInt16(0)
	x=Float64(0)
	y=Float64(0)
	z=Float64(0)

	nn=1
###BLOQUE 1++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	while c1 != "END" && nn<=limite
		c1=input[nn,1]
		c2=String("")

		if c1!="TER"&&c1!="END"
			j=parse(UInt16,input[nn,2])
			c2=input[nn,3]
			c3=input[nn,4]
			k=parse(UInt16,input[nn,5])
			x=parse(Float64,input[nn,6])
			y=parse(Float64,input[nn,7])
			z=parse(Float64,input[nn,8])
		end
	##Segunda y Tercera SC, se pueden juntar en un sólo código
	##eliminando a iflag
		if dtipus(c2)
			 c2=="H1" ? atom[n]="H" : atom[n]=c2	
			 res[n]=c3
			 ind2[n]=k
			 r[n,1]=x
			 r[n,2]=y
			 r[n,3]=z
			 imol[n]=ims
			 n=n+1
		end
	

		if c1=="TER"
			iseg=1
			ims=ims+1
		end
##FIN BLOQUE 1++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		
		##Bloque 2##############################################
		if iseg==1
			natom=n-1

			##Primer SC se omite ya que solo sirve
			##para inicializan los valores de rb
			##para una molécula	
			i=natom0+1	
			while i<natom
			
				m=caso(atom[i])



				if m==1

					j=i
					rij=0.
					while   j<limite
						j+=1	
						rij1=r[i,1]-r[j,1]
						rij2=r[i,2]-r[j,2]
						rij3=r[i,3]-r[j,3]	
						rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)

						if atom[j] =="C"  rb[i,j]=rij
						elseif atom[j] =="O"  rb[i,j]=rij
						elseif atom[j] =="N"  rb[i,j]=rij
						elseif atom[j] =="H"  rb[i,j]=rij
						elseif atom[j] =="CB"  rb[i,j]=rij
						elseif atom[j] =="CG"  rb[i,j]=rij
						elseif atom[j] =="CG1"  rb[i,j]=rij
						elseif atom[j] =="CG2"  rb[i,j]=rij
						elseif atom[j] =="OG"  rb[i,j]=rij
						elseif atom[j] =="OG1"  rb[i,j]=rij
						elseif atom[j] =="OG2"  rb[i,j]=rij
						elseif atom[j] =="SG"  rb[i,j]=rij
						elseif atom[j] =="CA" || atom[j] =="OXT"
						break
						end
					end
						rb[i,j]=rij

	
				elseif m==2
					j=i
					rij=0.

					while   j<limite
						j+=1	
						rij1=r[i,1]-r[j,1]
						rij2=r[i,2]-r[j,2]
						rij3=r[i,3]-r[j,3]	
						rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)

						if atom[j]=="C"
							break
						elseif atom[j]=="H" 
							rb[i,j]=rij 
						elseif atom[j]=="CD"&&res[j]=="PRO" rb[i,j]=rij
						elseif atom[j]=="CG"&&res[j]=="PRO" rb[i,j]=rij
						elseif atom[j]=="CA" rb[i,j]=rij
						elseif atom[j]=="CB" rb[i,j]=rij
						end
					end

					rb[i,j]=rij



				elseif m==3


						j=i
						rij=0.
						while   j<limite
							
							j+=1	
							rij1=r[i,1]-r[j,1]
							rij2=r[i,2]-r[j,2]
							rij3=r[i,3]-r[j,3]	
							rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)	

							if atom[j]=="C"
							break
							elseif atom[j]=="CG" rb[i,j]=rij
							elseif atom[j]=="CG1" rb[i,j]=rij
							elseif atom[j]=="CG2" rb[i,j]=rij
							elseif atom[j]=="OG" rb[i,j]=rij
							elseif atom[j]=="OG1" rb[i,j]=rij
							elseif atom[j]=="OG2" rb[i,j]=rij
							elseif atom[j]=="CD" rb[i,j]=rij
							elseif atom[j]=="CD1" rb[i,j]=rij
							elseif atom[j]=="CD2" rb[i,j]=rij
							elseif atom[j]=="OD" rb[i,j]=rij
							elseif atom[j]=="OD1" rb[i,j]=rij
							elseif atom[j]=="OD2" rb[i,j]=rij
							elseif atom[j]=="ND" rb[i,j]=rij
							elseif atom[j]=="ND1" rb[i,j]=rij
							elseif atom[j]=="ND2" rb[i,j]=rij
							elseif atom[j]=="SG" rb[i,j]=rij
							elseif atom[j]=="SD" rb[i,j]=rij
							elseif atom[j]=="OH" rb[i,j]=rij
							elseif atom[j]=="CE2" rb[i,j]=rij
							elseif atom[j]=="CE3" rb[i,j]=rij
							elseif atom[j]=="CH2" rb[i,j]=rij
							elseif atom[j]=="CA" rb[i,j]=rij
							elseif res[j][1:2]=="HI" &&(atom[j]=="NE2"||atom[j]=="CE1" ) rb[i,j]=rij 
							elseif atom[j][1:2]=="CZ" && res[j]!="ARG" rb[i,j]=rij
							end 
						end
						rb[i,j]=rij

				elseif m==4
						j=i
							while   j<limite
								j+=1	
								rij1=r[i,1]-r[j,1]
								rij2=r[i,2]-r[j,2]
								rij3=r[i,3]-r[j,3]	
								rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)	

								if atom[j]=="C"
									break      
								elseif atom[j]=="CA" rb[i,j]=rij
								elseif atom[j]=="CB" rb[i,j]=rij
								elseif atom[j]=="CD" rb[i,j]=rij
								elseif atom[j]=="CD1" rb[i,j]=rij
								elseif atom[j]=="CD2" rb[i,j]=rij
								elseif atom[j]=="OD" rb[i,j]=rij
								elseif atom[j]=="OD1" rb[i,j]=rij
								elseif atom[j]=="OD2" rb[i,j]=rij
								elseif atom[j]=="ND" rb[i,j]=rij
								elseif atom[j]=="ND1" rb[i,j]=rij
								elseif atom[j]=="ND2" rb[i,j]=rij
								elseif atom[j]=="CE" rb[i,j]=rij
								elseif atom[j]=="CE1" rb[i,j]=rij
								elseif atom[j]=="CE2" rb[i,j]=rij
								elseif atom[j]=="OE" rb[i,j]=rij
								elseif atom[j]=="OE1" rb[i,j]=rij
								elseif atom[j]=="OE2" rb[i,j]=rij
								elseif atom[j]=="NE" rb[i,j]=rij
								elseif atom[j]=="NE1" rb[i,j]=rij
								elseif atom[j]=="NE2" rb[i,j]=rij
								elseif atom[j]=="SD" rb[i,j]=rij
								elseif atom[j]=="CE2" rb[i,j]=rij
								elseif atom[j]=="CE3" rb[i,j]=rij
								elseif atom[j]=="CH2" rb[i,j]=rij
								elseif atom[j][1:2]=="CZ"&&res[j]!="ARG" rb[i,j]=rij
								end
							end


				elseif m==5

						j=i
							while   j<limite
								j+=1	
								rij1=r[i,1]-r[j,1]
								rij2=r[i,2]-r[j,2]
								rij3=r[i,3]-r[j,3]	
								rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)

								if atom[j]=="C"
									break
								elseif atom[j]=="OG2" rb[i,j]=rij
								elseif atom[j]=="CG2" rb[i,j]=rij
								elseif atom[j]=="NG2" rb[i,j]=rij
								elseif atom[j]=="CD1" rb[i,j]=rij
								elseif atom[j]=="OD1" rb[i,j]=rij
								elseif atom[j]=="ND1" rb[i,j]=rij
								elseif atom[j]=="CE1" rb[i,j]=rij
								elseif atom[j]=="OE1" rb[i,j]=rij
								elseif atom[j]=="NE1" rb[i,j]=rij
								
								end
							end

				elseif m==6
						j=i
							while   j<limite
								j+=1	
								rij1=r[i,1]-r[j,1]
								rij2=r[i,2]-r[j,2]
								rij3=r[i,3]-r[j,3]	
								rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)

								if atom[j]=="C"
								break
								elseif atom[j]=="OG1" rb[i,j]=rij
								elseif atom[j]=="CG1" rb[i,j]=rij
								elseif atom[j]=="NG1" rb[i,j]=rij
								elseif atom[j]=="CD2" rb[i,j]=rij
								elseif atom[j]=="OD2" rb[i,j]=rij
								elseif atom[j]=="ND2" rb[i,j]=rij
								elseif atom[j]=="CE2" rb[i,j]=rij
								elseif atom[j]=="OE2" rb[i,j]=rij
								elseif atom[j]=="NE2" rb[i,j]=rij
								end
							end


				elseif m==7
						j=i
							while   j<limite
								j+=1	
								rij1=r[i,1]-r[j,1]
								rij2=r[i,2]-r[j,2]
								rij3=r[i,3]-r[j,3]	
								rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)

								if atom[j]=="C"
								break
								elseif atom[j]=="CE"  rb[i,j]=rij
								elseif atom[j]=="CE1" rb[i,j]=rij
								elseif atom[j]=="OE1" rb[i,j]=rij
								elseif atom[j]=="NE1" rb[i,j]=rij
								elseif atom[j]=="CE2" rb[i,j]=rij
								elseif atom[j]=="OE2" rb[i,j]=rij
								elseif atom[j]=="NE2" rb[i,j]=rij
								elseif atom[j]=="NE" rb[i,j]=rij
								elseif atom[j]=="NZ" rb[i,j]=rij
								elseif atom[j]=="CZ" rb[i,j]=rij
								elseif atom[j]=="CG" rb[i,j]=rij
								elseif atom[j]=="CB" rb[i,j]=rij
								elseif atom[j]=="CA" rb[i,j]=rij
								end
							end


				elseif m==8

						j=i
							while   j<limite
								j+=1	
								rij1=r[i,1]-r[j,1]
								rij2=r[i,2]-r[j,2]
								rij3=r[i,3]-r[j,3]	
								rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)

								if atom[j]=="C"
								break
								elseif     atom[j]=="CD1"||atom[j]=="OD1"||atom[j]=="ND1" rb[i,j]=rij
								elseif atom[j]=="CD2"||atom[j]=="OD2"||atom[j]=="ND2" rb[i,j]=rij
								elseif atom[j]=="CE1"||atom[j]=="OE1"||atom[j]=="NE1" rb[i,j]=rij
								elseif atom[j]=="CE2"||atom[j]=="OE2"||atom[j]=="NE2" rb[i,j]=rij
								elseif atom[j][1:2]=="CZ" rb[i,j]=rij
								elseif atom[j]=="CE3" rb[i,j]=rij
								elseif atom[j]=="CH2" rb[i,j]=rij
								end
							end

				elseif m==9

						j=i
							while   j<limite
								j+=1	
								rij1=r[i,1]-r[j,1]
								rij2=r[i,2]-r[j,2]
								rij3=r[i,3]-r[j,3]	
								rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
								
								if atom[j]=="C"
								break
								elseif atom[j]=="NH1" rb[i,j]=rij
								elseif atom[j]=="NH2" rb[i,j]=rij
								end
							end

				elseif m==10

					j=i
						while   j<limite
							j+=1	
							rij1=r[i,1]-r[j,1]
							rij2=r[i,2]-r[j,2]
							rij3=r[i,3]-r[j,3]	
							rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)

							if atom[j]=="C"
							break
							elseif atom[j]=="NZ" rb[i,j]=rij
							elseif atom[j]=="CE1" rb[i,j]=rij
							elseif atom[j]=="CE2" rb[i,j]=rij
							elseif atom[j]=="OE1" rb[i,j]=rij
							elseif atom[j]=="OE2" rb[i,j]=rij
							elseif atom[j]=="NE1" rb[i,j]=rij
							elseif atom[j]=="NE2" rb[i,j]=rij
							elseif atom[j]=="CD1" rb[i,j]=rij
							elseif atom[j]=="CD2" rb[i,j]=rij
							elseif atom[j]=="OD1" rb[i,j]=rij
							elseif atom[j]=="OD2" rb[i,j]=rij
							elseif atom[j]=="ND1" rb[i,j]=rij
							elseif atom[j]=="ND2" rb[i,j]=rij
							elseif atom[j][1:2]=="CZ" rb[i,j]=rij
							elseif atom[j]=="CH2" rb[i,j]=rij
							elseif atom[j]=="CE3" rb[i,j]=rij
							end
						end


				elseif m==11
					j=i
						while   j<limite
							j+=1	
							rij1=r[i,1]-r[j,1]
							rij2=r[i,2]-r[j,2]
							rij3=r[i,3]-r[j,3]	
							rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)

							if atom[j]=="C"
							break
							elseif atom[j]=="NH1" rb[i,j]=rij
							elseif atom[j]=="NH2" rb[i,j]=rij
							elseif atom[j]=="OH" rb[i,j]=rij
							elseif atom[j]=="CE2" rb[i,j]=rij
							elseif atom[j]=="CD2" rb[i,j]=rij
							elseif atom[j]=="CZ" rb[i,j]=rij
							end
						end

				elseif m==12
					j=i
						while   j<limite
							j+=1	
							rij1=r[i,1]-r[j,1]
							rij2=r[i,2]-r[j,2]
							rij3=r[i,3]-r[j,3]	
							rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)

							if atom[j]=="C"
							break
							elseif atom[j]=="CZ2"||atom[j]=="CZ3"||atom[j]=="CH2"||atom[j]=="CE3"||atom[j]=="CD2" 
							rb[i,j]=rij	
							end
						end

				elseif m==13
						j=i
						rij=0.
						while   j<limite
							j+=1	
							rij1=r[i,1]-r[j,1]
							rij2=r[i,2]-r[j,2]
							rij3=r[i,3]-r[j,3]	
							rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)

							if atom[j]=="N" rb[i,j]=rij
							elseif atom[j]=="CA" || atom[j]=="OXT"
							break
							end
						end
							rb[i,j]=rij

				elseif m==14

						j=i
						rij=0.
						while   j<limite
							j+=1	
							rij1=r[i,1]-r[j,1]
							rij2=r[i,2]-r[j,2]
							rij3=r[i,3]-r[j,3]	
							rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)

							if atom[j]=="O" rb[i,j]=rij
							elseif atom[j]=="N" rb[i,j]=rij
							elseif atom[j]=="H" rb[i,j]=rij
							elseif atom[j]=="CA"||atom[j]=="OXT" 
							break
							end
						end
							rb[i,j]=rij
			
				
				elseif m==15
				j=i
				band=false ##Para salir del for
					while   j<limite
						j+=1	
						rij1=r[i,1]-r[j,1]
						rij2=r[i,2]-r[j,2]
						rij3=r[i,3]-r[j,3]	
						rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
						if atom[j]=="CA" rb[i,j]=rij
						break
						elseif atom[j]=="OXT"
						band=true
						break
						end

					end  ##FIN while

					if band
					break
					end

				elseif m==0
					##print("\nTERMINO IGNORADO\n")
				        ##print("\nERROR INESPERADO\n")
				end ##FIN IF

				i+=1
			end  ##FIN WHILE


		##Se ha acabado la proteína 
			rb[i,j]=rij
			natom=j

		##Hay una conexión O-OXT o un término X
			i=j-1
			rij1=r[i,1]-r[j,1]
			rij2=r[i,2]-r[j,2]
			rij3=r[i,3]-r[j,3]	
			rij=sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
			natom0=natom		
			##nmoc=nmoc+1 inecesario
		iseg=0
		end	     ##FIN IF
		##FIN Bloque 2###############################################
	#print("$c1 $n $(i::Int64nput[nn+1,1:8]) \n")
	nn+=1
	end    ###Fin del while


	if c1!="END"
		print("\nERROR EL ARCHIVO NO ES FORMATO AMBER\n")
	else
	
	##Impresión del archivo structure.pdb

		if n>1

		#Escritura de structure.pdb
			n=n-1
			atm=fill("ATOM",1,n)		
		
			io=open("structureJulia.pdb","w")
		
			for i=1:n
			 	write(io,"$(atm[i])   $(lpad(i,4))  $(rpad(atom[i],3)) $(lpad(res[i],3)) $(cadind(imol[i]))$(lpad(ind2[i],4))     $(lpad(string(r[i,1]),7)) $(lpad(r[i,2],7)) $(lpad(r[i,3],7))\n")
			end
			
			close(io)
		#
		#Escritura de topology.dat



		 io=open("topologyJulia.dat","w")

			for i=1:n-1
				for j=i+1:n
					if imol[i] == imol[j]
						if rb[i,j] > (10^-10)
						   write(io,"$(lpad(i,12)) $(lpad(j,11))   $(SubString(lpad(rb[i,j],10),1,10))    \n")	
						end	
					end
				end
			end

		 close(io)



		else
			#print("\nError archivo de entrada vacío\n")

		end    #Fin IF

	end            #FIN ELSE 
	##IMPRIMIR ARCHIVOS CON RESULTADOS

end 


topol(["END" ""])
topol(["END" ""])
topol(["END" ""])

@time topol(input0)


