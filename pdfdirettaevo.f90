program pdfevo1_1

implicit none 

integer i,j,k,l,e,d
integer,parameter::n=4096,m=100000,ntot=409600,np=1000
real*8 deltau(100,n),freq(ntot),u(ntot),x(100,n),step(n),sqm(100),mean(100),minimo,massimo
real*8 y,gauss,dx,dn
real*8,parameter::pi=2.0d0*dasin(1.0d0)

open(unit=25,file="realeudiretta.dat",status='old')
open(unit=21,file="pdfdirettaevo.dat")
open(unit=34,file="deltau.dat")
open(unit=41,file="gaussiana.dat")
!!leggo!!

do i=1,ntot
	read(25,*)u(i)
end do
!suddivido in gruppi!!
k=1
do j=1,100
	do i=1,n
		x(j,i)=u(k)
	k=k+1
	end do
end do

!k=1
!do j=1,5
!	do i=1,n
!		x(j,i)=u(k)
!		k=k+1
!	end do
!end do

!!calcoli su ogni gruppo!!

write(*,*)'Inserisci la lunghezza di scala'
read(*,*) d

do j=1,100
	do i=1,n
		if((i+d).le.n) then
			deltau(j,i)=x(j,i+d)-x(j,i)
		else if((i+d).gt.n) then
			deltau(j,i)=x(j,i+d-n)-x(j,i)
		end if
	end do
end do

do j=1,100
	mean(j)=0.d0
end do

do j=1,100
	do i=1,n
		mean(j)=mean(j)+deltau(j,i)/dfloat(n)
	end do
end do

write(*,*)'media = ',mean

do j=1,100
	sqm=0.d0
end do

do j=1,100
	do i=1,n			
		sqm(j)=sqm(j)+(deltau(j,i)-mean(j))**2/dfloat(n)
	end do
end do
!write(*,*)'varianza= ',sqm

do j=1,100
	do i=1,n
		deltau(j,i)=(deltau(j,i)-mean(j))/dsqrt(sqm(j))
	end do
end do

minimo=0.d0
massimo=0.d0
massimo=deltau(1,1)
minimo=deltau(1,1)
do j=1,100
	do i=1,n
		if(deltau(j,i).gt.massimo) massimo=deltau(j,i)
	end do

	do i=1,n
		if(deltau(j,i).lt.minimo) minimo=deltau(j,i)
	end do
end do
write(*,*)minimo,massimo

dx=0.4d0
e=(massimo-minimo)/dx

do i=1,e+1
	freq(i)=0.d0
end do

do j=1,e+2
step(j)=minimo+0.2d0+dfloat(j-1)*dx
end do

do j=1,100
      do i=1,n
         if(deltau(j,i).le.step(1)) then
           freq(1)=freq(1)+1.d0
	   else
           if(deltau(j,i).gt.step(e+2)) then
           freq(e+2)=freq(e+2)+1.d0
	   else
         if(deltau(j,i).gt.step(1).and.deltau(j,i).le.step(e+2)) then
           l=int((deltau(j,i)-minimo)/dx)
	   l=l+1
	   freq(l)=freq(l)+1.d0
	  end if
	  end if
	  end if
      end do
end do

do j=1,e+1
	freq(j)=freq(j)/(dfloat(ntot)*dx)
	!freq(i)=freq(i)/(dfloat(n)*dx)
end do

dn=10.d0/dfloat(np)
do j=1,np
	y=-5.d0+dn*(dfloat(j))
	gauss=dexp(-y**2/2.d0)/(dsqrt(2.d0*pi))
	write(41,*) y,gauss
 end do

!do i=1,e+1      
!   write(21,*) step(i),freq(i)
!end do
!	write(21,*) 

!end do

do j=1,e+1      
   write(21,*) step(j),freq(j)!/dfloat(ntot)
end do

 close(41)
 close(25)
 close(21)
 close(34)
end program pdfevo1_1
