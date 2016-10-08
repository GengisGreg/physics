!!Programma che calcola la cascata diretta mmt, hamiltoniana, incrementi e istogramma pdf degli incrementi!!
program pdfd
implicit none

integer i,j,ier,inc,lenc,lensav,lenwrk,seed
integer, parameter::n=4096,m=100000,nprints=1000
real*8, parameter::pi=2.0d0*dasin(1.0d0),l=2.d0*pi
real*8 lin,nlin
real*8 dt,dk,k(n),h,esp(n),ham,dx,quart(n),f(n),numeno,nupiu,nukm(n),nukp(n),psiquadro(n),deltau(n)
complex*16 u(n),psi(n),derfraz(n),u0(n),w(n)
complex*16,parameter::im=(0.d0,1.d0)
real*8,allocatable,dimension(:):: work
real*8,allocatable,dimension(:):: wsave

!!cose per l'uso della libreria fftpack5!!
lensav=2*n+int(log(real(n))/log(2.0d0))+4
lenwrk=2*n
allocate(work(1:lenwrk))
allocate(wsave(1:lensav))

dt=0.015d0
h=l/dfloat(n)
dk=2.0d0*pi/l

open(unit=31,file="realeudiretta.dat")
open(unit=46,file="cascatadiretta.dat")
open(unit=53,file="ham.dat")
open(unit=63,file="kappa.dat")
open(unit=34,file="rapportoham.dat")


!!SUBROUTINE DI INIZIALIZZAZIONE FFTPACK5!!
 call cfft1i(n,wsave,lensav,ier)
 inc=1
 lenc=n

!!MODI K WRAP AROUND FFT!!
do i=1,n

	if (i.le.n/2+1) then
		k(i)=dk*(i-1)
	else
		k(i)=dk*(i-n-1)
	end if
write(63,*)k(i)
end do

do i=1,n
f(i)=0.d0

	if((i.eq.3).or.(i.eq.4).or.(i.eq.4094).or.(i.eq.4095)) then!(i.eq.16382).or.(i.eq.16383)) then!

		f(i)=0.1d0

	end if

	if((i.eq.2).or.(i.eq.4096)) then!(i.eq.16384)) then!
		nukm(i)=10d0
	end if

	if((i.gt.2031).and.(i.lt.2067)) then
	!if((i.gt.8173).and.(i.lt.8212)) then	
		nukp(i)=10d0
	end if

end do

do i=1,n
u0(i)=f(i)
end do

!!CALCOLO ESPONENZIALE OPERATORE LINEARE!!
do i=1,n
esp(i)=1.0d0*((dabs(k(i)))**(0.5d0))
end do
!!LOOP SUL TEMPO!!

do j=1,m

!!CALCOLO PEZZO OPERATORE NON LINEARE!!
do i=1,n
	u(i)=cdexp(im*dt*(cdabs(u0(i)))**2)*u0(i)
end do
!!FFT DELLA PARTE NON LINEARE!!

 call cfft1f(n,inc,u,lenc,wsave,lensav,work,lenwrk,ier)

!!MOLTIPLICO PER PARTE LINEARE!!
do i=1,n
u(i)=cdexp(im*dt*esp(i))*u(i)
end do
!!DE-ALIASING!
do i=n/4,n-n/4
	u(i)=0.d0	
end do

!!MOLTIPLICO PER FATTORE FORZANTE E DISSIPATIVO!!
do i=1,n
u(i)=dexp((f(i)-nukm(i)-nukp(i))*dt)*u(i)
end do

do i=1,n
psiquadro(i)=psiquadro(i)+dt*((cdabs(u(i)))**2)
end do
!!ANTITRASFORMO CON FFT!!

 call cfft1b(n,inc,u,lenc,wsave,lensav,work,lenwrk,ier)

!!CALCOLO HAMILTONIANA!!
!!#################################################!!
do i=1,n
	psi(i)=u(i)
end do

do i=1,n
	quart(i)=0.5d0*(cdabs(psi(i))**4)
end do

call cfft1f(n,inc,psi,lenc,wsave,lensav,work,lenwrk,ier)

do i=1,n
	psi(i)=(cdabs(im*k(i))**(0.25d0))*psi(i)
end do

call cfft1b(n,inc,psi,lenc,wsave,lensav,work,lenwrk,ier)

do i=1,n
	derfraz(i)=cdabs(psi(i))*cdabs(psi(i))
end do

ham=0.0d0
do i=1,n
	ham=ham+h*(derfraz(i)+quart(i))
	nlin=nlin+quart(i)
	lin=lin+derfraz(i)
end do
write(53,*)j,ham
write(34,*)j,nlin/lin
!!###########################################################!!

!if(mod(j,nprints).eq.0) then
if(j.eq.80000) then
do i=1,n
	write(31,*)dreal(u(i))
end do
	write(31,*) 
psiquadro=psiquadro/dfloat(j)

	do i=1,n
		write(46,*)k(i),psiquadro(i)
	end do
		write(46,*) 
end if

do i=1,n
	u0(i)=u(i)
end do

end do

 close(34)
 close(31)
 close(46)
 close(53)
 close(63)
 
deallocate(wsave)
deallocate(work)
end program pdfd

