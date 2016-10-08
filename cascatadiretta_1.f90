!!!!!ALGORITMO RISOLUZIONE EQUAZIONE MMT CON alpha=1/2 e beta=0, UTILIZZANDO SPLIT-STEP METHOD!!!!!
!!UTILIZZO LIBRERIE FFTPACK5 gfortran nomefile -lfftpack5.1d!!

program mmt
implicit none

integer i,j,ier,inc,lenc,lensav,lenwrk,seed,q
integer, parameter::n=4096,m=100000,nprints=10
real*8, parameter::pi=2.0d0*dasin(1.0d0),l=2.0d0*pi
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
h=l/real(n)
dk=2.0d0*pi/l

open(unit=46,file="cascatadiretta.dat")
open(unit=63,file="kappa.dat")
open(unit=68,file="realeudiretta.dat")
open(unit=53,file="ham.dat")
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
!write(63,*)k(i)
end do

!numeno=163.84d0
!nupiu=5D-48

numeno=100d0
nupiu=1D-48

do i=1,n
f(i)=0.

nukm(i)=numeno*((dabs(k(i)))**(-8.))
nukp(i)=nupiu*((dabs(k(i)))**(16.))

if((i.eq.3).or.(i.eq.4).or.(i.eq.5).or.(i.eq.4091).or.(i.eq.4092).or.(i.eq.4093)) then

		f(i)=0.2d0

	end if
end do

u0=f

!!CALCOLO ESPONENZIALE OPERATORE LINEARE!!

esp=1.0d0*((dabs(k))**(0.5d0))

!!LOOP SUL TEMPO!!

do j=1,m

!!CALCOLO PEZZO OPERATORE NON LINEARE!!

u=cdexp(im*dt*(cdabs(u0))**2)*u0

!!FFT DELLA PARTE NON LINEARE!!

 call cfft1f(n,inc,u,lenc,wsave,lensav,work,lenwrk,ier)

!!MOLTIPLICO PER PARTE LINEARE!!

u=cdexp(im*dt*esp)*u

do i=n/4,n-n/4

	u(i)=0.d0	

end do

!!MOLTIPLICO PER FATTORE FORZANTE E DISSIPATIVO!!

u=dexp((f-nukm-nukp)*dt)*u

psiquadro=psiquadro+dt*((cdabs(u))**2)

!!ANTITRASFORMO CON FFT!!

 call cfft1b(n,inc,u,lenc,wsave,lensav,work,lenwrk,ier)

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
write(53,*)j*dt,ham
write(34,*)j*dt,nlin/lin

!!###########################################################!!

!if((j.eq.100000)) then
if(mod(j,nprints).eq.0) then

psiquadro=psiquadro/dfloat(j)

	do i=1,n
		write(46,*)k(i),psiquadro(i)
	end do
		write(46,*) 
	do i=1,n
		write(68,*)dreal(u(i))
	end do
		write(68,*) 
end if

u0=u

end do

 close(46)
 close(63)
 close(61)
 close(68)
 close(53)
 close(34)
deallocate(wsave)
deallocate(work)
end program mmt
