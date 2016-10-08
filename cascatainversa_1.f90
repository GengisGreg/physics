!!!!!ALGORITMO RISOLUZIONE EQUAZIONE MMT CON alpha=1/2 e beta=0, UTILIZZANDO SPLIT-STEP METHOD!!!!!
!!UTILIZZO LIBRERIE FFTPACK5 gfortran nomefile -lfftpack5.1d!!

program mmt
implicit none

integer i,j,ier,inc,lenc,lensav,lenwrk,seed
integer, parameter::n=4096,m=100000,nprints=10
real*8, parameter::pi=2.0d0*dasin(1.0d0),l=2.0d0*pi
real*8 dt,dk,k(n),h,esp(n),f(n),numeno,nupiu,nukm(n),nukp(n),psiquadro(n),deltau(n),dx,xmin,xmax
complex*16 u(n),psi(n),u0(n)
complex*16,parameter::im=(0.d0,1.d0)
real*8,allocatable,dimension(:):: work
real*8,allocatable,dimension(:):: wsave

!!cose per l'uso della libreria fftpack5!!
lensav=2*n+int(log(real(n))/log(2.0d0))+4
lenwrk=2*n
allocate(work(1:lenwrk))
allocate(wsave(1:lensav))

dt=0.05d0
h=l/real(n)
dk=2.0d0*pi/(l)

open(unit=46,file="cascatainv.dat")
open(unit=63,file="kappa.dat")
open(unit=61,file="realeuinversa.dat")

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

end do
numeno=1D4
nupiu=1D-27

do i=1,n
f(i)=0.
nukm(i)=numeno*((dabs(k(i)))**(-16.))
nukp(i)=nupiu*((dabs(k(i)))**(8.))

if((i.eq.1000).or.(i.eq.1001).or.(i.eq.1002).or.(i.eq.3095).or.(i.eq.3096).or.(i.eq.3097)) then
		
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

!!MOLTIPLICO PER FATTORE FORZANTE E DISSIPATIVO!!

u=dexp((f-nukm-nukp)*dt)*u

psiquadro=psiquadro+dt*((cdabs(u))**2)

!!ANTITRASFORMO CON FFT!!

 call cfft1b(n,inc,u,lenc,wsave,lensav,work,lenwrk,ier)
!if((j.eq.100000)) then

if(mod(j,nprints).eq.0) then

psiquadro=psiquadro/dfloat(j)

do i=1,n
	write(46,*)k(i),psiquadro(i)
!write(63,*)i,k(i)
end do
	write(46,*) 
do i=1,n
		write(61,*)dreal(u(i))
	end do
		write(61,*)  
end if

u0=u

end do

 close(61)
 close(46)
 close(63)
deallocate(wsave)
deallocate(work)
end program mmt

