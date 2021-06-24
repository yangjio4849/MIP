MODULE params 
  implicit none
  integer,PARAMETER :: dp=kind(1.0d0)
  real(dp),PARAMETER :: tol8= 0.00000001_dp
  real(dp),PARAMETER :: pi=3.141592653589793238462643383279502884197_dp
  real(dp),PARAMETER :: zero=0._dp
  REAL(dp),PARAMETER   :: kb  = 8.6174103E-5_dp     ! eV/K
  REAL(dp),PARAMETER   :: hbar=1.05457148E-34_dp    ! Js
  REAL(dp),PARAMETER   :: e   = 1.602176462E-19_dp  ! ELECTRON CHARGE [C]
  REAL(dp),PARAMETER   :: me     = 9.1093826E-31_dp    ! ELECTRON MASS   [kg]
  REAL(dp),PARAMETER   :: A0_SI        = 0.5291772083E-10_dp ! Bohr radius [m]
  REAL(dp),PARAMETER   :: auprAA       = 0.5291772083_dp     ! Bohr radius [A]
  REAL(dp),PARAMETER   :: auprkg       = 1.66053886E-27_dp   ! atomic mass unit [kg] 
END MODULE params


MODULE Fermifunction
  use params
  implicit none

  CONTAINS 
REAL(dp) FUNCTION DFERMI(E,EF,T,spi)
  REAL(dp) :: e,ef,t
  REAL(dp) :: factor,efact
  integer spi
  factor=(e-ef)/(kb*T)
  if(factor>40) then
    DFERMI=0.0
  else if (factor<-40) then
    DFERMI=0.0
  else
    efact=EXP(factor)
    DFERMI=efact/(kb*T*(1.0+efact)**2)*2/spi
  endif
END FUNCTION DFERMI

REAL(dp) FUNCTION FERMI(E,EF,T,spi)
  REAL(dp) :: e,ef,t
  REAL(dp) :: factor,efact
  integer spi
  factor=(e-ef)/(kb*T)
  if(factor>40) then
    FERMI=0.0
  else if (factor<-40) then
    FERMI=2.0/spi
  else
    efact=EXP(factor)
    FERMI=2.0/spi/(efact+1)
  endif
END FUNCTION FERMI  
END MODULE Fermifunction 

MODULE variable 
  use params 
  implicit none
  real(dp),save:: totale,T,efermi
  real(dp),save,allocatable::bandenergy(:,:,:),kwpt(:)
  real(dp),save:: scales,vec(3,3),volume
  integer i,j,k,l,eof
  integer ik,iband,ispin,ivec
  integer,save:: nk,nband,nspin
end MODULE variable

module ef
use variable
use Fermifunction
implicit none
interface efdetermin
  module procedure getef 
end interface

CONTAINS
SUBROUTINE getef(ferminorigin,fermaxorigin)
implicit none
real(kind=dp),intent(in),optional::ferminorigin,fermaxorigin
real(kind=dp) fermin,fermax
real(kind=dp) de
real(kind=dp) ne(nspin),netemp(nspin)
character*1 nomeaning

de=1.0
fermin=2.0
fermax=20.0
if (present(ferminorigin).and.present(fermaxorigin)) then
fermin=ferminorigin
fermax=fermaxorigin
endif
do while(de.GE.tol8)
  efermi=(fermin+fermax)/2.0
  ne=0.0
  do ik=1,nk
    netemp=0.0
    do iband=1,nband
      do ispin=1,nspin
        netemp(ispin)=netemp(ispin)+FERMI(bandenergy(ispin,iband,ik),efermi,T,nspin)
      enddo
    enddo
    do ispin=1,nspin
      ne(ispin)=ne(ispin)+netemp(ispin)*kwpt(ik)
    enddo
  enddo
  de=abs(sum(ne)-totale)
  if(sum(ne)>totale)fermax=efermi
  if(sum(ne)<totale)fermin=efermi
enddo
return
END SUBROUTINE getef 
END MODULE ef


