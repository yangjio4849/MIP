MODULE variables 
  use params
  implicit none
  real(dp),allocatable::parsimp(:,:,:,:,:)   ! simplified edition for par
  integer nion,iion,nlm,ilm  ! number of ions; number of lm-orbit; counting number
  integer ionstart,ionend
  integer ncol,icol                           ! ncol=4 or 5, (including total) depends on f orbitals
END MODULE variables

program dospar 
  use ef
  use variables
  implicit none
  integer nemu
  real(dp) emu,estart,de,zeropoint
  integer sisband
  real(dp) sis
  character*150 nomeaning,outputformat
  real(dp),allocatable::par(:)            ! nlm
  logical yon,splitd,SOC                  ! yes or no  ,  d eletrons split to t2g & eg, spin-orbit interaction

open(unit=38,file="EIGENVAL",status="old",iostat=eof)
if (eof/=0) stop "EIGENVAL not found!"
read(38,*)nomeaning,nomeaning,nomeaning,nspin
do i=1,4
  read(38,*) nomeaning
enddo
read(38,*) totale,nk,nband

allocate(bandenergy(nspin,nband,nk),kwpt(nk))
do ik=1,nk
  read(38,*)nomeaning,nomeaning,nomeaning,kwpt(ik)
  do iband=1,nband
    read(38,*)nomeaning,(bandenergy(ispin,iband,ik),ispin=1,nspin)
  enddo
enddo
close(38)

! get the input info iteratively 
write(*,"(A)")"The energy for zero point is:"
read(*,*)zeropoint
write(*,"(A)")"The energy starts from (with respect to zero point):"
read(*,*)estart
write(*,"(A)")"The energy step is:"
read(*,*)de
write(*,"(A)")"How many points in this energy step?"
read(*,*)nemu
write(*,"(A)")"Tell the temperature entering the Fermi-Dirac function:"
read(*,*)T
write(*,"(A)")"What scissors for band gap?(T/F)"
read(*,*)yon
sisband=1
sis=0.0
if (yon) then
  write(*,"(A)")"Tell me the band index, including and below which the band will be addded a constant value (i.e. VBM)"
  read(*,*)sisband
  write(*,"(A)")"Tell me the value of scissor difference, bandenergy(1:sisband)=bandenergy(1:sisband)+sis"
  read(*,*)sis
endif
write(*,"(A)")"Does the calculation has spin-orbit interaction (T/F)? Just read the first component with SOC"
read(*,*)SOC

bandenergy(:,1:sisband,:)=bandenergy(:,1:sisband,:)+sis
  
call gettotaldos(zeropoint,estart,de,nemu)

open(unit=38,file="PROCAR",status="old",iostat=eof)
if (eof/=0) stop "PROCAR not found!"

splitd=.FALSE.
write(*,"(A)") "Want to split the d orbitals(i.e. t2g & eg)? (T/F)"
read(*,*) splitd

read(38,*)nomeaning
read(38,"(A14,A5,A20,A4,A19,I4)")nomeaning,nomeaning,nomeaning,nomeaning,nomeaning,nion
read(38,*)nomeaning
read(38,*)nomeaning
read(38,"(A)")nomeaning   !  The empty line
read(38,"(A122)")nomeaning
if (index(nomeaning,"f")==0) then
  nlm=9
  ncol=4           ! s,p,d+total
else
  nlm=16
  ncol=5           ! spdf+total
endif
if (splitd) ncol=ncol+1

allocate(par(nlm))
allocate(parsimp(ncol,nion,nband,nk,nspin))
rewind(38)
read(38,*)nomeaning   ! PROCAR lm decomposed
do ispin=1,nspin
  read(38,*)nomeaning !# of k-points:    4         # of bands: 140         # of ions:  33
  do ik=1,nk
    read(38,*)nomeaning  ! k-point    1 :    0.00000000 0.00000000 0.00000000     weight = ..
    do iband=1,nband
      read(38,*)nomeaning  !band   1 # energy  -33.64471171 # occ.  2.00000000
      read(38,*)nomeaning  !ion      s     py     pz     px  
      do iion=1,nion
        read(38,*)nomeaning,(par(ilm),ilm=1,nlm)
        do icol=1,ncol-1
          if (splitd) then
            select case(icol)
            case(1)
              parsimp(icol,iion,iband,ik,ispin)=par(1)
            case(2)
              parsimp(icol,iion,iband,ik,ispin)=sum(par(2:4))
            case(3)   ! d_t2g
              parsimp(icol,iion,iband,ik,ispin)=par(5)+par(6)+par(8)
            case(4)   ! d_eg
              parsimp(icol,iion,iband,ik,ispin)=par(7)+par(9)
            case(5)
              parsimp(icol,iion,iband,ik,ispin)=sum(par(10:16))
            end select
          else
            if (icol==2) then
              parsimp(2,iion,iband,ik,ispin)=par(2) 
            else
              parsimp(icol,iion,iband,ik,ispin)=sum(par(icol*icol-2*icol+2:icol*icol))  ! very smart!!
            endif
          endif
        enddo
        parsimp(ncol,iion,iband,ik,ispin)=sum(parsimp(1:ncol-1,iion,iband,ik,ispin))
      enddo
      if (nion/=1) read(38,*)nomeaning   ! tot
      if (SOC) then
        do iion=1,3*(nion+1)
          read(38,*) nomeaning
        enddo
      endif
    enddo
  enddo
enddo
close(38)

write(*,"(A)") "Want partial dos for special atoms?(T/F)"
read(*,*)yon
do while (yon)
  write(*,"(A)") "Please give two numbers n,m, & I'll sum ion number from n to m"
  read(*,*)ionstart,ionend
  if (1<=ionstart.and.ionstart<=nion.and.1<=ionend.and.ionend<=nion) then
    call getspecialpdos(zeropoint,estart,de,nemu)
  else
    write(*,"(A)")"No kidding! Out of the ion range!"
    stop
  endif
  write(*,"(A)")"One more time? (T/F)"
  read(*,*)yon
enddo
  
write(*,"(A)") "Want partial dos for each type of atoms? (T/F)"
read(*,*) yon
if (yon) call getpardos(zeropoint,estart,de,nemu)

end

SUBROUTINE gettotaldos(zeropoint,estart,de,nemu)
use variable
use variables
use Fermifunction
implicit none
real(dp) dos(nspin),dostemp(nspin),intdos(nspin),intdostemp(nspin)
real(dp) emu,zeropoint,estart,de,realef
character*80 nomeaning,dosformat
integer nemu

open(unit=38,file="tdos.txt",status="replace")
write(dosformat,"(A,I1,A)")"(F10.6,",2*nspin,"G17.8)"
do i=1,nemu
  realef=zeropoint+(i-1)*de+estart
  dos=0.0
  intdos=0.0
  do ik=1,nk
    dostemp=0.0
    intdostemp=0.0
    do iband=1,nband
      do ispin=1,nspin
        dostemp(ispin)=dostemp(ispin)+DFERMI(bandenergy(ispin,iband,ik),realef,T,nspin)
        intdostemp(ispin)=intdostemp(ispin)+FERMI(bandenergy(ispin,iband,ik),realef,T,nspin)
      enddo
    enddo
    do ispin=1,nspin
      dos(ispin)=dos(ispin)+dostemp(ispin)*kwpt(ik)
      intdos(ispin)=intdos(ispin)+intdostemp(ispin)*kwpt(ik)
    enddo
  enddo
  write(38,dosformat) (i-1)*de+estart,dos,intdos
enddo
close(38)  
END SUBROUTINE gettotaldos

SUBROUTINE getpardos(zeropoint,estart,de,nemu)
use variable
use variables
use Fermifunction
implicit none
real(dp),allocatable::dos(:,:) ! dos(ncol,nspin)
real(dp) dostemp
real(dp) emu,zeropoint,estart,de,realef
character*80 nomeaning,dosformat
integer nemu

integer nol,species,sstart              ! sstart is the start number of species
integer,allocatable::num(:)             ! store the number of each type
character(len=4),allocatable::symbol(:) ! contain the symbol
character(len=2)::symname               ! contain the temp symbol
character(len=2)::symnum                ! contain the No. of same symbol

open(unit=40,file="POTCAR",status="old",iostat=i)
if(i/=0) stop "POTCAR failed!"
nol=0
species=0
do while(.TRUE.)
  read(40,*,iostat=i) nomeaning
  if(i/=0)exit
  nol=nol+1
  if (nomeaning=="TITEL") species=species+1
end do

allocate(symbol(species))
allocate(num(species))
rewind(40)
j=1
symbol="    "
do i=1,nol
  read(40,*)nomeaning
  if(nomeaning=="TITEL")then
    backspace(40)
    read(40,*)nomeaning,nomeaning,nomeaning,symname
      l=0
      do k=1,j-1
        if (symbol(k)(1:2)==symname) l=l+1
      end do
      if (l==0) then
        symbol(j)=symname
      else
        write(symnum,"(I2.2)") l
        symbol(j)=symname//symnum
      end if
    j=j+1
  end if
end do
close(40)

open(unit=39,file="POSCAR",status="old",iostat=i)
if(i/=0) stop "POSCAR failed!"
do i=1,5
  read(39,*) nomeaning
enddo
read(39,*,iostat=eof) num
if (eof/=0)  read(39,*) num
close(39)

write(dosformat,"(A,I2,A)")"(F12.6,",ncol*nspin,"G17.8)"

allocate(dos(ncol,nspin))

sstart=1
do i=1,species
  open(unit=38,file=trim(symbol(i))//".txt",status="replace")
  do j=1,nemu
    realef=zeropoint+(j-1)*de+estart
    dos=0
    do ispin=1,nspin
      do icol=1,ncol
        do ik=1,nk
          dostemp=0.0
          do iband=1,nband
            dostemp=dostemp+sum(parsimp(icol,sstart:sstart+num(i)-1,iband,ik,ispin))*DFERMI(bandenergy(ispin,iband,ik),realef,T,nspin)
          enddo   ! band
          dos(icol,ispin)=dos(icol,ispin)+dostemp*kwpt(ik)
        enddo     ! k-points
      enddo       ! column
    enddo         ! spin
    write(38,dosformat) (j-1)*de+estart,dos
  enddo           ! energy
  close(38)  
  sstart=sstart+num(i)
enddo             ! species
END SUBROUTINE getpardos



SUBROUTINE getspecialpdos(zeropoint,estart,de,nemu)
use variable
use variables
use Fermifunction
implicit none
real(dp),allocatable::dos(:,:) ! dos(ncol,nspin)
real(dp) dostemp
real(dp) emu,zeropoint,estart,de,realef
character*80 nomeaning,dosformat,istart,iend
integer nemu

write(dosformat,"(A,I2,A)")"(F12.6,",ncol*nspin,"G17.8)"

allocate(dos(ncol,nspin))

write(istart,"(I4)")ionstart
write(iend,"(I4)")ionend
do i=1,79
  if(istart(i:i)==" ")  istart(i:79)=istart(i+1:80)
  if(iend(i:i)==" ")  iend(i:79)=iend(i+1:80)
enddo

open(unit=38,file=trim(istart)//"-"//trim(iend)//"_pdos.txt",status="replace")
do j=1,nemu
  realef=zeropoint+(j-1)*de+estart
  dos=0
  do ispin=1,nspin
    do icol=1,ncol
      do ik=1,nk
        dostemp=0.0
        do iband=1,nband
          dostemp=dostemp+sum(parsimp(icol,ionstart:ionend,iband,ik,ispin))*DFERMI(bandenergy(ispin,iband,ik),realef,T,nspin)
        enddo   ! band
        dos(icol,ispin)=dos(icol,ispin)+dostemp*kwpt(ik)
      enddo     ! k-points
    enddo       ! column
  enddo         ! spin
  write(38,dosformat) (j-1)*de+estart,dos
enddo           ! energy
close(38)  
END SUBROUTINE getspecialpdos
