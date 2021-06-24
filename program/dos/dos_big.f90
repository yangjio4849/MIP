program dos_big 
! Calculate the dos 
use ef
implicit none
integer nemu,nfiles,ifiles, itotalk
integer,allocatable::nkpoints(:)
real(dp) emu,estart,de,zeropoint,fermi0
real(dp),allocatable::fermi_file(:)
integer control 
integer sisband
real(dp) sis,markscale,dfdeT,dfdescale
character*80 nomeaning,outputformat

open(unit=38,file="EIGENVAL_COMBINED",status="old",iostat=eof)
if (eof/=0) stop "EIGENVAL_COMBINED not found!"
read(38,*) fermi0, nband, nspin, nfiles  
allocate(nkpoints(nfiles), fermi_file(nfiles))
do ifiles=1,nfiles
  read(38,*) nkpoints(ifiles), fermi_file(ifiles)
enddo

nk=sum(nkpoints)
allocate(bandenergy(nspin,nband,nk),kwpt(nk))
kwpt=1.0/nk

!read the data from EIGENVAL_COMBINED, and align the fermi level band after band
itotalk=1
do ifiles=1,nfiles
  do ik=1,nkpoints(ifiles)
    read(38,*)nomeaning
    do iband=1,nband
      read(38,*) nomeaning,(bandenergy(ispin,iband,itotalk),ispin=1,nspin)
      bandenergy(:,iband,itotalk)=bandenergy(:,iband,itotalk)-fermi_file(ifiles)+fermi0
    enddo
    itotalk=itotalk+1
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
write(*,"(A)")"You want scissors for band gaps?"
write(*,"(A)")"Type 0 for no, 1 for yes"
read(*,*)control
if(control<0.or.control>1) call error() 
sis=0.0
sisband=1
if(control==1) then
  write(*,"(A)")"Tell me the band index, including and below which the band will be addded a constant value"
  read(*,*)sisband
  write(*,"(A)")"Tell me the value of scissor difference"
  read(*,*)sis
endif

bandenergy(:,1:sisband,:)=bandenergy(:,1:sisband,:)+sis
  
call getdos(zeropoint,estart,de,nemu)
end

SUBROUTINE error()
implicit none
pause "No kidding.  Please any key to quit..."
stop
end subroutine error

SUBROUTINE getdos(zeropoint,estart,de,nemu)
use variable
use Fermifunction
implicit none
real(dp) dos(nspin),dostemp(nspin),intdos(nspin),intdostemp(nspin),dosvk(nspin),dosvktemp(nspin)
real(dp) vk2(nspin,nband,nk),vk(3,nspin)
real(dp) emu,zeropoint,estart,de,realef
character*80 nomeaning,dosformat
integer nemu
logical putvk

putvk=.FALSE.
if(putvk) then
open(unit=38,file="GROUPVEC",status="old",iostat=eof)
if(eof==0) then
  vk2=0.0
!  putvk=.TRUE.
  do ik=1,nk
    read(38,*)nomeaning
    do iband=1,nband
      read(38,*)nomeaning,((vk(i,ispin),ispin=1,nspin),i=1,3)
      do ispin=1,nspin
        vk2(ispin,:,:)=sum(vk(:,ispin))
      enddo  
    enddo
  enddo
  close(38)
  vk2=vk2*(1.0E-10*hbar/e/me)**2
endif
endif

open(unit=38,file="dos.txt",status="replace")
if(putvk) then
  write(dosformat,"(A,I1,A)")"(F10.6,",2*nspin+1,"G17.8)"
else
  write(dosformat,"(A,I1,A)")"(F10.6,",2*nspin,"G17.8)"
endif
do i=1,nemu
  realef=zeropoint+(i-1)*de+estart
  dos=0.0
  intdos=0.0
  if(putvk)dosvk=0.0
  do ik=1,nk
    dostemp=0.0
    intdostemp=0.0
    if(putvk)dosvktemp=0.0
    do iband=1,nband
      do ispin=1,nspin
        dostemp(ispin)=dostemp(ispin)+DFERMI(bandenergy(ispin,iband,ik),realef,T,nspin)
        intdostemp(ispin)=intdostemp(ispin)+FERMI(bandenergy(ispin,iband,ik),realef,T,nspin)
        if(putvk)dosvktemp(ispin)=dosvktemp(ispin)+DFERMI(bandenergy(ispin,iband,ik),realef,T,nspin)*sqrt(vk2(ispin,iband,ik))
      enddo
    enddo
    do ispin=1,nspin
      dos(ispin)=dos(ispin)+dostemp(ispin)*kwpt(ik)
      intdos(ispin)=intdos(ispin)+intdostemp(ispin)*kwpt(ik)
      if(putvk)dosvk(ispin)=dosvk(ispin)+dosvktemp(ispin)*kwpt(ik)
    enddo
  enddo
  if(putvk)then
    write(38,dosformat) (i-1)*de+estart,dos,intdos,sum(dosvk)**2.0/sum(dos)
  else
    write(38,dosformat) (i-1)*de+estart,dos,intdos
  endif
enddo
close(38)  
END SUBROUTINE getdos

