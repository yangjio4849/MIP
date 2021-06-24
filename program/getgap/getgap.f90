program getgap2017
implicit none
character*80::nonsense
logical::alive
real(kind=8),allocatable::bandenergy(:,:,:)
integer nspin,nband,nk,ispin,iband,ik,ispincbm,ibandcbm,ikcbm,ispinvbm,ibandvbm,ikvbm
integer i,j,k,eof
real(kind=8)::efermi, totale , cbm , vbm 
real(kind=8),allocatable::kpt(:,:)

real(kind=8)::emin,emax,gapgama,gaptotal,emintemp,emaxtemp,bandshift
real(kind=8)::vbmx,vbmy,vbmz,cbmx,cbmy,cbmz,tbmx,tbmy,tbmz


inquire(file="fermilevel",exist=alive)

if(alive) then
  open(unit=38,file="fermilevel",status="old")
  read(38,*)efermi
  close(38)
else
  open(unit=38,file="OUTCAR",status="old",iostat=eof)
  if(eof/=0) stop "OUTCAR not found!"
  out:do while(.TRUE.)
    read(38,*)nonsense
    if(nonsense=="E-fermi") then
      backspace(38)
      read(38,*) nonsense,nonsense,efermi
      exit out
    endif
  end do out
close(38)
endif

open(unit=38,file="EIGENVAL",status="old",iostat=eof)
if (eof/=0) stop "EIGENVAL not found!"
read(38,*)nonsense,nonsense,nonsense,nspin
do i=1,4
  read(38,*) nonsense
enddo
read(38,*,iostat=eof) totale,nk,nband
if (eof/=0) stop "Too many kpoints. Number of electrons and kpoints should be seperated in EIGENVAL (line 6)! "
allocate(kpt(3,nk),bandenergy(nspin,nband,nk))

do ik=1,nk
  read(38,*)kpt(:,ik)
  do iband=1,nband
    read(38,*)nonsense,(bandenergy(ispin,iband,ik),ispin=1,nspin)
  enddo
enddo
close(38)

cbm=1E6
vbm=-1E6

ispincbm=1
ibandcbm=1
ikcbm=1
ispinvbm=1
ibandvbm=1
ikvbm=1

do ik=1,nk
  do iband=1,nband
    do ispin=1,nspin
      if (bandenergy(ispin,iband,ik)>efermi.and.bandenergy(ispin,iband,ik)<cbm) then
        cbm=bandenergy(ispin,iband,ik)
        ispincbm=ispin
        ibandcbm=iband
        ikcbm=ik
      endif
      if (bandenergy(ispin,iband,ik)<efermi.and.bandenergy(ispin,iband,ik)>vbm) then
        vbm=bandenergy(ispin,iband,ik)
        ispinvbm=ispin
        ibandvbm=iband
        ikvbm=ik
      endif
    enddo
  enddo
enddo

open(unit=40,file="gap.txt",status="replace")
write(40,"('The k-point of VBM is (',F9.7,',',F9.7,',',F9.7,'), with value ',F12.8,' , band no. ',I4)")kpt(:,ikvbm),vbm,ibandvbm
write(40,"('The k-point of CBM is (',F9.7,',',F9.7,',',F9.7,'), with value ',F12.8,' , band no. ',I4)")kpt(:,ikcbm),cbm,ibandcbm
write(40,"('The Fermilevel is ',F12.8)") efermi
write(40,"('The distance between VBM and fermilevel is ',F12.8)")efermi-vbm
write(40,"('The distance between CBM and fermilevel is ',F12.8)")cbm-efermi
write(40,"('The gap is ',F12.8)")cbm-vbm
write(*,"(F12.8)")cbm-vbm
close(40)

end
