program getgapfinal
implicit none
character*80::nonsense
logical::alive
integer nelec,nkp,nba,i,j,nomin,nomax,notemp,shift,eof             ! Electrons, kpoints and bands in one cell, respectively.
real(kind=8)::non1,non2,non3,efermi
real(kind=8)::emin,emax,gapgama,gaptotal,emintemp,emaxtemp,bandshift
real(kind=8)::vbmx,vbmy,vbmz,cbmx,cbmy,cbmz,tbmx,tbmy,tbmz


inquire(file="fermilevel",exist=alive)

shift=1                                                  !  if shift==1, no band shift
if(alive) then
  open(unit=38,file="fermilevel",status="old")
  read(38,"(A80)")nonsense
  read(nonsense,*,iostat=shift)efermi,bandshift
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

open(unit=39,file="EIGENVAL",status="old")
emax=100.0
emin=-100.0
read(39,*)nonsense
read(39,*)nonsense
read(39,*)nonsense
read(39,*)nonsense
read(39,*)nonsense
read(39,*) nelec,nkp,nba

if(shift==0) then                   !  shift the EIGENVAL below the efermi
  rewind(39)
  open(unit=40,file="EIGENVAL1",status="replace")
  do i=1,7
    read(39,"(A80)")nonsense
    write(40,"(A80)")nonsense
  enddo
  do i=1,nkp
    read(39,"(A80)")nonsense
    write(40,"(A80)")nonsense
    do j=1,nba
      read(39,*)notemp,non1
      if(non1.LT.efermi)non1=non1+bandshift
      write(40,"(I4,F14.4)")notemp,non1
    enddo
    if(i==nkp)cycle
    read(39,"(A80)")nonsense
    write(40,"(A80)")nonsense
  enddo 
stop "EIGENVAL1 generated!"
endif


do i=1,nkp
  read(39,*)tbmx,tbmy,tbmz
  read(39,*)notemp,emaxtemp
  do while(emaxtemp<efermi)
    emintemp=emaxtemp
    read(39,*)notemp,emaxtemp
  end do
  do j=1,nba-notemp
    read(39,*)nonsense
  end do
  if(emaxtemp<emax) then
    emax=emaxtemp
    cbmx=tbmx
    cbmy=tbmy
    cbmz=tbmz
  end if
  if(emintemp>emin) then
    emin=emintemp
    vbmx=tbmx
    vbmy=tbmy
    vbmz=tbmz
  end if
end do
close(39)

open(unit=40,file="gap.txt",status="replace")
write(40,"('The k-point of VBM is (',F9.7,',',F9.7,',',F9.7,'), with value ',F12.8,' , band no. ',I4)")vbmx,vbmy,vbmz,emin,notemp-1
write(40,"('The k-point of CBM is (',F9.7,',',F9.7,',',F9.7,'), with value ',F12.8,' , band no. ',I4)")cbmx,cbmy,cbmz,emax,notemp
write(40,"('The Fermilevel is ',F12.8)") efermi
write(40,"('The distance between VBM and fermilevel is ',F12.8)")efermi-emin
write(40,"('The distance between CBM and fermilevel is ',F12.8)")emax-efermi
write(40,"('The gap is ',F12.8)")emax-emin
write(*,*)emax-emin
close(40)

end
