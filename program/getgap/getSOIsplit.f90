program getSOIsplit
implicit none
character*80::nonsense
logical::alive,split
integer nelec,nkp,nba,i,j,nomin,nomax,notemp,eof,nVBM              
real(kind=8)::non1,non2,non3,efermi
real(kind=8)::emin,emax,gapgama,gaptotal,emintemp,emaxtemp
real(kind=8)::vbmx,vbmy,vbmz,cbmx,cbmy,cbmz,tbmx,tbmy,tbmz,tbmxtemp,tbmytemp,tbmztemp
real(kind=8)::gapmin, splitmin, EVBM, EVBMm1, ECBM, ECBMp1, splittemp

gapmin=0.01   !  The minimum gap value that go check the splitting of VBM-(VBM-1) and CBM+1-CBM
splitmin=0.07 !  The minimum splitting value that is worth reporting

inquire(file="fermilevel",exist=alive)

if(alive) then
  open(unit=38,file="fermilevel",status="old")
  read(38,*) efermi
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
close (39)

open(unit=40,file="gap.txt",status="replace")
write(40,"('The k-point of VBM is (',F10.7,',',F10.7,',',F10.7,'), with value ',F12.8,', band no. ',I4)")vbmx,vbmy,vbmz,emin,notemp-1
write(40,"('The k-point of CBM is (',F10.7,',',F10.7,',',F10.7,'), with value ',F12.8,', band no. ',I4)")cbmx,cbmy,cbmz,emax,notemp
write(40,"('The Fermilevel is ',F12.8)") efermi
write(40,"('The distance between VBM and fermilevel is ',F12.8)")efermi-emin
write(40,"('The distance between CBM and fermilevel is ',F12.8)")emax-efermi
write(40,"('The gap is ',F12.8)")emax-emin
close(40)

if (emax-emin>gapmin) then
  nVBM=notemp-1
  open(unit=39,file="EIGENVAL",status="old")
  read(39,*)nonsense
  read(39,*)nonsense
  read(39,*)nonsense
  read(39,*)nonsense
  read(39,*)nonsense
  read(39,*) nelec,nkp,nba
  split=.false.
  do i=1,nkp
    read(39,*)tbmxtemp,tbmytemp,tbmztemp
    do j=1,nVBM-2
      read(39,*)nonsense
    enddo
    read(39,*)notemp,EVBMm1
    read(39,*)notemp,EVBM
    read(39,*)notemp,ECBM
    read(39,*)notemp,ECBMp1
    do j=nVBM+3,nba
      read(39,*)nonsense
    enddo
    if (max(EVBM-EVBMm1,ECBMp1-ECBM)>splitmin) then
      split=.true.
      splitmin=max(EVBM-EVBMm1,ECBMp1-ECBM)
      tbmx=tbmxtemp
      tbmy=tbmytemp
      tbmz=tbmztemp
    endif
  end do
  if (split) then
    write(*,"('The maximum splitting is at (',F10.7,',',F10.7,',',F10.7,'), with value ',F12.8)")tbmx,tbmy,tbmz,splitmin
  endif
  close (39)
endif

end
