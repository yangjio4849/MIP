program getgap_simple
use ef
implicit none 
character*80 nomeaning
real(dp) temp,vbm,cbb
real(dp),allocatable::vb(:),cb(:),band(:)
integer nvb,ncb

open(unit=38,file="EIGENVAL",status="old",iostat=eof)
if (eof/=0) stop "EIGENVAL not found!"
read(38,*)nomeaning,nomeaning,nomeaning,nspin
do i=1,4
  read(38,"(A80)") nomeaning
enddo
read(38,*) totale,nk,nband

allocate(bandenergy(nspin,nband,nk),band(nspin*nband*nk))

do ik=1,nk
  read(38,*)nomeaning
  do iband=1,nband
    read(38,*)nomeaning,(bandenergy(ispin,iband,ik),ispin=1,nspin)
  enddo
enddo
close(38)

open(unit=38,file="OUTCAR",status="old",iostat=eof)
if(eof/=0) stop "OUTCAR not found!"
out:do while(.TRUE.)
  read(38,*)nomeaning
  if(nomeaning=="E-fermi") then
    backspace(38)
    read(38,*) nomeaning,nomeaning,efermi
    exit out
  endif
end do out
close(38)

vbm=minval(bandenergy)
cbb=maxval(bandenergy)

do ik=1,nk
  do iband=1,nband
    do ispin=1,nspin
      temp=bandenergy(ispin,iband,ik)
      if (temp<efermi.and.temp>vbm) vbm=temp
      if (temp>efermi.and.temp<cbb) cbb=temp
    enddo
  enddo
enddo

write(*,"(F12.8,A)") vbm,"  vbm"
write(*,"(F12.8,A)") cbb,"  cbb"
write(*,"(F12.8,A)") cbb-vbm,"  gap"

!band=pack(bandenergy,bandenergy>-1000000_dp)
!nvb=count(band<0_dp)
!ncb=count(band>0_dp)
!allocate(vb(nvb),cb(ncb))
!vb=pack(band,band<0_dp)
!cb=pack(band,band>0_dp)

!write(*,"(A,F10.6)") minval(cb)-maxval(vb)
end
