PROGRAM chkpt 
  implicit none
  character(len=100) nomeaning
  integer nirkpt,niraddkpt,nkpt,nkptori,nsym,nbands,eof
  integer ik,iband,i,j,k,iirk,isym
  integer,allocatable::corresponding(:)  ! corresponding(:) for each k-point out of the whole pool of nkpt, the origin one from irkpt
  real(kind=8),allocatable::sym(:,:,:),eig(:,:),occ(:)    ! eig(nbands,nirkpt),occupation(nbands)
  real(kind=8),allocatable::irkpt(:,:),kpt(:,:)
  real(kind=8) kpttemp(3),kpttempijk(3),err,temp(3,3)
  logical new

! read SYMMETRY
open(unit=38,file="SYMMETRY")
read(38,*)nsym
allocate(sym(3,3,nsym*2))
do isym=1,nsym
  do j=1,3
    read(38,*)(sym(j,k,2*isym-1),k=1,3)
  enddo
    sym(:,:,2*isym-1)=transpose(sym(:,:,2*isym-1))
    sym(:,:,2*isym)=sym(:,:,2*isym-1)*(-1.0)
enddo
nsym=nsym*2
close(38)

! read EIGENVAL
open(unit=38,file="EIGENVAL",status="old")
read(38,*)nomeaning
read(38,*)nomeaning
read(38,*)nomeaning
read(38,*)nomeaning
read(38,*)nomeaning
read(38,"(A100)")nomeaning
read(nomeaning,*,iostat=eof)ik,niraddkpt,nbands
if (eof/=0) stop "Too many kpoints. Number of electrons and kpoints should be seperated in EIGENVAL (line 6)! "
allocate(irkpt(4,niraddkpt),eig(nbands,niraddkpt),kpt(3,nsym*niraddkpt*2),corresponding(nsym*niraddkpt*2)) 
irkpt=0
do iirk=1,niraddkpt
  read(38,*)(irkpt(j,iirk),j=1,4)   
  do iband=1,nbands
    read(38,*)nomeaning,eig(iband,iirk)
  enddo
enddo
close(38)

forall (i=1:3,j=1:niraddkpt,abs(irkpt(i,j))<1E-4)
  irkpt(i,j)=0
end forall

if (abs(irkpt(1,1))+abs(irkpt(2,1))+abs(irkpt(3,1))>1E-4) stop "The KPOINTS should be Gamma centered (G instead of M)"

nirkpt=count(irkpt(4,:)>1E-8)

! main loop  
nkpt=0
kpt=0.0
corresponding=0
do iirk=1,nirkpt
  do isym=1,nsym
    kpttemp(1)=irkpt(1,iirk)*sym(1,1,isym)+irkpt(2,iirk)*sym(1,2,isym)+irkpt(3,iirk)*sym(1,3,isym)
    kpttemp(2)=irkpt(1,iirk)*sym(2,1,isym)+irkpt(2,iirk)*sym(2,2,isym)+irkpt(3,iirk)*sym(2,3,isym)
    kpttemp(3)=irkpt(1,iirk)*sym(3,1,isym)+irkpt(2,iirk)*sym(3,2,isym)+irkpt(3,iirk)*sym(3,3,isym)
! stupid 0,1 problem Jiong Yang
    do i=0,1
      do j=0,1
        do k=0,1
          new=.TRUE.
          kpttempijk(1)=kpttemp(1)+i
          kpttempijk(2)=kpttemp(2)+j
          kpttempijk(3)=kpttemp(3)+k
          do while (kpttempijk(1)<0.0) 
            kpttempijk(1)=kpttempijk(1)+1
          enddo
          do while (kpttempijk(1)-1.0>0.0)
            kpttempijk(1)=kpttempijk(1)-1
          enddo
          do while (kpttempijk(2)<0.0) 
            kpttempijk(2)=kpttempijk(2)+1
          enddo
          do while (kpttempijk(2)-1.0>0.0)
            kpttempijk(2)=kpttempijk(2)-1
          enddo
          do while (kpttempijk(3)<0.0) 
            kpttempijk(3)=kpttempijk(3)+1
          enddo
          do while (kpttempijk(3)-1.0>0.0)
            kpttempijk(3)=kpttempijk(3)-1
          enddo
!          kpttempijk(4)=1000000*kpttempijk(1)+1000*kpttempijk(2)+kpttempijk(3)
          do ik=1,nkpt
            err=abs(kpttempijk(1)-kpt(1,ik))+abs(kpttempijk(2)-kpt(2,ik))+abs(kpttempijk(3)-kpt(3,ik))   ! if use kpttempijk(4), the nkpt will be overestimated
            if (err<=1E-4) then
              new=.FALSE.
              exit
            endif
          enddo
          if (new) then
            nkpt=nkpt+1
            kpt(:,nkpt)=kpttempijk(:)
            corresponding(nkpt)=iirk
          endif
        enddo   ! k
      enddo   ! j
    enddo   !i
  enddo  ! isym
enddo  ! iirk

! for additional k-points when you use KPOINTS.recommend from checkresults and also check the duplications! I guess you won't need it!
do iirk=nirkpt+1,niraddkpt
  do ik=1,nkpt
    new=.true.
    err=abs(irkpt(1,iirk)-kpt(1,ik))+abs(irkpt(2,iirk)-kpt(2,ik))+abs(irkpt(3,iirk)-kpt(3,ik)) 
    if (err<=1E-4) then
      new=.FALSE.
      exit
    endif
  enddo
  if (new) then
    nkpt=nkpt+1
    kpt(:,nkpt)=irkpt(1:3,iirk)
    corresponding(nkpt)=iirk
  endif
enddo

! sort the kpt and corresponding in ascending order, bubbling method
do i=1,nkpt-1
  do j=i+1,nkpt
    if (kpt(1,i)-kpt(1,j)>1E-4.or.abs(kpt(1,i)-kpt(1,j))<1E-4.and.kpt(2,i)-kpt(2,j)>1E-4.or.abs(kpt(1,i)-kpt(1,j))<1E-4.and.abs(kpt(2,i)-kpt(2,j))<1E-4.and.kpt(3,i)-kpt(3,j)>1E-4) then
      kpttemp=kpt(:,i)
      k=corresponding(i)
      kpt(:,i)=kpt(:,j)
      corresponding(i)=corresponding(j)
      kpt(:,j)=kpttemp
      corresponding(j)=k
    endif
  enddo
enddo


open(unit=40,file="full_eig",status="replace")
write(40,"(A,I10,A,I10)")"NKPTS = ",nkpt,"  ,  NBANDS = ",nbands
open(unit=41,file="OUTCAR",status="old")
do while (.TRUE.) 
  read(41,"(A100)")nomeaning
  if(nomeaning(1:71)=="      direct lattice vectors                 reciprocal lattice vectors") then
    write(40,"(A71)") nomeaning
    do i=1,3
      read(41,"(A84)") nomeaning
      write(40,"(A84)") nomeaning
    enddo
  endif
  if (nomeaning(1:8)==" E-fermi") then
    write(40,"(A65)")nomeaning
    exit
  endif
enddo
close(41)
open(unit=42,file="klist",status="replace")
do ik=1,nkpt
  write(40,"(3F20.14)") kpt(1:3,ik)
  write(42,"(3F20.14)") kpt(1:3,ik)
  do iband=1,nbands
    write(40,"(I5,5X,F30.10,F20.5)")iband,eig(iband,corresponding(ik)),1.00
  enddo
  write(40,*)""
enddo

close(40)
close(42)

END PROGRAM 

subroutine inverse(a)
real(kind=8)::a(3,3),b(3,3)
real(kind=8)::debt
debt=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(2,1)*a(3,2)*a(1,3)-a(1,3)*a(2,2)*a(3,1)-a(1,2)*a(2,1)*a(3,3)-a(1,1)*a(2,3)*a(3,2)
b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
b(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
b(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
b(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
b(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
b(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
b(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2)
b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
b=b/debt
a=b
end subroutine

