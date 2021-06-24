program lattice_parameters 

!The program gets the a,b,c,alpha,beta,gamma from the line 2 to line 5 of POSCAR
!file.

!annoucement begin
implicit none
real*8::vec(3,3)          !lattice vector in POSCAR
real*8::a,b,c                         !The length of three vectors
real*8::alph,beta,gama                !The three angles of every two vectors
real*8::scales                        !The scale in POSCAR
integer,allocatable::num(:)           !The number of atoms;
integer i,j,k                         !The counting device as usual
integer species                       !number of kinds
integer nol                           !the total line of POTCAR
integer numb                          !number of all atoms
character(len=15) remark              !The remark which will be written in the first line
character(len=5) symmetry             !The symmetry of the lattice
character(len=2),allocatable::symbol(:) !The elemental symbol
character(len=20)::non1,non2,non3     
character(len=1)::non
character(len=20) forma
logical dirt  ! respective coordinates,  if false, means cartisian

!announcement end
real*8,parameter::pi=3.14159265
!Begin reading data from poscar.input


open(unit=38,file="POSCAR",status="old")
read(38,*)remark
read(38,*)scales
do i=1,3
  read(38,*)vec(i,1),vec(i,2),vec(i,3)
end do

!End reading, then begin converting something
a=sqrt(vec(1,1)**2+vec(1,2)**2+vec(1,3)**2)
b=sqrt(vec(2,1)**2+vec(2,2)**2+vec(2,3)**2)
c=sqrt(vec(3,1)**2+vec(3,2)**2+vec(3,3)**2)

gama=acos((vec(1,1)*vec(2,1)+vec(1,2)*vec(2,2)+vec(1,3)*vec(2,3))/a/b)
beta=acos((vec(1,1)*vec(3,1)+vec(1,2)*vec(3,2)+vec(1,3)*vec(3,3))/a/c)
alph=acos((vec(3,1)*vec(2,1)+vec(3,2)*vec(2,2)+vec(3,3)*vec(2,3))/c/b)


a=a*scales
b=b*scales
c=c*scales
alph=alph*180/pi
beta=beta*180/pi
gama=gama*180/pi


write(*,60) a,b,c,alph,beta,gama
60 format(F14.7,F14.7,F14.7,F11.5,F11.5,F11.5)

end

