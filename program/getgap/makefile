#F90=pgf90
#F90=gfortran
F90=ifort
getgap:getgap.f90
	$(F90) -o getgap getgap.f90
getgap_simple:paramfer.o getgap_simple.o
	$(F90) -o getgap_simple getgap_simple.o paramfer.o
paramfer.o:paramfer.f90
	$(F90) -c paramfer.f90
getgap_simple.o:getgap_simple.f90
	$(F90) -c getgap_simple.f90
clean:
	rm -f getgap getgap_simple *.mod *.o
