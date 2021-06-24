#F90=pgf90
F90=ifort
OBJECTd= paramfer.o dos.o
OBJECTdb= paramfer.o dos_big.o
OBJECTp= paramfer.o dospar.o
dos:$(OBJECTd)
	$(F90) -o dos $(OBJECTd)
dos_big:$(OBJECTdb)
	$(F90) -o dos_big $(OBJECTdb)
dospar:$(OBJECTp)
	$(F90) -o dospar $(OBJECTp)
dos.o:dos.f90
	$(F90) -c -g dos.f90
dos_big.o:dos_big.f90
	$(F90) -c -g dos_big.f90
dospar.o:dospar.f90
	$(F90) -c -g dospar.f90
paramfer.o:paramfer.f90
	$(F90) -c -g paramfer.f90
clean:
	rm -f dos dospar *.o *.mod 
