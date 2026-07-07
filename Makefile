CC	= gfortran
OPS = -ftrace=full
OBJ = \
obj/mod_const.o \
obj/mod_rfpmLLGS.o \
obj/mod_funcLLGS.o  \
obj/mod_mLLGS.o \
obj/mod_cmLLGS.o
OBJ2 = \
obj/mod_readfile.o \
obj/mod_fft.o  
.SUFFIXES: .f95

main: $(OBJ)
	$(CC) -o main.exe main.f95 $^

read: $(OBJ2)
	$(CC) -o read.exe cmplx.f95 $^

$(OBJ): obj/%.o: ./module/%.f95
	$(CC) -c $< -o $@

$(OBJ2): obj/%.o: ../module/%.f95
	$(CC) -c $< -o $@

clean:
	rm -f ../module/*.o *~
	rm -f ./obj/*.o *~
	rm -f ./obj/*.mod *~
	rm -f *.mod
	rm -f *.exe


