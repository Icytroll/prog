CFLAGS = -Wall -std=gnu1x
LDLIBS = -lm

default: Output/A.txt Output/B.txt

Output/A.txt Output/B.txt: A.txt B.txt
	mv A.txt Output/A.txt
	mv B.txt Output/B.txt

A.txt B.txt: main
	./$<

main: main.o A_decomp.o A_solve.o B_inverse.o qr_gs_decomp.o qr_gs_solve.o qr_gs_inverse.o vector.o vector.h matrix.o matrix.h

A_solve.o vector.o qr_gs_solve: vector.h
A_solve.o A_decomp.o B_inverse.o matrix.o qr_gs_decomp qr_gs_solve: matrix.h

.PHONEY:clean
clean:
	$(RM) main *.o Output/*.txt
