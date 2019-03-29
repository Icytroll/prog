CFLAGS = -Wall -std=gnu1x
LDLIBS = -lm

default: Output/A.txt Output/B.txt Output/C.txt

Output/A.txt Output/B.txt Output/C.txt: A.txt B.txt C.txt
	mv A.txt Output/A.txt
	mv B.txt Output/B.txt
	mv C.txt Output/C.txt

A.txt B.txt C.txt: main
	./$<

main: main.o A_decomp.o A_solve.o B_inverse.o bi_inverse.o mv_mult.o C_bidiag.o bi_solve.o bi_det.o qr_gs_decomp.o qr_gs_solve.o qr_gs_inverse.o bi_decomp.o vector.o vector.h matrix.o matrix.h

A_solve.o C_bidiag.o vector.o qr_gs_solve mv_mult.o bi_inverse.o: vector.h
A_solve.o A_decomp.o B_inverse.o C_bidiag.o matrix.o mv_mult.o bi_inverse.o qr_gs_decomp bi_det.o bi_solve.o qr_gs_solve bi_decomp.o: matrix.h

.PHONEY:clean
clean:
	$(RM) main *.o Output/*.txt
