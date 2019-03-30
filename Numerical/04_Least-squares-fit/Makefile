CFLAGS = -Wall -std=gnu1x 
LDLIBS = -lm

default: A.svg

A.svg: A.gpi A.txt
	gnuplot $< > $@

A.txt: main
	./$<

main: main.o lsfit.o inverse.o qr_gs_inverse.o qr_gs_decomp.o qr_gs_solve.o mv_mult.o matrix.o vector.o matrix.h vector.h

lsfit.o inverse.o qr_gs_inverse.o qr_gs_decomp.o qr_gs_solve.o mv_mult.o matrix.o: matrix.h
lsfit.o qr_gs_solve.o mv_mult.o vector.o: vector.h

.PHONEY:clean
clean:
	$(RM) main *.o *.txt
