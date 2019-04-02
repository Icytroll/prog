CFLAGS = -Wall -std=gnu1x 
LDLIBS = -lm

default: Output/fit.svg Output/out.txt

Output/fit.svg Output/out.txt: fit.svg out.txt
	mv fit.svg Output/fit.svg
	mv out.txt Output/out.txt

fit.svg: fit.gpi fit.txt
	gnuplot $< > $@

out.txt fit.txt: main
	./$< > out.txt

main: main.o lsfit.o inverse.o qr_gs_inverse.o qr_gs_decomp.o qr_gs_solve.o mv_mult.o matrix.o vector.o matrix.h vector.h

lsfit.o inverse.o qr_gs_inverse.o qr_gs_decomp.o qr_gs_solve.o mv_mult.o matrix.o: matrix.h
lsfit.o qr_gs_solve.o mv_mult.o vector.o: vector.h

.PHONEY:clean
clean:
	$(RM) main *.o *.txt Output/*.svg
