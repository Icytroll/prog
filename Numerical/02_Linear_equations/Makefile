CFLAGS = -Wall -std=gnu1x
LDLIBS = -lm

defaiult: A.txt
	cat $<

A.txt: main
	./$<

main: main.o A_decomp.o A_solve.o qr_gs_decomp.o qr_gs_solve.o vector.o vector.h matrix.o matrix.h

main.o A_solve.o vector.o qr_gs_solve: vector.h
main.o A_solve.o A_decomp.o matrix.o qr_gs_decomp qr_gs_solve: matrix.h

.PHONEY:clean
clean:
	$(RM) main *.o *.txt
