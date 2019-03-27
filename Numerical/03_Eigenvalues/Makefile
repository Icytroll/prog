CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

default: A.svg out.txt

A.svg: A.gpi A_time.txt
	gnuplot $< > $@

out.txt A_time.txt: main
	./$< 1> out.txt

main: main.o jacobi_cyclic.o matrix_print.o

.PHONEY:clean
clean:
	$(RM) main *.o *.txt *.svg
