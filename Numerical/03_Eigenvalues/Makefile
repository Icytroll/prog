CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

default: cyclic_sweep_time.svg out.txt

cyclic_sweep_time.svg: cyclic_sweep_time.gpi cyclic_sweep_time.txt
	gnuplot $< > $@

out.txt cyclic_sweep_time.txt: main
	./$< 1> out.txt

main: main.o jacobi_cyclic.o jacobi_single_row.o matrix_print.o

.PHONEY:clean
clean:
	$(RM) main *.o *.txt *.svg
