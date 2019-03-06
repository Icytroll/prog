CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

default: rosen.svg fit.svg

rosen.svg: rosen.gpi rosen.txt 
	gnuplot $< > $@

fit.svg: fit.gpi fit.txt input.data
	gnuplot $< > $@

rosen.txt fit.txt: main Makefile 
	./$< -4 3 $$(wc --lines input.data) < input.data 1> rosen.txt 2> fit.txt

main: main.o rosenbrock.o LS_fit.o

.PHONEY:clean
clean:
	$(RM) main *.o *.txt *.svg
