CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

default: rosen.txt log.txt rosen.svg

rosen.svg: rosen.txt rosen.gpi
	gnuplot rosen.gpi > $@

log.txt rosen.txt: main
	./$< 2> log.txt 1> rosen.txt

main: main.o rosen_extremum.o

.PHONEY:clean
clean:
	$(RM) main *.o *.txt *.svg
