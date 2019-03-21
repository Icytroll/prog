CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

default: linear_spline.svg

linear_spline.svg: linear.gpi data.txt
	gnuplot $< > $@

data.txt: main
	./$<

main: main.o linterp.o linterp_integ.o binary_search.o qspline.o

main.o qspline.o: qspline.h

.PHONEY:clean
clean:
	$(RM) main *.o *.txt
