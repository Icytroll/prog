CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

default: linear_spline.svg quadratic_spline.svg cubic_spline.svg

linear_spline.svg: linear.gpi data.txt
	gnuplot $< > Output/linear_spline.svg

quadratic_spline.svg: quadratic.gpi data.txt
	gnuplot $< > Output/quadratic_spline.svg

cubic_spline.svg: cubic.gpi data.txt
	gnuplot $< > Output/cubic_spline.svg

data.txt: main
	./$<

main: main.o linterp.o linterp_integ.o binary_search.o qspline.o cspline.o

main.o qspline.o: qspline.h

main.o cspline.o: cspline.h

.PHONEY:clean
clean:
	$(RM) main *.o *.txt *.svg
