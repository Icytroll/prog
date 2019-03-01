CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

default: logSqrt.txt E.txt E.svg

E.svg: E.gpi Makefile
	gnuplot $< > $@

logSqrt.txt E.txt: main
	./$< 1> logSqrt.txt 2> E.txt
	cat logSqrt.txt

main: main.o logSqrt.o E_numerical.o

.PHONEY:clean
clean:
	$(RM) main *.o *.txt *.svg
