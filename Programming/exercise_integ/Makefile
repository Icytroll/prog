CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

default: logSqrt.txt E.txt E.svg

E.svg: E.txt Makefile
	gnuplot -e "\
	set terminal svg background 'white' enhanced; \
	set yrange [0:5]; \
	set xlabel '{/Symbol a}'; \
	set ylabel 'E({/Symbol a})'; \
	plot \
		0.5 title 'E({/Symbol a}) = 0.5' \
		,'$<' using 1:2 with lines title 'Expectation value E({/Symbold a})' \
	" > $@

logSqrt.txt E.txt: main
	./$< 1> logSqrt.txt 2> E.txt

main: main.o logSqrt.o expecValue.o

.PHONEY:clean
clean:
	$(RM) main *.o *.txt *.svg
