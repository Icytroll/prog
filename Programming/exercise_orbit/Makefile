CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY:default clean
default: logistic.data logistic.svg equatorial.data equatorial.svg

equatorial.svg: equatorial.gpi equatorial.data Makefile
	gnuplot $< > $@

logistic.svg: logistic.data Makefile
	gnuplot -e "\
	set terminal svg background 'white'; \
	set xlabel 'x';\
	set ylabel 'y';\
	plot '$<' using 1:2 with lines title 'mylogistic' \
	, '$<' using 1:3 with lines linewidth 5 dashtype (10,20) title 'Analytic' \
	" > $@

logistic.data equatorial.data: main
	./$< 1> logistic.data 2> equatorial.data

main: main.o mylogistic.o myequatorial.o

clean:
	$(RM) main *.o *.data *.svg
