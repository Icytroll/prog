CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY:default clean
default: logistic.data logistic.svg equatorial.data equatorial.svg

equatorial.svg: equatorial.data Makefile
	gnuplot -e "\
	set terminal svg background 'white'; \
	plot $< using (1/$2)*(sin($1)):(1/$2)*(cos($1)) with lines 'Normal circle'\
	" > $@

logistic.svg: logistic.data Makefile
	gnuplot -e "\
	set terminal svg background 'white'; \
	plot '$<' using 1:2 with lines title 'mylogistic' \
	" > $@

logistic.data equatorial.data: main
	./$< 1> logistic.data 2> equatorial.data

main: main.o mylogistic.o myequatorial.o

clean:
	$(RM) main *.o *.data *.svg
