CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

default: plot.svg

plot.svg: output.txt Makefile
	gnuplot -e "\
	set term svg background 'white'; \
	set xlabel 'N'; \
	set grid xtics; \
	set grid ytics; \
	plot '$<' u 1:2 w l title 'myarctan' \
	, '$<' u 1:3 w l lw 3 dashtype (10,20) title 'atan' \
	" > $@

output.txt: main
	./$< > $@

main: main.o myarctan.o

.PHONEY:clean
clean:
	$(RM) main *.o
