CFLAGS = -fopenmp -Wall -std=gnu11
LDLIBS = -lm

default: Pi.svg

Pi.svg: PiData.txt Makefile
	gnuplot -e "\
	set term svg background 'white'; \
	set xlabel 'N'; \
	plot '$<' index 0 using 1:2 with lines title 'Pthreads' \
	, '$<' index 1 using 1:2 with lines title 'OpenMP' \
	, 3.14159265359 lc 7 lw 2 title 'Pi'\
	" > $@

PiData.txt: main
	./$< > $@

.PHONEY:clean
clean:
	$(RM) main *.o *.txt *.svg
