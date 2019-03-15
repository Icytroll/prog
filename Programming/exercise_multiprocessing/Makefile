CFLAGS = -fopenmp -Wall -std=gnu11
LDLIBS = -lm

default: Pi.svg

Pi.svg: PiData.txt Makefile
	gnuplot -e "\
	set term svg background 'white'; \
	set xlabel 'N'; \
	set log x;\
	set log y;\
	a=1;\
	f(x)=a/sqrt(x);\
	pi=4*atan(1);\
	fit f(x) '$<' index 0 using 1:(abs((\$$2)-pi)) via a;\
	plot '$<' index 0 using 1:(abs((\$$2)-pi)) with lines title 'Pthreads' \
	, '$<' index 1 using 1:(abs((\$$2)-pi))  with lines title 'OpenMP' \
	, f(x) with lines\
	" > $@

PiData.txt: main
	./$< 12 30 > $@

.PHONEY:clean
clean:
	$(RM) main *.o *.txt *.svg
