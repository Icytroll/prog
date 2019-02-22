CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY:all
all: airy.svg out.txt data.txt

out.txt data.txt: main
	./$< 1> out.txt

main: main.o

airy.svg: airy.gpi data.txt
	gnuplot $< > $@

.PHONEY:clean
clean:
	$(RM) *.o *.svg *.txt main
