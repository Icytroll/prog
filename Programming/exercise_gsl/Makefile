CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = `gsl-config --libs`

all: airy.svg data.txt out.txt

main: main.o

airy.svg: airy.gpi data.txt
	gnuplot $< > $@

out.txt data.txt: main
	./$< > $@

.PHONEY:clean
clean:
	$(RM) *.o *.svg *.txt main
