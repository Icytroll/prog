CFLAGS = -Wall -std=gnu99 $$(gls-config --clfags)
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
