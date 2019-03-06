CFLAGS = -Wall -std=gnu11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

default: rosen.txt log.txt rosen.svg s_wave.svg

rosen.svg: rosen.txt rosen.gpi
	gnuplot rosen.gpi > $@

s_wave.svg: s_wave.gpi sData.txt
	gnuplot $< > $@

log.txt rosen.txt sData.txt: main
	./$< 2> log.txt 1> rosen.txt

main: main.o rosen_extremum.o s_wave_root.o fe.o

.PHONEY:clean
clean:
	$(RM) main *.o *.txt *.svg
