CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

default: Report.pdf

Report.pdf: Report.tex atan.pdf tan.pdf
	pdflatex $<
	pdflatex $<

atan.pdf: output.txt Makefile
	gnuplot -e "\
	set term pdf background 'white'; \
	set xlabel 'x'; \
	set ylabel 'a'; \
	set grid xtics; \
	set grid ytics; \
	set key top left; \
	plot '$<' u 1:2 w l lw 2 title 'myarctan' \
	, '$<' u 1:3 w l lw 4 dashtype (20,30) title 'atan' \
	" > $@

tan.pdf: Makefile
	gnuplot -e "\
	set term pdf background 'white'; \
	set samples 1000; \
	set size ratio -1; \
	set xlabel 'a'; \
	set ylabel 'x'; \
	set grid xtics; \
	set grid ytics; \
	set arrow from 0,-7 to 0,0 nohead lc 7 dashtype (20,12) lw 2; \
	set label 1 at 0,-7.5 '{a_0}' center; \
	set arrow from 0,0 to -7,-7 nohead lc 7 lw 2; \
	set label 2 at -7,-7.5 '{a_1}' center; \
	set arrow from -7,-7 to -7,-0.87 nohead lc 7 dashtype (20,12) lw 2; \
	set arrow from -7,-0.87 to -10,-6.5 nohead lc 7 lw 2; \
	plot [-10:5] [-10:5] tan(x) lc 3 title 'tan(a)'\
	, -7 lc 8 lw 2 title 'x = 7'\
	" > $@


output.txt: main
	./$< > $@

main: main.o myarctan.o

.PHONEY:clean
clean:
	$(RM) main *.o *.txt *.pdf
