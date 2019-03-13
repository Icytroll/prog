CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

default: Report.pdf

Report.pdf: Report.tex plot.pdf
	pdflatex $<
	pdflatex $<

plot.pdf: plot.txt Makefile
	gnuplot -e "\
	set term pdf background 'white'; \
	set xlabel 'x'; \
	set ylabel 'y'; \
	set grid ytics; \
	set grid xtics; \
	plot '$<' u 1:2 w l title 'myerf(x)' \
	, '$<' u 1:3 w l lw 3 dashtype (10,20) title 'gsl erf(x)' \
	" > $@

plot.txt: main
	./$< -3 3 0.01 > $@

main: main.o myerf.o

.PHONEY:clean
clean:
	$(RM) main *.o *.txt *.pdf
