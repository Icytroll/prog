CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

default: Report.pdf

Report.pdf: Report.tex plot.pdf
	pdflatex $<
	pdflatex $<

plot.pdf: output.txt Makefile
	gnuplot -e "\
	set term pdf background 'white'; \
	set xlabel 'x'; \
	set ylabel 'a'; \
	set grid xtics; \
	set grid ytics; \
	set key top left; \
	plot '$<' u 1:2 w l title 'myarctan' \
	, '$<' u 1:3 w l lw 3 dashtype (10,20) title 'atan' \
	" > $@

output.txt: main
	./$< > $@

main: main.o myarctan.o

.PHONEY:clean
clean:
	$(RM) main *.o
