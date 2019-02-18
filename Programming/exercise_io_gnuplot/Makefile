CFLAGS = -Wall -std=gnu99
LDLIBS = -lm
CC = gcc

default: test.io.out.txt test.cmd.out.txt out.io.txt out.cmd.txt plot.svg
test.io.out.txt: main-stdio
	echo 1 2 3 4 5 | ./$< > $@
	seq 1 0.5 4 | ./$< > $@

test.cmd.out.txt: main-cmdline
	./$< 1 2 3 4 5 6 > $@
	./$< `seq 0 0.5 2` > $@

main-stdio: main-stdio.o

main-cmdline: main-cmdline.o

out.io.txt: main-stdio input.data.txt
	cat input.data.txt | ./main-stdio > $@

out.cmd.txt: main-cmdline input.data.txt
	./main-cmdline `cat input.data.txt` > $@

input.data.txt:
	seq 0 0.2 5 | sed -e 's/,/./g' > $@

plot.svg:
	gnuplot plot.gpi

.PHONEY:clean
clean:
	$(RM) main-stdio main-stdio.o test.io.out.txt main-cmdline main-cmdline.o test.cmd.out.txt ou.io.txt out.cmd.txt input.data.txt plot.svg

