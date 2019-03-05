CFLAGS = -Wall -std=gnu99 -O
LDLIBS = -lm

.PHONEY: default clean
default: out.txt
	cat $<

out.txt: main
	./$< > $@

main: main.o komplex.o

main.o komplex.o: komplex.h

clean:
	$(RM) main main.o komplex.o out.txt
