CFLAGS = -Wall -std=gnu99
LDLIBS = -lm

.PHONEY: default clean
default: out.txt
	cat $<

out.txt: main
	./$< > $@

main: main.o nvector.o double_equal.o

main.o nvector.o: nvector.h

clean:
	$(RM) main main.o nvector.o double_equal.o out.txt 
