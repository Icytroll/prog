CFLAGS = -Wall -std=gnu99 -O1
LDLIBS = -lm
CC = gcc

default: out.txt

main: main.o equal.o

out.txt: main
	./$< 1 1 0.001 0.002 > $@

.PHONEY:clean
clean:
	$(RM) *.o out.txt main
