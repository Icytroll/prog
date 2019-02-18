CFLAGS = -Wall -std=gnu99 -O1
LDLIBS = -lm
CC = gcc

default: out.txt

main: main.o equal.o
	gcc -o main main.o equal.o

main.o: main.c
	gcc -c main.c

equal.o: equal.c
	gcc -c equal.c

out.txt: main
	./$< 1 1 0.001 0.002 > $@

.PHONEY:clean
clean:
	$(RM) *.o out*
	find -type f -executable -delete
