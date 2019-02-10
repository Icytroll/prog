CFLAGS = -Wall -std=gnu99 -O1
LDLIBS = -lm
CC = gcc

main: main.o equal.o
	gcc -o main main.o equal.o

main.o: main.c
	gcc -c main.c

equal.o: equal.c
	gcc -c equal.c

.PHONEY:clean
clean:
	$(RM) *.o
	find -type f -executable -delete
