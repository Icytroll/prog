CFLAGS = -Wall -std=gnu1x
LDLIBS = -lm

default: Output/A.txt

Output/A.txt: A.txt
	mv A.txt Output/A.txt

A.txt: main
	./$<

main: main.o matrix.o vector.h matrix.o matrix.h

.PHONEY:clean
clean:
	$(RM) main *.o Output/*.txt
