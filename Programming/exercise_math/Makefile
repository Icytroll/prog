CFLAGS = -Wall -std=gnu99 -O1
LDLIBS = -lm

default: out.txt

out.txt: main
	./$< > $@

main: main.o

.PHONEY:clean
clean:
	$(RM) main main.o out.txt
