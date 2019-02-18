CFLAGS = -Wall -std=gnu99

default: out.txt 

main: main.o

out.txt: main
	./$< > $@

.PHONEY:clean
clean:
	$(RM) main main.o out.txt
