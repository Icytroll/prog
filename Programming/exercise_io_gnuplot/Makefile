CFLAGS = -Wall -std=gnu99
LDLIBS = -lm

.PHONEY: test.io test.cmdclean
test.io          : test.io.out.txt   ; cat $<
test.io.out.txt  : main-stdio        ; ./$< > $@
main-stdio       : main-stdio.o	     # built-in linking rule is used here

test.cmd         : test.cmd.out.txt  ; cat $<
test.cmd.out.txt : main-cmdline      ; ./$< > $@
main-cmdline     : main-cmdline.o

clean            :                   ; $(RM) main-stdio main-stdio.o test.io.out.txt main-cmdline main-cmdline.o test.cmd.out.txt

