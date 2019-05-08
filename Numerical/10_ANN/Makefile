HDRDIR = ./hdr
CC     = gcc
CFLAGS = -Wall -std=gnu1x -I$(HDRDIR)
LDLIBS = -lm

SRCDIR = ./src
OBJDIR = ./obj

SRCS := $(shell find $(SRCDIR) -name *.c)
OBJS := $(SRCS:%:=$(OBJDIR)/%.o)
HDRS := $(shell find $(HDRDIR) -name *.h)

default: output/1Dplot.svg output/log.txt

output/log.txt: log.txt
	mv log.txt output/log.txt

output/1Dplot.svg: 1Dplot.svg
	mv 1Dplot.svg output/1Dplot.svg
	mv 1Ddata.txt output/1Ddata.txt

1Dplot.svg: src/1Dplot.gpi 1Ddata.txt
	gnuplot $< > $@

log.txt 1Ddata.txt: main
	./$< 1> log.txt 2> 1Ddata.txt

main: $(OBJS)
	$(CC) $(CFLAGS) $(HDRS) -o main $(OBJS) $(LDLIBS)

$(OJBDIR)/%.o: $(SRCDIR)/%.c $(HDRS)
	$(CC) $(CFLAGS) $(HDRS) -c $< -o $@

show_files:
	@echo 'Source files: ' $(SRCS)
	@echo 'Object files: ' $(OBJS)
	@echo 'Header files: ' $(HDRS)

.PHONEY:clean
clean:
	$(RM) main $(OJBDIR)/*.o output/*.txt output/*.svg output/*.data
