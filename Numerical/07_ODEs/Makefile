HDRDIR = ./hdr
CC     = gcc
CFLAGS = -Wall -std=gnu1x -I$(HDRDIR)
LDLIBS = -lm

SRCDIR = ./src
OBJDIR = ./obj

SRCS := $(shell find $(SRCDIR) -name *.c)
OBJS := $(SRCS:%:=$(OBJDIR)/%.o)
HDRS := $(shell find $(HDRDIR) -name *.h)

default: output/log.txt output/iterations.txt output/fit.svg
	
output/log.txt output/iterations.txt output/fit.svg: log.txt iterations.txt fit.svg
	mv log.txt output/log.txt
	mv iterations.txt output/iterations.txt
	mv fit.svg output/fit.svg

fit.svg: $(SRCDIR)/fit.gpi fit_data.txt
	gnuplot $< > $@
	mv fit_data.txt output/fit_data.txt

log.txt iterations.txt fit_data.txt: main
	./$< 1> log.txt 2> iterations.txt

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
	$(RM) main $(OJBDIR)/*.o output/*.txt output/*.svg
