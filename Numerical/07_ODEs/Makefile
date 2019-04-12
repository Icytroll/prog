HDRDIR = ./hdr
CC     = gcc
CFLAGS = -Wall -std=gnu1x -I$(HDRDIR)
LDLIBS = -lm

SRCDIR = ./src
OBJDIR = ./obj

SRCS := $(shell find $(SRCDIR) -name *.c)
OBJS := $(SRCS:%:=$(OBJDIR)/%.o)
HDRS := $(shell find $(HDRDIR) -name *.h)

default: output/ODE.svg output/log.txt
	
output/ODE.svg output/log.txt: ODE.svg log.txt
	mv ODE.svg output/ODE.svg
	mv log.txt output/log.txt

ODE.svg: $(SRCDIR)/ODE.gpi ODE_data.txt
	gnuplot $< > $@
	mv ODE_data.txt output/ODE_data.txt

log.txt ODE_data.txt: main
	./$< 1> log.txt 2> ODE_data.txt

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
