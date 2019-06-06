HDRDIR = ./hdr
CC     = gcc
CFLAGS = -Wall -std=gnu1x -I$(HDRDIR) $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

SRCDIR = ./src
OBJDIR = ./obj

SRCS := $(shell find $(SRCDIR) -name *.c)
OBJS := $(SRCS:%:=$(OBJDIR)/%.o)
HDRS := $(shell find $(HDRDIR) -name *.h)

default: output/out.txt output/subdivisions.webm
	
output/out.txt: out.txt
	mv out.txt output/out.txt

output/subdivisions.webm: subdivisions.webm
	mv subdivisions.webm output/subdivisions.webm
	mv subdivisions.txt output/subdivisions.txt

subdivisions.webm: $(SRCDIR)/subdivisions.gpi subdivisions.txt
	gnuplot $<

out.txt subdivisions.txt: main
	./$< 1> out.txt 2> subdivisions.txt

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
	$(RM) main $(OJBDIR)/*.o output/*.txt output/*.svg output/*.webm
