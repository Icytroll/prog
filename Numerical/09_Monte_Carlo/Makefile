HDRDIR = ./hdr
CC     = gcc
CFLAGS = -Wall -std=gnu1x -I$(HDRDIR)
LDLIBS = -lm

SRCDIR = ./src
OBJDIR = ./obj

SRCS := $(shell find $(SRCDIR) -name *.c)
OBJS := $(SRCS:%:=$(OBJDIR)/%.o)
HDRS := $(shell find $(HDRDIR) -name *.h)

default: output/out.txt output/error.svg
	
output/out.txt: out.txt
	mv out.txt output/out.txt

output/error.svg: error.svg
	mv error.svg output/error.svg
	mv error.data output/error.data

error.svg: src/error.gpi error.data
	gnuplot $< > $@

out.txt error.data: main
	./$< 1> out.txt 2> error.data

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
