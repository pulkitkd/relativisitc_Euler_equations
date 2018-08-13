# Set the C compiler
CC            = gcc -std=c11 -O2 -g \
		-Wall \
		-Wextra \
		-Wformat-nonliteral \
		-Wcast-align \
		-Wpointer-arith \
		-Wbad-function-cast  \
		-Winline -Wundef \
		-Wno-unused-parameter \
		-Wfloat-equal \
		-pedantic \
		-Wcast-qual\
		-Wshadow \
		-Wwrite-strings
# 		-Wmissing-declarations
# 		-Wnested-externs  \
# 		-Wmissing-prototypes  -Wstrict-prototypes


HDRS	        = dg.h dg1d.h

LD            = $(CC)

MAKEFILE      = makefile

OBJS          = main.o adaptation.o atgauss.o cerk.o div.o error.o flux.o gauss.o \
		init.o meshvel.o movegrid.o predictor.o project.o result.o   \
		shape.o taylor.o  time.o update.o utils.o con2prim.o

PROGRAM       = dg

all:            $(PROGRAM)

.c.o: $(HDRS)
		$(CC) -c $*.c

$(PROGRAM):     $(OBJS)
		$(LD) -o $(PROGRAM) $(OBJS) -lm

clean:;         rm -f $(OBJS) core $(PROGRAM)
