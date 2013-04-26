# Variable definition
CC = g++
LD = g++
CFLAGS =-g #-DDEBUG
#CFLAGS =-DNDEBUG
LDFLAGS = -lm 
PROG = twoSpk
SRCS = modelAdapt.cc acoust.cc acoustDynamic.cc tool.cc gmm.cc filter.cc gammaTone.cc frontEnd.cc twoSpk.cc
OBJS = modelAdapt.o acoust.o acoustDynamic.o tool.o gmm.o filter.o gammaTone.o frontEnd.o twoSpk.o
RM = rm
# end

all:$(PROG) clean

# Link
 $(PROG):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(PROG)

# Compile
%.o:%.cc
	$(CC) -c $(CFLAGS) $<

clean:
	$(RM) $(OBJS)

