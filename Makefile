CC = g++
SRCS = main.cpp linearalgebra.cpp keyFunction.cpp inout.cpp 
OBJS = $(SRCS:.c=.o)
all:
	$(CC) -o $(SRCS)
.PHONY : clean
clean :
	rm *.o main
