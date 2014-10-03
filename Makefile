#CC = musl-gcc
CC = gcc
FLAGS = -pedantic -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes -lm -std=c11 -O3 -s
# add this later!
# -Werror

all: dense ccsA hardcodedCCSA clean

dense: denseA.o matrix.o utils.o ccsMatrix.o
	$(CC) $(FLAGS) -o denseA denseA.o matrix.o utils.o

denseA.o: denseA.c matrix.h
	$(CC) $(FLAGS) -c denseA.c

ccsA: ccsA.o utils.o ccsMatrix.o
	$(CC) $(FLAGS) -o ccsA ccsA.o ccsMatrix.o utils.o

ccsA.o: ccsA.c matrix.h
	$(CC) $(FLAGS) -c ccsA.c

hardcodedCCSA: hardcodedCCSA.o utils.o ccsMatrix.o
	$(CC) $(FLAGS) -o hardcodedCCSA hardcodedCCSA.o ccsMatrix.o utils.o

hardcodedCCSA.o: hardcodedCCSA.c
	$(CC) $(FLAGS) -c hardcodedCCSA.c

matrix.o: matrix.c matrix.h
	$(CC) $(FLAGS) -c matrix.c

ccsMatrix.o: ccsMatrix.c matrix.h
	$(CC) $(FLAGS) -c ccsMatrix.c

utils.o: utils.c utils.h
	$(CC) $(FLAGS) -c utils.c

.PHONY: clean
clean:
	rm -f *.o
