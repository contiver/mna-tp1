CC = gcc
FLAGS = -pedantic -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes -O3 -lm -std=c11
# add this later!
# -Werror

all: default clean

default: main.o matrix.o ccsMatrix.o utils.o
	$(CC) $(FLAGS) -o A_matrix main.o matrix.o ccsMatrix.o utils.o

main.o: main.c matrix.o
	$(CC) $(FLAGS) -c main.c

matrix.o: matrix.c matrix.h
	$(CC) $(FLAGS) -c matrix.c

ccsMatrix.o: ccsMatrix.c matrix.h
	$(CC) $(FLAGS) -c ccsMatrix.c

utils.o: utils.c utils.h
	$(CC) $(FLAGS) -c utils.c

.PHONY: clean
clean:
	rm -f *.o
