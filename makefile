CC = gcc
FLAGS = -Wall -O3 -lm

default: main.o matrix.o
	$(CC) $(FLAGS) -o A_matrix main.o matrix.o

main.o: main.c matrix.o
	$(CC) $(FLAGS) -c main.c

matrix.o: matrix.c matrix.h
	$(CC) $(FLAGS) -c matrix.c

.PHONY: clean
clean:
	rm -f *.o
