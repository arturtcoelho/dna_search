FLAGS=-O3 -Wall -fopenmp

CC=gcc

RM=rm -f

EXEC=dna3

all: $(EXEC)

$(EXEC): dna3.c
	$(CC) $(FLAGS) dna3.c -c -o dna3.o
	$(CC) $(FLAGS) dna3.o -o $(EXEC)

run:
	./$(EXEC)

clean:
	$(RM) dna3.o $(EXEC)