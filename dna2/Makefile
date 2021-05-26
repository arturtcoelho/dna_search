FLAGS=-O3 -Wall -fopenmp

CC=gcc

RM=rm -f

EXEC=dna2

all: $(EXEC)

$(EXEC): dna2.c
	$(CC) $(FLAGS) dna2.c -c -o dna2.o
	$(CC) $(FLAGS) dna2.o -o $(EXEC)

run:
	./$(EXEC)

clean:
	$(RM) dna2.o $(EXEC)