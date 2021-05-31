FLAGS=-O3 -Wall -fopenmp

CC=gcc

RM=rm -f

EXEC=dna4

all: $(EXEC)

$(EXEC): dna4.c
	$(CC) $(FLAGS) dna4.c -c -o dna4.o
	$(CC) $(FLAGS) dna4.o -o $(EXEC)

run:
	./$(EXEC)

clean:
	$(RM) dna4.o $(EXEC)