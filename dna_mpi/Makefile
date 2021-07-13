FLAGS=-O3 -Wall -fopenmp

CC=mpicc

RM=rm -f

EXEC=dna

all: $(EXEC)

$(EXEC): dna.c
	$(CC) $(FLAGS) dna.c -c -o dna.o
	$(CC) $(FLAGS) dna.o -o $(EXEC)

run:
	mpirun ./$(EXEC)

clean:
	$(RM) *.o

purge: clean
	$(RM) $(EXEC)