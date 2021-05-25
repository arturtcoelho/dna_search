FLAGS=-O3 -Wall

CC=gcc

RM=rm -f

EXEC=dna

all: $(EXEC)

$(EXEC): dna.c
	$(CC) $(FLAGS) dna.c -c -o dna.o
	$(CC) $(FLAGS) dna.o -o $(EXEC)

run:
	./$(EXEC)

clean:
	$(RM) dna.o $(EXEC)