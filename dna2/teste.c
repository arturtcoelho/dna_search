/* For the size of the file. */
#include <sys/stat.h>
/* This contains the mmap calls. */
#include <sys/mman.h> 
/* These are for error printing. */
#include <errno.h>
#include <string.h>
#include <stdarg.h>
/* This is for open and close. */
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
/* For exit. */
#include <stdlib.h>
/* For the final part of the example. */
#include <ctype.h>
#include <time.h>

char *dna, *query;
long dna_s, query_s;
char **query_map, **dna_map;
long query_map_s, dna_map_s;
char **dna_remap;
int num_sectors;

/* "check" checks "test" and prints an error and exits if it is
   true. */
void check (int test, const char * message, ...)
{
    if (test) {
        va_list args;
        va_start (args, message);
        vfprintf (stderr, message, args);
        va_end (args);
        fprintf (stderr, "\n");
        exit (EXIT_FAILURE);
    }
}

void map_files()
{
    /* The file descriptor. */
    int fd;
    /* Information about the file. */
    struct stat s;
    int status;
    /* The file name to open. */
    char * file_name = "dna.in";

    /* Open the file for reading. */
    fd = open (file_name, O_RDONLY);
    check (fd < 0, "open %s failed: %s", file_name, strerror (errno));

    /* Get the size of the file. */
    status = fstat (fd, & s);
    check (status < 0, "stat %s failed: %s", file_name, strerror (errno));
    dna_s = s.st_size;

    /* Memory-map the file. */
    dna = mmap (0, dna_s, PROT_READ, MAP_PRIVATE, fd, 0);
    check (dna == MAP_FAILED, "mmap %s failed: %s", file_name, strerror (errno));

    file_name = "query.in";

    fd = open (file_name, O_RDONLY);
    check (fd < 0, "open %s failed: %s", file_name, strerror (errno));

    /* Get the size of the file. */
    status = fstat (fd, & s);
    check (status < 0, "stat %s failed: %s", file_name, strerror (errno));
    query_s = s.st_size;

    /* Memory-map the file. */
    query = mmap (0, query_s, PROT_READ, MAP_PRIVATE, fd, 0);
    check (query == MAP_FAILED, "mmap %s failed: %s", file_name, strerror (errno));
}

#define MAX_QUERY 200000
#define MAX_DNA 100000
#define QUERY_SIZE 100000
#define DNA_SIZE 100000
#define MAX_SECTORS 100

void remove_eol(char *line, int len) 
{
	int i = len - 1;
	while (line[i] == '\n' || line[i] == '\r') {
		line[i] = 0;
		i--;
	}
}

void alloc_string_map()
{
    dna_map = malloc(MAX_DNA*sizeof(char*));
    query_map = malloc(MAX_QUERY*sizeof(char*));
    for (int i = 0; i < MAX_DNA; i++) {
        dna_map[i] = malloc(DNA_SIZE);
        dna_map[i][0] = 0;
    }
    for (int i = 0; i < MAX_QUERY; i++) {
        query_map[i] = malloc(QUERY_SIZE);
    }

}

void map_strings()
{
	long map_index = 0;
	int last_pos = 0;

    int sectors[MAX_SECTORS];
    num_sectors = 0;

	for (long i = 0; i < dna_s; i++){
		if (dna[i] == '\n') {
            if (*(dna+last_pos) == '>' || *(dna+last_pos+1) == '>') {
                last_pos = i;
                sectors[num_sectors++] = map_index;
                continue;
            }
            int len = i-last_pos;
			memcpy(*(&dna_map[map_index]), dna+last_pos+1, len);
            remove_eol(dna_map[map_index], len);
            map_index++;
            last_pos = i;
		}
	}
    
    dna_map_s = map_index;

    dna_remap = malloc((num_sectors+1) * sizeof(char*));
    for (int i = 0; i < num_sectors; i++){
        dna_remap[i] = malloc(DNA_SIZE * sizeof(char));
    }
    sectors[num_sectors] = dna_map_s;
    
    for (int i = 0; i < num_sectors; i++){
        int k = 0;
        for (int j = sectors[i]; j < sectors[i+1]; j++) {
            int len = strlen(dna_map[j]);   
            memcpy(dna_remap[i]+k, dna_map[j], len);
            k += len;
        }        
    }

    map_index = 0;
	last_pos = 0;

    for (long i = 0; i < query_s; i++){
		if (query[i] == '\n') {
            if (*(query+last_pos) == '>') {
                last_pos = i;
                continue;
            }
            int len = i-last_pos;
			memcpy(query_map[map_index], query+last_pos+1, len);
            remove_eol(query_map[map_index], len);
            map_index++;
			last_pos = i;
		}
	}

    query_map_s = map_index;
}

int main ()
{
    map_files();
    alloc_string_map();
    map_strings();

    return 0;
}