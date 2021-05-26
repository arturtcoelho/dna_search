#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h> 

#include <errno.h>
#include <stdarg.h>

// MAX char table (ASCII)
#define MAX 256
#define QUERY_LINE_SIZE 1000
#define DESC_LINE_SIZE 100

#define MAX_QUERY 200000
#define MAX_DNA 100000
#define QUERY_SIZE 100000
#define DNA_SIZE 100000
#define MAX_SECTORS 100

FILE *fout;
char *dna, *query;
long dna_s, query_s;
char **query_map, **dna_map;
long query_map_s, dna_map_s;
char **dna_remap;
int num_sectors;

void alloc_string_map();
void check (int test, const char * message, ...);
void openfiles();
void map_files();
void map_strings();
int bmhs(char *string, int n, char *substr, int m);
void closefiles();
void remove_eol(char *line, int len);
void free_all();

int main()
{
    map_files();
    alloc_string_map();
    map_strings();
	openfiles();

	int found = 0;
	int result;

	for (long i = 0; i < query_map_s; i++) {
		fprintf(fout, ">Query string #%ld\n", i);
		
		// read database and search
		found = 0;
		for (long j = 0; j < num_sectors; j++){

			result = bmhs(dna_remap[j], strlen(dna_remap[j]), query_map[i], strlen(query_map[i]));
			if (result > 0) {
				fprintf(fout, "> Escherichia coli K-12 MG1655 section %ld of 400 of the complete genome\n%d\n", j+1, result);
				found++;
			}
		}

		if (!found)
			fprintf(fout, "NOT FOUND\n");
	}

	closefiles();
	free_all();

	return EXIT_SUCCESS;
}

void alloc_string_map()
{
    dna_map = malloc(MAX_DNA*sizeof(char*));
    query_map = malloc(MAX_QUERY*sizeof(char*));

	if (!dna_map || !query_map) exit(EXIT_FAILURE);

    for (int i = 0; i < MAX_DNA; i++) {
        dna_map[i] = malloc(DNA_SIZE);
		if (!dna_map[i]) exit(EXIT_FAILURE);
    }
    for (int i = 0; i < MAX_QUERY; i++) {
        query_map[i] = malloc(QUERY_SIZE);
		if (!query_map[i]) exit(EXIT_FAILURE);
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
	if (!dna_remap) exit(EXIT_FAILURE);
    for (int i = 0; i < num_sectors; i++){
        dna_remap[i] = malloc(DNA_SIZE * sizeof(char));
		if (!dna_remap[i]) exit(EXIT_FAILURE);
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
            if (*(query+last_pos) == '>' || *(query+last_pos+1) == '>') {
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
    // File descriptor.
    int fd;
    // Information about the file opened
    struct stat s;
    int status;

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

	close (fd);

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

	close (fd);
}

// Boyers-Moore-Hospool-Sunday algorithm for string matching
int bmhs(char *string, int n, char *substr, int m) {

	int d[MAX];
	int i, j, k;

	// pre-processing
	for (j = 0; j < MAX; j++)
		d[j] = m + 1;
	for (j = 0; j < m; j++)
		d[(int) substr[j]] = m - j;

	// searching
	i = m - 1;
	while (i < n) {
		k = i;
		j = m - 1;
		while ((j >= 0) && (string[k] == substr[j])) {
			j--;
			k--;
		}
		if (j < 0)
			return k + 1;
		i = i + d[(int) string[i + 1]];
	}

	return -1;
}

void remove_eol(char *line, int len) 
{
	int i = len - 1;
	while (line[i] == '\n' || line[i] == '\r') {
		line[i] = 0;
		i--;
	}
}

void openfiles() 
{
	fout = fopen("dna.out", "w");
	if (fout == NULL) {
		perror("fout");
		exit(EXIT_FAILURE);
	}
}


void closefiles() 
{
	fflush(fout);
	fclose(fout);
}

void free_all()
{
	munmap(dna, dna_s);
	munmap(query, query_s);
	for (long i = 0; i < query_map_s; i++) free(query_map[i]);
	for (long i = 0; i < dna_map_s; i++) free(dna_map[i]);
	free(query_map);
	free(dna_map);
	for (long i = 0; i < num_sectors; i++) free(dna_remap[i]);
	free(dna_remap);	
}