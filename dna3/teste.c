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

char *fout;
long fout_s;

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
    fout_s = query_s/30; // Heuristic

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

void map_output()
{
    int i;
    int fd;
    int result;
    char *map;  /* mmapped array of int's */
    char *filename = "dna.out";  /* mmapped array of int's */

    /* Open a file for writing.
     *  - Creating the file if it doesn't exist.
     *  - Truncating it to 0 size if it already exists. (not really needed)
     *
     * Note: "O_WRONLY" mode is not sufficient when mmaping.
     */
    fd = open(filename, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
    if (fd == -1) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }

    /* Stretch the file size to the size of the (mmapped) array of ints
     */
    result = lseek(fd, fout_s-1, SEEK_SET);
    if (result == -1) {
        close(fd);
        perror("Error calling lseek() to 'stretch' the file");
        exit(EXIT_FAILURE);
    }
    
    /* Something needs to be written at the end of the file to
     * have the file actually have the new size.
     * Just writing an empty string at the current file position will do.
     *
     * Note:
     *  - The current position in the file is at the end of the stretched 
     *    file due to the call to lseek().
     *  - An empty string is actually a single '\0' character, so a zero-byte
     *    will be written at the last byte of the file.
     */
    result = write(fd, "", 1);
        if (result != 1) {
        close(fd);
        perror("Error writing last byte of the file");
        exit(EXIT_FAILURE);
    }

    /* Now the file is ready to be mmapped.
     */
    fout = mmap(0, fout_s, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (map == MAP_FAILED) {
        close(fd);
        perror("Error mmapping the file");
        exit(EXIT_FAILURE);
    }
}

int main ()
{
    map_files();
    alloc_string_map();
    map_strings();
    map_output();

    for (int i = 0; i < 30; i++)
    {
        fout[i] = i+97;
    }
    fout[11] = '\n';
    

    return 0;
}