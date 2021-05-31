#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Mapping and system header functions
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h> 

// Error traceback
#include <errno.h>

// openMP lib
#include <omp.h>

// MAX char table (ASCII)
#define MAX 256

// File line max size
#define QUERY_LINE_SIZE 1000
#define DESC_LINE_SIZE 200

// Max file structure size
#define MAX_QUERY 200000
#define MAX_DNA 100000
#define QUERY_SIZE 100000
#define DNA_SIZE 100000
#define MAX_SECTORS 100

// based on the total found strings to be writen
#define OUT_MAP_MULTIPLIER 2.1

// File maps and lengths
// these files are maped using mmap, and can be accessed as strings
// The lengths are obtained with the file descriptor
#define OUTPUT_FILE_NAME "dna.out"
#define DNA_FILE_NAME "dna.in"
#define QUERY_FILE_NAME "query.in"
char *dna, *query, *fout;
long dna_s, query_s, fout_s;

// String maps to read, and their lengths
// they are strings of strings, each line representing a file line
// Dna remap is used to join the dna map lines by sector
char **query_map, **dna_map, **dna_remap;
long query_map_s, dna_map_s, num_sectors;

// this array of strings is where each thread write 
// the correspondent line to be put on the final file
char **output_map;
long out_map_s;

// Open and map read files
void map_files();

// Aux error checking function 
void check(int, char *);

// Allocate memory to read data
void alloc_string_map();

// Map file data to program data
void map_strings();

// Aux remove end of line char
void remove_eol(char *, int);

// The main seach algorithm
int bmhs(char *, int, char *, int);

// Open and map the output file
void map_output();

// Free memory and close files
void free_all();


// Entry point
int main()
{
	clock_t begin = clock();

    // Open files and allocate memory
    map_files();
    alloc_string_map();
    map_strings();

    fout_s = 0; // output file size

    int num_threads = omp_get_max_threads();
    
    // to know where to write independently
    long *line_number = calloc(num_threads, 64);
    if (!line_number) exit(EXIT_FAILURE);
    
    long total_writen = 0;
    
    #pragma omp parallel reduction (+:total_writen) 
    {
        int init = 0;
        #pragma omp for schedule(static)
        for (long i = 0; i < query_map_s; i++) { // for each query in block

            int id = omp_get_thread_num();
            int found = 0;

            if (!init) {
                line_number[id] = i * OUT_MAP_MULTIPLIER;
                init++;
            }

            total_writen += sprintf(output_map[line_number[id]++], 
                                    ">Query string #%ld\n", i);

            for (long j = 0; j < num_sectors; j++){ // for each database sector
                
                // actualy search with bmhs
                int result = bmhs(dna_remap[j], strlen(dna_remap[j]), 
                                    query_map[i], strlen(query_map[i]));

                // if found store the value on output array
                if (result > 0) {
                    #pragma omp critical
                    {
                        total_writen += sprintf(output_map[line_number[id]], 
                                                "> Escherichia coli K-12 MG1655 section %ld of 400 of the complete genome\n%d\n"
                                                , j+1, result);
                    }
                    line_number[id]+=2;
                    found++;
                }
            }
            
            if (!found) {
                #pragma omp critical
                {
                    total_writen += sprintf(output_map[line_number[id]++], 
                                            "NOT FOUND\n");
                }
            }
        }
    }

    fout_s = total_writen;

    // Create and mat output file, now with know size
	map_output();

    // white to file
	for (long i = 0, len = 0; i < out_map_s-1; i++){
        if (output_map[i][0]) {
            len += sprintf(fout+len, "%s", output_map[i]);
        }
	}

	free_all();

	clock_t end = clock();
	printf("%f\n", (double)(end - begin) / CLOCKS_PER_SEC);

	return 0;
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

    // for each dna line, copy it
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

    // join dna lines    
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

    // copy query lines
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
	out_map_s = (query_map_s)*OUT_MAP_MULTIPLIER; // heuristic of how many found objects

    // get a output map, knowing the input size
	output_map = malloc(out_map_s * sizeof(char*));
	if (!output_map) exit(EXIT_FAILURE);
	for (long i = 0; i < out_map_s-1; i++) {
		output_map[i] = malloc(DESC_LINE_SIZE * sizeof(char));
		if (!output_map[i]) exit(EXIT_FAILURE);
        output_map[i][0] = '\0';
	}
	
}

void map_output()
{
    int fd;
    int result;

    // Open file descriptor
    fd = open(OUTPUT_FILE_NAME, O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
    check(fd == -1, "Error opening file for writing");

    // Go to end of file with set distance
    result = lseek(fd, fout_s-1, SEEK_SET);
    check(result == -1, "Error calling lseek() to 'stretch' the file");
    
    // write byte to end of file to make sure we have space
    result = write(fd, "", 1);  
    check(result != 1, "Error writing last byte of the file");

    // map entire virtual space to string
    fout = mmap(0, fout_s, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    check(fout == MAP_FAILED, "Error mmapping the file");
}

void check(int test, char * message)
{
    if (test) {
        fprintf (stderr, "%s\n", message);
        exit (EXIT_FAILURE);
    }
}

void map_files()
{
    int fd;
    struct stat s;
    int status;

    // Dna file
    // Open file descriptor
    fd = open (DNA_FILE_NAME, O_RDONLY);
    check(fd < 0, "open dna failed");

    // Get system header
    status = fstat (fd, &s);
    check(status < 0, "stat dna failed");
    dna_s = s.st_size;

    // Map virtual memory to string
    dna = mmap (0, dna_s, PROT_READ, MAP_PRIVATE, fd, 0);
    check(dna == MAP_FAILED, "mmap dna failed");

	close (fd);

    // Query file
    fd = open (QUERY_FILE_NAME, O_RDONLY);
    check(fd < 0, "open query failed");

    status = fstat (fd, & s);
    check(status < 0, "stat query failed");
    query_s = s.st_size;

    query = mmap (0, query_s, PROT_READ, MAP_PRIVATE, fd, 0);
    check(query == MAP_FAILED, "mmap query failed");

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

void free_all()
{
	munmap(dna, dna_s);
	munmap(query, query_s);
	munmap(fout, fout_s);
	for (long i = 0; i < query_map_s; i++) free(query_map[i]);
	for (long i = 0; i < dna_map_s; i++) free(dna_map[i]);
	for (long i = 0; i < num_sectors; i++) free(dna_remap[i]);
	free(query_map);
	free(dna_map);
	free(dna_remap);

	free(output_map);
}