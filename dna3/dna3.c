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

#include <omp.h>

// MAX char table (ASCII)
#define MAX 256
#define QUERY_LINE_SIZE 1000
#define DESC_LINE_SIZE 200

#define MAX_QUERY 200000
#define MAX_DNA 100000
#define QUERY_SIZE 100000
#define DNA_SIZE 100000
#define MAX_SECTORS 100

#define OUT_MAP_MULTIPLIER 2.1

// File maps and lengths
// these files are maped using mmap, and can be accessed as strings
// The lengths are obtained with the file descriptor
char *dna, *query, *fout;
long dna_s, query_s, fout_s;

// String maps to read, and their lengths
// they are strings of strings, each line representing a file line
// Dna remap is used to join the dna map lines by sector
char **query_map, **dna_map, **dna_remap;
long query_map_s, dna_map_s, num_sectors;

// 
char **output_map;
long out_map_s;

void alloc_string_map();
void check (int test, const char * message, ...);
void openfiles();
void map_files();
void map_strings();
void map_output();
int bmhs(char *string, int n, char *substr, int m);
void closefiles();
void remove_eol(char *line, int len);
void free_all();

int main()
{
	// clock_t begin_0 = clock();

    map_files();
    alloc_string_map();
    map_strings();

	// clock_t end = clock();
	// printf("time reading and allocating memory = %f\n", (double)(end - begin) / CLOCKS_PER_SEC);

	int num_threads;
    fout_s = 0;
    
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		if (!id) num_threads = omp_get_num_threads();

		// to know number of threads
		#pragma omp barrier

		long chunk = query_map_s / num_threads;
		long init_for = id*chunk;
		long end_for = (id+1)*chunk;
		long line_number = init_for*OUT_MAP_MULTIPLIER;
        long total_writen = 0;

		// for each query in block
		for (long i = init_for; i < end_for; i++) {
            total_writen += sprintf(output_map[line_number++], ">Query string #%ld\n", i);
			int found = 0;

			// for each database sector
			for (long j = 0; j < num_sectors; j++){
				
				// actualy search
				int result = bmhs(dna_remap[j], strlen(dna_remap[j]), query_map[i], strlen(query_map[i]));

				// if found store the value on output array
				if (result > 0) {
                    total_writen += sprintf(output_map[line_number], "> Escherichia coli K-12 MG1655 section %ld of 400 of the complete genome\n%d\n", j+1, result);
                    line_number+=2;
					found++;
				}
			}
            if (!found) {
                total_writen += sprintf(output_map[line_number++], "NOT FOUND\n");
            }
		}

        #pragma omp critical
        {
            fout_s += total_writen;
        }

	}

	map_output();

	for (long i = 0, len = 0; i < out_map_s-1; i++){
        if (output_map[i][0]) {
            len += sprintf(fout+len, "%s", output_map[i]);
            // printf("%s", output_map[i]);
        }
	}

	free_all();

	// clock_t end = clock();
	// printf("Total internal time = %f\n", (double)(end - begin_0) / CLOCKS_PER_SEC);
	// printf("Total time = %f\n", ((double)(end - begin_0) / CLOCKS_PER_SEC)/num_threads);

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
	out_map_s = (query_map_s)*OUT_MAP_MULTIPLIER; // heuristic of how many found objects

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
    if (fout == MAP_FAILED) {
        close(fd);
        perror("Error mmapping the file");
        exit(EXIT_FAILURE);
    }
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