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
#include <mpi.h>
#define MASTER 0

// MAX char table (ASCII)
#define MAX 256

// File line max size
#define QUERY_LINE_SIZE 1024
#define DESC_LINE_SIZE 256

// Max file structure size
#define MAX_QUERY 262144
#define MAX_DNA 16384
#define QUERY_SIZE 1024
#define DNA_SIZE 32768
#define MAX_SECTORS 128

// based on the total found strings to be writen
#define OUT_MAP_MULTIPLIER 2.1

typedef struct {
    long init, end, size;
} message_t;

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
char **query_map, **dna_map, **dna_remap, **query_remap;
long query_map_s, dna_map_s, num_sectors, num_queries;

// this array of strings is where each thread write 
// the correspondent line to be put on the final file
char *output_map;
long out_map_s;

message_t proc_data;

// Open and map read files
void map_files();

// Aux error checking function 
void check(int, char *);

// Allocate memory to read data
void alloc_string_map();

// Map file data to program data
void map_query();
void map_dna();
void map_output_string(int);

// Aux remove end of line char
void remove_eol(char *, int);

// (Master) Send query and dna data to slaves 
void send_data(int);

// (Slaves) Get data from master 
void get_data(int);

// (Slaves) Send results to master
void send_results(int, long);

// (Master) Get results from slaves
void get_results(int, int, long);

// The main seach algorithm
int bmhs(char *, int, char *, int);

// Open and map the output file
void map_output();

// Free memory and close files
void free_all();

// Entry point
int main(int argc, char **argv)
{   
    MPI_Init(&argc, &argv);

    int id, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    alloc_string_map();

    if (id == MASTER) {
    
        // Open files and allocate memory
        map_files();

        // Send data via MPI to other procecesses
        send_data(num_procs);

    } else { // slave

        get_data(id);

    }

    map_query();
    map_dna();

    long query_chunk = num_queries/num_procs;
    proc_data.init = id * query_chunk;
    proc_data.end = (id+1) * query_chunk+1;
    proc_data.size = query_chunk;

    if (id == num_procs-1) {
        proc_data.end = num_queries+1;
    }

    map_output_string(num_procs);

    long total_writen = 0;
    fout_s = 0; // output file size

    for (long i = proc_data.init; i < proc_data.end-1; i++) { // for each query in block
        
        // printf("TESTE BAH%d\n", id);

        int found = 0;

        total_writen += sprintf(output_map+total_writen, 
                                ">Query string #%ld\n", i);

        for (long j = 0; j < num_sectors; j++){ // for each database sector
            // actualy search with bmhs

            int result = bmhs(dna_remap[j], strlen(dna_remap[j]), 
                                query_remap[i], strlen(query_remap[i]));

            // if found store the value on output array
            if (result > 0) {
                total_writen += sprintf(output_map+total_writen, 
                                        "> Escherichia coli K-12 MG1655 section %ld of 400 of the complete genome\n%d\n"
                                        , j+1, result);
                found++;
            }
        }
        
        if (!found) {
            total_writen += sprintf(output_map+total_writen, 
                                    "NOT FOUND\n");
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (id == MASTER) {
        get_results(id, num_procs, total_writen);

        // Create and mat output file, now with know size
        map_output();

        // white to file    
        memcpy(fout, output_map, fout_s);

    } else {
        send_results(id, total_writen);
    }

    free_all();

    MPI_Finalize();

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

void map_dna()
{
    // Dna

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

}

void map_query()
{
	// Queries

	long map_index = 0;
	int last_pos = 0;

	int queries_disp[MAX_QUERY];
    num_queries = 0;

    for (long i = 0; i < query_s; i++){
		if (query[i] == '\n') {
            if (*(query+last_pos) == '>' || *(query+last_pos+1) == '>') {
                last_pos = i;
				queries_disp[num_queries++] = map_index;
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

    query_remap = malloc((num_queries+1) * sizeof(char*));
	if (!query_remap) exit(EXIT_FAILURE);
    for (int i = 0; i < num_queries; i++){
        query_remap[i] = malloc(QUERY_SIZE * sizeof(char));
		if (!query_remap[i]) exit(EXIT_FAILURE);
    }
    queries_disp[num_queries] = query_map_s;
    
    for (int i = 0; i < num_queries; i++){
        int k = 0;
        for (int j = queries_disp[i]; j < queries_disp[i+1]; j++) {
            int len = strlen(query_map[j]);   
            memcpy(query_remap[i]+k, query_map[j], len);
            k += len;
        }        
    }

}

void map_output_string(int num_procs)
{
    out_map_s = proc_data.size*num_procs*OUT_MAP_MULTIPLIER*100; // heuristic of how many found objects

    output_map = malloc(out_map_s);
    if (!output_map) exit(EXIT_FAILURE);
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

// (Master) Send data to slaves 
void send_data(int num_procs) 
{
    for (int i = 1; i < num_procs; i++) {
        MPI_Send(query, query_s, MPI_CHAR, i, i, MPI_COMM_WORLD);
        MPI_Send(dna, dna_s, MPI_CHAR, i, i*10, MPI_COMM_WORLD);
    }
}

// (Slaves) Get data from master 
void get_data(int id)
{
    MPI_Status st;
    query = malloc(MAX_QUERY*QUERY_SIZE);
    dna = malloc(MAX_DNA*DNA_SIZE);
    if (!query || !dna) exit(EXIT_FAILURE);

    MPI_Recv(query, MAX_QUERY*QUERY_SIZE, MPI_CHAR, 0, id, MPI_COMM_WORLD, &st); 
    query_s = st._ucount;
    MPI_Recv(dna, MAX_DNA*DNA_SIZE, MPI_CHAR, 0, id*10, MPI_COMM_WORLD, &st); 
    dna_s = st._ucount;

}

// (Slaves) Send results to master
void send_results(int id, long writen)
{
    MPI_Send(output_map, writen, MPI_CHAR, 0, id, MPI_COMM_WORLD);
}

// (Master) Get results from slaves
void get_results(int id, int num_procs, long writen)
{
    void *buff = malloc((MAX_QUERY*QUERY_SIZE));
    if (!buff) exit(EXIT_FAILURE);
    MPI_Status st;

    fout_s = writen;

    for (int i = 1; i < num_procs; i++) {
        MPI_Recv(buff, (MAX_QUERY*QUERY_SIZE), MPI_CHAR, i, i, MPI_COMM_WORLD, &st);
        long len = st._ucount;
        memcpy(output_map+fout_s, buff, len);
        fout_s += len;
    }

    free(buff);
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
