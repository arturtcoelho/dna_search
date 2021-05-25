#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h> 

#include <errno.h>
#include <stdarg.h>

// MAX char table (ASCII)
#define MAX 256
#define QUERY_LINE_SIZE 1000
#define DESC_LINE_SIZE 100

// Global files
FILE *fdatabase, *fquery, *fout;
char *dna, *query;
long dna_s, query_s;

// Global arrays of bases 
char *bases, *str;

void check (int test, const char * message, ...);
void map_files();
void openfiles();
int bmhs(char *string, int n, char *substr, int m);
void closefiles();
void remove_eol(char *line);
void alloc_arrays();

int main()
{

	alloc_arrays();
	openfiles();
	map_files();

	char desc_dna[DESC_LINE_SIZE], desc_query[DESC_LINE_SIZE];
	char line[QUERY_LINE_SIZE];
	int i, found, result;

	// double time_spent0 = 0;
	double time_spent1 = 0;
	clock_t begin, end;
	
	fgets(desc_query, DESC_LINE_SIZE, fquery);
	remove_eol(desc_query);

	printf("%ld\n", dna_s);
	printf("%ld\n", query_s);

	while (!feof(fquery)) {
		fprintf(fout, "%s\n", desc_query);
		
		// read query string
		fgets(line, QUERY_LINE_SIZE, fquery);
		remove_eol(line);
		str[0] = 0;
		i = 0;
		do {
			strcat(str + i, line);
			if (!fgets(line, QUERY_LINE_SIZE, fquery))
				break;
			remove_eol(line);
			i += QUERY_LINE_SIZE-1;
		} while (line[0] != '>');
		strcpy(desc_query, line);

		// read database and search
		found = 0;
		fseek(fdatabase, 0, SEEK_SET);
		fgets(line, DESC_LINE_SIZE, fdatabase);
		remove_eol(line);
		while (!feof(fdatabase)) {
			// read dna section
			begin = clock();
			
			strcpy(desc_dna, line);
			bases[0] = 0;
			i = 0;
			fgets(line, 100, fdatabase);
			remove_eol(line);
			do {
				strcat(bases + i, line);
				if (!fgets(line, 100, fdatabase))
					break;
				remove_eol(line);
				i += 80;
			} while (line[0] != '>');

			end = clock();
			time_spent1 += (double)(end - begin) / CLOCKS_PER_SEC;

			// search with str in bases

			result = bmhs(bases, strlen(bases), str, strlen(str));
			if (result > 0) {
				fprintf(fout, "%s\n%d\n", desc_dna, result);
				found++;
			}
		}

		if (!found)
			fprintf(fout, "NOT FOUND\n");

	}

	closefiles();
	free(str);
	free(bases);

	// printf("time 0 = %f\n", time_spent0);
	printf("time 1 = %f\n", time_spent1);

	return EXIT_SUCCESS;
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

void openfiles() 
{

	fdatabase = fopen("dna.in", "r+");
	if (fdatabase == NULL) {
		perror("dna.in");
		exit(EXIT_FAILURE);
	}

	fquery = fopen("query.in", "r");
	if (fquery == NULL) {
		perror("query.in");
		exit(EXIT_FAILURE);
	}

	fout = fopen("dna.out", "w");
	if (fout == NULL) {
		perror("fout");
		exit(EXIT_FAILURE);
	}

}

void closefiles() 
{
	fflush(fdatabase);
	fclose(fdatabase);

	fflush(fquery);
	fclose(fquery);

	fflush(fout);
	fclose(fout);
}

void remove_eol(char *line) 
{
	int i = strlen(line) - 1;
	while (line[i] == '\n' || line[i] == '\r') {
		line[i] = 0;
		i--;
	}
}

void alloc_arrays()
{
	bases = malloc(sizeof(char) * 1000001);
	if (bases == NULL) {
		perror("malloc");
		exit(EXIT_FAILURE);
	}

	str = malloc(sizeof(char) * 1000001);
	if (str == NULL) {
		perror("malloc str");
		exit(EXIT_FAILURE);
	}
}
