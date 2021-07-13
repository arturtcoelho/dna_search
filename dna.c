#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// MAX char table (ASCII)
#define MAX 256

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

// Global files
FILE *fdatabase, *fquery, *fout;
void openfiles() {

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

void closefiles() {
	fflush(fdatabase);
	fclose(fdatabase);

	fflush(fquery);
	fclose(fquery);

	fflush(fout);
	fclose(fout);
}

void remove_eol(char *line) {
	int i = strlen(line) - 1;
	while (line[i] == '\n' || line[i] == '\r') {
		line[i] = 0;
		i--;
	}
}

// Global arrays of bases 
char *bases;
char *str;
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

int main(void) {

	clock_t beg0 = clock();

	alloc_arrays();
	openfiles();

	char desc_dna[100], desc_query[100];
	char line[100];
	int i, found, result;

	printf("alloc: %f\n", (double)(clock() - beg0) / CLOCKS_PER_SEC);

	clock_t beg1 = clock();
	double total = 0;
	double tot_bhms = 0;

	fgets(desc_query, 100, fquery);
	remove_eol(desc_query);
	while (!feof(fquery)) {
		fprintf(fout, "%s\n", desc_query);

		clock_t beg2 = clock();

		// read query string
		fgets(line, 100, fquery);
		remove_eol(line);
		str[0] = 0;
		i = 0;
		do {
			strcat(str + i, line);
			if (fgets(line, 100, fquery) == NULL)
				break;
			remove_eol(line);
			i += 80;
		} while (line[0] != '>');
		strcpy(desc_query, line);

		total += (double)(clock() - beg2) / CLOCKS_PER_SEC;

		// read database and search
		found = 0;
		fseek(fdatabase, 0, SEEK_SET);
		fgets(line, 100, fdatabase);
		remove_eol(line);
		while (!feof(fdatabase)) {

			clock_t beg3 = clock();

			strcpy(desc_dna, line);
			bases[0] = 0;
			i = 0;
			fgets(line, 100, fdatabase);
			remove_eol(line);
			do {
				strcat(bases + i, line);
				if (fgets(line, 100, fdatabase) == NULL)
					break;
				remove_eol(line);
				i += 80;
			} while (line[0] != '>');

			total += (double)(clock() - beg3) / CLOCKS_PER_SEC;

			clock_t beg5 = clock();

			result = bmhs(bases, strlen(bases), str, strlen(str));

			tot_bhms +=(double)(clock() - beg5) / CLOCKS_PER_SEC;

			if (result > 0) {
				fprintf(fout, "%s\n%d\n", desc_dna, result);
				found++;
			}
		}

		if (!found)
			fprintf(fout, "NOT FOUND\n");
	}

	printf("total read: %f\n", total);
	printf("total bhms: %f\n", tot_bhms);
	printf("total loop: %f\n", (double)(clock() - beg1) / CLOCKS_PER_SEC);

	clock_t beg4 = clock();

	closefiles();
	free(str);
	free(bases);

	printf("free: %f\n", (double)(clock() - beg4) / CLOCKS_PER_SEC);

	clock_t end = clock();
	printf("total: %f\n", (double)(end - beg0) / CLOCKS_PER_SEC);

	return EXIT_SUCCESS;
}
