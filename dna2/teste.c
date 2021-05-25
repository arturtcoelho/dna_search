/* For the size of the file. */
#include <sys/stat.h>
/* This contains the mmap calls. */
#include <sys/mman.h> 
/* These are for error printing. */
#include <errno.h>
#include <string.h>
#include <stdarg.h>
/* This is for open. */
#include <fcntl.h>
#include <stdio.h>
/* For exit. */
#include <stdlib.h>
/* For the final part of the example. */
#include <ctype.h>

char *dna, *query;
long dna_s, query_s;

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
    size_t size;
    /* The file name to open. */
    char * file_name = "dna.in";
    /* The memory-mapped thing itself. */
    char * mapped;
    int i;

    /* Open the file for reading. */
    fd = open (file_name, O_RDONLY);
    check (fd < 0, "open %s failed: %s", file_name, strerror (errno));

    /* Get the size of the file. */
    status = fstat (fd, & s);
    check (status < 0, "stat %s failed: %s", file_name, strerror (errno));
    dna_s = s.st_size;

    /* Memory-map the file. */
    dna = mmap (0, size, PROT_READ, MAP_PRIVATE, fd, 0);
    check (mapped == MAP_FAILED, "mmap %s failed: %s",
           file_name, strerror (errno));

    file_name = "query.in";

    fd = open (file_name, O_RDONLY);
    check (fd < 0, "open %s failed: %s", file_name, strerror (errno));

    /* Get the size of the file. */
    status = fstat (fd, & s);
    check (status < 0, "stat %s failed: %s", file_name, strerror (errno));
    query_s = s.st_size;

    /* Memory-map the file. */
    query = mmap (0, size, PROT_READ, MAP_PRIVATE, fd, 0);
    check (mapped == MAP_FAILED, "mmap %s failed: %s",
           file_name, strerror (errno));
}

int main ()
{
    map_files();
    int i = 0;
    while (dna[i++] != '\n') putc(dna[i], stdout);
    i = 0;
    while (query[i++] != '\n') putc(query[i], stdout);
    return 0;
}