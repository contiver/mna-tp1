#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "utils.h"

void *
xrealloc(void *p, size_t len) {
    if( !(p = realloc(p, len)) )
        die("Out of memory\n");
    return p;
}

void *
xmalloc(size_t len) {
    void *p = malloc(len);

    if(!p)
        die("Out of memory.\n");
    return p;
}

void *
xcalloc(size_t nmemb, size_t size) {
    void *p = calloc(nmemb, size);

    if(!p)
        die("Out of memory.\n");
    return p;
}

void
die(const char *errstr, ...) {
    va_list ap;

    va_start(ap, errstr);
    vfprintf(stderr, errstr, ap);
    va_end(ap);
    exit(EXIT_FAILURE);
}
