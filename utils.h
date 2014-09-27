#ifndef UTILS_H
#define UTILS_H

void  *xmalloc(size_t len);
void  *xcalloc(size_t nmemb, size_t size);
void  *xrealloc(void *p, size_t len);
void   die(const char *errstr, ...);

#endif
