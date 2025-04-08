#ifndef __BENCH
#define __BENCH

#ifdef __cplusplus
extern "C" {
#endif
#include <openssl/rand.h>
#include <math.h>
#ifdef __cplusplus
}
#endif

double bench_pcg(size_t n, size_t c, size_t t,int party, int port);
#endif