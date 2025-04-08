#ifndef __TESTF4
#define __TESTF4

#ifdef __cplusplus
extern "C" {
#endif
#include <openssl/rand.h>
#include <math.h>
#ifdef __cplusplus
}
#endif

void test_pcg(const size_t n,const size_t c,const size_t t ,int party,int port);
void sample_a_and_a2(uint8_t *fft_a, uint32_t *fft_a2, size_t poly_size, size_t c);

#endif
