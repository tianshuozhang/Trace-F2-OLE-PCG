#ifndef __TESTF2
#define __TESTF2

#ifdef __cplusplus
extern "C" {
#endif
#include <openssl/rand.h>
#include <math.h>
#ifdef __cplusplus
}
#endif

void sample_a_and_tensor(uint8_t *fft_a, uint32_t *fft_a2, uint32_t *fft_a2a,size_t poly_size);
void test_pcg_F2(const size_t n,const size_t t ,int party,int port);
void squre_poly(const size_t n, const size_t t,const size_t block_size,
    std::vector<uint8_t> &err_poly_coeffs_A, std::vector<size_t> &err_poly_positions_A,
    std::vector<uint8_t> &err_polys_A2,std::vector<uint8_t> &err_poly_coeffs_A2, std::vector<size_t> &err_poly_positions_A2);
void sample_local_err(const size_t n, const size_t t,const size_t block_size,std::vector<uint8_t> &err_polys_A,
    std::vector<uint8_t> &err_poly_coeffs_A, std::vector<size_t> &err_poly_positions_A);
#endif
