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
void sample_error_polys_A(std::vector<uint8_t> &err_polys_A, std::vector<uint8_t> &err_poly_coeffs_A, std::vector<size_t> &err_poly_positions_A,
    const  size_t c,const  size_t t, const size_t block_size);
void FFT_on_error(std::vector<uint8_t> &err_polys_A, std::vector<uint8_t> &fft_eA,const size_t n,const size_t c, const size_t poly_size);


void KeyGen(PRFKeys* prf_keys,NetIO *io,int party, const size_t n, const size_t c, const size_t t, const size_t dpf_domain_bits, 
    std::vector<DPFParty> &dpf_party, std::vector<uint8_t> &share,
    std::vector<size_t> &err_poly_positions_A);

void Fullevaluation(std::vector<DPFParty> &dpf_party, const size_t c, const size_t t, 
    const size_t packed_block_size, const size_t packed_poly_size, const size_t block_size, 
    std::vector<uint128_t> &packed_polys_A);

void Evaluation_on_FFT(std::vector<uint128_t> &packed_polys_A, std::vector<uint32_t> &fft_uA, const size_t n, const size_t c, 
    const size_t poly_size,const size_t block_size, const size_t packed_block_size, const size_t packed_poly_size);
#endif
