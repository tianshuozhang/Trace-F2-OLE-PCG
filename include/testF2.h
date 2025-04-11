#pragma once

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
class Beaver_Triple_Seed {
    public:
        
        explicit Beaver_Triple_Seed(size_t t);  // 声明
    
        // 成员变量保持public（根据原始代码设计）
        std::vector<uint8_t> err_poly_coeffs_A;
        std::vector<size_t> err_poly_positions_A;
        std::vector<uint8_t> s_poly_coeffs_A;
        std::vector<size_t> s_poly_positions_A;
    };

void Beaver_Triple_test(const size_t n,const size_t t ,int party,int port);

void sample_a_and_tensor(uint8_t *fft_a, uint32_t *fft_a2, uint32_t *fft_a2a,size_t poly_size);

void test_pcg_F2(const size_t n,const size_t t ,int party,int port,
    std::vector<uint8_t> &fft_a,std::vector<uint8_t>&fft_a2,std::vector<uint8_t>&fft_a2a,
    Beaver_Triple_Seed &beaver_triple_seed,
    std::vector<uint8_t> &z_poly_A,std::vector<uint8_t> &x_poly_A);

void squre_poly(const size_t n, const size_t t,const size_t block_size,
    std::vector<uint8_t> &err_poly_coeffs_A, std::vector<size_t> &err_poly_positions_A,
    std::vector<uint8_t> &err_polys_A2,std::vector<uint8_t> &err_poly_coeffs_A2, std::vector<size_t> &err_poly_positions_A2);

void get_full_poly(const size_t n, const size_t t,const size_t block_size,std::vector<uint8_t> &err_polys_A,
    std::vector<uint8_t> &err_poly_coeffs_A, std::vector<size_t> &err_poly_positions_A);

void sample_local_err(const size_t n, const size_t t,const size_t block_size,
    std::vector<uint8_t> &err_poly_coeffs_A, std::vector<size_t> &err_poly_positions_A);

#endif
