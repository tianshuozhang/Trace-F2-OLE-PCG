
#ifdef __cplusplus
extern "C" {
#endif
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "prf.h"
#include "fft.h"
#include "utils.h"
#include "f4ops.h"
#ifdef __cplusplus

}
#endif
#include"bench.h"
#include"otdpf.h"
#include"testF4.h"


#define MAX(a, b) ((a) > (b) ? (a) : (b))

// Benchmarks are less documented compared to test.c; see test.c to
// better understand what is done here for timing purposes.

#define DPF_MSG_SIZE 1

double bench_pcg(size_t n, size_t c, size_t t,int party, int port)
{
    if (c > 4)
    {
        printf("ERROR: currently only implemented for c <= 4");
        exit(0);
    }

    const size_t poly_size = ipow(3, n);

    //************************************************************************
    // Step 0: Sample the global (1, a1 ... a_c-1) polynomials
    //************************************************************************
    NetIO* io = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port);
    std::vector<uint8_t> fft_a(poly_size);
    
    std::vector<uint32_t>fft_a2(poly_size);
    
    if (party==1){
        sample_a_and_a2(fft_a.data(), fft_a2.data(), poly_size, c);
        io->send_data(fft_a.data(),poly_size);
        io->flush();
        io->send_data(fft_a2.data(),poly_size*4);
        io->flush();
    }
    else{
        io->recv_data(fft_a.data(),poly_size);
        io->recv_data(fft_a2.data(),poly_size*4);
    }

    //************************************************************************
    // Step 1: Sample DPF keys for the cross product.
    // For benchmarking purposes, we sample random DPF functions for a
    // sufficiently large domain size to express a block of coefficients.
    //************************************************************************
    size_t dpf_domain_bits = ceil(log_base(poly_size / (t * DPF_MSG_SIZE), 3));
    printf("dpf_domain_bits = %zu \n", dpf_domain_bits);


    size_t dpf_block_size = DPF_MSG_SIZE * ipow(3, dpf_domain_bits);
    size_t block_size = ceil(poly_size / t);

    printf("block_size = %zu \n", block_size);

    
    std::vector<DPFParty> dpf_party(c * c * t * t);
    

    // Sample PRF keys for the DPFs
    struct PRFKeys *prf_keys = new PRFKeys;
    PRFKeyGen(prf_keys);
    std::vector<uint8_t> trit_decomp_A(log_base(t,3));
    std::vector<uint8_t> trit_decomp_B(log_base(t,3));
    
    
    // Sample DPF keys for each of the t errors in the t blocks
    for (size_t i = 0; i < c; i++)
    {
        for (size_t j = 0; j < c; j++)
        {
            std::vector<uint8_t> next_idx(t,0);
            for (size_t k = 0; k < t; k++)
            {
                for (size_t l = 0; l < t; l++)
                {
                    int_to_trits(k, trit_decomp_A.data(), log_base(t,3));
                    int_to_trits(l, trit_decomp_B.data(), log_base(t,3));
                    for (size_t k = 0; k < trit_decomp_A.size(); k++) {
                        // printf("[DEBUG]: trits_A[%zu]=%i, trits_B[%zu]=%i\n",
                        //    k, trit_decomp_A[k], k, trit_decomp_B[k]);
                        trit_decomp_A[k] = (trit_decomp_A[k] + trit_decomp_B[k]) % 3;
                    }
                    size_t pos = trits_to_int(trit_decomp_A.data(), trit_decomp_A.size());
                    size_t idx = next_idx[pos];
                    next_idx[pos]++;

                    size_t index = i * c * t * t + j * t * t + k * t + l;

                    // Pick a random index for benchmarking purposes
                    size_t alpha = random_index(block_size);

                    // Pick a random output message for benchmarking purposes
                    // uint128_t beta[DPF_MSG_SIZE] = {0};
                    std::vector<uint128_t> beta(DPF_MSG_SIZE);

                    RAND_bytes((uint8_t *)beta.data(), DPF_MSG_SIZE * sizeof(uint128_t));

                    

                    // Message (beta) is of size 8 blocks of 128 bits
                    DPFParty dpf(prf_keys,dpf_domain_bits,alpha,beta,party);
                    dpf.generate(io);
                    
                    // printf("current index = %zu\n", index);
                    const uint32_t offset =i * c * t * t + j * t * t + pos * t + idx;
                    dpf_party[offset] = dpf;   
                }
            }
        }
    }

    //************************************************
    printf("Benchmarking PCG evaluation \n");
    //************************************************

    // Allocate memory for the DPF outputs (this is reused for each evaluation)
    
    

    // Allocate memory for the concatenated DPF outputs
    const size_t packed_block_size = ceil(block_size / 64.0);
    const size_t packed_poly_size = t * packed_block_size;

    std::vector<uint128_t>packed_polys(c * c * packed_poly_size,0);

    // Allocate memory for the output FFT
    std::vector<uint32_t> fft_u(poly_size,0);

    // Allocate memory for the final inner product
    std::vector<uint8_t> z_poly(poly_size,0);
    std::vector<uint32_t> res_poly_mat(poly_size);

    //************************************************************************
    // Step 3: Evaluate all the DPFs to recover shares of the c*c polynomials.
    //************************************************************************

    clock_t time;
    time = clock();
    
    
    size_t key_index;
    uint128_t *poly_block;
    size_t i, j, k, l, w;
    
    for (i = 0; i < c; i++)
    {
        for (j = 0; j < c; j++)
        {
            const size_t poly_index = i * c + j;
            uint128_t *packed_poly = &packed_polys[poly_index * packed_poly_size];

            for (k = 0; k < t; k++)
            {
                std::vector<uint128_t> totalsharesA(ipow(3,dpf_domain_bits));
                poly_block = &packed_poly[k * packed_block_size];

                for (l = 0; l < t; l++)
                {
                    key_index = i * c * t * t + j * t * t + k * t + l;
                    auto dpf = dpf_party[key_index];
                    std::vector<uint128_t> shares;
                    dpf.fulldomainevaluation(shares);
                    
                    for(size_t w = 0; w < shares.size(); w++) {
                        totalsharesA[w] ^= shares[w];
                    }

                }
                for (size_t w = 0; w < totalsharesA.size(); w++) {
                    // 提取低2位
                    uint128_t lsb = totalsharesA[w] & 0b11;
                
                    // 计算存储位置
                    size_t block_index = w / 64;
                    size_t bit_offset = (63-w % 64) * 2;
                
                    // 边界检查（防止越界）
                    if (block_index >= packed_block_size) {
                        break; // 或处理错误
                    }
                    
                    // 写入到目标位置
                    poly_block[block_index] |= (lsb << bit_offset);
                    
                }
            }
        }
    }
    
    //************************************************************************
    // Step 3: Compute the transpose of the polynomials to pack them into
    // the parallel FFT format.
    //
    // TODO: this is the bottleneck of the computation and can be improved
    // using SIMD operations for performing matrix transposes (see TODO in test.c).
    //************************************************************************
    for (size_t i = 0; i < c * c; i++)
    {
        size_t poly_index = i * packed_poly_size;
        const uint128_t *poly = &packed_polys[poly_index];
        __builtin_prefetch(&poly[0], 0, 3);

        size_t block_idx, packed_coeff_idx, coeff_idx;
        uint8_t packed_bit_idx;
        uint128_t packed_coeff;

        block_idx = 0;
        packed_coeff_idx = 0;
        coeff_idx = 0;

        for (size_t k = 0; k < poly_size - 64; k += 64)
        {
            packed_coeff = poly[block_idx * packed_block_size + packed_coeff_idx];
            __builtin_prefetch(&fft_u[k], 0, 0);
            __builtin_prefetch(&fft_u[k], 1, 0);

            for (size_t l = 0; l < 64; l++)
            {
                packed_coeff = packed_coeff >> 2;
                fft_u[k + l] |= packed_coeff & 0b11;
                fft_u[k + l] = fft_u[k + l] << 2;
            }

            packed_coeff_idx++;
            coeff_idx += 64;

            if (coeff_idx > block_size)
            {
                coeff_idx = 0;
                block_idx++;
                packed_coeff_idx = 0;
                __builtin_prefetch(&poly[block_idx * packed_block_size], 0, 2);
            }
        }

        packed_coeff = poly[block_idx * packed_block_size + packed_coeff_idx];
        for (size_t k = poly_size - 64 + 1; k < poly_size; k++)
        {
            packed_coeff = packed_coeff >> 2;
            fft_u[k] |= packed_coeff & 0b11;
            fft_u[k] = fft_u[k] << 2;
        }
    }
    
    fft_recursive_uint32(fft_u.data(), n, poly_size / 3);
    multiply_fft_32(fft_a2.data(), fft_u.data(), res_poly_mat.data(), poly_size);

    // Perform column-wise XORs to get the result
    for (size_t i = 0; i < poly_size; i++)
    {
        // XOR the (packed) columns into the accumulator
        for (size_t j = 0; j < c * c; j++)
        {
            z_poly[i] ^= res_poly_mat[i] & 0b11;
            res_poly_mat[i] = res_poly_mat[i] >> 2;
        }
    }
    

    time = clock() - time;
    double time_taken = ((double)time) / (CLOCKS_PER_SEC / 1000.0); // ms

    printf("Eval time (total) %f ms\n", time_taken);
    printf("DONE\n\n");

    DestroyPRFKey(prf_keys);
   

    return time_taken;
}
