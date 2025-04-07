
#ifdef __cplusplus
extern "C" {
#endif
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
// #include "test.h"
#include "prf.h"
#include "fft.h"
#include "utils.h"
#include "f4ops.h"
#ifdef __cplusplus

}
#endif

#include"otdpf.h"



#define MAX(a, b) ((a) > (b) ? (a) : (b))

// Benchmarks are less documented compared to test.c; see test.c to
// better understand what is done here for timing purposes.

#define DPF_MSG_SIZE 8

// samples the a polynomials and axa polynomials
void sample_a_and_a2(uint8_t *fft_a, uint32_t *fft_a2, size_t poly_size, size_t c) {
    RAND_bytes((uint8_t *)fft_a, sizeof(uint8_t) * poly_size);

    // make a_0 the identity polynomial (in FFT space) i.e., all 1s
    for (size_t i = 0; i < poly_size; i++) {
        fft_a[i] = fft_a[i] >> 2;
        fft_a[i] = fft_a[i] << 2;
        fft_a[i] |= 1;
    }
    // FOR DEBUGGING: set fft_a to the identity
    // for (size_t i = 0; i < poly_size; i++)
    // {
    //     fft_a[i] = (0xaaaa >> 1);
    // }

    uint32_t prod;
    for (size_t j = 0; j < c; j++) {
        for (size_t k = 0; k < c; k++) {
            for (size_t i = 0; i < poly_size; i++) {
                prod = mult_f4((fft_a[i] >> (2 * j)) & 0b11, (fft_a[i] >> (2 * k)) & 0b11);
                size_t slot = j * c + k;
                fft_a2[i] |= prod << (2 * slot);
            }
        }
    }
}

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
    size_t dpf_domain_bits = ceil(log_base(poly_size / (t * DPF_MSG_SIZE * 64), 3));
    printf("dpf_domain_bits = %zu \n", dpf_domain_bits);

    size_t seed_size_bits = (128 * (dpf_domain_bits * 3 + 1)*DPF_MSG_SIZE) * c * c * t * t;
    printf("PCG seed size: %.2f MB\n", seed_size_bits / 8000000.0);

    size_t dpf_block_size = DPF_MSG_SIZE * ipow(3, dpf_domain_bits);
    size_t block_size = ceil(poly_size / t);

    printf("block_size = %zu \n", block_size);

    
    std::vector<DPFParty> dpf_party;
    

    // Sample PRF keys for the DPFs
    struct PRFKeys *prf_keys = new PRFKeys;
    PRFKeyGen(prf_keys);
    
    
    
    // Sample DPF keys for each of the t errors in the t blocks
    for (size_t i = 0; i < c; i++)
    {
        for (size_t j = 0; j < c; j++)
        {
            for (size_t k = 0; k < t; k++)
            {
                for (size_t l = 0; l < t; l++)
                {
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
                    dpf_party.push_back(dpf);
                    // printf("current index = %zu\n", index);
                    
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
                poly_block = &packed_poly[k * packed_block_size];

                for (l = 0; l < t; l++)
                {
                    key_index = i * c * t * t + j * t * t + k * t + l;
                    auto dpf = dpf_party[key_index];
                    std::vector<uint128_t> shares;
                    dpf.fulldomainevaluation(shares);
                    
                    for (w = 0; w < packed_block_size; w++)
                        poly_block[w] ^= shares[w];

                }
            }
        }
    }
    if(party==1){
        int nums = 0;
        std::vector<uint128_t>receive(packed_polys.size());
        io->recv_data(receive.data(),packed_polys.size()*16);
        for(int w = 0; w < packed_polys.size(); ++w){
            if((packed_polys[w]^receive[w])!= 0 ) ++nums;
        }
        std::cout<<"非零个数：\t"<<nums<<"\n";
    }
    else{
        io->send_data(packed_polys.data(),packed_polys.size()*16);
        io->flush();
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

int main(int argc , char** argv){
    
    int party, port;
    parse_party_and_port(argv, &party, &port);
    size_t c = 4;
    size_t t = 9;
    size_t n = 10;
    std::cout<<"party:\t"<<party<<std::endl;
    bench_pcg(n,c,t,party,port);
    
}