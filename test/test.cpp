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

// This test case implements Figure 1 from https://eprint.iacr.org/2024/429.pdf.
// It uses /libs/fft and libs/tri-dpf extensively.
// Several simplifying design choices are made:
// 1. We assume that c*c <= 16 so that we can use a parallel FFT packing of F4
// elements using a uint32_t type.
// 2. We assume that t is a power of 3 so that the block size of each error
// vector divides the size of the polynomial. This makes the code significantly
// more readable and easier to understand.

// TODO[feature]: The current implementation assumes that C*C <= 16 in order
// to parallelize the FFTs and other components. Making the code work with
// arbitrary values of C is left for future work.

// TODO[feature]: modularize the different components of the test case and
// design more unit tests.

#define N 16 // 3^N number of OLEs generated in total

// The C and T parameters are computed using the SageMath script that can be
// found in https://github.com/mbombar/estimator_folding

#define C 4  // compression factor
#define T 27 // noise weight

// This test evaluates the full PCG.Expand for both parties and
// checks correctness of the resulting OLE correlation.
void test_pcg()
{
    clock_t time;
    time = clock();

    const size_t n = N;
    const size_t c = C;
    const size_t t = T;
    const size_t poly_size = ipow(3, n);

    //************************************************************************
    // Step 0: Sample the global (1, a1 ... a_c-1) polynomials
    //************************************************************************

    std::vector<uint8_t> fft_a(poly_size);
    
    std::vector<uint32_t>fft_a2(poly_size);
    sample_a_and_a2(fft_a.data(), fft_a2.data(), poly_size, c); // Each is of poly_size

    //************************************************************************
    // Here, we figure out a good block size for the error vectors such that
    // t*block_size = 3^n and block_size/L*128 is close to a power of 3.
    // We pack L=256 coefficients of F4 into each DPF output (note that larger
    // packing values are also okay, but they will do increase key size).
    //************************************************************************
    size_t dpf_domain_bits = n;
    if (dpf_domain_bits == 0)
        dpf_domain_bits = 1;

    printf("DPF domain bits %zu \n", dpf_domain_bits);

    // 4*128 ==> 256 coefficients in F4
    size_t dpf_block_size = poly_size;

    printf("dpf_block_size = %zu\n", dpf_block_size);

    // Note: We assume that t is a power of 3 and so it divides poly_size
    size_t block_size = poly_size / t;

    printf("block_size = %zu \n", block_size);

    printf("[       ]Done with Step 0 (sampling the public values)\n");

    //************************************************************************
    // Step 1: Sample error polynomials eA and eB (c polynomials in total)
    // each polynomial is t-sparse and has degree (t * block_size) = poly_size.
    //************************************************************************
    
    std::vector<uint8_t> err_polys_A(c * poly_size);

    // coefficients associated with each error vector
    
    std::vector<uint8_t> err_poly_coeffs_A(c * t);

    // positions of the T errors in each error vector
   

    std::vector<uint8_t> err_poly_positions_A(c * t);

    for (size_t i = 0; i < c; i++)
    {
        for (size_t j = 0; j < t; j++)
        {
            size_t offset = i * t + j;

            // random *non-zero* coefficients in F4
            uint8_t a = rand_f4x();
            
            err_poly_coeffs_A[offset] = a;
            

            // random index within the block
            size_t pos_A = random_index(block_size - 1);
            

            if (pos_A >= block_size)
            {
                printf("FAIL: position > block_size: %zu\n", pos_A);
                exit(0);
            }

            err_poly_positions_A[offset] = pos_A;
            

            // set the coefficient at the error position to the error value
            err_polys_A[i * poly_size + j * block_size + pos_A] = a;
            
        }
    }

    // Compute FFT of eA and eB in packed form.
    // Note that because c = 4, we can pack 4 FFTs into a uint8_t
    
    std::vector<uint8_t>fft_eA(poly_size);
    
    uint8_t coeff_A;

    // This loop essentially computes a transpose to pack the coefficients
    // of each polynomial into one "row" of the parallel FFT matrix
    for (size_t j = 0; j < c; j++)
    {
        for (size_t i = 0; i < poly_size; i++)
        {
            // extract the i-th coefficient of the j-th error polynomial
            coeff_A = err_polys_A[j * poly_size + i];
            

            // pack the extracted coefficient into the j-th FFT slot
            fft_eA[i] |= (coeff_A << (2 * j));
           
        }
    }

    // Evaluate the FFTs on the error polynomials eA and eB
    fft_recursive_uint8(fft_eA.data(), n, poly_size / 3);
    

    printf("[.      ]Done with Step 1 (sampling error vectors)\n");

    //************************************************************************
    // Step 2: compute the inner product xA = <a, eA> and xB = <a, eB>
    //************************************************************************

    // Initialize polynomials to zero (accumulators for inner product)
    
    std::vector<uint8_t> x_poly_A(poly_size);

    // Compute the coordinate-wise multiplication over the packed FFT result
    
    std::vector<uint8_t> res_poly_A(poly_size);
    multiply_fft_8(fft_a.data(), fft_eA.data(), res_poly_A.data(), poly_size); // a*eA
    

    // XOR the result into the accumulator.
    // Specifically, we XOR all the columns of the FFT result to get a
    // vector of size poly_size.
    for (size_t j = 0; j < c; j++)
    {
        for (size_t i = 0; i < poly_size; i++)
        {
            x_poly_A[i] ^= (res_poly_A[i] >> (2 * j)) & 0b11;
           
        }
    }

    printf("[..     ]Done with Step 2 (computing the local vectors)\n");

    //************************************************************************
    // Step 3: Compute cross product (eA x eB) using the position vectors
    //************************************************************************
    
    printf("[...    ]Done with Step 3 (computing the cross product)\n");

    //************************************************************************
    // Step 4: Sample the DPF keys for the cross product (eA x eB)
    //************************************************************************
    std::vector<DPFParty> dpf_party;
    

    // Sample PRF keys for the DPFs
    struct PRFKeys *prf_keys = new PRFKeys;
    PRFKeyGen(prf_keys);
    int party,port;
    NetIO* io = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port);

    // Sample DPF keys for each of the t errors in the t blocks
    for (size_t i = 0; i < c; i++) {
        for (size_t j = 0; j < c; j++) {
            for (size_t k = 0; k < t; k++) {
                for (size_t l = 0; l < t; l++) {
                    

                    

                    // Message (beta) is of size 4 blocks of 128 bits
                    
                    DPFParty dpf(prf_keys,dpf_domain_bits,err_poly_positions_A[i*t+k],err_poly_coeffs_A[i*t+k],party);

                    dpf.generate(io);
                    dpf_party.push_back(dpf);
                    
                }
            }
        }
    }

    printf("[....   ]Done with Step 4 (sampling DPF keys)\n");

    //************************************************************************
    // Step 5: Evaluate the DPFs to compute shares of (eA x eB)
    //************************************************************************

    // Allocate memory for the DPF outputs (this is reused for each evaluation)
    uint128_t *shares_A = malloc(sizeof(uint128_t) * dpf_block_size);
    uint128_t *shares_B = malloc(sizeof(uint128_t) * dpf_block_size);
    uint128_t *cache = malloc(sizeof(uint128_t) * dpf_block_size);

    // Allocate memory for the concatenated DPF outputs
    size_t packed_block_size = ceil(block_size / 64.0);
    size_t packed_poly_size = t * packed_block_size;

    // printf("[DEBUG]: packed_block_size = %zu\n", packed_block_size);
    // printf("[DEBUG]: packed_poly_size = %zu\n", packed_poly_size);

    uint128_t *packed_polys_A = calloc(c * c * packed_poly_size, sizeof(uint128_t));
    uint128_t *packed_polys_B = calloc(c * c * packed_poly_size, sizeof(uint128_t));

    // Allocate memory for the output FFT
    uint32_t *fft_uA = calloc(poly_size, sizeof(uint32_t));
    uint32_t *fft_uB = calloc(poly_size, sizeof(uint32_t));

    // Allocate memory for the final inner product
    uint8_t *z_poly_A = calloc(poly_size, sizeof(uint8_t));
    uint8_t *z_poly_B = calloc(poly_size, sizeof(uint8_t));
    uint32_t *res_poly_mat_A = malloc(sizeof(uint32_t) * poly_size);
    uint32_t *res_poly_mat_B = malloc(sizeof(uint32_t) * poly_size);

    for (size_t i = 0; i < c; i++) {
        for (size_t j = 0; j < c; j++) {
            const size_t poly_index = i * c + j;
            // each entry is of length packed_poly_size
            uint128_t *packed_polyA = &packed_polys_A[poly_index * packed_poly_size];
            uint128_t *packed_polyB = &packed_polys_B[poly_index * packed_poly_size];

            for (size_t k = 0; k < t; k++) {
                // each entry is of length packed_block_size
                uint128_t *poly_blockA = &packed_polyA[k * packed_block_size];
                uint128_t *poly_blockB = &packed_polyB[k * packed_block_size];

                for (size_t l = 0; l < t; l++) {
                    size_t index = i * c * t * t + j * t * t + k * t + l;
                    struct DPFKey *dpf_keyA = dpf_keys_A[index];
                    struct DPFKey *dpf_keyB = dpf_keys_B[index];

                    DPFFullDomainEval(dpf_keyA, cache, shares_A);
                    DPFFullDomainEval(dpf_keyB, cache, shares_B);

                    // Sum all the DPFs for the current block together
                    // note that there is some extra "garbage" in the last
                    // block of uint128_t since 64 does not divide block_size.
                    // We deal with this slack later when packing the outputs
                    // into the parallel FFT matrix.
                    for (size_t w = 0; w < packed_block_size; w++) {
                        poly_blockA[w] ^= shares_A[w];
                        poly_blockB[w] ^= shares_B[w];
                    }
                }
            }
        }
    }

    // Here, we test to make sure all polynomials have at most t^2 errors
    // and fail the test otherwise.
    for (size_t i = 0; i < c; i++) {
        for (size_t j = 0; j < c; j++) {
            size_t err_count = 0;
            size_t poly_index = i * c + j;
            uint128_t *poly_A = &packed_polys_A[poly_index * packed_poly_size];
            uint128_t *poly_B = &packed_polys_B[poly_index * packed_poly_size];

            for (size_t p = 0; p < packed_poly_size; p++)
            {
                uint128_t res = poly_A[p] ^ poly_B[p];
                for (size_t l = 0; l < 64; l++)
                {
                    if (((res >> (2 * (63 - l))) & 0b11) != 0)
                        err_count++;
                }
            }

            // printf("[DEBUG]: Number of non-zero coefficients in poly (%zu,%zu) is %zu\n", i, j, err_count);

            if (err_count > t * t) {
                printf("FAIL: Number of non-zero coefficients is %zu > t*t\n", err_count);
                exit(0);
            } else if (err_count == 0) {
                printf("FAIL: Number of non-zero coefficients in poly (%zu,%zu) is %zu\n", i, j, err_count);
                exit(0);
            }
        }
    }

    printf("[.....  ]Done with Step 5 (evaluating all DPFs)\n");

    //************************************************************************
    // Step 6: Compute an FFT over the shares of (eA x eB)
    //************************************************************************

    // Pack the coefficients into FFT blocks
    //
    // TODO[optimization]: use AVX and fast matrix transposition algorithms.
    // The transpose is the bottleneck of the current implementation and
    // therefore improving this step can result in significant performance gains.
    size_t poly_index_jk, poly_index_kj, poly_index_jj, ctr;
    uint128_t packedA_jk, packedA_kj, packedA_jj, packedB_jk, packedB_kj, packedB_jj;
    uint32_t coeffA_jk, coeffA_kj, coeffA_jj, coeffB_jk, coeffB_kj, coeffB_jj;

    for (size_t j = 0; j < c; j++) {
        for (size_t k = 0; k < c; k++) {
            uint8_t *test_poly_A = calloc(poly_size, sizeof(uint8_t));
            uint8_t *test_poly_B = calloc(poly_size, sizeof(uint8_t));

            size_t poly_index = j * c + k;
            uint128_t *poly_A = &packed_polys_A[poly_index * packed_poly_size];
            uint128_t *poly_B = &packed_polys_B[poly_index * packed_poly_size];

            size_t block_idx = 0;
            size_t bit_idx = 0;
            for (size_t i = 0; i < poly_size; i++)
            {
                if (i % block_size == 0 && i != 0)
                {
                    block_idx++;
                    bit_idx = 0;
                }

                size_t packed_idx = block_idx * packed_block_size + floor(bit_idx / 64.0);
                size_t packed_bit = (63 - bit_idx % 64);
                test_poly_A[i] = (poly_A[packed_idx] >> (2 * packed_bit)) & 0b11;
                test_poly_B[i] = (poly_B[packed_idx] >> (2 * packed_bit)) & 0b11;

                bit_idx++;
            }

        }
    }

    // TODO[optimization]: for arbitrary values of C, we only need to perform
    // C*(C+1)/2 FFTs which can lead to a more efficient implementation.
    // Because we assume C=4, we have C*C = 16 which fits perfectly into a
    // uint32 packing.

    for (size_t j = 0; j < c; j++)
    {
        for (size_t k = 0; k < c; k++)
        {
            size_t poly_index = (j * c + k) * packed_poly_size;
            const uint128_t *polyA = &packed_polys_A[poly_index];
            const uint128_t *polyB = &packed_polys_B[poly_index];

            size_t block_idx = 0;
            size_t bit_idx = 0;
            for (size_t i = 0; i < poly_size; i++)
            {
                if (i % block_size == 0 && i != 0)
                {
                    block_idx++;
                    bit_idx = 0;
                }

                size_t packed_idx = block_idx * packed_block_size + floor(bit_idx / 64.0);
                size_t packed_bit = (63 - bit_idx % 64);

                // Extract the i-th (packed) coefficient of the (j,k)-th polynomial
                uint128_t packedA = polyA[packed_idx];
                uint128_t packedB = polyB[packed_idx];

                // Extract the i-th coefficient from the packed coefficient
                uint32_t coeffA = (packedA >> (2 * packed_bit)) & 0b11;
                uint32_t coeffB = (packedB >> (2 * packed_bit)) & 0b11;

                // Pack the extracted coefficient into the ctr-th FFT slot
                size_t idx = j * c + k;
                fft_uA[i] |= coeffA << (2 * idx);
                fft_uB[i] |= coeffB << (2 * idx);

                bit_idx++;
            }
        }
    }

    fft_recursive_uint32(fft_uA, n, poly_size / 3);
    fft_recursive_uint32(fft_uB, n, poly_size / 3);

    printf("[...... ]Done with Step 6 (computing FFTs)\n");

    //************************************************************************
    // Step 7: Compute shares of z = <axa, u>
    //************************************************************************
    multiply_fft_32(fft_a2.data(), fft_uA, res_poly_mat_A, poly_size);
    multiply_fft_32(fft_a2.data(), fft_uB, res_poly_mat_B, poly_size);

    size_t num_ffts = c * c;

    // XOR the (packed) columns into the accumulator.
    // Specifically, we perform column-wise XORs to get the result.
    for (size_t j = 0; j < c * c; j++)
    {
        for (size_t i = 0; i < poly_size; i++)
        {
            z_poly_A[i] ^= (res_poly_mat_A[i] >> (2 * j)) & 0b11;
            z_poly_B[i] ^= (res_poly_mat_B[i] >> (2 * j)) & 0b11;
        }
    }

    // Now we check that we got the correct OLE correlations and fail
    // the test otherwise.
    for (size_t i = 0; i < poly_size; i++)
    {
        uint8_t res = z_poly_A[i] ^ z_poly_B[i];
        uint8_t exp = mult_f4(x_poly_A[i], x_poly_B[i]);

        // printf("[DEBUG]: Got: (%i,%i), Expected: (%i, %i)\n",
        //        (res >> 1) & 1, res & 1, (exp >> 1) & 1, exp & 1);

        if (res != exp)
        {
            printf("FAIL: Incorrect correlation output at index %zu\n", i);
            printf("Got: (%i,%i), Expected: (%i, %i)\n",
                   (res >> 1) & 1, res & 1, (exp >> 1) & 1, exp & 1);
            exit(0);
        }
    }

    time = clock() - time;
    double time_taken = ((double)time) / (CLOCKS_PER_SEC / 1000.0); // ms

    printf("[.......]Done with Step 7 (recovering shares)\n\n");

    printf("Time elapsed %f ms\n", time_taken);

    //************************************************************************
    // Step 8: Cleanup
    //************************************************************************
    DestroyPRFKey(prf_keys);
    
}

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
