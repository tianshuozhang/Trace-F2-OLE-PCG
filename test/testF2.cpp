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

#include"baseole.h"
#include"otdpf.h"

void sample_a_and_tensor(uint8_t *fft_a, uint32_t *fft_a2, uint32_t *fft_a2a,size_t poly_size, size_t c) {

   

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
    for (size_t j = 0; j < c; j++) {
        for (size_t k = 0; k < c; k++) {
            for (size_t i = 0; i < poly_size; i++) {
                uint8_t u = (fft_a[i] >> (2 * j)) & 0b11;
                uint8_t v = (fft_a[i] >> (2 * k)) & 0b11;
                uint32_t w = mult_f4(u,v);
                uint32_t prod = mult_f4(u, w);
                size_t slot = j * c + k;
                fft_a2[i] |= w<<(2*slot);
                fft_a2a[i] |= prod <<(2*slot);
            }
        }
    }
    printf("Done with sampling the public values\n");
}
void test_pcg_F2(const size_t n,const size_t c,const size_t t ,int party,int port){
    clock_t time;
    time = clock();
    clock_t start = clock();
    double time_taken;
    const size_t poly_size = ipow(3, n);

    //************************************************************************
    // Step 0: Sample the global (1, a1 ... a_c-1) polynomials
    //************************************************************************

    std::vector<uint8_t> fft_a(poly_size);
    
    std::vector<uint32_t>fft_a2(poly_size);

    std::vector<uint32_t>fft_a2a(poly_size);
    NetIO* io = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port);
    
    
    if (party==1){
        sample_a_and_tensor(fft_a.data(), fft_a2.data(), fft_a2a.data(), poly_size, c);
        io->send_data(fft_a.data(),poly_size);
        io->flush();
        io->send_data(fft_a2.data(),poly_size*4);
        io->flush();
        io->send_data(fft_a2a.data(),poly_size*4);
        io->flush();
    }
    else{
        io->recv_data(fft_a.data(),poly_size);
        io->recv_data(fft_a2.data(),poly_size*4);
        io->recv_data(fft_a2a.data(),poly_size*4);
    }

    //************************************************************************
    // Here, we figure out a good block size for the error vectors such that
    // t*block_size = 3^n and block_size/L*128 is close to a power of 3.
    // We pack L=256 coefficients of F4 into each DPF output (note that larger
    // packing values are also okay, but they will do increase key size).
    //************************************************************************
    const size_t dpf_domain_bits = n-log_base(t,3);
    

    printf("DPF domain bits %zu \n", dpf_domain_bits);

    
    const size_t dpf_block_size = pow(3,dpf_domain_bits);

    printf("dpf_block_size = %zu\n", dpf_block_size);

    // Note: We assume that t is a power of 3 and so it divides poly_size
    const size_t block_size = poly_size / t;

    printf("block_size = %zu \n", block_size);

    printf("[       ]Done with Step 0 (sampling the public values)\n");



    //************************************************************************
    std::vector<uint8_t> err_polys_A(c * poly_size);

    // coefficients associated with each error vector
    
    std::vector<uint8_t> err_poly_coeffs_A(c * t);

    // positions of the T errors in each error vector
   
    std::vector<size_t> err_poly_positions_A(c * t);


    std::vector<uint8_t> err_polys_A2(c * poly_size);

    // coefficients associated with each error vector
    
    std::vector<uint8_t> err_poly_coeffs_A2(c * t);

    // positions of the T errors in each error vector
   
    std::vector<size_t> err_poly_positions_A2(c * t);

    

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

            
            a = mult_f4(a,a);

            err_poly_positions_A2[offset] = pos_A;

            err_poly_coeffs_A2[offset] = a;

            err_polys_A2[i * poly_size + j * block_size + pos_A] = a;

            
        }
        
    }
    std::vector<uint8_t> shareA(c*c*t*t);
    std::vector<uint8_t> shareA2(c*c*t*t);
    if(party==1){
        sender_OT(io,err_poly_coeffs_A,c,t,shareA);
        sender_OT(io,err_poly_coeffs_A,c,t,shareA2);
       
    }
    else{
        receiver_OT(io,err_poly_coeffs_A,c,t,shareA);
        receiver_OT(io,err_poly_coeffs_A2,c,t,shareA2);
    }
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
    fft_recursive_uint8(fft_eA.data(), n, poly_size / 3);
    if(party==1){
        printf("[.      ]Done with Step 1 (sampling error vectors)\n");

        time_taken = ((double)clock() - time) / (CLOCKS_PER_SEC / 1000.0); // ms
        printf("Step 1 time elapsed %f ms\n", time_taken);
        time = clock();
    }

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
    // Apply Trace Function
    for (size_t i = 0; i < poly_size; i++)
    {
        x_poly_A[i] = (x_poly_A[i]^mult_f4(x_poly_A[i],x_poly_A[i])) & 0b01;
        
    }

    if(party==1){
        printf("[..     ]Done with Step 2 (computing the local vectors)\n");
        time_taken = ((double)clock() - time) / (CLOCKS_PER_SEC / 1000.0); // ms
        printf("Step 2 time elapsed %f ms\n", time_taken);
        time = clock();
    }
    // Sample PRF keys for the DPFs
    struct PRFKeys *prf_keys = new PRFKeys;
    PRFKeyGen(prf_keys);
    std::vector<uint8_t> trit_decomp_A(log_base(t,3));
    std::vector<uint8_t> trit_decomp_B(log_base(t,3));
    std::vector<DPFParty> dpf_party(c * c * t * t);
    if(party==1){
    
        for (size_t i = 0; i < c; i++) {
            for (size_t j = 0; j < c; j++) {
                std::vector<size_t> next_idx(t,0);
    
                for (size_t k = 0; k < t; k++) {
                    for (size_t l = 0; l < t; l++) {
                        int_to_trits(k, trit_decomp_A.data(), n-dpf_domain_bits);
                        int_to_trits(l, trit_decomp_B.data(), n-dpf_domain_bits);
                        for (size_t k = 0; k < trit_decomp_A.size(); k++) {
                            // printf("[DEBUG]: trits_A[%zu]=%i, trits_B[%zu]=%i\n",
                            //    k, trit_decomp_A[k], k, trit_decomp_B[k]);
                            trit_decomp_A[k] = (trit_decomp_A[k] + trit_decomp_B[k]) % 3;
                        }
                        const size_t pos = trits_to_int(trit_decomp_A.data(), trit_decomp_A.size());
                        const size_t idx = next_idx[pos];
                        next_idx[pos]++;
    
                        const size_t index = i * c * t * t + j * t * t + k * t + l;
                        // Message (beta) is of size 4 blocks of 128 bits
                        std::vector<uint128_t>beta;
                        
                        const uint128_t share_value = static_cast<__uint128_t>(shareA[index]&0b11);
                        
                        beta.push_back(share_value);
                        DPFParty dpf(prf_keys,dpf_domain_bits,err_poly_positions_A[i*t+k],beta,party);
                        const uint128_t offset =i * c * t * t + j * t * t + pos * t + idx;
                        dpf.generate(io);
    
                        dpf_party[offset] = dpf;                 
                    }
                }
            }
        }
    
        }
        else{
            for (size_t i = 0; i < c; i++) {
                for (size_t j = 0; j < c; j++) {
                    
                    std::vector<size_t> next_idx(t,0);
    
                    for (size_t k = 0; k < t; k++) {
                        for (size_t l = 0; l < t; l++) {
                            int_to_trits(k, trit_decomp_A.data(), n-dpf_domain_bits);
                            int_to_trits(l, trit_decomp_B.data(), n-dpf_domain_bits);
                            for (size_t k = 0; k < trit_decomp_A.size(); k++) {
                                // printf("[DEBUG]: trits_A[%zu]=%i, trits_B[%zu]=%i\n",
                                //    k, trit_decomp_A[k], k, trit_decomp_B[k]);
                                trit_decomp_A[k] = (trit_decomp_A[k] + trit_decomp_B[k]) % 3;
                            }
                            const size_t pos = trits_to_int(trit_decomp_A.data(), trit_decomp_A.size());
                            const size_t idx = next_idx[pos];
                            next_idx[pos]++;
                            const size_t index = i * c * t * t + j * t * t + k * t + l;
                            // Message (beta) is of size 4 blocks of 128 bits
                            std::vector<uint128_t>beta;
                            const uint128_t share_value = static_cast<__uint128_t>(shareA[index]&0b11);
                        
                            beta.push_back(share_value);
                            DPFParty dpf(prf_keys,dpf_domain_bits,err_poly_positions_A[j*t+l],beta,party);
        
                            dpf.generate(io);
    
                            const uint128_t offset =i * c * t * t + j * t * t + pos * t + idx;
                            dpf_party[offset] = dpf;      
    
                            
                            
                        }
                    }
                }
            }
        }
    std::vector<DPFParty> dpf_partyA2(c * c * t * t);
    if(party==1){
    
        for (size_t i = 0; i < c; i++) {
            for (size_t j = 0; j < c; j++) {
                std::vector<size_t> next_idx(t,0);
    
                for (size_t k = 0; k < t; k++) {
                    for (size_t l = 0; l < t; l++) {
                        int_to_trits(k, trit_decomp_A.data(), n-dpf_domain_bits);
                        int_to_trits(l, trit_decomp_B.data(), n-dpf_domain_bits);
                        for (size_t k = 0; k < trit_decomp_A.size(); k++) {
                            // printf("[DEBUG]: trits_A[%zu]=%i, trits_B[%zu]=%i\n",
                            //    k, trit_decomp_A[k], k, trit_decomp_B[k]);
                            trit_decomp_A[k] = (trit_decomp_A[k] + trit_decomp_B[k]) % 3;
                        }
                        const size_t pos = trits_to_int(trit_decomp_A.data(), trit_decomp_A.size());
                        const size_t idx = next_idx[pos];
                        next_idx[pos]++;
    
                        const size_t index = i * c * t * t + j * t * t + k * t + l;
                        // Message (beta) is of size 4 blocks of 128 bits
                        std::vector<uint128_t>beta;
                        
                        const uint128_t share_value = static_cast<__uint128_t>(shareA2[index]&0b11);
                        
                        beta.push_back(share_value);
                        DPFParty dpf(prf_keys,dpf_domain_bits,err_poly_positions_A[i*t+k],beta,party);
                        const uint128_t offset =i * c * t * t + j * t * t + pos * t + idx;
                        dpf.generate(io);
    
                        dpf_partyA2[offset] = dpf;                 
                    }
                }
            }
        }
    
        }
        else{
            for (size_t i = 0; i < c; i++) {
                for (size_t j = 0; j < c; j++) {
                    
                    std::vector<size_t> next_idx(t,0);
    
                    for (size_t k = 0; k < t; k++) {
                        for (size_t l = 0; l < t; l++) {
                            int_to_trits(k, trit_decomp_A.data(), n-dpf_domain_bits);
                            int_to_trits(l, trit_decomp_B.data(), n-dpf_domain_bits);
                            for (size_t k = 0; k < trit_decomp_A.size(); k++) {
                                // printf("[DEBUG]: trits_A[%zu]=%i, trits_B[%zu]=%i\n",
                                //    k, trit_decomp_A[k], k, trit_decomp_B[k]);
                                trit_decomp_A[k] = (trit_decomp_A[k] + trit_decomp_B[k]) % 3;
                            }
                            const size_t pos = trits_to_int(trit_decomp_A.data(), trit_decomp_A.size());
                            const size_t idx = next_idx[pos];
                            next_idx[pos]++;
                            const size_t index = i * c * t * t + j * t * t + k * t + l;
                            // Message (beta) is of size 4 blocks of 128 bits
                            std::vector<uint128_t>beta;
                            const uint128_t share_value = static_cast<__uint128_t>(shareA2[index]&0b11);
                        
                            beta.push_back(share_value);
                            DPFParty dpf(prf_keys,dpf_domain_bits,err_poly_positions_A2[j*t+l],beta,party);
        
                            dpf.generate(io);
    
                            const uint128_t offset =i * c * t * t + j * t * t + pos * t + idx;
                            dpf_partyA2[offset] = dpf;      
    
                            
                            
                        }
                    }
                }
            }
        }
}

