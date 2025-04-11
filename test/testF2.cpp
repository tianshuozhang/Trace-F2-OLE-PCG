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
#include"testF4.h"
#include"testF2.h"

Beaver_Triple_Seed::Beaver_Triple_Seed(size_t t) {
    err_poly_coeffs_A.resize(t);
    err_poly_positions_A.resize(t);
    s_poly_coeffs_A.resize(t);
    s_poly_positions_A.resize(t);
}

void sample_a_and_tensor(uint8_t *fft_a, uint8_t *fft_a2, uint8_t *fft_a2a,size_t poly_size) {
    // RAND_bytes((uint8_t *)fft_a, sizeof(uint8_t) * poly_size);

    // make a_0 the identity polynomial (in FFT space) i.e., all 1s
    for (size_t i = 0; i < poly_size; i++) {
        fft_a[i] = rand_f4x();

        
    }
    for (size_t i = 0; i < poly_size; i++) {
        uint8_t u = fft_a[i]& 0b11;
        uint8_t w = mult_f4(u,u);
        uint8_t prod = mult_f4(u, w);
        
        fft_a2[i] = w& 0b11;
        fft_a2a[i] = prod& 0b11;
    }
        
    printf("Done with sampling the public values\n");
}
void test_pcg_F2(const size_t n,const size_t t ,int party,int port,
     std::vector<uint8_t> &fft_a,std::vector<uint8_t>&fft_a2,std::vector<uint8_t>&fft_a2a,
     Beaver_Triple_Seed &beaver_triple_seed,
     std::vector<uint8_t> &z_poly_A,std::vector<uint8_t> &x_poly_A){
    clock_t time;
    time = clock();
    clock_t start = clock();
    double time_taken;
    const size_t poly_size = ipow(3, n);

    //************************************************************************
    // Step 0: Sample the global (1, a1 ... a_c-1) polynomials
    //************************************************************************

    NetIO* io = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port);
    if(party==1){
        
        io->send_data(fft_a.data(),poly_size);
        io->flush();
        io->send_data(fft_a2.data(),poly_size);
        io->flush();
        io->send_data(fft_a2a.data(),poly_size);
        io->flush();
    }
    else{
        io->recv_data(fft_a.data(),poly_size);
        io->recv_data(fft_a2.data(),poly_size);
        io->recv_data(fft_a2a.data(),poly_size);
    }

    assert(fft_a.size() == poly_size);
    assert(fft_a2.size() == poly_size);
    assert(fft_a2a.size() == poly_size);

    
    

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
    std::vector<uint8_t> err_polys_A(poly_size);

    // coefficients associated with each error vector
    
    auto err_poly_coeffs_A = beaver_triple_seed.err_poly_coeffs_A;

    // positions of the T errors in each error vector
   
    
    auto err_poly_positions_A = beaver_triple_seed.err_poly_positions_A;

    std::vector<uint8_t> s_polys_A(poly_size);

    // coefficients associated with each error vector
    
    auto s_poly_coeffs_A = beaver_triple_seed.s_poly_coeffs_A;

    // positions of the T errors in each error vector
   
    
    auto s_poly_positions_A = beaver_triple_seed.s_poly_positions_A;

    std::vector<uint8_t> err_polys_A2(poly_size);

    // coefficients associated with each error vector
    
    std::vector<uint8_t> err_poly_coeffs_A2(t);

    // positions of the T errors in each error vector
   
    std::vector<size_t> err_poly_positions_A2(t);

    std::vector<uint8_t> s_polys_A2(poly_size);

    // coefficients associated with each error vector
    
    std::vector<uint8_t> s_poly_coeffs_A2(t);

    // positions of the T errors in each error vector
   
    std::vector<size_t> s_poly_positions_A2(t);

    get_full_poly(n,t,block_size,err_polys_A,err_poly_coeffs_A, err_poly_positions_A);
    get_full_poly(n,t,block_size,s_polys_A,s_poly_coeffs_A, s_poly_positions_A);

    if(party!=1){
        squre_poly(n,t,block_size,err_poly_coeffs_A,err_poly_positions_A,err_polys_A2,err_poly_coeffs_A2, err_poly_positions_A2);
        squre_poly(n,t,block_size,s_poly_coeffs_A,s_poly_positions_A,s_polys_A2,s_poly_coeffs_A2, s_poly_positions_A2);
    }
    
    std::vector<std::vector<uint8_t>> shareA(8, std::vector<uint8_t>(t*t));
    
    if(party==1){
        sender_OT(io,s_poly_coeffs_A,1,t,shareA[0]);
        sender_OT(io,err_poly_coeffs_A,1,t,shareA[1]);
        sender_OT(io,s_poly_coeffs_A,1,t,shareA[2]);
        sender_OT(io,err_poly_coeffs_A,1,t,shareA[3]);
        sender_OT(io,s_poly_coeffs_A,1,t,shareA[4]);
        sender_OT(io,err_poly_coeffs_A,1,t,shareA[5]);
        sender_OT(io,s_poly_coeffs_A,1,t,shareA[6]);
        sender_OT(io,err_poly_coeffs_A,1,t,shareA[7]);
       
    }
    else{
        receiver_OT(io,s_poly_coeffs_A,1,t,shareA[0]);
        receiver_OT(io,s_poly_coeffs_A,1,t,shareA[1]);
        receiver_OT(io,err_poly_coeffs_A,1,t,shareA[2]);
        receiver_OT(io,err_poly_coeffs_A,1,t,shareA[3]);
        receiver_OT(io,s_poly_coeffs_A2,1,t,shareA[4]);
        receiver_OT(io,s_poly_coeffs_A2,1,t,shareA[5]);
        receiver_OT(io,err_poly_coeffs_A2,1,t,shareA[6]);
        receiver_OT(io,err_poly_coeffs_A2,1,t,shareA[7]);
    }


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
    std::vector<uint8_t> fft_eA(s_polys_A.data(), s_polys_A.data() + poly_size);

    fft_recursive_uint8(fft_eA.data(), n, poly_size / 3);

    

    // Compute the coordinate-wise multiplication over the packed FFT result
    
    std::vector<uint8_t> res_poly_e(poly_size,0);
    std::vector<uint8_t> res_poly_s(poly_size,0);
    multiply_fft_8(fft_a.data(), fft_eA.data(), res_poly_e.data(), poly_size); // a*eA

    
    memcpy(fft_eA.data(), err_polys_A.data(), poly_size);

    fft_recursive_uint8(fft_eA.data(), n, poly_size / 3);

    multiply_fft_8(fft_a2a.data(), fft_eA.data(), res_poly_s.data(), poly_size); // a*eA

    // XOR the result into the accumulator.
    // Specifically, we XOR all the columns of the FFT result to get a
    // vector of size poly_size.
    
    for (size_t i = 0; i < poly_size; i++)
    {
        x_poly_A[i] = (res_poly_e[i]^res_poly_s[i])& 0b11;
        
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
    
    std::vector<std::vector<DPFParty>>dpf_partyA(8,std::vector<DPFParty>(t*t));
    
    if(party==1){
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[0], shareA[0], 
        s_poly_positions_A);
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[1], shareA[1], 
        err_poly_positions_A);
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[2], shareA[2], 
        s_poly_positions_A);
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[3], shareA[3], 
        err_poly_positions_A);
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[4], shareA[4], 
        s_poly_positions_A);
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[5], shareA[5], 
        err_poly_positions_A);
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[6], shareA[6], 
        s_poly_positions_A);
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[7], shareA[7], 
        err_poly_positions_A);
    }
    else{
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[0], shareA[0], 
            s_poly_positions_A);
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[1], shareA[1], 
        s_poly_positions_A);
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[2], shareA[2], 
        err_poly_positions_A);
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[3], shareA[3], 
        err_poly_positions_A);
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[4], shareA[4], 
        s_poly_positions_A2);
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[5], shareA[5], 
        s_poly_positions_A2);
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[6], shareA[6], 
        err_poly_positions_A2);
        KeyGen(prf_keys,io,party,n,1,t, dpf_domain_bits, 
            dpf_partyA[7], shareA[7], 
        err_poly_positions_A2);

    }
    
    
    if(party==1){
        printf("[....   ]Done with Step 3 (sampling DPF keys)\n");
        time_taken = ((double)clock() - time) / (CLOCKS_PER_SEC / 1000.0); // ms
        printf("Step 3 time elapsed %f ms\n", time_taken);
        time = clock();
    }

    // Allocate memory for the concatenated DPF outputs
    const size_t packed_block_size = ceil(block_size / 64.0);
    const size_t packed_poly_size = t * packed_block_size;

    // printf("[DEBUG]: packed_block_size = %zu\n", packed_block_size);
    // printf("[DEBUG]: packed_poly_size = %zu\n", packed_poly_size);

    // std::vector<uint128_t>packed_polys_A(packed_poly_size, 0);
    std::vector<std::vector<uint128_t>>packed_polys_A(8,std::vector<uint128_t>(packed_poly_size, 0));

    // Allocate memory for the output FFT
    // std::vector<uint32_t>fft_uA(poly_size, 0);
    std::vector<std::vector<uint8_t>>fft_uA(8,std::vector<uint8_t>(poly_size, 0));

    // Allocate memory for the final inner product
    
    
    std::vector<std::vector<uint8_t>>res_poly_mat_A(8,std::vector<uint8_t>(poly_size, 0));
    

    for(uint8_t i=0;i<8;i++){
        Fullevaluation(dpf_partyA[i], 1, t, packed_block_size, packed_poly_size,block_size, packed_polys_A[i]);
    }
    
    if(party==1){

        printf("[.....  ]Done with Step 4 (evaluating all DPFs)\n");
        time_taken = ((double)clock() - time) / (CLOCKS_PER_SEC / 1000.0); // ms
        printf("Step 4 time elapsed %f ms\n", time_taken);
        time = clock();
    }

    // TODO[optimization]: for arbitrary values of C, we only need to perform
    // C*(C+1)/2 FFTs which can lead to a more efficient implementation.
    // Because we assume C=4, we have C*C = 16 which fits perfectly into a
    // uint32 packing.
    for(uint8_t j=0;j<8;j++){
        auto polyA = packed_polys_A[j];
        
        size_t block_idx = 0;
        size_t bit_idx = 0;
        for (size_t i = 0; i < poly_size; i++)
        {
            if (i % block_size == 0 && i != 0)
            {
                block_idx++;
                bit_idx = 0;
            }
            const size_t packed_idx = block_idx * packed_block_size + floor(bit_idx / 64.0);
            const size_t packed_bit = (63 - bit_idx % 64);
            
            // Extract the i-th (packed) coefficient of the (j,k)-th polynomial
            const uint128_t packedA = polyA[packed_idx];
            

            // Extract the i-th coefficient from the packed coefficient
            const uint8_t coeffA =((packedA >> (2 * packed_bit)) & 0b11);
            
            // Pack the extracted coefficient into the ctr-th FFT slot
        
            fft_uA[j][i] = coeffA;
            

            bit_idx++;
        }
    
        fft_recursive_uint8(fft_uA[j].data(), n, poly_size / 3);
    }

   
    if(party==1){
    printf("[...... ]Done with Step 5 (computing FFTs)\n");
    time_taken = ((double)clock() - time) / (CLOCKS_PER_SEC / 1000.0); // ms
    printf("Step 5 time elapsed %f ms\n", time_taken);
    time = clock();
    }

    //************************************************************************
    // Step 6: Compute shares of z = <axa, u>
    //************************************************************************
   
        
    multiply_fft_8(fft_a2.data(), fft_uA[0].data(), res_poly_mat_A[0].data(), poly_size);
    multiply_fft_8(fft_a.data(), fft_uA[1].data(), res_poly_mat_A[1].data(), poly_size);
    multiply_fft_8(fft_a.data(), fft_uA[2].data(), res_poly_mat_A[2].data(), poly_size);
    multiply_fft_8(fft_a2a.data(), fft_uA[3].data(), res_poly_mat_A[3].data(), poly_size);
    multiply_fft_8(fft_a2a.data(), fft_uA[4].data(), res_poly_mat_A[4].data(), poly_size);
    multiply_fft_8(fft_a2.data(), fft_uA[5].data(), res_poly_mat_A[5].data(), poly_size);
    multiply_fft_8(fft_a.data(), fft_uA[6].data(), res_poly_mat_A[6].data(), poly_size);
    multiply_fft_8(fft_a2a.data(), fft_uA[7].data(), res_poly_mat_A[7].data(), poly_size);

    

    // XOR the (packed) columns into the accumulator.
    // Specifically, we perform column-wise XORs to get the result.
    for(size_t w = 0; w < 8; w++)
    {
        auto res_poly = res_poly_mat_A[w];
        for (size_t i = 0; i < poly_size; i++)
        {
        z_poly_A[i] ^= (res_poly[i]& 0b11);
        }
       
    }
    
    

    

    // Apply Trace Function
    for (size_t i = 0; i < poly_size; i++)
    {
        z_poly_A[i] = (z_poly_A[i]^mult_f4(z_poly_A[i],z_poly_A[i])) & 0b01;
    }

    if(party==1){
        printf("[.......]Done with Step 6 (recovering shares)\n\n");
        time_taken = ((double)clock() - time) / (CLOCKS_PER_SEC / 1000.0); // ms
        printf("Step 6 time elapsed %f ms\n", time_taken);

        
        time_taken = ((double)clock() - start) / (CLOCKS_PER_SEC / 1000.0); // ms
        printf("Time elapsed %f ms\n", time_taken);
    }
    if (party==1){
        int nums = 0;
        std::vector<uint8_t> z_poly_B(poly_size);
        std::vector<uint8_t> x_poly_B(poly_size);
        io->recv_data(z_poly_B.data(),poly_size);
        io->recv_data(x_poly_B.data(),poly_size);
        for (size_t i = 0; i < poly_size; i++)
        {
            const uint8_t res = z_poly_A[i] ^ z_poly_B[i];
            const uint8_t exp = x_poly_A[i]&x_poly_B[i];
            
            // printf("[DEBUG]: Got: (%i,%i), Expected: (%i, %i)\n",
            //        (res >> 1) & 1, res & 1, (exp >> 1) & 1, exp & 1);
            if (res != exp)
            {
                nums++;
                // printf("[DEBUG]: i=%zu, res=%i, exp=%i\n", i, res, exp);
                // printf("FAIL: Incorrect correlation output at index %zu\n", i);
                // printf("Got: (%i,%i), Expected: (%i, %i)\n",
                    // (res >> 1) & 1, res & 1, (exp >> 1) & 1, exp & 1);
                // exit(0);
            }
            // if(i%200==0) printf("Got: (%i,%i), Expected: (%i, %i)\n",z_poly_A[i] & 1, z_poly_B[i] & 1, x_poly_A[i] & 1, x_poly_B[i] & 1);
        }
        std::cout<<"出错的次数为："<<nums<<std::endl;

    }
    else{
        io->send_data(z_poly_A.data(),poly_size);
        io->flush();
        io->send_data(x_poly_A.data(),poly_size);
        io->flush();
        
    }
    delete io;
    DestroyPRFKey(prf_keys);
}

void get_full_poly(const size_t n, const size_t t,const size_t block_size,std::vector<uint8_t> &err_polys_A,
    std::vector<uint8_t> &err_poly_coeffs_A, std::vector<size_t> &err_poly_positions_A){
    
    for (size_t j = 0; j < t; j++)
    {
        
        // random *non-zero* coefficients in F4
        uint8_t a = err_poly_coeffs_A[j];
        // random index within the block
        size_t pos_A =  err_poly_positions_A[j];

        if (pos_A >= block_size)
        {
            printf("FAIL: position > block_size: %zu\n", pos_A);
            exit(0);
        }
        // set the coefficient at the error position to the error value
        err_polys_A[j * block_size + pos_A] = a;

    }
}
void sample_local_err(const size_t n, const size_t t,const size_t block_size,
    std::vector<uint8_t> &err_poly_coeffs_A, std::vector<size_t> &err_poly_positions_A){
    for (size_t j = 0; j < t; j++)
    {
        // random *non-zero* coefficients in F4
       
        
        uint8_t a = rand_f4x();
        
        err_poly_coeffs_A[j] = a;

        // random index within the block
        size_t pos_A = random_index(block_size - 1);
        
        err_poly_positions_A[j] = pos_A; 

        if (pos_A >= block_size)
        {
            printf("FAIL: position > block_size: %zu\n", pos_A);
            exit(0);
        }
        
    }
}
void squre_poly(const size_t n, const size_t t,const size_t block_size,
    std::vector<uint8_t> &err_poly_coeffs_A, std::vector<size_t> &err_poly_positions_A,
    std::vector<uint8_t> &err_polys_A2,std::vector<uint8_t> &err_poly_coeffs_A2, std::vector<size_t> &err_poly_positions_A2)
{
    std::vector<uint8_t> tritsA(log_base(t,3));
    std::vector<uint8_t> tritsB(n-log_base(t,3));
    for (size_t j = 0; j < t; j++)
    {
        

        // random *non-zero* coefficients in F4
        
        
        uint8_t a = err_poly_coeffs_A[j];
        

        // random index within the block
        size_t pos_A = err_poly_positions_A[j];
        
        int_to_trits(j, tritsA.data(), tritsA.size());
        int_to_trits(pos_A, tritsB.data(), tritsB.size());
        for(size_t k = 0; k < tritsA.size(); k++)
        {
            tritsA[k] = (tritsA[k] + tritsA[k]) % 3;
        }
        for(size_t k = 0; k < tritsB.size(); k++)
        {
            tritsB[k] = (tritsB[k] + tritsB[k]) % 3;
        }

        pos_A = trits_to_int(tritsB.data(), tritsB.size());
        size_t idx = trits_to_int(tritsA.data(), tritsA.size());

        if (pos_A >= block_size)
        {
            printf("FAIL: position > block_size: %zu\n", pos_A);
            exit(0);
        }

        
        

        // set the coefficient at the error position to the error value

        a = mult_f4(a,a);

        err_poly_positions_A2[idx] = pos_A;

        err_poly_coeffs_A2[idx] = a;

        err_polys_A2[idx * block_size + pos_A] = a;
    }
}

void Beaver_Triple_test(const size_t n,const size_t t ,int party,int port)
{
    const size_t poly_size = ipow(3, n);
    const size_t block_size = poly_size / t;
    const size_t dpf_domain_bits = n-log_base(t,3);
    Beaver_Triple_Seed beaver_triple_seed(t);
    sample_local_err(n,t,block_size,
        beaver_triple_seed.err_poly_coeffs_A,beaver_triple_seed.err_poly_positions_A);
    sample_local_err(n,t,block_size,
        beaver_triple_seed.s_poly_coeffs_A,beaver_triple_seed.s_poly_positions_A);
    Beaver_Triple_Seed beaver_triple_seed2(t);
    sample_local_err(n,t,block_size,
        beaver_triple_seed2.err_poly_coeffs_A,beaver_triple_seed2.err_poly_positions_A);
    sample_local_err(n,t,block_size,
        beaver_triple_seed2.s_poly_coeffs_A,beaver_triple_seed2.s_poly_positions_A);
    
    std::vector<uint8_t> fft_a(poly_size);
    std::vector<uint8_t> fft_a2(poly_size);
    std::vector<uint8_t> fft_a2a(poly_size);
    if(party==1)
        sample_a_and_tensor(fft_a.data(), fft_a2.data(), fft_a2a.data(), poly_size);
    std::vector<uint8_t> x_poly_A(poly_size,0);
    std::vector<uint8_t> z_poly_A(poly_size,0);
    std::vector<uint8_t> x_poly_A2(poly_size,0);
    std::vector<uint8_t> z_poly_A2(poly_size,0);
    test_pcg_F2(n,t,party,port,fft_a,fft_a2,fft_a2a,
        beaver_triple_seed,z_poly_A,x_poly_A);
    test_pcg_F2(n,t,party,port,fft_a,fft_a2,fft_a2a,
        beaver_triple_seed2,z_poly_A2,x_poly_A2);
    for(size_t w = 0 ; w < poly_size ;w++){
        z_poly_A[w] = (z_poly_A[w]^z_poly_A2[w]^(x_poly_A2[w]&x_poly_A[w]))& 0b1;
    }
    NetIO* io = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port);
    if(party==1){
        io->send_data(z_poly_A.data(),poly_size);
        io->flush();
        io->send_data(x_poly_A.data(),poly_size);
        io->flush();
        io->send_data(x_poly_A2.data(),poly_size);
        io->flush();
    }
    else{
        
        std::vector<uint8_t> receive_x(poly_size);
        std::vector<uint8_t> receive_z(poly_size);
        std::vector<uint8_t> receive_x2(poly_size);
        io->recv_data(receive_z.data(),poly_size);
        io->recv_data(receive_x.data(),poly_size);
        io->recv_data(receive_x2.data(),poly_size);
        int nums = 0;
        for(size_t w = 0 ; w < poly_size ;w++){
            uint8_t res = (z_poly_A[w]^receive_z[w]);
            uint8_t exp = (x_poly_A[w]^receive_x2[w])&(x_poly_A2[w]^receive_x[w]);
            if(res!=exp){
                nums++;
                // printf("FAIL: Incorrect correlation output at index %zu\n", i);
                // printf("Got: (%i,%i), Expected: (%i, %i)\n",
                //        (res >> 1) & 1, res & 1, (exp >> 1) & 1, exp & 1);
            }
            
        }
        std::cout<<"出错的次数为："<<nums<<std::endl;
    }
    
}

