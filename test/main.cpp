#include"otdpf.h"
#include"testF4.h"
#include"bench.h"
#include"testF2.h"
void testOutputCorrectness(
    uint128_t *shares,
    size_t num_outputs,
    size_t secret_index,
    uint128_t *secret_msg,
    size_t msg_len)
{
    for (size_t i = 0; i < msg_len; i++)
    {
       
        uint128_t res = shares[secret_index * msg_len + i];

        std::cout<<static_cast<int>(secret_index * msg_len + i)<<"\t"<<static_cast<int>(secret_msg[i])<<"\n";
        std::cout<<static_cast<int>(secret_index)<<"\t"<<static_cast<int>(res)<<"\n";

        if (res != secret_msg[i])
        {
            
            printf("FAIL (wrong message)\n");
            exit(0);
        }
    }

    for (size_t i = 0; i < num_outputs; i++)
    {
        if (i == secret_index)
            continue;

        for (size_t j = 0; j < msg_len; j++)
        {
            uint128_t res = shares[i * msg_len + j];


            if (res != 0)
            {
                
                printf("FAIL (non-zero) %zu\n", i);

                exit(0);
            }
        }
    }
}
size_t randIndex(size_t MAXRANDINDEX = 81)
{
    uint8_t index;
    RAND_bytes((uint8_t *)&index, sizeof(uint8_t));
    return index % ((size_t)MAXRANDINDEX);
}

uint128_t randMsg()
{
    uint128_t msg;
    RAND_bytes((uint8_t *)&msg, sizeof(uint128_t));
    return msg;
}

// 三进制按位异或（定义：对应位相加后取模3）
size_t ternary_xor(size_t a, size_t b) {
    size_t result = 0;
    size_t power = 1; // 当前位的权值 (3^0, 3^1, ...)
    
    while (a > 0 || b > 0) {
        // 取出当前位的三进制值
        int digit_a = a % 3;
        int digit_b = b % 3;
        
        // 计算异或结果（相加后取模3）
        int xor_digit = (digit_a + digit_b) % 3;
        
        // 将结果累加到最终值
        result += xor_digit * power;
        
        // 更新权值和剩余数值
        power *= 3;
        a /= 3;
        b /= 3;
    }
    return result;
}

void test_DPF(int party, int port,const size_t size,const size_t msg_len){
    // 初始化网络
    NetIO* io = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port);

    
    size_t MAXRANDINDEX = pow(3,size);
    size_t secret_index = randIndex(MAXRANDINDEX);
    
    // sample a random message of size msg_len
    std::vector<uint128_t> secret_msg(msg_len) ;
    for (size_t i = 0; i < msg_len; i++)
        secret_msg[i] = randMsg();

    
    struct PRFKeys *prf_keys = new PRFKeys ;
    if(party==ALICE) PRFKeyGen(prf_keys);
    clock_t time;
    clock_t start = clock();
    // sample the DPF keys
    DPFParty dpf(prf_keys,size,secret_index,secret_msg,party);
    std::cout<<static_cast<int>(secret_index)<<"\n";
    dpf.generate(io);
    time = clock();
    double time_taken = ((double)time-start) / (CLOCKS_PER_SEC / 1000.0); // ms
    printf("Key generation time %f ms\n", time_taken);

    std::vector<uint128_t> shares;
    dpf.fulldomainevaluation(shares);
    std::vector<uint128_t> receive(shares.size());
    uint128_t len = shares.size()*sizeof(uint128_t);
    
    time_taken = ((double)clock()-time) / (CLOCKS_PER_SEC / 1000.0); // ms
    printf("Full Evaluation time %f ms\n", time_taken);

    
    time_taken = ((double)clock()-start) / (CLOCKS_PER_SEC / 1000.0); // ms
    printf("Eval time (total) %f ms\n", time_taken);
    printf("DONE\n\n");
    
    if(party==1){
        io->send_data(shares.data(),len);
        io->flush();
        io->recv_data(receive.data(),len);
        for (size_t k = 0 ; k < shares.size() ; ++k)
        {
            shares[k]=shares[k] ^ receive[k];
           
        }
        std::vector<uint128_t> receive_msg(msg_len) ;
        io->recv_data(receive_msg.data(),msg_len*16);
        for(size_t k = 0 ; k < msg_len ; ++k) {
            secret_msg[k]^=receive_msg[k];
            std::cout<<static_cast<int>(k)<<"\t"<<static_cast<int>(secret_msg[k])<<std::endl;
            
        }
        size_t index = 0;
        io->recv_data(&index,sizeof(size_t));
        secret_index = ternary_xor(secret_index,index);
        std::cout<<static_cast<int> (secret_index)<<std::endl;
        testOutputCorrectness(shares.data(),MAXRANDINDEX,secret_index,secret_msg.data(),msg_len);
        std::cout<<"test pass\n";

    }
    else{
        io->recv_data(receive.data(),len);
        io->send_data(shares.data(),len);
        io->flush();
        io->send_data(secret_msg.data(),msg_len*16);
        io->flush();
        io->send_data(&secret_index,sizeof(size_t));
        io->flush();
    }
}

int main(int argc , char** argv){
    int party, port;
    size_t c = 4;
    size_t t = 9;
    size_t n = 6;
    size_t msg_len = 1;
    parse_party_and_port(argv, &party, &port);
    c = atoi(argv[3]);
    t = atoi(argv[4]);
    n = atoi(argv[5]);
    msg_len = atoi(argv[6]);
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--bench") == 0) {
            bench_pcg(n,c,t,party,port);
        } 
        else if (strcmp(argv[i], "--test_F4") == 0) {
            test_pcg(n,c,t,party,port);
        } else if (strcmp(argv[i], "--test_dpf")==0) {
            test_DPF(party,port,n,msg_len);
        }else if (strcmp(argv[i], "--test_F2")==0) {
            Beaver_Triple_test(n,t,party,port);
        }
    }

    return 0;
}

