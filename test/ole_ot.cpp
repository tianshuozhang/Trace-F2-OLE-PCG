#include <iostream>
#include <vector>
#include <random>
#include <cstdint>
#include <thread>
#include <mutex>
#include <emp-tool/emp-tool.h>
#include <emp-ot/emp-ot.h>
typedef unsigned __int128 uint128_t;
using namespace std;
using namespace emp;

// 互斥锁，防止多线程 `cout` 混乱
mutex mtx;

// GF(4) 乘法表
const uint8_t GF4_MULT_TABLE[4][4] = {
    {0, 0, 0, 0},
    {0, 1, 2, 3},
    {0, 2, 3, 1},
    {0, 3, 1, 2}
};

uint8_t gf4_multiply(const uint8_t a, const uint8_t b) {
    return GF4_MULT_TABLE[a][b];
}

// 生成随机 GF(4) 元素
uint8_t random_GF4_element() {
    static random_device rd;
    static mt19937 gen(rd());
    static uniform_int_distribution<int> dist(0, 3);
    return static_cast<uint8_t>(dist(gen));
}

// OT 发送方 (P1)
void sender_OT(NetIO* io, vector<uint8_t>& A, int c, int t, vector<uint8_t>& R) {
    const uint32_t length = t * c;
    const uint32_t size  = length*length;
    block* b0 = new block[size<<1]();
    block* b1 = new block[size<<1]();

    for (int jA = 0; jA < c; ++jA) {
        for (int jB = 0; jB < c ; ++jB) {
            for(int iA = 0 ; iA < t ; ++iA){
                uint8_t a = A[jA * t + iA];
                for(int iB = 0 ; iB < t ; ++iB){
                    uint32_t offset = (jA*c+jB)*t*t+iA*t+iB;
                    uint8_t r = random_GF4_element();
                    R[offset] = r;  // 存储 R 值
                    b0[offset] = makeBlock(0, r);
                    b1[offset] = makeBlock(0, r ^ gf4_multiply(a, 1));
                    b0[offset + size] = makeBlock(0, r ^ gf4_multiply(a, 2));
                    b1[offset + size] = makeBlock(0, r ^ gf4_multiply(a, 3));
                }
            }
           
        }
    }
    
           
            

            

    // 执行OT发送
    OTNEW<NetIO> ot(io);
    ot.send(b0, b1, size);

    delete[] b0;
    delete[] b1;

    lock_guard<mutex> lock(mtx);
    cout << "Sender: OT 发送完成" << endl;
}

// OT 接收方 (P2)
void receiver_OT(NetIO* io, vector<uint8_t>& B, int c, int t, vector<uint8_t>& S) {
    uint32_t length = t * c;
    const uint32_t size  = length*length;
    block* r = new block[size];
    bool* b = new bool[size<<1];

    for (int jA = 0; jA < c; ++jA) {
        for (int jB = 0; jB < c ; ++jB) {
            for(int iA = 0 ; iA < t ; ++iA){
                for(int iB = 0 ; iB < t ; ++iB){
                    uint8_t choice = B[jB*t+iB];
                    uint32_t offset = (jA*c+jB)*t*t+iA*t+iB;
                    b[offset] = choice & 1;
                    b[offset + size] = (choice >> 1) & 1;
                }
            }
           
        }
    }
    OTNEW<NetIO> ot(io);
    ot.recv(r, b, length * length);

    for (int jA = 0; jA < c; ++jA) {
        for (int jB = 0; jB < c ; ++jB) {
            for(int iA = 0 ; iA < t ; ++iA){
                for(int iB = 0 ; iB < t ; ++iB){
                   uint32_t offset = (jA*c+jB)*t*t+iA*t+iB;
                   S[offset] = (uint8_t)_mm_extract_epi64(r[offset], 0);
                }
            }
           
        }
    }
    delete[] r;
    delete[] b;

    lock_guard<mutex> lock(mtx);
    cout << "Receiver: OT 接收完成" << endl;
}

// 线程函数
void run_sender(int c, int t,vector<uint8_t>&A,vector<uint8_t> &R) {
    NetIO io(nullptr, 12345);  // 发送方监听端口
    

    for (int i = 0; i < c * t; ++i) A[i] = random_GF4_element();
    sender_OT(&io, A, c, t, R);

    lock_guard<mutex> lock(mtx);
    cout << "Sender A values: ";
    for (auto a : A) cout << (int)a << " ";
    cout << endl;

    
}

void run_receiver(int c, int t,vector<uint8_t>&B,vector<uint8_t>&S) {
    NetIO io("127.0.0.1", 12345);  // 接收方连接服务器
    

    for (int i = 0; i < c * t; ++i) B[i] = random_GF4_element();
    receiver_OT(&io, B, c, t, S);

    lock_guard<mutex> lock(mtx);
    cout << "Receiver B values: ";
    for (auto b : B) cout << (int)b << " ";
    cout << endl;

    
}

int main(int argc ,char**argv) {
    int c = 3, t = 2; // 示例参数
    if(argc==2) c = atoi(argv[1]);
    if(argc==3) c = atoi(argv[1]), t = atoi(argv[2]);

    
    vector<uint8_t> A(c * t);
    vector<uint8_t> R(c * c * t * t);
    vector<uint8_t> B(c * t);
    vector<uint8_t> S(c * c * t * t);
    // 启动发送方和接收方线程
    thread sender_thread(run_sender, c, t,ref(A),ref(R));
    thread receiver_thread(run_receiver, c, t,ref(B),ref(S));

    // 等待两个线程完成
    sender_thread.join();
    receiver_thread.join();
    cout<<"乘法运算结果：\n";
    for (int jA = 0; jA < c; ++jA) {
        for (int jB = 0; jB < c ; ++jB) {
            for(int iA = 0 ; iA < t ; ++iA){
                for(int iB = 0 ; iB < t ; ++iB){
                   cout<<(int)gf4_multiply(A[jA * t + iA], B[jB * t + iB])<<" ";
                }
            }
            cout<<endl;
           
        }
    }
    cout<<"加法share结果：\n";
    for (int jA = 0; jA < c; ++jA) {
        for (int jB = 0; jB < c ; ++jB) {
            for(int iA = 0 ; iA < t ; ++iA){
                for(int iB = 0 ; iB < t ; ++iB){
                    uint32_t offset = (jA*c+jB)*t*t+iA*t+iB;
                    cout<<(int)(R[offset]^S[offset])<<" ";
                }
            }
            cout<<endl;
           
        }
    }
    cout<<"OLE结果：\n";
    for (int jA = 0; jA < c; ++jA) {
        for (int jB = 0; jB < c ; ++jB) {
            for(int iA = 0 ; iA < t ; ++iA){
                for(int iB = 0 ; iB < t ; ++iB){
                    uint32_t offset = (jA*c+jB)*t*t+iA*t+iB;
                    cout<<((R[offset]^S[offset])==gf4_multiply(A[jA * t + iA], B[jB * t + iB]))<<" ";
                }
            }
            cout<<endl;
           
        }
    }

    return 0;
}