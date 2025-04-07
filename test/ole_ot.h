#ifndef OLE_OT
#define OLE_OT
#include <emp-tool/emp-tool.h>
#include <emp-ot/emp-ot.h>
typedef unsigned __int128 uint128_t;

void receiver_OT(NetIO* io, vector<uint8_t>& B, const size_t c,const size_t t, vector<uint8_t>& S);
void sender_OT(NetIO* io, vector<uint8_t>& A, const size_t c,const size_t t, vector<uint8_t>& R);
#endif