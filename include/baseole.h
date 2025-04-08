#ifndef _BASEOLE
#define _BASEOLE
#include <emp-tool/emp-tool.h>
#include <emp-ot/emp-ot.h>

void receiver_OT(NetIO* io, vector<uint8_t>& B, const size_t c,const size_t t, vector<uint8_t>& S);
void sender_OT(NetIO* io, vector<uint8_t>& A, const size_t c,const size_t t, vector<uint8_t>& R);
#endif