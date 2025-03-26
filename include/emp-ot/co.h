#ifndef EMP_OTCO_H__
#define EMP_OTCO_H__
#include <emp-tool/emp-tool.h>
#include "emp-ot/ot.h"
namespace emp {

template<typename IO>
/*
 * Chou Orlandi OT
 * [REF] Implementation of "The Simplest Protocol for Oblivious Transfer"
 * https://eprint.iacr.org/2015/267.pdf
 */

class OTCO: public OT<IO> { public:
	IO* io;
	Group *G = nullptr;
	bool delete_G = true;
	OTCO(IO* io, Group * _G = nullptr) {
		this->io = io;
		if (_G == nullptr)
			G = new Group();
		else {
			G = _G;
			delete_G = false;
		}
	}
	~OTCO() {
		if (delete_G)
			delete G;
	}

	void send(const block* data0, const block* data1, int64_t length) override {
		BigInt a;
		Point A, AaInv;
		block res[2];
		Point * B = new Point[length];
		Point * BA = new Point[length];

		G->get_rand_bn(a);
		A = G->mul_gen(a);
		io->send_pt(&A);
		AaInv = A.mul(a);
		AaInv = AaInv.inv();

		for(int64_t i = 0; i < length; ++i) {
			io->recv_pt(G, &B[i]);
			B[i] = B[i].mul(a);
			BA[i] = B[i].add(AaInv);
		}
		io->flush();

		for(int64_t i = 0; i < length; ++i) {
			res[0] = Hash::KDF(B[i], i) ^ data0[i];
			res[1] = Hash::KDF(BA[i], i) ^ data1[i];
			io->send_data(res, 2*sizeof(block));
		}

		delete[] BA;
		delete[] B;
	}

	void recv(block* data, const bool* b, int64_t length) override {
		BigInt * bb = new BigInt[length];
		Point * B = new Point[length],
				* As = new Point[length],
				A;

		for(int64_t i = 0; i < length; ++i)
			G->get_rand_bn(bb[i]);

		io->recv_pt(G, &A);

		for(int64_t i = 0; i < length; ++i) {
			B[i] = G->mul_gen(bb[i]);
			if (b[i]) 
				B[i] = B[i].add(A);
			io->send_pt(&B[i]);
		}
		io->flush();

		for(int64_t i = 0; i < length; ++i)
			As[i] = A.mul(bb[i]);

		block res[2];
		for(int64_t i = 0; i < length; ++i) {
			io->recv_data(res, 2*sizeof(block));
			data[i] = Hash::KDF(As[i], i);
			if(b[i])
				data[i] = data[i] ^ res[1];
			else
				data[i] = data[i] ^ res[0];
		}
		
		delete[] bb;
		delete[] B;
		delete[] As;
	}
};

template<typename IO>
/*
 * Chou Orlandi OT
 * [REF] Implementation of "The Simplest Protocol for Oblivious Transfer"
 * https://eprint.iacr.org/2015/267.pdf
 */

class OTNEW: public OT<IO> { public:
	IO* io;
	Group *G = nullptr;
	bool delete_G = true;
	OTNEW(IO* io, Group * _G = nullptr) {
		this->io = io;
		if (_G == nullptr)
			G = new Group();
		else {
			G = _G;
			delete_G = false;
		}
	}
	~OTNEW() {
		if (delete_G)
			delete G;
	}

	void send(const block* data0, const block* data1,int64_t length) override {
		
		BigInt a;
		Point A, AaInv;
        BigInt a1;
		Point A1, AaInv1;
		block res[4];
		Point * B = new Point[length];
		Point * BA = new Point[length];
        Point * B1 = new Point[length];
		Point * BA1 = new Point[length];

		G->get_rand_bn(a);
		A = G->mul_gen(a);
		io->send_pt(&A);
		AaInv = A.mul(a);
		AaInv = AaInv.inv();

		for(int64_t i = 0; i < length; ++i) {
			io->recv_pt(G, &B[i]);
			B[i] = B[i].mul(a);
			BA[i] = B[i].add(AaInv);
		}
		io->flush();

        G->get_rand_bn(a1);
		A1 = G->mul_gen(a1);
		io->send_pt(&A1);
		AaInv1 = A1.mul(a1);
		AaInv1 = AaInv1.inv();

        for(int64_t i = 0; i < length; ++i) {
			io->recv_pt(G, &B1[i]);
			B1[i] = B1[i].mul(a1);
			BA1[i] = B1[i].add(AaInv1);
		}
		io->flush();

		for(int64_t i = 0; i < length; ++i) {
			res[0] = Hash::KDF(B[i], i) ^Hash::KDF(B1[i], i)^ data0[i];
			res[1] = Hash::KDF(BA[i], i)^ Hash::KDF(B1[i], i)^ data1[i];
            res[2] = Hash::KDF(B[i], i)^Hash::KDF(BA1[i], i)^ data0[length+i];
			res[3] = Hash::KDF(BA[i], i)^Hash::KDF(BA1[i], i) ^ data1[length+i];
			io->send_data(res, 4*sizeof(block));
		}

		delete[] BA;
		delete[] B;
        delete[] BA1;
		delete[] B1;
	}

	void recv(block* data, const bool* b, int64_t length) override {
		BigInt * bb = new BigInt[length];
		Point * B = new Point[length],
				* As = new Point[length],
				A;
        BigInt * bb1 = new BigInt[length];
		Point * B1 = new Point[length],
				* As1 = new Point[length],
				A1;


		for(int64_t i = 0; i < length; ++i)
			G->get_rand_bn(bb[i]);

		io->recv_pt(G, &A);

		for(int64_t i = 0; i < length; ++i) {
			B[i] = G->mul_gen(bb[i]);
			if (b[i]) 
				B[i] = B[i].add(A);
			io->send_pt(&B[i]);
		}
		io->flush();

        for(int64_t i = 0; i < length; ++i)
			G->get_rand_bn(bb1[i]);

		io->recv_pt(G, &A1);

		for(int64_t i = 0; i < length; ++i) {
			B1[i] = G->mul_gen(bb1[i]);
			if (b[length+i]) 
				B1[i] = B1[i].add(A1);
			io->send_pt(&B1[i]);
		}
		io->flush();

		for(int64_t i = 0; i < length; ++i)
			As[i] = A.mul(bb[i]);

        for(int64_t i = 0; i < length; ++i)
			As1[i] = A1.mul(bb1[i]);

		block res[4];
		for(int64_t i = 0; i < length; ++i) {
			io->recv_data(res, 4*sizeof(block));
			data[i] = Hash::KDF(As[i], i)^Hash::KDF(As1[i], i);
			short ind = b[length+i] ? 1 : 0;
            ind = (ind<<1) + (b[i] ? 1 : 0) ;
            data[i] = data[i] ^ res[ind];
		}
		
		delete[] bb;
		delete[] B;
		delete[] As;
        delete[] bb1;
		delete[] B1;
		delete[] As1;
	}
};


}//namespace
#endif// OT_CO_H__
