/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include <iomanip>
#include "EncGWAS.h"
#include <StringUtils.h>
#include <NTL/BasicThreadPool.h>

EncGWAS::EncGWAS(Scheme& scheme, SecretKey& secretKey) : scheme(scheme), secretKey(secretKey) {
}

void EncGWAS::addGWASKeys(long n2) {

	scheme.addLeftRotKeys(secretKey);
	scheme.addLeftRotKey(secretKey, 3);
	scheme.addLeftRotKey(secretKey, 2 * n2 - 8);
	scheme.addLeftRotKey(secretKey, n2 - 4);
	scheme.addLeftRotKey(secretKey, (n2 - 4) * 2);
	scheme.addLeftRotKey(secretKey, (n2 - 4) * 3);
	scheme.addLeftRotKey(secretKey, 12);
	scheme.addLeftRotKey(secretKey, n2);
	scheme.addLeftRotKey(secretKey, 2 * n2);
	scheme.addLeftRotKey(secretKey, 2 * n2 * 2);
	scheme.addLeftRotKey(secretKey, 2 * n2 * 3);
	scheme.addLeftRotKey(secretKey, (2 * n2 - 1) * 1);
	scheme.addLeftRotKey(secretKey, (2 * n2 - 1) * 2);
	scheme.addLeftRotKey(secretKey, (2 * n2 - 1) * 3);
	scheme.addLeftRotKey(secretKey, (2 * n2 - 4) * 1);
	scheme.addLeftRotKey(secretKey, (2 * n2 - 4) * 2);
	scheme.addLeftRotKey(secretKey, (2 * n2 - 4) * 3);

	// long logsl_thread = (log2(n2) - 2) / 3;
	// for(long j = 0; j < logsl_thread; j++){
	// 	for(int i = 0; i < 7; i++){
	// 		scheme.addLeftRotKey(secretKey, (1 << (3 * j + 2)) * (i + 1));
	// 	}
	// }
	// scheme.addLeftRotKey(secretKey, 2 * n2 + 1);
	// scheme.addLeftRotKey(secretKey, 2 * n2 + 2);
	// scheme.addLeftRotKey(secretKey, 2 * n2 + 3);
	// scheme.addLeftRotKey(secretKey, 4 * n2 + 1);
	// scheme.addLeftRotKey(secretKey, 4 * n2 + 2);
	// scheme.addLeftRotKey(secretKey, 4 * n2 + 3);
	// scheme.addLeftRotKey(secretKey, 6 * n2 + 1);
	// scheme.addLeftRotKey(secretKey, 6 * n2 + 2);
	// scheme.addLeftRotKey(secretKey, 6 * n2 + 3);


	scheme.addRightRotKeys(secretKey);
	scheme.addRightRotKey(secretKey, 3);
	scheme.addRightRotKey(secretKey, n2 - 1);
	scheme.addRightRotKey(secretKey, n2 - 2);
	scheme.addRightRotKey(secretKey, n2 - 3);
	scheme.addRightRotKey(secretKey, 12);

	// logsl_thread = (log2(8 * n2) - 2) / 3;
	// for(long j = 0; j < logsl_thread; j++){
	// 	for(int i = 0; i < 7; i++){
	// 		scheme.addRightRotKey(secretKey, (1 << (3 * j + 2)) * (i + 1));
	// 	}
	// }

}

void EncGWAS::encryptX(Ciphertext& encX, double** X, long n2, long logp, long logq) {
	long slots = n2 << 3;
	double* tmp = new double[slots]();
	for (long j = 0; j < 4; ++j) {
		for (long i = 0; i < (n2 << 1); ++i) {
			tmp[i + (j * (n2 << 1))] = X[i % n2][j];
		}
		encX = scheme.encrypt(tmp, slots, logp, logq);
	}
	delete[] tmp;
}

void EncGWAS::encryptY(Ciphertext& encY, double* Y, long n2, long logp, long logq) {
	encY = scheme.encrypt(Y, n2, logp, logq);
}

void EncGWAS::encryptS(Ciphertext* encS, double** S, long n2, long m, long ScolNum, long sNum, long logp, long logq){
	double* tmp = new double[n2 * ScolNum]();
	for (long k = 0; k < sNum - 1; ++k) {
		for (long i = 0; i < n2; ++i) {
			for (long j = 0; j < ScolNum; ++j) {
				tmp[n2 * j + i] = S[i][ScolNum * k + j];
			}
		}
		encS[k] = scheme.encrypt(tmp, n2 * ScolNum, logp, logq);
		for (long i = 0; i < n2 * ScolNum; ++i) {
			tmp[i] = 0.0;
		}
	}
	long rest = m - ScolNum * (sNum - 1);
	for (long i = 0; i < n2; ++i) {
		for (long j = 0; j < rest; ++j) {
			tmp[n2 * j + i] = S[i][ScolNum * (sNum - 1) + j];
		}
	}
	encS[sNum - 1] = scheme.encrypt(tmp, n2 * ScolNum, logp, logq);

	delete[] tmp;
}

void EncGWAS::encryptDX(Ciphertext* encDX, double** X, long n2, long logp, long logq) {
	long slots = n2 << 2;
	double* tmp = new double[slots]();
	for (long j = 0; j < 4; ++j) {
		for (long i = 0; i < n2; ++i) tmp[i + (((4 + i - j) % 4) * n2)] = X[i][((4 + i - j) % 4)];
		encDX[j] = scheme.encrypt(tmp, slots, logp, logq);
		for (int i = 0; i < slots; ++i) tmp[i] = 0.0;
	}
	delete[] tmp;
}

void EncGWAS::expandDX(Ciphertext* encDXExp, Ciphertext* encDX, long n2) {
	for (long i = 0; i < 4; ++i) {
		encDXExp[i] = encDX[i];
		for (long j = 1; j < 4; j <<= 1) {
			Ciphertext rot = scheme.leftRotateFast(encDXExp[i], n2 * j);
			scheme.addAndEqual(encDXExp[i], rot);
			rot.kill();
		}
	}
}

void EncGWAS::expandDXnew(Ciphertext* encDXExp, Ciphertext* encDX, long n2) {
	for (long i = 0; i < 4; ++i) {
		encDXExp[i] = encDX[i];
		for (long j = 1; j < 4; j <<= 1) {
			// Ciphertext rot = scheme.leftRotateFast(encDXExp[i], n2 * j); // original
			Ciphertext rot = scheme.leftRotateFast(encDXExp[i], j);  //modified
			scheme.addAndEqual(encDXExp[i], rot);
			rot.kill();
		}
	}
}

ZZ* EncGWAS::maskAdjU(long n2, long logc) {
	ZZ* maskAdjU = new ZZ[scheme.ring.N];
	long slots = n2 << 3;
	double* vals = new double[slots]();
	for (long i = 0; i < 4; i++) {
		for (long j = 0; j < (n2 << 1); j++) vals[i * (n2 << 1) + j] = 2 * (i % 2 == j % 2) - 1;
	}
	scheme.ring.encode(maskAdjU, vals, slots, logc);
	delete[] vals;
	return maskAdjU;
}

ZZ* EncGWAS::maskDetU(long n2, long logc) {
	ZZ* maskDetU = new ZZ[scheme.ring.N];
	long slots = n2 << 3;
	double* vals = new double[slots]();
	for (long i = 0; i < 4; i++) {
		for (long j = 0; j < (n2 >> 1); j++) vals[i * (n2 << 1) + j] = 1.0;
		for (long j = (n2 >> 1); j < (n2 << 1); j++) vals[i * (n2 << 1) + j] = 0;
	}
	scheme.ring.encode(maskDetU, vals,  slots, logc);
	delete[] vals;
	return maskDetU;
}

ZZ** EncGWAS::maskBeta(long n2, long logc) {
	ZZ** maskBeta = new ZZ*[4];
	long slots = n2 << 3;
	double* vals = new double[slots]();
	for (long i = 0; i < 4; ++i) {
		maskBeta[i] = new ZZ[scheme.ring.N];
		for (long j = 0; j < (n2 << 3); j++) vals[j] = (j == (i * (n2 << 1)));
		scheme.ring.encode(maskBeta[i], vals, slots, logc);
		for (long j = 0; j < slots; ++j) vals[j] = 0.0;
	}
	delete[] vals;
	return maskBeta;
}

ZZ** EncGWAS::maskRU(long n2, long logc) {
	ZZ** maskRU = new ZZ*[4];
	long slots = n2 << 3;
	double* vals = new double[slots]();
	for (long i = 0; i < 4; ++i){
		maskRU[i] = new ZZ[scheme.ring.N];
		for (long j = 0; j < slots; ++j) vals[j] = 0.0;
		for (long j = 0; j < 4; ++j) vals[i * (n2 << 1) + j] = 1.0;
		scheme.ring.encode(maskRU[i], vals, slots, logc);
	}
	delete[] vals;
	return maskRU;
}

ZZ** EncGWAS::maskXWS(long n2, long logc) {
	ZZ** maskXWS = new ZZ*[6];
	double* vals = new double[n2]();
	for (long i = 0; i < 3; ++i) {
		maskXWS[i] = new ZZ[scheme.ring.N];
		for (long j = 0; j < n2 - i - 1; ++j) vals[j] = 1.0;
		scheme.ring.encode(maskXWS[i], vals, n2, logc);
		for (long j = 0; j < n2; ++j) vals[j] = 0.0;
	}
	for (long i = 0; i < 3; ++i) {
		maskXWS[i + 3] = new ZZ[scheme.ring.N];
		for (long j = n2 - i - 1; j < n2; ++j) vals[j] = 1.0;
		scheme.ring.encode(maskXWS[i + 3], vals, n2, logc);
		for (long j = 0; j < n2; ++j) vals[j] = 0.0;
	}
	delete[] vals;
	return maskXWS;
}

ZZ** EncGWAS::maskDU(long logc) {
	ZZ** maskDU = new ZZ*[4];
	double* vals = new double[16]();
	for (long i = 0; i < 4; ++i) {
		maskDU[i] = new ZZ[scheme.ring.N];
		for (long j = 0; j < 4; ++j) vals[(i + j) % 4 + j * 4] = 1.0;
		scheme.ring.encode(maskDU[i], vals, 16, logc);
		for (long j = 0; j < 16; ++j) vals[j] = 0.0;
	}
	delete[] vals;
	return maskDU;
}

ZZ** EncGWAS::maskDVUV(long logc) {
	ZZ** maskDVUV = new ZZ*[6];
	double* vals = new double[4]();
	for (long i = 0; i < 3; ++i) {
		maskDVUV[i] = new ZZ[scheme.ring.N];
		vals[i] = 1.0;
		scheme.ring.encode(maskDVUV[i], vals, 4, logc);
	}
	vals[3] = 1.0;
	for (long i = 0; i < 3; ++i) {
		maskDVUV[i + 3] = new ZZ[scheme.ring.N];
		vals[i] = 0.0;
		scheme.ring.encode(maskDVUV[i + 3], vals, 4, logc);
	}
	delete[] vals;
	return maskDVUV;
}

ZZ** EncGWAS::maskPack(long n2, long packNum, long logc) {
	ZZ** maskPack = new ZZ*[packNum];
	double * vals = new double[n2]();
	for (long i = 0; i < packNum; ++i) {
		maskPack[i] = new ZZ[scheme.ring.N];
		vals[4 * i] = 1.0; vals[4 * i + 1] = 1.0; vals[4 * i + 2] = 1.0; vals[4 * i + 3] = 1.0;
		scheme.ring.encode(maskPack[i], vals, n2, logc);
		vals[4 * i] = 0.0; vals[4 * i + 1] = 0.0; vals[4 * i + 2] = 0.0; vals[4 * i + 3] = 0.0;
	}
	delete[] vals;
	return maskPack;
}
ZZ* EncGWAS::maskIPS(long n2, long logc) {
	ZZ* maskIPS = new ZZ[scheme.ring.N];
	double* vals = new double[scheme.ring.N >> 1]();
	for(int i = 0; i < scheme.ring.N >> 1; i++){
		if(i % 4 == 0){
			vals[i] = 1;
		}
	}
	scheme.ring.encode(maskIPS, vals, scheme.ring.N >> 1, logc);
	return maskIPS;
}


void EncGWAS::evalUinit(Ciphertext& encU, Ciphertext& encX, Ciphertext* encDXExp, long n2, long logp) {
	Ciphertext encWX = scheme.divByPo2(encX, 2);
	Ciphertext* ctmp = new Ciphertext[4];

	NTL_EXEC_RANGE(4, first, last);
	for (long i = first; i < last; ++i) {
		ctmp[i] = scheme.modDownTo(encDXExp[i], encWX.logq);
		scheme.multAndEqual(ctmp[i], encWX);
		if(i != 0) {
			scheme.leftRotateFastAndEqual(ctmp[i], i);
		}
	}
	NTL_EXEC_RANGE_END;

	encU = ctmp[0];
	for (long i = 1; i < 4; ++i) {
		scheme.addAndEqual(encU, ctmp[i]);
	}

	delete[] ctmp;

	long logsl = log2(n2);
	//Normal version
	for (long j = 2; j < logsl; ++j) {
		Ciphertext rot = scheme.leftRotateFast(encU, (1 << j));
		scheme.addAndEqual(encU, rot);
	}

	//Multithreading version
	// Ciphertext* rot = new Ciphertext[7];
	// long logsl_thread = (logsl - 2) / 3;
	// for(long j = 0; j < logsl_thread; j++){
	// 	NTL_EXEC_RANGE(7, first, last);
	// 	for(long i = first; i < last; i++){
	// 		rot[i] = scheme.leftRotateFast(encU, (1 << (3 * j + 2)) * (i + 1));
	// 	}
	// 	NTL_EXEC_RANGE_END;
	// 	for(long i = 0; i < 7; i++){
	// 		scheme.addAndEqual(encU, rot[i]);
	// 	}
	// }
	// for(long j = 3 * logsl_thread; j < logsl - 2; j++){
	// 	Ciphertext rotat = scheme.leftRotateFast(encU, (1 << (j + 2)));
	// 	scheme.addAndEqual(encU, rotat);
	// }
	// delete[] rot;

	scheme.reScaleByAndEqual(encU, logp);
	encWX.kill();
}

void EncGWAS::evalU(Ciphertext& encU, Ciphertext& encX, Ciphertext& encW, Ciphertext* encDXExp, long n2, long logp) {
	Ciphertext encWX = scheme.modDownTo(encX, encW.logq);
	scheme.multAndEqual(encWX, encW);

	Ciphertext* ctmp = new Ciphertext[4];

	NTL_EXEC_RANGE(4, first, last);
	for (long i = first; i < last; ++i) {
		ctmp[i] = scheme.modDownTo(encDXExp[i], encWX.logq);
		scheme.multAndEqual(ctmp[i], encWX);
		if(i != 0) {
			scheme.leftRotateAndEqual(ctmp[i], i);
		}
	}
	NTL_EXEC_RANGE_END;

	encU = ctmp[0];
	for (long i = 1; i < 4; ++i) {
		scheme.addAndEqual(encU, ctmp[i]);
	}

	delete[] ctmp;

	long logsl = log2(n2);

	//Normal version
	for (long j = 2; j < logsl; ++j) {
		Ciphertext rot = scheme.leftRotateFast(encU, (1 << j));
		scheme.addAndEqual(encU, rot);
	}

	// //Multithreading version
	// Ciphertext* rot = new Ciphertext[7];
	// long logsl_thread = (logsl - 2) / 3;
	// for(long j = 0; j < logsl_thread; j++){
	// 	NTL_EXEC_RANGE(7, first, last);
	// 	for(long i = first; i < last; i++){
	// 		rot[i] = scheme.leftRotateFast(encU, (1 << (3 * j + 2)) * (i + 1));
	// 	}
	// 	NTL_EXEC_RANGE_END;
	// 	for(long i = 0; i < 7; i++){
	// 		scheme.addAndEqual(encU, rot[i]);
	// 	}
	// }
	// for(long j = 3 * logsl_thread; j < logsl - 2; j++){
	// 	Ciphertext rotat = scheme.leftRotateFast(encU, (1 << (j + 2)));
	// 	scheme.addAndEqual(encU, rotat);
	// }
	// delete[] rot;


	scheme.reScaleByAndEqual(encU, 2 * logp);
	encWX.kill();
}

void EncGWAS::evalAdjU(Ciphertext& encAdjU, Ciphertext& encU, ZZ* maskAdjU, long n2, long logp, long logc) {
	Ciphertext* ctmp = new Ciphertext[9];

	NTL_EXEC_RANGE(3, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext rot = scheme.leftRotateFast(encU, i + 1);
		for (long j = 0; j < 3; ++j) {
			ctmp[3 * i + j] = scheme.leftRotateFast(rot, 2 * n2 * (j + 1));
		}
		rot.kill();
	}
	NTL_EXEC_RANGE_END;
	// NTL_EXEC_RANGE(8, first, last);
	// for (long i = first; i < last; ++i) {
	// 	ctmp[i] = scheme.leftRotateFast(encU, 2 * n2 * ((i % 3) + 1) + ((i / 3) + 1));
	// }
	// NTL_EXEC_RANGE_END;
	// ctmp[8] = scheme.leftRotateFast(encU, 6 * n2 + 3); 


	NTL_EXEC_RANGE(3, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext diag1 = scheme.mult(ctmp[(3 * (i + 1) + 1) % 9], ctmp[(3 * (i + 2) + 2) % 9]);
		Ciphertext diag2 = scheme.mult(ctmp[(3 * (i + 2) + 1) % 9], ctmp[(3 * (i + 1) + 2) % 9]);
		scheme.subAndEqual(diag1, diag2);
		scheme.multAndEqual(ctmp[3 * i], diag1);
		diag1.kill();
		diag2.kill();
	}
	NTL_EXEC_RANGE_END;

	encAdjU = scheme.add(ctmp[0], ctmp[3]);
	scheme.addAndEqual(encAdjU, ctmp[6]);

	delete[] ctmp;

	scheme.multByPolyAndEqual(encAdjU, maskAdjU, logc);

	scheme.reScaleByAndEqual(encAdjU, 2 * logp + logc);
}

void EncGWAS::evalDetU(Ciphertext& encDetU, Ciphertext& encU, Ciphertext& encAdjU, ZZ* maskDetU, long n2, long logp, long logc) {
	encDetU = scheme.modDownTo(encU, encAdjU.logq);
	scheme.multAndEqual(encDetU, encAdjU);

	for (long i = 1; i <= 2; ++i) {
		Ciphertext rot = scheme.leftRotateFast(encDetU, 2 * n2 * i);
		scheme.addAndEqual(encDetU, rot);
		rot.kill();
	}

	scheme.multByPolyAndEqual(encDetU, maskDetU, logc);

	for (long i = 1; i <= 2; ++i){
		Ciphertext rot = scheme.rightRotateFast(encDetU, i * (n2 / 2));
		scheme.addAndEqual(encDetU, rot);
		rot.kill();
	}

	scheme.reScaleByAndEqual(encDetU, logp + logc);
}

void EncGWAS::rearrange(Ciphertext& encU, ZZ** maskRU, long n2, long logp, long logc) {
	Ciphertext* cvec = new Ciphertext[4];

	NTL_EXEC_RANGE(4, first, last);
	for (long i = first; i < last; ++i) {
		cvec[i] = scheme.multByPoly(encU, maskRU[i], logc);
		if (i != 0) {
			scheme.leftRotateFastAndEqual(cvec[i], (2 * n2 - 4) * i);
		}
	}
	NTL_EXEC_RANGE_END;

	encU = cvec[0];
	for (long i = 1; i < 4; ++i) {
		scheme.addAndEqual(encU, cvec[i]);
	}

	delete[] cvec;

	long logsl = log2(n2 * 8);
	for (long i = 4; i < logsl; ++i){
		Ciphertext rot = scheme.leftRotateFast(encU, (1 << i));
		scheme.addAndEqual(encU, rot);
		rot.kill();
	}

	scheme.reScaleByAndEqual(encU, logc);
}

void EncGWAS::evalXv(Ciphertext& encw, Ciphertext& encv, Ciphertext* encDX, long n2, long logp) {
	Ciphertext* ctmp = new Ciphertext[4];

	NTL_EXEC_RANGE(4, first, last);
	for (long i = first; i < last; ++i) {
		ctmp[i] = i == 0 ? encv : scheme.leftRotate(encv, 4 - i);
		Ciphertext tmp = scheme.modDownTo(encDX[i], ctmp[i].logq);
		scheme.multAndEqual(ctmp[i], tmp);
		tmp.kill();
	}
	NTL_EXEC_RANGE_END;

	encw = ctmp[0];
	for (long i = 1; i < 4; ++i) {
		scheme.addAndEqual(encw, ctmp[i]);
	}

	delete[] ctmp;

	for (long j = 0; j < 2; ++j) {
		Ciphertext rot = scheme.leftRotateFast(encw, n2 * (1 << j));
		scheme.addAndEqual(encw, rot);
		rot.kill();
	}

	scheme.reScaleByAndEqual(encw, logp);
}

void EncGWAS::evalXtv(Ciphertext& encw, Ciphertext& encv, Ciphertext* encDXExp, long n2, long logp) {
	Ciphertext* ctmp = new Ciphertext[4];

	NTL_EXEC_RANGE(4, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext tmp = scheme.modDownTo(encDXExp[i], encv.logq);
		ctmp[i] = scheme.mult(tmp, encv);
		if(i != 0) {
			scheme.leftRotateAndEqual(ctmp[i], i);
		}
		tmp.kill();
	}
	NTL_EXEC_RANGE_END;

	encw = ctmp[0];
	for (long i = 1; i < 4; ++i) {
		scheme.addAndEqual(encw, ctmp[i]);
	}

	delete[] ctmp;

	long logsl = log2(n2);
	for (long j = 2; j < logsl; ++j){
		Ciphertext rot = scheme.leftRotateFast(encw, (1 << j));
		scheme.addAndEqual(encw, rot);
		rot.kill();
	}

	scheme.reScaleByAndEqual(encw, logp);
}

void EncGWAS::evalBeta(Ciphertext& encBeta, Ciphertext& encAdjU, Ciphertext& encY, Ciphertext* encDXExp, ZZ** maskBeta, long n2, long logp, long logc) {
	Ciphertext encYP = scheme.addConst(encY, -0.5, logp);
	Ciphertext tmp;
	evalXtv(tmp, encYP, encDXExp, n2, logp);
	scheme.modDownToAndEqual(tmp, encAdjU.logq);
	scheme.multAndEqual(tmp, encAdjU);

	for (long i = 1; i < 3; ++i){
		Ciphertext rot = scheme.leftRotateFast(tmp, i);
		scheme.addAndEqual(tmp, rot);
	}

	Ciphertext* ctmp = new Ciphertext[4];

	NTL_EXEC_RANGE(4, first, last);
	for (long i = first; i < last; ++i) {
		ctmp[i] = scheme.multByPoly(tmp, maskBeta[i], logc);
		if(i != 0) {
			scheme.leftRotateFastAndEqual(ctmp[i], (2 * n2 - 1) * i);
		}
	}
	NTL_EXEC_RANGE_END;

	encBeta = ctmp[0];
	for (long i = 1; i < 4; ++i) {
		scheme.addAndEqual(encBeta, ctmp[i]);
	}

//	scheme.divByPo2AndEqual(encBeta, 7); // for arbitrary det!

	delete[] ctmp;

	long logsl = log2(n2 * 8);
	//Normal version	
	for (long i = 2; i < logsl; ++i) {
		Ciphertext rot = scheme.rightRotateFast(encBeta, 1 << i);
		scheme.addAndEqual(encBeta, rot);
	}

	// Multithreading version
	// Ciphertext* rot = new Ciphertext[7];
	// long logsl_thread = (logsl - 2) / 3;
	// for(long j = 0; j < logsl_thread; j++){
	// 	NTL_EXEC_RANGE(7, first, last);
	// 	for(long i = first; i < last; i++){
	// 		rot[i] = scheme.rightRotateFast(encBeta, (1 << (3 * j + 2)) * (i + 1));
	// 	}
	// 	NTL_EXEC_RANGE_END;
	// 	for(long i = 0; i < 7; i++){
	// 		scheme.addAndEqual(encBeta, rot[i]);
	// 	}
	// }
	// for(long j = 3 * logsl_thread; j < logsl - 2; j++){
	// 	Ciphertext rotat = scheme.rightRotateFast(encBeta, (1 << (j + 2)));
	// 	scheme.addAndEqual(encBeta, rotat);
	// }
	// delete[] rot;


	scheme.reScaleByAndEqual(encBeta, logp + logc);
//	scheme.leftShiftAndEqual(encBeta, 3);				//	alpha = 8
//	scheme.divByPo2AndEqual(encBeta, 7);				//	alpha = 1/128
	scheme.multByConstAndEqual(encBeta, 0.01, logc);	// 	alpha = 0.01
	scheme.reScaleByAndEqual(encBeta, logc);
}

void EncGWAS::evaldBeta(Ciphertext& encdBeta, Ciphertext& encAdjU, Ciphertext& encY, Ciphertext& encP, Ciphertext* encDXExp, ZZ** maskBeta,  long n2, long logp, long logc) {
	Ciphertext encYP = scheme.modDownTo(encY, encP.logq);
	scheme.subAndEqual(encYP, encP);

	Ciphertext tmp;
	evalXtv(tmp, encYP, encDXExp, n2, logp);
	scheme.modDownToAndEqual(tmp, encAdjU.logq);
	scheme.multAndEqual(tmp, encAdjU);


	// scheme.multByConstAndEqual(tmp, 0.01, logp); // for arbitrary det!

//	scheme.divByPo2AndEqual(tmp, 7); // for arbitrary det!

	for (long i = 1; i < 3; ++i){
		Ciphertext rot = scheme.leftRotateFast(tmp, i);
		scheme.addAndEqual(tmp, rot);
	}

	Ciphertext* ctmp = new Ciphertext[4];

	NTL_EXEC_RANGE(4, first, last);
	for (long i = first; i < last; ++i) {
		ctmp[i] = scheme.multByPoly(tmp, maskBeta[i], logc);
		if(i != 0) {
			scheme.leftRotateFastAndEqual(ctmp[i], (2 * n2 - 1) * i);
		}
	}
	NTL_EXEC_RANGE_END;

	encdBeta = ctmp[0];
	for (long i = 1; i < 4; i++) {
		scheme.addAndEqual(encdBeta, ctmp[i]);
	}

	delete[] ctmp;

	long logsl = log2(n2 * 8);
	//Normal version
	for (long i = 2; i < logsl; ++i){
		Ciphertext rot = scheme.rightRotateFast(encdBeta, 1 << i);
		scheme.addAndEqual(encdBeta, rot);
	}

	// //Multithreading version
	// Ciphertext* rot = new Ciphertext[7];
	// long logsl_thread = (logsl - 2) / 3;
	// for(long j = 0; j < logsl_thread; j++){
	// 	NTL_EXEC_RANGE(7, first, last);
	// 	for(long i = first; i < last; i++){
	// 		rot[i] = scheme.rightRotateFast(encdBeta, (1 << (3 * j + 2)) * (i + 1));
	// 	}
	// 	NTL_EXEC_RANGE_END;
	// 	for(long i = 0; i < 7; i++){
	// 		scheme.addAndEqual(encdBeta, rot[i]);
	// 	}
	// }
	// for(long j = 3 * logsl_thread; j < logsl - 2; j++){
	// 	Ciphertext rotat = scheme.rightRotateFast(encdBeta, (1 << (j + 2)));
	// 	scheme.addAndEqual(encdBeta, rotat);
	// }
	// delete[] rot;

	scheme.reScaleByAndEqual(encdBeta, logp + logc);
}

void EncGWAS::evalP(Ciphertext& encP, Ciphertext& encXBeta, long logp, long loga, long kdeg) {
	scheme.reScaleByAndEqual(encXBeta, loga);
	Ciphertext encXBeta2 = scheme.square(encXBeta);
	scheme.reScaleByAndEqual(encXBeta2, logp);
	if(kdeg == 3) {
		scheme.addConstAndEqual(encXBeta2, degree3[1] / degree3[2], logp - 2 * loga);
		encP = scheme.multByConst(encXBeta, degree3[2], logp + 3 * loga);
		scheme.reScaleByAndEqual(encP, logp);
		scheme.modDownToAndEqual(encP, encXBeta2.logq);
		scheme.multAndEqual(encP, encXBeta2);
		scheme.reScaleByAndEqual(encP, logp);
		scheme.addConstAndEqual(encP, degree3[0], logp);
	} else if (kdeg == 5) {
		Ciphertext encXBeta4 = scheme.square(encXBeta2);
		scheme.reScaleByAndEqual(encXBeta4, logp);
		scheme.multByConstAndEqual(encXBeta2, degree5[2] / degree5[3], logp - 2 * loga);
		scheme.reScaleByAndEqual(encXBeta2, logp);
		scheme.addAndEqual(encXBeta4, encXBeta2);
		scheme.addConstAndEqual(encXBeta4, degree5[1] / degree5[3], logp - 4 * loga);
		encP = scheme.multByConst(encXBeta, degree5[3], logp + 5 * loga);
		scheme.reScaleByAndEqual(encP, logp);
		scheme.modDownToAndEqual(encP, encXBeta4.logq);
		scheme.multAndEqual(encP, encXBeta4);
		scheme.reScaleByAndEqual(encP, logp);
		scheme.addConstAndEqual(encP, degree5[0], logp);
	} else {
		Ciphertext encXBeta4 = scheme.square(encXBeta2);
		scheme.reScaleByAndEqual(encXBeta4, logp);
		Ciphertext encXBeta2c = scheme.multByConst(encXBeta2, degree7[3] / degree7[4], logp - 2 * loga);
		scheme.reScaleByAndEqual(encXBeta2c, logp);
		scheme.addAndEqual(encXBeta4, encXBeta2c);
		scheme.addConstAndEqual(encXBeta4, degree7[2] / degree7[4], logp - 4 * loga);
		encP = scheme.multByConst(encXBeta, degree7[4], logp + 7 * loga);
		scheme.reScaleByAndEqual(encP, logp);
		scheme.modDownToAndEqual(encP, encXBeta2.logq);
		scheme.multAndEqual(encP, encXBeta2);
		scheme.reScaleByAndEqual(encP, logp);
		scheme.modDownToAndEqual(encP, encXBeta4.logq);
		scheme.multAndEqual(encP, encXBeta4);
		scheme.reScaleByAndEqual(encP, logp);
		Ciphertext tmp = scheme.multByConst(encXBeta, degree7[1], logp + loga);
		scheme.reScaleByAndEqual(tmp, logp);
		scheme.modDownToAndEqual(tmp, encP.logq);
		scheme.addAndEqual(encP, tmp);
		scheme.addConstAndEqual(encP, degree7[0], logp);
	}
}

void EncGWAS::evalW(Ciphertext& encW, Ciphertext& encP, long logp) {
	Ciphertext encP2 = scheme.addConst(encP, -1.0, logp);
	encW = scheme.mult(encP, encP2);
	scheme.reScaleByAndEqual(encW, logp);
	scheme.negateAndEqual(encW);
}

void EncGWAS::evalWS(Ciphertext& encW, Ciphertext* encS, long sNum, long logp) {
	NTL_EXEC_RANGE(sNum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modDownToAndEqual(encS[i], encW.logq);
		scheme.multAndEqual(encS[i], encW);
		scheme.reScaleByAndEqual(encS[i], logp);
	}
	NTL_EXEC_RANGE_END;
}

void EncGWAS::evalV(Ciphertext* encV, Ciphertext* encWS, Ciphertext* encDXExp, ZZ** maskXWS, ZZ** maskPack, long n2, long packNum, long lastpackNum, long packctxNum, long logp, long logc) {
	long tmppackNum;
	for (long i = 0; i < packctxNum; ++i) {
		if (i == packctxNum - 1) {
			tmppackNum = lastpackNum;
		} else {
			tmppackNum = packNum;
		}

		// Ciphertext* ctmp = new Ciphertext[tmppackNum];

		// NTL_EXEC_RANGE(tmppackNum, first, last);
		// for (long j = first; j < last; ++j) {
		// 	evalVi(ctmp[j], encWS[i * packNum + j], encDXExp, maskXWS, n2, logp, logc, 4 * j);
		// 	scheme.multByPolyAndEqual(ctmp[j], maskPack[j], logc);
		// 	scheme.reScaleByAndEqual(ctmp[j], logc);
		// }
		// NTL_EXEC_RANGE_END;

		// encV[i] = ctmp[0];
		// for (long j = 1; j < tmppackNum; ++j) {
		// 	scheme.addAndEqual(encV[i], ctmp[j]);
		// }
		long quo = tmppackNum / 8;
		long rem = tmppackNum - quo * 8;
		for(long k = 0; k < quo; k++) {
			Ciphertext* ctmp = new Ciphertext[8];
			NTL_EXEC_RANGE(8, first, last);
			for(long j = first; j < last; ++j) {
				evalVi(ctmp[j], encWS[i * packNum + j + 8 * k], encDXExp, maskXWS, n2, logp, logc, 4 * (j + 8 * k));
				scheme.multByPolyAndEqual(ctmp[j], maskPack[j + 8 * k], logc);
				scheme.reScaleByAndEqual(ctmp[j], logc);
			}
			NTL_EXEC_RANGE_END;
			if(k == 0){
				encV[i] = ctmp[0];
				for (long j = 1; j < 8; ++j) {
					scheme.addAndEqual(encV[i], ctmp[j]);
				}
			}
			else{
				for (long j = 0; j < 8; ++j) {
					scheme.addAndEqual(encV[i], ctmp[j]);
				}
			}
			delete[] ctmp;
		}
		Ciphertext* ctmp = new Ciphertext[rem];
		NTL_EXEC_RANGE(rem, first, last);
		for(long j = first; j < last; ++j) {
			evalVi(ctmp[j], encWS[i * packNum + j + 8 * quo], encDXExp, maskXWS, n2, logp, logc, 4 * (j + 8 * quo));
			scheme.multByPolyAndEqual(ctmp[j], maskPack[j + 8 * quo], logc);
			scheme.reScaleByAndEqual(ctmp[j], logc);
		}
		NTL_EXEC_RANGE_END;
		for (long j = 0; j < rem; ++j) {
			scheme.addAndEqual(encV[i], ctmp[j]);
		}
		delete[] ctmp;
	}
}

void EncGWAS::evalVi(Ciphertext& encVi, Ciphertext& encWS, Ciphertext* encDXExp, ZZ** maskXWS, long n2, long logp, long logc, long pos) {
	Ciphertext* ctmp = new Ciphertext[4];
	NTL_EXEC_RANGE(4, first, last);
	for (long i = first; i < last; ++i) {
		ctmp[i] = scheme.modDownTo(encDXExp[i], encWS.logq);
		scheme.multAndEqual(ctmp[i], encWS);
		if(i != 0) {
			Ciphertext rot = scheme.leftRotateFast(ctmp[i], i);
			scheme.multByPolyAndEqual(rot, maskXWS[i - 1], logc);
			scheme.rightRotateFastAndEqual(ctmp[i], n2 - i);
			scheme.multByPolyAndEqual(ctmp[i], maskXWS[2 + i], logc);
			scheme.addAndEqual(ctmp[i], rot);
			scheme.reScaleByAndEqual(ctmp[i], logc);
			rot.kill();
		} else {
			scheme.modDownByAndEqual(ctmp[i], logc);
		}
	}
	NTL_EXEC_RANGE_END;

	encVi = ctmp[0];
	for (long i = 1; i < 4; ++i) {
		scheme.addAndEqual(encVi, ctmp[i]);
	}
	scheme.reScaleByAndEqual(encVi, logp);

	delete[] ctmp;

	for (long j = n2 / 2; j >= 4; j >>= 1) {
		if (pos < j) {
			Ciphertext rot = scheme.leftRotateFast(encVi, j);
			scheme.addAndEqual(encVi, rot);
			rot.kill();
		} else {
			Ciphertext rot = scheme.rightRotateFast(encVi, j);
			scheme.addAndEqual(encVi, rot);
			pos -= j;
			rot.kill();
		}
	}
}

void EncGWAS::evalDSWS(Ciphertext* encDSWS, Ciphertext* encWS, ZZ** maskPack, long n2, long packNum, long lastpackNum, long packctxNum, long logp, long logc) {
	long tmppackNum;
	for (long i = 0; i < packctxNum; ++i) {
		if (i == packctxNum - 1) {
			tmppackNum = lastpackNum;
		} else {
			tmppackNum = packNum;
		}

		// Ciphertext* ctmp = new Ciphertext[tmppackNum];

		// NTL_EXEC_RANGE(tmppackNum, first, last);
		// for (long j = first; j < last; ++j) {
		// 	ctmp[j] = encWS[i * packNum + j];
		// 	long pos = 4 * j;
		// 	for (long k = n2 / 2; k >= 1; k >>= 1) {
		// 		if (pos < k) {
		// 			Ciphertext rot = scheme.leftRotateFast(ctmp[j], k);
		// 			scheme.addAndEqual(ctmp[j], rot);
		// 		} else {
		// 			Ciphertext rot = scheme.rightRotateFast(ctmp[j], k);
		// 			scheme.addAndEqual(ctmp[j], rot);
		// 			pos -= k;
		// 		}
		// 	}
		// 	scheme.multByPolyAndEqual(ctmp[j], maskPack[j], logc);
		// }
		// NTL_EXEC_RANGE_END;

		// encDSWS[i] = ctmp[0];
		// for (long j = 1; j < tmppackNum; ++j) {
		// 	scheme.addAndEqual(encDSWS[i], ctmp[j]);
		// }
		// scheme.reScaleByAndEqual(encDSWS[i], logc);

		long quo = tmppackNum / 8;
		long rem = tmppackNum - quo * 8;
		for(int l = 0; l < quo; l++){
			Ciphertext* ctmp = new Ciphertext[8];
			NTL_EXEC_RANGE(8, first, last);
			for (long j = first; j < last; ++j) {
				ctmp[j] = encWS[i * packNum + j + 8 * l];
				long pos = 4 * (j + 8 * l);
				for (long k = n2 / 2; k >= 1; k >>= 1) {
					if (pos < k) {
						Ciphertext rot = scheme.leftRotateFast(ctmp[j], k);
						scheme.addAndEqual(ctmp[j], rot);
						rot.kill();
					} else {
						Ciphertext rot = scheme.rightRotateFast(ctmp[j], k);
						scheme.addAndEqual(ctmp[j], rot);
						pos -= k;
						rot.kill();
					}
				}
				scheme.multByPolyAndEqual(ctmp[j], maskPack[j + 8 * l], logc);
			}
			NTL_EXEC_RANGE_END;
			if(l == 0){
				encDSWS[i] = ctmp[0];
				for (long j = 1; j < 8; ++j) {
					scheme.addAndEqual(encDSWS[i], ctmp[j]);
				}
			}
			else{
				for (long j = 0; j < 8; ++j) {
					scheme.addAndEqual(encDSWS[i], ctmp[j]);
				}
			}
			delete[] ctmp;
		}
		Ciphertext* ctmp = new Ciphertext[rem];
		NTL_EXEC_RANGE(rem, first, last);
		for (long j = first; j < last; ++j) {
				ctmp[j] = encWS[i * packNum + j + 8 * quo];
				long pos = 4 * (j + 8 * quo);
				for (long k = n2 / 2; k >= 1; k >>= 1) {
					if (pos < k) {
						Ciphertext rot = scheme.leftRotateFast(ctmp[j], k);
						scheme.addAndEqual(ctmp[j], rot);
						rot.kill();
					} else {
						Ciphertext rot = scheme.rightRotateFast(ctmp[j], k);
						scheme.addAndEqual(ctmp[j], rot);
						pos -= k;
						rot.kill();
					}
				}
				scheme.multByPolyAndEqual(ctmp[j], maskPack[j + 8 * quo], logc);
			}
			NTL_EXEC_RANGE_END;
			for (long j = 0; j < rem; ++j) {
				scheme.addAndEqual(encDSWS[i], ctmp[j]);
			}
			delete[] ctmp;
		scheme.reScaleByAndEqual(encDSWS[i], logc);
	}
}

void EncGWAS::evalDAdjUExp(Ciphertext* encDUExp, Ciphertext& encAdjU, ZZ** maskDU, long logc) {
	Ciphertext* ctmp = new Ciphertext[4];

	NTL_EXEC_RANGE(4, first, last);
	for (long i = first; i < last; ++i) {
		ctmp[i] = scheme.multByPoly(encAdjU, maskDU[i], logc);
		scheme.reScaleByAndEqual(ctmp[i], logc);
	}
	NTL_EXEC_RANGE_END;

	expandDX(encDUExp, ctmp, 4);
	delete[] ctmp;
}

void EncGWAS::evalDVAdjUV(Ciphertext* encDVUV, Ciphertext* encV, Ciphertext* encDUExp, ZZ** maskDVUV, long packctxNum, long logp, long logc) {
	for (long i = 0; i < packctxNum; ++i) {
		Ciphertext* ctmp = new Ciphertext[4];
		
		NTL_EXEC_RANGE(4, first, last);
		for (long j = first; j < last; ++j) {
			encDUExp[j].n = encV[i].n;				// slots
			ctmp[j] = scheme.modDownTo(encDUExp[j], encV[i].logq);				//logp, logq
			scheme.multAndEqual(ctmp[j], encV[i]);								//2logp, logq
			if(j == 0) {
				scheme.multAndEqual(ctmp[j], encV[i]);							//3logp, logq
				scheme.reScaleByAndEqual(ctmp[j], logp);						//2logp, logq - logp

				////output V^T * (adj U)_j * V_j
				//double* vuv_part = new double[2048 * 32]();
				//complex<double>* tmp1 = scheme.decrypt(secretKey, ctmp[j]);
				//int r = 0;
				//for (int w = 0; w < 4; w++) {
				//	for (int s = 0; s < 64; s++) {
				//		for (int x = 0; x < 256; x++) {
				//			vuv_part[r] = tmp1[w + 256 * x + 4 * s].real();
				//			r++;
				//		}
				//	}
				//}
				//ofstream VUVpart("outVUV_part.csv");
				//for (int i = 0; i < 2048 * 32; i++) {
				//	VUVpart << *vuv_part;
				//	if (i < 2048 * 32 - 1)
				//		VUVpart << "\n";
				//	vuv_part++;
				//}

			}
			else {
				scheme.reScaleByAndEqual(ctmp[j], logp);						//logp, logq - logp
				Ciphertext rot1 = scheme.leftRotateFast(encV[i], 4 - j);		// logp, logq
				Ciphertext rot2 = scheme.rightRotateFast(encV[i], j);			// logp, logq
				scheme.multByPolyAndEqual(rot1, maskDVUV[j - 1], logc);			// logp + logc, logq
				scheme.multByPolyAndEqual(rot2, maskDVUV[j + 2], logc);			// logp + logc, logq
				scheme.addAndEqual(rot2, rot1);									// logp + logc, logq
				scheme.reScaleByAndEqual(rot2, logc);							// logp, logq - logc
				scheme.modDownToAndEqual(rot2, ctmp[j].logq);					// logp, logq - logp			
				scheme.multAndEqual(ctmp[j], rot2);								// 2logp, logq - logp

			////output rot2(rotation of encV[i])
			//	double* rot2plain = new double[2048 * 32]();
			//	rot2.n = 2048 * 32;
			//	complex<double>* tmp = scheme.decrypt(secretKey, rot2);
			//	int r = 0;
			//	for (int w = 0; w < 4; w++) {
			//		for (int s = 0; s < 64; s++) {
			//			for (int x = 0; x < 256; x++) {
			//				rot2plain[r] = tmp[w + 256 * x + 4 * s].real();
			//				r++;
			//			}
			//		}
			//	}
			//	ofstream Vrot("outVtrans_rot.csv" + to_string(j));
			//	for (int i = 0; i < 2048 * 32; i++) {
			//		Vrot << *rot2plain;
			//		if (i < 2048 * 32 - 1)
			//			Vrot << "\n";
			//		rot2plain++;
			//	}

			//	//output V^T * (adj U)_j * V_j
			//	double* vuv_part = new double[2048 * 32]();
			//	complex<double>* tmp1 = scheme.decrypt(secretKey, ctmp[j]);
			//	r = 0;
			//	for (int w = 0; w < 4; w++) {
			//		for (int s = 0; s < 64; s++) {
			//			for (int x = 0; x < 256; x++) {
			//				vuv_part[r] = tmp1[w + 256 * x + 4 * s].real();
			//				r++;
			//			}
			//		}
			//	}
			//	ofstream VUVpart("outVUV_part.csv" + to_string(j));
			//	for (int i = 0; i < 2048 * 32; i++) {
			//		VUVpart << *vuv_part;
			//		if (i < 2048 * 32 - 1)
			//			VUVpart << "\n";
			//		vuv_part++;
			//	}

			}


		}
		NTL_EXEC_RANGE_END;

		encDVUV[i] = ctmp[0];
		for (long j = 1; j < 4; ++j) {
			scheme.addAndEqual(encDVUV[i], ctmp[j]);							// 2logp, logq - logp
		}

		// scheme.modDownByAndEqual(encV[i], logc);								// logp, logq - logc
		// scheme.multAndEqual(encDVUV[i], encV[i]);								// 3logp, logq - logc

		delete[] ctmp;

		for (long j = 0; j < 2; ++j) {
			Ciphertext rot = scheme.leftRotateFast(encDVUV[i], (1 << j));
			// rot.n = encDVUV[i].n;
			scheme.addAndEqual(encDVUV[i], rot);
			rot.kill();
		}
		scheme.reScaleByAndEqual(encDVUV[i], logp);							// logp, logq - 2logp

	}

	//////log p = log c 인 경우!!, rescale 한번으로 줄임

	// for (long i = 0; i < packctxNum; ++i) {
	// 	Ciphertext* ctmp = new Ciphertext[4];
		
	// 	NTL_EXEC_RANGE(4, first, last);
	// 	for (long j = first; j < last; ++j) {
	// 		ctmp[j] = scheme.modDownTo(encDUExp[j], encV[i].logq);				//logp, logq
	// 		scheme.multAndEqual(ctmp[j], encV[i]);								//2logp, logq
	// 		if(j == 0) {
	// 			scheme.multAndEqual(ctmp[j], encV[i]);							//3logp, logq
	// 			scheme.multByConstAndEqual(ctmp[j], 1.0, logp);					//4logp, logq
	// 		} else {
	// 			Ciphertext rot1 = scheme.leftRotateFast(encV[i], 4 - j);		// logp, logq
	// 			Ciphertext rot2 = scheme.rightRotateFast(encV[i], j);			// logp, logq
	// 			scheme.multByPolyAndEqual(rot1, maskDVUV[j - 1], logc);			// 2logp, logq
	// 			scheme.multByPolyAndEqual(rot2, maskDVUV[j + 2], logc);			// 2logp, logq
	// 			scheme.addAndEqual(rot2, rot1);									// 2logp, logq			
	// 			scheme.multAndEqual(ctmp[j], rot2);								// 4logp, logq
	// 		}
	// 	}
	// 	NTL_EXEC_RANGE_END;

	// 	encDVUV[i] = ctmp[0];
	// 	for (long j = 1; j < 4; ++j) {
	// 		scheme.addAndEqual(encDVUV[i], ctmp[j]);							// 4logp, logq
	// 	}

	// 	delete[] ctmp;

	// 	for (long j = 0; j < 2; ++j) {
	// 		Ciphertext rot = scheme.leftRotateFast(encDVUV[i], (1 << j));
	// 		rot.n = encDVUV[i].n;
	// 		scheme.addAndEqual(encDVUV[i], rot);
	// 		rot.kill();
	// 	}
	// 	scheme.reScaleByAndEqual(encDVUV[i], 3 * logp);							// logp, logq - 3logp

	// }
}

void EncGWAS::evalIPS(Ciphertext* encIPS, Ciphertext* encDSWS, Ciphertext* encDVAdjUV, Ciphertext& encDetU, ZZ* maskIPS, long packctxNum, long logp, long logc) {
//	Ciphertext* rot = new Ciphertext[packctxNum];
	NTL_EXEC_RANGE(packctxNum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modDownToAndEqual(encDSWS[i], encDetU.logq);
		scheme.multAndEqual(encDSWS[i], encDetU);
		scheme.reScaleByAndEqual(encDSWS[i], logp);
		scheme.modDownToAndEqual(encDSWS[i], encDVAdjUV[i].logq);
		encIPS[i] = scheme.sub(encDSWS[i], encDVAdjUV[i]);	

//		complex<double>* ddips = scheme.decrypt(secretKey, encIPS[i]);
//		cout << "encIPS " << endl;
//		for (int j = 0; j < 256; j++){
//			cout << ddips[j].real() <<", ";
//			if (j%4 == 3){
//				cout <<"["<<j/4 + 1<<"] ";
//			}
//		}
//		cout << endl;

		scheme.multByPolyAndEqual(encIPS[i], maskIPS, logc);
		scheme.reScaleByAndEqual(encIPS[i], logc);

		// ddips = scheme.decrypt(secretKey, encIPS[i]);
		// cout << "encIPS slots: " << encIPS[i].n << endl;
  //               cout << "ips_after_repac k" << i << " = ";
  //               for(int j = 0; j < 256; j++){
		// 	cout << ddips[j].real() << ", ";
		// 	if (j%4 == 3){
		// 		cout <<"["<<j/4 + 1<<"] ";
		// 	}
  //               }


//		complex<double>* dIPS_dw = scheme.decrypt(secretKey, encIPS[i]);
//              double max = 0;
//		cout << "IPS_enc" << i << " = ";
//                for(int j = 0; j < 128; i++){
//			if (max < dIPS_dw[j].real()){
//				max = dIPS_dw[j].real();
//			}
//		}
//		cout << "max " << i << " is " << max << endl;

		// Ciphertext rot = scheme.rightRotateFast(encIPS[i], 1);
		// scheme.addAndEqual(encIPS[i], rot);
		// rot = scheme.rightRotateFast(encIPS[i], 2);
		// scheme.addAndEqual(encIPS[i], rot);
	}
	NTL_EXEC_RANGE_END;
}

void EncGWAS::evalSWz(Ciphertext* Sz, Ciphertext* encS, Ciphertext& encY, Ciphertext& encP, ZZ** maskPack, long n2, long packNum, long lastpackNum, long packctxNum, long logp, long logc) {
	Ciphertext encYP = scheme.modDownTo(encY, encP.logq);
	scheme.subAndEqual(encYP, encP);
	long tmppackNum;
	for (long i = 0; i < packctxNum; ++i) {
		if (i == packctxNum - 1) {
			tmppackNum = lastpackNum;
		} else {
			tmppackNum = packNum;
		}

		// Ciphertext* ctmp = new Ciphertext[tmppackNum];

		// NTL_EXEC_RANGE(tmppackNum, first, last);
		// for (long j = first; j < last; ++j) {
		// 	ctmp[j] = scheme.modDownTo(encS[i * packNum + j], encYP.logq);
		// 	scheme.multAndEqual(ctmp[j], encYP);
		// 	long pos = 4 * j;
		// 	for (long k = n2 / 2; k >= 1; k >>= 1) {
		// 		if (pos < k) {
		// 			Ciphertext rot = scheme.leftRotateFast(ctmp[j], k);
		// 			scheme.addAndEqual(ctmp[j], rot);
		// 		} else {
		// 			Ciphertext rot = scheme.rightRotateFast(ctmp[j], k);
		// 			scheme.addAndEqual(ctmp[j], rot);
		// 			pos -= k;
		// 		}
		// 	}
		// 	scheme.multByPolyAndEqual(ctmp[j], maskPack[j], logc);
		// }
		// NTL_EXEC_RANGE_END;

		// Sz[i] = ctmp[0];
		// for (long j = 1; j < tmppackNum; ++j) {
		// 	scheme.addAndEqual(Sz[i], ctmp[j]);
		// }
		// scheme.reScaleByAndEqual(Sz[i], logp + logc);

		long quo = tmppackNum / 8;
		long rem = tmppackNum - quo * 8;
		for(int l = 0; l < quo; l++){
			Ciphertext* ctmp = new Ciphertext[8];
			NTL_EXEC_RANGE(8, first, last);
			for (long j = first; j < last; ++j) {
				ctmp[j] = scheme.modDownTo(encS[i * packNum + j + 8 * l], encYP.logq);
				scheme.multAndEqual(ctmp[j], encYP);
				long pos = 4 * (j + 8 * l);
				for (long k = n2 / 2; k >= 1; k >>= 1) {
					if (pos < k) {
						Ciphertext rot = scheme.leftRotateFast(ctmp[j], k);
						scheme.addAndEqual(ctmp[j], rot);
						rot.kill();
					} else {
						Ciphertext rot = scheme.rightRotateFast(ctmp[j], k);
						scheme.addAndEqual(ctmp[j], rot);
						pos -= k;
						rot.kill();
					}
				}
				scheme.multByPolyAndEqual(ctmp[j], maskPack[j + 8 * l], logc);
			}
			NTL_EXEC_RANGE_END;
			if(l == 0){
				Sz[i] = ctmp[0];
				for (long j = 1; j < 8; ++j) {
					scheme.addAndEqual(Sz[i], ctmp[j]);
				}
			}
			else{
				for (long j = 0; j < 8; ++j) {
					scheme.addAndEqual(Sz[i], ctmp[j]);
				}
			}
			delete[] ctmp;
		}
		Ciphertext* ctmp = new Ciphertext[rem];
		NTL_EXEC_RANGE(rem, first, last);
		for (long j = first; j < last; ++j) {
			ctmp[j] = scheme.modDownTo(encS[i * packNum + j + 8 * quo], encYP.logq);
			scheme.multAndEqual(ctmp[j], encYP);
			long pos = 4 * (j + 8 * quo);
			for (long k = n2 / 2; k >= 1; k >>= 1) {
				if (pos < k) {
					Ciphertext rot = scheme.leftRotateFast(ctmp[j], k);
					scheme.addAndEqual(ctmp[j], rot);
					rot.kill();
				} else {
					Ciphertext rot = scheme.rightRotateFast(ctmp[j], k);
					scheme.addAndEqual(ctmp[j], rot);
					pos -= k;
					rot.kill();
				}
			}
			scheme.multByPolyAndEqual(ctmp[j], maskPack[j + 8 * quo], logc);
		}
		NTL_EXEC_RANGE_END;
		for (long j = 0; j < rem; ++j) {
			scheme.addAndEqual(Sz[i], ctmp[j]);
		}
		delete[] ctmp;
		scheme.reScaleByAndEqual(Sz[i], logp + logc);
	}
	encYP.kill();
}

void EncGWAS::evalInvIPS(Ciphertext* encInvIPS, Ciphertext* encIPS, long iter, long logscale, long packctxNum, long logp, long logc){
	NTL_EXEC_RANGE(packctxNum, first, last);
	for (long i = first; i < last; ++i) {

		// Initialization of x, y for iteration
//		// deg 2 approx at x = 1 : 1/x =  1 - (x-1) + (x-1)^2
//		Ciphertext x = scheme.addConst(encIPS[i], -1, logp);
//		Ciphertext x2 = scheme.square(x);
//		scheme.reScaleByAndEqual(x2, logp);
//		scheme.modDownToAndEqual(x, x2.logq);
//		scheme.subAndEqual(x2, x);
//		scheme.addConstAndEqual(x2, 1, logp);
//		x = x2;
//		Ciphertext _x = scheme.negate(x);

		// deg 2 approx at x = 2 : 1/x =  1/2 - (x-2)/4 + (x-2)^2 /8
		Ciphertext x = encIPS[i];
		scheme.divByPo2AndEqual(x, logscale);
		scheme.addConstAndEqual(x, -2, logp);
		Ciphertext x2 = scheme.square(x);
		scheme.reScaleByAndEqual(x2, logp);
		scheme.modDownToAndEqual(x, x2.logq);
		scheme.divByPo2AndEqual(x, 2);
		scheme.divByPo2AndEqual(x2, 3);
		scheme.subAndEqual(x2, x);
		scheme.addConstAndEqual(x2, 0.5, logp);
		x = x2;
		Ciphertext _x = scheme.negate(x);

		//simple approx as 1/2


//		deg 1 approx
//		Ciphertext x = scheme.multByConst(encIPS[i], -0.25, logp);
//		scheme.addConstAndEqual(x, 1, logp + logp);
//		scheme.reScaleByAndEqual(x, logp);
//		Ciphertext _x = scheme.negate(x);

		// y = 1 - x * (approx of 1/x)
		Ciphertext y = encIPS[i];
		scheme.divByPo2AndEqual(y, logscale); 
		scheme.modDownToAndEqual(y, _x.logq);
		scheme.multAndEqual(y, _x);
		scheme.reScaleByAndEqual(y, logp);
		scheme.addConstAndEqual(y, 1, logp);
		scheme.modDownToAndEqual(x, y.logq);

		// Iteration;  x_{i+1} = x_i (y_i + 1),  y_{i+1} = (y_i)^2
		Ciphertext tmp;
		for (long j = 0; j < iter; ++j){
//			clearCpx(y, logp);
//			clearCpx(x, logp);
			scheme.modDownToAndEqual(x, y.logq);
			tmp = scheme.mult(x, y);
			scheme.reScaleByAndEqual(tmp, logp);
			scheme.modDownToAndEqual(x, tmp.logq);
			scheme.addAndEqual(x, tmp);
			scheme.squareAndEqual(y);
			scheme.reScaleByAndEqual(y, logp);

//			cout << "x.logq = " << x.logq << endl;
//			complex<double>* dx = scheme.decrypt(secretKey, x);
//			cout << "x_" << j << " = ";
//			for(int i = 0; i < 128; i++){
//				cout << dx[i].real() << ", ";
//				if (i%4 == 3){
//					cout << "["<< i/4 + 1 << "] ";	
//				}
//			}
//			cout << endl;
//			cout << "y.logq = " << y.logq << endl;
//			complex<double>* dy = scheme.decrypt(secretKey, y);
//			cout << "y_"<< j <<" = ";
//			for(int i = 0; i < 128; i++){
//				cout << dy[i].real() << ", ";
//				if (i%4 == 3){
//					cout << "[" << i/4 + 1 << "] ";
//				}
//			}
//			cout << endl;

		}

		scheme.divByPo2AndEqual(x, logscale);
		encInvIPS[i] = x;
		x.kill();
		x2.kill();
		_x.kill();
		y.kill();
		tmp.kill();
	}
	NTL_EXEC_RANGE_END;
}

void EncGWAS::showMat(complex<double>* dvec, long n, long k) {
	cout << "-------" << endl;
	for (long i = 0; i < n; ++i) {
		cout << "[";
		for (int j = 0; j < k; ++j) {
			cout << setw(10) << dvec[i + j * n].real() << ",";
		}
		cout << "]" << endl;
	}
	cout << "-------" << endl;
}


void EncGWAS::showByVec(complex<double>* dvec, long n) {
	cout << "-------" << endl;
	cout << "{";
	for (int i = 0; i < n; ++i) {
		cout << dvec[i].real() << ",";
	}
	cout << "}";
}
