/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef ENCGWAS_H_
#define ENCGWAS_H_

#include <complex>
#include <Scheme.h>
#include <SecretKey.h>
#include <Ciphertext.h>

using namespace std;
using namespace NTL;

static double degree3[3] = {0.5,0.15012,-0.001593};
static double degree5[4] = {0.5,0.19131,-0.0045963, 0.0000412332};
static double degree7[5] = {0.5,0.216884,-0.00819276,0.000165861,-0.00000119581};

class EncGWAS {
public:

	Scheme& scheme;
	SecretKey& secretKey;

	EncGWAS(Scheme& scheme, SecretKey& secretKey);

	void addGWASKeys(long n2);

	void encryptX(Ciphertext& encX, double** X, long n2, long logp, long logq);
	void encryptY(Ciphertext& encY, double* Y, long n2, long logp, long logq);
	void encryptS(Ciphertext* encS, double** S, long n2, long m, long ScolNum, long sNum, long logp, long logq);
	void encryptDX(Ciphertext* encDX, double** X, long n2, long logp, long logq);
	void expandDX(Ciphertext* encDXExp, Ciphertext* encDX, long n2);
	void expandDXnew(Ciphertext* encDXExp, Ciphertext* encDX, long n2);

	ZZ* maskAdjU(long n2, long logc);
	ZZ* maskDetU(long n2, long logc);
	ZZ** maskBeta(long n2, long logc);
	ZZ** maskRU(long n2, long logc);
	ZZ** maskXWS(long n2, long logc);
	ZZ** maskDU(long logc);
	ZZ** maskDVUV(long logc);
	ZZ** maskPack(long n2, long packNum, long logc);
	ZZ* maskIPS(long n2, long logc);

	void evalUinit(Ciphertext& encU, Ciphertext& encX, Ciphertext* encDXExp, long n2, long logp);
	void evalU(Ciphertext& encU, Ciphertext& encX, Ciphertext& encW, Ciphertext* encDXExp, long n2, long logp);
	void evalAdjU(Ciphertext& encAdjU, Ciphertext& encU, ZZ* maskAdjU, long n2, long logp, long logc);
	void evalDetU(Ciphertext& encDetU, Ciphertext& encU, Ciphertext& encAdjU, ZZ* maskDetU, long n2, long logp, long logc);
	void rearrange(Ciphertext& encU, ZZ** maskRU, long n2, long logp, long logc);

	void evalXv(Ciphertext& encw, Ciphertext& encv, Ciphertext* encDX, long n2, long logp);
	void evalXtv(Ciphertext& encw, Ciphertext& encv, Ciphertext* encDXExp, long n2, long logp);

	void evalBeta(Ciphertext& encBeta, Ciphertext& encAdjU, Ciphertext& encY, Ciphertext* encDXExp, ZZ** maskBeta, long n2, long logp, long logc);
	void evaldBeta(Ciphertext& encdBeta, Ciphertext& encAdjU, Ciphertext& encY, Ciphertext& encP, Ciphertext* encDXExp, ZZ** maskBeta, long n2, long logp, long logc);
	void evalP(Ciphertext& encP, Ciphertext& encXBeta, long logp, long loga, long kdeg);
	void evalW(Ciphertext& encW, Ciphertext& encP, long logp);

	void evalWS(Ciphertext& encW, Ciphertext* encS, long sNum, long logp);
	void evalV(Ciphertext* encV, Ciphertext* encWS, Ciphertext* encDXExp, ZZ** maskXWS, ZZ** maskPacking, long n2, long packNum, long lastpackNum, long packctxNum, long logp, long logc);
	void evalVi(Ciphertext& encU, Ciphertext& encWS, Ciphertext* encDXExp, ZZ** maskXWS, long n2, long logp, long logc, long pos);
	void evalDSWS(Ciphertext* encdiagSWS, Ciphertext* encWS, ZZ** maskPack, long n2, long packNum, long lastpackNum, long packctxNum, long logp, long logc);
	void evalDAdjUExp(Ciphertext* encDUExp, Ciphertext& encAdjU, ZZ** maskDU, long logc);
	void evalDVAdjUV(Ciphertext* encdiagVadjUV, Ciphertext* encV, Ciphertext* encadjUdiagExp, ZZ** maskCyclicRot, long packctxNum, long logp, long logc);
	void evalIPS(Ciphertext* encIPS, Ciphertext* encDSWS, Ciphertext* encDVAdjUV, Ciphertext& encDetU, ZZ* maskIPS, long packctxNum, long logp, long logc);
	void evalSWz(Ciphertext* Sz_star, Ciphertext* encS, Ciphertext& encY, Ciphertext& encP, ZZ** maskPack, long n2, long packNum, long lastpackNum, long packctxNum, long logp, long logc);
	void evalInvIPS(Ciphertext* encInvIPS, Ciphertext* encIPS, long iter, long logscale, long packctxNum, long logp, long logc);


	/* Check Functions */
	void showMat(complex<double>* dvec, long n, long k);
	void showByVec(complex<double>* dvec, long n);

};

#endif /* ENCGWAS_H_ */
