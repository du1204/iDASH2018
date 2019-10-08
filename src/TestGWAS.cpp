/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "TestGWAS.h"

#include <iostream>
#include "Extractor.h"
#include "EncGWAS.h"
#include "MemoryUsage.h"
#include <StringUtils.h>
#include <TimeUtils.h>
#include <SerializationUtils.h>
#include <NTL/BasicThreadPool.h>
#include <math.h>

using namespace std;

void TestGWAS::testEnc(double*& Statsq, double*& IPS, double*& InvIPS, double*& Sz, double** X, double* y, double** S, long n2, long m, long logN,
	long logQ, long logp, long logc, long logbp, long logbt, long logbT, long logbI, long fnum, long kdeg, long logscale) {
	cout << "========================" << endl;
	cout << "Start Test Enc" << endl;
	cout << "========================" << endl;

	long loga = 2;

	long ScolNum = (1 << (logN - 1)) / n2;
	long sNum = ceil(float(m) / float(ScolNum));

	long packNum = min(n2 >> 2, sNum);
	long packctxNum = ceil(float(sNum) / float(packNum));
	long lastpackNum = sNum - packNum * (packctxNum - 1);

	srand(time(NULL));
	SetNumThreads(8);
	//SetNumThreads(4);
	TimeUtils timeutils;

	size_t currentAfterSchemeSize = getCurrentRSS() >> 20;
	size_t peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB" << endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB" << endl;
	//-----------------------------------------
	timeutils.start("Scheme Generation");
	Ring ring(logN, logQ, 3.2, 78);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring, true);
	EncGWAS encGWAS(scheme, secretKey);
	encGWAS.addGWASKeys(n2);
	scheme.addBootKey(secretKey, 2, logbt + logbI); // for boot beta
	scheme.addBootKey(secretKey, 4, logbt + logbI); // for boot AdjU
	scheme.addBootKey(secretKey, 0, logbt + logbI); // for boot DetU
	timeutils.stop("Scheme Generation");
	//-----------------------------------------
	currentAfterSchemeSize = getCurrentRSS() >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB" << endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB" << endl;
	//-----------------------------------------
	timeutils.start("Mask Generation");
	ZZ* maskAdjU = encGWAS.maskAdjU(n2, logc);
	ZZ* maskDetU = encGWAS.maskDetU(n2, logc);
	ZZ** maskBeta = encGWAS.maskBeta(n2, logc);
	ZZ** maskRU = encGWAS.maskRU(n2, logc);
	ZZ** maskXWS = encGWAS.maskXWS(n2, logc);
	ZZ** maskDU = encGWAS.maskDU(logc);
	ZZ** maskDVUV = encGWAS.maskDVUV(logc);
	ZZ** maskPack = encGWAS.maskPack(n2, packNum, logc);
	ZZ* maskIPS = encGWAS.maskIPS(n2, logc);
	timeutils.stop("Mask Generation");
	//-----------------------------------------
	currentAfterSchemeSize = getCurrentRSS() >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB" << endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB" << endl;
	//-----------------------------------------
	Ciphertext encX, encY, encW, encP, encBeta, encdBeta, encXBeta, encU, encAdjU, encDetU;
	Ciphertext *encSvec, *encDXvec, *encDXExpvec, *encSzvec, *encDSWSvec, *encVvec,
		*encDUExpvec, *encDVUVvec, *encIPSvec, *encInvIPSvec, *encSzsqrvec, *encStatvec;

	timeutils.start("encrypt");
	encGWAS.encryptX(encX, X, n2, logp, logQ);
	encGWAS.encryptY(encY, y, n2, logp, logQ);
	encDXvec = new Ciphertext[4];
	encGWAS.encryptDX(encDXvec, X, n2, logp, logQ);
	encDXExpvec = new Ciphertext[4];
	encGWAS.expandDX(encDXExpvec, encDXvec, n2);
	encSvec = new Ciphertext[sNum];
	encGWAS.encryptS(encSvec, S, n2, m, ScolNum, sNum, logp, logQ);
	timeutils.stop("encrypt");

	string folder = "../sercipher/";

	for (long i = 0; i < sNum; ++i) {
		SerializationUtils::writeCiphertext(&encSvec[i], folder + "encSvec" + to_string(i) + ".txt");
		encSvec[i].kill();
	}
	delete[] encSvec;

	//-----------------------------------------
	currentAfterSchemeSize = getCurrentRSS() >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB" << endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB" << endl;
	//-----------------------------------------
	// FISHER SCORING
	//-----------------------------------------
	timeutils.start("step [1] - Fisher Scoring");
	encGWAS.evalUinit(encU, encX, encDXExpvec, n2, logp);
	encGWAS.evalAdjU(encAdjU, encU, maskAdjU, n2, logp, logc);
	encGWAS.evalBeta(encBeta, encAdjU, encY, encDXExpvec, maskBeta, n2, logp, logc);
	complex<double>* dbeta = scheme.decrypt(secretKey, encBeta);
	encGWAS.showMat(dbeta, 1, 4);
	delete[] dbeta;
	for (long i = 0; i < fnum; ++i) {
		encGWAS.evalXv(encXBeta, encBeta, encDXvec, n2, logp);
		encGWAS.evalP(encP, encXBeta, logp, loga, kdeg);
		encGWAS.evalW(encW, encP, logp);
		encGWAS.evalU(encU, encX, encW, encDXExpvec, n2, logp);
		encGWAS.evalAdjU(encAdjU, encU, maskAdjU, n2, logp, logc);
		encGWAS.evaldBeta(encdBeta, encAdjU, encY, encP, encDXExpvec, maskBeta, n2, logp, logc);
		encdBeta.n = 4;
		scheme.reScaleByAndEqual(encdBeta, logp - logbp);
		scheme.bootstrapAndEqual(encdBeta, logbt, logQ, logbT, logbI);
		scheme.leftShiftAndEqual(encdBeta, 3 + logp - logbp);
		encdBeta.n = n2;
		encdBeta.logp = logp;
		scheme.modDownToAndEqual(encBeta, encdBeta.logq);
		scheme.addAndEqual(encBeta, encdBeta);
		complex<double>* dbeta = scheme.decrypt(secretKey, encBeta);
		encGWAS.showMat(dbeta, 1, 4);
		delete[] dbeta;
	}
	timeutils.stop("step [1]");


	SerializationUtils::writeCiphertext(&encX, folder + "encX.txt");
	SerializationUtils::writeCiphertext(&encY, folder + "encY.txt");
	SerializationUtils::writeCiphertext(&encdBeta, folder + "encdBeta.txt");
	SerializationUtils::writeCiphertext(&encU, folder + "encU.txt");
	SerializationUtils::writeCiphertext(&encAdjU, folder + "encAdjU.txt");


//	delete &encX; delete &encY; delete &encdBeta; delete &encU; delete &encAdjU;
	encX.kill();
	encY.kill();
	encdBeta.kill();
	encU.kill();
	encAdjU.kill();

	for (long i = 0; i < 4; ++i) {
		SerializationUtils::writeCiphertext(&encDXExpvec[i], folder + "encDXExpvec" + to_string(i) + ".txt");
		encDXExpvec[i].kill();
	}
	delete[] encDXExpvec;

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [2] - Eval W");
	encGWAS.evalXv(encXBeta, encBeta, encDXvec, n2, logp);
	encGWAS.evalP(encP, encXBeta, logp, loga, kdeg);
	encGWAS.evalW(encW, encP, logp);
	timeutils.stop("step [2]");

	SerializationUtils::writeCiphertext(&encXBeta, folder + "encXBeta.txt");
	SerializationUtils::writeCiphertext(&encBeta, folder + "encBeta.txt");
	SerializationUtils::writeCiphertext(&encP, folder + "encP.txt");
//	delete &encXBeta; delete &encBeta; delete &encP; 
	encXBeta.kill();
	encBeta.kill();
	encP.kill();
	for(long i = 0; i < 4; i++){
		encDXvec[i].kill();
	}
	delete[] encDXvec;

	encU = *SerializationUtils::readCiphertext(folder + "encU.txt");
	encX = *SerializationUtils::readCiphertext(folder + "encX.txt");
	encAdjU = *SerializationUtils::readCiphertext(folder + "encAdjU.txt");

	encDXExpvec = new Ciphertext[4];
	for (long i = 0; i < 4; ++i) {
		encDXExpvec[i] = *SerializationUtils::readCiphertext(folder + "encDXExpvec" + to_string(i) + ".txt");
	}

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [3] & [4] - Eval U, detU, AdjU");

	encGWAS.evalU(encU, encX, encW, encDXExpvec, n2, logp);
	encGWAS.evalAdjU(encAdjU, encU, maskAdjU, n2, logp, logc);
	encGWAS.evalDetU(encDetU, encU, encAdjU, maskDetU, n2, logp, logc);
	encGWAS.rearrange(encU, maskRU, n2, logp, logc);
	encGWAS.rearrange(encAdjU, maskRU, n2, logp, logc);

	//encAdjU.n = 16;
	//scheme.reScaleByAndEqual(encAdjU, logp - logbp);
	//scheme.bootstrapAndEqual(encAdjU, logbt, logQ, logbT, logbI);
	//scheme.leftShiftAndEqual(encAdjU, logp - logbp);
	//encAdjU.n = n2;
	//encAdjU.logp = logp;

	//encDetU.n = 1;
	//scheme.reScaleByAndEqual(encDetU, logp - logbp);
	//scheme.bootstrapAndEqual(encDetU, logbt, logQ, logbT, logbI);
	//scheme.leftShiftAndEqual(encDetU, logp - logbp);
	//encDetU.n = n2;
	//encDetU.logp = logp;
	timeutils.stop("step [3] & [4]");

	SerializationUtils::writeCiphertext(&encU, folder + "encU.txt");
	SerializationUtils::writeCiphertext(&encX, folder + "encX.txt");
	SerializationUtils::writeCiphertext(&encW, folder + "encW.txt");
	SerializationUtils::writeCiphertext(&encAdjU, folder + "encAdjU.txt");
	SerializationUtils::writeCiphertext(&encDetU, folder + "encDetU.txt");
//	delete &encU; delete &encX; delete &encW; delete &encAdjU; delete &encDetU;
	encU.kill();
	encX.kill();
	encW.kill();
	encAdjU.kill();
	encDetU.kill();

	for (long i = 0; i < 4; ++i) {
		SerializationUtils::writeCiphertext(&encDXExpvec[i], folder + "encDXExpvec" + to_string(i) + ".txt");
		encDXExpvec[i].kill();
	}
	delete[] encDXExpvec;

	encY = *SerializationUtils::readCiphertext(folder + "encY.txt");
	encP = *SerializationUtils::readCiphertext(folder + "encP.txt");
	encSzvec = new Ciphertext[packctxNum];
	encSvec = new Ciphertext[sNum];
	for (long i = 0; i < sNum; ++i) {
		encSvec[i] = *SerializationUtils::readCiphertext(folder + "encSvec" + to_string(i) + ".txt");
	}

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [5] : Compute S^t (y-p)");
	encGWAS.evalSWz(encSzvec, encSvec, encY, encP, maskPack, n2, packNum, lastpackNum, packctxNum, logp, logc);
	timeutils.stop("step [5]");

	// complex<double>* dectest;
	// cout << "encSz " << endl;
	// for (int i = 0; i < packctxNum; i++) {
	// 	dectest = scheme.decrypt(secretKey, encSzvec[i]);
	// 	for (int j = 0; j < (n2 >> 2); j++) {
	// 		for (int k = 0; k < ScolNum; k++) {
	// 			cout << (1 << (logN-3))*i + j*ScolNum + k << ": " << dectest[n2*k + 4*j].real() <<", ";			
	// 		}
	// 		cout << endl;
	// 	}
	// }

	SerializationUtils::writeCiphertext(&encY, folder + "encY.txt");
	SerializationUtils::writeCiphertext(&encP, folder + "encP.txt");
//	delete &encY; delete &encP;
	encY.kill();
	encP.kill();

	for (long i = 0; i < packctxNum; ++i) {
		SerializationUtils::writeCiphertext(&encSzvec[i], folder + "encSzvec" + to_string(i) + ".txt");
		encSzvec[i].kill();
	}
	delete[] encSzvec;

	encW = *SerializationUtils::readCiphertext(folder + "encW.txt");

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [6] : Convert encS into encWS");
	encGWAS.evalWS(encW, encSvec, sNum, logp);
	timeutils.stop("step [6]");

	SerializationUtils::writeCiphertext(&encW, folder + "encW.txt");
	//delete &encW;
	encW.kill();

	encDSWSvec = new Ciphertext[packctxNum];

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [7] : Compute diag(S_i^t W S_i) and Repack");
	encGWAS.evalDSWS(encDSWSvec, encSvec, maskPack, n2, packNum, lastpackNum, packctxNum, logp, logc);
	timeutils.stop("step [7]");


	for (long i = 0; i < packctxNum; ++i) {
		SerializationUtils::writeCiphertext(&encDSWSvec[i], folder + "encDSWSvec" + to_string(i) + ".txt");
	}
	delete[] encDSWSvec;

	encVvec = new Ciphertext[packctxNum];
	encDXExpvec = new Ciphertext[4];
	for (long i = 0; i < 4; ++i) {
		encDXExpvec[i] = *SerializationUtils::readCiphertext(folder + "encDXExpvec" + to_string(i) + ".txt");
	}

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [8] : Compute V = X^t W S_i and Repack into X^t W S");
	encGWAS.evalV(encVvec, encSvec, encDXExpvec, maskXWS, maskPack, n2, packNum, lastpackNum, packctxNum, logp, logc);
	timeutils.stop("step [8]");

	for (long i = 0; i < 4; ++i) {
		SerializationUtils::writeCiphertext(&encDXExpvec[i], folder + "encDXExpvec" + to_string(i) + ".txt");
		encDXExpvec[i].kill();
	}
	delete[] encDXExpvec;
	for(long i = 0; i < sNum; i++){
		encSvec[i].kill();
	}
	delete[] encSvec;

	encAdjU = *SerializationUtils::readCiphertext(folder + "encAdjU.txt");

	encDUExpvec = new Ciphertext[4];
	encDVUVvec = new Ciphertext[packctxNum];

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [9] : Compute diag(V^t adjU V)");
	encGWAS.evalDAdjUExp(encDUExpvec, encAdjU, maskDU, logc);
	encGWAS.evalDVAdjUV(encDVUVvec, encVvec, encDUExpvec, maskDVUV, packctxNum, logp, logc);
	timeutils.stop("step [9]");


//	delete &encAdjU;
	encAdjU.kill();
	for(long i = 0; i < packctxNum; i++){
		encVvec[i].kill();
	}
	delete[] encVvec;
	for(long i = 0; i < 4; i++){
		encDUExpvec[i].kill();
	}
	delete[] encDUExpvec;

	encDetU = *SerializationUtils::readCiphertext(folder + "encDetU.txt");
	encIPSvec = new Ciphertext[packctxNum];
	encDSWSvec = new Ciphertext[packctxNum];

	for (long i = 0; i < packctxNum; ++i) {
		encDSWSvec[i] = *SerializationUtils::readCiphertext(folder + "encDSWSvec" + to_string(i) + ".txt");
	}



	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [10] : (detU * [7]) - [9] = det * diag(S*^t W S*) ");
	encGWAS.evalIPS(encIPSvec, encDSWSvec, encDVUVvec, encDetU, maskIPS, packctxNum, logp, logc);
	timeutils.stop("step [10]");

	cout << "final logQ: " << encIPSvec[0].logq << endl;

	// cout << "encIPS logq: " << encIPSvec[0].logq << endl;

	// cout << "encIPS " << endl;
	// for (int i = 0; i < packctxNum; i++) {
	// 	dectest = scheme.decrypt(secretKey, encIPSvec[i]);
	// 	for (int j = 0; j < (n2 >> 2); j++) {
	// 		for (int k = 0; k < ScolNum; k++) {
	// 			cout << (1 << (logN-3))*i + j*ScolNum + k << ": " << dectest[n2*k + 4*j].real() <<", ";			
	// 		}
	// 		cout << endl;
	// 	}
	// }

	SerializationUtils::writeCiphertext(&encDetU, folder + "encDetU.txt");
//	delete &encDetU;
	encDetU.kill();
	for(long i = 0; i < packctxNum; i++){
		encDSWSvec[i].kill();
		encDVUVvec[i].kill();
	}
	delete[] encDSWSvec;
	delete[] encDVUVvec;

	encInvIPSvec = new Ciphertext[packctxNum];

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;
	cout << "IPS" << encIPSvec[0].logq << endl;

	timeutils.start("step [11] : Calculate the square of Statistics = det * [5]^2 / [10] ");
	encGWAS.evalInvIPS(encInvIPSvec, encIPSvec, 5, logscale, packctxNum, logp, logc);
	cout << "invIPS" << encInvIPSvec[0].logq << endl;

	// cout << "encInvIPS.logq = " << encInvIPSvec[0].logq;
	// for(long i = 0; i < packctxNum; i++){
	// 	encIPSvec[i].kill();
	// }

	// cout << "encInvIPS " << endl;
	// for (int i = 0; i < packctxNum; i++) {
	// 	dectest =scheme.decrypt(secretKey, encInvIPSvec[i]);
	// 	for (int j = 0; j < (n2 >> 2); j++) {
	// 		for (int k = 0; k < ScolNum; k++) {
	// 			cout << (1 << (logN-3))*i + j*ScolNum + k << ": " << dectest[n2*k + 4*j].real() <<", ";			
	// 		}
	// 		cout << endl;
	// 	}
	// }

	delete[] encIPSvec;

	encDetU = *SerializationUtils::readCiphertext(folder + "encDetU.txt");
	encSzsqrvec = new Ciphertext[packctxNum];
	encStatvec = new Ciphertext[packctxNum];
	encSzvec = new Ciphertext[packctxNum];
	for (long i = 0; i < packctxNum; ++i) {
		encSzvec[i] = *SerializationUtils::readCiphertext(folder + "encSzvec" + to_string(i) + ".txt");
	}

	for(int i = 0; i < packctxNum; i++){
		encSzsqrvec[i] = scheme.mult(encSzvec[i], encSzvec[i]);
		scheme.reScaleByAndEqual(encSzsqrvec[i], logp);
		scheme.modDownToAndEqual(encDetU, encSzsqrvec[i].logq);
		scheme.multAndEqual(encSzsqrvec[i], encDetU);
		scheme.reScaleByAndEqual(encSzsqrvec[i], logp);
		scheme.modDownToAndEqual(encSzsqrvec[i], encInvIPSvec[i].logq);
		encStatvec[i] = scheme.mult(encSzsqrvec[i], encInvIPSvec[i]);
		scheme.reScaleByAndEqual(encStatvec[i], logp);
	}
	timeutils.stop("step [11]");
	cout << "final logq: " << encStatvec[0].logq << endl;

	// cout << "encStatvec.logq = " << encStatvec[0].logq;

//	delete &encDetU;
	encDetU.kill();
	for(long i = 0; i < packctxNum; i++){
		encInvIPSvec[i].kill();
		encSzsqrvec[i].kill();
		encSzvec[i].kill();
	}
	delete[] encInvIPSvec; delete[] encSzsqrvec; delete[] encSzvec;

	cout <<"End of the Encrpyted Semi-Parallel GWAS Computations!" << endl;
	//-----------------------------------------
	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;
	//-----------------------------------------
	long m2 = sNum * ScolNum;
	Statsq = new double[m2]();
	timeutils.start("decryption");
	for (long k = 0; k < packctxNum; ++k) {
		long tmppackNum = (k == packctxNum - 1) ? lastpackNum : packNum;
		encStatvec[k].n = n2 * ScolNum;
		complex<double>* dStatsq = scheme.decrypt(secretKey, encStatvec[k]);
		for (long i = 0; i < tmppackNum; ++i) {
			for (long j = 0; j < ScolNum; ++j) {
				Statsq[k * ScolNum * packNum + i * ScolNum + j] = dStatsq[(i << 2) + j * n2].real();
			}
		}
		delete[] dStatsq;
	}
	timeutils.stop("decryption");
	//-----------------------------------------
	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;
	//-----------------------------------------
	cout << "========================" << endl;
	cout << "End Test Enc" << endl;
	cout << "========================" << endl;
}


void TestGWAS::testEnc_seperate(double*& numerator, double*& IPS, double** X, double* y, double** S, long n2, long m, long logN,
			long logQ, long logp, long logc, long logbp, long logbt, long logbT, long logbI, long fnum, long kdeg){
	cout << "========================" << endl;
	cout << "Start Test Enc Separate" << endl;
	cout << "========================" << endl;

	long loga = 2;

	long ScolNum = (1 << (logN - 1)) / n2;
	long sNum = ceil(float(m) / float(ScolNum));

	long packNum = min(n2 >> 2, sNum);
	long packctxNum = ceil(float(sNum) / float(packNum));
	long lastpackNum = sNum - packNum * (packctxNum - 1);

	srand(time(NULL));
	//SetNumThreads(8);
	SetNumThreads(8);
	TimeUtils timeutils;

	size_t currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	size_t peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;
	//-----------------------------------------
	timeutils.start("Scheme Generation");
	Ring ring(logN, logQ, 3.2, 56); // 91
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring, true);
	EncGWAS encGWAS(scheme, secretKey);
	encGWAS.addGWASKeys(n2);
	scheme.addBootKey(secretKey, 2, logbt + logbI); // for boot beta
	scheme.addBootKey(secretKey, 4, logbt + logbI); // for boot AdjU
	scheme.addBootKey(secretKey, 0, logbt + logbI); // for boot DetU
	timeutils.stop("Scheme Generation");
	//-----------------------------------------
	timeutils.start("Mask Generation");
	ZZ* maskAdjU = encGWAS.maskAdjU(n2, logc);
	ZZ* maskDetU = encGWAS.maskDetU(n2, logc);
	ZZ** maskBeta = encGWAS.maskBeta(n2, logc);
	ZZ** maskRU = encGWAS.maskRU(n2, logc);
	ZZ** maskXWS = encGWAS.maskXWS(n2, logc);
	ZZ** maskDU = encGWAS.maskDU(logc);
	ZZ** maskDVUV = encGWAS.maskDVUV(logc);
	ZZ** maskPack = encGWAS.maskPack(n2, packNum, logc);
	ZZ* maskIPS = encGWAS.maskIPS(n2, logc);
	timeutils.stop("Mask Generation");
	//-----------------------------------------
	Ciphertext encX, encY, encW, encP, encBeta, encdBeta, encXBeta, encU, encAdjU, encDetU;
	Ciphertext *encSvec, *encDXvec, *encDXExpvec, *encSzvec, *encDSWSvec, *encVvec, *encDUExpvec, *encDVUVvec, *encIPSvec;

	timeutils.start("encrypt");
	encGWAS.encryptX(encX, X, n2, logp, logQ);
	encGWAS.encryptY(encY, y, n2, logp, logQ);
	encDXvec = new Ciphertext[4];
	encGWAS.encryptDX(encDXvec, X, n2, logp, logQ);
	encDXExpvec = new Ciphertext[4];
	encGWAS.expandDX(encDXExpvec, encDXvec, n2);
	encSvec = new Ciphertext[sNum];
	encGWAS.encryptS(encSvec, S, n2, m, ScolNum, sNum, logp, logQ);
	timeutils.stop("encrypt");

	string folder = "../sercipher/";

	for (long i = 0; i < sNum; ++i) {
		SerializationUtils::writeCiphertext(&encSvec[i], folder + "encSvec" + to_string(i) + ".txt");
		encSvec[i].kill();
	}
	delete[] encSvec;

	//-----------------------------------------
	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;
	//-----------------------------------------
	// FISHER SCORING
	//-----------------------------------------
	timeutils.start("step [1] - Fisher Scoring");
	encGWAS.evalUinit(encU, encX, encDXExpvec, n2, logp);
	encGWAS.evalAdjU(encAdjU, encU, maskAdjU, n2, logp, logc);
	encGWAS.evalBeta(encBeta, encAdjU, encY, encDXExpvec, maskBeta, n2, logp, logc);
	readcsv(scheme, secretKey, encAdjU, "adjU_initial");
//	encAdjU.n = 16;
	
	cout << "Beta logQ initial: " << encBeta.logq << endl;
	for (long i = 0; i < fnum; ++i) {
		encGWAS.evalXv(encXBeta, encBeta, encDXvec, n2, logp);
		encGWAS.evalP(encP, encXBeta, logp, loga, kdeg);
		encGWAS.evalW(encW, encP, logp);
		encGWAS.evalU(encU, encX, encW, encDXExpvec, n2, logp);
		encGWAS.evalAdjU(encAdjU, encU, maskAdjU, n2, logp, logc);
		encGWAS.evaldBeta(encdBeta, encAdjU, encY, encP, encDXExpvec, maskBeta, n2, logp, logc);

		cout << "dBeta logQ: " << encdBeta.logq << endl;

		encdBeta.n = 4;
		scheme.reScaleByAndEqual(encdBeta, logp - logbp);
		scheme.bootstrapAndEqual(encdBeta, logbt, logQ, logbT, logbI);
//		scheme.leftShiftAndEqual(encdBeta, 3 + logp - logbp);				//	alpha = 8
//		scheme.leftShiftAndEqual(encdBeta, logp - logbp - 7);				// 	alpha = 1/128
		scheme.leftShiftAndEqual(encdBeta, logp - logbp);
		scheme.multByConstAndEqual(encdBeta, 0.01, logc);					// 	alpha = 0.01
		scheme.reScaleByAndEqual(encdBeta, logc);

		encdBeta.n = n2;
		encdBeta.logp = logp;
		scheme.modDownToAndEqual(encBeta, encdBeta.logq);
		scheme.addAndEqual(encBeta, encdBeta);

		readcsv(scheme, secretKey, encAdjU, "adjU_Fisheriter" + to_string(i));
//		encAdjU.n = 16;
	
		cout << "dBeta logQ after " << i <<"-th iteration, and boot: " << encdBeta.logq << endl;
	}
	timeutils.stop("step [1]");

	//double true_beta[4] = { 0.6742541, 4.0439744, 5.5658792, -3.7415836 };
	//encBeta = scheme.encrypt(true_beta, 4, encBeta.logp, encBeta.logq);


	cout << "Beta logQ after: "<< encBeta.logq << endl;

	SerializationUtils::writeCiphertext(&encX, folder + "encX.txt");
	SerializationUtils::writeCiphertext(&encY, folder + "encY.txt");
	SerializationUtils::writeCiphertext(&encdBeta, folder + "encdBeta.txt");
	SerializationUtils::writeCiphertext(&encU, folder + "encU.txt");
	SerializationUtils::writeCiphertext(&encAdjU, folder + "encAdjU.txt");

	//	delete &encX; delete &encY; delete &encdBeta; delete &encU; delete &encAdjU;
	encX.kill();
	encY.kill();
	encdBeta.kill();
	encU.kill();
	encAdjU.kill();

	for (long i = 0; i < 4; ++i) {
		SerializationUtils::writeCiphertext(&encDXExpvec[i], folder + "encDXExpvec" + to_string(i) + ".txt");
		encDXExpvec[i].kill();
	}
	delete[] encDXExpvec;

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [2] - Eval W");
	encGWAS.evalXv(encXBeta, encBeta, encDXvec, n2, logp);
	encGWAS.evalP(encP, encXBeta, logp, loga, kdeg);
	encGWAS.evalW(encW, encP, logp);
	timeutils.stop("step [2]");

	
	// readcsv(scheme, secretKey, encW, "outW");

	
	SerializationUtils::writeCiphertext(&encXBeta, folder + "encXBeta.txt");
	SerializationUtils::writeCiphertext(&encBeta, folder + "encBeta.txt");
	SerializationUtils::writeCiphertext(&encP, folder + "encP.txt");
//	delete &encXBeta; delete &encBeta; delete &encP; 
	encXBeta.kill();
	encBeta.kill();
	encP.kill();
	for(long i = 0; i < 4; i++){
		encDXvec[i].kill();
	}
	delete[] encDXvec;

	encU = *SerializationUtils::readCiphertext(folder + "encU.txt");
	encX = *SerializationUtils::readCiphertext(folder + "encX.txt");
	encAdjU = *SerializationUtils::readCiphertext(folder + "encAdjU.txt");

	encDXExpvec = new Ciphertext[4];
	for (long i = 0; i < 4; ++i) {
		encDXExpvec[i] = *SerializationUtils::readCiphertext(folder + "encDXExpvec" + to_string(i) + ".txt");
	}

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [3] & [4] - Eval U, detU, AdjU");
	encGWAS.evalU(encU, encX, encW, encDXExpvec, n2, logp);
	encGWAS.evalAdjU(encAdjU, encU, maskAdjU, n2, logp, logc);
	encGWAS.evalDetU(encDetU, encU, encAdjU, maskDetU, n2, logp, logc);
	encGWAS.rearrange(encU, maskRU, n2, logp, logc);
	encGWAS.rearrange(encAdjU, maskRU, n2, logp, logc);

	readcsv(scheme, secretKey, encAdjU, "adjU_before_BT");

	encAdjU.n = 16;
	scheme.reScaleByAndEqual(encAdjU, logp - logbp);
	scheme.bootstrapAndEqual(encAdjU, logbt, logQ, logbT, logbI);
	scheme.leftShiftAndEqual(encAdjU, logp - logbp);
	encAdjU.n = n2;
	encAdjU.logp = logp;

	readcsv(scheme, secretKey, encAdjU, "adjU_after_BT");

	encAdjU.n = (1 << (logN - 1));
	for (long i = 4; i < logN - 1; ++i) {
		Ciphertext rot = scheme.leftRotateFast(encAdjU, (1 << i));
		scheme.addAndEqual(encAdjU, rot);
		rot.kill();
	}
	encAdjU.n = n2;
	scheme.divByPo2AndEqual(encAdjU, logN - 5);

	// readcsv(scheme, secretKey, encAdjU, "adjU_after_BT+average");


	encDetU.n = 1;
	scheme.reScaleByAndEqual(encDetU, logp - logbp);
	scheme.bootstrapAndEqual(encDetU, logbt, logQ, logbT, logbI);
	scheme.leftShiftAndEqual(encDetU, logp - logbp);
	encDetU.n = n2;
	encDetU.logp = logp;

	readcsv(scheme, secretKey, encDetU, "DetU");
	encDetU.n = n2;
	encDetU.logp = logp;

	timeutils.stop("step [3] & [4]");

	SerializationUtils::writeCiphertext(&encU, folder + "encU.txt");
	SerializationUtils::writeCiphertext(&encX, folder + "encX.txt");
	SerializationUtils::writeCiphertext(&encW, folder + "encW.txt");
	SerializationUtils::writeCiphertext(&encAdjU, folder + "encAdjU.txt");
	SerializationUtils::writeCiphertext(&encDetU, folder + "encDetU.txt");
//	delete &encU; delete &encX; delete &encW; delete &encAdjU; delete &encDetU;
	encU.kill();
	encX.kill();
	encW.kill();
	encAdjU.kill();
	encDetU.kill();

	for (long i = 0; i < 4; ++i) {
		SerializationUtils::writeCiphertext(&encDXExpvec[i], folder + "encDXExpvec" + to_string(i) + ".txt");
		encDXExpvec[i].kill();
	}
	delete[] encDXExpvec;

	encY = *SerializationUtils::readCiphertext(folder + "encY.txt");
	encP = *SerializationUtils::readCiphertext(folder + "encP.txt");
	encSzvec = new Ciphertext[packctxNum];
	encSvec = new Ciphertext[sNum];
	for (long i = 0; i < sNum; ++i) {
		encSvec[i] = *SerializationUtils::readCiphertext(folder + "encSvec" + to_string(i) + ".txt");
	}
	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [5] : Compute S^t (y-p)");
	encGWAS.evalSWz(encSzvec, encSvec, encY, encP, maskPack, n2, packNum, lastpackNum, packctxNum, logp, logc);
	timeutils.stop("step [5]");

	// complex<double>* dectest;
	// cout << "encSz " << endl;
	// for (int i = 0; i < packctxNum; i++) {
	// 	dectest = scheme.decrypt(secretKey, encSzvec[i]);
	// 	for (int j = 0; j < (n2 >> 2); j++) {
	// 		for (int k = 0; k < ScolNum; k++) {
	// 			cout << (1 << (logN-3))*i + j*ScolNum + k << ": " << dectest[n2*k + 4*j].real() <<", ";			
	// 		}
	// 		cout << endl;
	// 	}
	// }

	SerializationUtils::writeCiphertext(&encY, folder + "encY.txt");
	SerializationUtils::writeCiphertext(&encP, folder + "encP.txt");
//	delete &encY; delete &encP;
	encY.kill();
	encP.kill();

	for (long i = 0; i < packctxNum; ++i) {
		SerializationUtils::writeCiphertext(&encSzvec[i], folder + "encSzvec" + to_string(i) + ".txt");
		encSzvec[i].kill();
	}
	delete[] encSzvec;


	encW = *SerializationUtils::readCiphertext(folder + "encW.txt");

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [6] : Convert encS into encWS");
	encGWAS.evalWS(encW, encSvec, sNum, logp);
	timeutils.stop("step [6]");

	SerializationUtils::writeCiphertext(&encW, folder + "encW.txt");
	//delete &encW;
	encW.kill();

	encDSWSvec = new Ciphertext[packctxNum];

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [7] : Compute diag(S_i^t W S_i) and Repack");
	encGWAS.evalDSWS(encDSWSvec, encSvec, maskPack, n2, packNum, lastpackNum, packctxNum, logp, logc);
	timeutils.stop("step [7]");

	for (long i = 0; i < packctxNum; ++i) {
		SerializationUtils::writeCiphertext(&encDSWSvec[i], folder + "encDSWSvec" + to_string(i) + ".txt");
	}
	delete[] encDSWSvec;

	encVvec = new Ciphertext[packctxNum];
	encDXExpvec = new Ciphertext[4];
	for (long i = 0; i < 4; ++i) {
		encDXExpvec[i] = *SerializationUtils::readCiphertext(folder + "encDXExpvec" + to_string(i) + ".txt");
	}

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [8] : Compute V = X^t W S_i and Repack into X^t W S");
	encGWAS.evalV(encVvec, encSvec, encDXExpvec, maskXWS, maskPack, n2, packNum, lastpackNum, packctxNum, logp, logc);
	timeutils.stop("step [8]");

	/*double* vtrans = new double[2048 * 32]();
	encVvec[0].n = 2048 * 32;
	complex<double>* tmp = scheme.decrypt(secretKey, encVvec[0]);
	int r = 0;
	for (int i = 0; i < 4; i++) {
		for (int s = 0; s < 64; s++) {
			for (int j = 0; j < 256; j++) {
				vtrans[r] = tmp[i + 256 * j + 4 * s].real();
				r++;
			}
		}
	}
	ofstream Vtrans("outVtrans.csv");
	for (int i = 0; i < 2048 * 32; i++) {
		Vtrans << *vtrans;
		if (i < 2048 * 32 - 1)
			Vtrans << "\n";
		vtrans++;
	}*/
	
	/*scheme.multAndEqual(encVvec[0], encVvec[0]);		
	scheme.reScaleByAndEqual(encVvec[0], logp);
	cout << "out V squares to outV2trans.csv." << endl;


	encVvec[0].n = 2048 * 32;
	tmp = scheme.decrypt(secretKey, encVvec[0]);
	r = 0;
	for (int i = 0; i < 4; i++) {
		for (int s = 0; s < 64; s++) {
			for (int j = 0; j < 256; j++) {
				vtrans[r] = tmp[i + 256 * j + 4 * s].real();
					r++;
			}
		}
	}
	ofstream Vtrans2("outV2trans.csv");
	for (int i = 0; i < 2048 * 32; i++) {
		Vtrans2 << *vtrans;
		if (i < 2048 * 32 - 1)
			Vtrans2 << "\n";
		vtrans++;
	}*/

	/*double* dvvec = new double[m]();
	int	j = 0;
	for (int i = 0; i < packctxNum; i++) {
		complex<double>* tmpvvec = scheme.decrypt(secretKey, encVvec[i]);
		for (int s = 0; s < 64; s++) {
			for (int r = 0; r < 8; r++) {
				dvvec[j] = tmpvvec[256 * r + 4 * s].real();
				j++;
			}
		}
	}
	ofstream DVV("outVVec.csv");
	for (int i = 0; i < 2048; ++i) {
		DVV << *dvvec;
		if (i < 2048 - 1)
			DVV << "\n";
		dvvec++;
	}*/
	
	for (long i = 0; i < 4; ++i) {
		SerializationUtils::writeCiphertext(&encDXExpvec[i], folder + "encDXExpvec" + to_string(i) + ".txt");
		encDXExpvec[i].kill();
	}
	delete[] encDXExpvec;
	for(long i = 0; i < sNum; i++){
		encSvec[i].kill();
	}
	delete[] encSvec;

	encAdjU = *SerializationUtils::readCiphertext(folder + "encAdjU.txt");

	encDUExpvec = new Ciphertext[4];
	encDVUVvec = new Ciphertext[packctxNum];

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [9] : Compute diag(V^t adjU V)");
	encGWAS.evalDAdjUExp(encDUExpvec, encAdjU, maskDU, logc);
	encGWAS.evalDVAdjUV(encDVUVvec, encVvec, encDUExpvec, maskDVUV, packctxNum, logp, logc);
	timeutils.stop("step [9]");

	//for(int i = 0; i < 4; i++){
	//	//readcsv(scheme, secretKey, encDUExpvec[i], "AdjU_diagexpand" + to_string(i));	
	//	
	//	double* vuv_part = new double[2048 * 32]();
	//	complex<double>* tmp1 = scheme.decrypt(secretKey, encDUExpvec[i]);
	//	int r = 0;
	//	for (int w = 0; w < 4; w++) {
	//		for (int s = 0; s < 64; s++) {
	//			for (int x = 0; x < 256; x++) {
	//				vuv_part[r] = tmp1[w + 256 * x + 4 * s].real();
	//				r++;
	//			}
	//		}
	//	}
	//	ofstream VUVpart("AdjU_diaexp.csv" + to_string(i));
	//	for (int i = 0; i < 2048 * 32; i++) {
	//		VUVpart << *vuv_part;
	//		if (i < 2048 * 32 - 1)
	//			VUVpart << "\n";
	//		vuv_part++;
	//	}
	//}
	

	
//	readcsv(scheme, secretKey, encDVUVvec[0], "outDVUV");


	/*double* dvuv = new double[m*32]();
	int	j = 0;
	for (int i = 0; i < packctxNum; i++) {
		encDVUVvec[i].n = 2048*32;
		complex<double>* tmpdvuv = scheme.decrypt(secretKey, encDVUVvec[i]);
		for (int s = 0; s < 64; s++) {
			for (int r = 0; r < 8*32; r++) {
				dvuv[j] = tmpdvuv[256 * r + 4 * s].real();
				j++;
			}
		}
	}
	ofstream DVUV("outVUV.csv");
	for (int i = 0; i < 2048*32; ++i) {
		DVUV << *dvuv;
		if (i < 2048*32 - 1)
			DVUV << "\n";
		dvuv++;
	}*/

	//	delete &encAdjU;
	encAdjU.kill();
	for(long i = 0; i < packctxNum; i++){
		encVvec[i].kill();
	}
	delete[] encVvec;
	for(long i = 0; i < 4; i++){
		encDUExpvec[i].kill();
	}
	delete[] encDUExpvec;

	encDetU = *SerializationUtils::readCiphertext(folder + "encDetU.txt");
	encIPSvec = new Ciphertext[packctxNum];
	encDSWSvec = new Ciphertext[packctxNum];

	for (long i = 0; i < packctxNum; ++i) {
		encDSWSvec[i] = *SerializationUtils::readCiphertext(folder + "encDSWSvec" + to_string(i) + ".txt");
	}


	/*double* dsws = new double[m*32]();
	j = 0;
	for (int i = 0; i < packctxNum; i++) {
		complex<double>* tmpdsws = scheme.decrypt(secretKey, encDSWSvec[i]);	
		for (int s = 0; s < 64; s++) {
			for (int r = 0; r < 8*32; r++) {
				dsws[j] = tmpdsws[256 * r + 4 * s].real();
				j++;
			}
		}
	}
	ofstream DSWS("outSWS.csv");
	for (int i = 0; i < 2048*32; ++i) {
	DSWS << *dsws;
	if (i < 2048*32 - 1)
	DSWS << "\n";
	dsws++;
	}*/


//	readcsv(scheme, secretKey, encDSWSvec[0], "outDSWS");


	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [10] : (detU * [7]) - [9] = det * diag(S*^t W S*) ");
	encGWAS.evalIPS(encIPSvec, encDSWSvec, encDVUVvec, encDetU, maskIPS, packctxNum, logp, logc);
	timeutils.stop("step [10]");

	// cout << "encIPS " << endl;
	// for (int i = 0; i < packctxNum; i++) {
	// 	dectest = scheme.decrypt(secretKey, encIPSvec[i]);
	// 	for (int j = 0; j < (n2 >> 2); j++) {
	// 		for (int k = 0; k < ScolNum; k++) {
	// 			cout << (1 << (logN-3))*i + j*ScolNum + k << ": " << dectest[n2*k + 4*j].real() <<", ";			
	// 		}
	// 		cout << endl;
	// 	}
	// }


// 	SerializationUtils::writeCiphertext(&encDetU, folder + "encDetU.txt");
// //	delete &encDetU;
// 	encDetU.kill();
	for(long i = 0; i < packctxNum; i++){
		encDSWSvec[i].kill();
		encDVUVvec[i].kill();
	}
	delete[] encDSWSvec;
	delete[] encDVUVvec;

	encSzvec = new Ciphertext[packctxNum];
	for (long i = 0; i < packctxNum; ++i) {
		encSzvec[i] = *SerializationUtils::readCiphertext(folder + "encSzvec" + to_string(i) + ".txt");
	}

	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	timeutils.start("step [11] : Calculate det * [5]^2");
	Ciphertext* encnumerator = new Ciphertext[packctxNum];
	for(int i = 0; i < packctxNum; i++){
		encnumerator[i] = scheme.mult(encSzvec[i], encSzvec[i]);
		scheme.reScaleByAndEqual(encnumerator[i], logp);
		scheme.modDownToAndEqual(encDetU, encnumerator[i].logq);
		scheme.multAndEqual(encnumerator[i], encDetU);
		scheme.reScaleByAndEqual(encnumerator[i], logp);
		}
	timeutils.stop("step [11]");

	cout <<"End of the Encrpyted Semi-Parallel GWAS Computations!" << endl;
	cout <<"Final logQ: " << encnumerator[0].logq << endl;
	cout <<"Final logQ of IPS: " << encIPSvec[0].logq << endl;


	
	//-----------------------------------------
	currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	peakAfterSchemeSize = getPeakRSS() >> 20;
	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;
	//-----------------------------------------

	long m2 = sNum * ScolNum;
	timeutils.start("decryption");
	IPS = new double[m2]();
	for (long k = 0; k < packctxNum; ++k) {
		long tmppackNum = (k == packctxNum - 1) ? lastpackNum : packNum;
		encIPSvec[k].n = n2 * ScolNum;
		complex<double>* dIPS = scheme.decrypt(secretKey, encIPSvec[k]);
		for (long i = 0; i < tmppackNum; ++i) {
			for (long j = 0; j < ScolNum; ++j) {
				IPS[k * ScolNum * packNum + i * ScolNum + j] = dIPS[(i << 2) + j * n2].real();
			}
		}
		delete[] dIPS;
	}

	numerator = new double[m2]();
	for (long k = 0; k < packctxNum; ++k) {
		long tmppackNum = (k == packctxNum - 1) ? lastpackNum : packNum;
		encnumerator[k].n = n2 * ScolNum;
		complex<double>* dnum = scheme.decrypt(secretKey, encnumerator[k]);
		for (long i = 0; i < tmppackNum; ++i) {
			for (long j = 0; j < ScolNum; ++j) {
				numerator[k * ScolNum * packNum + i * ScolNum + j] = dnum[(i << 2) + j * n2].real();
			}
		}
		delete[] dnum;
	}
	timeutils.stop("decryption");

	cout <<"NOTE: the square of Statistics = numerator / IPS !!" << endl;

	cout << "========================" << endl;
	cout << "End Test Enc Separate" << endl;
	cout << "========================" << endl;
}

void TestGWAS::readcsv(Scheme& scheme, SecretKey& secretKey, Ciphertext& cipher, string path) {
	complex<double>* check = scheme.decrypt(secretKey, cipher);
	//double* checkd = new double[cipher.slots]();

	ofstream outfile(path + ".csv");
	for (int i = 0; i < cipher.n; ++i) {
		outfile << (*check).real();
		if (i < cipher.n - 1)
			outfile << "\n";
		check++;
	}
}