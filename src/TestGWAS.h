/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef TESTGWAS_H_
#define TESTGWAS_H_
#include "EncGWAS.h"
#include <iostream>

using namespace std;

class TestGWAS {
public:
	static void readcsv(Scheme& scheme, SecretKey& secretKey, Ciphertext& cipher, string path);

	static void testEnc(double*& Statsq, double*& IPS, double*& InvIPS, double*& Sz, double** X, double* y, double** S, long n2, long m, long logN,
			long logQ, long logp, long logc, long logbp, long logbt, long logbT, long logbI, long fnum, long kdeg, long logscale);

	static void testEnc_seperate(double*& numerator, double*& IPS, double** X, double* y, double** S, long n2, long m, long logN,
			long logQ, long logp, long logc, long logbp, long logbt, long logbT, long logbI, long fnum, long kdeg);
	
};

#endif /* TESTGWAS_H_ */
