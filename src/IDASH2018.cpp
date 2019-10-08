/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
:*/

#include <iostream>
#include <math.h>
#include "TestGWAS.h"
#include "Extractor.h"

using namespace std;

// ./IDASH2018 ../data/covariates.csv ../data/snpMat.txt 17 1500 45 20 30 40 3 4 4 7 1
int main(int argc, char **argv) {

	string filenameXY(argv[1]);
	string filenameS(argv[2]);

	//scheme parameters
	long logN = atoi(argv[3]); // 17
	long logQ = atoi(argv[4]); // 1500
	long logp = atoi(argv[5]); // 45
	long logc = atoi(argv[6]); // 20

	//scheme parameters for bootstrapping
	long logbp = atoi(argv[7]); // 30
	long logbt = atoi(argv[8]); // 40
	long logbT = atoi(argv[9]); // 3
	long logbI = atoi(argv[10]); // 4

	//other parameters
	long fnum = atoi(argv[11]); // 4
	long kdeg = atoi(argv[12]); // 7
	long logscale = atoi(argv[13]); // 1

	long n, n2, m;
	double **X, **S;
	double *y, *dStatsq, *dIPS, *dInvIPS, *dSz;

	Extractor::readXY(filenameXY, X, y, n, n2);
//	m = 2048;
//	Extractor::readSfixm(filenameS, S, n, m, n2);
	Extractor::readS(filenameS, S, n, m, n2);
	Extractor::normalizeGlobal(X, n);

	cout << "extractor complete" << endl;
	cout << "n:" << n << ", n2:" << n2 << ", m:" << m << endl; 
// TestGWAS::testEnc(dStatsq, dIPS, dInvIPS, dSz, X, y, S, n2, m, logN, logQ, logp, logc, logbp, logbt, logbT, logbI, fnum, kdeg, logscale);

// cout << "enc Statistics" << endl;
// double* dStat = new double[m];
// for(long i = 0; i < m; ++i){
// 	if(dStatsq[i] > 0)
// 		dStat[i] = sqrt(dStatsq[i]);
// 	else
// 		dStat[i] = 0;
// }
// for (long i = 0; i < m; ++i) {
// 	cout << i << " " << dStat[i] << ";";
// 	if((i+1) % 8 == 0) cout << endl;
// }
// cout << endl;

// ofstream output("Output.csv");
// for(int i = 0; i < m; i++){
// 	output << *dStat;
//     	if (i < m - 1)
//	       	output << "\n";
//   	dStat++;
//    }


	// If the result of testEnc explodes, use following function which provides numerator and denominator respectively.	

	double *dnumerator;
	TestGWAS::testEnc_seperate(dnumerator, dIPS, X, y, S, n2, m, logN, logQ, logp, logc, logbp, logbt, logbT, logbI, fnum, kdeg);
	
	cout << "enc Statistics" << endl;
	double* dStat = new double[m];
	for(long i = 0; i < m; ++i){
		double stat = dnumerator[i] / dIPS[i];
	 	if(stat > 0)
			dStat[i] = sqrt(stat);
		else
			dStat[i] = 0;
	}

	// for (long i = 0; i < m; ++i) {
	// 	cout << i << " " << dnumerator[i] << ";";
	// 	if ((i+1) % 8 == 0) cout << endl;
	// }
	
	// for (long i = 0; i < m; ++i) {
	// 	cout << i << " " << dIPS[i] << ";";
	// 	if ((i+1) % 8 == 0) cout << endl;
	// }

	//	for (long i = 0; i < m; ++i) {
	//		cout << i << " " << dStat[i] << ";";
	//		if((i+1) % 8 == 0) cout << endl;
	//	}
	//	cout << endl;	
	 ofstream output("Output.csv");
	 for(int i = 0; i < m; i++){
		output << *dStat;
	      	if (i < m - 1)
	        	output << "\n";
	      	dStat++;
	}

	ofstream outnum("Outnum.csv");
	for (int i = 0; i < m; ++i){
		outnum << *dnumerator;
		if (i < m - 1)
			outnum << "\n";
		dnumerator++;
	}

	ofstream outde("Outdenom.csv");
	for (int i = 0; i < m; ++i){
		outde << *dIPS;
		if (i < m - 1)
			outde << "\n";
		dIPS++;
	}


	return 0;
}

