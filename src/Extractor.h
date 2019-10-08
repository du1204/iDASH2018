/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef EXTRACTOR_H_
#define EXTRACTOR_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

class Extractor {

public:
	static void readXY(string filename, double**& X, double*& y, long& n, long& n_pw2);

	static void normalize(double** X, long n);

	static void normalizeGlobal(double** X, long n);

	static void normalizeGlobalModified(double** X, long n);

	static void readS(string filename, double**& S, long& n, long& m, long& n2);

	static void readSfixm(string filename, double**& S, long& n, long& m, long& n2);

};

#endif /* EXTRACTOR_H_ */
