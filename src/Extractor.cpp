/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Extractor.h"

#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>    // std::shuffle
#include <array>        // std::array
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

void Extractor::readXY(string filename, double**& X, double*& y, long& n, long& n_pw2) {
	ifstream openFile(filename.data());
	vector< vector<double> > Xydata;
	n = 0;
	int* notNA;
	if(openFile.is_open()) {
		string line, temp;
		getline(openFile, line);
		notNA = new int[4]();
		size_t start, end;
		while(getline(openFile, line)){
			vector<double> Xline;
			end = line.find_first_of(',');
			for (long i = 0; i < 4; ++i) {
				start = end + 1;
				end = line.find_first_of (',', start);
				temp = line.substr(start, end - start);
				if(temp == "NA") {
					Xline.push_back(-1);
				} else {
					Xline.push_back(atof(temp.c_str()));
					notNA[i] += 1;
				}
			}
			Xydata.push_back(Xline);
			n++;
		}
	} else {
		cout << "Error: cannot read file" << endl;
	}
	long n_bits = (long) ceil(log2(n));
	n_pw2 = 1 << n_bits;
	X = new double*[n_pw2];
	for (long i = 0; i < n; ++i) {
		X[i] = new double[4]();
		X[i][0] = 1;
	}
	for (long i = n; i < n_pw2; ++i){
		X[i] = new double[4]();
	}

	y = new double[n_pw2]();


	if(notNA[0] < n) {
		double mean = 0.0;
		for (long j = 0; j < n; ++j) {
			if(Xydata[0][j] > -1) {
				mean += Xydata[j][0];
				y[j] = Xydata[j][0];
			}
		}
		mean /= notNA[0];
		for (long j = 0; j < n; ++j) {
			if(Xydata[j][0] < 0) {
				y[j] = mean;
			}
		}
	} else {
		for (long j = 0; j < n; ++j) {
			y[j] = Xydata[j][0];
		}
	}
	for (long i = 1; i < 4; ++i) {
		if(notNA[i] < n) {
			double sum = 0;
			for (long j = 0; j < n; ++j) {
				if(Xydata[j][i] >= 0) {
					sum += Xydata[j][i];
					X[j][i] = Xydata[j][i];
				}
			}
			sum /= notNA[i];
			for (long j = 0; j < n; ++j) {
				if(Xydata[j][i] < 0) {
					X[j][i] = sum;
				}
			}
		} else {
			for (long j = 0; j < n; ++j) {
				X[j][i] = Xydata[j][i];
			}
		}
	}

	delete[] notNA;
}

void Extractor::normalize(double** X, long n) {
	long i, j;
	double m;

	for (i = 1; i < 4; ++i) {
		m = X[0][1];
		for (j = 1; j < n; ++j) {
			m = max(m, X[j][i]);
		}
		for (j = 0; j < n; ++j) {
			X[j][i] = X[j][i] / m;
		}
	}
}

void Extractor::normalizeGlobal(double** X, long n) {
	long i, j;
	double m, s;
	m = X[0][1];
	s = X[0][1];
	for (i = 1; i < 4; ++i) {
		for (j = 0; j < n; ++j) {
			m = max(m, X[j][i]);
			s = min(s, X[j][i]);
		}
	}
	for (long i = 1; i < 4; ++i) {
		for (j = 0; j < n; ++j) {
			X[j][i] = (X[j][i] - s) / (m-s);
		}
	}
}

void Extractor::normalizeGlobalModified(double** X, long n) {
	long i, j;
	double m, s;

	for (i = 1; i < 4; ++i) {
		m = X[0][i];
		s = X[0][i];
		for (j = 0; j < n; ++j) {
			m = max(m, X[j][i]);
			s = min(s, X[j][i]);
		}
		for (j = 0; j < n; ++j) {
			X[j][i] = (X[j][i] - s) / (m - s);
		}
	}
}



void Extractor::readS(string filename, double**& S, long& n, long& m, long& n2) {
	vector<vector<double>> Sdata;
	m = 0;
	ifstream openFile(filename.data());
	if(openFile.is_open()) {
		string line, temp;
		getline(openFile, line);
		long i;
		size_t start, end;
		for(i = 0; i < line.length(); ++i) if(line[i] == ' ' ) m++;
		while(getline(openFile, line)){
			vector<double> vecline;
			do {
				end = line.find_first_of (' ', start);
				temp = line.substr(start,end - start);
				vecline.push_back(atof(temp.c_str()));
				start = end + 1;
			} while(start);
			Sdata.push_back(vecline);
		}
	} else {
		cout << "Error: cannot read file" << endl;
	}
	S = new double*[n2];
	for(long j = 0; j < n; ++j){
		S[j] = new double[m]();
		for (long i = 0; i < m; ++i) {
			S[j][i] = Sdata[j][i];
		}
	}
	for(long j = n; j < n2; ++j){
		S[j] = new double[m]();
	}
}

void Extractor::readSfixm(string filename, double**& S, long& n, long& m, long& n2) {
	vector<vector<double>> Sdata;
	ifstream openFile(filename.data());
	if(openFile.is_open()) {
		string line, temp;
		getline(openFile, line);
		size_t start, end;
		while(getline(openFile, line)){
			vector<double> vecline;
			do {
				end = line.find_first_of (' ', start);
				temp = line.substr(start,end - start);
				vecline.push_back(atof(temp.c_str()));
				start = end + 1;
			} while(start);
			Sdata.push_back(vecline);
		}
	} else {
		cout << "Error: cannot read file" << endl;
	}

	S = new double*[n2];
	for(long j = 0; j < n; ++j){
		S[j] = new double[m]();
		for (long i = 0; i < m; ++i) {
			S[j][i] = Sdata[j][i];
		}
	}
	for(long j = n; j < n2; ++j){
		S[j] = new double[m]();
	}
}
