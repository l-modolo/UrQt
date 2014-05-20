/*
UrQt quality and poly nucleotide trimming tool
Copyright (C) 2013  Laurent Modolo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DEF_Segmentation
#define DEF_Segmentation

#include "Read.hpp"
using namespace std;

class Read;

class Segmentation
{
public:
	Segmentation(Read* const read, bool estimate);
	Segmentation(Read* const read, double &pG, double &pC, double &pA, double &pT, bool trimmed, int cut_begin, int cut_end);
	~Segmentation();

	void polyNtrim();
	void polyNtrimEstimate();

private:
	bool m_init;
	Read* m_read;
	double* m_proba;
	int m_cut_begin = -1;
	int m_cut_end = -1;
	double m_log_read;
	double m_log_polyN;
	char m_N;
	int m_min_QC_phred;
	double m_min_QC_length;
	double m_G_probability;
	double m_C_probability;
	double m_A_probability;
	double m_T_probability;

	// compute the probability to have the right base
	void phred();
	
	// in a read the probability of observing one of the four bases is roughly the same
	// for each base we use the probability of being the right base
	double read(int begin, int end);
	double readEstimate(int begin, int end, double p_G, double p_C, double p_A, double p_T);
	double reverseRead(int begin, int end);
	double reverseReadEstimate(int begin, int end, double p_G, double p_C, double p_A, double p_T);

	// in a polyN segment of size n we expect to observe N with with a probability 1/n
	// for each N we use the probability of being the right base
	// for each non-N we use 1 - the probability of being the right base
	double polyN(int begin, int end);
	double reversePolyN(int begin, int end);

	double baseProba(double &pG, double &pC, double &pA, double &pT);
	void baseProbaTotal(double &pG, double &pC, double &pA, double &pT, bool trimmed, int cut_begin, int cut_end);

	double probaBaseDict(char base, double proba, double pG, double pC, double pA, double pT);
	void numberBaseDict(char base, double proba, double &G_number, double &C_number, double &A_number, double &T_number);

	// perform a quality check in the read, if it was constructed with the correct argument 
	bool QC_check();

};

#endif