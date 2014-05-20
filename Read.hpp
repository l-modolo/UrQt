/*
polyNtrimmer poly nucleotide trimming tool
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

#ifndef DEF_Read
#define DEF_Read

#include <zlib.h>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <string>
#include <queue>
#include <iostream>
#include <fstream>
#include <ostream>
#include <stdlib.h> 
#include <cctype>
#include <mutex>

using namespace std;

class Read
{
	public:
	Read(gzFile* fin, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number, bool remove_empty_reads, bool fill_empty_reads, char* fill_empty_reads_with, int min_QC_phred, double min_QC_length);
	Read(gzFile* fin, char* out, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number, bool remove_empty_reads, bool fill_empty_reads, char* fill_empty_reads_with, int min_QC_phred, double min_QC_length, int paired);
	Read(gzFile* fin, char* out, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number, bool remove_empty_reads, bool fill_empty_reads, char* fill_empty_reads_with, int min_QC_phred, double min_QC_length, bool estimation, int paired);
	void constructor(gzFile* fin,char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number, bool remove_empty_reads, bool fill_empty_reads, char* fill_empty_reads_with, int min_QC_phred, double min_QC_length);
	~Read();
	Read& operator=(Read const& readbis);
	void operator()();
	void init(gzFile* fin);
	void polyNtrim();
	void polyNtrimEstimate();
	
	void writeRead(ostream &stream, int* base_trimmed);
	string writeRead(int* base_trimmed);
	void writeRead();
	int writeSize();
	void run();
	void done();
	
	static double G_content();
	static double C_content();
	static double A_content();
	static double T_content();
	static bool sampling_done();

	static int empty_reads();
	static int trimmed_reads();
	static void reset();

	static void close();

	static void remove_empty_reads_paired(char* out, char* outpaired);

	private:
	bool m_init;
	int m_paired;
	bool m_estimation;
	bool m_sampled;
	bool m_trimmed;
	bool m_remove_empty_reads;
	bool m_fill_empty_reads;
	bool m_QC_check;
	int m_size;
	int m_cut;
	int m_read_number;
	int m_start;
	int m_stop;
	int m_min_QC_phred;
	string m_name;
	string m_seq;
	string m_phred;
	double* m_proba;
	double m_log_read;
	double m_log_polyN;
	double m_min_QC_length;
	char m_N;
	char m_fill_empty_reads_with;
	

	static double m_G_number;
	static double m_C_number;
	static double m_A_number;
	static double m_T_number;

	static double m_G_probability;
	static double m_C_probability;
	static double m_A_probability;
	static double m_T_probability;

	static queue<int> m_paired_pos_1;
	static queue<int> m_paired_pos_2;

	static bool m_out_open;
	static ofstream m_out;
	static bool m_phred_score_set;
	static int m_phred_score;
	static int m_empty_reads;
	static int m_trimmed_reads;
	
	static mutex m_read;
	
	// compute the probability to have the right base
	void phred();
	
	// in a read the probability of observing one of the four bases is roughly the same
	// for each base we use the probability of being the right base
	double read(int begin, int end);
	double readEstimate(int begin, int end, double p_G, double p_C, double p_A, double p_T);

	// in a polyN segment of size n we expect to observe N with with a probability 1/n
	// for each N we use the probability of being the right base
	// for each non-N we use 1 - the probability of being the right base
	double polyN(int begin, int end);

	double baseProba(double &pG, double &pC, double &pA, double &pT);
	void baseProbaTotal();

	double probaBaseDict(char base, double proba, double pG, double pC, double pA, double pT);
	void numberBaseDict(char base, double proba, double &G_number, double &C_number, double &A_number, double &T_number);

	// perform a quality check in the read, if it was constructed with the correct argument 
	bool QC_check();
};

#endif
