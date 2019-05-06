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

#ifndef DEF_Read
#define DEF_Read

# define BUFFER_LENGTH 1024000

#include <zlib.h>
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <string>
#include <queue>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <stdlib.h>
#include <cctype>
#include <mutex>
#include "ezRateProgressBar.hpp"
#include "gzstream.hpp"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "Segmentation.hpp"

using namespace std;

class Segmentation;

class Read
{
public:
	Read(igzstream &fin, char* N, int phred_score, int threshold, int max_head_trim, int max_tail_trim, int min_read_size, int read_number, bool remove_empty_reads, int min_QC_phred, double min_QC_length, int strand_bit);
	Read(igzstream &fin, char* out, bool gziped, char* N, int phred_score, int threshold, int max_head_trim, int max_tail_trim, int min_read_size, int read_number, bool remove_empty_reads, int min_QC_phred, double min_QC_length, int paired, int strand_bit);
	Read(igzstream &fin, char* out, bool gziped, char* N, int phred_score, int threshold, int max_head_trim, int max_tail_trim, int min_read_size, int read_number, bool remove_empty_reads, int min_QC_phred, double min_QC_length, bool estimation, int paired, int strand_bit);
	inline void constructor(igzstream &fin,char* N, int phred_score, int threshold, int max_head_trim, int max_tail_trim, int min_read_size, int read_number, bool remove_empty_reads, int min_QC_phred, double min_QC_length, int strand_bit);
	~Read();
	Read& operator=(Read const& readbis);
	void operator()();
	void init(igzstream &fin);

	inline void writeRead(ostream &stream);
	inline void writeRead();
	int writeSize();
	void run();
	void done();

	char base();
	int size();
	int start();
	int stop();
	int max_head_trim();
	int max_tail_trim();
	int min_read_size();
	char seq(int i);
	char phred(int i);
	int min_QC_phred();
	int threshold();
	double min_QC_length();
	static double G_probability();
	static double C_probability();
	static double A_probability();
	static double T_probability();
	int strand();
	void set_trim(bool trim, int cut_begin, int cut_end);
	bool QC_check();

	static double G_content();
	static double C_content();
	static double A_content();
	static double T_content();
	static bool sampling_done();
	static string base_trimmed();


	static int empty_reads();
	static int trimmed_reads();
	static void reset();

	static void open(char* out);
	static void close();

	static void remove_empty_reads_paired(char* out, char* outpaired, ez::ezRateProgressBar<int>* p, bool v);

private:
	bool m_init;
	bool m_estimation;
	bool m_sampled;
	bool m_trimmed;
	bool m_remove_empty_reads;
	bool m_QC_check;
	int m_size;
	int m_cut_begin;
	int m_cut_end;
	int m_read_number;
	int m_start;
	int m_stop;
	int m_max_head_trim;
	int m_max_tail_trim;
	int m_min_read_size;
	int m_min_QC_phred;
	int m_strand;
	int m_threshold;
	string m_name;
	string m_seq;
	string m_phred;
	double m_log_read;
	double m_log_polyN;
	double m_min_QC_length;
	char m_N;

	static double m_G_number;
	static double m_C_number;
	static double m_A_number;
	static double m_T_number;
	static long long m_base_number;
	static long long m_base_trimmed;

	static double m_G_probability;
	static double m_C_probability;
	static double m_A_probability;
	static double m_T_probability;

	static queue<int> m_paired_pos_1;
	static queue<int> m_paired_pos_2;

	static int m_paired;
	static bool m_out_open;
	static bool m_gziped;
	static ofstream m_out;
	static ogzstream m_out_gz;
	static bool m_phred_score_set;
	static char m_buffer[BUFFER_LENGTH];
	static int m_phred_score;
	static long long m_empty_reads;
	static long long m_trimmed_reads;

	static mutex m_read;
	static void writeFinal(queue<int> m_paired_pos, char* out, ez::ezRateProgressBar<int>* p, bool v);

};

#endif
