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

#include "Read.hpp"

// Definition of the static stuff
double Read::m_G_number = 0.0;
double Read::m_C_number = 0.0;
double Read::m_A_number = 0.0;
double Read::m_T_number = 0.0;
long long Read::m_base_number = 0;
long long Read::m_base_trimmed = 0;

double Read::m_G_probability = 0.25;
double Read::m_C_probability = 0.25;
double Read::m_A_probability = 0.25;
double Read::m_T_probability = 0.25;

queue<int> Read::m_paired_pos_1;
queue<int> Read::m_paired_pos_2;

mutex Read::m_read;
int Read::m_paired = -1;
bool Read::m_out_open = false;
bool Read::m_gziped = false;
ofstream Read::m_out;
ogzstream Read::m_out_gz;
char Read::m_buffer[BUFFER_LENGTH];
bool Read::m_phred_score_set = false;
int Read::m_phred_score = 33.0;
long long Read::m_empty_reads = 0;
long long Read::m_trimmed_reads = 0;

bool Read::sampling_done()
{
	if(m_G_number != 0.0 && m_C_number != 0.0 && m_A_number != 0.0 && m_T_number != 0.0)
		return true;
	return false;
}

double Read::G_content()
{
	m_G_probability = m_G_number/(m_G_number + m_C_number + m_A_number + m_T_number);
	return m_G_probability;
}
double Read::C_content()
{
	m_C_probability = m_C_number/(m_G_number + m_C_number + m_A_number + m_T_number);
	return m_C_probability;
}
double Read::A_content()
{
	m_A_probability = m_A_number/(m_G_number + m_C_number + m_A_number + m_T_number);
	return m_A_probability;
}
double Read::T_content()
{
	m_T_probability = m_T_number/(m_G_number + m_C_number + m_A_number + m_T_number);
	return m_T_probability;
}

int Read::empty_reads()
{
	return m_empty_reads;
}

string Read::base_trimmed()
{
	ostringstream os;
	os << m_base_trimmed << "/" << m_base_number << " (" << (double)m_base_trimmed / (double)m_base_number * 100.0 << "%)" << endl;
	return os.str();
}

int Read::trimmed_reads()
{
	return m_trimmed_reads;
}

void Read::reset()
{
	m_empty_reads = 0;
	m_trimmed_reads = 0;
	close();
	m_G_number = 0.0;
	m_C_number = 0.0;
	m_A_number = 0.0;
	m_T_number = 0.0;
	m_base_number = 0;
	m_base_trimmed = 0;

	m_G_probability = 0.25;
	m_C_probability = 0.25;
	m_A_probability = 0.25;
	m_T_probability = 0.25;
}

void Read::open(char* out)
{
	if(!m_out_open)
	{
		if(m_paired >= 1)
		{
			char *out_tmp = new char[strlen(out) + 5];
			strcpy(out_tmp, out);
			strcat(out_tmp, ".tmp");
			if(m_gziped)
				m_out_gz.open(out_tmp);
			else
				m_out.open(out_tmp);
			delete[] out_tmp;
		}
		else
		{
			if(m_gziped)
				m_out_gz.open(out);
			else
				m_out.open(out);
		}
		if(m_gziped)
		{
			if(!m_out_gz.good())
				throw logic_error("ERROR: while opening file");
			m_out_gz.rdbuf()->pubsetbuf(m_buffer, BUFFER_LENGTH);
		}
		else
		{
			if(!m_out.good())
				throw logic_error("ERROR: while opening file");
			m_out.rdbuf()->pubsetbuf(m_buffer, BUFFER_LENGTH);
		}
		m_out_open = true;
	}
}

void Read::close()
{
	if(m_out_open)
	{
		if(m_gziped)
		{
			m_out_gz.close();
			m_out_gz.clear();
		}
		else
		{
			m_out.close();
			m_out.clear();
		}
	}
	m_out_open = false;
}

void Read::remove_empty_reads_paired(char* out, char* outpair, ez::ezRateProgressBar<int>* p, bool v)
{
	close();
	// we only keep the position present in the two list of removed reads
	queue<int> m_paired_pos;
	while(!m_paired_pos_1.empty() && !m_paired_pos_2.empty())
	{
		while(!m_paired_pos_1.empty() && !m_paired_pos_2.empty() && m_paired_pos_1.front() < m_paired_pos_2.front())
		{
			m_paired_pos.push(m_paired_pos_1.front());
			m_paired_pos_1.pop();
		}
		while(!m_paired_pos_1.empty() && !m_paired_pos_2.empty() && m_paired_pos_2.front() < m_paired_pos_1.front())
		{
			m_paired_pos.push(m_paired_pos_2.front());
			m_paired_pos_2.pop();
		}
		while(!m_paired_pos_1.empty() && !m_paired_pos_2.empty() && m_paired_pos_1.front() == m_paired_pos_2.front())
		{
			m_paired_pos.push(m_paired_pos_1.front());
			m_paired_pos_1.pop();
			m_paired_pos_2.pop();
		}
	}
	writeFinal(m_paired_pos, out, p, v);
	writeFinal(m_paired_pos, outpair, p, v);
}

void Read::writeFinal(queue<int> m_paired_pos, char* out, ez::ezRateProgressBar<int>* p, bool v)
{
	// we process the file
	m_paired_pos_1 = m_paired_pos;
	char *out_tmp = new char[strlen(out) + 4];
	strcpy(out_tmp, out);
	strcat(out_tmp, ".tmp");
	ofstream fout;
	ogzstream fout_gz;
	char out_buffer[BUFFER_LENGTH];
	if(m_gziped)
	{
		fout_gz.open(out);
		if(!fout_gz.good())
			throw logic_error("ERROR: while opening output file");
		fout_gz.rdbuf()->pubsetbuf(out_buffer, BUFFER_LENGTH);
	}
	else
	{
		fout.open(out);
		if(!fout.good())
			throw logic_error("ERROR: while opening output file");
		fout.rdbuf()->pubsetbuf(out_buffer, BUFFER_LENGTH);
	}
	igzstream fin;
	char in_buffer[BUFFER_LENGTH];
	fin.open(out_tmp);
	if(!fin.good())
		throw logic_error("1 ERROR: while opening input tmp file");
	fin.rdbuf()->pubsetbuf(in_buffer, BUFFER_LENGTH);
	int line = 0;
	int read_number = 0;
	int next_read_to_skip = 0;
	if(!m_paired_pos_1.empty())
	{
		next_read_to_skip = m_paired_pos_1.front();
		m_paired_pos_1.pop();
	}
	string name, seq, phred;
	name.reserve(2048);
	seq.reserve(2048);
	phred.reserve(2048);
	
	while(!fin.eof() && fin.good())
	{
		// reading of read
		line = 0;
		while(!fin.eof() && fin.good() && line <= 4)
		{
			line++;
			switch(line)
			{
				case 1:
					getline(fin, name);
				break;
				case 2:
					getline(fin, seq);
				break;
				case 3:
					getline(fin, phred);
				break;
				case 4:
					getline(fin, phred);
					read_number++;
					if (v && (read_number+1)%10000 == 0)
						//p->update(read_number);
				break;
			}
		}
		if((read_number < next_read_to_skip || next_read_to_skip == 0) && line >= 4)
		{
			if(m_gziped)
				fout_gz << name << endl << seq << endl << "+" << endl << phred << endl;
			else
				fout << name << endl << seq << endl << "+" << endl << phred << endl;
		}
		else
		{
			if(m_paired_pos_1.empty())
			{
				next_read_to_skip = 0;
			}
			else
			{
				next_read_to_skip = m_paired_pos_1.front();
				m_paired_pos_1.pop();
			}
		}
	}
	if(m_gziped)
	{
		fout_gz.close();
		fout_gz.clear();
	}
	else
	{
		fout.close();
		fout.clear();
	}
	fin.close();
	fin.clear();
	remove(out_tmp);
	delete[] out_tmp;
}
// Definition of the non static stuff

// estimation of base probability for each read
Read::Read(igzstream &fin, char* out, bool gziped, char* N, int phred_score, int threshold, int max_head_trim, int max_tail_trim, int min_read_size, int read_number, bool remove_empty_reads, int min_QC_phred, double min_QC_length, bool estimate, int paired, int strand_bit)
{
	unique_lock<mutex> lk(m_read);
	constructor(fin, N, phred_score, threshold, max_head_trim, max_tail_trim, min_read_size, read_number, remove_empty_reads, min_QC_phred, min_QC_length, strand_bit);
	m_estimation = true;
	m_sampled = false;
	m_paired = paired;
	m_gziped = gziped;
	open(out);
}

// static base probability
Read::Read(igzstream &fin, char* out, bool gziped, char* N, int phred_score, int threshold, int max_head_trim, int max_tail_trim, int min_read_size, int read_number, bool remove_empty_reads, int min_QC_phred, double min_QC_length, int paired, int strand_bit)
{
	unique_lock<mutex> lk(m_read);
	constructor(fin, N, phred_score, threshold, max_head_trim, max_tail_trim, min_read_size, read_number, remove_empty_reads, min_QC_phred, min_QC_length, strand_bit);
	m_sampled = false;
	m_estimation = false;
	m_paired = paired;
	m_gziped = gziped;
	open(out);
}

// sampling of base probability
Read::Read(igzstream &fin, char* N, int phred_score, int threshold, int max_head_trim, int max_tail_trim, int min_read_size, int read_number, bool remove_empty_reads, int min_QC_phred, double min_QC_length, int strand_bit)
{
	unique_lock<mutex> lk(m_read);
	constructor(fin, N, phred_score, threshold, max_head_trim, max_tail_trim, min_read_size, read_number, remove_empty_reads, min_QC_phred, min_QC_length, strand_bit);
	m_sampled = true;
	m_estimation = false;
	m_paired = 0;
}

void Read::constructor(igzstream &fin,char* N, int phred_score, int threshold, int max_head_trim, int max_tail_trim, int min_read_size, int read_number, bool remove_empty_reads, int min_QC_phred, double min_QC_length, int strand_bit)
{
	m_trimmed = false;
	m_size = 0;
	m_cut_begin = 0;
	m_cut_end = 0;
	m_read_number = read_number;
	m_log_read = log(0.0);
	m_log_polyN = log(0.0);
	int line = 1;
	try{
		while(!fin.eof() && fin.good() && line <= 4)
		{
			switch(line)
			{
				case 1:
					getline(fin, m_name);
				break;
				case 2:
					getline(fin, m_seq);
				break;
				case 3:
					getline(fin, m_phred);
				break;
				case 4:
					getline(fin, m_phred);
				break;
			}
			line++;
		}
	}
	catch(exception const& e)
	{
		cerr << "FILE FORMAT ERROR : " << e.what() << endl;
		exit(-1);
	}
	m_size = m_seq.length();
	if(!m_phred_score_set)
		m_phred_score = phred_score;
	
	m_threshold = threshold;
	m_start = 0;
	m_stop = m_size;
	if(max_head_trim >= 0 && max_head_trim <= m_size)
		m_max_head_trim = max_head_trim;
	else
		m_max_head_trim = m_size;
	if(max_tail_trim >= 0 && max_tail_trim <= m_size)
		m_max_tail_trim = m_size - max_tail_trim;
	else
		m_max_tail_trim = 0;
	m_min_read_size = min_read_size;
	if(m_min_read_size > m_size || m_min_read_size < 0)
		m_min_read_size = 0;
	m_N = toupper(*N);
	m_remove_empty_reads = remove_empty_reads;
	m_strand = strand_bit;
	m_init = true;

	if (min_QC_length > 0.0)
	{
		m_QC_check = true;
		m_min_QC_phred = min_QC_phred;
		m_min_QC_length = min_QC_length;
	}
	else
		m_QC_check = false;

	if(line < 4)
		m_init = false;
}

Read::~Read()
{
}

Read& Read::operator=(Read const& readbis)
{
	if(this != &readbis)
	{
		m_trimmed = readbis.m_trimmed;
		m_size = readbis.m_size;
		m_cut_begin = readbis.m_cut_begin;
		m_cut_end = readbis.m_cut_end;
		m_read_number = readbis.m_read_number;
		m_name = readbis.m_name;
		m_seq = readbis.m_seq;
		m_phred = readbis.m_phred;
		m_log_read = readbis.m_log_read;
		m_log_polyN = readbis.m_log_polyN;
		m_N = readbis.m_N;
		m_phred_score = readbis.m_phred_score;
		m_start = readbis.m_start;
		m_stop = readbis.m_stop;
		m_max_head_trim = readbis.m_max_head_trim;
		m_max_tail_trim = readbis.m_max_tail_trim;
		m_min_read_size = readbis.m_min_read_size;
		m_sampled = readbis.m_sampled;
		m_estimation = readbis.m_estimation;
		m_init = readbis.m_init;
	}
	return *this;
}

void Read::run()
{
	if(m_init)
	{
		if(m_sampled || m_estimation)
			Segmentation work(this, true);
		else
			Segmentation work(this, true);
	}
}

void Read::done()
{
	if(m_init)
	{
		if(m_sampled && !m_estimation)
			Segmentation work(this, m_G_probability, m_C_probability, m_A_probability, m_T_probability, m_trimmed, m_cut_begin, m_cut_end);
		if(!m_sampled && !m_estimation)
			writeRead();
		if(m_estimation)
		{
			Segmentation work(this, m_G_probability, m_C_probability, m_A_probability, m_T_probability, m_trimmed, m_cut_begin, m_cut_end);
			writeRead();
		}
	}
}

char Read::base()
{
	return m_N;
}
int Read::size()
{
	return m_size;
}
int Read::start()
{
	return m_start;
}
int Read::stop()
{
	return m_stop;
}
int Read::max_head_trim()
{
	return m_max_head_trim;
}
int Read::max_tail_trim()
{
	return m_max_tail_trim;
}
int Read::min_read_size()
{
	return m_min_read_size;
}
char Read::seq(int i)
{
	return m_seq.at(i);
}
char Read::phred(int i)
{
	return m_phred.at(i) - m_phred_score;
}
int Read::min_QC_phred()
{
	return m_min_QC_phred;
}
int Read::threshold()
{
	return m_threshold;
}
double Read::min_QC_length()
{
	return m_min_QC_length;
}
double Read::G_probability()
{
	return m_G_probability;
}
double Read::C_probability()
{
	return m_C_probability;
}
double Read::A_probability()
{
	return m_A_probability;
}
double Read::T_probability()
{
	return m_T_probability;
}
int Read::strand()
{
	return m_strand;
}
void Read::set_trim(bool trim, int cut_begin, int cut_end)
{
	m_trimmed = trim;
	m_cut_begin = cut_begin;
	m_cut_end = cut_end;
}
bool Read::QC_check()
{
	return m_QC_check;
}

void Read::init(igzstream &fin)
{
	int line = 1;
	string line_tmp;
	while( !fin.eof() && fin.good() && line <= 4)
	{
			switch(line)
			{
				case 1:
					getline(fin, m_name);
				break;
				case 2:
					getline(fin, m_seq);
				break;
				case 3:
					getline(fin, line_tmp);
				break;
				case 4:
					getline(fin, m_phred);
				break;
			}
			line++;
	}
	m_size = m_seq.length();
	m_init = true;
}

inline void Read::writeRead()
{
	// we write the trimmed read if it's not just bad
	if(m_trimmed && m_cut_end - m_cut_begin > 0)
	{
		if(m_gziped)
			writeRead(m_out_gz);
		else
			writeRead(m_out);
	}
	else
	{
		m_empty_reads++;
		if(m_paired > 0)
		{
			if(m_gziped)
				writeRead(m_out_gz);
			else
				writeRead(m_out);
			if (m_paired == 1)
				m_paired_pos_1.push(m_read_number + 1);
			if (m_paired == 2)
				m_paired_pos_2.push(m_read_number + 1);
		}
		else
		{
			if(!m_remove_empty_reads)
			{
				if(m_gziped)
					writeRead(m_out_gz);
				else
					writeRead(m_out);
			}
		}
	}
	if(m_cut_end+1 < m_size || m_cut_begin > 0)
		m_trimmed_reads++;
	m_base_number += m_size;
	m_base_trimmed += m_size - (m_cut_end - m_cut_begin + 1);
}

inline void Read::writeRead(ostream &stream)
{
	if (m_name.at(0) == '@')
		stream << m_name << endl;
	else
		stream << "@" << m_name << endl;
	for (int i = m_cut_begin; i <= m_cut_end; i++)
		stream << m_seq.at(i);
	stream << endl << "+" << endl;
	for (int i = m_cut_begin; i <= m_cut_end; i++)
		stream << m_phred[i];
	stream << endl;
}

int Read::writeSize()
{
	return 1 + m_name.size() + 1 + m_size + 3 + m_size + 1;
}