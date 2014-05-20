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
	m_out_open = false;
}

void Read::remove_empty_reads_paired(char* out, char* outpair)
{
	try
	{
		close();
		// we only keep the position present in the two list of removed reads
		queue<int> m_paired_pos;
		while(!m_paired_pos_1.empty() && !m_paired_pos_2.empty())
		{
			while(!m_paired_pos_1.empty() && !m_paired_pos_2.empty() && m_paired_pos_1.front() < m_paired_pos_2.front())
				m_paired_pos_1.pop();
			while(!m_paired_pos_1.empty() && !m_paired_pos_2.empty() && m_paired_pos_2.front() < m_paired_pos_1.front())
				m_paired_pos_2.pop();
			while(!m_paired_pos_1.empty() && !m_paired_pos_2.empty() && m_paired_pos_1.front() == m_paired_pos_2.front())
			{
				m_paired_pos.push(m_paired_pos_1.front());
				m_paired_pos_1.pop();
				m_paired_pos_2.pop();
			}
		}
		// we process the first file
		m_paired_pos_1 = m_paired_pos;
		char *out_tmp = new char[strlen(out) + 5];
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

		string line_tmp;
		line_tmp.reserve(2048);
		int read_number = m_paired_pos_1.front()* 4 - 3;
		int line = 0;
		while(fin.good())
		{
			getline(fin, line_tmp);
			line++;
			if(m_paired_pos_1.empty())
			{
				if(m_gziped)
					fout_gz << line_tmp << endl;
				else
					fout << line_tmp << endl;
			}
			else
			{
				if(line < read_number)
				{
					if(m_gziped)
						fout_gz << line_tmp << endl;
					else
						fout << line_tmp << endl;
				}
				else
				{
					for(int i = 0; i < 3; i++)
					{
						getline(fin, line_tmp);
						line++;
					}
					m_paired_pos_1.pop();
					if(!m_paired_pos_1.empty())
						read_number = m_paired_pos_1.front() * 4 - 3;
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
		// remove(out_tmp);
		delete[] out_tmp;

		// we process the second file
		m_paired_pos_1 = m_paired_pos;
		out_tmp = new char[strlen(outpair) + 5];
		strcpy(out_tmp, outpair);
		strcat(out_tmp, ".tmp");

		if(m_gziped)
		{
			fout_gz.open(outpair);
			if(!fout_gz.good())
				throw logic_error("ERROR: while opening output file");
			fout_gz.rdbuf()->pubsetbuf(out_buffer, BUFFER_LENGTH);
		}
		else
		{
			fout.open(outpair);
			if(!fout.good())
				throw logic_error("ERROR: while opening output file");
			fout.rdbuf()->pubsetbuf(out_buffer, BUFFER_LENGTH);
		}
		fin.open(out_tmp);
		if(!fin.good())
			throw logic_error("2 ERROR: while opening input tmp file");
		fin.rdbuf()->pubsetbuf(in_buffer, BUFFER_LENGTH);

		read_number = m_paired_pos_1.front()* 4 - 3;
		line = 0;
		while(fin.good())
		{
			getline(fin, line_tmp);
			line++;
			if(m_paired_pos_1.empty())
			{
				if(m_gziped)
					fout_gz << line_tmp << endl;
				else
					fout << line_tmp << endl;
			}
			else
			{
				if(line < read_number)
				{
					if(m_gziped)
						fout_gz << line_tmp << endl;
					else
						fout << line_tmp << endl;
				}
				else
				{
					for(int i = 0; i < 3; i++)
					{
						getline(fin, line_tmp);
						line++;
					}
					m_paired_pos_1.pop();
					if(!m_paired_pos_1.empty())
						read_number = m_paired_pos_1.front() * 4 - 3;
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
		// remove(out_tmp);
		delete[] out_tmp;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " " << out << " in : void Read::remove_empty_reads_paired(char* out, char* outpair)" << endl;
		exit(-1);
	}
}

// Definition of the non static stuff

// estimation of base probability for each read
Read::Read(igzstream &fin, char* out, bool gziped, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number, bool remove_empty_reads, int min_QC_phred, double min_QC_length, bool estimate, int paired, int strand_bit)
{
	try
	{
		unique_lock<mutex> lk(m_read);
		constructor(fin, N, phred_score, min_read_size, min_polyN_size, read_number, remove_empty_reads, min_QC_phred, min_QC_length, strand_bit);
		m_estimation = true;
		m_sampled = false;
		m_paired = paired;
		m_gziped = gziped;
		open(out);
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " " << out << " in : Read::Read(igzstream &fin, Buffer* buffer, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number)" << endl;
		exit(-1);
	}
}

// static base probability
Read::Read(igzstream &fin, char* out, bool gziped, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number, bool remove_empty_reads, int min_QC_phred, double min_QC_length, int paired, int strand_bit)
{
	try
	{
		unique_lock<mutex> lk(m_read);
		constructor(fin, N, phred_score, min_read_size, min_polyN_size, read_number, remove_empty_reads, min_QC_phred, min_QC_length, strand_bit);
		m_sampled = false;
		m_estimation = false;
		m_paired = paired;
		m_gziped = gziped;
		open(out);
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " " << out << " in : Read::Read(igzstream &fin, Buffer* buffer, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number)" << endl;
		exit(-1);
	}
}

// sampling of base probability
Read::Read(igzstream &fin, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number, bool remove_empty_reads, int min_QC_phred, double min_QC_length, int strand_bit)
{
	try
	{
		unique_lock<mutex> lk(m_read);
		constructor(fin, N, phred_score, min_read_size, min_polyN_size, read_number, remove_empty_reads, min_QC_phred, min_QC_length, strand_bit);
		m_sampled = true;
		m_estimation = false;
		m_paired = 0;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : Read::Read(igzstream &fin, Buffer* buffer, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number)" << endl;
		exit(-1);
	}
}

void Read::constructor(igzstream &fin,char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number, bool remove_empty_reads, int min_QC_phred, double min_QC_length, int strand_bit)
{
		m_trimmed = false;
		m_size = 0;
		m_cut_begin = 0;
		m_cut_end = 0;
		m_read_number = read_number;
		m_proba = nullptr;
		m_log_read = log(0.0);
		m_log_polyN = log(0.0);
		int line = 1;
		bool new_line = true;
		string line_tmp;
		while( fin.good() && line <= 4)
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
		if(!m_phred_score_set)
		{
			switch(phred_score)
			{
				case 1: m_phred_score = 33.0;
				break;
				case 2: m_phred_score = 64.0;
				break;
				case 3: m_phred_score = 59.0;
				break;
				default: m_phred_score = 33.0;
			}
		}
		m_start = min_read_size;
		if(m_start >= m_size)
			m_start = 0;
		m_stop = m_size - min_polyN_size;
		if(m_stop <= 0)
			m_stop = m_size;
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
	try
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
			m_proba = nullptr;
			m_log_read = readbis.m_log_read;
			m_log_polyN = readbis.m_log_polyN;
			m_N = readbis.m_N;
			m_phred_score = readbis.m_phred_score;
			m_start = readbis.m_start;
			m_stop = readbis.m_stop;
			m_sampled = readbis.m_sampled;
			m_estimation = readbis.m_estimation;
			m_init = readbis.m_init;
		}
		return *this;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : Read& Read::operator=(Read const& readbis)" << endl;
		exit(-1);
	}
}

void Read::run()
{
	try
	{
		if(m_init)
		{
			if(m_sampled || m_estimation)
				polyNtrimEstimate();
			else
				polyNtrim();
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void Read::run()" << endl;
		exit(-1);
	}
}

void Read::done()
{
	try
	{
		if(m_init)
		{
			if(m_sampled && !m_estimation)
				baseProbaTotal();
			if(!m_sampled && !m_estimation)
				writeRead();
			if(m_estimation)
			{
				baseProbaTotal();
				writeRead();
			}
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void Read::done()" << endl;
		exit(-1);
	}
}

void Read::init(igzstream &fin)
{
	try
	{
		int line = 1;
		bool new_line = true;
		string line_tmp;
		while( fin.good() && line <= 4)
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
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void Read::init(igzstream &fin)" << endl;
		exit(-1);
	}
}

void Read::polyNtrim()
{
	try
	{
		if(m_init)
		{
			phred();
			
			// we find the maximum likelihood cut point between a read segment and a poly N segment
			double logL = log(0.0);
			double newlogL = 0.0;
			// if i = 0 we have a polyN read, if i = m_size we have a read without polyN
			if(m_strand == 0 || m_strand == 2) // we find the cut point for a poly N tail
			{
				for(int i = m_start; i <= m_stop; i++) // the last possible cut point is outside the read
				{
					newlogL = read(m_start, i) + polyN(i, m_stop);
					if(newlogL > logL)
					{
						logL = newlogL;
						m_cut_end = i;
					}
				}
			}
			if(m_strand == 1 || m_strand == 2) // we find the cut point for a poly N head
			{
				for(int i = m_start; i <= m_stop; i++) // the last possible cut point is outside the read
				{
					newlogL = reversePolyN(m_start, i) + reverseRead(i, m_stop);
					if(newlogL > logL)
					{
						logL = newlogL;
						m_cut_begin = i;
					}
				}
			}
			delete[] m_proba;
			if(QC_check()) // if the read pass the quality check
				m_trimmed = true;
			else
				m_trimmed = false;
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void Read::polyNtrim()" << endl;
		exit(-1);
	}
}

void Read::polyNtrimEstimate()
{
	try
	{
		if(m_init)
		{
			phred();
			
			double pG = 0.25;
			double pC = 0.25;
			double pA = 0.25;
			double pT = 0.25;
			double old_pG, old_pC, old_pA, old_pT;
			int iter = 0;

			// we find the maximum likelihood cut point between a read segment and a polyN segment 
			double logL = log(0.0);
			m_log_read = log(0.0);
			m_log_polyN = log(0.0);
			double newlogL = 1.0;
			double oldlogL = 0.0;

			if(m_strand == 0 || m_strand == 2) // we find the cut point for a poly N tail
			{
				while( labs(newlogL - oldlogL) > 0.01 && iter < 100 )
				{
					oldlogL = newlogL;
					logL = log(0.0);
					// if i = 0 we have a polyN read, if i = m_size we have a read without polyN
					for(int i = m_start; i <= m_stop; i++) // the last possible cut point is outside the read
					{
						newlogL = readEstimate(m_start, i, pG, pC, pA, pT) + polyN(i, m_stop);
						if(newlogL > logL)
						{
							logL = newlogL;
							m_cut_end = i;
						}
					}
					newlogL = logL;
					old_pG = pG;
					old_pC = pC;
					old_pA = pA;
					old_pT = pT;
					baseProba(pG, pC, pA, pT);
					iter++;
				}
			}
			if( m_cut_end == 0)
				m_cut_end = m_size;
			
			if(m_strand == 1 || m_strand == 2) // we find the cut point for a poly N head
			{
				iter = 0;
				logL = log(0.0);
				m_log_read = log(0.0);
				m_log_polyN = log(0.0);
				newlogL = 1.0;
				oldlogL = 0.0;
				while( labs(newlogL - oldlogL) > 0.01 && iter < 100 )
				{
					oldlogL = newlogL;
					logL = log(0.0);
					// if i = 0 we have a polyN read, if i = m_size we have a read without polyN
					for(int i = m_start; i <= m_cut_end; i++) // the last possible cut point is outside the read
					{
						newlogL = reversePolyN(m_start, i) + reverseReadEstimate(i, m_stop, pG, pC, pA, pT);
						if(newlogL > logL)
						{
							logL = newlogL;
							m_cut_begin = i;
						}
					}
					newlogL = logL;
					old_pG = pG;
					old_pC = pC;
					old_pA = pA;
					old_pT = pT;
					baseProba(pG, pC, pA, pT);
					iter++;
				}
			}
			switch(m_N)
			{
				case 'G':
					if(pG >= 0.99) m_cut_end = m_cut_begin;
				break;
				case 'C':
					if(pC >= 0.99) m_cut_end = m_cut_begin;
				break;
				case 'A':
					if(pA >= 0.99) m_cut_end = m_cut_begin;
				break;
				case 'T':
					if(pT >= 0.99) m_cut_end = m_cut_begin;
				break;
			}
			if(m_cut_end != m_cut_begin)
			{
				if(m_cut_end < m_size) // we keep the last A
					m_cut_end++;
				if(m_cut_begin > 0) // we keep the first A
					m_cut_begin--;
			}
			delete[] m_proba;
			if(QC_check()) // if the read pass the quality check
				m_trimmed = true;
			else
				m_trimmed = false;
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void Read::polyNtrim()" << endl;
		exit(-1);
	}
}

double Read::read(int begin, int end)
{
	try
	{
		double logL = 0.0;
		if(begin < end)
		{
			if(end == m_start) // if this is the first iteration we have to compute the full logL
			{
				logL = (end - begin ) * log(1.0 / (end - begin));
				for (int i= begin; i < end; i++)
				{
					logL += probaBaseDict(m_seq.at(i), m_proba[i], m_G_probability, m_C_probability, m_A_probability, m_T_probability) ;
				}
			}
			else // else we only have to add the new base
			{
				if(end != m_start+1)
					logL = m_log_read - (end - begin -1) * log(1.0 / (end - begin - 1)); // we remove the old uniform probability
				logL += (end - begin) * log(1.0 / (end - begin)); // we add the new uniform probability
				logL += probaBaseDict(m_seq.at(end-1), m_proba[end-1], m_G_probability, m_C_probability, m_A_probability, m_T_probability); // we add the new base probability
			}
		}
		m_log_read = logL; // we record the new logL
		return logL;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Read::read(int begin, int end)" << endl;
		exit(-1);
	}
}

double Read::readEstimate(int begin, int end, double pG, double pC, double pA, double pT)
{
	try
	{
		double logL = 0.0;
		if(begin < end)
		{
			if(end == m_start) // if this is the first iteration we have to compute the full logL
			{
				logL = (end - begin ) * log(1.0 / (end - begin));
				for (int i= begin; i < end; i++)
				{
					logL += probaBaseDict(m_seq.at(i), m_proba[i], pG, pC, pA, pT);
				}
			}
			else // else we only have to add the new base
			{
				if(end != m_start+1)
					logL = m_log_read - (end - begin -1) * log(1.0 / (end - begin -1)); // we remove the old uniform probability
				logL += (end - begin) * log(1.0 / (end - begin)); // we add the new uniform probability
				logL += probaBaseDict(m_seq.at(end-1), m_proba[end-1], pG, pC, pA, pT); // we add the new base probability
			}
		}
		m_log_read = logL; // we record the new logL
		return logL;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Read::readEstimate(int begin, int end, double pG, double pC, double pA, double pT)" << endl;
		exit(-1);
	}
}

double Read::polyN(int begin, int end)
{
	try
	{
		if(begin >= end)
			return 0.0;
		double logL = 0.0;
		if(begin < end)
		{
			if(begin == m_start) // if this is the first iteration we have to compute the full logL
			{
				double unifN = log(1.0 / (end - begin));
				double unif = log(1.0 / (end - begin)) + log(1.0/4.0);
				for (int i= begin; i < end; i++)
				{
					if ((char)toupper(m_seq.at(i)) == m_N)
						logL += unifN + log(m_proba[i]);
					else
						logL += unif + log(1-m_proba[i]);
				}
			}
			else
			{
				logL = m_log_polyN - (end - begin +1) * log(1.0 / (end - begin +1)); // we remove the old uniform probability
				logL += (end - begin) * log(1.0 / (end - begin)); // we add the new uniform probability
				if ((char)toupper(m_seq.at(begin-1)) == m_N) // we remove the old first base
					logL -= log(m_proba[begin-1]);
				else
					logL -= log(1.0/4.0) + log(1-m_proba[begin-1]);
			}
		}
		m_log_polyN = logL;
		return logL;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Read::polyN(int begin, int end)" << endl;
		exit(-1);
	}
}

double Read::reverseRead(int begin, int end)
{
	try
	{
		if(begin >= end)
			return 0.0;
		double logL = 0.0;
		if(begin < end)
		{
			if(begin == m_start) // if this is the first iteration we have to compute the full logL
			{
				logL = (end - begin ) * log(1.0 / (end - begin));
				for (int i= begin; i < end; i++)
				{
					logL += probaBaseDict(m_seq.at(i), m_proba[i], m_G_probability, m_C_probability, m_A_probability, m_T_probability) ;
				}
			}
			else // else we only have to add the new base
			{
				logL = m_log_read - (end - begin +1) * log(1.0 / (end - begin +1));  // we remove the old uniform probability
				logL += (end - begin) * log(1.0 / (end - begin)); // we add the new uniform probability
				logL -= probaBaseDict(m_seq.at(begin-1), m_proba[begin-1], m_G_probability, m_C_probability, m_A_probability, m_T_probability); // we remove the old first base probability
			}
		}
		m_log_read = logL; // we record the new logL
		return logL;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Read::reverseRead(int begin, int end)" << endl;
		exit(-1);
	}
}

double Read::reverseReadEstimate(int begin, int end, double pG, double pC, double pA, double pT)
{
	try
	{
		if(begin >= end)
			return 0.0;
		double logL = 0.0;
		if(begin < end)
		{
			if(begin == m_start) // if this is the first iteration we have to compute the full logL
			{
				logL = (end - begin ) * log(1.0 / (end - begin));
				for (int i= begin; i < end; i++)
				{
					logL += probaBaseDict(m_seq.at(i), m_proba[i], pG, pC, pA, pT);
				}
			}
			else // else we only have to add the new base
			{
				logL = m_log_read - (end - begin +1) * log(1.0 / (end - begin +1)); // we remove the old uniform probability
				logL += (end - begin) * log(1.0 / (end - begin)); // we add the new uniform probability
				logL -= probaBaseDict(m_seq.at(begin-1), m_proba[begin-1], pG, pC, pA, pT); // we remove the old first base probability
			}
		}
		m_log_read = logL; // we record the new logL
		return logL;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Read::reverseReadEstimate(int begin, int end, double pG, double pC, double pA, double pT)" << endl;
		exit(-1);
	}
}

double Read::reversePolyN(int begin, int end)
{
	try
	{
		if(begin >= end)
			return 0.0;
		double logL = 0.0;
		if(begin < end)
		{
			if(end == m_start) // if this is the first iteration we have to compute the full logL
			{
				double unifN = log(1.0 / (end - begin));
				double unif = log(1.0 / (end - begin)) + log(1.0/4.0);
				for (int i= begin; i < end; i++)
				{
					if ((char)toupper(m_seq.at(i)) == m_N)
						logL += unifN + log(m_proba[i]);
					else
						logL += unif + log(1-m_proba[i]);
				}
			}
			else
			{
				if(end != m_start+1)
					logL = m_log_polyN - (end - begin -1) * log(1.0 / (end - begin -1)); // we remove the old uniform probability
				logL += (end - begin) * log(1.0 / (end - begin)); // we add the new uniform probability
				if ((char)toupper(m_seq.at(end-1)) == m_N) // we add the old first base
					logL += log(m_proba[end-1]);
				else
					logL += log(1.0/4.0) + log(1-m_proba[end-1]);
			}
		}
		m_log_polyN = logL;
		return logL;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Read::reversePolyN(int begin, int end)" << endl;
		exit(-1);
	}
}


double Read::baseProba(double &pG, double &pC, double &pA, double &pT)
{
	try
	{
		if(m_cut_end > 0)
		{
			double G_number = 0.0;
			double C_number = 0.0;
			double A_number = 0.0;
			double T_number = 0.0;
			for (int i = m_cut_begin; i < m_cut_end; i++)
			{
				numberBaseDict(m_seq.at(i), m_proba[i], G_number, C_number, A_number, T_number);
			}
			pG = G_number / (G_number + C_number + A_number + T_number);
			pC = C_number / (G_number + C_number + A_number + T_number);
			pA = A_number / (G_number + C_number + A_number + T_number);
			pT = T_number / (G_number + C_number + A_number + T_number);
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Read::baseProba(double &pG, double &pC, double &pA, double &pT)" << endl;
		exit(-1);
	}
}

void Read::baseProbaTotal()
{
	try
	{
		if(m_trimmed)
		{
			if(m_cut_end > 0)
			{
				double G_number = 0.0;
				double C_number = 0.0;
				double A_number = 0.0;
				double T_number = 0.0;
				phred();
				for (int i = m_cut_begin; i < m_cut_end; i++)
				{
					numberBaseDict(m_seq.at(i), m_proba[i], G_number, C_number, A_number, T_number);
				}
				m_G_number += G_number;
				m_C_number += C_number;
				m_A_number += A_number;
				m_T_number += T_number;
				delete[] m_proba;
			}
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Read::baseProbaTotal(double &pG, double &pC, double &pA, double &pT)" << endl;
		exit(-1);
	}
}

void Read::writeRead()
{
	try
	{
		if(m_trimmed)
		{
			// we write the trimmed read if it's not just poly N
			if( m_cut_end > m_start +1 && m_cut_end != m_cut_begin)
			{
				if(m_gziped)
					writeRead(m_out_gz);
				else
					writeRead(m_out);
				if(m_cut_end < m_size || m_cut_begin > 0)
					m_trimmed_reads++;
			}
			else
			{
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
					m_trimmed_reads++;
				}
				if(!m_remove_empty_reads)
				{
					if(m_gziped)
						writeRead(m_out_gz);
					else
						writeRead(m_out);
					m_trimmed_reads++;
				}
				m_empty_reads++;
			}
			m_base_number += m_size;
			m_base_trimmed += (m_size - m_cut_end) + m_cut_begin;
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void Read::writeRead()" << endl;
		exit(-1);
	}
}

void Read::writeRead(ostream &stream)
{
	try
	{
		if (m_name.at(0) == '@')
			stream << m_name << endl;
		else
			stream << "@" << m_name << endl;
		for (int i = m_cut_begin; i < m_cut_end; i++)
			stream << m_seq.at(i);
		stream << endl << "+" << endl;
		for (int i = m_cut_begin; i < m_cut_end; i++)
			stream << m_phred[i];
		stream << endl;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void Read::writeRead(ostream &stream, int* base_trimmed)" << endl;
		exit(-1);
	}
}

int Read::writeSize()
{
	return 1 + m_name.size() + 1 + m_size + 3 + m_size + 1;
}

void Read::phred()
{
	try
	{
		m_proba = new double[m_size];
		for (int i = 0; i < m_size; i++)
		{
			m_proba[i] = 1.0 - pow(10.0,- ( (double)(m_phred[i] - m_phred_score) )/10.0);
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void Read::phred(int phred_score)" << endl;
		exit(-1);
	}
}

double Read::probaBaseDict(char base, double proba, double pG, double pC, double pA, double pT)
{
	double logL = 0.0;
	switch((char)toupper(base))
	{
		case 'G' : logL += log(proba) + log(pG);
		break;
		case 'C' : logL += log(proba) + log(pC);
		break;
		case 'A' : logL += log(proba) + log(pA);
		break;
		case 'T' : logL += log(proba) + log(pT);
		break;
		case 'U' : logL += log(proba) + log(pT);
		break;
		case 'R' : logL += log(proba) + log( (pA + pG)/2.0 ); // A or G 	puRine
		break;
		case 'Y' : logL += log(proba) + log( (pC + pT)/2.0 ); // C, T or U 	pYrimidines
		break;
		case 'K' : logL += log(proba) + log( (pG + pT)/2.0 ); // G, T or U 	bases which are Ketones
		break;
		case 'M' : logL += log(proba) + log( (pA + pC)/2.0 ); // A or C 	bases with aMino groups
		break;
		case 'S' : logL += log(proba) + log( (pC + pG)/2.0 ); // C or G 	Strong interaction
		break;
		case 'W' : logL += log(proba) + log( (pA + pT)/2.0 ); // A, T or U 	Weak interaction
		break;
		case 'B' : logL += log(proba) + log( 1.0 - pA ); // not A (i.e. C, G, T or U) 	B comes after A
		break;
		case 'D' : logL += log(proba) + log( 1.0 - pC ); // not C (i.e. A, G, T or U) 	D comes after C
		break;
		case 'H' : logL += log(proba) + log( 1.0 - pG ); // not G (i.e., A, C, T or U) 	H comes after G
		break;
		case 'V' : logL += log(proba) + log( 1.0 - pT ); // neither T nor U (i.e. A, C or G) 	V comes after U
		break;
		case 'N' : logL += log(proba) + log(1.0/4.0); // A C G T U 	aNy
		break;
		case 'X' : logL -= log(proba) + log(1.0/4.0); // masked
		break;
		case '-' : logL -= log(proba) + log(1.0/4.0); // gap of indeterminate length
		break;
		default:
			logL += log(1.0/4.0) + log(proba);
	}
	return logL;
}

void Read::numberBaseDict(char base, double proba, double &G_number, double &C_number, double &A_number, double &T_number)
{
	switch((char)toupper(base))
	{
		case 'G' : 
			G_number += proba;
			C_number += (1-proba)/3.0;
			A_number += (1-proba)/3.0;
			T_number += (1-proba)/3.0;
		break;
		case 'C' : 
			C_number += proba;
			G_number += (1-proba)/3.0;
			A_number += (1-proba)/3.0;
			T_number += (1-proba)/3.0;
		break;
		case 'A' : 
			A_number += proba;
			G_number += (1-proba)/3.0;
			C_number += (1-proba)/3.0;
			T_number += (1-proba)/3.0;
		break;
		case 'T' : 
			T_number += proba;
			G_number += (1-proba)/3.0;
			C_number += (1-proba)/3.0;
			A_number += (1-proba)/3.0;
		break;
		case 'R' :  // A or G 	puRine
			G_number += proba/2.0;
			C_number += (1-proba)/2.0;
			A_number += proba/2.0;
			T_number += (1-proba)/2.0;
		break;
		case 'Y' : // C, T or U 	pYrimidines
			G_number += (1-proba)/2.0;
			C_number += proba/2.0;
			A_number += (1-proba)/2.0;
			T_number += proba/2.0;
		break;
		case 'K' : // G, T or U 	bases which are Ketones
			G_number += proba/2.0;
			C_number += (1-proba)/2.0;
			A_number += (1-proba)/2.0;
			T_number += proba/2.0;
		break;
		case 'M' : // A or C 	bases with aMino groups
			G_number += (1-proba)/2.0;
			C_number += proba/2.0;
			A_number += proba/2.0;
			T_number += (1-proba)/2.0;
		break;
		case 'S' : // C or G 	Strong interaction
			G_number += proba/2.0;
			C_number += proba/2.0;
			A_number += (1-proba)/2.0;
			T_number += (1-proba)/2.0;
		break;
		case 'W' : // A, T or U 	Weak interaction
			G_number += (1-proba)/2.0;
			C_number += (1-proba)/2.0;
			A_number += proba/2.0;
			T_number += proba/2.0;
		break;
		case 'B' : // not A (i.e. C, G, T or U) 	B comes after A
			G_number += proba/3.0;
			C_number += proba/3.0;
			A_number += 1-proba;
			T_number += proba/3.0;
		break;
		case 'D' : // not C (i.e. A, G, T or U) 	D comes after C
			G_number += proba/3.0;
			C_number += 1-proba;
			A_number += proba/3.0;
			T_number += proba/3.0;
		break;
		case 'H' : // not G (i.e., A, C, T or U) 	H comes after G
			G_number += 1-proba;
			C_number += proba/3.0;
			A_number += proba/3.0;
			T_number += proba/3.0;
		break;
		case 'V' : // neither T nor U (i.e. A, C or G) 	V comes after U
			G_number += proba/3.0;
			C_number += proba/3.0;
			A_number += proba/3.0;
			T_number += 1-proba;
		break;
		case 'N' : // A C G T U 	aNy
			G_number += 1/4.0;
			C_number += 1/4.0;
			A_number += 1/4.0;
			T_number += 1/4.0;
		break;
		case 'X' : // masked
			G_number += 0.0;
			C_number += 0.0;
			A_number += 0.0;
			T_number += 0.0;
		break;
		case '-' : // gap of indeterminate length
			G_number += 0.0;
			C_number += 0.0;
			A_number += 0.0;
			T_number += 0.0;
		break;
		default:
			G_number += 1/4.0;
			C_number += 1/4.0;
			A_number += 1/4.0;
			T_number += 1/4.0;
	}
}

bool Read::QC_check()
{
	if(m_QC_check)
	{
		int base_checked = 0;
		for (int i = m_cut_begin; i < m_cut_end; i++)
		{
			if ( ((int)m_phred[i]-m_phred_score) >= m_min_QC_phred)
				base_checked++;
		}
		if( (((double)base_checked) / ((double)m_cut_end - (double)m_cut_begin) * 100.0) >= m_min_QC_length)
			return true;
		else
			return false;
	}
	else
		return true;
}