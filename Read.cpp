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

double Read::m_G_probability = 0.25;
double Read::m_C_probability = 0.25;
double Read::m_A_probability = 0.25;
double Read::m_T_probability = 0.25;

mutex Read::m_read;
bool Read::m_out_open = false;
ofstream Read::m_out;
bool Read::m_phred_score_set = false;
int Read::m_phred_score = 33.0;
int Read::m_empty_reads = 0;
int Read::m_trimmed_reads = 0;

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

int Read::trimmed_reads()
{
	return m_trimmed_reads;
}

void Read::close()
{
	if(m_out_open)
		m_out.close();
}

// Definition of the non static stuff

Read::Read(gzFile* fin, char* out, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number, bool remove_empty_reads, char* fill_empty_reads, int min_QC_phred, double min_QC_length, bool estimate)
{
	try
	{
		unique_lock<mutex> lk(m_read);
		constructor(fin, N, phred_score, min_read_size, min_polyN_size, read_number, remove_empty_reads, fill_empty_reads, min_QC_phred, min_QC_length);
		if(!m_out_open)
		{
			m_out.open(out, ofstream::out);
			if(!m_out)
				throw logic_error("ERROR: while opening file");
			m_out_open = true;
		}
		m_estimation = true;
		m_sampled = false;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : Read::Read(gzFile* fin, Buffer* buffer, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number)" << endl;
		exit(-1);
	}
}

Read::Read(gzFile* fin, char* out, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number, bool remove_empty_reads, char* fill_empty_reads, int min_QC_phred, double min_QC_length)
{
	try
	{
		unique_lock<mutex> lk(m_read);
		constructor(fin, N, phred_score, min_read_size, min_polyN_size, read_number, remove_empty_reads, fill_empty_reads, min_QC_phred, min_QC_length);
		if(!m_out_open)
		{
			m_out.open(out, ofstream::out);
			if(!m_out)
				throw logic_error("ERROR: while opening file");
			m_out_open = true;
		}
		m_sampled = false;
		m_estimation = false;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : Read::Read(gzFile* fin, Buffer* buffer, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number)" << endl;
		exit(-1);
	}
}

Read::Read(gzFile* fin, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number, bool remove_empty_reads, char* fill_empty_reads, int min_QC_phred, double min_QC_length)
{
	try
	{
		unique_lock<mutex> lk(m_read);
		constructor(fin, N, phred_score, min_read_size, min_polyN_size, read_number, remove_empty_reads, fill_empty_reads, min_QC_phred, min_QC_length);
		m_sampled = true;
		m_estimation = false;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : Read::Read(gzFile* fin, Buffer* buffer, char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number)" << endl;
		exit(-1);
	}
}

void Read::constructor(gzFile* fin,char* N, int phred_score, int min_read_size, int min_polyN_size, int read_number, bool remove_empty_reads, char* fill_empty_reads, int min_QC_phred, double min_QC_length)
{
		m_trimmed = false;
		m_size = 0;
		m_cut = 0;
		m_read_number = read_number;
		m_proba = NULL;
		m_log_read = log(0.0);
		m_log_polyN = log(0.0);
		int line = 1;
		bool new_line = true;
		char c;
		while( (c = gzgetc(*fin)) != EOF && line <= 4)
		{
			if(c == '\n' || c == '\r')
				line++;
			else
			{
				switch(line)
				{
					case 1:
						m_name.push_back(c);
					break;
					case 2:
						m_seq.push_back(c);
					break;
					case 4:
						m_phred.push_back(c);
				}
			}
		}
		m_size = m_seq.length();
		if(!m_phred_score_set)
		{
			switch(phred_score)
			{
				case 1: m_phred_score = 33;
				break;
				case 2: m_phred_score = 64;
				break;
				case 3: m_phred_score = 59;
				break;
				default: m_phred_score = 33;
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
		m_fill_empty_reads = *fill_empty_reads;
		m_init = true;

		if (min_QC_length > 0.0)
		{
			m_QC_check = true;
			m_min_QC_phred = min_QC_phred;
			m_min_QC_length = min_QC_length;
		}
		else
			m_QC_check = false;
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
			m_cut = readbis.m_cut;
			m_read_number = readbis.m_read_number;
			m_name = readbis.m_name;
			m_seq = readbis.m_seq;
			m_phred = readbis.m_phred;
			m_proba = NULL;
			m_log_read = readbis.m_log_read;
			m_log_polyN = readbis.m_log_polyN;
			m_N = readbis.m_N;
			m_phred_score = readbis.m_phred_score;
			m_start = readbis.m_start;
			m_stop = readbis.m_stop;
			m_sampled = readbis.m_sampled;
			m_estimation = readbis.m_estimation;
			m_init = true;
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
		if(m_sampled || m_estimation)
			polyNtrimEstimate();
		else
			polyNtrim();
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
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void Read::done()" << endl;
		exit(-1);
	}
}

void Read::init(gzFile* fin)
{
	try
	{
		int line = 1;
		bool new_line = true;
		char c;
		while( (c = gzgetc(*fin)) != EOF && line <= 4)
		{
			if(c == '\n' || c == '\r')
			{
				line++;
			}
			else
			{
				switch(line)
				{
					case 1:
						m_name.push_back(c);
					break;
					case 2:
						m_seq.push_back(c);
					break;
					case 4:
						m_phred.push_back(c);
				}
			}
		}
		m_size = m_seq.length();
		m_init = true;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void Read::init(gzFile* fin)" << endl;
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
			for(int i = m_start; i <= m_stop; i++) // the last possible cut point is outside the read
			{
				newlogL = read(m_start, i) + polyN(i, m_stop);
				if(newlogL > logL)
				{
					logL = newlogL;
					m_cut = i;
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
			double newlogL = 1.0;
			double oldlogL = 0.0;
			while( labs(newlogL - oldlogL) > 0.01 && iter < 100 )
			{
				oldlogL = newlogL;
				// m_log_read = 0.0;
				// m_log_polyN = 0.0;
				logL = log(0.0);
				// if i = 0 we have a polyN read, if i = m_size we have a read without polyN
				for(int i = m_start; i <= m_stop; i++) // the last possible cut point is outside the read
				{
					newlogL = readEstimate(m_start, i, pG, pC, pA, pT) + polyN(i, m_stop);
					if(newlogL > logL)
					{
						logL = newlogL;
						m_cut = i;
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
			// if(end == m_start) // if this is the first iteration we have to compute the full logL
			// {
				logL = (end - begin ) * log(1.0 / (end - begin));
				for (int i= begin; i < end; i++)
				{
					logL += probaBaseDict(m_seq.at(i), m_proba[i], m_G_probability, m_C_probability, m_A_probability, m_T_probability) ;
				}
			// }
			// else // else we only have to add the new base
			// {
			// 	logL = m_log_read - (end - begin -1) * log(1.0 / (end - begin)); // we remove the old uniform probability
			// 	logL += (end - begin) * log(1.0 / (end - begin)) + probaBaseDict(m_seq.at(end-1), m_proba[end-1], m_G_probability, m_C_probability, m_A_probability, m_T_probability); // we add the new base plus the new uniform probability
			// }
		}
		// m_log_read = logL;
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
			// if(end - m_start == 1) // if this is the first iteration we have to compute the full logL
			// {
				logL = (end - begin ) * log(1.0 / (end - begin));
				for (int i= begin; i < end; i++)
				{
					logL += probaBaseDict(m_seq.at(i), m_proba[i], pG, pC, pA, pT);
				}
			// }
			// else // else we only have to add the new base
			// {
			// 	logL = m_log_read - (end - begin -1) * log(1.0 / (end - begin)); // we remove the old uniform probability
			// 	logL += (end - begin) * log(1.0 / (end - begin)) + probaBaseDict(m_seq.at(end-1), m_proba[end-1], pG, pC, pA, pT); // we add the new base plus the new uniform probability
			// }
		}
		// m_log_read = logL;
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
			// if( begin == m_start) // if this is the first iteration we have to compute the full logL
			// {
				double unifN = log(1.0 / (end - begin));
				double unif = log(1.0 / (end - begin)) + log(1.0/4.0);
				for (int i= begin; i < end; i++)
				{
					if ((char)toupper(m_seq.at(i)) == m_N)
						logL += unifN + log(m_proba[i]);
					else
						logL += unif + log(1-m_proba[i]);
				}
			// }
			// else
			// {
			// 	logL = m_log_polyN - (end - begin +1) * log(1.0 / (end - begin)); // we remove the old uniform probability
			// 	logL += (end - begin) * log(1.0 / (end - begin)); // we add the new uniform probability
			// 	if ((char)toupper(m_seq.at(begin-1)) == m_N) // we remove the old first base
			// 		logL -= log(m_proba[begin-1]);
			// 	else
			// 		logL -= log(1.0/4.0) + log(1-m_proba[begin-1]);
			// }
		}
		// m_log_polyN = logL;
		return logL;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Read::polyN(int begin, int end)" << endl;
		exit(-1);
	}
}


double Read::baseProba(double &pG, double &pC, double &pA, double &pT)
{
	try
	{
		if(m_cut > 0)
		{
			double G_number = 0.0;
			double C_number = 0.0;
			double A_number = 0.0;
			double T_number = 0.0;
			for (int i = 0; i < m_cut; i++)
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
			if(m_cut > 0)
			{
				double G_number = 0.0;
				double C_number = 0.0;
				double A_number = 0.0;
				double T_number = 0.0;
				phred();
				for (int i = 0; i < m_cut; i++)
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
			if( m_cut > m_start +1)
			{
				if (m_name.at(0) == '@')
					m_out << m_name << endl;
				else
					m_out << "@" << m_name << endl;
				for (int i= 0; i < m_cut; i++)
					m_out << m_seq.at(i);
				m_out << endl << "+" << endl;
				for (int i= 0; i < m_cut; i++)
					m_out << m_phred[i];
				m_out << endl;
				if(m_cut < m_size)
					m_trimmed_reads++;
			}
			else
			{
				if(!m_remove_empty_reads)
				{
					if (m_name.at(0) == '@')
						m_out << m_name << endl;
					else
						m_out << "@" << m_name << endl;
					for (int i= 0; i < m_cut; i++)
						m_out << m_fill_empty_reads;
					m_out << endl << "+" << endl;
					for (int i= 0; i < m_cut; i++)
						m_out << m_phred[i];
					m_out << endl;
					if(m_cut < m_size)
						m_trimmed_reads++;
				}
				m_empty_reads++;
			}
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void Read::writeRead()" << endl;
		exit(-1);
	}
}

void Read::writeRead(ostream &stream, int* base_trimmed)
{
	try
	{
		if(m_trimmed)
		{
			// we write the trimmed read if it's not just poly N
			if( m_cut > 0)
			{
				stream << "@" << m_name << endl;
				for (int i= 0; i < m_cut; i++)
					stream << m_seq.at(i);
				stream << endl << "+" << endl;
				for (int i= 0; i < m_cut; i++)
					stream << m_phred[i];
				stream << endl;
				if(m_cut < m_size)
					m_trimmed_reads++;
			}
			*base_trimmed += m_size - m_cut;
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void Read::writeRead(ostream &stream, int* base_trimmed)" << endl;
		exit(-1);
	}
}

string Read::writeRead(int* base_trimmed)
{
	try
	{
		string output;
		if(m_trimmed)
		{
			// we write the trimmed read if it's not just poly N
			if( m_cut > 0)
			{
				output.append("@");
				output.append(m_name);
				output.append("\n");
				for (int i= 0; i < m_cut; i++)
					output.push_back(m_seq.at(i));
				output.append("\n+\n");
				for (int i= 0; i < m_cut; i++)
					output.push_back(m_phred[i]);
				output.append("\n");
				if(m_cut < m_size)
					m_trimmed_reads++;
			}
			*base_trimmed += m_size - m_cut;
		}
		return output;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : string Read::writeRead(int* base_trimmed)" << endl;
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
		for (int i = 0; i < m_cut; i++)
		{
			if ( ((int)m_phred[i]-m_phred_score) >= m_min_QC_phred)
				base_checked++;
		}
		if( (((double)base_checked) / ((double)m_cut) * 100.0) >= m_min_QC_length)
			return true;
		else
			return false;
}
	else
		return true;
}