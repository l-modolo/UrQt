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

#include "Segmentation.hpp"

Segmentation::Segmentation(Read* const read, bool estimate)
{
	m_init = true;
	m_read = read;
	m_proba = nullptr;
	m_cut_begin = 0;
	m_cut_end = 0;
	m_min_QC_phred = m_read->min_QC_phred();
	m_min_QC_length = m_read->min_QC_length();
	m_N = (char)toupper(m_read->base());
	m_G_probability = m_read->G_probability();
	m_C_probability = m_read->C_probability();
	m_A_probability = m_read->A_probability();
	m_T_probability = m_read->T_probability();
	phred();
	if(estimate)
		polyNtrimEstimate();
	else
		polyNtrim();
}

Segmentation::Segmentation(Read* read, double &pG, double &pC, double &pA, double &pT, bool trimmed, int cut_begin, int cut_end)
{
	m_init = true;
	m_read = read;
	m_proba = nullptr;
	phred();
	baseProbaTotal(pG, pC, pA, pT, trimmed, cut_begin, cut_end);
}

Segmentation::~Segmentation()
{
	m_init = false;
	m_read = nullptr;
	if(m_proba != nullptr)
		delete[] m_proba;
}

void Segmentation::phred()
{
	if(m_init)
	{
		m_proba = new double[m_read->size()];
		for (int i = 0; i < m_read->size(); i++)
			m_proba[i] = 1.0 - pow(10.0,- ( (double)(m_read->phred(i)) )/10.0);
	}
}

void Segmentation::polyNtrim()
{
	try
	{
		if(m_init)
		{
			// we find the maximum likelihood cut point between a read segment and a poly N segment
			double logL = log(0.0);
			double newlogL = 0.0;
			// if i = 0 we have a polyN read, if i = m_size we have a read without polyN
			if(m_read->strand() == 0 || m_read->strand() == 2) // we find the cut point for a poly N tail
			{
				for(int i = m_read->start(); i <= m_read->stop(); i++) // the last possible cut point is outside the read
				{
					newlogL = read(m_read->start(), i) + polyN(i, m_read->stop());
					if(newlogL > logL)
					{
						logL = newlogL;
						m_cut_end = i;
					}
				}
			}
			if(m_read->strand() == 1 || m_read->strand() == 2) // we find the cut point for a poly N head
			{
				for(int i = m_read->start(); i <= m_read->stop(); i++) // the last possible cut point is outside the read
				{
					newlogL = reversePolyN(m_read->start(), i) + reverseRead(i, m_read->stop());
					if(newlogL > logL)
					{
						logL = newlogL;
						m_cut_begin = i;
					}
				}
			}
			if(QC_check()) // if the read pass the quality check
				m_read->set_trim(true, m_cut_begin, m_cut_end);
			else
				m_read->set_trim(false, m_cut_begin, m_cut_end);
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void Segmentation::polyNtrim()" << endl;
		exit(-1);
	}
}

void Segmentation::polyNtrimEstimate()
{
	try
	{
		if(m_init)
		{
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

			if(m_read->strand() == 0 || m_read->strand() == 2) // we find the cut point for a poly N tail
			{
				while( labs(newlogL - oldlogL) > 0.01 && iter < 100 )
				{
					oldlogL = newlogL;
					logL = log(0.0);
					// if i = 0 we have a polyN read, if i = m_size we have a read without polyN
					for(int i = m_read->start(); i <= m_read->stop(); i++) // the last possible cut point is outside the read
					{
						newlogL = readEstimate(m_read->start(), i, pG, pC, pA, pT) + polyN(i, m_read->stop());
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
				m_cut_end = m_read->size();
			
			if(m_read->strand() == 1 || m_read->strand() == 2) // we find the cut point for a poly N head
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
					for(int i = m_read->start(); i <= m_cut_end; i++) // the last possible cut point is outside the read
					{
						newlogL = reversePolyN(m_read->start(), i) + reverseReadEstimate(i, m_read->stop(), pG, pC, pA, pT);
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
				if(m_N == 'G' || m_N == 'C' || m_N == 'A' || m_N == 'T')
				{
					if(m_cut_end < m_read->size()) // we keep the last N if we are not doing a QC
						m_cut_end++;
					if(m_cut_begin > 0) // we keep the first N if we are not doing a QC
						m_cut_begin--;
				}
			}
			if(QC_check()) // if the read pass the quality check
				m_read->set_trim(true, m_cut_begin, m_cut_end);
			else
				m_read->set_trim(false, m_cut_begin, m_cut_end);
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void Segmentation::polyNtrim()" << endl;
		exit(-1);
	}
}

double Segmentation::read(int begin, int end)
{
	try
	{
		double logL = 0.0;
		if(begin < end)
		{
			if(end == m_read->start()) // if this is the first iteration we have to compute the full logL
			{
				logL = (end - begin ) * log(1.0 / (end - begin));
				for (int i= begin; i < end; i++)
				{
					logL += probaBaseDict(m_read->seq(i), m_proba[i], m_G_probability, m_C_probability, m_A_probability, m_T_probability) ;
				}
			}
			else // else we only have to add the new base
			{
				if(end != m_read->start()+1)
					logL = m_log_read - (end - begin -1) * log(1.0 / (end - begin - 1)); // we remove the old uniform probability
				logL += (end - begin) * log(1.0 / (end - begin)); // we add the new uniform probability
				logL += probaBaseDict(m_read->seq(end-1), m_proba[end-1], m_G_probability, m_C_probability, m_A_probability, m_T_probability); // we add the new base probability
			}
		}
		m_log_read = logL; // we record the new logL
		return logL;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Segmentation::read(int begin, int end)" << endl;
		exit(-1);
	}
}

double Segmentation::readEstimate(int begin, int end, double pG, double pC, double pA, double pT)
{
	try
	{
		double logL = 0.0;
		if(begin < end)
		{
			if(end == m_read->start()) // if this is the first iteration we have to compute the full logL
			{
				logL = (end - begin ) * log(1.0 / (end - begin));
				for (int i= begin; i < end; i++)
				{
					logL += probaBaseDict(m_read->seq(i), m_proba[i], pG, pC, pA, pT);
				}
			}
			else // else we only have to add the new base
			{
				if(end != m_read->start()+1)
					logL = m_log_read - (end - begin -1) * log(1.0 / (end - begin -1)); // we remove the old uniform probability
				logL += (end - begin) * log(1.0 / (end - begin)); // we add the new uniform probability
				logL += probaBaseDict(m_read->seq(end-1), m_proba[end-1], pG, pC, pA, pT); // we add the new base probability
			}
		}
		m_log_read = logL; // we record the new logL
		return logL;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Segmentation::readEstimate(int begin, int end, double pG, double pC, double pA, double pT)" << endl;
		exit(-1);
	}
}

double Segmentation::polyN(int begin, int end)
{
	try
	{
		if(begin >= end)
			return 0.0;
		double logL = 0.0;
		if(begin < end)
		{
			if(begin == m_read->start()) // if this is the first iteration we have to compute the full logL
			{
				double unifN = log(1.0 / (end - begin));
				double unif = log(1.0 / (end - begin)) + log(1.0/4.0);
				for (int i= begin; i < end; i++)
				{
					if (m_read->seq(i) == m_N)
						logL += unifN + log(m_proba[i]);
					else
						logL += unif + log(1-m_proba[i]);
				}
			}
			else
			{
				logL = m_log_polyN - (end - begin +1) * log(1.0 / (end - begin +1)); // we remove the old uniform probability
				logL += (end - begin) * log(1.0 / (end - begin)); // we add the new uniform probability
				if (m_read->seq(begin-1) == m_N) // we remove the old first base
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
		cerr << "ERROR : " << e.what() << " in : double Segmentation::polyN(int begin, int end)" << endl;
		exit(-1);
	}
}

double Segmentation::reverseRead(int begin, int end)
{
	try
	{
		if(begin >= end)
			return 0.0;
		double logL = 0.0;
		if(begin < end)
		{
			if(begin == m_read->start()) // if this is the first iteration we have to compute the full logL
			{
				logL = (end - begin ) * log(1.0 / (end - begin));
				for (int i= begin; i < end; i++)
				{
					logL += probaBaseDict(m_read->seq(i), m_proba[i], m_G_probability, m_C_probability, m_A_probability, m_T_probability) ;
				}
			}
			else // else we only have to add the new base
			{
				logL = m_log_read - (end - begin +1) * log(1.0 / (end - begin +1));  // we remove the old uniform probability
				logL += (end - begin) * log(1.0 / (end - begin)); // we add the new uniform probability
				logL -= probaBaseDict(m_read->seq(begin-1), m_proba[begin-1], m_G_probability, m_C_probability, m_A_probability, m_T_probability); // we remove the old first base probability
			}
		}
		m_log_read = logL; // we record the new logL
		return logL;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Segmentation::reverseSegmentation(int begin, int end)" << endl;
		exit(-1);
	}
}

double Segmentation::reverseReadEstimate(int begin, int end, double pG, double pC, double pA, double pT)
{
	try
	{
		if(begin >= end)
			return 0.0;
		double logL = 0.0;
		if(begin < end)
		{
			if(begin == m_read->start()) // if this is the first iteration we have to compute the full logL
			{
				logL = (end - begin ) * log(1.0 / (end - begin));
				for (int i= begin; i < end; i++)
				{
					logL += probaBaseDict(m_read->seq(i), m_proba[i], pG, pC, pA, pT);
				}
			}
			else // else we only have to add the new base
			{
				logL = m_log_read - (end - begin +1) * log(1.0 / (end - begin +1)); // we remove the old uniform probability
				logL += (end - begin) * log(1.0 / (end - begin)); // we add the new uniform probability
				logL -= probaBaseDict(m_read->seq(begin-1), m_proba[begin-1], pG, pC, pA, pT); // we remove the old first base probability
			}
		}
		m_log_read = logL; // we record the new logL
		return logL;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Segmentation::reverseSegmentationEstimate(int begin, int end, double pG, double pC, double pA, double pT)" << endl;
		exit(-1);
	}
}

double Segmentation::reversePolyN(int begin, int end)
{
	try
	{
		if(begin >= end)
			return 0.0;
		double logL = 0.0;
		if(begin < end)
		{
			if(end == m_read->start()) // if this is the first iteration we have to compute the full logL
			{
				double unifN = log(1.0 / (end - begin));
				double unif = log(1.0 / (end - begin)) + log(1.0/4.0);
				for (int i= begin; i < end; i++)
				{
					if (m_read->seq(i) == m_N)
						logL += unifN + log(m_proba[i]);
					else
						logL += unif + log(1-m_proba[i]);
				}
			}
			else
			{
				if(end != m_read->start()+1)
					logL = m_log_polyN - (end - begin -1) * log(1.0 / (end - begin -1)); // we remove the old uniform probability
				logL += (end - begin) * log(1.0 / (end - begin)); // we add the new uniform probability
				if (m_read->seq(end-1) == m_N) // we add the old first base
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
		cerr << "ERROR : " << e.what() << " in : double Segmentation::reversePolyN(int begin, int end)" << endl;
		exit(-1);
	}
}

double Segmentation::baseProba(double &pG, double &pC, double &pA, double &pT)
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
				numberBaseDict(m_read->seq(i), m_proba[i], G_number, C_number, A_number, T_number);
			}
			pG = G_number / (G_number + C_number + A_number + T_number);
			pC = C_number / (G_number + C_number + A_number + T_number);
			pA = A_number / (G_number + C_number + A_number + T_number);
			pT = T_number / (G_number + C_number + A_number + T_number);
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Segmentation::baseProba(double &pG, double &pC, double &pA, double &pT)" << endl;
		exit(-1);
	}
}

void Segmentation::baseProbaTotal(double &pG, double &pC, double &pA, double &pT, bool trimmed, int cut_begin, int cut_end)
{
	try
	{
		if(trimmed)
		{
			if(cut_end > 0)
			{
				double G_number = 0.0;
				double C_number = 0.0;
				double A_number = 0.0;
				double T_number = 0.0;
				for (int i = cut_begin; i < cut_end; i++)
				{
					numberBaseDict(m_read->seq(i), m_proba[i], G_number, C_number, A_number, T_number);
				}
				pG += G_number;
				pG += C_number;
				pG += A_number;
				pG += T_number;
			}
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : double Read::baseProbaTotal(double &pG, double &pC, double &pA, double &pT)" << endl;
		exit(-1);
	}
}

bool Segmentation::QC_check()
{
	if(m_read->QC_check())
	{
		int base_checked = 0;
		for (int i = m_cut_begin; i < m_cut_end; i++)
		{
			if ( ((int)m_read->phred(i)) >= m_min_QC_phred)
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

double Segmentation::probaBaseDict(char base, double proba, double pG, double pC, double pA, double pT)
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

void Segmentation::numberBaseDict(char base, double proba, double &G_number, double &C_number, double &A_number, double &T_number)
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