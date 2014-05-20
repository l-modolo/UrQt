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

# define BUFFER_LENGTH 1024000

#include <zlib.h>
#include <cstdlib>
#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <cmath>
#include <vector>
#include <set>
#include <random>
#include <algorithm>
#include <ctime>
#include "mThread/mThread.hpp"
#include "ezRateProgressBar.hpp"
#include "Read.hpp"

#include "gzstream.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

int main(int argc, char **argv)
{
	// getting the parameters
	int c;
	int option_index = 0;
	char* in = nullptr;
	char* inpair = nullptr;
	char* out = nullptr;
	char* outpair = nullptr;
	char* strand = nullptr;
	char N = 'A';
	int phred_score = 1;
	int min_read_size = 0;
	int min_polyN_size = 0;
	int thread_number = 1;
	int sampling = 0;
	bool estimation = true;
	bool remove_empty_reads = true;
	bool fill_empty_reads = false;
	int paired = 0;
	int strand_bit = 0;
	char fill_empty_reads_with = '-';
	double min_QC_length = -1.0;
	double min_QC_phred = -1.0;
	bool v = false;
	bool gziped = false;
	
	static struct option long_options[] =
	{
		{"in", required_argument, 0, 'a'},
		{"inpair", required_argument, 0, 'o'},
		{"out", required_argument, 0, 'b'},
		{"outpair", required_argument, 0, 'p'},
		{"N", required_argument, 0, 'c'},
		{"phred", required_argument, 0, 'd'},
		{"min_read_size", required_argument, 0, 'e'},
		{"min_polyN_size", required_argument, 0, 'f'},
		{"m", required_argument, 0, 'g'},
		{"v", no_argument, 0, 'h'},
		{"s", required_argument, 0, 'i'},
		{"S", no_argument, 0, 'j'},
		{"r", required_argument, 0, 'n'},
		{"R", required_argument, 0, 'k'},
		{"min_QC_length", required_argument, 0, 'l'},
		{"min_QC_phred", required_argument, 0, 'm'},
		{"gz", no_argument, 0, 'q'},
		{"pos", required_argument, 0, 'r'},
		{nullptr, 0, 0, 0}
	};
	
	while ((c = getopt_long_only(argc, argv, "a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r", long_options, &option_index)) != -1) {
		switch(c)
		{
			case 'a': in = optarg;
			break;
			case 'o': inpair = optarg; paired = 1;
			break;
			case 'b': out = optarg;
			break;
			case 'p': outpair = optarg; paired = 1;
			break;
			case 'c': N = *optarg;
			break;
			case 'd': phred_score = atoi(optarg);
			break;
			case 'e': min_read_size = atoi(optarg);
			break;
			case 'f': min_polyN_size = atoi(optarg);
			break;
			case 'g': thread_number = atoi(optarg);
			break;
			case 'h': v = true;
			break;
			case 'i': sampling = atoi(optarg); estimation = false;
			break;
			case 'j': estimation = false;
			break;
			case 'k': remove_empty_reads = false; fill_empty_reads = *optarg;
			break;
			case 'l': min_QC_length = atof(optarg);
			break;
			case 'm': min_QC_phred = atof(optarg);
			break;
			case 'n': remove_empty_reads = false;
			break;
			case 'q': gziped = true;
			break;
			case 'r': strand = optarg;
			break;
			default : cout <<  "without argument : " << optopt << endl;
		}
	}
	if (min_read_size < 0)
		min_read_size = 0;
	if (min_polyN_size < 0)
		min_polyN_size = 0;
	if ( phred_score != 1 && phred_score != 2 && phred_score != 3)
		phred_score = 1;
	if (thread_number < 0)
		thread_number == 1;
	
	// the only required argument is the blast file
	if(in == nullptr || out == nullptr){
		cout <<  "Argument must be defined." << endl;
		cout <<  "Usage: " << argv[0] <<" --in <input.fastq> --out <output.fastq>" << endl;
		cout <<  "       --in input fastq file" << endl;
		cout <<  "       --out output fastq file" << endl;
		cout <<  "Optional:" << endl;
		cout <<  "       --inpair input fastq file for paired end data, only read full of N in both input files will be removed" << endl;
		cout <<  "       --m <number> number of thread to use" << endl;
		cout <<  "       --N <character>  polyN to remove (default: A)" << endl;
		cout <<  "       --phred <number> [1 = Sanger (ASCII 33 to 126), 2 = Illumina 1.3 (ASCII 64 to 126), 3 = Solexa/Illumina 1.0 (ASCII 59 to 126)] (default: 1)" << endl;
		cout <<  "       --min_read_size <number> (default: 0)" << endl;
		cout <<  "       --min_polyN_size <number> (default: 0)" << endl;
		cout <<  "       --min_polyN_size <number> (default: 0)" << endl;
		cout <<  "       --pos <head|tail|both> (expected position of polyN in the sequences)" << endl;
		cout <<  "       --s <number> number of reads to sample to compute the fixe proportion of the 4 different nucleotides instead of being computed in the partitioning of each reads" << endl;
		cout <<  "       --S if present the proportion of the 4 different nucleotides is set to 1/4 instead of being computed in the partitioning of each reads" << endl;
		cout <<  "       --r no removing of empty reads (100% polyN) (default: the empty reads are removed from the output)" << endl;
		cout <<  "       --R <character> if present fill the empty reads (100% polyN) with this letter (default: the empty reads are removed from the output)" << endl;
		cout <<  "       --min_QC_length <double> if present with --min_QC_phred the minimum percentage of base with min_QC_phred necessary to keep a read" << endl;
		cout <<  "       --min_QC_phred <int> if present with --min_QC_length, the minimum phred score on min_QC_length percent of the base necessary to keep a read" << endl;
		cout <<  "       --v verbose" << endl;
		cout <<  "       --gz gziped output" << endl;
		return 1;
	}
	else
	{
		if(v)
		{
			cout << "input:                   " << in << endl;
			if(paired > 0)
				cout << "input pair:              " << inpair << endl;
			cout << "output:                  " << out << endl;
			if(paired > 0)
				cout << "output pair:             " << outpair << endl;
			cout << "removing:                poly" << N << endl;
			
			switch(phred_score)
			{
				case 1 :
					cout << "phred score:             Sanger (ASCII 33 to 126)" << endl;
				break;
				case 2 :
					cout << "phred score:             Illumina 1.3 (ASCII 64 to 126)" << endl;
				break;
				case 3 :
					cout << "phred score:             Solexa/Illumina 1.0 (ASCII 59 to 126)" << endl;
				break;
				default :
					cout << "phred score:             Sanger (ASCII 33 to 126)" << endl;
					phred_score = 1;
			}
			cout << "min read size:           " << min_read_size << endl;
			cout << "min polyN size:          " << min_polyN_size << endl;
			cout << "thread number:           " << thread_number << endl;
			if(strand == nullptr)
			{
				strand_bit = 0;
				cout << "remove ployN at:         end" << endl;
			}
			else
			{
				if(strcmp(strand,"head") == 0)
				{
					strand_bit = 1;
					cout << "remove ployN at:         head" << endl;
				}
				else
				{
					if(strcmp(strand,"both") == 0)
					{
						strand_bit = 2;
						cout << "remove ployN at:         head and tail" << endl;
					}
					else
					{
						strand_bit = 0;
						cout << "remove ployN at:         tail" << endl;
					}
				}
			}
				if(sampling > 0)
				cout << "read to sample:          " << sampling << endl;
			cout << "probability estimation:  ";
			if(estimation)
				cout << "for each reads" << endl;
			else
				cout << "global" << endl;
			if(remove_empty_reads)
				cout << "removing empty reads:    true" << endl;
			else
				if(fill_empty_reads)
				{
					cout << "filling empty reads with:" << fill_empty_reads_with << endl;
					if(min_QC_length > 0.0 && min_QC_phred > 0)
					{
						cout << "removing reads with :" << endl;
						cout << " -phred score          >=" << min_QC_phred << endl;
						cout << " -on:                    " << min_QC_length << "% of their sequences"<< endl;
					}
				}
				else
					cout << "removing empty reads:    false" << endl;
			if(gziped)
				cout << "compressed output" << endl;
		}
	}
	
	char *in_tmp = nullptr;
	char *out_tmp = nullptr;
	while(paired >=0 && paired <= 2)
	{
		igzstream fin;
		char buffer[BUFFER_LENGTH];
		if(paired <= 1)
		{
			in_tmp = new char[strlen(in) + 1];
			out_tmp = new char[strlen(out) + 1];
			strcpy(in_tmp, in);
			strcpy(out_tmp, out);
		}
		else
		{
			in_tmp = new char[strlen(inpair) + 1];
			out_tmp = new char[strlen(outpair) + 1];
			strcpy(in_tmp, inpair);
			strcpy(out_tmp, outpair);
			
		}
		fin.open(in_tmp);
		fin.rdbuf()->pubsetbuf(buffer, BUFFER_LENGTH);
		if(!fin.good())
		{
			cerr << "open of " << in_tmp << " failed" << endl;
			exit(-1);
		}
		
		cout << "processing:              " << in_tmp << endl;
		// counting number of reads
		int number_of_lines = 0;

		string line_tmp;
		line_tmp.reserve(2048);
		while(getline(fin, line_tmp))
			number_of_lines++;
		line_tmp.clear();
		fin.close();
		fin.clear();

		fin.open(in_tmp);
		fin.rdbuf()->pubsetbuf(buffer, BUFFER_LENGTH);
		if(!fin.good())
		{
			cerr << "open of " << in_tmp << " failed" << endl;
			exit(-1);
		}

		number_of_lines = number_of_lines / 4;
		cout <<  "reads number:            " << number_of_lines << endl;
		
		int t_number;
		if(thread_number < 1)
			thread_number = 1;
		if(thread_number > number_of_lines)
		{
			if(number_of_lines > 1)
				t_number = number_of_lines-1;
			else
				t_number = 1;
		}
		else
		{
			t_number = thread_number;
		}
		
		// estimation of the ATGC probability in the reads
		if(estimation || sampling == number_of_lines)
		{
			//trimming of the reads
			mThread<Read> readTrim(t_number, true);
			
			ez::ezRateProgressBar<int> p(number_of_lines);
			p.units = "reads";
			if(v) p.start();
			for ( int i = 0; i < number_of_lines; i++)
			{
				if (v && (i+1)%1000 == 0)
					p.update(i);
				readTrim.add(new Read(fin, out_tmp, gziped, &N, phred_score, min_read_size, min_polyN_size, i, remove_empty_reads, fill_empty_reads, &fill_empty_reads_with, min_QC_phred, min_QC_length, true, paired, strand_bit));
			}
			if(v) p.update(number_of_lines);
			readTrim.stop();

			cout << "number of empty reads:   " << Read::empty_reads() << endl;
			cout << "number of trimmed reads: " << Read::trimmed_reads() << endl;
			cout << "number of trimmer bases: " << Read::base_trimmed() << endl;
		}
		else 
		{
			if(sampling > 0 && sampling != number_of_lines) //estimation of the ATGC probability in the reads from a sample
			{
				cout << "sampling..." << endl;
				
				set<int> reads_to_sample;
				mThread<Read> readSample(t_number, true);
				while(!Read::sampling_done())
				{
					random_device rd;
					default_random_engine generator( rd() );
					// default_random_engine generator;
					uniform_int_distribution<int> distribution(0, number_of_lines-1);
					int i = 0;
					int number_of_lines_to_sample = sampling;
					// we select a random set of reads
					while(i < number_of_lines_to_sample)
					{
						reads_to_sample.insert(distribution(generator));
						i++;
					}
					
					// we compute the ATCG for those reads
					i = 0;
					set<int>::iterator it = reads_to_sample.begin();
					
					while( i < number_of_lines_to_sample && it != reads_to_sample.end())
					{
						if(*it == i)
						{
							readSample.add(new Read(fin, &N, phred_score, min_read_size, min_polyN_size, i, remove_empty_reads, fill_empty_reads, &fill_empty_reads_with, min_QC_phred, min_QC_length, strand_bit));
							it++;
						}
						i++;
					}
				}
				readSample.stop();
				
				cout << "probability of G:        " << Read::G_content() << endl;
				cout << "probability of C:        " << Read::C_content() << endl;
				cout << "probability of A:        " << Read::A_content() << endl;
				cout << "probability of T:        " << Read::T_content() << endl;			}
			
			//trimming of the reads			
			mThread<Read> readTrim(t_number, true);
			
			ez::ezRateProgressBar<int> p(number_of_lines);
			p.units = "reads";
			if(v) p.start();
			for ( int i = 0; i < number_of_lines; i++)
			{
				if (v && (i+1)%1000 == 0)
					p.update(i);
				readTrim.add(new Read(fin, out_tmp, gziped,  &N, phred_score, min_read_size, min_polyN_size, i, remove_empty_reads, fill_empty_reads, &fill_empty_reads_with, min_QC_phred, min_QC_length, paired, strand_bit));
			}
			if(v) p.update(number_of_lines);
			readTrim.stop();
			
			cout << "number of empty reads:   " << Read::empty_reads() << endl;
			cout << "number of trimmed reads: " << Read::trimmed_reads() << endl;
		}
		if(paired > 0)
			paired++;
		else
			paired--;
		Read::reset();
		fin.close();
		fin.clear();
		delete[] in_tmp;
		delete[] out_tmp;
	}
	if(paired > 0)
	{
		cout << "syncronisation of the paired files..." << endl;
		Read::remove_empty_reads_paired(out, outpair);
	}
	return 0;
}

