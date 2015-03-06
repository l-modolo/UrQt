# UrQt
UrQt an efficient software for NGS data quality trimming

# Licence

UrQt is licensed under the General Public License v3 (GPLv3).
The last version of this document is available in the UrQt website:
https://lbbe.univ-lyon1.fr/-UrQt-.html

# Installation Procedure :

You can directly clone the UrQt git repository:
```sh
git clone https://github.com/l-modolo/UrQt
```
or download the last version of UrQt: UrQt.1.0.17.tar.gz
```sh
wget ftp://pbil.univ-lyon1.fr/pub/logiciel/UrQt/UrQt.1.0.17.tar.gz
mkdir UrQt
tar xvzf UrQt.1.0.17.tar.gz -C UrQt
```

then compile it
```sh
cd UrQt
make
```
You can compile a static binary with the following commands:
```sh
make clean
make static
```

You should have a UrQt binary in your folder.
Precompiled binary are also available [here](ftp://pbil.univ-lyon1.fr/pub/logiciel/UrQt/bin/).

You may need to install zlib for UrQt to work/compile.
For Ubuntu :
```sh
sudo apt-get install zlib1g-dev
```


# Galaxy Installation Procedure :

With `GALAYX_PATH` the path to your [galaxy distribution](https://wiki.galaxyproject.org/Admin/GetGalaxy) and after the compilation of UrQt.

You can install UrQt as a galaxy tools with the following commands :
```sh
mkdir GALAXY_PATH/tools/UrQt
cp UrQt UrQt.xml GALAXY_PATH/tools/UrQt/
```

and by appening the following line in a relevent section of the file `GALAXY_PATH/config/tool_conf.xml` :
```xml
<tool file="UrQt/UrQt.xml" />
```

Then restart the galaxy server
You can edit the line 5 of the file UrQt.xml to adjust the number of core to use. Any modification to this file requiere a restart of the galaxy server.

# Documentation

The last version of the documentation is available in the UrQt website:
[https://lbbe.univ-lyon1.fr/-UrQt-.html](https://lbbe.univ-lyon1.fr/-UrQt-.html)

UrQt (Unsupervised read Quality trimming) is a fast C++ software to trim nucleotides of unreliable quality from NGS data in fastq or fastq.gz format (automatically detected).
For the phred score encoding, the default is `33` = Sanger (ASCII 33 to 126), but this can be modified with the option
`--phred` to set for example `64` = Illumina 1.3 or `59` = Solexa/Illumina 1.0.

## Single-end

To use UrQt on a single-end fastq of fastq.gz file simply run the following command:
```sh
UrQt --in file.fastq --out file_trimmed.fastq
```
Both input and output files must be accessible and writeable to UrQt to prevent errors.

## Paired-end

To use UrQt on a paired-end fastq of fastq.gz file simply run the following command:
```sh
UrQt --in file_R1.fastq --inpair file_R2.fastq --out file_R2_trimmed.fastq --outpair file_R2_trimmed.fastq
```
By default UrQt remove empty reads (i.e. reads with zero nucleotides of good quality), and keep the correspondence between the paired-end files.
Note that we recommend to use the option `--gz` and output your file in fastq.gz for significant gains of disk space.

## Quality threshold

The quality threshold parameter `--t` threshold define the minimum phred score above which a phred score is considered as "good quality".
By default UrQt use a phred of 5 but this can be changed with the option `--t` threshold.
The classical definition of the quality threshold is obtained with `--t 3.0103`.
Note that UrQt won’t remove every base with a phred score below --t, but will find the best segmentation between two segments of "bad quality" framing a segment of "good quality".
This parameter is independent to the data and must be chosen according to the goal of the analysis.


Example to set a threshold of 10 :
```sh
UrQt --in file.fastq --out file_trimmed.fastq --t 10
```

## Homopolymer trimming

With the option `--N` letter you can define the poly-nucleotide to trim at the head or tail of the sequences.
For letters not present in the standard IUB/IUPAC dictionary, UrQt will perform QC trimming instead of poly-nucleotide trimming.


Example to trim polyA at the head and tail of the reads :
```sh
UrQt --in file.fastq --out file_trimmed.fastq --N A
```

## Verbose mode

By default UrQt display a minimal number of information.
If you want you can use the option `--v` to display all the options used, progress bars and time left.
```sh
UrQt --in file.fastq --out file_trimmed.fastq --v
```

## Multi-threading

By default UrQt use 3 thread (main plus two sub-threads) for a total CPU usage of 100% of one processing unit.
You can use the option `--m` thread_number to use more than one processing unit.
Each additional thread will use a new processing unit.


To run UrQt on 10 processing units:
```sh
UrQt --in file.fastq --out file_trimmed.fastq --m 10
```

## Head or Tail

By default, UrQt start by finding the best cut-point $k_1 \in [1,l]$ with $l$ the size of the read, between a segment of "good quality" and a segment of "bad quality", and then find the best cut-point $k_2 \in [1,k_1]$ between a segment of "bad quality" quality and a segment of "good quality".
Instead of `--pos` both, one can use the parameter `--pos head` to only trim the head of the reads or `--pos tail` to only trim the tail of the reads.


Example to trim the head and tail of the reads :
```sh
UrQt --in file.fastq --out file_trimmed.fastq UrQt --in file.fastq --out file_trimmed.fastq --pos both
```
Example to trim only the head of the reads :
```sh
UrQt --in file.fastq --out file_trimmed.fastq --pos head
```
Example to trim only the tail of the reads :
```sh
UrQt --in file.fastq --out file_trimmed.fastq --pos tail
```

## Minimum trimmed read size

You can tell UrQt to only report reads with a size superior to $n$ nucleotides with the option `--min_read_size n`
Example to report only reads with a size superior to 15 nucleotides :
```sh
UrQt --in file.fastq --out file_trimmed.fastq --min_read_size 15
```

## Maximum number of nucleotides trimmed

You can constrain UrQt to remove no more than n nucleotides at the tail of the reads with the option `--max_tail_trim n`.
The complementary option for the head of the reads is `--max_head_trim n`.


Example to trim at maximum 10 nucleotides at the head of the reads :
```sh
UrQt --in file.fastq --out file_trimmed.fastq --max_head_trim 10
```
Example to trim at maximum 10 nucleotides at the tail of the reads :
```sh
UrQt --in file.fastq --out file_trimmed.fastq --max_tail_trim 10
```
Example to trim at maximum 10 nucleotides at the head and at the tail of the reads :
```sh
UrQt --in file.fastq --out file_trimmed.fastq --max_head_trim 10 --max_tail_trim 10
```
By default empty reads are removed from the output, you can keep them with the option `--r`.


## Classical filter

You can tell UrQt to only keep reads with a minimum of $x$ percent of their length above y phred with the two following options: `--min_QC_length x` and `--min_QC_phred y`.
This filter will be applied after the trimming procedure of UrQt.


For example to retain only reads with a phred of more than 20 on 80% of their length after trimming:
```sh
UrQt --in file.fastq --out file_trimmed.fastq --min_QC_length 80.0 --min_QC_phred 20
```

## Nucleotides probability computation

By default, UrQt use the EM algorithm to compute the proportion of the 4 different nucleotides in a read and estimate the different cut-point.
You can tell UrQt to use fixed proportion of $1/4$ for each nucleotides with the option `--S`.
To compute the proportion of each nucleotide on a sample of size n reads, you can use the option `--s n`.
These two option speed-up the computation but we recommend to use the default parameters (no parameter) for a better estimate of the cut-points in the reads.