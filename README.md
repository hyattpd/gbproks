## GBProks: Download, update, and process prokaryotic genomes from NCBI

GBProks is a set of tools to create and update your own prokaryotic genome
repository containing all the bacterial and archaeal genomes (both draft and finished)
from NCBI's GenBank database.  These scripts will automatically download
(in an incremental fashion) the microbial genomes from NCBI, calculate and record additional
metadata such as MD5 checksums and GC content, and assist the user in running
any tool they wish on each new genome in the repository.

####Perl requirements

These scripts have been tested on Perl 5.8+ and use Digest::MD5,
File::Copy, File::Stat, IO::Uncompress, and LWP::Simple, all of which
should be part of a standard Perl installation.  The scripts expect perl
to be in /usr/bin/perl, but you can bypass this by running them with

```
perl downloadProks.pl
perl processProks.pl
```

or edit the scripts to reflect the correct path to your Perl binary.

####Setting up the root directory

The scripts require that you set an environment variable (**$NCBI_GENOME_ROOT**)
or specify the root directory at the script command line.  The download script
will create the repository in whichever root directory you specify.

To set the variable in bash, do:
```
export NCBI_GENOME_ROOT=/home/me/repository
```

In csh, do:
```
setenv NCBI_GENOME_ROOT /home/me/repository
```

Or, when you call each script, you can specify the root directory manually:
```
./downloadProks.pl --root /home/me/repository
./processProks.pl --root /home/me/repository
./runProdigal.pl --root /home/me/repository
```
####Downloading the genomes

To download the genomes, do the following:

```
./downloadProks.pl --root /home/me/repository
```

This will create the genome repository in **/home/me/repository/**.
Each genome is represented by its assembly ID, a unique identifier created
by NCBI for each genome assembly.  (This means that a single strain may have
multiple assemblies associated with it.)  

Genomes from Refseq and INSDC are placed in the **complete** directories,
i.e. **complete_genome_fasta** and **complete_genome_gbk**.  Genomes from WGS
projects (whole genome shotgun assemblies) are placed in the **wgs** directories,
i.e. **wgs_download**, **wgs_genome_fasta**, and **wgs_genome_gbk**.

You can choose to download only genomes of one type (for example, only the
complete genomes).

```
./downloadProks.pl --root /home/me/repository --complete
./downloadProks.pl --root /home/me/repository --wgsonly
```

For more details on the downloading process, see the Wiki.  You can also do

```
./downloadProks.pl --help
```

for a complete list of options.

####Processing the genomes

Once the genomes have been downloaded, the processing script can be used
to record metadata for each genome and update the summary file with information
for each genome.  When proteins are available, the script also extracts the
protein sequences from the Genbank flat files into the **complete_genome_fasta** 
and **wgs_genome_fasta** subdirectories.

To run the processing script, do:

```
./processProks.pl --root /home/me/repository
```

where the root directory provided is the same as the one you used in the 
download script.  The processing script creates (or updates) a tab-delimited file
**genbank.summary.txt** in the **genbank** directory.  This file contains
a great deal of information about each genome project.  For a detailed description
of each column in this file, see the Wiki.

In addition to the summary file, the script creates a tab-delimited **.txt** file for each
genome in the **complete_metadata** or **wgs_metadata** subdirectories, depending
if the genome is Refseq/INSDC or WGS, respectively.  These files contain metadata
for each contig in the assembly.  For a description of the content of these files,
see the Wiki.

For a detailed list of options, do:

```
./processProks.pl --help
```

####Running Prodigal on each genome

To run Prodigal on each genome, do the following:

```
./runProdigal.pl --root /home/me/repository --prodigal /usr/bin/prodigal
```

providing the path to the repository as well as to the Prodigal binary.  If the
**--prodigal** option is omitted, the script calls the binary 'prodigal' and
assumes it exists in the user's path.  The script writes the various Prodigal output
files to the **complete_prodigal** and **wgs_prodigal** subdirectories, including
**.faa** files for proteins, **.dna** files for gene nucleotide sequences, **.stt** files for
detailed information about each potential start site, **.gff** files for the lists
of gene coordinates, and **.smm** files for summary statistics for each genome.

The **runProdigal.pl** script can be run serially, but it is also designed to
be run in parallel.  For more details on this, do:

```
./runProdigal.pl --help
```

or see the Wiki.

####Running Other Tools on each Genome

The **runGenericTool.pl** script provides a template for the user
to run any tool on the genomes in the repository in an incremental fashion, 
provided the tool uses only one or more of the following as its inputs:

1. The nucleotide FASTA genome sequence
2. The FASTA file of the protein sequences extracted from the Genbank flat file.
3. The FASTA file of the protein sequences calculated by Prodigal.

Example usage:

```
./runGenericTool.pl --config /home/me/aragorn.cfg
```

where the file **aragorn.cfg** contains instructions on which tool to run,
which options to use, and which input and output files to read and write, respectively.
For a detailed description of this configuration file format and script options,
see the Wiki.

####Author

Direct questions and comments to the author, Doug Hyatt, at https://github.com/hyattpd.
