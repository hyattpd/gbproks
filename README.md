## NCBI Prokaryotic Genome Downloading and Processing Scripts

A set of tools to create your own genome repository that contains
all the bacterial and archaeal genomes (both draft and finished)
from NCBI's Genbank.

####Perl Requirements

These scripts have been tested on Perl 5.8+ and use Digest::MD5,
File::Copy, File::Stat, IO::Uncompress, and LWP::Simple, all of which
should be part of a standard Perl installation.  The scripts expect perl
to be in /usr/bin/perl, but you can bypass this by running them with

```
perl downloadGenbank.pl
perl processGenbank.pl
```

or edit the scripts to reflect the correct path to your Perl binary.

####Setting up the Root Directory

The scripts require that you set an environment variable (**$NCBI_GENOME_ROOT**)
or specify the root directory at the script command line.  The download script
will create a **genbank** directory in whichever root directory you specify.

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
./downloadGenbank.pl --root /home/me/repository
./processGenbank.pl --root /home/me/repository
./runProdigalGenbank.pl --root /home/me/repository
```
####Downloading the genomes

To download the genomes, simply do

```
./downloadGenbank.pl --root /home/me/repository
```

This will create the genome repository in **/home/me/repository/genbank**.
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
./downloadGenbank.pl --root /home/me/repository --complete
./downloadGenbank.pl --root /home/me/repository --wgsonly
```

For more details on the downloading process, see the Wiki.  You can also do

```
./downloadGenbank.pl --help
```

for a complete list of options.




