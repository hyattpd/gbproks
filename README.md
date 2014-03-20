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
####Download the FASTA and Genbank files for every NCBI bacterial and archaeal genome

text

