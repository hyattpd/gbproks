###############################################################################
# NCBI Prokaryotic Genome Downloading and Processing Scripts
#
# This is a simple set of scripts to download the entire set of bacterial and
# archaeal genomes from NCBI and do some parsing/processing on them.  The
# script reads the file located at the following location:
#
# ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt
# 
# Using this file, the download script downloads all Refseq, INSDC, and WGS
# genomes for which it can find data.  The processing script then parses the
# Genbank flat files and creates a summary file "genbank.summary.txt" with
# additional metadata.  The script does not currently download data from the
# SRA.
#
# Current scripts included are:
#
# downloadGenbank.pl:     Downloads the Genbank flat files and FASTA sequences
#                         for every bacterial and archaeal genome.
# processGenbank.pl:      Parses all the FASTA and Genbank files and records
#                         additional metadata such as MD5 checksums, GC
#                         content, individual contig lengths, full definition
#                         and taxonomy lines, number of gaps and nonstandard
#                         bases, etc.
# runProdigalGenbank.pl   Runs Prodigal gene predictions on every genome in
#                         the repository.
#
# ------------
# Requirements
# ------------
#
# Tested with Perl 5.8+
# Dependencies include LWP::Simple, IO::Uncompress, File::Stat, File::Copy,
# and Digest::MD5.
# These should all be included in a standard Perl installation.
# You will need Prodigal (https://github.com/hyattpd/prodigal) to run gene
# predictions.
#
# Each script must know where the "root directory" is located.  The Genbank
# repository is placed in <root dir>/genbank.  You can specify this directory
# using an environment variable (NCBI_GENOME_ROOT), or you can specify the
# directory on the command line of each script with the "--root" option.
#
# --------------
# Example usage:
# --------------

# ./downloadGenbank.pl --root /home/me/data
# ./processGenbank.pl --root /home/me/data
# ./runProdigalGenbank.pl --root /home/me/data --prodigal /home/me/prodigal
#
# ------------------------
# genbank.summary.txt file
# ------------------------
#
# The processing script creates a summary file which adds metadata to the
# original file from the ftp site.  Only genomes for which FASTA and Genbank
# flat files could be downloaded are retained in the summary file.  Columns
# in the file are as follows:
# Id:	A (hopefully) unique ID.  For WGS genomes, this ID is the 4-letter
#	project code.  For other genomes, the ID is of the form 
#	<bioproject id>.<taxon id>.  Note that these IDs sometimes change
#	over time, so existing genomes might suddenly change ID.
# MD5:	An MD5 checksum for this sequence.  See the section below on how these
#	are calculated.
# Organism:	The text name for this organism.
# TaxonID:	NCBI taxon ID for this organism.
# BioProject Acc:	NCBI Bioproject Acc (i.e. PRJNA some number).
# BioProject ID:	NCBI Bioproject ID (numerical)
# Group:	Group from NCBI file, i.e. "Proteobacteria"
# Subgroup:	Subgroup, i.e. "delta/epsilon divisions"
# SeqLen:	Sum of all the bases of all the contigs for the genome
# GC%:	GC content.
# NumSeq:	Number of sequences for this genome.
# Refseq Chr:	Comma-separated list of Refseq accessions for chromosomes.
# Refseq Plasmid:	Comma-separated list of Refseq accs for plasmids.
# INSDC Chr:	Comma-separated list of INSDC accessions for chromosomes.
# INSDC Plasmid:	Comma-separated list of INSDC accs for plasmids.
# WGS:	WGS 4-letter project ID with version.
# Full Taxonomy:	Full taxonomy string from Genbank file.
# Nonst Bases:	Number of non-ACTG bases among all the contigs.
# Num Gaps:	Number of runs of 10 or more N's from all the contigs.
# Passes Checks:	Y or N if it passes some simple quality checks or
#      			not.  Ways it can fail include having different
#			source lines for contigs in the multiple FASTA
#			file, or different taxon IDs, indicating multiple
#			genomes might have been put in one file.
# Gbk Genes:	# of genes in the Genbank flat files
# Gbk Proteins:	# of proteins in the Genbank flat files
# Release Date:	date genome was released
# Modify Date:	date genome was last modified
# Status:	genome status (finished, in scaffolds, SRA, etc.)
# Center:	who sequenced the genome
# Biosample Acc:	Biosample accession number
# Assembly Acc:	Assembly accession number
# Reference:	Journal reference
# Prod Protein File:	Relative path to the Prodigal FASTA file of proteins
# Fasta File:	Relative path to the genome FASTA file.
# 
# Note that this file can sometimes become lost/corrupted if you interrupt
# the processing script, or it exits with an error.  For this reason, a copy
# of the file is saved to <root>/genbank/backup/genbank.summary.txt when the
# script first executes.  The processing script is relatively fast, so you
# can always reconstruct this file by running it if worse comes to worse.

# ------------------------
# MD5 Checksum Calculation
# ------------------------
#
# Checksums are calculated by calling Perl's "md5_hex" function on each
# individual contig.  The resulting checksums are then sorted, joined into a
# single comma-delimited string (i.e. join ',', @checksums), and a final
# md5_hex is called on the comma-delimited string.
#
# ----------------
# The --doall Flag
# ----------------
#
# The scripts are designed to incrementally update, and remove any genomes
# which are no longer present in NCBI's summary file.  If you need to redo
# the entire download/analysis, you can call each script with the "--doall"
# flag.  This will force the scripts to redo everything.
#
# Questions or comments should be directed to doug.hyatt@gmail.com
###############################################################################
