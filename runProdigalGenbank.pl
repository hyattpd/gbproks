#!/usr/bin/perl

use strict;
use File::stat;

$| = 1;

###############################################################################
# Script to run Prodigal on all Genbank genomes (both complete and WGS)
#
# Note that this script can be run in parallel by passing it a task id
# and the number of tasks (i.e. "I am processor #5 of 100").  We don't 
# dictate how you run this script many times; that's up to you.  (We use qsub
# in a PBS queuing environment with the -t directive).  If, for example, the
# script is "taskid=5,numtask=100", this instance of the script will only run
# Prodigal on every 5th genome out of 100.  IDs can run from 0 to 
# "numtask minus 1" or from 1 to numtask; either is supported.
#
# In order for parallel execution to work, you need to launch one process for
# every task ID, i.e. for 0-9 with 10 tasks, 
#   runProdigalGenbank.pl -t 0 -n 10
#   runProdigalGenbank.pl -t 1 -n 10
#   runProdigalGenbank.pl -t 2 -n 10
#   runProdigalGenbank.pl -t 3 -n 10
#   runProdigalGenbank.pl -t 4 -n 10
#   runProdigalGenbank.pl -t 5 -n 10
#   runProdigalGenbank.pl -t 6 -n 10
#   runProdigalGenbank.pl -t 7 -n 10
#   runProdigalGenbank.pl -t 8 -n 10
#   runProdigalGenbank.pl -t 9 -n 10
#
# By default, the script just runs in serial mode on one core.
###############################################################################

# Parse command line arguments
my $version = "v0.6.0";
my $numTask = -1;
my $taskId = -1;
my $doAll = 0;
my $rootDir = $ENV{'NCBI_GENOME_ROOT'};
my $prodigal = "prodigal";
my $numArg = scalar @ARGV;
for(my $i = 0; $i < $numArg; $i++) {
  my $lastArg = ($numArg-1==$i?1:0);
  my $arg = $ARGV[$i];
  my $nextArg = $ARGV[$i+1];
  if($arg =~ /^\-/) { $arg =~ tr/A-Z/a-z/; }
  if($arg eq "-r" || $arg eq "--root") {
    if($lastArg == 1) { die "Error: $arg option requires parameter.\n"; }
    $rootDir = $nextArg;
    $i++;
  }
  elsif($arg eq "-p" || $arg eq "--prodigal") {
    if($lastArg == 1) { die "Error: $arg option requires parameter.\n"; }
    $prodigal = $nextArg;
    $i++;
  }
  elsif($arg eq "-n" || $arg eq "--numtask") {
    if($lastArg == 1) { die "Error: $arg option requires parameter.\n"; }
    $numTask = $nextArg;
    $i++;
  }
  elsif($arg eq "-t" || $arg eq "--taskid") {
    if($lastArg == 1) { die "Error: $arg option requires parameter.\n"; }
    $taskId = $nextArg;
    $i++;
  }
  elsif($arg eq "-d" || $arg eq "--doall") { $doAll = 1; }
  elsif($arg eq "-h" || $arg eq "--help") { help($version); }
  elsif($arg eq "-v" || $arg eq "--version") { version($version); }
  else { die "Error: unrecognized option $arg.  Do $0 --help for commands.\n"; }
}
if(!defined($rootDir)) { die "Error: no root directory specified.\n"; }
if(!(-e $rootDir)) { die "Error: $rootDir does not exist!\n"; }

# Check task ID for validity
if($taskId == -1 && $numTask != -1) {
  die "Error: You must specify both -t and -n (do $0 -h for options).\n";
}
elsif($taskId != -1 && $numTask == -1) {
  die "Error: You must specify both -t and -n (do $0 -h for options).\n";
}
elsif($taskId != -1 && $numTask != -1 && ($taskId =~ /[^0-9]/ ||
      $numTask =~ /[^0-9]/)) { 
  die "Error: Task ID and number of tasks must be numerical values.\n";
}
elsif(($taskId < 0 && $numTask != -1) || $taskId > $numTask) { 
  die "Error: Task ID must be between 0 and $numTask!\n";
}
if($numTask == -1) { $numTask = 1; }

# Print version information
print STDERR "##################################################\n";
print STDERR "NCBI Genbank Run_Prodigal Script $version [Mar 2014]\n";
print STDERR "##################################################\n";

# Set up local directories and remote ftp site
my $gbkDir = "$rootDir/genbank";
my $compFnaDir = "$gbkDir/complete_genome_fasta";
my $compProdDir = "$gbkDir/complete_prodigal";
my $wgsFnaDir = "$gbkDir/wgs_genome_fasta";
my $wgsProdDir = "$gbkDir/wgs_prodigal";
my $genomeList = "$gbkDir/genbank.summary.txt";

# If prodigal directory doesn't exist, then create it.
if(!(-e "$compProdDir") && -e $compFnaDir) { 
  print STDERR "Creating $compProdDir since it doesn't exist...\n";
  mkdir "$compProdDir" or die "...Error creating $compProdDir\n";
}
if(!(-e "$wgsProdDir") && -e $wgsFnaDir) { 
  print STDERR "Creating $wgsProdDir since it doesn't exist...\n";
  mkdir "$wgsProdDir" or die "...Error creating $wgsProdDir\n";
}

# Open the list of genomes and process one genome at a time.
print STDERR "Checking for genomes that need to be processed...\n";
my %sawGenome;
my $jobCtr = 0;
my $numSucc = 0, my $numFail = 0;
open FH, $genomeList or die "couldn't open $genomeList for reading\n";
while(my $line = <FH>) {
  next if($line =~ /^#/);
  my @genomeInfo = split /[\t\r\n]+/, $line;
  my $numField = @genomeInfo;
  my $id = $genomeInfo[0];
  my $prodFaa = "$gbkDir/$genomeInfo[$numField-2]";
  my $inputFna = "$gbkDir/$genomeInfo[$numField-1]";
  $sawGenome{$id} = 1;
  if($numTask == 1 || ($jobCtr%$numTask) == ($taskId%$numTask)) {
    if($doAll == 1 || needsProcessing($inputFna, $prodFaa) == 1) {
      print STDERR "...running Prodigal on genome $id...\n";
      my $retVal = runProdigal($prodigal, $inputFna, $prodFaa);
      if($retVal == -1) { 
        warn "...error processing genome $id\n";
        $numFail++;
      }
      else { $numSucc++; }
    }
  }
  $jobCtr++;
}
close FH;
print STDERR "...$numSucc Prodigal jobs run with $numFail failures.\n";

# Delete prodigal output files if no longer needed (only one node
# performs this task).
if($numTask == 1 || ($taskId%$numTask == 0)) {
  print STDERR "Deleting obsolete Prodigal output files...\n";
  my $numDel = 0;
  if(-e $compFnaDir) { # user has complete genomes
    if(opendir(DH, $compProdDir) != 0) {
      while(my $dline = readdir(DH)) {
        next if($dline !~ /\.[a-z][a-z][a-z]$/);
        my $stub = $dline;
        $stub =~ s/\.[a-z][a-z][a-z]$//g;
        if($sawGenome{$stub} eq "") {
          print STDERR "...deleting file $dline...\n";
          unlink("$compProdDir/$dline");
          $numDel++;
        }
      }
      closedir DH;
    }
    else {
      warn "...error opening $compProdDir for reading,";
      warn " unable to delete files...\n";
    }
  }
  if(-e $wgsFnaDir) { # user has wgs genomes
    if(opendir(DH, $wgsProdDir) != 0) {
      while(my $dline = readdir(DH)) {
        next if($dline !~ /\.[a-z][a-z][a-z]$/);
        my $stub = $dline;
        $stub =~ s/\.[a-z][a-z][a-z]$//g;
        if($sawGenome{$stub} eq "") {
          print STDERR "...deleting file $dline...\n";
          unlink("$wgsProdDir/$dline");
          $numDel++;
        }
      }
      closedir DH;
    }
    else {
      warn "...error opening $wgsProdDir for reading,";
      warn " unable to delete files...\n";
    }
  }
  print STDERR "...deleted $numDel files.\n";
}
exit 0;

# Use File::stat to determine if any files are newer than ones we've
# already processed.  If so, we process them.
sub needsProcessing() {
  my $inputFna = shift @_;
  my $prodStub = shift @_;
  my $prodFh, my $curFh;
  my $prodTime, my $curTime;

  $prodStub =~ s/\.faa$//g;
  if(!(-e "$prodStub.stt") || (-s "$prodStub.stt") == 0) { return 1; }
  if(!(-e "$prodStub.faa") || (-s "$prodStub.faa") == 0) { return 1; }
  if(!(-e "$prodStub.dna") || (-s "$prodStub.dna") == 0) { return 1; }
  if(!(-e "$prodStub.gff") || (-s "$prodStub.gff") == 0) { return 1; }
#  if(!(-e "$prodStub.smm") || (-s "$prodStub.smm") == 0) { return 1; }

  open $prodFh, "$prodStub.faa" or return 1;
  $prodTime = stat($prodFh)->mtime;
  close $prodFh;
  if(open($curFh, $inputFna) == 0) {
    warn "couldn't open $inputFna for reading, skipping...\n";
    return 0;
  }
  $curTime = stat($curFh)->mtime;
  close $curFh;
  if($curTime >= $prodTime) { return 1; }
  return 0; 
}

# Subroutine to run Prodigal on a genome. Partial genes are not allowed
# in complete genomes, but are allowed in WGS genomes.
sub runProdigal() {
  my $prodigal = shift @_;
  my $prodInput = shift @_;
  my $prodAmino = shift @_;

  my $anonFlag = "", my $cmd = "";
  my $prodStub = $prodAmino;
  $prodStub =~ s/\.faa$//g;
  my $prodOutput = "$prodStub.gff";
  my $prodAmino = "$prodStub.faa";
  my $prodMessage = "$prodStub.dna";
  my $prodStarts = "$prodStub.stt";
  my $prodSumm = "$prodStub.smm";

  if(!(-e $prodInput)) { 
    warn "...$prodInput file does not exist!...\n";
    return -1; 
  }
  if((-s $prodInput) < 100000) { $anonFlag = "-m anon"; }

  $cmd = "$prodigal -i $prodInput -o $prodOutput -s $prodStarts -a $prodAmino";
  $cmd .= " -d $prodMessage -w $prodSumm $anonFlag -q";
  if($prodInput =~ /^complete/) { $cmd .= " -c"; }
  return (system $cmd) >> 8;
}

# Print help message and exit
sub help() {
  my $version = shift @_;
  print STDERR "\nNCBI Prokaryotic Genome Downloader $version\n\n";
  print STDERR "runProdigalGenbank.pl [-d] [-h] [-n <numtask>] ";
  print STDERR "[-p <prod path>]\n";
  print STDERR "                      [-r <root dir>] [-t <taskID>] [-v]\n\n";
  print STDERR "-d,--doall   :   Do a full download (default: ";
  print STDERR "only do incremental update.\n";
  print STDERR "-h,--help    :   Print help information and exit.\n";
  print STDERR "-r,--root    :   Specify root directory (data is placed in";
  print STDERR "  <root>/genbank).\n";
  print STDERR "-p,--prodigal:   Specify path to the prodigal binary.";
  print STDERR " (Default 'prodigal').\n";
  print STDERR "-t,--taskid  :   Specify a task ID for this process.\n";
  print STDERR "                 (For parallel jobs only).\n";
  print STDERR "-n,--numtask :   Specify the total number of tasks.\n";
  print STDERR "                 (For parallel jobs only).\n";
  print STDERR "-v,--version :   Print version information and exit.\n\n";
  exit(0);
}

# Print version and exit
sub version() {
  my $version = shift @_;
  print STDERR "NCBI Prokaryotic Genome Downloader $version\n";
  exit(0);
}

