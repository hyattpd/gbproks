#!/usr/bin/perl

use strict;
use File::stat;
use LWP::Simple;
use LWP::UserAgent;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

###############################################################################
# Script to Download Microbial Genomes from NCBI
###############################################################################

# Parse command line arguments
my $version = "v0.5.0";
my $doAll = 0;
my $rootDir = $ENV{'NCBI_GENOME_ROOT'};
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
  elsif($arg eq "-d" || $arg eq "--doall") { $doAll = 1; }
  elsif($arg eq "-h" || $arg eq "--help") { help($version); }
  elsif($arg eq "-v" || $arg eq "--version") { version($version); }
  else { die "Error: unrecognized option $arg.  Do $0 --help for commands.\n"; }
}
if(!defined($rootDir)) { die "Error: no root directory specified.\n"; }
if(!(-e $rootDir)) { die "Error: $rootDir does not exist!\n"; }

# Date hash
my @dateText = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
                "Oct", "Nov", "Dec");
my %monthNum;
for(my $i = 0; $i < @dateText; $i++) {
  $monthNum{$dateText[$i]} = $i+1;
}
# Get current date
my ($curDay, $curMonth, $curYear) = (localtime)[3..5];
$curMonth++;
$curYear+= 1900;

# Print header
print STDERR "############################################################\n";
print STDERR "Prokaryotic Genome NCBI Downloading Script $version [Feb 2014]\n";
print STDERR "############################################################\n";

# Set up local directories and remote ftp site
my $gbkDir = "$rootDir/genbank"; 
my $compFnaDir = "$gbkDir/complete_genome_fasta";
my $compGbkDir = "$gbkDir/complete_genome_gbk";
my $wgsDlDir = "$gbkDir/wgs_download";
my $wgsFnaDir = "$gbkDir/wgs_genome_fasta";
my $wgsGbkDir = "$gbkDir/wgs_genome_gbk";
my $prokFtp = "ftp://ftp.ncbi.nlm.nih.gov/";
$prokFtp .= "genomes/GENOME_REPORTS/prokaryotes.txt";
my $wgsFtp = "ftp://ftp.ncbi.nlm.nih.gov/genbank/wgs";
my $ncbiTxt = "$gbkDir/ncbi.genomes.txt";
my $readFh, my $writeFh;

# Set up directories if they don't already exist.
if(!(-e "$gbkDir")) {
  print STDERR "Creating $gbkDir since it doesn't exist...\n";
  mkdir "$gbkDir" or die "...Error creating $gbkDir\n";
}
if(!(-e "$compFnaDir")) {
  print STDERR "Creating $compFnaDir since it doesn't exist...\n";
  mkdir "$compFnaDir" or die "...Error creating $compFnaDir\n";
}
if(!(-e "$compGbkDir")) {
  print STDERR "Creating $compGbkDir since it doesn't exist...\n";
  mkdir "$compGbkDir" or die "...Error creating $compGbkDir\n";
}
if(!(-e "$wgsDlDir")) {
  print STDERR "Creating $wgsDlDir since it doesn't exist...\n";
  mkdir "$wgsDlDir" or die "...Error creating $wgsDlDir\n";
}
if(!(-e "$wgsFnaDir")) {
  print STDERR "Creating $wgsFnaDir since it doesn't exist...\n";
  mkdir "$wgsFnaDir" or die "...Error creating $wgsFnaDir\n";
}
if(!(-e "$wgsGbkDir")) {
  print STDERR "Creating $wgsGbkDir since it doesn't exist...\n";
  mkdir "$wgsGbkDir" or die "...Error creating $wgsGbkDir\n";
}

# Get prokaryotes.txt genome reports file from the ftp site
print STDERR "Grabbing list of prokaryotic genomes from NCBI FTP site...\n";
my $retVal = getstore($prokFtp, $ncbiTxt);
while($retVal != 200) {
  print STDERR "...Error talking to NCBI, retrying in 30 seconds...\n";
  sleep(30);
  $retVal = getstore($prokFtp, $ncbiTxt);
}

# Grab list of WGS projects from the ftp site
print STDERR "Grabbing contents of WGS directory from NCBI FTP site...\n";
my $wgsData = get($wgsFtp);
while(!defined($wgsData)) {
  print STDERR "...Error talking to NCBI, retrying in 30 seconds...\n";
  sleep(30);
  $wgsData = get($wgsFtp);
}

# Go through the list of WGS projects and store the date and file size of each
# fasta/gbk file we see in a hash.
print STDERR "Compiling WGS file sizes and locations into hash...\n";
my %wgsDates, my %wgsSizes, my %wgsFnaFiles, my %wgsGbkFiles;
my @wgsRemote = split /[\n\r]+/, $wgsData;
foreach(@wgsRemote) {
  my @wgsFileInfo = split /[ \t]+/, $_;
  my $file = $wgsFileInfo[8];
  next if($file =~ /^\./ || $file =~ /\.mstr\./);
  next if($file !~ /fsa\_nt/ && $file !~ /gbff/);
  next if($file !~ /\.[A-Z][A-Z][A-Z][A-Z]\./);
  my ($rmMonth, $rmDay, $rmYear) = @wgsFileInfo[5..7];
  $rmMonth = $monthNum{$rmMonth}; 
  if($rmYear =~ /\:/) { 
    $rmYear = $curYear; 
    if($rmMonth > $curMonth || ($rmMonth == $curMonth && $rmDay > $curDay)) {
      $rmYear--;
    }
  }
  my $wgsId = ($file =~ /\.([A-Z][A-Z][A-Z][A-Z])\./)[0];
  if($file =~ /gbff/) { $wgsGbkFiles{$wgsId} = $file; }
  if($file =~ /fsa\_nt/) { $wgsFnaFiles{$wgsId} = $file; }
  $wgsDates{$file} = "$rmMonth/$rmDay/$rmYear";
  $wgsSizes{$file} = $wgsFileInfo[4];
}

# Go through and download any new genomes in the list
print STDERR "Downloading new or modified genome projects...\n";
my %sawWgs, my %sawComp, my $wgsAcc, my $compAcc, my $retVal;
my $localFnaFile, my $localGbkFile;
my $remoteFnaFile, my $remoteGbkFile;
my $numFail = 0, my $numSucc = 0; my $prokFh;
open $prokFh, $ncbiTxt or die "...couldn't open NCBI genome list...\n";
while(my $genome = <$prokFh>) {
  next if($genome =~ /^#/);
  my @genomeData = split /\t/, $genome;

  # WGS
  if($genomeData[12] ne "-" && $genomeData[12] ne "Unplaced") {
    my $wgsAcc = substr($genomeData[12], 0, 4);
    next if($wgsAcc !~ /[A-Z][A-Z][A-Z][A-Z]/); 
    next if(!defined($wgsFnaFiles{$wgsAcc}));
    next if(!defined($wgsGbkFiles{$wgsAcc}));
    $sawWgs{$wgsAcc} = 1;
    $localFnaFile = "$wgsDlDir/$wgsFnaFiles{$wgsAcc}";
    $localGbkFile = "$wgsDlDir/$wgsGbkFiles{$wgsAcc}";
    # Only skip this genome if both zip files exist and are the same size as 
    # at the ftp site.
    # Date information proved unreliable (too many ftp files have newer date
    # but size remains # unchanged), so we don't rely on date.
    next if(-e $localFnaFile && -s $localFnaFile == 
            $wgsSizes{$wgsFnaFiles{$wgsAcc}} && -e $localGbkFile && 
            -s $localGbkFile == $wgsSizes{$wgsGbkFiles{$wgsAcc}} &&
            $doAll == 0);
    print STDERR "...downloading files for WGS project $wgsAcc...\n";
    $retVal = getstore("$wgsFtp/$wgsFnaFiles{$wgsAcc}", $localFnaFile);
    if($retVal != 200) {
      print STDERR "......download failed on $localFnaFile\n"; $numFail++;
    }
    else { $numSucc++; }
    $retVal = getstore("$wgsFtp/$wgsGbkFiles{$wgsAcc}", $localGbkFile);
    if($retVal != 200) {
      print STDERR "......download failed on $localGbkFile\n"; $numFail++;
    }
    else { $numSucc++; }
    next;
  }

  # Skip this genome unless it has a chromosome sequence.
  # Then check the date and compare it to modify date.  Only download
  # if the modify date is newer than what we have.
  next if($genomeData[8] eq "-" && $genomeData[9] eq "-");
  my $modifyDate = ($genomeData[17] ne "-"?$genomeData[17]:$genomeData[16]);
  if($modifyDate eq "-") { print STDERR "......date error $genome"; }
  next if($modifyDate eq "-" || $genomeData[1] eq "-" || $genomeData[3] eq "-");
  $compAcc = "$genomeData[3].$genomeData[1]";
  if(defined($sawComp{$compAcc})) {
    print "error! duplicate entry $compAcc!\n";
  }
  else { $sawComp{$compAcc} = 1; }
  $localGbkFile = "$compGbkDir/$compAcc.gbk";
  $localFnaFile = "$compFnaDir/$compAcc.fna";
  next if(-e $localFnaFile && upToDate($modifyDate, $localGbkFile,
          $localFnaFile) == 1 && -e $localGbkFile && $doAll == 0);

  # Refseq
  my @accessionList = ();
  if($genomeData[8] ne "-") {
    push @accessionList,  split /\,/, $genomeData[8];
    if($genomeData[10] ne "-") {
      push @accessionList,  split /\,/, $genomeData[10];
    }
  }
  # INSDC
  elsif($genomeData[9] ne "-") {
    push @accessionList,  split /\,/, $genomeData[9];
    if($genomeData[10] ne "-") {
      push @accessionList,  split /\,/, $genomeData[11];
    }
  }
  my $accString = join ',', @accessionList;

  print STDERR "...downloading files for complete genome $compAcc...\n";
  # Fetch the fasta sequences for a given list of accession numbers.
  my $efetch = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?';
  $efetch .= 'db=nuccore&rettype=fasta&retmode=text&id='.$accString;
  my $tmpFile = "$gbkDir/sequence.tmp";
  $retVal = getstore($efetch, $tmpFile);
  if($retVal != 200) {
    print STDERR "...download failed on $accString\n"; 
    $numFail++;
    next;
  }
  # Check the sequence for correct number of contigs and remove the
  # blank lines returned by NCBI.
  my $numContig = 0;
  open $readFh, $tmpFile or die "......error opening $tmpFile...\n";
  open $writeFh, ">$localFnaFile" or die "......error opening $localFnaFile\n";
  while(my $line = <$readFh>) {
    next if($line =~ /^\s+$/);
    if($line =~ /^>/) { $numContig++; }
    print $writeFh $line;
  }
  close $writeFh;
  close $readFh;
  if($numContig != @accessionList) {
    print STDERR "......incorrect number of contigs for $compAcc...\n";
    $numFail++;
    unlink($localFnaFile);
    next;
  }
  # Fetch the genbank sequences for a given list of accession numbers.
  my $efetch = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?';
  $efetch .= 'db=nuccore&rettype=gbwithparts&retmode=text&id='.$accString;
  my $tmpFile = "$gbkDir/sequence.tmp";
  $retVal = getstore($efetch, $tmpFile);
  if($retVal != 200) {
    print STDERR "...download failed on $accString\n"; 
    $numFail++;
    next;
  }
  # Check the sequence for correct number of contigs and remove the
  # blank lines returned by NCBI.
  my $numContig = 0;
  open $readFh, $tmpFile or die "......error opening $tmpFile...\n";
  open $writeFh, ">$localGbkFile" or die "......error opening $localGbkFile\n";
  while(my $line = <$readFh>) {
    next if($line =~ /^\s+$/);
    if($line =~ /^\/\//) { $numContig++; }
    print $writeFh $line;
  }
  close $writeFh;
  close $readFh;
  if($numContig != @accessionList) {
    print STDERR "......incorrect number of contigs for $compAcc...\n";
    $numFail++;
    unlink($localFnaFile);
    unlink($localGbkFile);
    next;
  }
  $numSucc++;
  unlink($tmpFile);
}
close $prokFh;
print STDERR "...downloaded $numSucc files with $numFail failures.\n";

# Delete any WGS zip files which are no longer at the ftp site.
print STDERR "Deleting local WGS zip files that are no longer needed...\n";
my $numDel = 0;
opendir DH, $wgsDlDir or die "...open failed on $wgsDlDir\n";
while(my $dline = readdir(DH)) {
  next if($dline !~ /\.gz$/);
  my $wgsId = ($dline =~ /\.([A-Z][A-Z][A-Z][A-Z])\./)[0];
  next if(defined($sawWgs{$wgsId}));
  print STDERR "...deleting file $dline...\n";
  unlink("$wgsDlDir/$dline");
  $numDel++;
}
closedir DH;
print STDERR "...deleted $numDel obsolete files.\n";

# Uncompress any WGS files that need to be uncompressed
print STDERR "Uncompressing WGS zip files...\n";
my $srcFile, my $destFile, my $numUnzip = 0, my $retVal;
foreach(keys %sawWgs) {
  next if(!defined($wgsFnaFiles{$_}) || !defined($wgsGbkFiles{$_}));
  $srcFile = "$wgsDlDir/$wgsFnaFiles{$_}";
  $destFile = "$wgsFnaDir/$_.fna";
  if(needsUnzipping($srcFile, $destFile) == 1) {
    print STDERR "...uncompressing $destFile...\n";
    $numUnzip++;
    $retVal = gunzip $srcFile => $destFile;
    while($retVal == 0) {
      print STDERR "...error unzipping file $destFile $GunzipError,";
      print STDERR " retrying in 30 seconds...\n";
      sleep(30);
      $retVal = gunzip $srcFile => $destFile;
    }
  }
  $srcFile = "$wgsDlDir/$wgsGbkFiles{$_}";
  $destFile = "$wgsGbkDir/$_.gbk";
  if(needsUnzipping($srcFile, $destFile) == 1) {
    print STDERR "...uncompressing $destFile...\n";
    $retVal = gunzip $srcFile => $destFile;
    $numUnzip++;
    while($retVal == 0) {
      print STDERR "...error unzipping file $destFile $GunzipError,";
      print STDERR " retrying in 30 seconds...\n";
      sleep(30);
      $retVal = gunzip $srcFile => $destFile;
    }
  }
}
print STDERR "...unzipped $numUnzip files.\n";

# Delete any WGS gbk/fna files which are no longer at the ftp site.
print STDERR "Deleting local WGS gbk/fna files that are no longer needed...\n";
$numDel = 0;
opendir DH, $wgsFnaDir or die "...open failed on $wgsFnaDir\n";
while(my $dline = readdir(DH)) {
  next if($dline !~ /\.fna$/);
  my $wgsId = substr($dline, 0, 4);
  next if(defined($sawWgs{$wgsId}));
  print STDERR "...deleting file $dline...\n";
  unlink("$wgsFnaDir/$dline");
  $numDel++;
}
closedir DH;
opendir DH, $wgsGbkDir or die "...open failed on $wgsGbkDir\n";
while(my $dline = readdir(DH)) {
  next if($dline !~ /\.gbk$/);
  my $wgsId = substr($dline, 0, 4);
  next if(defined($sawWgs{$wgsId}));
  print STDERR "...deleting file $dline...\n";
  unlink("$wgsGbkDir/$dline");
  $numDel++;
}
closedir DH;
print STDERR "...deleted $numDel obsolete files.\n";

# Now delete refseq/insdc genomes that we no longer need
print STDERR "Deleting local finished gbk/fna files that are";
print STDERR " no longer needed...\n";
$numDel = 0;
opendir DH, $compFnaDir or die "...open failed on $compFnaDir\n";
while(my $dline = readdir(DH)) {
  next if($dline !~ /\.fna$/);
  my $compId = $dline;
  $compId =~ s/\.fna$//g;
  next if(defined($sawComp{$compId}));
  print STDERR "...deleting file $dline...\n";
  unlink("$compFnaDir/$dline");
  $numDel++;
}
closedir DH;
opendir DH, $compGbkDir or die "...open failed on $compGbkDir\n";
while(my $dline = readdir(DH)) {
  next if($dline !~ /\.gbk$/);
  my $compId = $dline;
  $compId =~ s/\.gbk$//g;
  next if(defined($sawComp{$compId}));
  print STDERR "...deleting file $dline...\n";
  unlink("$compGbkDir/$dline");
  $numDel++;
}
closedir DH;
print STDERR "...deleted $numDel obsolete files.\n";
exit 0;

# Return 1 if the dest unzipped file is less recent than the source zipped file.
sub needsUnzipping() {
  my $srcFile = shift @_;
  my $destFile = shift @_;
  my $srcTime, my $destTime;
  my $srcFH, my $destFH;
  if(!(-e $destFile) || (-s $destFile) == 0) { return 1; }
  open $srcFH, $srcFile or return 0;
  $srcTime = stat($srcFH)->mtime;
  close $srcFH;
  open $destFH, $destFile or return 0;
  $destTime = stat($destFH)->mtime;
  close $destFH;
  if($srcTime >= $destTime) { return 1; }
  return 0;
}

# Check if genomes are up to date. Returns 1 if yes, 0 if no.
sub upToDate() {
  my $date = shift @_;
  my $file1 = shift @_;
  my $file2 = shift @_;
  my $tmpFh;

  my ($destYear, $destMonth, $destDay) = split /\//, $date;

  open $tmpFh, $file1 or return 0;
  my ($srcDay, $srcMonth, $srcYear) = (localtime(stat($tmpFh)->mtime))[3..5]; 
  close $tmpFh;
  $srcMonth++; $srcYear += 1900;
  if($destYear > $srcYear) { return 0; }
  elsif($destYear == $srcYear && $destMonth > $srcMonth) { return 0; }
  elsif($destYear == $srcYear && $destMonth == $srcMonth && 
        $destDay >= $srcDay) { return 0; }

  open $tmpFh, $file2 or return 0;
  my ($srcDay, $srcMonth, $srcYear) = (localtime(stat($tmpFh)->mtime))[3..5]; 
  close $tmpFh;
  $srcMonth++; $srcYear += 1900;
  if($destYear > $srcYear) { return 0; }
  elsif($destYear == $srcYear && $destMonth > $srcMonth) { return 0; }
  elsif($destYear == $srcYear && $destMonth == $srcMonth && 
        $destDay >= $srcDay) { return 0; }

  return 1;
}

# Print help message and exit
sub help() {
  my $version = shift @_;
  print STDERR "\nNCBI Prokaryotic Genome Downloader $version\n\n";
  print STDERR "downloadGenbank.pl [-d] [-h] [-r <root dir>] [-v]\n\n";
  print STDERR "-d,--doall   :   Do a full download (default: ";
  print STDERR "only do incremental update.\n";
  print STDERR "-h,--help    :   Print help information and exit.\n";
  print STDERR "-r,--root    :   Specify root directory (data is placed in";
  print STDERR "  <root>/genbank).\n";
  print STDERR "-v,--version :   Print version information and exit.\n\n";
  exit(0);
}

# Print version and exit
sub version() {
  my $version = shift @_;
  print STDERR "NCBI Prokaryotic Genome Downloader $version\n";
  exit(0);
}
