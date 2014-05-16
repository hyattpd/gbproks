#!/usr/bin/perl

use strict;
use File::stat;
use LWP::Simple;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

sub needsUnzipping($$);
sub upToDate($$$);
sub help($);
sub version($);

###############################################################################
# Script to Download Microbial Genomes from NCBI
###############################################################################

# Parse command line arguments
my $version = "v0.6.0";
my $doAll = 0;
my $doDelete = 0;
my $completeOnly = 0;
my $wgsOnly = 0;
my $unzipOnly = 0;
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
  elsif($arg eq "-c" || $arg eq "--complete") { $completeOnly = 1; }
  elsif($arg eq "-w" || $arg eq "--wgsonly") { $wgsOnly = 1; }
  elsif($arg eq "-u" || $arg eq "--unzip") { $unzipOnly = 1; }
  elsif($arg eq "-a" || $arg eq "--doall") { $doAll = 1; }
  elsif($arg eq "-d" || $arg eq "--dodel") { $doDelete = 1; }
  elsif($arg eq "-h" || $arg eq "--help") { help($version); }
  elsif($arg eq "-v" || $arg eq "--version") { version($version); }
  else { die "Error: unrecognized option $arg.  Do $0 --help for commands.\n"; }
}
if(!defined($rootDir)) { die "Error: no root directory specified.\n"; }

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
print STDERR "Prokaryotic Genome NCBI Downloading Script $version [Mar 2014]\n";
print STDERR "############################################################\n";

# Set up local directories and remote ftp site
my $compFnaDir = "$rootDir/complete_genome_fasta";
my $compGbkDir = "$rootDir/complete_genome_gbk";
my $wgsDlDir = "$rootDir/wgs_download";
my $wgsFnaDir = "$rootDir/wgs_genome_fasta";
my $wgsGbkDir = "$rootDir/wgs_genome_gbk";
my $prokFtp = "ftp://ftp.ncbi.nlm.nih.gov/";
$prokFtp .= "genomes/GENOME_REPORTS/prokaryotes.txt";
my $wgsFtp = "ftp://ftp.ncbi.nlm.nih.gov/genbank/wgs";
my $ncbiTxt = "$rootDir/ncbi.genomes.txt";
my $readFh, my $writeFh;

# Set up directories if they don't already exist.
if(!(-e "$rootDir")) {
  print STDERR "Creating $rootDir since it doesn't exist...\n";
  mkdir "$rootDir" or die "...Error creating $rootDir\n";
}
if(!(-e "$compFnaDir") && $wgsOnly == 0) {
  print STDERR "Creating $compFnaDir since it doesn't exist...\n";
  mkdir "$compFnaDir" or die "...Error creating $compFnaDir\n";
}
if(!(-e "$compGbkDir") && $wgsOnly == 0) {
  print STDERR "Creating $compGbkDir since it doesn't exist...\n";
  mkdir "$compGbkDir" or die "...Error creating $compGbkDir\n";
}
if(!(-e "$wgsDlDir") && $completeOnly == 0) {
  print STDERR "Creating $wgsDlDir since it doesn't exist...\n";
  mkdir "$wgsDlDir" or die "...Error creating $wgsDlDir\n";
}
if(!(-e "$wgsFnaDir") && $completeOnly == 0) {
  print STDERR "Creating $wgsFnaDir since it doesn't exist...\n";
  mkdir "$wgsFnaDir" or die "...Error creating $wgsFnaDir\n";
}
if(!(-e "$wgsGbkDir") && $completeOnly == 0) {
  print STDERR "Creating $wgsGbkDir since it doesn't exist...\n";
  mkdir "$wgsGbkDir" or die "...Error creating $wgsGbkDir\n";
}

# Get prokaryotes.txt genome reports file from the ftp site
if($unzipOnly == 0) {
  print STDERR "Grabbing list of prokaryotic genomes from NCBI FTP site...\n";
  my $retVal = getstore($prokFtp, $ncbiTxt);
  while($retVal != 200) {
    print STDERR "...Error talking to NCBI, retrying in 30 seconds...\n";
    sleep(30);
    $retVal = getstore($prokFtp, $ncbiTxt);
  }
}

# Grab list of WGS projects from the ftp site
my $wgsData;
if($unzipOnly == 0 && $completeOnly == 0) {
  print STDERR "Grabbing contents of WGS directory from NCBI FTP site...\n";
  $wgsData = get($wgsFtp);
  while(!defined($wgsData)) {
    print STDERR "...Error talking to NCBI, retrying in 30 seconds...\n";
    sleep(30);
    $wgsData = get($wgsFtp);
  }
}

# Go through the list of WGS projects and store the date and file size of each
# fasta/gbk file we see in a hash.
my %wgsDates, my %wgsSizes, my %wgsFnaFiles, my %wgsGbkFiles;
if($unzipOnly == 0 && $completeOnly == 0) {
  print STDERR "Compiling WGS file sizes and locations into hash...\n";
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
}

# Go through and download any new genomes in the list
my %sawAssembly, my %newWgs, my $wgsAcc, my $assemblyId;
my $localFnaFile, my $localGbkFile;
my $remoteFnaFile, my $remoteGbkFile;
my $numFail = 0, my $numSucc = 0; my $prokFh;
if($unzipOnly == 0) {
  print STDERR "Downloading new or modified genome projects...\n";
  open($prokFh, $ncbiTxt) or die "...couldn't open NCBI genome list...\n";
  while(my $genome = <$prokFh>) {
    next if($genome =~ /^#/);
    my @genomeData = split /\t/, $genome;
  
    # Skip if no assembly ID
    $assemblyId = $genomeData[21];
    next if($assemblyId eq "-");
    $assemblyId =~ s/\.[^\.]+$//;
    if(defined($sawAssembly{$assemblyId})) {
      print STDERR "...warning duplicate assembly in";
      print STDERR " NCBI file [$assemblyId]...\n";
      next;
    }
  
    # WGS
    if($completeOnly == 0 && $genomeData[12] ne "-" && 
       $genomeData[12] ne "Unplaced") {
      my $wgsAcc = substr($genomeData[12], 0, 4);
      next if($wgsAcc !~ /[A-Z][A-Z][A-Z][A-Z]/); 
      next if(!defined($wgsFnaFiles{$wgsAcc}));
      next if(!defined($wgsGbkFiles{$wgsAcc}));
      $sawAssembly{$assemblyId} = 1;
      $localFnaFile = "$wgsDlDir/$assemblyId.fna.gz";
      $localGbkFile = "$wgsDlDir/$assemblyId.gbk.gz";
      $remoteFnaFile = "$wgsFtp/$wgsFnaFiles{$wgsAcc}";
      $remoteGbkFile = "$wgsFtp/$wgsGbkFiles{$wgsAcc}";
      # Only skip this genome if both zip files exist and are the same size as 
      # at the ftp site.
      # Date information proved unreliable (too many ftp files have newer date
      # but size remains # unchanged), so we don't rely on date.
      next if(-e $localFnaFile && -s $localFnaFile == 
              $wgsSizes{$wgsFnaFiles{$wgsAcc}} && -e $localGbkFile && 
              -s $localGbkFile == $wgsSizes{$wgsGbkFiles{$wgsAcc}} &&
              $doAll == 0);
      print STDERR "...downloading files for WGS project $wgsAcc ";
      print STDERR "[assembly $assemblyId]...\n";
      $newWgs{$assemblyId} = 1;
      my $retVal = getstore($remoteFnaFile, $localFnaFile);
      if($retVal != 200) {
        print STDERR "......download failed on $localFnaFile\n"; $numFail++;
      }
      else { $numSucc++; }
      $retVal = getstore($remoteGbkFile, $localGbkFile);
      if($retVal != 200) {
        print STDERR "......download failed on $localGbkFile\n"; $numFail++;
      }
      else { $numSucc++; }
      next;
    }
    next if($genomeData[12] =~ /^[A-Z][A-Z][A-Z][A-Z]/);
    next if($wgsOnly == 1);
  
    # Skip this genome unless it has a Refseq/INSDC sequence.
    # Then check the date and compare it to modify date.  Only download
    # if the modify date is newer than what we have.
    next if($genomeData[8] eq "-" && $genomeData[9] eq "-" &&
            $genomeData[10] eq "-" && $genomeData[11] eq "-");
    my $modifyDate = ($genomeData[17] ne "-"?$genomeData[17]:$genomeData[16]);
    if($modifyDate eq "-") { print STDERR "......date error $genome"; }
    next if($modifyDate eq "-" || $genomeData[1] eq "-" || 
            $genomeData[3] eq "-");
    $sawAssembly{$assemblyId} = 1;
    $localGbkFile = "$compGbkDir/$assemblyId.gbk";
    $localFnaFile = "$compFnaDir/$assemblyId.fna";
    next if(-e $localFnaFile && upToDate($modifyDate, $localGbkFile,
            $localFnaFile) == 1 && -e $localGbkFile && $doAll == 0);
  
    my @accessionList = (), my @uniqueAccessionList = ();
    # Refseq
    if($genomeData[8] ne "-" || $genomeData[10] ne "-") {
      if($genomeData[8] ne "-") {
        push @accessionList,  split /\,/, $genomeData[8];
      }
      if($genomeData[10] ne "-") {
        push @accessionList,  split /\,/, $genomeData[10];
      }
    }
    # INSDC
    elsif($genomeData[9] ne "-" || $genomeData[11] ne "-") {
      if($genomeData[9] ne "-") {
        push @accessionList,  split /\,/, $genomeData[9];
      }
      if($genomeData[11] ne "-") {
        push @accessionList,  split /\,/, $genomeData[11];
      }
    }
    # Remove redundant entries
    my %sawElement;
    foreach my $element (@accessionList) {
      if(!defined($sawElement{$element})) {
        push @uniqueAccessionList, $element;
        $sawElement{$element} = 1;
      }
    }
    my $accString = join ',', @uniqueAccessionList;
    my $numAcc = scalar @uniqueAccessionList;
  
    print STDERR "...downloading files for complete genome $assemblyId...\n";
    # Fetch the fasta sequences for a given list of accession numbers.
    my $efetch = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?';
    $efetch .= 'db=nuccore&rettype=fasta&retmode=text&id='.$accString;
    my $tmpFile = "$rootDir/sequence.tmp";
    my $retVal = getstore($efetch, $tmpFile);
    if($retVal != 200) {
      print STDERR "...download failed on $accString\n"; 
      $numFail++;
      next;
    }
    # Check the sequence for correct number of contigs and remove the
    # blank lines returned by NCBI.
    my $numContig = 0;
    if(open($readFh, $tmpFile) == 0) { 
      warn "......error opening $tmpFile...\n";
      $numFail++;
      next;
    }
    if(open($writeFh, ">$localFnaFile") == 0) {
      warn "......error opening $localFnaFile\n";
      $numFail++;
      next;
    }
    while(my $line = <$readFh>) {
      next if($line =~ /^\s+$/);
      if($line =~ /^>/) { $numContig++; }
      print $writeFh $line;
    }
    close $writeFh;
    close $readFh;
    if($numContig != $numAcc) {
      print STDERR "......incorrect number of FASTA sequences for ";
      print STDERR "$assemblyId ($numContig vs. $numAcc)...\n";
      $numFail++;
      unlink($localFnaFile);
      next;
    }
    # Fetch the genbank sequences for a given list of accession numbers.
    $efetch = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?';
    $efetch .= 'db=nuccore&rettype=gbwithparts&retmode=text&id='.$accString;
    $tmpFile = "$rootDir/sequence.tmp";
    $retVal = getstore($efetch, $tmpFile);
    if($retVal != 200) {
      print STDERR "...download failed on $accString\n"; 
      $numFail++;
      next;
    }
    # Check the sequence for correct number of contigs and remove the
    # blank lines returned by NCBI.
    $numContig = 0;
    if(open($readFh, $tmpFile) == 0) {
      warn "......error opening $tmpFile...\n";
      $numFail++;
      unlink($localFnaFile);
      next;
    }
    if(open($writeFh, ">$localGbkFile") == 0) {
      warn "......error opening $localGbkFile\n";
      $numFail++;
      unlink($localFnaFile);
      next;
    }
    while(my $line = <$readFh>) {
      next if($line =~ /^\s+$/);
      if($line =~ /^\/\//) { $numContig++; }
      print $writeFh $line;
    }
    close $writeFh;
    close $readFh;
    if($numContig != $numAcc) {
      print STDERR "......incorrect number of Genbank sequences for ";
      print STDERR "$assemblyId ($numContig vs. $numAcc)...\n";
      $numFail++;
      unlink($localFnaFile);
      unlink($localGbkFile);
      next;
    }
    $numSucc++;
    unlink($tmpFile);
  }
  close $prokFh;
  print STDERR "...downloaded $numSucc files successfully with ";
  print STDERR "$numFail failures.\n";
  
  # Delete any obsolete files which are no longer in NCBI's genome list.
  if($doDelete == 1) {
    my @directories = ( $wgsDlDir, $wgsFnaDir, $wgsGbkDir, $compFnaDir, 
                        $compGbkDir );
    my $numDel = 0;
    print STDERR "Deleting local files that are no longer needed...\n";
    foreach my $dir (@directories) {
      next if($completeOnly == 1 && $dir =~ /\/wgs\_/);
      next if($wgsOnly == 1 && $dir =~ /\/complete\_/);
      if(opendir(DH, $dir) == 0) {
        warn "...open failed on $dir, unable to delete files...\n";
        next;
      }
      while(my $dline = readdir(DH)) {
        next if($dline eq "." || $dline eq "..");
        my $id = ($dline =~ /^([^\.]+)\./)[0];
        next if(defined($id) && defined($sawAssembly{$id}));
        print STDERR "...deleting file $dir/$dline...\n";
        unlink("$dir/$dline");
        $numDel++;
      }
      closedir DH;
    }
    print STDERR "...deleted $numDel obsolete files.\n";
  }
}

# Uncompress any WGS files that need to be uncompressed
if($completeOnly == 0) {
  print STDERR "Checking for WGS zip files that need uncompressing...\n";
  my $srcFile, my $destFile, my $numUnzip = 0, my $retVal;
  my $wgsDh, my $dline, my $assId, my $numFail = 0;
  if(opendir($wgsDh, $wgsDlDir) == 1) {
    while($dline = readdir($wgsDh)) {
      next if($dline !~ /\.gz$/);
      $assId = ($dline =~ /^([^\.]+)\./)[0];
      next if(!defined($assId));
      next if($unzipOnly == 0 && !defined($newWgs{$assId})); 
      $srcFile = "$wgsDlDir/$dline";
      if($dline =~ /\.fna\.gz$/) { $destFile = "$wgsFnaDir/$assId.fna"; }
      elsif($dline =~ /\.gbk\.gz$/) { $destFile = "$wgsGbkDir/$assId.gbk"; }
      else { next; }
      if($doAll == 1 || needsUnzipping($srcFile, $destFile) == 1) {
        print STDERR "...uncompressing $destFile...\n";
        $retVal = gunzip $srcFile => $destFile;
        if($retVal == 0) { 
          perror("...error unzipping file $destFile $GunzipError...");
          $numFail++;
        } 
        else { $numUnzip++; }
      }
    }
    closedir($wgsDh);
  }
  print STDERR "...unzipped $numUnzip files with $numFail failures.\n";
}
exit 0;

# Return 1 if the dest unzipped file is less recent than the source zipped file.
sub needsUnzipping($$) {
  my ($srcFile, $destFile) = @_;
  my $srcTime, my $destTime;
  my $srcFH, my $destFH;
  if(!(-e $destFile) || (-s $destFile) == 0) { return 1; }
  open($srcFH, $srcFile) or return 0;
  $srcTime = stat($srcFH)->mtime;
  close $srcFH;
  open($destFH, $destFile) or return 0;
  $destTime = stat($destFH)->mtime;
  close $destFH;
  if($srcTime >= $destTime) { return 1; }
  return 0;
}

# Check if genomes are up to date. Returns 1 if yes, 0 if no.
sub upToDate($$$) {
  my ($date, $file1, $file2) = @_;
  my $tmpFh;

  my ($destYear, $destMonth, $destDay) = split /\//, $date;

  open($tmpFh, $file1) or return 0;
  my ($srcDay, $srcMonth, $srcYear) = (localtime(stat($tmpFh)->mtime))[3..5]; 
  close $tmpFh;
  $srcMonth++; $srcYear += 1900;
  if($destYear > $srcYear) { return 0; }
  elsif($destYear == $srcYear && $destMonth > $srcMonth) { return 0; }
  elsif($destYear == $srcYear && $destMonth == $srcMonth && 
        $destDay >= $srcDay) { return 0; }

  open($tmpFh, $file2) or return 0;
  ($srcDay, $srcMonth, $srcYear) = (localtime(stat($tmpFh)->mtime))[3..5]; 
  close $tmpFh;
  $srcMonth++; $srcYear += 1900;
  if($destYear > $srcYear) { return 0; }
  elsif($destYear == $srcYear && $destMonth > $srcMonth) { return 0; }
  elsif($destYear == $srcYear && $destMonth == $srcMonth && 
        $destDay >= $srcDay) { return 0; }

  return 1;
}

# Print help message and exit
sub help($) {
  my $version = @_;
  print STDERR "\nNCBI Prokaryotic Genome Downloader $version\n\n";
  print STDERR "$0 [-a] [-c] [-d] [-h] [-r <root dir>] [-u] [-v]";
  print STDERR " [-w]\n\n";
  print STDERR "-a,--doall   :   Do a full download/decompression (default: ";
  print STDERR "only do incremental update.\n";
  print STDERR "-d,--dodel   :   Delete files not in NCBI's list (default";
  print STDERR " false).\n";
  print STDERR "-c,--complete:   Only download complete genomes.\n";
  print STDERR "-h,--help    :   Print help information and exit.\n";
  print STDERR "-r,--root    :   Specify root directory.\n";
  print STDERR "-u,--unzip   :   Only decompress files (don't do downloads).\n";
  print STDERR "-v,--version :   Print version information and exit.\n";
  print STDERR "-w,--wgsonly :   Only download WGS genomes.\n\n";
  exit(0);
}

# Print version and exit
sub version($) {
  my $version = @_;
  print STDERR "NCBI Prokaryotic Genome Downloader $version\n";
  exit(0);
}
