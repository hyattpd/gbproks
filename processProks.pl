#!/usr/bin/perl

###############################################################################
# Script to Process NCBI WGS Draft Genomes and Write Metadata
###############################################################################

use strict;
use File::Copy qw(copy);
use File::stat;
use Digest::MD5 qw(md5_hex);

$| = 1;

sub processGenome($$$$);
sub upToDate($@);
sub help($);
sub version($);

# Parse command line arguments
my $version = "v0.6.0";
my $doAll = 0;
my $completeOnly = 0;
my $wgsOnly = 0;
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
  elsif($arg eq "-d" || $arg eq "--doall") { $doAll = 1; }
  elsif($arg eq "-h" || $arg eq "--help") { help($version); }
  elsif($arg eq "-v" || $arg eq "--version") { version($version); }
  else { die "Error: unrecognized option $arg.  Do $0 --help for commands.\n"; }
}
if(!defined($rootDir)) { die "Error: no root directory specified.\n"; }
if(!(-e $rootDir)) { die "Error: $rootDir does not exist!\n"; }

# Print version information
print STDERR "###########################################################\n";
print STDERR "Prokaryotic Genome NCBI Processing Script $version [Mar 2014]\n";
print STDERR "###########################################################\n";

# Set up local directories and remote ftp site
my $compFnaDir = "$rootDir/complete_genome_fasta";
my $compGbkDir = "$rootDir/complete_genome_gbk";
my $compFaaDir = "$rootDir/complete_genome_faa";
my $compProdDir = "$rootDir/complete_prodigal";
my $compMetaDir = "$rootDir/complete_metadata";
my $wgsFnaDir = "$rootDir/wgs_genome_fasta";
my $wgsGbkDir = "$rootDir/wgs_genome_gbk";
my $wgsFaaDir = "$rootDir/wgs_genome_faa";
my $wgsProdDir = "$rootDir/wgs_prodigal";
my $wgsMetaDir = "$rootDir/wgs_metadata";
my $backupDir = "$rootDir/backup";
my $ncbiText = "$rootDir/ncbi.genomes.txt";
our $summary = "$rootDir/genbank.summary.txt";
our $backupSummary = "$backupDir/genbank.summary.txt";

# If backup/faa/metadata directories don't exist, then create them.
if(!(-e "$backupDir")) { 
  print STDERR "Creating $backupDir since it doesn't exist...\n";
  mkdir "$backupDir" or die "...Error creating $backupDir\n";
}
if(!(-e "$compMetaDir") && $wgsOnly == 0) { 
  print STDERR "Creating $compMetaDir since it doesn't exist...\n";
  mkdir "$compMetaDir" or die "...Error creating $compMetaDir\n";
}
if(!(-e "$compFaaDir") && $wgsOnly == 0) { 
  print STDERR "Creating $compFaaDir since it doesn't exist...\n";
  mkdir "$compFaaDir" or die "...Error creating $compFaaDir\n";
}
if(!(-e "$wgsMetaDir") && $completeOnly == 0) { 
  print STDERR "Creating $wgsMetaDir since it doesn't exist...\n";
  mkdir "$wgsMetaDir" or die "...Error creating $wgsMetaDir\n";
}
if(!(-e "$wgsFaaDir") && $completeOnly == 0) { 
  print STDERR "Creating $wgsFaaDir since it doesn't exist...\n";
  mkdir "$wgsFaaDir" or die "...Error creating $wgsFaaDir\n";
}

# If summary file doesn't exist, then create it
if(!(-e "$summary")) {
  open FH, ">$summary" or die "...error creating $summary\n";
  close FH;
}

# Copy the current file list to a backup file in case something
# goes wrong with our script.
copy($summary,$backupSummary) or die "copy failed on file list\n";

# Read in existing summary information so we can reproduce it.
print STDERR "Reading in summary files...\n";
my $summaryFh, my %genomeData;
open $summaryFh, "$summary" or 
  die "...couldn't open $summary file for reading\n";
while(my $line = <$summaryFh>) {
  my @inf = split /[\n\r\t]+/, $line;
  $genomeData{$inf[0]} = $line; 
}
close FH;

# Open summary file for writing
open $summaryFh, ">$summary" or die "...couldn't create $summary file\n";
print $summaryFh "#Id\tMD5\tOrganism\tTaxonID\tBioProject Acc\tBioProject ID\t";
print $summaryFh "Group\tSubgroup\tSeqLen\tGC%\tNumSeq\tRefseq Chr\tRefseq";
print $summaryFh " Plasmid\tINSDC Chr\tINSDC Plasmid\tWGS\tFull Taxonomy\t";
print $summaryFh "Nonst Bases\tNum Gaps\tPasses Checks\tGbk Genes\tGbk ";
print $summaryFh "Proteins\tRelease Date\tModify Date\tStatus\tCenter\t";
print $summaryFh "Biosample Acc\tProd Protein File\tFasta File\n";

# Cycle through NCBI genomes file and process anything with fasta or
# gbk files more recent than the metadata files.
print "Processing genomes...\n";
my $ncbiFh, my $id, my $fnaFile, my $gbkFile, my $metaFile, my $prodFile;
my %sawAssembly, my $numSucc = 0, my $numFail = 0;
my $metaData;
open $ncbiFh, $ncbiText or die "...couldn't open $ncbiText for reading\n";
while(my $line = <$ncbiFh>) {
  next if($line =~ /^#/);
  my @ncbiInfo = split /\t/, $line;

  # Skip if no assembly id for this genome
  $id = $ncbiInfo[21];
  next if($id eq "-");
  $id =~ s/\.[^\.]+$//;
  if(defined($sawAssembly{$id})) {
    print STDERR "...warning duplicate assembly in";
    print STDERR " NCBI file [$id]...\n";
    next;
  }

  if($ncbiInfo[12] =~ /^[A-Z][A-Z][A-Z][A-Z]/) {
    next if($completeOnly == 1);
    $sawAssembly{$id} = 1;
    $fnaFile = "$wgsFnaDir/$id.fna";
    $gbkFile = "$wgsGbkDir/$id.gbk";
    $prodFile = "$wgsProdDir/$id.faa";
    $metaFile = "$wgsMetaDir/$id.txt";
  }
  elsif($ncbiInfo[8] ne "-" || $ncbiInfo[9] ne "-" ||
        $ncbiInfo[10] ne "-" || $ncbiInfo[11] ne "-") {
    next if($wgsOnly == 1);
    $sawAssembly{$id} = 1;
    $fnaFile = "$compFnaDir/$id.fna";
    $gbkFile = "$compGbkDir/$id.gbk";
    $prodFile = "$compProdDir/$id.faa";
    $metaFile = "$compMetaDir/$id.txt";
  }
  else { next; }
  if(!(-e $fnaFile) || !(-e $gbkFile)) {
    print STDERR "...couldn't find fna/gbk files for $id, skipping...\n";
    $numFail++;
    next;
  }
  # Skip if the metadata file exists and is more recent than the gbk/fna files,
  # and we have info on this genome in the existing summary file.
  if($doAll == 0 && defined($genomeData{$id}) && 
     upToDate($metaFile, $fnaFile, $gbkFile) == 1) {
    print $summaryFh $genomeData{$id};
    next;
  }
  # Process the genome
  print STDERR "...processing genome $id...\n";
  $metaData = processGenome($id, $fnaFile, $gbkFile, $metaFile);
  if(!defined($metaData)) { 
    $numFail++;
    warn "...error processing genome $id...\n"; 
    next;
  }

  # Print information to the genbank summary file
  my @localInfo = split /[\t\r\n]+/, $metaData;
  print $summaryFh "$localInfo[0]\t$localInfo[1]\t";
  for(my $i = 0; $i < 6; $i++) { print $summaryFh "$ncbiInfo[$i]\t"; }
  print $summaryFh "$localInfo[5]\t$localInfo[6]\t$localInfo[2]\t";
  for(my $i = 8; $i < 13; $i++) { print $summaryFh "$ncbiInfo[$i]\t"; }
  print $summaryFh "$localInfo[4]\t$localInfo[7]\t";
  print $summaryFh "$localInfo[8]\t$localInfo[9]\t";
  for(my $i = 14; $i < 21; $i++) { print $summaryFh "$ncbiInfo[$i]\t"; }
  $prodFile = substr($prodFile, length($rootDir));
  $fnaFile = substr($fnaFile, length($rootDir));
  print $summaryFh "$prodFile\t$fnaFile\n";
  $numSucc++;
}
close $ncbiFh;
close $summaryFh;
print STDERR "...$numSucc genomes successfully processed with";
print STDERR " $numFail failures.\n";

# Now delete complete genome metadata files that we no longer need
print STDERR "Deleting local metadata files that are no longer needed...\n";
my $numDel = 0;
if($wgsOnly == 0) {
  opendir DH, $compMetaDir or die "...open failed on $compMetaDir\n";
  while(my $dline = readdir(DH)) {
    next if($dline !~ /\.txt$/);
    my $compId = $dline;
    $compId =~ s/\.txt$//g;
    next if(defined($sawAssembly{$compId}));
    print STDERR "...deleting file $dline...\n";
    unlink("$compMetaDir/$dline");
    $numDel++;
  }
  closedir(DH);
}
# Same for WGS metadata files
if($completeOnly == 0) {
  opendir DH, $wgsMetaDir or die "...open failed on $wgsMetaDir\n";
  while(my $dline = readdir(DH)) {
    next if($dline !~ /\.txt$/);
    my $wgsId = $dline;
    $wgsId =~ s/\.txt$//g;
    next if(defined($sawAssembly{$wgsId}));
    print STDERR "...deleting file $dline...\n";
    unlink("$wgsMetaDir/$dline");
    $numDel++;
  }
  closedir(DH);
}
print STDERR "...deleted $numDel obsolete files.\n";
exit 0;

# Subroutine to see if the metadata file is up to date.
# If so, return 1.  If not, return 0.
sub upToDate($@) {
  my $trgFile = shift @_;
  my @srcFileList = @_;
  my $trgFh, my $srcFh;
  my $trgTime, my $srcTime;

  open $trgFh, $trgFile or return 0;
  $trgTime = stat($trgFh)->mtime;
  close $trgFh;

  foreach my $srcFile (@srcFileList) {
    open $srcFh, $srcFile or return 1;
    $srcTime = stat($srcFh)->mtime;
    close $srcFh;
    if($srcTime >= $trgTime) { return 0; }
  } 

  return 1;
}

# Subroutine to process a genome.
sub processGenome($$$$) {
  my ($id, $fastaFile, $gbkFile, $metaFile) = @_;
  my $sfh, my $metaFh;

  if(open($metaFh, ">$metaFile") == 0) { 
    warn "...couldn't open $metaFile for writing\n";
    return;
  }

  # FASTA files first.  Each file has the following that needs to be
  # calculated: (1) MD5 Checksum, (2) Length, (3) GC Bases, (4) Nonstandard
  # Bases, (5) Number of Gaps.
  my $sequence = "", my @checksum = (), my @seqLen = ();
  my @gcBases = (), my @nonBases = (), my @numGaps = ();
  my @giNum = (), my @accNum = (); my $seqCtr = -1;

  if(open($sfh, $fastaFile) == 0) {
    warn "...failed to open file $fastaFile for reading\n";
    return;
  }
  while(my $line = <$sfh>) {
    if($line =~ /^>/ && $seqCtr != -1) {
      $seqLen[$seqCtr] = length($sequence);
      $checksum[$seqCtr] = md5_hex($sequence);
      $gcBases[$seqCtr] = ($sequence =~ tr/GC/GC/);
      $nonBases[$seqCtr] = $seqLen[$seqCtr] - ($sequence =~ tr/ACTG/ACTG/);
      $numGaps[$seqCtr] = scalar (my @tmp = 
                          ($sequence =~ /[^N]NNNNNNNNNNN*[^N]/g));
    }
    if($line =~ /^>/) {
      $seqCtr++;
      $sequence = "";
      my @info = split /[\|>]+/, $line;
      if(!defined($info[2])) { $giNum[$seqCtr] = "#N/A"; }
      else { $giNum[$seqCtr] = $info[2]; }
      if(!defined($info[4])) { $accNum[$seqCtr] = "#N/A"; }
      else { $accNum[$seqCtr] = $info[4]; }
      next;
    }
    chop $line;
    $line =~ tr/a-z/A-Z/;
    $sequence .= $line;
  }
  close $sfh;
  $seqLen[$seqCtr] = length($sequence);
  $checksum[$seqCtr] = md5_hex($sequence);
  $gcBases[$seqCtr] = ($sequence =~ tr/GC/GC/);
  $nonBases[$seqCtr] = $seqLen[$seqCtr] - ($sequence =~ tr/ACTG/ACTG/);
  $numGaps[$seqCtr] = scalar (my @tmp = ($sequence =~ /[^N]NNNNNNNNNNN*[^N]/g));

  # Genbank files second.
  my $gbkText = "", my @dateMod = (), my @sourceOrg = ();
  my @taxonomy = (), my @taxonId = (), my $reading = 0;
  my @definition = (), my @repliconType = (), my @repliconDesc = ();
  $seqCtr = -1;
  if(open($sfh, $gbkFile) == 0) {
    warn "...failed to open file $gbkFile for reading\n";
    return;
  }
  while(my $gbkLine = <$sfh>) {
    if($gbkLine =~ /^\/\//) {
      $dateMod[$seqCtr] = ($gbkText =~ /LOCUS.+?\s(\S+)\n/sg)[0];
      if(!defined($dateMod[$seqCtr])) { $dateMod[$seqCtr] = "#N/A"; }
      $sourceOrg[$seqCtr] = ($gbkText =~ /\nSOURCE\s+(\S.+?)\n  ORGA/sg)[0];
      if(!defined($sourceOrg[$seqCtr])) { $sourceOrg[$seqCtr] = "#N/A"; }
      else { $sourceOrg[$seqCtr] =~ s/\n\s*/ /g; }
      $taxonomy[$seqCtr] = ($gbkText =~ /\n  ORGANISM\s+(\S.+?)\n[A-Z]/sg)[0];
      if(!defined($taxonomy[$seqCtr])) { $taxonomy[$seqCtr] = "#N/A"; }
      else {
        $taxonomy[$seqCtr] =~ s/\n\s*/ /g;
        my $ndx = index($taxonomy[$seqCtr], $sourceOrg[$seqCtr]);
        if($ndx != -1) {
          $ndx += length($sourceOrg[$seqCtr]) + 1;
          $taxonomy[$seqCtr] = substr($taxonomy[$seqCtr], $ndx);
          $taxonomy[$seqCtr] =~ s/^\s+//g;
        }
        else {
          my $tmpSource = $sourceOrg[$seqCtr];
          $tmpSource =~ s/\s*\([^\)]+\)\s*$//g;
          $ndx = index($taxonomy[$seqCtr], $tmpSource);
          if($ndx != -1) {
            $ndx += length($tmpSource) + 1;
            $taxonomy[$seqCtr] = substr($taxonomy[$seqCtr], $ndx);
            $taxonomy[$seqCtr] =~ s/^\s+//g;
          }
        }
      }
      $taxonId[$seqCtr] = ($gbkText =~ /db\_xref=\"taxon\:(\d+)\"/)[0];
      if(!defined($taxonId[$seqCtr])) { $taxonId[$seqCtr] = "#N/A"; }
      $definition[$seqCtr] = ($gbkText =~ /\nDEFINITION\s+(\S.+?)\n[A-Z]/sg)[0];
      if(!defined($definition[$seqCtr])) { $definition[$seqCtr] = "#N/A"; }
      else { $definition[$seqCtr] =~ s/\n\s*/ /g; }

      print $metaFh "$id\t$checksum[$seqCtr]\t$accNum[$seqCtr]\t";
      print $metaFh "$giNum[$seqCtr]\t$sourceOrg[$seqCtr]\t";
      print $metaFh "$dateMod[$seqCtr]\t$taxonId[$seqCtr]\t$taxonomy[$seqCtr]";
      print $metaFh "\t$seqLen[$seqCtr]\t$gcBases[$seqCtr]\t$nonBases[$seqCtr]";
      print $metaFh "\t$numGaps[$seqCtr]\t$definition[$seqCtr]\n";
    }
    if($gbkLine =~ /^LOCUS /) {
      $seqCtr++;
      $gbkText = "";
      $reading = 1;
    }
    if($reading == 1 && $gbkLine =~ /^ORIGIN/) { $reading = 0; next; }
    if($reading == 1) {
      $gbkLine =~ s/\r//g;
      $gbkText .= $gbkLine;
    }
  }
  close $metaFh;
  $seqCtr++;

  # Output final information to summary file
  my @sortedChecksum = (), my $finalChecksum;
  my $finalSeqLen = 0, my $finalGC = 0.0, my $finalNonBases = 0;
  my $finalNumGaps = 0, my $catMD5; my $keep = "Y";

  if($taxonomy[0] !~ /^Bacteria/ && $taxonomy[0] !~ /^Archaea/) { $keep = "N"; }
  for(my $i = 0; $i < $seqCtr; $i++) {
    $finalSeqLen += $seqLen[$i];
    $finalGC += $gcBases[$i];
    $finalNonBases += $nonBases[$i];
    $finalNumGaps += $numGaps[$i];
    if($sourceOrg[$i] ne $sourceOrg[0] && index($sourceOrg[0], 
       $sourceOrg[$i]) == -1 && index($sourceOrg[$i], $sourceOrg[0]) == -1) { 
      $keep = "N"; 
    }
    if($taxonomy[$i] ne $taxonomy[0]) { $keep = "N"; }
    if($taxonId[$i] ne $taxonId[0] &&
      ($taxonId[$i] ne "#N/A" || $taxonId[0] ne "#N/A")) { $keep = "N"; }
  }

  $finalGC = sprintf "%.2f", 100.0*$finalGC/$finalSeqLen;
  @sortedChecksum = sort { $a cmp $b } @checksum;
  $catMD5 = join ',', @sortedChecksum;
  $finalChecksum = md5_hex($catMD5);
  my $retString = "$id\t$finalChecksum\t$seqCtr\t";
  $retString .= "$sourceOrg[0]\t$taxonomy[0]\t$finalSeqLen\t$finalGC\t";
  $retString .= "$finalNonBases\t$finalNumGaps\t$keep";
  return $retString;
}

# Print help message and exit
sub help($) {
  my $version = @_;
  print STDERR "\nNCBI Prokaryotic Genome Processing Script $version\n\n";
  print STDERR "processGenbank.pl [-c] [-d] [-h] [-r <root dir>] [-v] [-w]\n\n";
  print STDERR "-c,--complete:   Only process complete genomes.\n";
  print STDERR "-d,--doall   :   Reprocess everything (default: ";
  print STDERR "only do incremental update.\n";
  print STDERR "-h,--help    :   Print help information and exit.\n";
  print STDERR "-r,--root    :   Specify root directory.\n";
  print STDERR "-v,--version :   Print version information and exit.\n";
  print STDERR "-w,--wgsonly :   Only process WGS genomes.\n\n";
  exit(0);
}

# Print version and exit
sub version($) {
  my $version = @_;
  print STDERR "NCBI Prokaryotic Genome Downloader $version\n";
  exit(0);
}
