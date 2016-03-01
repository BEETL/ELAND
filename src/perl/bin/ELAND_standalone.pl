#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;

# PROJECT: ELAND (Efficient Large-Scale Alignment of Nucleotide Databases)
# MODULE: ELAND_standalone.pl
# AUTHOR: A. J. Cox
#
# This software is covered by the "Illumina Genome Analyzer Software
# License Agreement" and the "Illumina Source Code License Agreement",
# and certain third party copyright/licenses, and any user of this
# source file is bound by the terms therein (see accompanying files
# Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
# Illumina_Source_Code_License_Agreement.pdf and third party
# copyright/license notices).

# Driver script to run ELAND analysis as done in the version 0.3 pipeline
# in a standalone fashion

use Carp;
use Cwd qw(abs_path);
use Data::Dumper;
use File::Basename;
use File::Temp;
use File::Spec;
use File::Copy;
use File::Path;
use POSIX qw(strftime);

use FindBin;
use lib "$FindBin::Bin/../src/perl/lib/";
use Casava::Common::Qseq;
use Casava::Common::QseqRead;
use Casava::Common::Utils;

sub executeCommand($$);
sub errorExit($);
sub getTimestamp;
sub printLog($);

# --------------------
# command line parsing
# --------------------

my $banner =<< "END";
===========================================================================
ELAND_standalone.pl - runs the alignment pipeline on supplied read files
CASAVA 1.8.0                                        (C) 2011 Illumina, Inc.
===========================================================================

END

my $usage =<< "END";
Required options:
  -if [ --input-file ] arg              the file containing the reads. Use
                                        twice for paired-end data sets.
                                        
  -ref [ --ref-sequences ] arg          path to the reference sequence
                                        directory.
                                                                                
Options:
  --bam                                 enables BAM output
 
  -bq [ --base-quality] arg             the base quality when parsing
      (=30)                             FASTA read files.
      
  -cr [ --copy-references ]             copies the references to the output
                                        directory. Use this option if your
                                        reference sequence directory is
                                        write-protected.
                                        
  --force                               forces existing output files to be
                                        overwritten.
 
  -it [ --input-type ] arg              the type of read files provided:
                                        'fasta', 'fastq', 'export', 'qseq'
                                        
  -l [ --log ] arg                      the path to the log file.
     (=ELAND_standalone.log)
     
  -od [ --output-directory] arg         the output directory.
  
  -op [ --output-prefix ] arg           the prefix for all of the files
      (=reanalysis)                     created by ELAND_standalone.
 
  -ko [ --kagu-options ] arg            additional options to kagu when
                                        processing paired-end data sets.
                                        e.g. -ko "-c" enables circular
                                        reference sequence support.
                                        
  -rt [ --remove-temps ]                removes all files except exports,
                                        BAM files, and log files upon
                                        successful completion.
  
  -sl [ --seed-length ] arg             the ELAND seed length. Use twice
      (=min(read length, 32))           for paired-end data sets.
      
  -ub [ --use-bases ] arg               the use-bases masking string. Use
      (=Y*n)                            twice for paired-end data sets.
 
Help:
  -h [ --help ]                         shows this help text.
END

print $banner;
die $usage if (@ARGV==0);
local $| = 1;

# global variables
my $maxSeedLength=32;
my $timeStampFormat = "%Y-%m-%d %H:%M:%S";

my $pipelineBinDir     = "$FindBin::Bin";
my $pipelineLibexecDir = "$FindBin::Bin";
#my $pipelineBinDir     = '@CASAVA_FULL_BINDIR@';
#my $pipelineLibexecDir = '@CASAVA_FULL_LIBEXECDIR@';

my @tempFiles = ();

my %options=(
    'base-quality'     => 30,
    'output-directory' => '.',
    'output-prefix'    => 'reanalysis',
    'index'            => 'NoIndex',
    'lane'             => 1);

GetOptions(\%options,
           "output-directory|od=s",
           "input-file|if=s@",
           "input-type|it=s",
           "seed-length|sl=i@",
           "ref-sequences|ref=s",
           "output-prefix|op=s",
           "base-quality|bq=i",
           "kagu-options|ko=s",
           "use-bases|ub=s@",
           "copy-references|cr",
           "bam",
           "log|l=s",
           "test=s",
           "force",
           "remove-temps|rt",
           "help");

# display the help menu
if(defined($options{'help'})) {
    print $usage;
    exit(1);
}

my $forceOutput = (defined($options{'force'}) ? 1 : 0);

# check that ref-sequences is defined
unless(defined($options{'ref-sequences'})) {
    print STDERR "ERROR: The reference sequence directory was not provided. Please use the --ref-sequences parameter.\n";
    exit(-1);
}

# retrieve the absolute path for the output directory
$options{'ref-sequences'}    = abs_path($options{'ref-sequences'});
$options{'output-directory'} = abs_path($options{'output-directory'});

print "Output files will be stored in: " . $options{'output-directory'} . "\n\n";

# create the output directory
unless(-e $options{'output-directory'}) {
    mkdir($options{'output-directory'}) || die("ERROR: Unable to create the output directory (" . $options{'output-directory'} . ").");
}

# check if the reference directory exists
unless (-e $options{'ref-sequences'}) {
    print STDERR "ERROR: The reference sequence directory (" . $options{'ref-sequences'} . ") does not exist.\n";
    exit(-1);
}

my $tempReferenceDir = File::Spec->catdir($options{'output-directory'}, "references");
my $copyReferences = (defined($options{'copy-references'}) ? 1 : 0);

# open the log
unless(defined($options{'log'})) {
    $options{'log'} = File::Spec->catfile($options{'output-directory'}, "ELAND_standalone.log");
}

print "Logging output to: " . $options{'log'} . "\n\n";
open(LOG, ">>".$options{'log'}) || errorExit("Could not open the log file.");
$|++;

print LOG $banner;
print LOG getTimestamp() . "started.\n\n";

# dump the chosen options
print LOG "Command line options:\n";
print LOG Dumper(\%options) . "\n";

# -----------
# subroutines
# -----------

# sub autoSquashReferenceSequences()
# checks the supplied directory for the presence of squashed reference
# sequences. If none are found, the references is automatically squashed.

sub autoSquashReferenceSequences() {
    
    # retrieve the filenames from the references directory
    opendir(DIR, $options{'ref-sequences'});
    my @files = readdir(DIR);
    closedir(DIR);
    
    # grab the squashed filenames
    my @twoBitFiles = grep(/.2bpb$/, @files);
    my @idxFiles    = grep(/.idx$/,  @files);
    my @vldFiles    = grep(/.vld$/,  @files);
    
    my $num2bpbFiles = scalar(@twoBitFiles);
    my $numIdxFiles  = scalar(@idxFiles);
    my $numVldFiles  = scalar(@vldFiles);
    
    # autosquash if we're missing any of the file types or if we have unequal
    # numbers of squashed files
    my $needAutosquash = 0;
    if(($num2bpbFiles == 0) || ($num2bpbFiles == 0) || ($num2bpbFiles == 0)) {
        $needAutosquash = 1;
    }
    
    unless(($num2bpbFiles == $numIdxFiles) && ($numIdxFiles == $numVldFiles)) {
        $needAutosquash = 1;
    }
    
    # perform the autosquashing
    if($needAutosquash) {
        
        # check if the reference sequence directory is writable
        if(($copyReferences == 0) && !(-w $options{'ref-sequences'})) {
            printStatus("WARNING: The reference sequence directory is not writable.\n");
            printStatus("         References will be copied and autosquashed in:\n");
            printStatus("         " . $tempReferenceDir . "\n");
            $copyReferences = 1;
        }
        
        # look for the FASTA files
        my @fastaFilenames = grep(/(.fa$|.fasta$)/, @files);
        my @fastaPaths  = map { $options{'ref-sequences'}."/$_" } @fastaFilenames;
        
        errorExit("No FASTA files (*.fa or *.fasta) were found in the reference directory (" . $options{'ref-sequences'} ."). Unable to automatically squash the reference files.") if (scalar(@fastaPaths) == 0);
    
        # create the temporary references directory
        if($copyReferences == 1) {
            
            # create the temporary references directory
            unless(-e $tempReferenceDir) {
                mkdir($tempReferenceDir) || errorExit("Unable to create the temporary reference sequence directory (" . $tempReferenceDir . ").");
            }
            
            # populate the temporary directory with symbolic links to the FASTA files
            my @newPaths = map { "$tempReferenceDir/$_" } @fastaFilenames;
            
            for(my $j = 0; $j < scalar(@fastaPaths); $j++) {
                unless(-e $newPaths[$j]) {
                    symlink($fastaPaths[$j], $newPaths[$j]) || errorExit("Unable to create a symbolic link from (" . $fastaPaths[$j] . ") to (" . $newPaths[$j] . ").");
                }
            }
            
            @fastaPaths = @newPaths;
        
            # update the reference directory
            $options{'ref-sequences'} = $tempReferenceDir;
        }
        
        # run squashGenome on each FASTA file
        printStatus("automatically squashing reference sequences... ");
        my $program = File::Spec->catfile($pipelineLibexecDir, "squashGenome");
        
        foreach my $file (@fastaPaths) {
            my $command = "${program} --chrom-name-source contigName " . $options{'ref-sequences'} . " ${file} 2>&1";
            executeCommand("squashGenome", $command);
        }
        
        print "finished.\n";
    }
}

# sub convertReadsToFastq($$$)
# Takes input file in fasta/fastq/export format and makes
# FASTQ files.

sub convertReadsToFastq($$$$) {
    
    croak "ERROR: convertReadsToFastq\n" unless (@_ == 4);
    my ($inputFilename,
        $inputFileFormat,
        $fastqFilename,
        $readNum) = @_;

    # check if our input filename exists
    unless (-e $inputFilename) {
        print STDERR "ERROR: The specified input file (" . $inputFilename . ") does not exist.\n";
        exit(-1);
    }
    
    # check if our output filename exists
    if((! -e $fastqFilename) || $forceOutput) {
        
        printStatus("converting ${inputFilename}... ");
        my $program = File::Spec->catfile($pipelineBinDir, "FastqConverter");
        my $tempFilename = File::Spec->catfile($options{'output-directory'}, "temp.fastq.gz");
        my $command= "${program} --in ${inputFilename} --out ${tempFilename} --read ${readNum}";
        
        # add the input type if specified
        if($inputFileFormat) {
            $command = $command . " --it ${inputFileFormat}";
            
            # add the default base qualities
            if($inputFileFormat eq "fasta") {
                $command = $command . " --bq " . $options{'base-quality'}
            }
        }
        
        # run the FASTQ converter
        executeCommand("FastqConverter", $command . " 2>&1");
        
        # rename the output file
        rename($tempFilename, $fastqFilename);
        push(@tempFiles, $fastqFilename);
        
        print "finished.\n";
    
    } else {
        
        printStatus("skipping FASTQ conversion since ${fastqFilename} already exists.\n");
    }
    
    # derive the read length
    return getFastqReadLength($fastqFilename);
}

sub executeCommand($$) {
    croak "ERROR: executeCommand\n" unless (@_ == 2);
    my($programName, $command) = @_;
    
    print LOG getTimestamp() . "About to run $command\n";
    
    my $output = qx{$command};
    my $status = $? >> 8;
    
    print LOG $output;
    
    if($status == 0) {
        print LOG getTimestamp() . "${programName} succeeded\n";
    } else {
        errorExit("An error occurred while running ${programName} (error code: ${status})");
    }
}

sub errorExit($) {
    croak "ERROR: errorExit\n" unless ( @_ == 1 );
    my($message) = @_;
    
    use Term::ANSIColor qw (:constants);
    print STDERR BOLD RED;
    print STDERR "ERROR: ";
    print STDERR RESET;
    print STDERR "$message\n";
    exit(1);
}

# sub getFastqReadLength($)
# Opens a FASTQ file using the supplied FASTQ filename and returns the read
# length of the first read.

sub getFastqReadLength($) {
    croak "ERROR: getFastqReadLength\n" unless ( @_ == 1 );
    my($filename) = @_;
    
    # retrieve the second line
    open(FASTQ, "gunzip -fc $filename |" ) || errorExit("Could not open fastq file $filename");
    my $line = <FASTQ>;
    chomp($line = <FASTQ>);
    close(FASTQ);
    
    return length($line);
}

# sub getTimestamp
# returns the current timestamp

sub getTimestamp {
    return "[" . strftime($timeStampFormat, localtime) . "] ";
}

# sub inspectInputFiles
# retrieve the lane and index by inspecting the input files or by
# checking the output directory

sub inspectInputFiles {
    
    # initialize
    my @lines    = ();
    my @colSizes = ();
    my $seqFormat  = "unknown";

    # sanity check
    my $inputFilename = $options{'input-file'}->[0];
    errorExit("The input filename ($inputFilename) doesn't exist.") unless(-e $inputFilename);

    # -------------------------------------------------------------------------
    # figure out which input file format we have (uses the same algorithm as in
    # lib/common/FileConversion.cpp)
    # -------------------------------------------------------------------------
    
    # read the three lines from the input filename
    open(IN, "gunzip -fc $inputFilename |") || die("Unable to open the input file ($inputFilename) for reading.");
    
    for(my $i = 0; $i < 3; $i++) {
        my $line = <IN>;
        last unless defined $line;
        chomp($line);
        push(@lines, $line);
        my @columns = split('\t', $line);
        push(@colSizes, scalar(@columns));
    }
    
    close(IN);
    
    my $numLines = scalar(@lines);
    
    # if the file starts with a greater than symbol, it's probably FASTA
    if(($numLines >= 2) && ($lines[0] =~ /^>/) && ($lines[1] =~ /^[acgtnACGTN]/)) {
        $seqFormat = "fasta";
    }

    # if the file has an at symbol and a plus symbol, it's probably FASTQ
    if(($numLines >= 3) && ($lines[0] =~ /^@/) && ($lines[2] =~ /^\+/)) {
        $seqFormat = "fastq";
    }

    # check for a consistent number of columns
    my $isQseq   = 1;
    my $isExport = 1;
    
    foreach my $colSize (@colSizes) {
        
        if($isQseq   && ($colSize != 11)) {
            $isQseq = 0;
        }
        
        if($isExport && ($colSize != 22)) {
            $isExport = 0;
        }
    }

    # set the QSEQ and export formats
    if($isQseq) {
        $seqFormat = "qseq";
    }
    
    if($isExport) {
        $seqFormat = "export";
    }
    
    # ----------------------------------------------
    # extract information from our files if possible
    # ----------------------------------------------
    
    # extract from qseq or export
    if(($seqFormat eq "qseq") || ($seqFormat eq "export")) {
        my @columns = split('\t', $lines[0]);
        $options{'lane'}  = $columns[2] if ($columns[2]);
        $options{'index'} = $columns[6] unless ($columns[6] eq "0");
    }
       
    # extract from FASTQ (if using CASAVA 1.8 format)
    if($seqFormat eq "fastq") {

        my @fastqFields = ($lines[0] =~ /^@([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):(\S+)\s+([^:]+):([^:]+):([^:]+):(.*)$/);
        
        # the FASTQ file is in CASAVA 1.8 format
        if(scalar(@fastqFields) == 11) {
            my $lane  = $fastqFields[3];
            my $index = $fastqFields[10];
            $options{'lane'}  = $lane  if ($lane);
            $options{'index'} = $index if ($index);
        }
    }
}

# sub makeBamFiles($$\@)
# converts the export files into BAM files

sub makeBamFiles($$\@) {

    croak "ERROR: makeBamFiles\n" unless (@_ == 3);
    my ($filenameStub, $genomeSizeFilename, $exportFilesRef) = @_;
    my @exportFiles = @$exportFilesRef;
    my $numExportFiles = scalar(@exportFiles);
    
    # check if we have the right number of export files
    if(($numExportFiles < 1) || ($numExportFiles > 2)) {
        errorExit("Expected one or two export files in makeBamFiles, but received $numExportFiles");
    }
    
    # define our filenames
    my $samFilename         = $filenameStub . ".sam";
    my $unsortedBamFilename = $filenameStub . "_unsorted.bam";
    my $bamFilename         = $filenameStub . ".bam";

    if((! -e $bamFilename) || $forceOutput) {
        
        # -------------------
        # create the SAM file
        # -------------------
    
        # check if our SAM filename exists
        if((! -e $samFilename) || $forceOutput) {
    
            # create the SAM header
            printStatus("creating SAM file... ");
            my $command = "grep chromosome ${genomeSizeFilename} | awk 'BEGIN { FS=\"\\\"\" } { print \"\@SQ\tSN:\" \$4 \"\tLN:\" \$6 }' > ${samFilename}.tmp";
            executeCommand("SAM header creation", $command);
            
            # convert the export file
            my $program = File::Spec->catfile($pipelineBinDir, "export2sam.pl");
            $command = "${program} --read1=" . $exportFiles[0];
            $command = $command . " --read2=" . $exportFiles[1] if($numExportFiles == 2);
            $command = $command . " 1>> ${samFilename}.tmp 2> /dev/null";
            
            # create the genome size xml file
            executeCommand("illumina_export2sam.pl", $command);
            
            # rename the output file
            rename($samFilename . ".tmp", $samFilename);
            push(@tempFiles, $samFilename);
            
            print "finished.\n";
            
        } else {
            
            printStatus("skipping SAM file creation since ${samFilename} already exists.\n");
        }
        
        # ------------------------
        # create unsorted BAM file
        # ------------------------
        
        # check if our unsorted BAM filename exists
        if((! -e $unsortedBamFilename) || $forceOutput) {
            
            printStatus("creating BAM file... ");
            my $program = File::Spec->catfile($pipelineLibexecDir, "samtools");
            my $command = "${program} view -bS ${samFilename} 1> ${unsortedBamFilename}.tmp 2> /dev/null";
            
            # create the uncompressed file
            executeCommand("samtools", $command);
            
            # rename the output file
            rename($unsortedBamFilename . ".tmp", $unsortedBamFilename);
            push(@tempFiles, $unsortedBamFilename);
            
            print "finished.\n";
            
        } else {
            
            printStatus("skipping BAM file creation since ${unsortedBamFilename} already exists.\n");
        }

        # -----------------
        # sort the BAM file
        # -----------------
        
        printStatus("sorting BAM file... ");
        my $program = File::Spec->catfile($pipelineLibexecDir, "samtools");
        my $command = "${program} sort ${unsortedBamFilename} ${filenameStub}_tmp 2>&1";
        
        # create the uncompressed file
        executeCommand("samtools", $command);
        
        # rename the output file
        rename(${filenameStub} . "_tmp.bam", $bamFilename);
        
        print "finished.\n";
        
        # -------
        # cleanup
        # -------
        
        unlink($unsortedBamFilename);
        unlink($samFilename);
        
    } else {
        
        printStatus("skipping BAM sorting since ${bamFilename} already exists.\n");
    }
}

# sub makeElandExtended($$$)
# Takes the processing as far as the eland_extended file
# This bit is common to both single read and paired read code
# Assumes read length is already either specified or worked out by
# Works out alignment seed length if necessary

sub makeElandExtended($$$$$$$$$) {
    croak "ERROR: makeElandExtended\n" unless (@_ == 9);
    my ($pipelineDir, $seedLength, $readNum, $useBases, $sampleName, $fastqDir,
        $genomeDir, $elandExtendedFilename, $readLength) = @_;
    
    # derive the seed length
    if ($seedLength==-1) {
        $seedLength = ($readLength > $maxSeedLength) ? $maxSeedLength : $readLength;
    } else {
        $seedLength = ($seedLength > $maxSeedLength) ? $maxSeedLength : $seedLength;
        $seedLength = ($seedLength > $readLength)    ? $readLength    : $seedLength;
    }

    print LOG getTimestamp() . "Setting seed length for $elandExtendedFilename to $seedLength\n";

    # check if our output filename exists
    if((! -e $elandExtendedFilename) || $forceOutput) {
                
        printStatus("aligning read ${readNum}... ");
        my $program = File::Spec->catfile($pipelineBinDir, "eland_ms");
        my $command = "${program} --oligo-length=${seedLength} --data-format fastq --lane " . $options{'lane'} ." --read ${readNum} --qseq-mask ${useBases} --cluster-sets 001 --sample ${sampleName} --barcode " . $options{'index'} . " --base-calls-dir ${fastqDir} --genome-directory ${genomeDir} --output-file ${elandExtendedFilename} --multi 2>&1";
        
        # align our reads
        executeCommand("ELAND", $command);
        push(@tempFiles, $elandExtendedFilename);
        
        print "finished.\n";
        
    } else {
        
        printStatus("skipping alignment since ${elandExtendedFilename} already exists.\n");
    }
    
    # return the seed length
    return $seedLength;
}

sub makeExportPaired($$$$$$$$$$$$$) {

    croak "ERROR: makeExportPaired\n" unless (@_ == 13);
    my ($refSeqDir, $genomeSizeFilename, $mate1ElandExtendedFilename,
        $mate2ElandExtendedFilename, $mate1FastqFilename, $mate2FastqFilename,
        $mate1ExportFilename, $mate2ExportFilename, $pairXmlFilename,
        $mate1UseBases, $mate2UseBases, $mate1SeedLength,
        $mate2SeedLength) = @_;

    # ---------------------------
    # create genome size xml file
    # ---------------------------

    # check if our genome size filename exists
    if((! -e $genomeSizeFilename) || $forceOutput) {
        
        printStatus("creating genome size xml file... ");
        my $program = File::Spec->catfile($pipelineLibexecDir, "squashGenome");
        my $command = "${program} --chrom-name-source contigName ${refSeqDir} 1> ${genomeSizeFilename}.tmp 2> /dev/null";

        # create the genome size xml file
        executeCommand("squashGenome", $command);
        
        # rename the output file
        rename($genomeSizeFilename . ".tmp", $genomeSizeFilename);
        push(@tempFiles, $genomeSizeFilename);
        
        print "finished.\n";
        
    } else {
        
        printStatus("skipping genome size file creation since ${genomeSizeFilename} already exists.\n");
    }
    
    # ----------------------
    # run the orphan aligner
    # ----------------------
    
    my $orphanMate1ElandExtendedFilename = $mate1ElandExtendedFilename . ".oa";
    my $orphanMate2ElandExtendedFilename = $mate2ElandExtendedFilename . ".oa";
    
    if((! -e $orphanMate1ElandExtendedFilename) || (! -e $orphanMate2ElandExtendedFilename) || $forceOutput) {
        
        printStatus("running the orphan aligner... ");
        my $program = File::Spec->catfile($pipelineBinDir, "orphanAligner");
        my $command = "${program} ${mate1ElandExtendedFilename} ${mate2ElandExtendedFilename} ${refSeqDir} .oa.tmp 10 2>&1";
    
        # run the orphan aligner
        executeCommand("orphanAligner", $command);
            
        # rename the output files
        rename($orphanMate1ElandExtendedFilename . ".tmp", $orphanMate1ElandExtendedFilename);
        rename($orphanMate2ElandExtendedFilename . ".tmp", $orphanMate2ElandExtendedFilename);
        push(@tempFiles, $orphanMate1ElandExtendedFilename);
        push(@tempFiles, $orphanMate2ElandExtendedFilename);
            
        print "finished.\n";
        
    } else {
        
        printStatus("skipping orphan aligner since ${orphanMate1ElandExtendedFilename} and ${orphanMate2ElandExtendedFilename} already exist.\n");
    }

    # -----------------------
    # create our export files
    # -----------------------
    
    if((! -e $mate1ExportFilename) || (! -e $mate2ExportFilename) || $forceOutput) {
        
        printStatus("creating export files... ");
        my $program = File::Spec->catfile($pipelineLibexecDir, "kagu");
        my $command = "${program} --ie1 ${orphanMate1ElandExtendedFilename} --ie2 ${orphanMate2ElandExtendedFilename} --if1 ${mate1FastqFilename} --if2 ${mate2FastqFilename} --irs ${genomeSizeFilename} --oe1 ${mate1ExportFilename}.tmp.gz --oe2 ${mate2ExportFilename}.tmp.gz --os ${pairXmlFilename}.tmp --ub1 ${mate1UseBases} --ub2 ${mate2UseBases} --sl1 ${mate1SeedLength} --sl2 ${mate2SeedLength} --ucn ";
        
        # add some additional options
        if(defined($options{'kagu-options'})) {
            $command = $command . $options{'kagu-options'} . " ";
        }
        
        $command = $command . "2>&1";
        
        # create the export file
        executeCommand("kagu", $command);
        
        # rename the output files
        rename($mate1ExportFilename . ".tmp.gz", $mate1ExportFilename);
        rename($mate2ExportFilename . ".tmp.gz", $mate2ExportFilename);
        rename($pairXmlFilename     . ".tmp",    $pairXmlFilename);
        
        print "finished.\n";
        
    } else {
        
        printStatus("skipping export file creation since ${mate1ExportFilename} and ${mate2ExportFilename} already exist.\n");
    }
}

sub makeExportSingle($$$$$$$) {
    
    croak "ERROR: makeExportSingle\n" unless (@_ == 7);
    my ($refSeqDir, $genomeSizeFilename, $elandExtendedFilename, $fastqFilename,
        $exportFilename, $useBases, $seedLength)=@_;

    # ---------------------------
    # create genome size xml file
    # ---------------------------

    # check if our genome size filename exists
    if((! -e $genomeSizeFilename) || $forceOutput) {
        
        printStatus("creating genome size xml file... ");
        my $program = File::Spec->catfile($pipelineLibexecDir, "squashGenome");
        my $command = "${program} ${refSeqDir} 1> ${genomeSizeFilename}.tmp 2> /dev/null";

        # create the genome size xml file
        executeCommand("squashGenome", $command);
        
        # rename the output file
        rename($genomeSizeFilename . ".tmp", $genomeSizeFilename);
        push(@tempFiles, $genomeSizeFilename);
        
        print "finished.\n";
        
    } else {
        
        printStatus("skipping genome size file creation since ${genomeSizeFilename} already exists.\n");
    }
    
    # ----------------------
    # create our export file
    # ----------------------
    
    if((! -e $exportFilename) || $forceOutput) {
        
        printStatus("creating export file... ");
        my $program = File::Spec->catfile($pipelineLibexecDir, "kagu");
        my $command = "${program} --ie1 ${elandExtendedFilename} --if1 ${fastqFilename} --irs ${genomeSizeFilename} --oe1 ${exportFilename}.tmp.gz --ub1 ${useBases} --sl1 ${seedLength} --ucn ";
        
        # add some additional options
        if(defined($options{'kagu-options'})) {
            $command = $command . $options{'kagu-options'} . " ";
        }
        
        $command = $command . "2>&1";
        
        # create the export file
        executeCommand("kagu", $command);
        
        # rename the output file
        rename($exportFilename . ".tmp.gz", $exportFilename);
        
        print "finished.\n";
        
    } else {
        
        printStatus("skipping export file creation since ${exportFilename} already exists.\n");
    }
}

# sub printHeader($)
# prints text in another color in order to highlight the text

sub printHeader($) {
    croak "ERROR: printLog\n" unless ( @_ == 1 );
    my($message) = @_;
    use Term::ANSIColor qw (:constants);
    print BOLD;
    print $message;
    print RESET;
}

# sub printStatus($)
# prints a timestamped status update

sub printStatus($) {
    croak "ERROR: printLog\n" unless ( @_ == 1 );
    my($message) = @_;
    use Term::ANSIColor qw (:constants);
    print DARK GREEN;
    print getTimestamp();
    print RESET;
    print $message;
    print LOG getTimestamp() . $message;
}

# sub removeTempFiles
# removes temporary files after successful completion of ELAND_standalone

sub removeTempFiles() {
    foreach my $file (@tempFiles) {
        print LOG getTimestamp() . "Deleting temporary file: ${file}\n";
        unlink($file);
    }
}

# ==========================
# check for required options
# ==========================

# check the input files
my $foundErrors   = 0;
my $numInputFiles = @{$options{'input-file'}};

unless(defined($options{'input-file'})) {
    print STDERR "ERROR: An input file was not specified. Please use the -if parameter.\n";
    $foundErrors = 1;
} else {
    # TODO: check the number of input files - only 1 or 2
    for(my $j = 0; $j < $numInputFiles; $j++) {
        unless(-e $options{'input-file'}->[$j]) {
            print STDERR "ERROR: The specified input file (" . $options{'input-file'}->[$j] . ") does not exist.\n";
            $foundErrors = 1;
        }
    }
}

# check the ELAND genome directory
unless(defined($options{'ref-sequences'})) {
    print STDERR "ERROR: A reference sequence directory was not specified. Please use the -ref parameter.\n";
    $foundErrors = 1;
} else {
    unless(-e $options{'ref-sequences'}) {
        print STDERR "ERROR: The specified reference sequence directory (" . $options{'ref-sequences'} . ") does not exist.\n";
        $foundErrors = 1;
    }
}

# convert the input type to lower case
if(defined($options{'input-type'})) {
    my $inputType = lc($options{'input-type'});
    $options{'input-type'} = $inputType;
    
    unless(($inputType eq "fasta") || ($inputType eq "fastq") || ($inputType eq "export") || ($inputType eq "qseq")) {
        print STDERR "ERROR: An unknown input type was encountered ($inputType). Please specify either 'fasta', 'fastq', 'export', or 'qseq'.\n";
        $foundErrors = 1;
    }
}

# display the help text
if($foundErrors == 1) {
    print "\n" . $usage;
    exit(-1);
}

# ================================================
# run either the single-end or paired-end workflow
# ================================================

# set the default seed length and use-bases
$options{'seed-length'}[0] = -1    unless (defined($options{'seed-length'}));
$options{'use-bases'}[0]   = "Y*n" unless (defined($options{'use-bases'}));

# retrieve the lane and index by inspecting the input files or by checking the
# output directory
inspectInputFiles();
my $laneString = sprintf("L%03d", $options{'lane'});

if($numInputFiles == 1)  {

    # -------------------
    # single-end workflow
    # -------------------

    printHeader("Running the single-end workflow:\n");
    print LOG getTimestamp() . "single input-file specified, will do single read analysis\n";
    
    # define our filenames
    my $filenameStub = File::Spec->catfile($options{'output-directory'}, $options{'output-prefix'}) . "_" . $options{'index'} . "_" . $laneString;
    my $mate1FastqFilename         = $filenameStub . "_R1_001.fastq.gz";
    my $mate1ElandExtendedFilename = $filenameStub . "_R1_001_eland_extended.txt";
    my $mate1ExportFilename        = $filenameStub . "_R1_001_export.txt.gz";
    my $genomeSizeFilename         = $filenameStub . "_genomesize.xml";
        
    # single-end workflow
    autoSquashReferenceSequences();
    
    $options{'read-length'}[0] = convertReadsToFastq(
        $options{'input-file'}->[0],
        $options{'input-type'},
        $mate1FastqFilename,
        1);
    
    # expand the use-bases string
    $options{'use-bases'}[0] = expandUseBases($options{'use-bases'}[0], $options{'read-length'}[0]);
    $options{'use-bases'}[0] =~ tr/y/Y/;
    print LOG getTimestamp() . "Expanded read 1 use-base string: " . $options{'use-bases'}[0] . "\n";
    
    $options{'seed-length'}[0] = makeElandExtended(
        $options{'pipeline-dir'},
        $options{'seed-length'}[0],
        1,
        $options{'use-bases'}[0],
        $options{'output-prefix'},
        $options{'output-directory'},
        $options{'ref-sequences'},
        $mate1ElandExtendedFilename,
        $options{'read-length'}[0]);

    makeExportSingle(
        $options{'ref-sequences'},
        $genomeSizeFilename,
        $mate1ElandExtendedFilename,
        $mate1FastqFilename,
        $mate1ExportFilename,
        $options{'use-bases'}[0],
        $options{'seed-length'}[0]);
    
    # optionally convert the export files to BAM
    if(defined($options{'bam'})) {
        my @exportFiles = ($mate1ExportFilename);
        makeBamFiles(
            File::Spec->catfile($options{'output-directory'}, $options{'output-prefix'}),
            $genomeSizeFilename,
            @exportFiles);
    }
          
} elsif ($numInputFiles == 2) {
    
    # -----------------------------
    # paired-end/mate-pair workflow
    # -----------------------------

    printHeader("Running the paired-end workflow:\n");
    print LOG getTimestamp() . "two input-files specified, will do paired read analysis\n";
    
    # define our filenames
    my $filenameStub = File::Spec->catfile($options{'output-directory'}, $options{'output-prefix'}) . "_" . $options{'index'} . "_" . $laneString;
    my $mate1FastqFilename         = $filenameStub . "_R1_001.fastq.gz";
    my $mate2FastqFilename         = $filenameStub . "_R2_001.fastq.gz";
    my $mate1ElandExtendedFilename = $filenameStub . "_R1_001_eland_extended.txt";
    my $mate2ElandExtendedFilename = $filenameStub . "_R2_001_eland_extended.txt";
    my $mate1ExportFilename        = $filenameStub . "_R1_001_export.txt.gz";
    my $mate2ExportFilename        = $filenameStub . "_R2_001_export.txt.gz";
    my $genomeSizeFilename         = $filenameStub . "_genomesize.xml";
    my $pairXmlFilename            = $filenameStub . "_001_pair.xml";
    
    # paired-end workflow
    autoSquashReferenceSequences();
        
    $options{'read-length'}[0] = convertReadsToFastq(
        $options{'input-file'}->[0],
        $options{'input-type'},
        $mate1FastqFilename,
        1);
    
    $options{'read-length'}[1] = convertReadsToFastq(
        $options{'input-file'}->[1],
        $options{'input-type'},
        $mate2FastqFilename,
        2);
    
    # adjust the seed length and use-bases
    $options{'use-bases'}[1]   = $options{'use-bases'}[0]   if (@{$options{'use-bases'}} == 1);
    $options{'seed-length'}[1] = $options{'seed-length'}[0] if (@{$options{'seed-length'}} == 1);
    
    # expand the use-bases strings
    for(my $i = 0; $i < $numInputFiles; $i++) {        
        $options{'use-bases'}[$i] = expandUseBases($options{'use-bases'}[$i], $options{'read-length'}[$i]);
        $options{'use-bases'}[$i] =~ tr/y/Y/;
        print LOG getTimestamp() . "Expanded read " . ($i + 1) . " use-base string: " . $options{'use-bases'}[$i] . "\n";
    }
    
    $options{'seed-length'}[0] = makeElandExtended(
        $options{'pipeline-dir'},
        $options{'seed-length'}[0],
        1,
        $options{'use-bases'}[0],
        $options{'output-prefix'},
        $options{'output-directory'},
        $options{'ref-sequences'},
        $mate1ElandExtendedFilename,
        $options{'read-length'}[0]);
    
    $options{'seed-length'}[1] = makeElandExtended(
        $options{'pipeline-dir'},
        $options{'seed-length'}[1],
        2,
        $options{'use-bases'}[1],
        $options{'output-prefix'},
        $options{'output-directory'},
        $options{'ref-sequences'},
        $mate2ElandExtendedFilename,
        $options{'read-length'}[1]);

    makeExportPaired(
        $options{'ref-sequences'},
        $genomeSizeFilename,
        $mate1ElandExtendedFilename,
        $mate2ElandExtendedFilename,
        $mate1FastqFilename,
        $mate2FastqFilename,
        $mate1ExportFilename,
        $mate2ExportFilename,
        $pairXmlFilename,
        $options{'use-bases'}[0],
        $options{'use-bases'}[1],
        $options{'seed-length'}[0],
        $options{'seed-length'}[1]);
    
    # optionally convert the export files to BAM
    if(defined($options{'bam'})) {
        my @exportFiles = ($mate1ExportFilename, $mate2ExportFilename);
        makeBamFiles(
            File::Spec->catfile($options{'output-directory'}, $options{'output-prefix'}),
            $genomeSizeFilename,
            @exportFiles);
    }
}

print "\n";
printStatus("finished.\n");

removeTempFiles() if (defined($options{'remove-temps'}));

close(LOG);
