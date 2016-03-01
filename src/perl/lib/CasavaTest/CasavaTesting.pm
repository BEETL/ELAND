package CasavaTest::CasavaTesting;

BEGIN {
    use Exporter();
    @ISA       = qw(Exporter);
    @EXPORT    = qw(&runCasava $casavaScript);
    @EXPORT_OK = qw();
}
# PROJECT: Pipeline
# MODULE:  Validate.pm
# AUTHOR:  L. Szajkowski
#
# Copyright (c)  2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.

# Functions used to test to test Casava
#

=pod

=head1 NAME

GATest::Validate.pm - Functions used to test to test Casava

=head2 SYNOPSIS

use GATest::CasavaTesting.pm qw();  

=head2 AUTHORSHIP

Copyright (c) 2007 Solexa, 2008 Illumina
This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

=head1 DESCRIPTION

=head2 Overview

Messaging Functions 

=head2 Exports

Global variables:

=head2 Depends on

    warnings, strict, POSIX, XML::Simple, constant, Carp, 
    GATest::TestingTools

=cut

use warnings FATAL => 'all';
use strict;
use POSIX;
use XML::Simple;
use constant DEBUG => 0;    # set to 1 to get debug info
use Data::Dumper;
use Carp;

use File::Copy;
use CasavaTest::TestingTools;

sub runCasava($;$;$;$;$;$;$;$);

our $casavaScript = "bin/CASAVA.pl";

=pod

=head1 FUNCTIONS

=head2 General Functions

=over 4

=cut

=pod

=item runCasava($GA_TEST_DATA_PATH, $testId, $testSetId, $type, 
$testRef, $testSetRef)

The procedure runs CASAVA.pl script for given $testId, $testSetId

B<Parameters:>

    $testsConfRef          - HASH MAP Ref to struture with pipeline_tests.xml
    $GA_TEST_DATA_PATH     - Path to GA_TEST_DATA folder (from command line)
    $ENV_NAME_GA_DIR_PATH  - Path to pipeline installation
    $testId                - Id of test in pipeline_tests.xml
    $testSetId             - Id of test set in pipeline_tests.xml
    $type                  - Type of test [IPAR_1.01 | GA]
    $testRef               - HASH MAP Ref to test
    $testSetRef            - HASH MAP Ref to test set

B<Returns:> 

    status, message
    
=cut
## runGoat, runBustar and run runCasava should be integrated into one procedure
sub runCasava($;$;$;$;$;$;$;$) {
    croak("ERROR: runCasava wrong number of parameters")
      unless ( @_ == 8 );
    my ( $testsConfRef, $GA_TEST_DATA_PATH, $ENV_NAME_GA_DIR_PATH, $testId,
         $testSetId, $type, $testRef, $testSetRef ) = @_;

    my $testSetFolderPath =
      File::Spec->catdir( $ENV_NAME_GA_DIR_PATH, $tempSpace, $testSetId );

    my $currentDirectory = Cwd->getcwd();
    chdir $testSetFolderPath;

    my $geraldFolderRef = getGeraldFolder( $type, %$testRef, %$testSetRef );
    if ( !defined $geraldFolderRef ) {
        die "ERROR: the test $testId is not suported for $testSetId\n";
    }

    my $dirsRef = parseGeraldPath( $geraldFolderRef->{path} );

##  run CASAVA
    my $casavaConfigPath = "";
    my $runFolder        = $dirsRef->{runFolder};
    my $bustardPath      =
      File::Spec->catdir( $testSetFolderPath, $dirsRef->{BustardFull} );
    my $topFolderPath = $bustardPath;


   print "Running CASAVA in $bustardPath\n";
    chdir $topFolderPath;

    my @casavaConfigs = @{ $testSetRef->{CASAVA} };
    if ( scalar(@casavaConfigs) == 1 ) {
        $casavaConfigPath =
          File::Spec->catdir( $testSetFolderPath, $runFolder,
            $casavaConfigs[0] );
    }
    else {
        die "ERROR: test framework supports one and only one casava config\n";
    }
    my $casavaScriptPath =
      File::Spec->catfile( $ENV_NAME_GA_DIR_PATH, $casavaScript );

    my $addCmdOption = $testRef->{addCmdOption};

    my $scriptCmd = sprintf( "%s %s $addCmdOption --EXPT_DIR %s --FORCE",
        $casavaScriptPath, $casavaConfigPath, $topFolderPath );

    #print "Running $casavaCmd\n";
    my ( $content, $error ) = getCmdOutput($scriptCmd);
    my $makeRunDir = "";

    if ( $error =~ /OUT_DIR\s+(\S+)\n/ ) {
        $makeRunDir = $1;
        chomp($makeRunDir);
    };
    
    my $contentMake = ""; 
    my $errorMake = "";

    if ( $testRef->{runMake} ne "N" ) {
        my $isRecursive = ( $testRef->{runMake} eq "R" ) ? 1 : 0;
        ( $contentMake, $errorMake ) =
          runMake( $makeRunDir, 2 * getProcessorsCount(), $isRecursive )
          ;    
    }

    my ( $status, $msg ) = checkResults(
        %$testsConfRef, $testId, $testSetId, $type, $makeRunDir, $scriptCmd,
        $content . $contentMake,
        $error . $errorMake
    );

    if ( $status != 0 ) {
        my $timeTmp  = strftime "%d_%m_%y-%H_%M_%S", localtime;
        my $webadd   = "";
        my $errorDir = "$testId-$testSetId-err-$timeTmp";
        my $temp     =
          File::Spec->catdir( $ENV_NAME_GA_DIR_PATH, $tempSpace, $errorDir );
        my $errorLog = File::Spec->catdir( $topFolderPath, "$errorDir.txt" );
        chdir File::Spec->catdir($currentDirectory);
        write2file( $errorLog,
            "\n\n    -------- Messages --------\n\n" . $msg . "\n\n    -------- Content --------\n\n" . $content . 
            "\n\n    -------- Content Make --------\n\n" . $contentMake . "\n\n    -------- Errors --------\n\n" . $error . 
            "\n\n    -------- Errors Make --------\n\n" . $errorMake );
        my $httpAddr = $testsConfRef->{conf}[0]->{WEB_DIR_ROOT}[0]->{address};
        my $httpPath = $errorLog;
        $httpPath =~ s/data//;
#           print "</pre><A href=\"$httpAddr$httpPath\">$errorDir.txt</A><pre>\n";
        print "$httpAddr$httpPath\n";

        $msg =
"test $testId failed for $testSetId see for details $errorLog\n";
    }

    if ( $testRef->{tearDown} eq "Y" ) {
        if ( $makeRunDir ne "" ) {
            executeCmd( "rm -fr $makeRunDir", 1 );
        }
    }

    chdir File::Spec->catdir($currentDirectory);

    return $status, $msg;
}    # runCasava

1;   # says use was ok
__END__

