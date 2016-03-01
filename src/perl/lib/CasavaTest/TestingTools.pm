package CasavaTest::TestingTools;

BEGIN {
    use Exporter();
    @ISA    = qw(Exporter);
    @EXPORT = qw(&copyImageAnaFiles &copyImages &executeCmd
      &executeCmds &copyFirecrestAll
      &write2file $IPAR_DIR_NAME &getCmdOutput &reorderRunFolders
      &getProcessorsCount &readFile &copyIpar &parseGeraldPath &checkResults
      &getGeraldFolder &getCasavaGeraldFolder &initialiseTestRefStructure 
      &initialiseTestSetRefStructure &runMake $tempSpace);
    @EXPORT_OK = qw();
}

# PROJECT: GERALD
# MODULE:  TestingTools.pm
# AUTHOR:  L. Szajkowski
#
# Copyright (c)  2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.

#

=pod

=head1 NAME

CasavaTest::Validate.pm - Functions used to test the pipeline

=head2 SYNOPSIS

use CasavaTest::TestingTools.pm qw();  

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

=head2 Overviewmv

Messaging Functions 

=head2 Exports

Global variables:

=head2 Depends on

    warnings, strict, POSIX, XML::Simple, constant, Carp

=cut

use warnings FATAL => 'all';
use strict;
use POSIX;
use XML::Simple;
use constant DEBUG => 0;    # set to 1 to get debug info
use Data::Dumper;
use Carp;
use File::Copy;
use IO::File;
sub copyImageAnaFiles($;$;$;\%);
sub copyImages($;$;$;$;\%);
sub executeCmd ($;$);
sub write2file($;$);
sub getCmdOutput($;$);
sub readFile($);
sub reorderRunFolders(\@);
sub getProcessorsCount ();
sub copyIpar($;$;$;$;\%);
sub parseGeraldPath($);
sub checkResults(\%;$;$;$;$;$;$;$);
sub getGeraldFolder($;\%;\%;);
sub getCasavaGeraldFolder($;$;$);
sub initialiseTestRefStructure($;$;$;$);
sub initialiseTestSetRefStructure($;$;$;$);
sub runMake ($;$;$);
sub copyFirecrestAll(\%;$;$;$;$;$);
sub executeCmds (\@);
sub xmlReConvert($);

my @chanes = ( 'a', 'c', 'g', 't' );
our $IPAR_DIR_NAME = "IPAR_1.01";
our $tempSpace     = File::Spec->catdir( "build", "Testing", "Temporary" );

=pod

=head1 FUNCTIONS

=head2 General Functions

=over 4

=cut

=pod

=item copyImages($longRunFolderName, $readPairMode)

The procedure copy images (tif files) from run folder Image directory
to destination directory. 

B<Parameters:>

    $longRunFolderName   - long name of run folder
    $readPairMode        - [0,1 or 2] Optional readPairMode appends '-1' or '-2' to readID

B<Returns:> 

    Nothing

=cut

sub copyImages($;$;$;$;\%) {
    croak("ERROR: copyImages wrong number of parameters")
      unless ( @_ == 5 );
    my ( $runDir, $lanesStr, $cyclesStr, $tilesStr, $dirsRef ) = @_;
    my @selectedLanes  = split ",", $lanesStr;
    my @selectedTiles  = split ",", $tilesStr;
    my @selectedCycles = split ",", $cyclesStr;

    foreach my $lane (@selectedLanes) {
        my $lanePath = File::Spec->catdir( $runDir, "Images", "L00$lane" );
        if ( $cyclesStr eq 'auto' ) {
            opendir( RUNFOLDER, $lanePath )
              or croak "couldn't open $lanePath: $!";
            @selectedCycles = sort grep !/^\./, readdir(RUNFOLDER);
            closedir RUNFOLDER;
        }
        my $outDirPath =
          File::Spec->catdir( $dirsRef->{imagesPath}, "L00$lane" );
        if ( !-d $outDirPath ) {
            mkdir "$outDirPath";
            foreach my $cycle (@selectedCycles) {

                my $cyclePath = File::Spec->catdir( $lanePath, "$cycle" );
                my $outCycleDirPath =
                  File::Spec->catdir( $dirsRef->{imagesPath},
                    "L00$lane", "$cycle" );

                if ( !-d $outCycleDirPath ) {
                    mkdir "$outCycleDirPath";
                }
                my $fileList = "";
                foreach my $tile (@selectedTiles) {
                    foreach my $chanel (@chanes) {
                        my $tilePath =
                          File::Spec->catfile( $cyclePath,
                            "s_$lane\_$tile\_$chanel.tif" );
                        if ( !-e "$outCycleDirPath/$tilePath" ) {
                            $fileList .= "$tilePath ";
                        }
                    }
                }

                if ( $fileList ne "" ) {
                    my $cmd = "cp $fileList $outCycleDirPath";
                    if ( system($cmd) != 0 ) {
                        return -1;
                    }
                }
            }
        }
    }
}    # sub copyImages

=pod

=item copyImageAnaFiles($runDir, $isIpar, $outDir)

The procedure copy images (tif files) from run folder Image directory
to destination directory. 

B<Parameters:>

    $longRunFolderName   - long name of run folder
    $isIpar              - is it it ipar run
    $readPairMode        - [0,1 or 2] Optional readPairMode appends '-1' or '-2' to readID
    $dirsRef             - HASH MAP Ref where directory configuration will be stored

B<Returns:> 

    Nothing

=cut

sub copyImageAnaFiles($;$;$;\%) {
    croak("ERROR: copyImageAnaFiles wrong number of parameters")
      unless ( @_ == 4 );
    my ( $runDir, $outDir, $isIpar, $dirsRef ) = @_;
    my $path          = File::Spec->catdir($runDir);
    my @laneList      = ();
    my @tmp           = split "/", $runDir;
    my $runFolderName = $tmp[-1];
    print "Processing $runFolderName\n";
    $dirsRef->{ELAND_test} = File::Spec->catdir( $outDir, "ELAND_test" );
    $dirsRef->{outPath} = File::Spec->catdir( $outDir, $runFolderName );
    $dirsRef->{dataPath} = File::Spec->catdir( $dirsRef->{outPath}, "Data" );
    $dirsRef->{imagesPath} =
      File::Spec->catdir( $dirsRef->{outPath}, "Images" );
    $dirsRef->{iparPath} =
      File::Spec->catdir( $dirsRef->{dataPath}, $IPAR_DIR_NAME );
    $dirsRef->{iparFirecresttPath} =
      File::Spec->catdir( $dirsRef->{iparPath}, "Firecrest" );

    #    print "Creating $outDir\n";
    if ( !-d $dirsRef->{outPath} )    { mkdir "$dirsRef->{outPath}"; }
    if ( !-d $dirsRef->{dataPath} )   { mkdir "$dirsRef->{dataPath}"; }
    if ( !-d $dirsRef->{imagesPath} ) { mkdir "$dirsRef->{imagesPath}"; }

    if ( !-d $dirsRef->{ELAND_test} ) {
        my $srcFile = File::Spec->catfile( $runDir, "..", "ELAND_test" );
        executeCmd( "cp -r $srcFile $outDir", 0 );
    }

    my %files = (
        config        => File::Spec->catfile( ".", "config.txt" ),
        PairedEndInfo => File::Spec->catfile( ".", "PairedEndInfo.xml" ),
        params1       => File::Spec->catfile( ".", "$runFolderName.params" ),
        params2       => File::Spec->catfile( ".", "Data", ".params" ),
        default_offsets =>
          File::Spec->catfile( ".", "Data", "default_offsets.txt" )
    );

    foreach my $file ( keys %files ) {
        my $srcFile = File::Spec->catfile( $runDir, $files{$file} );
        my $dscFile = File::Spec->catfile( $dirsRef->{outPath}, $files{$file} );

        #print "copy $srcFile, $dscFile\n";
        copy( $srcFile, $dscFile );
    }

    my $iparDir = File::Spec->catdir( $runDir, "Data", $IPAR_DIR_NAME );

    if ( -d $iparDir ) {
        if ( !-d $dirsRef->{iparPath} ) { mkdir "$dirsRef->{iparPath}"; }
        if ( !-d $dirsRef->{iparFirecresttPath} ) {
            mkdir "$dirsRef->{iparFirecresttPath}";
        }
        executeCmd( "rm -f $dirsRef->{iparPath}/Makefile*", 0 );
    }
}    # sub copyGEFirecrestFiles

=pod

=item copyImages($longRunFolderName, $readPairMode)

The procedure copy images (tif files) from run folder Image directory
to destination directory. 

B<Parameters:>

    $longRunFolderName   - long name of run folder
    $readPairMode        - [0,1 or 2] Optional readPairMode appends '-1' or '-2' to readID

B<Returns:> 

    Nothing

=cut

sub copyIpar($;$;$;$;\%) {
    croak("ERROR: copyImages wrong number of parameters")
      unless ( @_ == 5 );
    my ( $runDir, $lanesStr, $cyclesStr, $tilesStr, $dirsRef ) = @_;
    my @selectedLanes  = split ",", $lanesStr;
    my @selectedTiles  = split ",", $tilesStr;
    my @selectedCycles = split ",", $cyclesStr;

    foreach my $lane (@selectedLanes) {
        my $lanePath1 =
          File::Spec->catdir( $runDir, "Data", "$IPAR_DIR_NAME", "Firecrest",
            "L00$lane" );
        my $lanePath2 = File::Spec->catdir( $runDir, "Data", "$IPAR_DIR_NAME" );
        if ( $cyclesStr eq 'auto' ) {
            opendir( RUNFOLDER, $lanePath1 )
              or croak "couldn't open $lanePath1: $!";
            @selectedCycles = sort grep !/^\./, readdir(RUNFOLDER);
            closedir RUNFOLDER;
        }
        my $outDirPath =
          File::Spec->catdir( $dirsRef->{iparFirecresttPath}, "L00$lane" );
        if ( !-d $outDirPath ) {
            mkdir "$outDirPath";
            foreach my $tile (@selectedTiles) {
                my $fileList = "";
                my $tileStr  = sprintf "%04d", $tile;
                my $tilePath =
                  File::Spec->catfile( $lanePath1,
                    "s_$lane\_$tileStr\_*_qcm.xml" );
                $fileList .= "$tilePath ";
                $fileList .=
                  File::Spec->catfile( $lanePath1,
                    "s_$lane\_$tileStr\_pos.txt" )
                  . " ";
                my $cmd = "cp $fileList $outDirPath";
                system $cmd;
                $fileList = "";
                $fileList .=
                  File::Spec->catfile( $lanePath2,
                    "s_$lane\_$tileStr\_int.txt.p.gz" )
                  . " ";
                $fileList .=
                  File::Spec->catfile( $lanePath2,
                    "s_$lane\_$tileStr\_nse.txt.p.gz" )
                  . " ";
                $cmd = "cp $fileList $dirsRef->{iparPath}";

                #print $cmd . "\n";
                system $cmd;
            }
        }
    }
}    # sub copyImages

=pod

=item copyImages($longRunFolderName, $readPairMode)

The procedure copy all Firecrest/Bustard/Gerald data from data set folder to
    build/Testing folder 

B<Parameters:>

    $testsConfRef          - HASH MAP Ref to struture with pipeline_tests.xml
    $GA_TEST_DATA_PATH     - Path to GA_TEST_DATA folder (from command line)
    $ENV_NAME_GA_DIR_PATH  - Path to pipeline installation
    $testId                - Id of test in pipeline_tests.xml
    $testSetId             - Id of test set in pipeline_tests.xml
    $type                  - Type of test [IPAR_1.01 | GA]

B<Returns:> 

    0 - no need for gerald folder
    2 if Firecrest folder already exist; 1 - when no Firecrest folder

=cut

sub copyFirecrestAll(\%;$;$;$;$;$) {
    croak("ERROR: copyFirecrestAll wrong number of parameters")
      unless ( @_ == 6 );
    my ( $testsConfRef, $GA_TEST_DATA_PATH, $ENV_NAME_GA_DIR_PATH, $testId,
        $testSetId, $type )
      = @_;

    my $testRef         = $testsConfRef->{test}->{$testId};
    my $testSetRef      = $testsConfRef->{set}->{$testSetId};
    my $geraldFolderRef = getGeraldFolder( $type, %$testRef, %$testSetRef );
    if ( defined $geraldFolderRef ) {
        my $dirsRef       = parseGeraldPath( $geraldFolderRef->{path} );
        my $runFolder     = $dirsRef->{runFolder};
        my $firecrestPath = $dirsRef->{FirecrestFull};

        my $srcRunFolder =
          File::Spec->catdir( $GA_TEST_DATA_PATH, $testSetId, $firecrestPath );
        my $dscRunFolder = File::Spec->catdir( $ENV_NAME_GA_DIR_PATH,
            $tempSpace, $testSetId, $runFolder, "Data" );
        my $dscFirectestPath =
          File::Spec->catdir( $ENV_NAME_GA_DIR_PATH, $tempSpace, $testSetId,
            $firecrestPath );

        if ( !-d $dscFirectestPath ) {
            executeCmd( "cp -r $srcRunFolder $dscRunFolder", 0 );
            return 1;
        }
        return 2;
    }
    return 0;
}    # sub copyFirecrestAll

=pod

=head1 The procedure executes the command.

=over 4

=item executeCmd($command, $verbose)

The procedure executes the command but exits with error code when 
command failed.

Parameters:
    command - command to be executed
    verbose - verbose level (optional)
Returns
    Nothing 
=back

=cut

sub executeCmd ($;$) {
    croak "ERROR: executeCmd " . join( ',', @_ ) unless ( @_ == 1 || @_ == 2 );
    my ( $command, $verbose ) = @_;
    print $command . "\n" if ( defined $verbose && $verbose > 0 );
    system($command);
    my $exit_value  = $? >> 8;
    my $signal_num  = $? & 127;
    my $dumped_core = $? & 128;
    if ( $signal_num > 0 ) {
        croak
          "ERROR: executeCmd() Process killed by SIG($signal_num): $command\n";
    }    # if
    if ( $dumped_core > 0 ) {
        croak "ERROR: executeCmd() $command dumped core: $command \n";
    }    # if
    if ( $exit_value != 0 ) {
        croak "ERROR: executeCmd() $0 EXIT_CODE($exit_value): $command  $!\n";
    }    # if
}

=pod

=head1 The procedure executes an arrayr of commans command.

=over 4

=item executeCmd($command, $verbose)

The procedure executes and array of commands but exits with error code when 
command failed.

Parameters:
    $commandsRef - command to be executed
    verbose      - verbose level (optional)
Returns
    Nothing 
=back

=cut

sub executeCmds (\@) {
    die "ERROR: executeCmds " . join( ',', @_ ) unless ( @_ == 1 || @_ == 2 );
    my ( $commandsRef, $verbose ) = @_;
    foreach my $cmdRef (@$commandsRef) {
        executeCmd( $cmdRef->{cmd}, $verbose );
    }
}

=pod

=head1 The procedure executes the command.

=over 4

=item runMake($command, $verbose)

The procedure runs make in selected directory

Parameters:
    $runDir    - directory where make will be run
    $jobNumber - number of jobs (processors)
    $recursive - if 1 then make recursive
Returns
    Commnand output, command errors
=back

=cut

sub runMake ($;$;$) {
    die "ERROR: runMake " . join( ',', @_ ) unless ( @_ == 3 );
    my ( $runDir, $jobNumber, $recursive ) = @_;
    if ( $runDir eq "" ) {
        return "", "runMake empty folder\n";
    }
    my $curCir = Cwd->getcwd();
    chdir $runDir;

    my $makeCmd = "make ";
    if ( $jobNumber > 1 ) {
        $makeCmd = "$makeCmd -j $jobNumber";
    }

    if ( $recursive > 0 ) {
        $makeCmd = "$makeCmd recursive";
    }

    my ( $content, $error ) = getCmdOutput( $makeCmd, 0 );
    chdir File::Spec->catdir($curCir);

    return $content, $error;
}

=pod

=head1 The procedure executes the command.

=over 4

=item getCmdOutput($command, $verbose)

The procedure executes the command and retunrs commands output.

Parameters:
    command - command to be executed
    verbose - verbose level (optional)
Returns
    Commnand output, command errors
=back

=cut

sub getCmdOutput ($;$) {
    die "ERROR: getCmdOutput " . join( ',', @_ ) unless ( @_ == 1 || @_ == 2 );
    my ( $command, $verbose ) = @_;
    print $command . "\n" if ( defined $verbose && $verbose > 0 );
    my $tmpFile    = "/tmp/getCmdOutput-$$.txt";
    my $commandTmp = $command . " 2> $tmpFile";

    my $content      = "";
    my $errorContent = "";
    my $count        = 0;
    open( CMD, "$commandTmp |" )
      || die "Could not $commandTmp $! \n";

    while (<CMD>) {
        $content .= $_;
    }
    close(CMD);
    if ( -e $tmpFile ) {
        $errorContent = readFile($tmpFile);
    }
    unlink $tmpFile;
    return $content, $errorContent;
}

=pod

=item write2file($filePath, $content)

The procedure saves the content to a file

B<Parameters:>

    $filePath   - full path to file
    $content    - content to be save to a file

B<Returns:> 

    Nothing

=cut

sub write2file($;$) {
    croak("ERROR: write2file wrong number of parameters")
      unless ( @_ == 2 );
    my ( $filePath, $content ) = @_;
    my $file = IO::File->new( ">" . $filePath )
      || croak "ERROR: Couldn't create/open file handle for $filePath $!\n";
    print $file $content;
    close $file;
}

=pod

=item readFile($filePath)

The procedure reads the content from a file

B<Parameters:>

    $filePath   - full path to file

B<Returns:> 

    File content

=cut

sub readFile($) {
    croak("ERROR: readFile wrong number of parameters")
      unless ( @_ == 1 );
    my ($filePath) = @_;
    my $content    = "";
    my $file       = IO::File->new($filePath)
      || die "ERROR: Couldn't create/open file handle for $filePath $!\n";
    while (<$file>) {
        $content .= $_;
    }
    close $file;
    return $content;
}

=pod

=item parse Gerald path

The procedure parses gerald path into Firecrest, Bustard and run folder path

B<Parameters:>

    $geraldPath   - full path to gerald folder

B<Returns:> 

    HASH MAP Ref to parsed data

=cut

sub parseGeraldPath($) {
    croak("ERROR: parseGeraldPath wrong number of parameters")
      unless ( @_ == 1 );
    my ($geraldPath) = @_;

    my @dirs = split( "/", $geraldPath );
    if ( scalar(@dirs) < 5 ) {
        die "ERROR: parseGeraldPath [$geraldPath] at least full path "
          . "from run folder to gerald should be provided\n";
    }
    my %runDirs = ();
    $runDirs{GERALD}        = $dirs[-1];
    $runDirs{Bustard}       = $dirs[-2];
    $runDirs{Firecrest}     = $dirs[-3];
    $runDirs{Data}          = $dirs[-4];
    $runDirs{runFolder}     = $dirs[-5];
    $runDirs{runFolderFull} = $dirs[-5];
    $runDirs{DataFull}      = $runDirs{runFolderFull} . "/" . $runDirs{Data};
    $runDirs{FirecrestFull} = $runDirs{DataFull} . "/" . $runDirs{Firecrest};
    $runDirs{BustardFull}   = $runDirs{FirecrestFull} . "/" . $runDirs{Bustard};
    $runDirs{GERALDFull}    = $geraldPath;
    return \%runDirs;
}    # sub parseGeraldPath

=pod

=item reorderRunFolders($runFoldesRef)

The procedure creates arraye with run folder name(s) 
where for PFPE R1 is always before R2

B<Parameters:>

    $runFoldesRef  - ARRAY REF to array of folder names

B<Returns:> 

    folder names string ( folderName or folderNameR1 folderNameR2)

=cut

sub reorderRunFolders(\@) {
    croak("ERROR: reorderRunFolders wrong number of parameters")
      unless ( @_ == 1 );
    my ($runFoldesRef) = @_;
    my @folders = ();
    if ( scalar( @{$runFoldesRef} ) == 1 ) {
        push @folders, $runFoldesRef->[0];
    }
    elsif ( scalar( @{$runFoldesRef} ) == 2 ) {
        my $runFolderR1 = "";
        my $runFolderR2 = "";
        if (   $runFoldesRef->[0] =~ /.+R1$/
            && $runFoldesRef->[1] =~ /.+R2$/ )
        {
            push @folders, $runFoldesRef->[0];
            push @folders, $runFoldesRef->[1];
        }
        elsif ($runFoldesRef->[1] =~ /.+R1$/
            && $runFoldesRef->[0] =~ /.+R2$/ )
        {
            push @folders, $runFoldesRef->[1];
            push @folders, $runFoldesRef->[0];
        }
        else {
            croak "ERROR: $0 run folder configuration "
              . join( ":", @{$runFoldesRef} ) . "\n";
        }
    }
    else {
        croak "ERROR: $0 unsupported number of run folders "
          . scalar( @{$runFoldesRef} ) . "\n";
    }
    return @folders;
}    # sub joinRunFolders

=pod

=item getProcessorsCount()

Procedure to read the number of processors 

B<Returns:> 

    number of processors

=cut

sub getProcessorsCount () {
    die "ERROR: getProcessorsCount wrong parameters \n"
      unless ( @_ == 0 );
    my $command = "cat /proc/cpuinfo | grep processor | wc -l 2> /dev/null";
    my $count   = 0;
    open( SGE, "$command |" )
      || die
"ERROR: getProcessorsCount() Couldn't run cat /proc/cpuinfo | grep processor | wc $! \n";
    while (<SGE>) {
        chomp($_);
        $count = $_;
    }    # while
    return $count;
}

=pod

=item checkResults(testsConfRef, $testId,  $testSetId, $type,    
        $currentDir,   $command, $content,   $error)

The procedure runs checks expectedOut and expectedError rules
 on the command output and runs testCmd in current directory
runs 

B<Parameters:>

    $testsConfRef          - HASH MAP Ref to struture with pipeline_tests.xml
    $testId                - Id of test in pipeline_tests.xml
    $testSetId             - Id of test set in pipeline_tests.xml
    $type                  - Type of test [IPAR_1.01 | GA]
    $currentDir            - directory where testCmd will be run
    $command               - command which created results
    $content               - std output of the command
    $error                 - error output of the command

B<Returns:> 

    status, message
    
=cut

sub checkResults(\%;$;$;$;$;$;$;$) {
    croak("ERROR: checkResults wrong number of parameters")
      unless ( @_ == 8 );
    my (
        $testsConfRef, $testId,  $testSetId, $type,
        $currentDir,   $command, $content,   $error
      )
      = @_;

    my $msg           = "";
    my $status        = 0;
    my @expectedOut   = ();
    my @expectedError = ();

    my $testRef    = $testsConfRef->{test}->{$testId};
    my $testSetRef = $testsConfRef->{set}->{$testSetId};

    my $curCir = Cwd->getcwd();
    if ( $currentDir ne "" ) {

        #print "checkResults - Changing folder to $currentDir\n";
        chdir $currentDir;
    }
    if ( defined $testRef->{expectedOut} ) {
        @expectedOut = @{ $testRef->{expectedOut} };
    }
    if ( defined $testRef->{expectedError} ) {
        @expectedError = @{ $testRef->{expectedError} };
    }

    $msg .= "Running $command\n";

    #print "$msg";
    foreach my $ruleRef (@expectedOut) {
        my $regEx = $ruleRef->{out};
        $ruleRef->{out} = xmlReConvert($ruleRef->{out});
        next
          if ( $ruleRef->{testSet} ne "all"
            && $ruleRef->{testSet} ne $testSetId );

        if ( $content =~ /$regEx/ ) {

        }
        else {
            $status = -1;
            $msg .=
              "ERROR: \"$regEx\" missing in test [$testId] output\n$content";
        }
    }
    foreach my $ruleRef (@expectedError) {
        $ruleRef->{out} = xmlReConvert($ruleRef->{out});            
        my $regEx = $ruleRef->{out};
        next
          if ( $ruleRef->{testSet} ne "all"
            && $ruleRef->{testSet} ne $testSetId );
        if ( $error !~ /$regEx/ ) {
            $status = -1;
            $msg .=
"ERROR: \"$regEx\" missing in test [$testId] error output\n$error";
        }
    }

    if ( defined $testRef->{testCmd} ) {
        my @testCmds = @{ $testRef->{testCmd} };
        foreach my $cmd (@testCmds) {
            next
              if ( $cmd->{testSet} ne "all"
                && $cmd->{testSet} ne $testSetId );
            $cmd->{cmd} = xmlReConvert($cmd->{cmd});
            my ( $content, $error ) = getCmdOutput( $cmd->{cmd} );
            chomp($content);
            $cmd->{out} = xmlReConvert($cmd->{out});            
            
            if ( $content !~ /$cmd->{out}/ && length($error) == 0 ) {
                $status = -1;
                $msg .=
"ERROR: test $testId cmd $cmd->{cmd} was expecting [$cmd->{out}] but it got [$content/$error]\n";
            }
        }
    }

    if ( defined $testRef->{xpath} ) {
        my @xpathCmds = @{ $testRef->{xpath} };
        my ( $statusTmp, $msgTmp ) =
          checkXpath( $testsConfRef, $testId, $testSetId );
        if ( defined $statusTmp ) {
            $status = $statusTmp;
        }
        if ( defined $statusTmp ) {
            $msg = $msgTmp;
        }
    }
    chdir $curCir;

    return $status, $msg;
}    # checkResults

=pod

=item getGeraldFolder($type, $testRef, $testSetRef)

The procedure runs xpath query on a xml file

B<Parameters:>

$testsConfRef          - HASH MAP Ref to struture with pipeline_tests.xml
    $testId                - Id of test in pipeline_tests.xml
    $testSetId             - Id of test set in pipeline_tests.xml
B<Returns:> 

    status and mesaage 
    
=cut

sub checkXpath(\%;$;$) {
    croak("ERROR: getGeraldFolder wrong number of parameters")
      unless ( @_ == 3 );
    my ( $testsConfRef, $testId, $testSetId ) = @_;
    my ( $status, $msg );
    my $testRef = $testsConfRef->{test}->{$testId};

    my @xpathCmds = @{ $testRef->{xpath} };
    foreach my $xpathRef (@xpathCmds) {
        next
          if ( $xpathRef->{testSet} ne "all"
            && $xpathRef->{testSet} ne $testSetId );
        $xpathRef->{query} = xmlReConvert($xpathRef->{query});
        my $xmlFilePath = $xpathRef->{file};

        if ( !-e $xmlFilePath ) {
            $status = -1;
            $msg .=
              "ERROR: test $testId file " . "$xmlFilePath cannot be found \n";
            next;
        }

        my $cmd = "xpath $xmlFilePath '$xpathRef->{query}'";
        my ( $content, $error ) = getCmdOutput($cmd);
        chomp($content);

        if (   $xpathRef->{count} ne ""
            && $error !~ /Found $xpathRef->{count} nodes:/ )
        {
            $status = -1;
            $msg .=
                "ERROR: test $testId xpath query $xpathRef->{query} was "
              . "expecting [Found $xpathRef->{count} nodes:] but it got [$error]\n";
        }
        if ( $xpathRef->{out} ne "" && $content !~ /$xpathRef->{out}/ ) {
            $status = -1;
            $msg .=
                "ERROR: test $testId xpath query $xpathRef->{query} was "
              . "expecting [$xpathRef->{out}] but it got [$content]\n";
        }
    }
    return $status, $msg;
}    # getGeraldFolder

=pod

=item getGeraldFolder($type, $testRef, $testSetRef)

The procedure selectes gerald folder from testSet configuration based on testId

B<Parameters:>

    $type                  - Type of test [IPAR_1.01 | GA]

B<Returns:> 

    gerald folder path or empty string when the test has 
    no valida gerald path 
    
=cut

sub getGeraldFolder($;\%;\%;) {
    croak("ERROR: getGeraldFolder wrong number of parameters")
      unless ( @_ == 3 );
    my ( $type, $testRef, $testSetRef ) = @_;

    my $geraldFolderRef   = undef;
    my $testCompatibility = "";
    if ( defined $testRef->{compatibility} ) {
        $testCompatibility = $testRef->{compatibility};
    }

    my @GERALDFolders = @{ $testSetRef->{GERALDFolder} };
    foreach my $folderRef (@GERALDFolders) {

        #        print "$testCompatibility $testRef->{compatibility} $type\n";
        if (   defined $folderRef->{compatibility}
            && $folderRef->{compatibility} eq $testCompatibility
            && $folderRef->{type}          eq $type )
        {
            $geraldFolderRef = $folderRef;
        }
    }

    return $geraldFolderRef;
}    # getGeraldFolder

=pod

=item

=cut

sub getCasavaGeraldFolder($;$;$)
{
    croak("ERROR: getCasavaFolder wrong number of parameters")
        unless(@_ == 3);
    my($testParameterData, $suiteRef, $type) = @_;

    my $casavaFolderRef = undef;
    my $testCompatability = "";

    if(defined $testParameterData->{$suiteRef."-compatibility"})
    {
        $testCompatability = $suiteRef->{compatibility};
        print "Test compatibility - $testCompatability\n";
    }

    my $geraldFolders = $testParameterData->{"Suite-BF00151_ecoli_SE_verification_suite-GERALDFolder"};
    return $geraldFolders;
}    # getCasavaGeraldFolder

=pod

=item initialiseTestRefStructure($testRef, $testSpecificationPrefix, $testSpecificationData, $testParameterData)

This procedure is used by the CASAVA test code to fill a testRef structure
used by the testing framework for Firecrest/Bustard/GERALD.

TODO
Eventually refactor test framework to use CASAVA procedure

=cut

sub initialiseTestRefStructure($;$;$;$)
{
    croak("ERROR: initialiseTestRefStructure wrong number of parameters")
        unless(@_ == 4);
    my($testRef, $testSpecificationPrefix, $testConfigurationData, $testParameterData) = @_;
    $testRef->{tag} = $testParameterData->{$testSpecificationPrefix."tag"};
    $testRef->{set} = $testParameterData->{$testSpecificationPrefix."set"};
    $testRef->{cycles} = $testParameterData->{$testSpecificationPrefix."cycles"};
    $testRef->{lanes} = $testParameterData->{$testSpecificationPrefix."lanes"};
    $testRef->{tiles} = $testParameterData->{$testSpecificationPrefix."tiles"};
    $testRef->{script} = $testParameterData->{$testSpecificationPrefix."script"};
    $testRef->{enabled} = $testParameterData->{$testSpecificationPrefix."enabled"};
#    $testRef->{type} = $testParameterData->{$testSpecificationPrefix."type"};
    $testRef->{xpath} = $testParameterData->{$testSpecificationPrefix."xpath"};
#    $testRef->{xpath}->{testSet} = $testParameterData->{$testSpecificationPrefix."testSet"};
#    $testRef->{xpath}->{query} = $testParameterData->{$testSpecificationPrefix."xpathQuery"};
#    $testRef->{xpath}->{out} = $testParameterData->{$testSpecificationPrefix."xpathOut"};
#    $testRef->{xpath}->{count} = $testParameterData->{$testSpecificationPrefix."xpathCount"};
    $testRef->{expectedOutput} = $testParameterData->{$testSpecificationPrefix."expectedOutput"};
#    $testRef->{expectedOutput}->{testSet} = $testParameterData->{$testSpecificationPrefix."expectedOutputTestSet"};
#    $testRef->{expectedOutput}->{out} = $testParameterData->{$testSpecificationPrefix."expectedOutputOut"};
    $testRef->{expectedError} = $testParameterData->{$testSpecificationPrefix."expectedError"};
#    $testRef->{expectedError}->{testSet} = $testParameterData->{$testSpecificationPrefix."expectedErrorTestSet"};
#    $testRef->{expectedError}->{out} = $testParameterData->{$testSpecificationPrefix."expectedErrorOut"};
    $testRef->{runMake} = $testParameterData->{$testSpecificationPrefix."runMake"};
    $testRef->{script} = $testParameterData->{$testSpecificationPrefix."script"};
    $testRef->{setUp} = $testParameterData->{$testSpecificationPrefix."setUp"};
    $testRef->{setUpPostCommand} = $testParameterData->{$testSpecificationPrefix."setUpPostCommand"};
    $testRef->{tearDown} = $testParameterData->{$testSpecificationPrefix."tearDown"};
    $testRef->{tearDownPostCommand} = $testParameterData->{$testSpecificationPrefix."tearDownPostCommand"};
    $testRef->{testCommand} = $testParameterData->{$testSpecificationPrefix."test_command"};
#    $testRef->{testCommand}->{cmd} = $testParameterData->{$testSpecificationPrefix."testCommandCmd"};
#    $testRef->{testCommand}->{testSet} = $testParameterData->{$testSpecificationPrefix."testCommandTestSet"};
#    $testRef->{testCommand}->{out} = $testParameterData->{$testSpecificationPrefix."testCommandOut"};
}

=pod

=item initialiseTestSetRefStructure($testSetRef, $testSpecificationPrefix, $testSpecificationData)

This procedure is used by the CASAVA test code to fill a testSetRef structure
used by the testing framework for Firecrest/Bustard/GERALD.

TODO
Eventually refactor test framework to use CASAVA procedure

=cut

sub initialiseTestSetRefStructure($;$;$;$)
{
    croak("ERROR: initialiseTestSetRefStructure wrong number of parameters")
        unless(@_ == 4);
    my($testSetRef, $campaignSpecificationPrefix, $suiteSpecificationPrefix, $testSpecificationData) = @_;

#    print "The campaign specification prefix supplied is $campaignSpecificationPrefix\nThe testSpecificationData contains -\n";
#    print Dumper($testSpecificationData);

    $testSetRef->{cycles} = $testSpecificationData->{$suiteSpecificationPrefix."cycles"};
    $testSetRef->{enabled} = $testSpecificationData->{$suiteSpecificationPrefix."enabled"};
    $testSetRef->{fasta} = $testSpecificationData->{$suiteSpecificationPrefix."fasta"};
    $testSetRef->{GERALD} = $testSpecificationData->{$suiteSpecificationPrefix."GERALD"};
    $testSetRef->{GERALDFolder} = $testSpecificationData->{$suiteSpecificationPrefix."GERALDFolder"};
#    $testSetRef->{GERALDFolder}->{compatibility} = $testSpecificationData->{$suiteSpecificationPrefix."compatibility"};
    $testSetRef->{readMode} = $testSpecificationData->{$suiteSpecificationPrefix."readMode"};
    $testSetRef->{runFolder} = $testSpecificationData->{$campaignSpecificationPrefix."runFolder"};
    $testSetRef->{type} = $testSpecificationData->{$suiteSpecificationPrefix."type"};
}

=pod

=item xmlReConvert($type, $testRef, $testSetRef)

The procedure converts XML equivalents to charaters 

    $string =~ s/&quot;/"/g;
    $string =~ s/&apos;/'/;
    $string =~ s/&lt;/</;
    $string =~ s/&gt;/>/;
    $string =~ s/&quot;/"/g;
    $string =~ s/&apos;/'/;
    $string =~ s/&lt;/</;
    $string =~ s/&gt;/>/; 
    $string =~ s/&gt;/&/g; 
B<Parameters:>

    $string                  - string to convert

B<Returns:> 

    converted string
    
=cut

sub xmlReConvert($) {
    croak("ERROR: xmlReConvert wrong number of parameters")
      unless ( @_ == 1 );
    my ($string) = @_;

    $string =~ s/&quot;/"/g;
    $string =~ s/&apos;/'/;
    $string =~ s/&lt;/</;
    $string =~ s/&gt;/>/;
    $string =~ s/&quot;/"/g;
    $string =~ s/&apos;/'/;
    $string =~ s/&lt;/</;
    $string =~ s/&gt;/>/;
    $string =~ s/&amp;/&/g;
    $string =~ s/&or/\|/g;
    $string =~ s/&gt;/>/g;
    return $string;    
}     # getGeraldFolder
1;    # says use was ok
__END__

