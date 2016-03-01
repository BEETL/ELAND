# =head1 LICENSE
#
# Copyright (c) 2007-2010 Illumina, Inc.
#
# This software is covered by the "Illumina Genome Analyzer Software
# License Agreement" and the "Illumina Source Code License Agreement",
# and certain third party copyright/licenses, and any user of this
# source file is bound by the terms therein (see accompanying files
# Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
# Illumina_Source_Code_License_Agreement.pdf and third party
# copyright/license notices).
#
# This file is part of the Consensus Assessment of Sequence And VAriation
# (CASAVA) software package.
#

=head1 SUMMARY

bam - aggregate all reads from a build into a single BAM file with
optional chromosome re-labeling.

=head1 DESCRIPTION

This plugin consolidates all reads from a CASAVA build into a single BAM
file, with optional chromosome renaming.

A set of auxillary files is created with each BAM file to facilitate use
in downstream packages such as SAMtools or the Broad IGV. These files are:

=over 4

=item 1. sorted.bam           - the bam file itself

=item 2. sorted.bam.bai       - index of the bam file

=item 3. sorted.bam.fa.gz     - gzipped fasta file containing the reference sequence(s) (this file may be optionally omitted)

=back

The output of the BAM plugin can be found in:

F<Project_Dir/genome/bam/>

B<IMPORTANT:> Note that in order for the I<bam> target to work, the
I<sort> target must have been completed or be part of the current
workflow.

=head1 REALIGNED READ HANDLING

If the optional realigned read output was selected in the
I<callSmallVariants> stage, the 'sorted.realigned.bam' files resulting
from this process will also be merged and renamed into a single genome
file in parallel with the 'sorted.bam' files. The result will be two
additional files:

F<Project_Dir/genome/bam/realigned/sorted.realigned.bam>

F<Project_Dir/genome/bam/realigned/sorted.realigned.bam.bai>

Note that these files only contain reads for which the alignment was
changed during variant calling, and that the original alignments still
exist in the 'sorted.bam' files.

=head1 OPTIONS

=over 4

=item --bamChangeChromLabels=OFF|NOFA|UCSC

Change chromosome labels in the bam plugin output. The available
behaviors are:

=over 8

=item OFF

Use unmodified CASAVA chromosome labels (default behavior).

=item NOFA

Remove any ".fa" suffix found on each chromosome label. For example
"c11.fa" is changed to "c11".

=item UCSC

Remove any ".fa" suffix found on each chromosome label and attempt to
map the result to the corresponding UCSC human chromosome label. For
example "c11.fa" is changed to "chr11".

=back

=item --bamSkipRefSeq

Do not generate a reference sequence file with each bam file. The
default behavior can be restored with "--no-bamSkipRefSeq".

=back

=cut


package Casava::PostAlignment::Plugins::Sorted2Bam;

use warnings FATAL => 'all';
use strict;

use File::Spec;
use Getopt::Long;
use Sys::Hostname;

use lib '@CASAVA_FULL_PERL_LIBDIR@';
use Casava::Common::Log;
use Casava::Common::IOLib qw(createDirs executeCmd);
use Casava::PostAlignment::Sequencing::Config qw(isSpliceJunctionChrom %CONF_APP);
use Casava::PostAlignment::Sequencing::SamLib qw(checkSamtoolsBin getChangeChrLabelType changeChrLabel);

use Casava::TaskManager qw(executeSingleTask checkPointInit
    executeArrayTaskEx executeCheckPoint addTask2checkPoint tryAddTask2CheckPoint);
use Casava::PostAlignment::Plugins::RefSeq qw(getDoneCheckPoint);

#------------------------------------------------------------------------------
#
#=pod
#
#=head2 getTarget()
#
#All plugins are required to implement this method. The framework will execute
#the C<executeOrGenerateWorkflow> if user has requested the target to be added
#to the build.
#
#=over 4
#
#=item *
#
#Parameters:
#
#  No parameters
#
#=item *
#
#Returns:
#
#  name of target for which the plugin can perform workflow configuration
#
#=item *
#
#Exceptions
#
#  None
#
#=back
#
#=cut

sub getTarget { return 'bam'; }

#------------------------------------------------------------------------------
#
#=pod
#
#=head2 getOptionsMapping
#
#Returns command-line options understood by plugin. Unless the option is one of
#the standard ones or supported by a plugin, the CASAVA command line will not
#accept it. All options supplied via command line are accessible through
#%CONF_PROJ hash
#
#I<Parameters:>
#
#=over 4
#
#=item *
#
#\%PARAMS
#
#reference to target hash
#
#=back
#
#I<Returns:>
#
#Mappings suitable for Getopt::Long::GetOptions
#
#I<Exceptions:>
#
#=over 4
#
#=item *
#
#None
#
#=back
#
#=cut

sub getOptionsMapping(\%) {
    my ( $PARAMS ) = @_;

     return ('bamSkipRefSeq!'  => \$PARAMS->{bamSkipRefSeq},
             'bamChangeChromLabels=s' => \$PARAMS->{bamChangeChromLabels});
}

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#
#=pod
#
#=head2 executeOrGenerateWorkflow
#
#All plugins are required to implement this method. The framework calls it when
#the target returned by L</getTarget> has been requested to be added to the
#CASAVA build.
#
#I<Parameters:>
#
#=over 4
#
#=item *
#
#$prevTarget
#
#if C<start>, the target is the first target in the workflow. Otherwise, plugin
#is expected to specify dependencies on the other target tasks (if needed).
#
#=item *
#
#$projectDir
#
#Folder containing chromosome folders.
#
#=item *
#
#$CONF_PROJ_Ref
#
#HASH MAP Ref to project configuration
#
#=item *
#
#$chromsBinSizesRef
#
#HASH MAP Ref to chromosome sizes
#
#=back
#
#I<Returns:>
#
#Nothing
#
#I<Exceptions:>
#
#=over 4
#
#=item *
#
#Incorrect analysis type - indels are possible only for paired DNA
#
#=back
#
#=cut

sub executeOrGenerateWorkflow($$\%\%)
{
    my ($prevTarget, $projectDir, $CONF_PROJ_Ref, $chromsBinSizesRef) = @_;

    my $libexecDir          = '@CASAVA_FULL_LIBEXECDIR@';
    my $buildDir            = $CONF_PROJ_Ref->{dirBuildParsed};
    my $scriptDir           = $libexecDir;
    my $bamDirName          = $CONF_APP{dirBam};

    my $isBamWholeGenome    = (defined $CONF_PROJ_Ref->{bamWholeGenome} ? $CONF_PROJ_Ref->{bamWholeGenome} : 0 );
    my $isBamSkipRefSeq     = (defined $CONF_PROJ_Ref->{bamSkipRefSeq} ? $CONF_PROJ_Ref->{bamSkipRefSeq} : 0 );

    #my $sorted2BamScript = File::Spec->catfile($scriptDir, 'sorted2bam.pl');

    # configuration-time check for samtools binary and logging the path of which version was found:
    checkSamtoolsBin(%{$CONF_PROJ_Ref});

    # configuration-time check that a valid change-chrom code was provided:
    my $bamChangeChrType = getChangeChrLabelType($CONF_PROJ_Ref->{bamChangeChromLabels});

    my @chroms = keys %{$chromsBinSizesRef};

    {
        # additional check that chromosome name remapping will be valid:
        my %testLabels;
        for my $chrom (@chroms) {
            my $newLabel = changeChrLabel($chrom,$bamChangeChrType);
            if(exists $testLabels{$newLabel}){
                my $firstChrom=$testLabels{$newLabel};
                my $errMsg="ERROR: Current bam target options cannot be used with the chromosome labels in this build because the CASAVA chromosome labels: '$firstChrom' and '$chrom' would both be changed to the same BAM chromosome label: '$newLabel'.";
                errorExit($errMsg);
            }
            $testLabels{$newLabel} = $chrom;
        }
    }

    # setup output directory:
    my $genomeDirName = 'genome';

    my @dirs = ();
    push @dirs, $genomeDirName;
    push @dirs, File::Spec->catdir($genomeDirName,$bamDirName);
    createDirs($projectDir, @dirs);

    my $taskIdRefSeq = Casava::PostAlignment::Plugins::RefSeq::getDoneCheckPoint();

    my $allChrTaskLabel      = "genomeBam";
    my $allChrCheckPointName = "GENOMEBAM_PREREQ";
    my $allChrCheckPointRef  = checkPointInit( $allChrCheckPointName,
                                               $allChrTaskLabel );

    tryAddTask2CheckPoint( "SYNCH_MERGE_DONE", %{$allChrCheckPointRef} );

    executeCheckPoint( %$allChrCheckPointRef );

    my $genomeBamScript = File::Spec->catfile($scriptDir, 'makeGenomeBam.pl');
    my $genomeBamCmd = "$genomeBamScript --projectDir=$projectDir";
    my $genomeBamTaskLabel = "mergeBam genome bam";
    my $genomeBamTaskId = "MERGE_BAM";
    executeSingleTask( $genomeBamCmd , $genomeBamTaskId, $genomeBamTaskLabel, $allChrCheckPointName );

    {
        my $allRealignedTaskLabel      = "realigned genomeBam";
        my $allRealignedCheckPointName = "REALIGNED_GENOMEBAM_PREREQ";
        my $allRealignedCheckPointRef  = checkPointInit( $allRealignedCheckPointName,
                                                         $allRealignedTaskLabel );

        tryAddTask2CheckPoint( "VARIANTS", %{$allRealignedCheckPointRef} );
        executeCheckPoint( %$allRealignedCheckPointRef );

        my $realignedBamCmd = "$genomeBamCmd --realignedBam";
        my $realignedBamTaskLabel = "mergeBam realigned genome bam";
        my $realignedBamTaskId = "MERGE_REALIGNED_BAM";
        executeSingleTask( $realignedBamCmd , $realignedBamTaskId, $realignedBamTaskLabel, $allRealignedCheckPointName);
    }

    return 1 if($isBamSkipRefSeq);

    my $genomeFastaScript = File::Spec->catfile($scriptDir, 'makeGenomeBamFasta.pl');
    my $genomeFastaCmd = "$genomeFastaScript --projectDir=$projectDir";
    my $genomeFastaTaskLabel = "mergeBam genome fasta";
    my $genomeFastaTaskId = "MERGE_BAM_FASTA";
    executeSingleTask( $genomeFastaCmd , $genomeFastaTaskId, $genomeFastaTaskLabel, $taskIdRefSeq );
}


#------------------------------------------------------------------------------

1;
__END__

# =pod
#
# =head1 CONFIGURATION AND ENVIRONMENT
#
# Uses CASAVA build configuration data.
#
# =head1 DEPENDENCIES
#
# =over 4
#
# =item Standard perl modules
#
# strict, warnings
#
# =item External perl modules
#
# =item Casava perl modules
#
# =head1 BUGS AND LIMITATIONS
#
# Please report problems to Illumina Technical Support (support@illumina.com)
# 
# Patches are welcome.
#
# =head1 AUTHOR
#
# Original version by Roman Petrovski,
# Made pluggable by Richard Shaw.
# Production plugin by Chris Saunders
#
# =cut

#------------------------------------------------------------------------------
