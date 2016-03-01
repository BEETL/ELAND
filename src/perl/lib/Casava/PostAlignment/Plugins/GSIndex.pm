# =head1 LICENSE

# Copyright (c) 2007-2009 Illumina, Inc.

# This software is covered by the "Illumina Genome Analyzer Software
# License Agreement" and the "Illumina Source Code License Agreement",
# and certain third party copyright/licenses, and any user of this
# source file is bound by the terms therein (see accompanying files
# Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
# Illumina_Source_Code_License_Agreement.pdf and third party
# copyright/license notices).

# This file is part of the Consensus Assessment of Sequence And VAriation
# (CASAVA) software package.

=head1 SUMMARY

gsIndex - Create Genome Studio linear index for BAM files.

=head1 DESCRIPTION

This module creates a custom linear index for the BAM files in a build
which can be used by Genome Studio to provide very fast lookup and
multi-resolution read depth summaries. Note that this index is made in
addition to the standard samtools BAM index, which is already created
by default with each BAM file in CASAVA. The linear index is not
required to view reads in Genome Studio.

B<IMPORTANT:> the I<sort> target must have been completed or be part
of the current workflow in order for the I<gsIndex> target to
work

=begin comment

=head1 OPTIONS

=over 4

=item --indelsSaveTempFiles

Add this flag to save intermediate output files from each stage of the
indel assembly process.

=back

=end comment

=cut

#
## Not used in latest version of the indel finder:
#
# =item --indelsSrasThreshold=i
#
# Single read alignment score threshold. Default 10.
#
#
#=pod
#
#=head1 SUBROUTINES
#
#=cut
#


package Casava::PostAlignment::Plugins::GSIndex;

use strict;
use warnings "all";

use Casava::Common::Log;
use Casava::PostAlignment::Sequencing::Config qw(%CONF_APP isSpliceJunctionChrom);
use Casava::PostAlignment::Plugins::RefSeq qw(getDoneCheckPoint);
use Casava::TaskManager qw(executeSingleTask executeArrayTaskEx checkPointInit
                           executeCheckPoint addTask2checkPoint tryAddTask2CheckPoint);

#
#=pod
#
#=head2 getTarget()
#
#All plugins are required to implement this method. The framework will execute the
#C<executeOrGenerateWorkflow> if user has requested the target to be added to the
#build.
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
sub getTarget{ return 'gsIndex';}

#
#=pod
#
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

    return ();
}



#
#=pod
#
#=head2 executeOrGenerateWorkflow
#
#All plugins are required to implement this method. The framework calls it when 
#the target returned by L</getTarget> has been requested to be added to the CASAVA 
#build.
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
#is expected to specify dependencies on the other target tasks (if needed) 
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

sub executeOrGenerateWorkflow ($$\%\%)
{
    my ( $prevTarget, $projectDir, $CONF_PROJ_Ref, $chromsBinSizesRef ) = @_;

    my $libexecDir      = '@CASAVA_FULL_LIBEXECDIR@';
    my $gsIndexer       = File::Spec->catfile($libexecDir, 'makeGStudioBAMIndex.pl');
    my $dirCurrentBuild = $CONF_PROJ_Ref->{dirBuildParsed};
    my $bamDirName      = $CONF_APP{dirBam};
    my $bam_f           = 'sorted.bam';

    my $taskIdRefSeq = Casava::PostAlignment::Plugins::RefSeq::getDoneCheckPoint();

    my $checkPointName = "GSINDEX";
    my $checkPointRef  = checkPointInit($checkPointName , "GS linear BAM index" );

    my @chroms = keys %$chromsBinSizesRef;
    for my $chrom ( @chroms ) {
        next if ( isSpliceJunctionChrom($chrom, %{$CONF_PROJ_Ref}) );

        my $checkPointNamePrereq = "PREREQ_GSINDEX_$chrom";
        my $checkPointPrereq  = checkPointInit( $checkPointNamePrereq, "GS linear BAM index prereq for chrom: $chrom" );

        addTask2checkPoint( $taskIdRefSeq , %{$checkPointPrereq} );

        # dependent on whole BAM file now:
        my $taskIdSplitSameBinPrereq = "SYNCH_MERGE_$chrom\_DONE";
        tryAddTask2CheckPoint( $taskIdSplitSameBinPrereq, %{$checkPointPrereq} );

        executeCheckPoint( %{$checkPointPrereq} );

        my $chromDir = File::Spec->catdir($dirCurrentBuild, $chrom);
        my $bamPath  = File::Spec->catfile($chromDir,$bamDirName,$bam_f);

        my $cmd = "$gsIndexer --bam=$bamPath";
        executeArrayTaskEx( $cmd, "GSINDEX_$chrom", "GS linear BAM index for $chrom", $checkPointNamePrereq , %$checkPointRef );
    }
    executeCheckPoint( %$checkPointRef );
}

1;
__END__

# =head1 CONFIGURATION AND ENVIRONMENT
# =head1 DEPENDENCIES
# =over 4
# =item Standard perl modules
# strict, warnings, 
# =item External perl modules
# =item Casava perl modules
# Casava::Common::Log, 
# Casava::PostAlignment::Sequencing::Config, 
# Casava::TaskManager, 
# Casava::Common::IOLib
# =back

# =pod

# =head1 BUGS AND LIMITATIONS

# There are no known bugs in this module.

# All documented features are fully implemented.

# Please report problems to Illumina Technical Support (support@illumina.com)

# Patches are welcome.

# =cut
