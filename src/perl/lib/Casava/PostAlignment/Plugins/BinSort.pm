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

sort - Bin alignments into genomic segments and sort by
position.

=head1 DESCRIPTION

The sort module takes all alignments provided as input for a
post-alignment CASAVA build, bins these into contiguous genomic
segments and sorts the alignments by start position. This operation is
a prerequisite for all subsequent variant discovery and analysis
targets in CASAVA.

For paired-end DNA-Seq builds, the sort module will also perform PCR
duplicate removal by default (this can be optionally disabled).

For RNA-Seq analyses the sort module is responsible for re-mapping
spliced read alignments produced by ELAND-RNA into genomic coordinates.

The sort operation takes alignments in Illumina 'export' format and
transforms these into BAM format during the sort process. Alignments
in BAM format are required for CASAVA's downstream variant-analysis
code. The sort module may be combined with the 'bam' module (which
aggregates all alignments in the build to a single bam file) to
efficiently produce a single sorted and indexed BAM file from a set of
export files.

The sort module has two modes, filtered and archival. Filtered mode is
run by default.

In filtered mode, only alignments or alignment pairs which are usable
for variant calling are written into the final sorted BAM files. Thus,
all alignments which fail the purity filter or are marked as PCR
duplicates are filtered out. Additionally any alignment which is not
mapped will be filtered out of single-end read data and any
alignment pairs where both reads are not mapped will be filtered out
of paired-end data.

In archival mode, no read filtration takes place. Any reads or read
pairs which are not mapped as described above are placed in a special
subdirectory of the build's parsed folder called 'notMapped'. All
reads which either fail purity filtration or are determined to be PCR
duplicates have the appropriate bit marked in the BAM FLAG field, but
are otherwise placed in the output BAM files following the same
procedure used for all other alignments: mapped data are placed in the
sorted BAM file for the appropriate genomic segment, and unmapped data
are placed into the appropriate BAM file in the 'notMapped' directory.

The sort module does not trim reads in either mode. Thus, an archival
build will contain all of the build's input reads in their
entirety. Note that the alignments themselves may be trimmed in
RNA-Seq builds via BAM's 'soft-clip' indicator, but this will not trim
the read.

Using archival mode should have no effect on downstream variant
calling results because the extra data included in the build are not
used by CASAVA's variant callers. This mode will increase the cpu, I/O
and storage requirements of the sort operation.

Note that if the 'bam' module is run following an archival sort
operation, it will produce a single BAM file which is also archival in
the sense that it contains all reads provided as input to the CASAVA
build.

=head1 OPTIONS

=over 4

=item --sortKeepAllReads

Run the sort module in archival mode instead of the default filtered
mode.

=item --rmDup=YES|NO

Turn On/Off PCR duplicate marking/removal for paired-end reads
(default YES).

=item --sortBufferSize=INTEGER

Buffer size used by the read sorting process, in megabytes (default:
1984).

=back

=begin comment

=head2 RD OPTIONS

=over 4

=item --sortNoCompressPair

Disable compression of paired export files between the first split and
duplicate removal steps.

=back

=end comment

=cut


package Casava::PostAlignment::Plugins::BinSort;

use warnings FATAL => 'all';
use strict;

use File::Spec;
# use Getopt::Long;
# use Sys::Hostname;

use lib '@CASAVA_FULL_PERL_LIBDIR@';
use Casava::Common::Log;
use Casava::Common::IOLib qw(executeCmd);
use Casava::PostAlignment::Sequencing::Config qw(%CONF_APP isSpliceJunctionChrom);
use Casava::PostAlignment::Sequencing::SamLib qw(getSamtoolsBin);
# use Casava::PostAlignment::Sequencing::SamLib qw(checkSamtoolsBin);

# use Casava::PostAlignment::Plugins::RefSeq qw(getDoneCheckPoint);
use Casava::TaskManager qw(executeSingleTask checkPointInit
    executeCheckPoint
    addTask2checkPoint);

sub executeOrGenerateWorkflow ($$\%\%);

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

sub getTarget { return 'sort'; }

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

    return (
            'sortKeepAllReads!' => \$PARAMS->{sortKeepAllReads},
            'rmDup=s' => \$PARAMS->{rmDup},

            # Performance Option:
            'sortBufferSize=i' => \$PARAMS->{sortBufferSize},

            # RD Option:
            'sortNoCompressPair!' => \$PARAMS->{sortNoCompressPair},

            # undocumented legacy option from configureBuild.pl -- seems to have
            # something to do with tile filtration :
            'buildFile=s' => \$PARAMS->{buildFilePath} );
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

    if (!$CONF_PROJ_Ref->{alignSortedBam})
    {
        my $cmd;
        my $libExecDir      = '@CASAVA_FULL_LIBEXECDIR@';
        my $runExportCmd    = File::Spec->catfile( $libExecDir , $CONF_APP{cmdRunExport} );
        my $runSplitCmd     = File::Spec->catfile( $libExecDir , $CONF_APP{cmdRunSplit} );
        my $isBuildFilePath = (defined $CONF_PROJ_Ref->{buildFilePath} and $CONF_PROJ_Ref->{buildFilePath});
        my $sortBufferSize  = $CONF_PROJ_Ref->{sortBufferSize};
    
        if( $sortBufferSize < 1 ) {
            errorExit("ERROR: invalid argument to sortBufferSize: '$sortBufferSize'\n");
        }
    
        # TODO - import these scripts into direct plugin code:
        $cmd = "$runExportCmd --projectDir=$projectDir";
        $cmd .= " --buildFile=" .$CONF_PROJ_Ref->{buildFilePath} if ($isBuildFilePath);
        executeCmd($cmd, 4);
    
        $cmd = "$runSplitCmd --projectDir=$projectDir --prevTarget=export";
        executeCmd($cmd, 4);
    }
    else
    {
        errorExit("ERROR: alignSortedBam requires --inSortedBam") if (!$CONF_PROJ_Ref->{inSortedBam});
        errorExit("ERROR: alignSortedBam requires --inSortedBamCmd") if(!$CONF_PROJ_Ref->{inSortedBamCmd});
        my $cmd = "$CONF_PROJ_Ref->{inSortedBamCmd} '$CONF_PROJ_Ref->{inSortedBam}' '$CONF_PROJ_Ref->{alignSortedBam}'";

        my $bamSourceIndexFilePath = "$CONF_PROJ_Ref->{inSortedBam}\.bai";
        if (!-e $bamSourceIndexFilePath)
        {
            my $samtoolsBin     = getSamtoolsBin(%$CONF_PROJ_Ref);
            
            $cmd .= " && $samtoolsBin index '$CONF_PROJ_Ref->{alignSortedBam}'";
        }
        else
        {
            my $bamTargetIndexFilePath = "$CONF_PROJ_Ref->{alignSortedBam}\.bai";
            $cmd .= " && $CONF_PROJ_Ref->{inSortedBamCmd} '$bamSourceIndexFilePath' '$bamTargetIndexFilePath'";
        }

        my @chromsBuild   = keys %$chromsBinSizesRef;

        my $taskLabel      = "sortMerge";
        my $importCheckPointName = "SYNCH_MERGE_DONE";
        my $checkPointRef  = checkPointInit( $importCheckPointName, $taskLabel );
        executeSingleTask( $cmd, "IMPORT_BAM", "Import external bam", 'N/A' );
        addTask2checkPoint( "IMPORT_BAM" , %$checkPointRef );
        executeCheckPoint( %$checkPointRef );

        for my $chrom (@chromsBuild) {
            next if (isSpliceJunctionChrom($chrom, %$CONF_PROJ_Ref));
            my $checkPointChromName = "SYNCH_MERGE_$chrom\_DONE";
            my $checkPointChromRef  = checkPointInit( $checkPointChromName, "sortMerge for $chrom" );
            
            addTask2checkPoint( $importCheckPointName , %$checkPointChromRef );
            executeCheckPoint( %$checkPointChromRef );
            
        }
        
    }
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
# Chris Saunders
#
# =cut

#------------------------------------------------------------------------------
