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

assembleIndels - search for candidate indels using local de-novo read assembly

=head1 DESCRIPTION

This module implements the GROUPER algorithm, which searches for
clusters of unaligned read pairs and poorly aligned reads in the same
region and attempts to de-novo assemble these into a contig to be
aligned back the genome.

B<IMPORTANT:> the I<sort> target must have been completed or be part
of the current workflow in order for the I<assembleIndels> target to
work

=head1 OPTIONS

=over 4

=item --indelsSpReadThresholdIndels=NUMBER

Spanning read score threshold. The higher the score, the more
unlikely it is to see this pattern of mismatches given the read's
quality values. Default value is 25. Drop this value to add more
reads into the indel finding process, at the possible expense of
introducing noise. For an alignment with no mismatches this
metric has a value of zero.

=item --indelsPrasThreshold=i

Paired read alignment score threshold. If a read has a paired read
alignment score of at least this, then it is used to update the base
quality stats for that sample prep.

=item --indelsAlignScoreThresh=NUMBER

If an alignment score for a read exceeds this threshold then
the output file is updated to incorporate this alignment.
Otherwise the read's entry remains as per the input file.
Default value is 120.
A low value will cause some reads to be wrongly placed
(albeit within a small interval).

=item --indelsNumLowInsertSizeSds

If the insert size of a read pair is more than this number of SDs
below the median the reads are considered as evidence for an insertion.
(default: 5.0).

=item --indelsNumHighInsertSizeSds

If the insert size of a read pair is more than this number of SDs
above the median the reads are considered as evidence for a deletion.
(default: 3.0).

=item --indelsSdFlankWeight=NUMBER

Number of standard deviations to use when defining the
genomic interval to align the read to (default: 1).

=item --indelsMinGroupSize=NUMBER

Only output clusters if they contain at least this many reads.

=item --indelsSpReadThresholdClusters=NUMBER

Spanning read score threshold. This is calculated in exactly the same
way as --indelsSpReadThresholdIndels. However it is used in the
opposite way. Here the point to find reads with few or no mismatches,
which are presumed to arise from repeats and not from indels, and
exclude them from the clustering process.

=item --indelsMaxDistAligned=NUMBER

Max distance between semi-aligned reads for them to be in same cluster
(default 5.0)

=item --indelsMaxDistInterStrand=NUMBER

Max distance between mean positions of clusters of semi-aligned (& shadow)
reads aligned to different strands for them to be merged
(default 12.0)

=item --indelsMinCoverage=NUMBER

min coverage to extend contig (default 3)

=item --indelsMinContext=NUMBER

Demand at least x exact matching bases either side of variant (default is 6).
The idea here is to ensure that an indel has a minimum number of exactly matching
bases on either side. Setting this to zero might be good for finding reads which
align to breakpoints.

=item --indelsMinTranslocPriAlignScore=NUMBER

Min alignment score for candidate translocation contig to primary chromosome
search region

=item --indelsMinTranslocSecAlignScore=NUMBER

Min alignment score for candidate translocation contig to secondary chromosome
search region

=item --indelsSaveTempFiles

Add this flag to save intermediate output files from each stage of the
indel assembly process.

=item --indelsReadStartDepthCutoff=NUMBER

Ignore regions where more than NUMBER of reads start at the same position.
Default is 0 (no cutoff).

=back

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


package Casava::PostAlignment::Plugins::IndelAssembler;

use strict;
use warnings "all";
#use Carp;
#use Exporter 'import';

use Casava::Common::Log;
use Casava::PostAlignment::Sequencing::Config qw(%CONF_APP %chrEnds isSpliceJunctionChrom);
use Casava::Common::IOLib qw(createDirs);
use Casava::PostAlignment::Plugins::RefSeq qw(getDoneCheckPoint);
use Casava::TaskManager qw(executeSingleTask checkPointInit
                           executeCheckPoint addTask2checkPoint tryAddTask2CheckPoint);
use File::Spec;
use FileHandle;

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
sub getTarget{ return 'assembleIndels';}

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

    $PARAMS->{indelsErrorRate} = "";

    return ("indelsSrasThreshold=i"             => \$PARAMS->{indelsSrasThreshold},
            "indelsPrasThreshold=i"             => \$PARAMS->{indelsPrasThreshold},
            "indelsSpReadThresholdIndels=s"     => \$PARAMS->{indelsSpReadThresholdIndels},
            "indelsSdFlankWeight=i"             => \$PARAMS->{indelsSdFlankWeight},
            "indelsAlignScoreThresh=i"          => \$PARAMS->{indelsAlignScoreThresh},
            "indelsNumLowInsertSizeSds=f"       => \$PARAMS->{indelsNumLowInsertSizeSds},
            "indelsNumHighInsertSizeSds=f"      => \$PARAMS->{indelsNumHighInsertSizeSds},
            "indelsReadStartDepthCutoff=i"      => \$PARAMS->{indelsReadStartDepthCutoff},
            "indelsMaxDistAligned=f"            => \$PARAMS->{indelsMaxDistAligned},
            "indelsMaxDistInterStrand=f"        => \$PARAMS->{indelsMaxDistInterStrand},
            "indelsMinGroupSize=i"              => \$PARAMS->{indelsMinGroupSize},
            "indelsSpReadThresholdClusters=i"   => \$PARAMS->{indelsSpReadThresholdClusters},
            "indelsMinCoverage=i"               => \$PARAMS->{indelsMinCoverage},
            "indelsMinContext=i"                => \$PARAMS->{indelsMinContext},
            "indelsMinTranslocPriAlignScore=i"  => \$PARAMS->{indelsMinTranslocPriAlignScore},
            "indelsMinTranslocSecAlignScore=i"  => \$PARAMS->{indelsMinTranslocSecAlignScore},
            "indelsSaveTempFiles!"              => \$PARAMS->{indelsSaveTempFiles},
            "CNVsegDisable!"                    => \$PARAMS->{CNVsegDisable});
}

#------------------------------------------------------------------------------


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

    my $libexecDir          = '@CASAVA_FULL_LIBEXECDIR@';
    my $preTradeGrouperBin  = File::Spec->catfile($libexecDir,
                                                 'runPreTradeGrouperBin.pl');
    my $grouperTrading      = File::Spec->catfile($libexecDir,
                                                  'runGrouperTrading.pl');
    my $postTradeGrouperBin = File::Spec->catfile($libexecDir,
                                                 'runPostTradeGrouperBin.pl');
    my $dataSetSuffix       = $CONF_PROJ_Ref->{dataSetSuffix};

    my @chroms = keys %$chromsBinSizesRef;

    my @dirs = ();
    for my $chrom ( @chroms ) {
        my $indelDir = File::Spec->catdir( $chrom, "Indel" );
        push @dirs, $indelDir;
    } # foreach

    createDirs($CONF_PROJ_Ref->{dirBuildParsed}, @dirs);

    my $taskIdRefSeq
      = Casava::PostAlignment::Plugins::RefSeq::getDoneCheckPoint();

    my $preTradeCheckPointName = "GROUPER_PRE_TRADE";
    my $preTradeCheckPointRef
      = checkPointInit($preTradeCheckPointName,
                       "Grouper pre-trade chrom-level processing");

    for my $chrom ( @chroms ) {
        next if ( isSpliceJunctionChrom($chrom, %{$CONF_PROJ_Ref}) );

        my $preTradeCheckPointChromName = "GROUPER_PRE_TRADE_$chrom";
        my $preTradeCheckPointChromRef
          = checkPointInit($preTradeCheckPointChromName,
                           "Grouper pre-trade bin-level processing");

        for ( my $i = 0 ; $i < $chromsBinSizesRef->{$chrom} ; $i++ ) {
            my $binId = sprintf "%04d", $i;

            my $cmd = "$preTradeGrouperBin --projectDir=$projectDir "
              . "--chrom=$chrom --binId=$binId";
            if ($CONF_PROJ_Ref->{alignSortedBam})
            {
                $cmd .= " --bamFile='$CONF_PROJ_Ref->{alignSortedBam}'";
            }

            # enumerate all prerequisite tasks for this bin:
            #
            my $checkPointNamePrereq = "PREREQ_GROUPER_$chrom/$binId";
            my $checkPointPrereq
              = checkPointInit( $checkPointNamePrereq,
                                "Grouper prereq for $chrom/$binId" );
            addTask2checkPoint( $taskIdRefSeq , %{$checkPointPrereq} );
            my $depId ="SYNCH_MERGE_$chrom\_DONE";
            tryAddTask2CheckPoint( $depId , %{$checkPointPrereq} );
            executeCheckPoint( %{$checkPointPrereq} );

            my $taskId = "GROUPER_PRE_TRADE_$chrom/$binId";

            executeSingleTask($cmd, $taskId,
                              "Grouper pre-trade $chrom $binId",
                              $checkPointNamePrereq);
            addTask2checkPoint($taskId, %$preTradeCheckPointChromRef);
        }

        executeCheckPoint(%$preTradeCheckPointChromRef);
        addTask2checkPoint($preTradeCheckPointChromName,
                           %$preTradeCheckPointRef);
    }

    executeCheckPoint(%$preTradeCheckPointRef);


    # Trading

    my $tradingCmd = "$grouperTrading --projectDir=${projectDir}";
    my $tradingTaskId = "INDEL_FINDER_READ_BROKER";
    executeSingleTask($tradingCmd, $tradingTaskId,
                      "GROUPER Trading", $preTradeCheckPointName);


    # Post-trading

    # Overall completion name used externally so do not qualify.
    my $postTradeCheckPointName = "GROUPER";
    my $postTradeCheckPointRef
      = checkPointInit($postTradeCheckPointName,
                       "Grouper post-trade chrom-level processing");

    for my $chrom ( @chroms ) {
        next if ( isSpliceJunctionChrom($chrom, %{$CONF_PROJ_Ref}) );

        my $postTradeCheckPointChromName = "GROUPER_$chrom";
        my $postTradeCheckPointChromRef
          = checkPointInit($postTradeCheckPointChromName,
                           "Grouper post-trade bin-level processing");

        for (my $binNum = 0; $binNum < $chromsBinSizesRef->{$chrom};
             ++$binNum) {
            my $binId = sprintf "%04d", $binNum;

            my $cmd = "$postTradeGrouperBin --projectDir=$projectDir "
              . "--chrom=$chrom --binId=$binId";
            if ($CONF_PROJ_Ref->{alignSortedBam})
            {
                $cmd .= " --bamFile='$CONF_PROJ_Ref->{alignSortedBam}'";
            }
            
            my $taskId = "GROUPER_$chrom/$binId";
            executeSingleTask($cmd, $taskId,
                              "Grouper post-trade $chrom $binId",
                              $tradingTaskId);
            addTask2checkPoint($taskId, %$postTradeCheckPointChromRef);
        }

        executeCheckPoint(%$postTradeCheckPointChromRef );
        addTask2checkPoint($postTradeCheckPointChromName,
                           %$postTradeCheckPointRef);
    }

    executeCheckPoint(%$postTradeCheckPointRef);
}

#------------------------------------------------------------------------------

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

=pod

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=head1 AUTHOR

Roman Petrovski

=cut
