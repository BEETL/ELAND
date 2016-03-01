=head1 LICENSE

Copyright (c) 2007-2009 Illumina, Inc.

This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

This file is part of the Consensus Assessment of Sequence And VAriation
(CASAVA) software package.

=head1 NAME

Casava::PostAlignment::Plugins::RnaSeqCounts- Plugin for RNA feature counting workflow.

=head1 OPTIONS

Plugin options:

=over 4

=item --rnaCountMethod=s

Defines which feature counting method to use. Valid methods are: B<readBases> (default).

=over 8

=item readBases

Number of read bases falling into each feature is counted except for splice
junction counts. Splice junctions are counted as number of reads overlapping the
junction point.

=back

=item --refFlatFile=s

Path to refFlat.txt.gz file containing the genome annotation in UCSC format

=item --seqGeneMdFile=s

Path to seq_gene_md.gz file containing the genome annotation in NCBI format

=item --seqGeneMdGroupLabel=s

Assembly group label to use with --seqGeneMdFile

=item --nonoverlappingExonCoordsParams=s

Command-line parameters to be passed to nonoverlappingExonCoords.pl

=item --seqGeneMd2RefFlatParams=s

Command-line parameters to be passed to seqGeneMd2RefFlat.pl

=back

=head1 DESCRIPTION

This script configures workflow for calculating the exon, gene and splice 
junction expression.

B<IMPORTANT:> the I<sort> and I<callSmallVariants> targets must have
been completed or be part of the current workflow

=cut
#=pod
#
#=head1 SUBROUTINES
#
#=cut



package Casava::PostAlignment::Plugins::RnaSeqCounts;

use strict;
use warnings 'all';

use Casava::Common::Log;
use Casava::Common::IOLib;
use Casava::PostAlignment::Sequencing::Config qw(%CONF_APP %CONF_PROJ %chrEnds
    isSpliceJunctionChrom getSpliceJunctionFileName);
use Casava::PostAlignment::Sequencing::GenomicIOLib qw(%readsIdxFields);
use Casava::PostAlignment::Sequencing::RnaSeqLib qw($tmpSpliceCountFileName);
use Casava::TaskManager qw(executeSingleTask checkPointInit
    executeArrayTaskEx executeCheckPoint addTask2checkPoint tryAddTask2CheckPoint);

my $taskLabel = "rnaCount";
my $checkPointName = "SYNCH_RNA_COUNT";

sub executeOrGenerateWorkflow ($$\%\%);

#
#=pod
#
#=head2 getTarget
#
#All plugins are required to implement this method. The framework will execute the
#L</executeOrGenerateWorkflow> if user has requested the target to be added to the
#build.
#
#I<Parameters:>
#
#=over 4
#
#=item *
#
#No parameters
#
#=back
#
#I<Returns:>
#
#Target name for which the plugin can perform workflow configuration
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
#
#=cut

sub getTarget() { return 'rnaCounts';}
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
#none
#
#=back
#
#=cut
sub getOptionsMapping(\%)
{
    my ( $PARAMS ) = @_;

    return ('rnaCountMethod|rcm=s'              => \$PARAMS->{rnaCountMethod},
            'nonoverlappingExonCoordsParams=s'  => \$PARAMS->{nonoverlappingExonCoordsParams},
            'seqGeneMd2RefFlatParams=s'         => \$PARAMS->{seqGeneMd2RefFlatParams},
            'refFlatFile=s'                     => \$PARAMS->{refFlatFile},
            'seqGeneMdFile=s'                   => \$PARAMS->{seqGeneMdFile},
            'seqGeneMdGroupLabel=s'             => \$PARAMS->{seqGeneMdGroupLabel});
}

#=head2 getPrintExonCoordsCmd
#
#Validates configuration and calculates command to produce exon_coords.txt file.
#
#I<Parameters:>
#
#None 
#
#I<Returns:>
#
#String that can be used as exon_coords file path
#
#I<Exceptions:>
#
#=over 4
#
#=item *
#
#Invalid combination of project settings
#
#=back
#
#=cut
sub getPrintExonCoordsCmd()
{
    my ( $PARAMS ) = @_;

    my $nonoverlappingExonCoords = File::Spec->catfile('@CASAVA_FULL_LIBEXECDIR@', "nonoverlappingExonCoords.pl");
    my $seqGeneMd2RefFlat = File::Spec->catfile('@CASAVA_FULL_LIBEXECDIR@', 'Conversion', 'seqGeneMd2RefFlat.pl');

    my $printExonCoordsCmd = 'cat ' . File::Spec->catfile($CONF_PROJ{genesListPath}, $CONF_PROJ{featureFileName});
    if (!$CONF_PROJ{featureFileName}) 
    {
        if (!$CONF_PROJ{refFlatFile})
        {
            if (!$CONF_PROJ{seqGeneMdFile})
            {
                errorExit "ERROR: Either --featureFileName or --refFlatFile or --seqGeneMdFile is required for RNA feature counting";
            }
            else
            {
                errorExit "ERROR: --seqGeneMdGroupLabel is required for $CONF_PROJ{seqGeneMdFile}" unless $CONF_PROJ{seqGeneMdGroupLabel};
                $printExonCoordsCmd = "$CONF_PROJ{cmdGunzip} -c '$CONF_PROJ{seqGeneMdFile}'"
                                . " | $seqGeneMd2RefFlat --input-file - --group-label '$CONF_PROJ{seqGeneMdGroupLabel}' $CONF_PROJ{seqGeneMd2RefFlatParams}"
                                . " | $nonoverlappingExonCoords --chromNameSource=$CONF_PROJ{chromNameSource} $CONF_PROJ{nonoverlappingExonCoordsParams}";
            }
        }
        else
        {
            $printExonCoordsCmd = "$CONF_PROJ{cmdGunzip} -c '$CONF_PROJ{refFlatFile}'"
                            . " | $nonoverlappingExonCoords --chromNameSource=$CONF_PROJ{chromNameSource} $CONF_PROJ{nonoverlappingExonCoordsParams}";
        }
    }
    else
    {
        if( $CONF_APP{spliceJunctionAuto} eq $CONF_PROJ{spliceJunction} )
        {
            errorExit "--spliceJunction must be specified if --featureFileName is set";
        }
        else
        {
            my $genesListPath = $CONF_PROJ{genesListPath};
            unless (-d $genesListPath)
            {
                errorExit "--genesListPath must point to an existing folder. Current: $genesListPath";
            }
        
            my $featureFileNamePath = File::Spec->catfile( $genesListPath, $CONF_PROJ{featureFileName} );
        
            if ( !-e $featureFileNamePath ) {
                errorExit("ERROR: $featureFileNamePath does not exist. Please check --featureFileName paremeter.");
            }
            
        }
    }
    return $printExonCoordsCmd;
}

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
#none
#
#=back
#
#    $prevTarget         - previous target (dependent - default empty)    
#    $projectDir         - HASH MAP Ref allele caller configuration    
#    $CONF_PROJ_Ref->_Ref      - HASH MAP Ref to project configuration
#    $chromsBinSizesRef  - HASH MAP Ref to chromosome sizes    
#    
#B<Returns:> 
#    Nothing
#    
#=cut

sub executeOrGenerateWorkflow ($$\%\%) {
    my ( $prevTarget, $projectDir, $CONF_PROJ_Ref, $chromsBinSizesRef ) = @_;

    if ($CONF_PROJ{applicationType} eq 'RNA')
    {
        errorExit "ERROR: rnaCounts cannot be added incrementally to an existing build" unless ( 'start' ne $prevTarget ); 
        errorExit "ERROR: Unsupported rnaCountMethod $CONF_PROJ_Ref->{rnaCountMethod}" unless ('readBases' eq $CONF_PROJ_Ref->{rnaCountMethod});
    
        my $rnaSpliceCountsBin = File::Spec->catfile('@CASAVA_FULL_LIBEXECDIR@', "DataAnalysis", "rnaSpliceCounts.pl");
        my $rnaExonCountsBin = File::Spec->catfile('@CASAVA_FULL_LIBEXECDIR@', "DataAnalysis", "rnaExonCounts.pl");
        my $exonCoordsProjectPath = File::Spec->catdir( $projectDir, 'features', 'exon_coords.txt' );
        my $readsIdxFilePath = File::Spec->catfile( $projectDir, 'stats', $CONF_APP{f_reads_indx} );
        
        my $dirCurrentBuild = $CONF_PROJ{dirBuildParsed};
        my $checkPointRef  = checkPointInit( $checkPointName, $taskLabel );
    
        my $printExonCoordsCmd = getPrintExonCoordsCmd();
        my $makeExonCoordsCmd = "$printExonCoordsCmd > '$exonCoordsProjectPath'";
        executeArrayTaskEx($makeExonCoordsCmd, "EXON_COORDS", "making exon_coords.txt", 'N/A', %$checkPointRef);
    
        my $taskIdSplit = "SYNCH_EXPORT3";
        tryAddTask2CheckPoint($taskIdSplit, %$checkPointRef);
    
        executeCheckPoint( %{$checkPointRef} );
    
        my $cutFieldReads = $readsIdxFields{goodReads} + 1;
        my $chromTotalUsedReadsExpr = '$(( `' . "grep -vE '^#' '$readsIdxFilePath' |cut -f$cutFieldReads |paste -sd+" . '` ))' ;
        my $cutFieldBases = $readsIdxFields{goodReadBases} + 1;
        my $chromTotalUsedBasesExpr = '$(( `' . "grep -vE '^#' '$readsIdxFilePath' |cut -f$cutFieldBases |paste -sd+" . '` ))' ;
    
        my $QVCutoffSplice = $CONF_PROJ{QVCutoffSingle};
    
        my @chroms = keys(%{$chromsBinSizesRef});
        for my $chrom (@chroms)
        {
            next if ( isSpliceJunctionChrom($chrom, %$CONF_PROJ_Ref));
    
            my $chromDir = File::Spec->catdir( $dirCurrentBuild, $chrom);
    
            my @sitesFiles;
            my $binCount = $chromsBinSizesRef->{$chrom};
            for ( my $i = 0 ; $i < $binCount; ++$i ) {
                my $binId = sprintf "%04d", $i;
                my $sitesFile = File::Spec->catfile( $chromDir, $binId, "sites.txt.gz" );
                push @sitesFiles, $sitesFile;
            }
    
            my $checkPointChromPrereqName = "PREREQ_RNA_COUNT_$chrom";
            my $checkPointChromPrereq  = checkPointInit( $checkPointChromPrereqName, "Prereq for rna-counting on chrom $chrom");
            addTask2checkPoint( $checkPointName, %{$checkPointChromPrereq} );
            my $taskIdVariants = "VARIANTS_$chrom";
            tryAddTask2CheckPoint( $taskIdVariants, %{$checkPointChromPrereq} );
    
            executeCheckPoint( %{$checkPointChromPrereq} );
    
            my $exonCountsFileName = File::Spec->catfile( $chromDir, "$chrom\_exon_count.txt");
            my $geneCountsFileName = File::Spec->catfile( $chromDir, "$chrom\_genes_count.txt");
            my $rnaExonCountsCmd = "$rnaExonCountsBin"
                ." --chrom='$chrom'"
                ." --featureFile='$exonCoordsProjectPath'"
                ." --geneCountsFile='$geneCountsFileName'"
                ." --totalBases=$chromTotalUsedBasesExpr"
                ." " . join(' ',@sitesFiles) ." > '$exonCountsFileName'";
    
            executeSingleTask( $rnaExonCountsCmd, "rnaExonCounts$chrom",
                               "Calculate exon and gene counts RPKM for $chrom", $checkPointChromPrereqName);
    
            my $tmpCountsPath = File::Spec->catfile( $chromDir, $tmpSpliceCountFileName);
            my $spliceCountsFileName = File::Spec->catfile( $chromDir, "$chrom\_splice_count.txt");
            my $rnaSpliceCountsCmd = "$rnaSpliceCountsBin"
                ." --chrom='$chrom'"
                ." --totalReads=$chromTotalUsedReadsExpr"
                ." <'$tmpCountsPath' > '$spliceCountsFileName' && rm '$tmpCountsPath'";
    
            executeSingleTask( $rnaSpliceCountsCmd, "rnaSpliceCounts$chrom",
                               "Calculate splice RPKM for $chrom", $checkPointName);
        }
    }
    else
    {
        return 'N/A';
    }
}
1;
__END__

=pod

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over 4

=item Standard perl modules

strict, warnings, File::Spec 

=item External perl modules

=item Casava perl modules

Casava::Common::Log, 
Casava::PostAlignment::Sequencing::Config, 
Casava::TaskManager, 
Casava::Common::IOLib

=back

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=head1 AUTHOR

Roman Petrovski

=cut
