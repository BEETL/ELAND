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

Casava::PostAlignment::Plugins::RefSeq - Plugin for reference sequence
cache optimization.

=head1 OPTIONS

Plugin options:

=over 4

=item *

--refSeqCacheDirName=s

Name of the folder that is created under the project folder to store the
reference fasta files

=back

=head1 DESCRIPTION

In case when --chromNameSource=contigName, the B<refSeq> target creates
a copy of the reference sequence in single-entry fasta file format and
upldates the project --refSequences and --maskRefSeqFiles to point to the 
project copy of the reference. 

=cut
#=pod
#
#=head1 SUBROUTINES
#
#=cut


package Casava::PostAlignment::Plugins::RefSeq;

use strict;
use warnings "all";
use Exporter 'import';
our @EXPORT_OK = qw(getRuntimeRefSeqFile getDoneCheckPoint);

use Casava::Common::Log;
use Casava::Common::IOLib  qw(createDir);
use Casava::PostAlignment::Plugins qw(getPluginScriptsDir);
use Casava::TaskManager qw(checkPointInit executeArrayTaskEx executeCheckPoint tryAddTask2CheckPoint);


sub getRuntimeRefSeqFile($\%$);

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
sub getTarget{ return 'refSeq';}

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
    return ('refSeqCacheDirName=s' => \$PARAMS->{refSeqCacheDirName});
}

sub getDoneCheckPoint()
{
    return getTarget() . '_DONE';
}

sub getRuntimeRefSeqFile($\%$)
{
    my ( $projectDir, $CONF_PROJ_Ref, $chrom) = @_;
    
    my $dir = ('contigName' eq $CONF_PROJ_Ref->{chromNameSource}) 
        ? File::Spec->catdir($projectDir, $CONF_PROJ_Ref->{refSeqCacheDirName})
        : $CONF_PROJ_Ref->{dirRefSeq};

    my $ret = File::Spec->catfile($dir, $chrom);
    return $ret;
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
#  Input reference files not found
#
#=item *
#
#  Failed to create folder for split reference
#
#=back
#
#=cut

sub executeOrGenerateWorkflow($$\%\%)
{
    my ( $prevTarget, $projectDir, $CONF_PROJ_Ref, $chromsBinSizesRef ) = @_;

    my $taskLabel      = "Split reference into single-entry fastas.";
    my $checkPointName = getDoneCheckPoint();
    my $checkPointRef  = checkPointInit( $checkPointName, $taskLabel );

    if (defined $checkPointRef)
    {
        if ('contigName' eq $CONF_PROJ_Ref->{chromNameSource})
        {
            my $cachedSeqFilesDir = File::Spec->catdir($projectDir, $CONF_PROJ_Ref->{refSeqCacheDirName});
            createDir($cachedSeqFilesDir);
            # checkpoint does not exist
            my $globMask = File::Spec->catfile( $CONF_PROJ_Ref->{dirRefSeq}, $CONF_PROJ_Ref->{maskRefSeqFiles});
            my @files = glob($globMask);
            errorExit "ERROR: Could not find a single reference in $globMask" if (!scalar(@files));
            my $filesString = join (' ', @files);

            my $filesString_dir = $CONF_PROJ_Ref->{dirRefSeq};
            my $copyAlignability_files_cmd = "cp -v $filesString_dir/*forward_seqability_* $cachedSeqFilesDir";
            my $transformFastaCmd = File::Spec->catfile('@CASAVA_FULL_LIBEXECDIR@', "transformFasta.pl");

            my $cmd = "$copyAlignability_files_cmd; $transformFastaCmd --chromNameSource=$CONF_PROJ_Ref->{chromNameSource} --targetPath=$cachedSeqFilesDir $filesString "
                    . "&& $0 --chromNameSource=fileName --refSequences=$cachedSeqFilesDir --maskRefSeqFiles='*' --outDir=$projectDir --targets continue";

            # SYNCH_EXPORT1 is the task that does some validation before export starts. In particular
            # it will cause the genomesize.xml file validation to fail after we switch the reference
            # sequence to fileName. Make sure we run after it is satisified.
            my $checkPointNamePrereq = "PREREQ_$checkPointName";
            my $checkPointPrereq  = checkPointInit( $checkPointNamePrereq, "RefSeq prerequisites" );
            tryAddTask2CheckPoint( 'SYNCH_EXPORT1', %{$checkPointPrereq} );
            executeCheckPoint( %{$checkPointPrereq} );

            executeArrayTaskEx( $cmd, 'SPLIT_REF_SEQ', 'Split reference', $checkPointNamePrereq, %$checkPointRef );
        }
        executeCheckPoint(%$checkPointRef);
    }
}

1;
__END__

=pod

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over 4

=item Standard perl modules

strict, warnings, Pod::Usage, File::Spec

=item External perl modules

=item Casava perl modules

Casava::Common::Log, 
Casava::PostAlignment::Sequencing::Config,
Casava::TaskManager, 

=back

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=head1 AUTHOR

Roman Petrovski

=cut
