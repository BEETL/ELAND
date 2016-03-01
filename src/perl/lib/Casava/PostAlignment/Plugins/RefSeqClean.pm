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

Casava::PostAlignment::Plugins::RefSeqClean - Plugin for reference sequence
cache optimization cleanup.

=head1 OPTIONS

Plugin options:

=over 4

=item *

--refSeqCacheDirName=s

Name of the folder that is created under the project folder to store the
reference fasta files

=item *

--keepCachedRefSeq

If set, the cached reference sequence copy is not destroyed at the workflow end

=back

=head1 DESCRIPTION

Frees up disk space allocated for reference sequence cache. Add as post-workflow
target to minimize data footprint. Especially useful for small builds.

=cut
#=pod
#
#=head1 SUBROUTINES
#
#=cut


package Casava::PostAlignment::Plugins::RefSeqClean;

use strict;
use warnings "all";

use Casava::Common::Log;
use Casava::TaskManager qw(checkPointInit executeSingleTask addEverything2CheckPoint executeCheckPoint);

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
sub getTarget{ return 'refSeqClean';}

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
    return ('refSeqCacheDirName=s' => \$PARAMS->{refSeqCacheDirName},
            'keepCachedRefSeq'     => \$PARAMS->{keepCachedRefSeq});
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


    if ('contigName' eq $CONF_PROJ_Ref->{chromNameSource} && !$CONF_PROJ_Ref->{keepCachedRefSeq})
    {
        my $taskLabel      = "Cleanup reference cache prerequisites.";
        my $checkPointName = getTarget().'_PREREQ';
        my $checkPointRef  = checkPointInit( $checkPointName, $taskLabel );
        if (defined $checkPointRef)
        {
            addEverything2CheckPoint(%$checkPointRef);
            executeCheckPoint(%$checkPointRef);

            my $cachedSeqFilesDir = File::Spec->catdir($projectDir, $CONF_PROJ_Ref->{refSeqCacheDirName});

            my $cmd = "$0 --outDir=$projectDir --targets continue"
                      . " --chromNameSource=contigName"
                      . " --refSequences=$CONF_PROJ_Ref->{dirRefSeq}"
                      . " --maskRefSeqFiles=$CONF_PROJ_Ref->{maskRefSeqFiles}"
                    . " && rm -rf $cachedSeqFilesDir ";
            executeSingleTask( $cmd, 'SPLIT_REF_SEQ_CLEANUP', 'Split reference cleanup', $checkPointName);
        }
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
