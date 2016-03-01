package Casava::PostAlignment::Sequencing::WorkflowManagement;

# PROJECT: NCBISubmission
# MODULE:  $RCSfile$
# AUTHOR:  Lukasz Szajkowski
#
# Copyright (c) 2008, 2009 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.
#
# The library contains procedures and variables usefull in operating
# TaskManager

=pod

=head1 NAME

Casava::PostAlignment::Sequencing::WorkflowManagement.pm - The library - TaskManager
	

=head1 SYNOPSIS

The library contains procedures and variables usefull in operating 
TaskManagers.
 
use Casava::PostAlignment::Sequencing::WorkflowManagement.pm qw();  

=head1 DESCRIPTION

Exports: 

		    
# Global variable


=head1 AUTHORSHIP

Copyright (c) 2008, 2009 Illumina
This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).


=cut

BEGIN {
    use Exporter();
    @ISA       = qw(Exporter);
    @EXPORT    = qw(&taskManagerTarget);
    @EXPORT_OK = qw();
}
use warnings FATAL => 'all';
use strict;
use POSIX qw(strftime);
use Carp;
use Sys::Hostname;

use File::Spec;

use Casava::Common::Log;
use Casava::Common::IOLib
  qw(executeCmd testIfSupports getProcessorsCount);
use Casava::TaskManager::Sge qw(isSgeQueueOk);    
=pod

=item taskManagerTarget($PARAMS_Ref, $CONF_APP_Ref, $CONF_PROJ_Ref)

Runs TaskManager target

B<Parameters:>

    $PARAMS_Ref     - 
    $CONF_APP_Ref   - 
    $CONF_PROJ_Ref  - 
    
B<Returns:> 
    nothing
    
=cut

sub taskManagerTarget (\%\%\%) {
    my ( $PARAMS_Ref, $CONF_APP_Ref, $CONF_PROJ_Ref ) = @_;

    my $taskManagerCmd = File::Spec->catfile('@CASAVA_FULL_BINDIR@', $CONF_APP_Ref->{cmdTaskManager});
    my $task2MakeCmd = File::Spec->catfile('@CASAVA_FULL_LIBEXECDIR@', 'PostAlignment', 'task2make.pl');

    my $jobsLimitLocal = (!defined $PARAMS_Ref->{jobsLimit}) ? 1 : $PARAMS_Ref->{jobsLimit};
    my $hostname = hostname();
    my $cmdLocal = undef;
    my $cmdSGE = undef;
    if ($PARAMS_Ref->{isWorkflowMake})
    {
        my $makeFile = $CONF_PROJ_Ref->{workflowFile};
        $makeFile .= 'mk' unless $makeFile =~ s/\.txt$/\.mk/;
        $makeFile = (File::Spec->splitpath($makeFile))[2];
        $makeFile = File::Spec->catfile( $CONF_PROJ_Ref->{dirBuild}, $makeFile);
        my $cmdConfigureMake = "$task2MakeCmd";
        $cmdConfigureMake .= " $PARAMS_Ref->{task2MakeParams}" if defined $PARAMS_Ref->{task2MakeParams};
        $cmdConfigureMake .= " <$CONF_PROJ_Ref->{workflowFile} >$makeFile";
        executeCmd($cmdConfigureMake);
        my $inputMakefile = File::Spec->catfile( '@CASAVA_FULL_DATADIR@', "makefiles", 'PostAlignment', 'Makefile');
        my $outputMakefile = File::Spec->catfile( $CONF_PROJ_Ref->{dirBuild}, 'Makefile');
        use File::Copy;
        copy($inputMakefile, $outputMakefile) or errorExit "ERROR: Failed to copy $inputMakefile to $outputMakefile: $!";
        executeCmd("echo 'include $makeFile' >> $outputMakefile");
        
        $cmdSGE = "qmake";
        $cmdSGE .= " -q $PARAMS_Ref->{sgeQueue}" if defined $PARAMS_Ref->{sgeQueue};
        $cmdSGE .= " $CONF_PROJ_Ref->{sgeQmakeFlags}" if defined $CONF_PROJ_Ref->{sgeQmakeFlags};
        $cmdSGE .= " $PARAMS_Ref->{sgeQmakeFlags}" if defined $PARAMS_Ref->{sgeQmakeFlags};
        $cmdSGE .= " -- ";
        $cmdSGE .= " -j $PARAMS_Ref->{jobsLimit}" if defined $PARAMS_Ref->{jobsLimit};
        $cmdSGE .= " -C $CONF_PROJ_Ref->{dirBuild}";

        $cmdLocal = "make -C $CONF_PROJ_Ref->{dirBuild}";
        $cmdSGE .= " $PARAMS_Ref->{makeFlags}" if defined $PARAMS_Ref->{makeFlags};
        $cmdLocal .= " -j $PARAMS_Ref->{jobsLimit}" if defined $PARAMS_Ref->{jobsLimit};
    }
    else
    {
        $cmdLocal = "$taskManagerCmd --tasksFile=$CONF_PROJ_Ref->{workflowFile} --host=localhost --jobsLimit=$jobsLimitLocal --mode=local ";
        $cmdSGE = "$taskManagerCmd --tasksFile=$CONF_PROJ_Ref->{workflowFile} --host=$hostname --mode=sge ";
        $cmdSGE .= " --sgeQueue=$PARAMS_Ref->{sgeQueue}" if defined $PARAMS_Ref->{sgeQueue};
        $cmdSGE .= " --jobsLimit=$PARAMS_Ref->{jobsLimit}" if defined $PARAMS_Ref->{jobsLimit};
        $cmdSGE .= " --sgeQsubFlags='$PARAMS_Ref->{sgeQsubFlags}'" if defined $PARAMS_Ref->{sgeQsubFlags};
    }
    
    if ( -e $CONF_PROJ_Ref->{workflowFile} ) 
    {
        printLog( "Task description in $CONF_PROJ_Ref->{workflowFile}.\n", 0 );
        if (!$PARAMS_Ref->{isWorkflowAuto} )
        {
            if ( testIfSupports('qstat') ) {
                logInfo( "Support for SGE (qstat) detected\n", 0 );
                logWarning "SGE queue check failed for '$PARAMS_Ref->{sgeQueue}'" 
                    unless (!$PARAMS_Ref->{sgeQueue} || isSgeQueueOk($PARAMS_Ref->{sgeQueue}));
            }
            else {
                logInfo( "Support for SGE (qstat) not detected\n", 0 );
                $cmdSGE = undef;
            }
        }
        if ( !getWarningsCount() )
        {
            if ( $PARAMS_Ref->{isWorkflowAuto} )
            {
                printLog( "Running $cmdLocal\n", 0 );
                system($cmdLocal);
                return;
            }
            elsif ( $PARAMS_Ref->{isWorkflowSGE} )
            {
                printLog( "Running $cmdSGE\n", 0 );
                system($cmdSGE);
                return;
            }
        }
        else
        {
            logWarning 'Warnings reported. --workflowAuto or --sgeAuto is not '
                .'allowed. Please review logged messages before attempting '
                .'the following command(s).' 
                unless !$PARAMS_Ref->{isWorkflowAuto} && !$PARAMS_Ref->{isWorkflowSGE};
        }

        if ( $cmdSGE ) {
            printLog( "To run on SGE, submit:\n$cmdSGE\n", 0 , 'message');
        }
        printLog( "To run locally on $jobsLimitLocal processors:\n$cmdLocal\n", 0 , 'message');

        system("echo '' >> /dev/null");
    }

}
1;    # says use was ok
__END__
