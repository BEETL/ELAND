# PROJECT: CASAVA
# MODULE:  $RCSfile$
# AUTHOR:  Lukasz Szajkowski
#
# Copyright (c) 2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.
#

=pod

=head1 NAME

Casava::TaskManager::Sge.pm - Perl wraper library for accessing Sun Grid Engine (SGE).

=head1 SYNOPSIS

# include what functions you need... 
use Casava::TaskManager::Sge qw();  

=head1 DESCRIPTION

    
# Global variable

=head1 AUTHORSHIP

Copyright (c) 2008 Illumina
This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

Created by Lukasz Szajkowski <lszajkowski@illumina.com>

This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

=head1 SOURCE CODE

The most current release of the Perl source code 
for this module is available through CVS at genie01 a.k.a. 10.44.0.81 
cvs co BullFrog

=cut

package Casava::TaskManager::Sge;

#
# Place functions/variables you want to *export*, ie be visible
# from the caller package into @EXPORT_OK
#
BEGIN {
	use Exporter();
	@ISA       = qw(Exporter);
	@EXPORT    = qw($VERSION);
	@EXPORT_OK = qw(submitJobs2SGE submit2SGE isSgeQueueOk);
}
use warnings FATAL => 'all';
use strict;
use POSIX qw(strftime);

use Carp;


use Casava::Common::Log;
use Casava::Common::IOLib qw(testIfSupports);
sub submit2SGE($);
sub isSgeQueueOk($);

=pod

=head1 The procedure submits SGE job to SGE.

=over 4

=item submitJobs2SGE

The procedure submits SGE job to SGE and returns the job id assigned to the job.

Parameters:
    $queueName        - name of sge queue
    $sgePrioritie     - sge prioritie (what procentage of tasks should be submited)
    $agentPrams       - command line parameters to pass to the agent
    $errorPath        - a path where out.txt and err.txt files will be stored
    $numberOfTaskToShedule         - number of agent jobs to launch
    $sgeQsubFlags                  - extra switches to pass to qsub
Returns:
    sge job id, -1 if $numberOfTaskToShedule is 0. 
=back

=cut 

sub submitJobs2SGE($$$$$$) {
	my (
		$queueName, 
		$sgePrioritie,   $agentPrams,
		$errorPath, $numberOfTaskToShedule, $sgeQsubFlags
	  )
	  = @_;

    if ($numberOfTaskToShedule)
    {
    	printLog( "Submit prepare\n", 6 );
    	my $queueparam = '';
    	if ( $queueName ) {
    		$queueparam = "-q '$queueName' ";
    	}    
    	my $cmd = "qsub $sgeQsubFlags -cwd -t 1-$numberOfTaskToShedule $queueparam"
    		  . "-o $errorPath.out.txt "
    		  . "-e $errorPath.err.txt "
              . "-p $sgePrioritie "
    		  . File::Spec->catfile('@CASAVA_FULL_LIBEXECDIR@', 'TaskManager', 'taskAgent.pl') . ' '
    		  . "--type=sge $agentPrams";

    	printLog( $cmd . "\n", 6 );
    	printLog( "Submit start\n", 6 );
    	my $sgeJobId = submit2SGE($cmd);
    	printLog( "Submit end\n", 6 );
        return $sgeJobId;
    }
    return -1;
}

=pod

=head1 The procedure submits SGE job to SGE.

=over 4

=item submit2SGE($command)

The procedure submits SGE job to SGE and returns the job id assigned to the job.

Parameters:
    $command - sge formated command to be submited to sge

Returns:
    sge job id
=back

=cut

sub submit2SGE($) {
	croak "ERROR: submit2SGE\n" unless ( @_ == 1 );
	my ($command) = @_;
	my $jobid;
	open( SGE, "$command 2> /dev/null |" )
	  || errorExit "Could not submit $command to SGE $! \n";
	my $respond = '';
	while (<SGE>) {
		print $_;
		$respond .= $_;
		if ( $_ =~ /your job-array (\d+)/i ) {
			$jobid = $1;
		}
		if ( $_ =~ /your job (\d+)/i ) {
			$jobid = $1;
		}
	}
	if ( !defined($jobid) ) {
		system("qstat -f");
		system("qstat -g c");
		errorExit "ERROR: Failed to qsub $command $!. Response = [$respond]\n";
	}
	else {
		logInfo( "Submitted with $jobid $command\n", 6 );
	}
    close (SGE);
	#system("sleep 1s");
	return $jobid;
}

=pod

=head1 verify if SGE queue is valid

=over 4

=item isSgeQueueOk($queue)

The procedure runs the qstat q <queue> and checks for the result code.

Returns:
    1 if test qas successful
=back

=cut

sub isSgeQueueOk($) {
    my ($queue) = @_;
    return testIfSupports("qstat -q $queue");
}
1;    # says use was ok
__END__

