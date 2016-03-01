#!/usr/bin/env perl

=head1 LICENSE

Copyright (c) 2009 Illumina, Inc.


This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

This file is part of the Consensus Assessment of Sequence And VAriation

CASAVA) software package.

=cut

=head1 NAME

Builder.pm

=cut

=head1 DESCRIPTION

Class to produce actual_build.xml from requested_build.xml and run configureBuild.pl.


=cut

=head1 AUTHOR

Richard Shaw

=cut

use warnings;
use strict;

#------------------------------------------------------------------------------

use File::Spec;

use lib '@CASAVA_FULL_PERL_LIBDIR@';

use Casava::Common::Log;
use Casava::PostAlignment::QC::Build_Spec;
use Casava::PostAlignment::QC::Rule_Set;
use Casava::PostAlignment::QC::Summary;

#------------------------------------------------------------------------------

package Casava::PostAlignment::QC::Builder;

#------------------------------------------------------------------------------

my $build_dir_path_key = 'build_dir_path';
my $requested_build_path_key = 'requested_build_path';
my $actual_build_path_key = 'actual_build_path';
my $run_cmd_key = 'run_cmd';

#------------------------------------------------------------------------------

sub new($;$;$;$;$)
{
    my ($class, $build_dir_path, $requested_build_path, $actual_build_path,
        $run_script_path) = @_;
    my $self = {};
    bless $self, $class;

    $self->{$build_dir_path_key} = $build_dir_path;
    $self->{$requested_build_path_key} = $requested_build_path;
    $self->{$actual_build_path_key} = $actual_build_path;
    $self->{$run_cmd_key} = join(' ', $run_script_path,
                                 "--projectDir=$build_dir_path");

    return $self;
}

#------------------------------------------------------------------------------

sub add_run_to_run_cmd($;$;$;$)
{
    my ($self, $run_num, $export_dir_path, $lanes_used_ref) = @_;

    my $run_id = "Run${run_num}";
    my $lanes_used_list = join(',', @{$lanes_used_ref});

    $self->{$run_cmd_key} = join(' ', $self->{$run_cmd_key},
                                 "--runId=${run_id}",
                                 "--exportDir=${export_dir_path}",
                                 "--lanes=${lanes_used_list}");
}

#------------------------------------------------------------------------------

sub build($)
{
    my ($self) = @_;

    my $actual_build_path = $self->{$actual_build_path_key};

    if (-f $actual_build_path) {
        my $msg = "Nothing to do - actual build file exists: "
          . $actual_build_path . "\n";
        Casava::Common::Log::printLog($msg, 0);
        return;
    }

    my $build_spec
      = Casava::PostAlignment::QC::Build_Spec->new($self->{$requested_build_path_key});

    my $tile_rule_set_path = $build_spec->get_tile_selection_config_path();
    my $tile_rule_set = Casava::PostAlignment::QC::Rule_Set->new($tile_rule_set_path);

    my $min_yield = $build_spec->get_min_yield();

    $self->{$run_cmd_key} = join(' ', $self->{$run_cmd_key},
                                 $build_spec->get_run_script_fixed_args());

    my $num_lanes_used = 0;
    my $total_yield = 0;
    my $run_num = 0;

    my $qc_log_path
          = File::Spec->catfile($self->{$build_dir_path_key}, 'QC.log');
    open(QC_LOG, ">$qc_log_path")
      or die("Failed to open QC log file $qc_log_path for writing.");

    my $summary_filename = 'Summary.xml';

    foreach my $export_dir ($build_spec->get_export_dirs()) {
        ++$run_num;
        my $export_dir_spec = $build_spec->get_export_dir_spec($export_dir);

        my $summary_path = File::Spec->catfile($export_dir, $summary_filename);

        if (! -e $summary_path) {
            my $msg = "Did not find $summary_filename "
              . "under $export_dir : skipping\n";
            Casava::Common::Log::printLog($msg, 0);
            next;
        }

        my $summary = Casava::PostAlignment::QC::Summary->new($summary_path);

        if (!$summary) {
            my $msg = "Failed to parse $summary_path : skipping $export_dir"; 
            Casava::Common::Log::printLog($msg, 0);
            next;
        }

        $summary->apply_rule_set($tile_rule_set, *QC_LOG);
        my $flowcell_yield_contrib = 0;

        my @lane_nums_used = ();

        Casava::Common::Log::printLog("Summary : $summary_path\n", 0);

        foreach my $lane_num (sort {$a <=> $b} $export_dir_spec->get_lane_nums()) {
            my @bad_tile_nums = ();
            my $lane_yield = 0;

            if ($summary->lane_is_valid($lane_num,
                                        \$lane_yield,
                                        \@bad_tile_nums)) {
                ++$num_lanes_used;
                $flowcell_yield_contrib += $lane_yield;
                $total_yield += $lane_yield;

                my $msg = "Lane $lane_num : Filtered yield $lane_yield\n";
                Casava::Common::Log::printLog($msg, 0);

                my $lane_spec = $export_dir_spec->get_lane_spec($lane_num);
                $lane_spec->add_bad_tiles(\@bad_tile_nums);

                push(@lane_nums_used, $lane_num);
            } else {
                my $msg = "Lane $lane_num : excluded (yield 0)\n";
                Casava::Common::Log::printLog($msg, 0);
                $export_dir_spec->delete_lane($lane_num);
            }
        }

        if (@lane_nums_used) {
            $self->add_run_to_run_cmd($run_num, $export_dir,
                                      \@lane_nums_used);
        }

        Casava::Common::Log::printLog("Flowcell yield contrib : "
                                        . $flowcell_yield_contrib . "\n", 0);
    }

    close(QC_LOG);

    Casava::Common::Log::printLog("Num lanes used    : $num_lanes_used\n",
				    0);
    Casava::Common::Log::printLog("Required yield    : $min_yield\n", 0);
    Casava::Common::Log::printLog("Accumulated yield : $total_yield\n", 0);

    if ($total_yield >= $min_yield) {
        $build_spec->write($actual_build_path);
        Casava::Common::Log::printLog("Wrote actual build file : "
                                        . $actual_build_path . "\n", 0);

        my $run_cmd = $self->{$run_cmd_key};

        # Add passing of actual build spec.
        $run_cmd = join(' ', $run_cmd,
                        "--buildFile=${actual_build_path}");

        # Add redirection of output.
        my $run_script_output_path
          = File::Spec->catfile($self->{$build_dir_path_key}, 'run.log');
        $run_cmd = join(' ', $run_cmd,
                        "> $run_script_output_path",
                        "2> $run_script_output_path");

        # Fork off configureBuild.pl, so caller can continue.
        Casava::Common::Log::printLog("Forking run script : $run_cmd\n", 0);
        system("${run_cmd} &");
    } else {
        Casava::Common::Log::printLog("Insufficient yield to proceed\n", 0);
    }
}

#------------------------------------------------------------------------------

return 1;

#------------------------------------------------------------------------------
