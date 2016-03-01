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

Summary.pm

=cut

=head1 DESCRIPTION

Class to encapsulate a Summary.xml file.

=cut

=head1 AUTHOR

Richard Shaw

=cut



use warnings;
use strict;

use XML::Simple;
use Data::Dumper;

#------------------------------------------------------------------------------

use lib '@CASAVA_FULL_PERL_LIBDIR@';

use Casava::PostAlignment::QC::Lane;
use Casava::PostAlignment::QC::Rule_Set;
use Casava::PostAlignment::QC::Rule;

#------------------------------------------------------------------------------

package Casava::PostAlignment::QC::Summary;

#------------------------------------------------------------------------------

my $lane_str = 'Lane';

my $summary_path_key = 'SummaryPath';
my $lanes_key = 'Lanes';
my $yield_key = 'yield';

my $lane_parameter_summary_tag = 'LaneParameterSummary';
my $lane_tag = 'Lane';
my $lane_number_tag = 'laneNumber';
my $read_tag = 'Read';
my $read_number_tag = 'readNumber';

#------------------------------------------------------------------------------

sub make_lane_key($)
{
    my $lane_num = shift;
    return join(':', $lane_str, $lane_num);
}

#------------------------------------------------------------------------------

sub new($;$)
{
    my ($class, $summary_path) = @_;
    my $self = {};
    bless $self, $class;

    $self->{$summary_path_key} = $summary_path; # For error messages
    $self->{$lanes_key} = {};

    if (! -f $summary_path) {
        print "File $summary_path not found\n";
        return 0;
    }

    # print "Loading Summary from $summary_path\n";

    my $summary_xml = new XML::Simple(suppressempty => '',
                                      XMLDecl => 1);

    my $summary_hash_ref
      = $summary_xml->XMLin($summary_path,
                            forcearray => [qw(Lane Read Tile)]);

    # print "Loaded Summary from $summary_path\n";

    if (!defined($summary_hash_ref)) {
        return 0;
    }

    $self->process_parameter_summary($summary_hash_ref,
                                     $lane_parameter_summary_tag);

    my @lane_table_tags = ('LaneResultsSummary', 'ExpandedLaneSummary');

    foreach my $lane_table_tag (@lane_table_tags) {
        $self->process_lane_table($summary_hash_ref, $lane_table_tag);
    }

    my @tile_table_tags = ('TileResultsByLane', 'TileErrorsByLane');

    foreach my $tile_table_tag (@tile_table_tags) {
        $self->process_tile_table($summary_hash_ref, $tile_table_tag);
    }

    return $self;
}

#------------------------------------------------------------------------------

sub process_parameter_summary($;$;$)
{
    my ($self, $summary_hash_ref, $table_tag) = @_;

    my $table_ref = $summary_hash_ref->{$table_tag};

    my $lanes_ref = $table_ref->{$lane_tag};

    foreach my $lane_ref (@$lanes_ref) {
        my $lane_number = $lane_ref->{$lane_number_tag};

        $self->{$lanes_key}->{$lane_number} = 1;

        my $lane_key = make_lane_key($lane_number);

        if (!defined($self->{$lane_key})) {
            $self->{$lane_key} = Casava::PostAlignment::QC::Lane->new($lane_number);
        }

        $self->{$lane_key}->add_params($table_tag, $lane_ref,
                                       $lane_number,
                                       $self->{$summary_path_key});
    }
}

#------------------------------------------------------------------------------
# <Read><Lane><stat>

sub process_lane_table($;$;$)
{
    my ($self, $summary_hash_ref, $lane_table_tag) = @_;

    my $lane_table_ref = $summary_hash_ref->{$lane_table_tag};

    my $lane_table_reads_ref = $lane_table_ref->{$read_tag};

    foreach my $lane_table_read_ref (@$lane_table_reads_ref) {
        my $read_number = $lane_table_read_ref->{$read_number_tag};
        my $lane_table_read_lanes_ref = $lane_table_read_ref->{$lane_tag};

        foreach my $read_lane_ref (@$lane_table_read_lanes_ref) {
            my $lane_number = $read_lane_ref->{$lane_number_tag};

            $self->{$lanes_key}->{$lane_number} = 1;
            my $lane_key = make_lane_key($lane_number);

            if (!defined($self->{$lane_key})) {
                $self->{$lane_key} = Casava::PostAlignment::QC::Lane->new($lane_number);
            }

            $self->{$lane_key}->add_read_lane_stats($lane_table_tag,
                                                    $read_number,
                                                    $read_lane_ref,
                                                    $self->{$summary_path_key});
        }
    }
}

#------------------------------------------------------------------------------
# <Lane><Read><Tile><stat>

sub process_tile_table($;$;$)
{
    my ($self, $summary_hash_ref, $tile_table_tag) = @_;

    my $tile_table_ref = $summary_hash_ref->{$tile_table_tag};

    my $tile_table_lanes_ref = $tile_table_ref->{$lane_tag};

    foreach my $tile_table_lane_ref (@$tile_table_lanes_ref) {
        my $lane_number = $tile_table_lane_ref->{$lane_number_tag};

        $self->{$lanes_key}->{$lane_number} = 1;
        my $lane_key = make_lane_key($lane_number);

        if (!defined($self->{$lane_key})) {
            $self->{$lane_key} = Casava::PostAlignment::QC::Lane->new($lane_number);
        }

        $self->{$lane_key}->add_tile_stats($tile_table_tag,
                                           $tile_table_lane_ref,
                                           $self->{$summary_path_key});
    }
}

#------------------------------------------------------------------------------

sub apply_rule($;$;$)
{
    my ($self, $rule, $qc_log_hndl) = @_;

    # DEBUG
    # print("Summary : apply_rule ", $rule->get_xpath(), "\n");

    my $match_lane_num = 0;
    my $any_lane_num = 0;

    $rule->reduce_xpath($lane_str, \$any_lane_num, \$match_lane_num);

    foreach my $lane_num (sort {$a <=> $b} keys %{$self->{$lanes_key}}) {
        if ($any_lane_num || ($lane_num == $match_lane_num)) {

            # DEBUG
            # print("Summary : apply_rule (Lane ", $lane_num, ") : ",
            #       $rule->get_xpath(), "\n");

            $self->{make_lane_key($lane_num)}->apply_rule($rule, $qc_log_hndl);
        }
    }

    $rule->undo_reduce_xpath();
}

#------------------------------------------------------------------------------

sub update_stats($)
{
    my ($self) = @_;
    my $summary_path = $self->{$summary_path_key};

    foreach my $lane_num (keys %{$self->{$lanes_key}}) {
        # DEBUG
        # print("update_stats : lane $lane_num\n");

        $self->{make_lane_key($lane_num)}->update_stats($summary_path);
    }
}

#------------------------------------------------------------------------------

sub apply_rule_set($;$;$)
{
    my ($self, $rule_set, $qc_log_hndl) = @_;

    print $qc_log_hndl 'Filtering : ', $self->{$summary_path_key}, "\n";

    my @rule_xpaths = $rule_set->xpaths();

    foreach my $rule_xpath (@rule_xpaths) {
        my $rule = $rule_set->rule($rule_xpath);
        $self->apply_rule($rule, $qc_log_hndl);
    }

    $self->update_stats();
}

#------------------------------------------------------------------------------

sub lane_is_valid($;$;$;$)
{
    my ($self, $lane_num, $lane_yield_ref, $bad_tile_nums_ref) = @_;

    my $lane_ref = $self->{make_lane_key($lane_num)};

    if (!defined($lane_ref)) {
        die("Attempt to access unknown lane $lane_num.");
    }

    $$lane_yield_ref = $lane_ref->get_yield();
    @{$bad_tile_nums_ref} = $lane_ref->get_bad_tile_nums();

    return $lane_ref->is_valid();
}

#------------------------------------------------------------------------------

sub dump_bad_tile_nums($)
{
    my $self = shift;

    print "BAD TILES\n";

    foreach my $lane_num (sort {$a <=> $b} keys %{$self->{$lanes_key}}) {
        print "Lane ", $lane_num, "\n";

        my $lane_yield = 0;
        my @bad_tile_nums = ();

        if ($self->lane_is_valid($lane_num, \$lane_yield, \@bad_tile_nums)) {
            print join(',', @bad_tile_nums), "\n";
        } else {
            print "LANE REJECTED\n";
        }
    }
}

#------------------------------------------------------------------------------

sub dump_yields($)
{
    my $self = shift;
    print "Yields\n";

    foreach my $lane_num (sort {$a <=> $b} keys %{$self->{$lanes_key}}) {
        print "Lane ", $lane_num, "\n";

        my $lane_ref = $self->{make_lane_key($lane_num)};

        if ($lane_ref->is_valid()) {
            print $lane_ref->get_yield(), "\n";
        } else {
            print "LANE REJECTED\n";
        }
    }
}

#------------------------------------------------------------------------------

sub dump($)
{
    my $self = shift;

    foreach my $lane_num (sort {$a <=> $b} keys %{$self->{$lanes_key}}) {
        print "Lane ", $lane_num, "\n";

        $self->{make_lane_key($lane_num)}->dump();

        print "-----\n";
    }
}

#------------------------------------------------------------------------------

return 1;

#------------------------------------------------------------------------------
