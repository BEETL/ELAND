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

Lane.pm

=cut

=head1 DESCRIPTION

Class to encapsulate a lane component of a Summary.xml file.

=cut

=head1 AUTHOR

Richard Shaw

=cut


use warnings;
use strict;

use Data::Dumper;

#------------------------------------------------------------------------------

use lib '@CASAVA_FULL_PERL_LIBDIR@';

use Casava::Common::Log;
use Casava::PostAlignment::QC::Tile;
use Casava::PostAlignment::QC::Rule;

#------------------------------------------------------------------------------

package Casava::PostAlignment::QC::Lane;

#------------------------------------------------------------------------------

my $lane_number_tag = 'laneNumber';
my $read_tag = 'Read';
my $read_number_tag = 'readNumber';
my $tile_tag = 'Tile';
my $tile_number_tag = 'tileNumber';

my $lengths_list_tag = 'lengthsList';

my $read_str = 'Read';
my $tile_str = 'Tile';

my $lane_number_key = 'lane_number';
my $tiles_key = 'Tiles';
my $lengths_key = 'Lengths';

my $read_nums_key = 'read_nums';
my $stats_key = 'stats';
my $tiles_checked_key = 'tiles_checked';
my $valid_key = 'valid';
my $bad_tile_nums_key = 'bad_tile_nums';
my $yield_key = 'yield';

my $mean_str = 'mean';
my $stdev_str = 'stdev';
my $sumsq_str = 'sumsq';
my @stat_field_names = ($mean_str, $stdev_str, $sumsq_str);

#------------------------------------------------------------------------------

sub make_tile_key($)
{
    my $tile_num = shift;
    return sprintf("%s%04d", $tile_str, $tile_num);
}

#------------------------------------------------------------------------------

sub make_read_key($)
{
    my ($read_num) = @_;
    return $read_str . $read_num;
}

#------------------------------------------------------------------------------

sub make_table_read_stat_key($;$;$)
{
    my ($table_name, $read_num, $stat_name) = @_;
    return '/' . join('/', $table_name, make_read_key($read_num), $stat_name);
}

#------------------------------------------------------------------------------

sub make_table_read_stat_field_key($;$;$;$)
{
    my ($table_name, $read_num, $stat_name, $field_name) = @_;
    return join('/',
                make_table_read_stat_key($table_name, $read_num, $stat_name),
                $field_name);
}

#------------------------------------------------------------------------------

sub new($;$)
{
    my ($class, $lane_number) = @_;
    my $self = {};
    bless $self, $class;

    $self->{$lane_number_key} = $lane_number;
    $self->{$read_nums_key} = {};
    $self->{$stats_key} = {};
    $self->{$valid_key} = 1;
    $self->{$tiles_checked_key} = 0;

    return $self;
}

#------------------------------------------------------------------------------
# Want just the lengths - for yield calculation.

sub add_params($;$;$;$;$)
{
    my ($self, $table_tag, $lane_ref, $lane_number, $summary_path) = @_;

    my $lengths_str = $lane_ref->{$lengths_list_tag};

    if (!defined($lengths_str)) {
        my $msg = "Failed to find $lengths_str for lane $lane_number "
          . "under $table_tag in $summary_path.";
        Casava::Common::Log::printLog($msg, 0);
        $lengths_str = '';
    }

    my @lengths = split(', ', $lengths_str);

    @{$self->{$lengths_key}} = @lengths;
}

#------------------------------------------------------------------------------

sub add_read_lane_stats($;$;$;$;$)
{
    my ($self, $table_name, $read_number, $read_lane_ref, $summary_path) = @_;
    $self->{$read_nums_key}->{$read_number} = 1;

    my @stat_names = keys %$read_lane_ref;

    foreach my $stat_name (@stat_names) {
        next if ($stat_name eq $lane_number_tag);

        my $stat_fields_ref = $read_lane_ref->{$stat_name};

        # Some statistics have just a value - not mean, stdev, etc.
        if (not ref($stat_fields_ref)) {
            $self->{$stats_key}{make_table_read_stat_key($table_name,
                                                         $read_number,
                                                         $stat_name)}
              = $stat_fields_ref;
            next;
        }

        foreach my $field_name (@stat_field_names) {
            my $field_val = $stat_fields_ref->{$field_name};

            $self->{$stats_key}{make_table_read_stat_field_key($table_name,
                                                               $read_number,
                                                               $stat_name,
                                                               $field_name)}
              = $field_val;
        }
    }
}

#------------------------------------------------------------------------------

sub add_tile_stats($;$;$;$)
{
    my ($self, $tile_table_tag, $tile_table_lane_ref, $summary_path) = @_;

    # Keep track of the tiles found.
    $self->{$tiles_key} = {};

    my $read_refs_ref = $tile_table_lane_ref->{$read_tag};

    if (!defined($read_refs_ref)) {
        my $msg = 'No tile results found for lane '
          . $self->{$lane_number_key} . ' in ' . $summary_path . "\n";
        Casava::Common::Log::printLog($msg, 0);

        # DEBUG
        # print("tile lane table ref\n",
        #      Data::Dumper->Dumper($tile_table_lane_ref), "\n");
        return;
    }

    my @reads_ref = @{$read_refs_ref};

    foreach my $read_ref (@reads_ref) {
        my $read_number = $read_ref->{$read_number_tag};
        $self->{$read_nums_key}->{$read_number} = 1;
        my $tile_set_ref = $read_ref->{$tile_tag};

        my @tiles_ref = @$tile_set_ref;

        foreach my $tile_xml_ref (@tiles_ref) {
            my $tile_number = $tile_xml_ref->{$tile_number_tag};
            my $tile_key = make_tile_key($tile_number);

            # Update of list of tiles without duplication.
            $self->{$tiles_key}->{$tile_number} = 1;

            # Create an instance of the tile in the first Read it is found in.
            if (!defined($self->{make_tile_key($tile_number)})) {
                $self->{$tile_key} = new Casava::PostAlignment::QC::Tile($tile_number);
            }

            $self->{$tile_key}->add_read_stats($tile_table_tag,
                                               $read_number,
                                               $tile_xml_ref);
        }
    }
}

#------------------------------------------------------------------------------

sub apply_lane_rule($;$;$)
{
    my ($self, $rule, $qc_log_hndl) = @_;

    # DEBUG
    # print("Lane : apply_lane_rule ", $rule->get_xpath(), "\n");

    foreach my $read_num (keys %{$self->{$read_nums_key}}) {
        my $rule_xpath = $rule->localised_xpath($read_str, $read_num);

        if ($rule_xpath) {
            my $stat_val = $self->{$stats_key}->{$rule_xpath};

            if (!$rule->val_ok($stat_val)) {
                $self->{$valid_key} = 0;

                print($qc_log_hndl 'Lane ',  $self ->{$lane_number_key},
                      ' (Read ', $read_num, ')',
                      ' : FAILED rule ', $rule->to_str(),
                      ' : Read ', $read_num,
                      ' val was ', $stat_val, "\n");
            }
        }
    }
}

#------------------------------------------------------------------------------

sub apply_rule($;$;$)
{
    my ($self, $rule, $qc_log_hndl) = @_;

    # DEBUG
    # print("Lane : apply_rule ", $rule->get_xpath(), "\n");

    my $match_tile_num = 0;
    my $any_tile_num = 0;

    my $did_reduce
      = $rule->reduce_xpath($tile_str, \$any_tile_num, \$match_tile_num);

    if ($did_reduce) {
        my $lane_num = $self->{$lane_number_key};

        # Tile stat rule
        foreach my $tile_num (sort {$a <=> $b} keys %{$self->{$tiles_key}}) {
            if ($any_tile_num || ($tile_num == $match_tile_num)) {
                $self->{make_tile_key($tile_num)}->apply_rule($rule,
                                                              $lane_num,
                                                              $qc_log_hndl);
            }
        }
    } else {
        # Lane stat rule
        $self->apply_lane_rule($rule, $qc_log_hndl);
    }

    $rule->undo_reduce_xpath();
}

#------------------------------------------------------------------------------
# Currently derives only bad tiles list, valid flag and yield.

sub update_stats($;$)
{
    my ($self, $summary_path) = @_;

    my @bad_tile_nums = ();
    my $num_good_tiles = 0;
    my $num_valid_tile_pf_clusters = 0;

    foreach my $tile_num (keys %{$self->{$tiles_key}}) {
        my $tile_ref = $self->{make_tile_key($tile_num)};

        if ($tile_ref->is_valid()) {
            ++$num_good_tiles;
            $num_valid_tile_pf_clusters += $tile_ref->num_pf_clusters();
        } else {
            push (@bad_tile_nums, $tile_num);
        }
    }


    # Calculate lane yield as :-
    # sum<reads>(length used) * sum<over valid tiles>(num PF clusters).
    my @lengths = @{$self->{$lengths_key}};
    my $total_length = 0;
    my $read_num = 0;

    foreach my $length (@lengths) {
        ++$read_num;

        if (!($length =~ /\d+/)) {
            my $lane_num = $self->{$lane_number_key};
            my $msg = "Read $read_num length $length for lane $lane_num in "
              . $summary_path . "\n";
            Casava::Common::Log::printLog($msg, 0);
            $length = 0;
        }

        $total_length += $length;
    }

    my $yield = $num_valid_tile_pf_clusters * $total_length;


    # Update members.

    $self->{$tiles_checked_key} = 1;

    if ($num_good_tiles == 0) {
        $self->{$valid_key} = 0;
    }

    @{$self->{$bad_tile_nums_key}} = @bad_tile_nums;
    $self->{$yield_key} = $yield;
}

#------------------------------------------------------------------------------

sub get_bad_tile_nums($)
{
    my ($self) = @_;

    if (!$self->{$tiles_checked_key}) {
        die("Lane::is_valid called before Lane::update_stats");
    }

    return @{$self->{$bad_tile_nums_key}};
}

#------------------------------------------------------------------------------

sub is_valid($)
{
    my $self = shift;

    if (!$self->{$tiles_checked_key}) {
        die("Lane::is_valid called before Lane::update_stats");
    }

    # FIXME : Add Lane stat checks?

    return $self->{$valid_key};
}

#------------------------------------------------------------------------------

sub get_yield($)
{
    my ($self) = @_;
    return $self->{$yield_key};
}

#------------------------------------------------------------------------------

sub dump($)
{
    my $self = shift;

    print "Lane ", $self->{$lane_number_key}, "\n";

    foreach my $stat_path (sort keys %{$self->{$stats_key}}) {
        print $stat_path, " : ", $self->{$stats_key}->{$stat_path}, "\n";
    }

    print "---\n";

    foreach my $tile_num (sort {$a <=> $b} keys %{$self->{$tiles_key}}) {
        print "Tile ", $tile_num, "\n";
        my $tile_key = make_tile_key($tile_num);
        $self->{$tile_key}->dump();
    }
}

#------------------------------------------------------------------------------

return 1;

#------------------------------------------------------------------------------
