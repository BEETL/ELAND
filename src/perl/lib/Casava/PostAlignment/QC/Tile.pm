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

Tile.pm

=cut

=head1 DESCRIPTION

Class to encapsulate a tile component of a Summary.xml file.

=cut

=head1 AUTHOR

Richard Shaw

=cut


use warnings;
use strict;

#------------------------------------------------------------------------------

use lib '@CASAVA_FULL_PERL_LIBDIR@';

use Casava::PostAlignment::QC::Rule;

#------------------------------------------------------------------------------

package Casava::PostAlignment::QC::Tile;

#------------------------------------------------------------------------------

my $read_str = 'Read';

my $tile_results_by_lane_tag = 'TileResultsByLane';
my $cluster_count_pf_tag = 'clusterCountPF';
my $tile_number_tag = 'tileNumber';

my $tile_number_key = 'tile_number';
my $stats_key = 'stats';
my $table_names_key = 'table_names';
my $read_nums_key = 'read_nums';
my $valid_key = 'valid';

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

sub new($;$)
{
    my ($class, $tile_num) = @_;
    my $self = {};

    # Keep the tile number for error reporting only.
    $self ->{$tile_number_key} = $tile_num;
    $self->{$read_nums_key} = {};
    $self->{$stats_key} = {};
    $self->{$valid_key} = 1;

    bless $self, $class;
    return $self;
}

#------------------------------------------------------------------------------

sub add_read_stats {
    my ($self, $table_name, $read_num, $tile_xml_ref) = @_;

    $self->{$read_nums_key}->{$read_num} = 1;

    # print("Tile ", $self ->{$tile_number_key}, " : adding table $table_name"
    #       . " read $read_num\n";

    if (!defined($self->{$table_name})) {
        $self ->{$table_name} = {};
    }

    my $read_key = make_read_key($read_num);

    if (defined($self->{$table_name}->{$read_key})) {
        die("Attempt to add ", $table_name,
            " stats for Read ", $read_num,
            " to tile ", $self ->{$tile_number_key},
            " multiple times.");
    }

    my @stat_names = keys %$tile_xml_ref;

    foreach my $stat_name (@stat_names) {
        next if ($stat_name eq $tile_number_tag);

        $self->{$stats_key}{make_table_read_stat_key($table_name,
                                                     $read_num,
                                                     $stat_name)}
          = $tile_xml_ref->{$stat_name};
    }
}

#------------------------------------------------------------------------------

sub set_val() {
    my ($self, $key, $val) = @_;

    $self->{$key} = $val;
}

#------------------------------------------------------------------------------

sub get_val {
    my ($self, $key) = @_;

    return $self->{$key};
}

#------------------------------------------------------------------------------

sub apply_rule($;$;$;$)
{
    my ($self, $rule, $lane_num, $qc_log_hndl) = @_;

    # DEBUG
    # print("Tile : apply_rule ", $rule->get_xpath(), "\n");

    foreach my $read_num (keys %{$self->{$read_nums_key}}) {
        my $rule_xpath = $rule->localised_xpath($read_str, $read_num);

        if ($rule_xpath) {
            my $stat_val = $self->{$stats_key}->{$rule_xpath};

            if (!$rule->val_ok($stat_val)) {
                $self->{$valid_key} = 0;

                print($qc_log_hndl 'Lane ', $lane_num,
                      ' : Tile ',  $self ->{$tile_number_key},
                      ' : FAILED rule ', $rule->to_str(),
                      ' : Read ', $read_num,
                      ' val was ', $stat_val, "\n");
            }
        }
    }
}

#------------------------------------------------------------------------------

sub is_valid($)
{
    my $self = shift;
    return $self->{$valid_key};
}

#------------------------------------------------------------------------------
# For yield

sub num_pf_clusters($)
{
    my $self = shift;
    my $read_num = 1;

    return $self->{$stats_key}
      ->{make_table_read_stat_key($tile_results_by_lane_tag, 
                                  $read_num,
                                  $cluster_count_pf_tag)};
}

#------------------------------------------------------------------------------

sub dump($)
{
    my $self = shift;

    print "Tile ", $self ->{$tile_number_key}, "\n";

    foreach my $stat_path (sort keys %{$self->{$stats_key}}) {
        print $stat_path, " : ", $self->{$stats_key}->{$stat_path}, "\n";
    }

    print "---\n";
}

#------------------------------------------------------------------------------

return 1;

#------------------------------------------------------------------------------
