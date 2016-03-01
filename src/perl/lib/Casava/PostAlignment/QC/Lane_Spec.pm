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

Lane_Spec.pm

=cut

=head1 DESCRIPTION

Class to encapsulate the lane component of a build specification.

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

#------------------------------------------------------------------------------

package Casava::PostAlignment::QC::Lane_Spec;

#------------------------------------------------------------------------------

my $number_tag = 'Number';
my $bad_tiles_tag = 'BadTiles';

my $xml_key = 'XML';
my $bad_tiles_key = 'bad_tiles';

#------------------------------------------------------------------------------

sub new($;$)
{
    my ($class, $lane_xml_ref) = @_;
    my $self = {};
    bless $self, $class;

    $self->{$xml_key} = $lane_xml_ref;
    $self->parse_spec();

    # DEBUG
    # print "Lane :-\n", Data::Dumper->Dumper($lane_xml_ref);

    return $self;
}

#------------------------------------------------------------------------------

sub parse_spec($)
{
    my ($self) = @_;
    $self->{$bad_tiles_key} = ();

    my $bad_tiles_list = $self->{$xml_key}->{$bad_tiles_tag};

    if (defined($bad_tiles_list)) {
        $bad_tiles_list =~ s/\s//g;
        @{$self->{$bad_tiles_key}} = split(',', $bad_tiles_list);
    }
}

#------------------------------------------------------------------------------

sub add_bad_tiles($;$)
{
    my ($self, $new_bad_tiles_arr_ref) = @_;

    my $old_bad_tiles_arr_ref = $self->{$bad_tiles_key};

    my @old_bad_tiles = (defined($old_bad_tiles_arr_ref)
                         ? @{$old_bad_tiles_arr_ref}
                         : ());
    my @new_bad_tiles = sort {$a <=> $b} @{$new_bad_tiles_arr_ref};

    my @bad_tiles = ();

    if (!@old_bad_tiles) {
        # No previously specified bad tiles.
        @bad_tiles = @new_bad_tiles;
    } else {
        my %is_bad_tile = map {$_ => 1} (@old_bad_tiles, @new_bad_tiles);
        @bad_tiles = sort {$a <=> $b} keys %is_bad_tile;
    }

    @{$self->{$bad_tiles_key}} = @bad_tiles;
    $self->{$xml_key}->{$bad_tiles_tag} = join(',', @bad_tiles);
}

#------------------------------------------------------------------------------

sub get_bad_tiles($)
{
    my ($self) = @_;
    return @{$self->{$bad_tiles_key}};
}

#------------------------------------------------------------------------------

return 1;

#------------------------------------------------------------------------------
