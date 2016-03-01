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

Export_Dir_Spec.pm

=cut

=head1 DESCRIPTION

Class to encapsulate the export directory component of a build specification.

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

use Casava::PostAlignment::QC::Lane_Spec;

#------------------------------------------------------------------------------

package Casava::PostAlignment::QC::Export_Dir_Spec;

#------------------------------------------------------------------------------

my $lanes_tag = 'Lanes';
my $lane_tag = 'Lane';

my $xml_key = 'XML';
my $lane_set_ref_key = 'lane_set_ref';
my $lanes_key = 'lanes';

#------------------------------------------------------------------------------

sub new($;$)
{
    my ($class, $export_dir_xml_ref) = @_;
    my $self = {};
    bless $self, $class;

    $self->{$xml_key} = $export_dir_xml_ref;
    $self->parse_spec();

    # DEBUG
    # print Data::Dumper->Dumper($export_dir_xml_ref);

    return $self;
}

#------------------------------------------------------------------------------

sub parse_spec($)
{
    my ($self) = @_;
    $self->{$lanes_key} = {};

    my $lanes_xml_ref = $self->{$xml_key}->{$lanes_tag};

    if (!defined($lanes_xml_ref)) {
        die("Failed to find $lanes_tag in export dir XML");
    }

    my $lane_set_ref = $lanes_xml_ref->{$lane_tag};

    # Need to store this in case lanes have to be deleted.
    $self->{$lane_set_ref_key} = $lane_set_ref;

    if (!defined($lane_set_ref)) {
        warn("Failed to find $lane_tag under $lanes_tag in export dir XML");
        return;
    }

    foreach my $lane_number (keys %$lane_set_ref) {
        # DEBUG
        # print "Lane : ", $lane_number, "\n";

        my $lane_xml_ref = $lane_set_ref->{$lane_number};
        $self->{$lanes_key}->{$lane_number}
          = new Casava::PostAlignment::QC::Lane_Spec($lane_xml_ref);
    }
}

#------------------------------------------------------------------------------

sub get_lane_nums($)
{
    my ($self) = @_;
    return keys %{$self->{$lanes_key}};
}

#------------------------------------------------------------------------------

sub get_lane_spec($;$)
{
    my ($self, $lane_num) = @_;
    return $self->{$lanes_key}->{$lane_num};
}

#------------------------------------------------------------------------------

sub delete_lane($;$)
{
    my ($self, $lane_num) = @_;

    # Delete the Lane object.
    delete($self->{$lanes_key}->{$lane_num});

    # Delete the XML for the lane.
    delete($self->{$lane_set_ref_key}->{$lane_num});
}

#------------------------------------------------------------------------------

sub dump($)
{
    my ($self) = @_;
    print Data::Dumper->Dumper($self->{$xml_key});
}

#------------------------------------------------------------------------------

