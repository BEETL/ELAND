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

Build_Spec.pm

=cut

=head1 DESCRIPTION

Class to encapsulate a build specification (before or after filtering)

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

use Casava::PostAlignment::QC::Export_Dir_Spec;

#------------------------------------------------------------------------------

package Casava::PostAlignment::QC::Build_Spec;

#------------------------------------------------------------------------------

my $tile_selection_config_tag = 'TileSelectionConfig';
my $min_yield_gb_tag = 'MinYieldForBuildGb';
my $run_script_fixed_args_tag = 'RunScriptFixedArgs';

# Min yield Gb is optional. If not supplied then build unconditionally.
my @required_val_tags = ($tile_selection_config_tag, $run_script_fixed_args_tag);

my $alignments_tag = 'Alignments';
my $alignment_tag = 'Alignment';
my $export_dir_tag = 'ExportDir';
my $lane_tag = 'Lane';
my $number_tag = 'Number';

my $xml_key = 'XML';
my $export_dirs_key = 'export_dirs';

#------------------------------------------------------------------------------

sub new($;$)
{
    my ($class, $build_spec_path) = @_;
    my $self = {};
    bless $self, $class;

    if (! -f $build_spec_path) {
        die("Request file $build_spec_path not found\n");
    }

    my $build_xml = new XML::Simple(suppressempty => '',
                                    XMLDecl => 1);

    my $build_hash = undef;

    eval {
        $build_hash = $build_xml->XMLin($build_spec_path,
                                        forcearray => [$alignment_tag,
                                                       $lane_tag],
                                        keyattr => {
                                                    $alignment_tag => $export_dir_tag,
                                                    $lane_tag => $number_tag});
    };

    if ($@) {
        die("Caught XML::Simple throw for $build_spec_path : $@\n");
    }

    $self->{$xml_key} = $build_hash;
    $self->parse_spec($build_spec_path);

    return $self;
}

#------------------------------------------------------------------------------

sub parse_spec($;$)
{
    my ($self, $build_spec_path) = @_;

    $self->check_required_vals($build_spec_path);

    $self->{$export_dirs_key} = {};

    my $aligns_ref = $self->{$xml_key}->{$alignments_tag};

    if (!defined($aligns_ref)) {
        die("Failed to find $alignments_tag in $build_spec_path"); 
    }

    my $align_set_ref = $aligns_ref->{$alignment_tag};

    if (!defined($align_set_ref)) {
        die("Failed to find $alignment_tag under $alignments_tag "
            . "in $build_spec_path");
    }

    foreach my $export_dir (keys %$align_set_ref) {
        # DEBUG
        # print "Export dir : ", $export_dir, "\n";

        my $export_dir_xml_ref = $align_set_ref->{$export_dir};
        $self->{$export_dirs_key}->{$export_dir}
          = new Casava::PostAlignment::QC::Export_Dir_Spec($export_dir_xml_ref);
    }
}

#------------------------------------------------------------------------------

sub check_required_vals($;$)
{
    my ($self, $build_spec_path) = @_;

    foreach my $required_val_tag (@required_val_tags) {
        my $val = $self->{$xml_key}->{$required_val_tag};

        if (!defined($val)) {
            die("Value must be specified for $required_val_tag "
                . "in $build_spec_path.");
        }
    }
}

#------------------------------------------------------------------------------

sub get_tile_selection_config_path($)
{
    my ($self) = @_;
    return $self->{$xml_key}->{$tile_selection_config_tag};
}

#------------------------------------------------------------------------------

sub get_min_yield($)
{
    my ($self) = @_;
    my $min_yield_gb = $self->{$xml_key}->{$min_yield_gb_tag};

    return (defined($min_yield_gb)
            ? $min_yield_gb * 1000 * 1000 * 1000
            : 0);
}

#------------------------------------------------------------------------------

sub get_run_script_fixed_args($)
{
    my ($self) = @_;
    return $self->{$xml_key}->{$run_script_fixed_args_tag};
}

#------------------------------------------------------------------------------

sub write($;$)
{
    my ($self, $output_path) = @_;

    open(OUT, ">$output_path")
      or die("Failed to open $output_path");

    my $xml_decl = '';
    my $xml_obj = new XML::Simple;

    print(OUT $xml_obj->XMLout($self->{$xml_key},
                               NoAttr => 1,
                               KeyAttr => {$alignment_tag => $export_dir_tag,
                                           $lane_tag => $number_tag},
                               XMLDecl => $xml_decl,
                               RootName => 'Build'));
}

#------------------------------------------------------------------------------

sub get_export_dirs($)
{
    my ($self) = @_;
    return keys %{$self->{$export_dirs_key}};
}

#------------------------------------------------------------------------------

sub get_export_dir_spec($;$)
{
    my ($self, $export_dir) = @_;
    return $self->{$export_dirs_key}->{$export_dir};
}

#------------------------------------------------------------------------------

sub dump($)
{
    my $self = shift;

    print "Build_Spec :\n", Data::Dumper->Dumper($self->{$xml_key});
}

#------------------------------------------------------------------------------

return 1;

#------------------------------------------------------------------------------

