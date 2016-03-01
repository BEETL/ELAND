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

Rule_Set.pm

=cut

=head1 DESCRIPTION

Container class for filter rules

=cut

=head1 AUTHOR

Richard Shaw

=cut


use warnings;
use strict;

use XML::Simple;

#------------------------------------------------------------------------------

use lib '@CASAVA_FULL_PERL_LIBDIR@';

use Casava::PostAlignment::QC::Rule;

#------------------------------------------------------------------------------

package Casava::PostAlignment::QC::Rule_Set;

#------------------------------------------------------------------------------

sub is_magic_tag($)
{
    my $tag = shift;

    if (($tag eq 'Lane')
        || ($tag eq 'Read')
        || ($tag eq 'Tile')) {
        return 1;
    }

    return 0;
}

#------------------------------------------------------------------------------

sub new {
    my ($class, $rule_set_path) = @_;
    my $self = {};
    bless $self, $class;

    if (! -f $rule_set_path) {
        print "File $rule_set_path not found\n";
        return 0;
    }

    my $rule_set_xml = new XML::Simple(suppressempty => '',
                                       XMLDecl => 1);

    my $rule_set_hash_ref
      = $rule_set_xml->XMLin($rule_set_path,
                             forcearray
                             => [qw(table Lane Read Tile column statistic)],
                             keyattr => [qw(name number)]);

    my $odd_level = 0;
    parse_level($self, '', $rule_set_hash_ref, $odd_level, '');

    return $self;
}

#------------------------------------------------------------------------------

sub parse_level($;$;$;$;$);

sub parse_level($;$;$;$;$)
{
    my ($self, $xpath, $hash_ref, $odd_level, $last_tag) = @_;

    my @tags = keys %{$hash_ref};

    $odd_level = ($odd_level == 0) ? 1 : 0;

    foreach my $tag (@tags) {
        if (($tag eq 'lbound') || ($tag eq 'ubound')) {
            # print $xpath, "\n";
            $self->{$xpath} = Casava::PostAlignment::QC::Rule->new($xpath, $hash_ref);
            return;
        }

        my $next_xpath = $xpath;

        if ($odd_level == 0) {
            $next_xpath = join('/', $xpath,
                               (is_magic_tag($last_tag)
                                ? $last_tag . $tag
                                : $tag));
        }

        my $next_hash_ref = $hash_ref->{$tag};

        parse_level($self, $next_xpath, $next_hash_ref, $odd_level, $tag);
    }
}

#------------------------------------------------------------------------------
# To iterate through the Rules, get the array of xpaths by calling this fn,
# then iterate through them passing each to the rule fn.

sub xpaths($)
{
    my ($self) = @_;
    return keys %{$self};
}

#------------------------------------------------------------------------------

sub rule($;$)
{
    my ($self, $xpath) = @_;
    return $self->{$xpath};
}

#------------------------------------------------------------------------------

sub dump($)
{
    my $self = shift;

    foreach my $xpath (keys %$self) {
        $self->{$xpath}->dump();
    }
}

#------------------------------------------------------------------------------

return 1;

#------------------------------------------------------------------------------
