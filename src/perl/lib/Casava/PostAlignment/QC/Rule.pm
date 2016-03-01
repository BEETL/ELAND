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

Rule.pm

=cut

=head1 DESCRIPTION

Filter rule class

=cut

=head1 AUTHOR

Richard Shaw

=cut


use warnings;
use strict;

#------------------------------------------------------------------------------

use lib '@CASAVA_FULL_PERL_LIBDIR@';

#------------------------------------------------------------------------------

package Casava::PostAlignment::QC::Rule;

#------------------------------------------------------------------------------

my $lbound_tag = 'lbound';
my $ubound_tag = 'ubound';
my $action_tag = 'action';
my $content_tag = 'content';

my $orig_xpath_key = 'orig_xpath';
my $curr_xpath_key = 'curr_xpath';
my $xpath_stack_key = 'xpath_stack';

my $bound_type_key = 'bound_type';
my $action_key = 'action';
my $value_key = 'value';

#------------------------------------------------------------------------------

sub new($;$;$)
{
    my ($class, $xpath, $rule_ref) = @_;
    my $self = {};
    bless $self, $class;

    my @tags = keys %{$rule_ref};
    my $bound_type = $tags[0];

    my $rule_info_ref = $rule_ref->{$bound_type};
    my $action = $rule_info_ref->{$action_tag};
    my $value = $rule_info_ref->{$content_tag};

    $self->{$orig_xpath_key} = $xpath;
    $self->{$curr_xpath_key} = $xpath;
    $self->{$xpath_stack_key} = ();

    $self->{$bound_type_key} = $bound_type;
    $self->{$action_key} = $action;
    $self->{$value_key} = $value;

    # DEBUG
    # print join(' : ', $xpath, $bound_type, $action, $value), "\n";

    return $self;
}

#------------------------------------------------------------------------------

sub get_xpath($)
{
    my $self = shift;
    return $self->{$curr_xpath_key};
}

#------------------------------------------------------------------------------
# If the current XPath is /TileErrorsByLane/Lane*/Read*/Tile*/errorPF
# and the $key is `Read', strip out the `/Read*' from the current xpath
# and tell the caller that Read is wild : $any_val is 1, $val is -1.
#
# If the current XPath is /TileErrorsByLane/Lane*/Read1/Tile*/errorPF
# and the $key is `Read', strip out the `/Read1' from the current xpath
# and tell the caller the Read number that matches : $any_val is 0, 
#
# If $key is matched in either of the above ways, return 1; else return 0.

sub reduce_xpath($;$;$;$)
{
    my ($self, $key, $any_val_ref, $val_ref) = @_;

    my $curr_xpath = $self->{$curr_xpath_key};

    # Push the current xpath onto the stack - whether or not it changes.
    push (@{$self->{$xpath_stack_key}}, $curr_xpath);

    my $did_reduce = 0;
    $$any_val_ref = 0;
    $$val_ref = -1;

    if ($curr_xpath =~ /(.*)(\/${key}\*)(\/.*)/) {
        # Wildcard
        $curr_xpath = $1 . $3;
        $did_reduce = 1;

        $$any_val_ref = 1;
    } elsif ($curr_xpath =~ /(.*)(\/${key})(\d+)(\/.*)/) {
        # Specific numerical value
        $curr_xpath = $1 . $4;
        $did_reduce = 1;

        $$val_ref = $3;
    }

    if ($did_reduce) {
        $self->{$curr_xpath_key} = $curr_xpath;
    }

    return $did_reduce;
}

#------------------------------------------------------------------------------

sub undo_reduce_xpath($)
{
    my $self = shift;
    $self->{$curr_xpath_key} = pop(@{$self->{$xpath_stack_key}});
}

#------------------------------------------------------------------------------

sub localised_xpath($;$;$)
{
    my ($self, $key, $val) = @_;

    my $curr_xpath = $self->{$curr_xpath_key};
    my $local_xpath = '';

    if ($curr_xpath =~ /(.*)(\/${key}\*)(\/.*)/) {
        $local_xpath = $curr_xpath;
        $local_xpath =~ s|\/${key}\*\/|\/${key}${val}\/|;
    }

    return $local_xpath;
}

#------------------------------------------------------------------------------

sub val_ok($;$)
{
    my ($self, $val) = @_;

    if (!defined($val)) {
        return 0;
    }

    my $bound_type = $self->{$bound_type_key};
    my $bound_val = $self->{$value_key};

    if ($bound_type eq $lbound_tag) {
        return ($val >= $bound_val);
    } elsif ($bound_type eq $ubound_tag) {
        return ($val <= $bound_val);
    }

    return 0;
}

#------------------------------------------------------------------------------

sub to_str($)
{
    my $self = shift;
    return join(' : ',
                $self->{$orig_xpath_key},
                $self->{$action_key},
                $self->{$bound_type_key},
                $self->{$value_key});
}

#------------------------------------------------------------------------------

sub dump($)
{
    my $self = shift;
    print $self->to_str(), "\n";
}

#------------------------------------------------------------------------------

