=head1 LICENSE

Copyright (c) 2007-2009 Illumina, Inc.

This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

This file is part of the Consensus Assessment of Sequence And VAriation
(CASAVA) software package.

=cut

package Casava::Io::SeqGeneMd::GeneFeature;

use strict;
use warnings "all";
use Carp;
use Exporter 'import';

use lib '@CASAVA_FULL_PERL_LIBDIR@'; # substituted by CMake during install

my $FALSE   = 0;
my $TRUE    = 1;

# intended to be base class of RNA, UTR, and CDS feature_types. Possibly PSEUDO too, though not currently used.

sub new {
    my $class = shift;
    my $args = shift;   # accept hash ref with args

    if (!defined $args) {
        die "Error - Gene_Feature has required arguments\n";

    # check for required args
    } elsif (! exists $args->{'chr'}) {
        die "Error - chr is required in Gene_Feature\n";
    } elsif (! exists $args->{'start'}) {
        die "Error - start is required in Gene_Feature\n";
    } elsif (! exists $args->{'end'}) {
        die "Error - end is required in Gene_Feature\n";
    } elsif (! exists $args->{'orient'}) {
        die "Error - orient is required in Gene_Feature\n";

    } else {
        bless $args, $class;
    }
}


sub chr {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'chr'} = $val;
    } else {
        return $self->{'chr'};
    }
}
sub start {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'start'} = $val;
    } else {
        return $self->{'start'};
    }
}
sub end {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'end'} = $val;
    } else {
        return $self->{'end'};
    }
}
sub orient {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'orient'} = $val;
    } else {
        return $self->{'orient'};
    }
}



# given a full package name, return it's basename (right-most name)
sub package_basename {
    my ($self, $pkg) = @_;
    my @elems = split /::/, $pkg;
    return pop @elems
}


1;
