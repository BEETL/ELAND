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

package Casava::Io::SeqGeneMd::CDS;

use strict;
use warnings "all";
use Carp;
use Exporter 'import';

use lib '@CASAVA_FULL_PERL_LIBDIR@'; # substituted by CMake during install

my $FALSE   = 0;
my $TRUE    = 1;

our @ISA = 'Casava::Io::SeqGeneMd::GeneFeature';

sub new {
    my $class = shift;
    my $args = shift;   # accept hash ref with args

    if (!defined $args) {
        die "Error - CDS has required arguments\n";

    # check for required args
    } elsif (! exists $args->{'protein_id'}) {
        die "Error - protein_id is required in CDS\n";

    } else {
        bless $args, $class;
    }
}

# from feature_name seq_gene.md field
sub protein_id {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'protein_id'} = $val;
    } else {
        return $self->{'protein_id'};
    }
}


1;
