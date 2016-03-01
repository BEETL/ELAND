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

package Casava::Io::SeqGeneMd::UTR;

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

    if (defined $args) {
        bless $args, $class;
    } else {
        bless {}, $class;
    }
}

# not recording feature_name, which should be the same as transcript_id


1;
