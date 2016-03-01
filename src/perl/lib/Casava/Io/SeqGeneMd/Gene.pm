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

package Casava::Io::SeqGeneMd::Gene;

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
        die "Error - Gene has required arguments\n";

    # check for required args
    } elsif (! exists $args->{'gene_id'}) {
        die "Error - gene_id is required in Gene\n";
    } elsif (! exists $args->{'symbol'}) {
        die "Error - symbol is required in Gene\n";
    } elsif (! exists $args->{'index'}) {
        die "Error - index is required in Gene\n";

    } else {
        # set defaults
        if (! exists $args->{'RNAs'}) {
            $args->{'RNAs'} = [];
        }
        bless $args, $class;
    }
}


# seq_gene.md feature_type GENE
sub gene_id {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'gene_id'} = $val;
    } else {
        return $self->{'gene_id'};
    }
}
sub symbol {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'symbol'} = $val;
    } else {
        return $self->{'symbol'};
    }
}
# index is the instance of this gene in the genome (0 for single genes (most))
sub index {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'index'} = $val;
    } else {
        return $self->{'index'};
    }
}
# array ref of RNA objects
sub RNAs {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'RNAs'} = $val;
    } else {
        return $self->{'RNAs'};
    }
}




# print this Gene to fh, and append transcript_id_suffix to all transcript_ids 
# (since the same gene can occur in mulitple places in the genome
sub print_refFlat_refGene {
    my ($self, $type, $fh, $refGene_cnt) = @_;

    # format index for making unique the transcript-id of multi-genes
    my $index = $self->index;
    if ($index == 0) {
        $index = "";
    } else {
        $index = "-" . $index;
    }

    foreach my $rna (@{$self->RNAs}) {

        if ($type =~ /refFlat/) {
            $rna->print_refFlat($fh, $self->symbol, $index);

        } elsif ($type =~ /refGene/) {
            $rna->print_refGene($fh, $self->symbol, $index, $refGene_cnt);
            $refGene_cnt++;

        } else {
            print STDERR "Warning - unrecognized print type $type\n";
        }
    }
    return $refGene_cnt;
}






# add RNA, UTR, or CDS 
sub add_feature {
    my ($self, $transcript_id, $feature) = @_;

    # check for RNA
    if ($self->package_basename(ref $feature) eq 'RNA') {
        $self->add_rna($feature);

    # find RNA with transcript_id of feature
    } else {
        my $found_rna = $FALSE;
        foreach my $rna (@{$self->RNAs}) {
            
            if ($rna->transcript_id eq $transcript_id) {
                $rna->add_feature($feature);
                $found_rna = $TRUE;
                last;
            }
        }
        if ($found_rna == $FALSE) {
            print STDERR "Warning - couldn't find RNA that matches transcript_id $transcript_id of current feature.\n";
        }
    }
}


# add an RNA object to RNAs
sub add_rna { 
    my ($self, $rna) = @_;
    my $RNAs = $self->RNAs;

    push @{$RNAs}, $rna;
}


1;
