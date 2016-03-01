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

package Casava::Io::SeqGeneMd::GeneCollection;

use strict;
use warnings "all";
use Carp;
use Exporter 'import';

use lib '@CASAVA_FULL_PERL_LIBDIR@'; # substituted by CMake during install

my $FALSE   = 0;
my $TRUE    = 1;


sub new {
    my $class = shift;
    my $args = shift;   # accept hash ref with args

    # set defaults for genes and refGene_cnt
    if (defined $args) {
        if (!exists $args->{'genes'}) {
            $args->{'genes'} = [];
        }
        if (!exists $args->{'refGene_cnt'}) {
            $args->{'refGene_cnt'} = 1;
        }
    } else {
        $args->{'genes'} = [];
        $args->{'refGene_cnt'} = 1;
    }
    bless $args, $class;
}


# array ref of Gene objects
sub genes {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'genes'} = $val;
    } else {
        return $self->{'genes'};
    }
}
sub refGene_cnt {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'refGene_cnt'} = $val;
    } else {
        return $self->{'refGene_cnt'};
    }
}






# print all genes that are 'finished'. That is, as we parse a seq_gene file,
# print genes that we've been adding features to where we're now on a diff 
# chromosome or beyond the end coord of that gene.
sub print_finished_genes {
    my ($self, $type, $fh, $cur_chr, $cur_start) = @_;

    my @cur_genes;
    foreach my $gene ( @{$self->genes} ) {
    
        # keep GENEs that we're still gathering features for
        if (($gene->chr eq $cur_chr) 
                && ($gene->end >= $cur_start)) {
            push @cur_genes, $gene;

        } else {
            my $new_cnt = $gene->print_refFlat_refGene($type, $fh, $self->refGene_cnt);
            $self->refGene_cnt($new_cnt);
        }
    }
    # update genes with GENEs in progress
    $self->genes(\@cur_genes);
}


# print all genes, intended for printing final genes at end of file parsing
sub print_all_genes {
    my ($self, $type, $fh) = @_;

    foreach my $gene ( @{$self->genes} ) {
        my $new_cnt = $gene->print_refFlat_refGene($type, $fh, $self->refGene_cnt);
        $self->refGene_cnt($new_cnt);
    }
    # clear genes
    $self->genes( [] );
}




# add a RNA, UTR, or CDS object
sub add_feature {
    my ($self, $new_gene_id, $new_transcript_id, $new_feature) = @_;
    my $found_gene = $FALSE;

    # check if we've already seen this gene_id
    foreach my $gene ( @{$self->genes} ) {

        if ($gene->gene_id eq $new_gene_id) {

            $gene->add_feature($new_transcript_id, $new_feature);
            $found_gene = $TRUE;
            last;
        }
    }
    # should've already seen this gene_id as a GENE -- warn if not
    if ($found_gene == $FALSE) {
        print STDERR "Warning- encountered GeneId $new_gene_id before GENE. Missing GENE for this feature? Ignoring this feature.\n";
    }
}




# add a Gene object
sub add_gene {
    my ($self, $gene) = @_;
    my $genes = $self->genes;
    
    push @{$genes}, $gene;
}

1;
