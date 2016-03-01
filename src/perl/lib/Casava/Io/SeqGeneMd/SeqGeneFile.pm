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

package Casava::Io::SeqGeneMd::SeqGeneFile;

use strict;
use warnings "all";
use Carp;
use Exporter 'import';

use lib '@CASAVA_FULL_PERL_LIBDIR@'; # substituted by CMake during install

use Casava::Io::SeqGeneMd::GeneFeature;
use Casava::Io::SeqGeneMd::GeneCollection;
use Casava::Io::SeqGeneMd::Gene;
use Casava::Io::SeqGeneMd::RNA;
use Casava::Io::SeqGeneMd::UTR;
use Casava::Io::SeqGeneMd::CDS;

my $FALSE   = 0;
my $TRUE    = 1;


sub new {
    my $class = shift;
    my $args = shift;   # accept hash ref with args

    if (!defined $args) {
        die "Error - Seq_Gene_File has required arguments\n";

    # check for required args
    } elsif (! exists $args->{'fh'}) {
        die "Error - fh is required in Seq_Gene_File\n";
    } elsif (! exists $args->{'group_label'}) {
        die "Error - group_label is required in Seq_Gene_File\n";

    } else {
        bless $args, $class;
    }
}


sub fh {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'fh'} = $val;
    } else {
        return $self->{'fh'};
    }
}
sub group_label {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'group_label'} = $val;
    } else {
        return $self->{'group_label'};
    }
}
# a Gene_Collection object
sub gene_collection {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'gene_collection'} = $val;
    } else {
        return $self->{'gene_collection'};
    }
}






#
# print seq_gene to given FH, or STDOUT if FH not given;
# only print entries with this group label
# if rename_group_label is also given, rename these group labels
#
# 
sub print_seq_gene {
    my ($self, $FH_OUTPUT, $rename_group_label) = @_;
    
    $FH_OUTPUT      = *STDOUT if (!defined $FH_OUTPUT);
    my $FH_INPUT    = $self->fh;
    my $group_label = $self->group_label;

    # print seq_gene.md header
    my $hdr = <$FH_INPUT>;
    print $FH_OUTPUT $hdr;

    my $match_group_label = $FALSE;

    while (<$FH_INPUT>) {
        chomp;
        my ($tax_id,     $chromosome,   $chr_start,  $chr_end,
            $chr_orient, $contig,       $ctg_start,  $ctg_end,
            $ctg_orient, $feature_name, $feature_id, $feature_type,
            $label,      $accver,       $evidence
        ) = split(/\t/);


        # restrict parsing to primary assembly
        next if ($label ne $group_label);
        $match_group_label = $TRUE;

        # rename group label, if desired
        if (defined $rename_group_label) {
            $label = $rename_group_label;
        }

        # print seq_gene line
        my $print_line = join ( "\t",
            (   $tax_id,     $chromosome,   $chr_start,  $chr_end,
                $chr_orient, $contig,       $ctg_start,  $ctg_end,
                $ctg_orient, $feature_name, $feature_id, $feature_type,
                $label,      $accver,       $evidence,
            ) );
        print $FH_OUTPUT "$print_line\n";
    }

    # throw error if group_label wasn't matched
    if ($match_group_label == $FALSE) {
        die "ERROR: group label " . $self->group_label . " didn't match any entry in seq_gene.md\n";
    }
}








# parse seq_gene file into objects, and print to fh; type is string that must
# contain refFlat or refGene
sub print_refFlat_refGene {
    my ($self, $type, $fh) = @_;
    $fh = *STDOUT if (!defined $fh);

    # check print type
    if ( (!defined $type) 
            || (($type !~ /refFlat/) && ($type !~ /refGene/))){
        print STDERR "ERROR: print type must be refFlat or refGene\n";
        return;
    }

    $self->gene_collection( Casava::Io::SeqGeneMd::GeneCollection->new() );
    
    # first sort seq_gene by chr, then chr_start, then feature type
    # so we always see GENE first, then RNAs, then UTRs and CDSs
    my @lines = $self->_sort_seq_gene();

    # next parse entire file for GENE entries, counting the number of times
    # each GeneId is found. Since a GeneId may occur more than once in a genome,
    # and we must produce a unique transcript_id for refFlat, we will append
    # an 'index' number to transcipt_ids of genes we observe multiple times.
    my $gene_cnt = $self->_count_multi_genes(\@lines);

    my $match_group_label = $FALSE;

    foreach my $line (@lines) {
        chomp $line;
        my (
            $tax_id,     $chr,          $chr_start,  $chr_end,
            $chr_orient, $contig,       $ctg_start,  $ctg_end,
            $ctg_orient, $feature_name, $feature_id, $feature_type,
            $label,      $transcript,   $evidence
        ) = split(/\t/, $line);

        # restrict parsing to primary assembly
        next if ($label ne $self->group_label);
        $match_group_label = $TRUE;

        # ignore features that don't have a transcript id
        next if $feature_type ne "GENE" and $transcript eq "-";

        if ( $feature_id =~ /^GeneID:(\d+)$/ ) {
            my $gene_id = $1;


            if ( $feature_type eq "GENE" ) {

                # empty and print past GENEs that we've passed up
                $self->gene_collection->print_finished_genes($type, $fh, $chr, $chr_start);

                # get index of multi-gene, if this gene was found more than once
                my $cnt = 0;
                if (exists $gene_cnt->{$gene_id}) {
                    $cnt = $gene_cnt->{$gene_id}--;
                }

                # a GENE 'feature_name' is actually the gene symbol
                my $gene = Casava::Io::SeqGeneMd::Gene->new( {
                    'gene_id'   => $gene_id,
                    'symbol'    => $feature_name,
                    'chr'       => $chr,
                    'start'     => $chr_start,
                    'end'       => $chr_end,
                    'orient'    => $chr_orient,
                    'index'     => $cnt,
                     } );
                $self->gene_collection->add_gene($gene);

            } elsif ( $feature_type eq "RNA" ) {

                # a RNA 'feature_name' is actually the transcript_id 
                # (which is redundant with $transcript)
                my $rna = Casava::Io::SeqGeneMd::RNA->new( {
                    'transcript_id' => $transcript,
                    'chr'           => $chr,
                    'start'         => $chr_start,
                    'end'           => $chr_end,
                    'orient'        => $chr_orient,
                     } );
                $self->gene_collection->add_feature($gene_id, $transcript, $rna);

            } elsif ( $feature_type eq "UTR" ) {

                # a RNA 'feature_name' is actually the transcript_id 
                # (which is redundant with $transcript)
                my $utr = Casava::Io::SeqGeneMd::UTR->new( {
                    'chr'           => $chr,
                    'start'         => $chr_start,
                    'end'           => $chr_end,
                    'orient'        => $chr_orient,
                     } );
                $self->gene_collection->add_feature($gene_id, $transcript, $utr);

            } elsif ( $feature_type eq "CDS" ) {

                # a CDS 'feature_name' is actually the protein_id
                my $cds = Casava::Io::SeqGeneMd::CDS->new( {
                    'chr'           => $chr,
                    'start'         => $chr_start,
                    'end'           => $chr_end,
                    'orient'        => $chr_orient,
                    'protein_id'    => $feature_name,
                     } );
                $self->gene_collection->add_feature($gene_id, $transcript, $cds);

            # unrecognized feature type
            } elsif ( $feature_type eq "PSEUDO" ) {

            } else {
                print STDERR "WARNING: Unrecognized feature_type $feature_type. Skipping.\n";
            }
    
    
        } else {
            print STDERR "WARNING: feature_id $feature_id didn't include GeneID\n\n";
        }
        
    }
    # print remaining genes
    $self->gene_collection->print_all_genes($type, $fh);


    # throw error if group_label wasn't matched
    if ($match_group_label == $FALSE) {
        die "ERROR: group label " . $self->group_label . " didn't match any entry in seq_gene.md\n";
    }
}




# return hash of number of times each GENE's GeneID was encountered
# but only for genes observed more than once
sub _count_multi_genes {
    my ($self, $lines) = @_;
    my %gene_cnt;

    foreach my $line (@$lines) {
        chomp $line;
        my (
            $tax_id,     $chr,          $chr_start,  $chr_end,
            $chr_orient, $contig,       $ctg_start,  $ctg_end,
            $ctg_orient, $feature_name, $feature_id, $feature_type,
            $label,      $transcript,   $evidence
        ) = split(/\t/, $line);

        # restrict parsing to primary assembly
        next if ($label ne $self->group_label);

        # if GENE, add to count of this GeneID
        if ( $feature_id =~ /^GeneID:(\d+)$/ ) {
            my $gene_id = $1;

            if ( $feature_type eq "GENE" ) {
                $gene_cnt{$gene_id}++;
            }
        }
    }

    # now remove genes that were only seen once
    foreach my $gene (keys %gene_cnt) {
        if ($gene_cnt{$gene} == 1) {
            delete $gene_cnt{$gene};   
        }
    }
    return \%gene_cnt;
}





# sort seq_gene lines -- see _sort_Seq_Gene method for criteria
sub _sort_seq_gene {
    my $self = shift;
    my $fh = $self->fh;
    my @lines;

    #discard header
    my $hdr = <$fh>;
    
    # slurp all lines into a list
    while(<$fh>) {
       push @lines, $_; 
    }
    return (sort _sort_Seq_Gene @lines);
}





# Sort Seq_Gene lines in order of chromosome, chrom start pos, then feature_type
# This sort allows parsing in a simple single pass, so we encounter GENE types
# first, then RNA types, and then their associated UTR and CDS types.
sub _sort_Seq_Gene {

    my @first  = split( '\t', $a);
    my @second = split( '\t', $b);
    my $compare;

    # first sort (alpha) by chromosome
    $compare = ($first[1] cmp $second[1] );
    if ($compare != 0) {
        return $compare;
    }
    # sort (numeric) by chr start position
    $compare = ($first[2] <=> $second[2] );
    if ($compare != 0) {
        return $compare;
    }
    # sort by feature_type in this order: GENE, RNA, CDS, UTR, PSEUDO
    my $first_type = $first[11];
    my $second_type = $second[11];
    my $first_rank;
    my $second_rank;

    my %feature_type = (
        "GENE"  => 1,
        "RNA"   => 2,
        "UTR"   => 3,
        "CDS"   => 4,
        "PSEUDO"=> 5,   );

    if (!exists $feature_type{$first_type}) {
        $first_rank = 10;
    } else {
        $first_rank = $feature_type{$first_type};
    }
    if (!exists $feature_type{$second_type}) {
        $second_rank = 10;
    } else {
        $second_rank = $feature_type{$second_type};
    }
    $compare = ($first_rank <=> $second_rank);

    return $compare;
}

1;
