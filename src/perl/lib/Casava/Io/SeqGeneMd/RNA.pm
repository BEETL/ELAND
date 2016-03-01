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

package Casava::Io::SeqGeneMd::RNA;

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
        die "Error - RNA has required arguments\n";

    # check for required args
    } elsif (! exists $args->{'transcript_id'}) {
        die "Error - transcript_id is required in RNA\n";

    } else {
        # set defaults
        if (! exists $args->{'UTRs'}) {
            $args->{'UTRs'} = [];
        }
        if (! exists $args->{'CDSs'}) {
            $args->{'CDSs'} = [];
        }
        bless $args, $class;
    }
}


# from transcript field of seq_gene.md, feature_type RNA
sub transcript_id {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'transcript_id'} = $val;
    } else {
        return $self->{'transcript_id'};
    }
}
# array ref of UTR objects
sub UTRs {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'UTRs'} = $val;
    } else {
        return $self->{'UTRs'};
    }
}
# array ref of CDS objects
sub CDSs {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'CDSs'} = $val;
    } else {
        return $self->{'CDSs'};
    }
}




# print refFlat to fh. RefFlat start coords are zero-based and end coords are 
# 1-based.
sub print_refFlat {
    my ($self, $fh, $gene_symbol, $transcript_id_suffix) = @_;

    my ($num_exons, $exon_starts, $exon_ends) = $self->_exon_stats();

    my $refFlat = join("\t",    $gene_symbol, 
                                $self->transcript_id . $transcript_id_suffix,
                                $self->chr,
                                $self->orient,
                                $self->start - 1,
                                $self->end,
                                $self->_cds_start() - 1,
                                $self->_cds_end(),
                                $num_exons,
                                $exon_starts,
                                $exon_ends,        );

    print $fh "$refFlat\n";
}



# print refGene to fh. RefGene start coords are zero-based and end coords are 
# 1-based.
sub print_refGene {
    my ($self, $fh, $gene_symbol, $transcript_id_suffix, $refGene_cnt) = @_;

    my ($num_exons, $exon_starts, $exon_ends) = $self->_exon_stats();

    my $refGene = join("\t",    $refGene_cnt,
                                $self->transcript_id . $transcript_id_suffix,
                                $self->chr,
                                $self->orient,
                                $self->start - 1,
                                $self->end,
                                $self->_cds_start() - 1,
                                $self->_cds_end(),
                                $num_exons,
                                $exon_starts,
                                $exon_ends,        );
    $refGene .= "\t\t\t\t\t";

    print $fh "$refGene\n";
}



# return number of exons, exon start coords, and exon end coords for refFlat
# or refGene printing. Coords are comma-delimited, start is zero-based, end 
# is one-based
sub _exon_stats {
    my ($self) = @_;

    # sort UTRs and CDSs by start coord
    my @sorted_features = sort { $a->start <=> $b->start } 
            ( @{$self->UTRs}, @{$self->CDSs} );

    # find exon starts & ends
    my @starts;
    my @ends;
    my $last_end = -1;
    foreach my $f (@sorted_features) {

        push @starts, $f->start - 1 if ($f->start != $last_end + 1);
        push @ends, $last_end     if (($last_end != -1) 
                                    && ($last_end != $f->start - 1));
        $last_end = $f->end;
    }
    push @ends, $last_end if ($last_end != -1);

    # convert lists to comma-delimited strings
    my $num_exons = @starts;
    my $starts = join (',', @starts);
    my $ends = join (',', @ends);

    # add trailing comma (UCSC example has this)
    $starts .= ',' if ($num_exons > 0);
    $ends .= ',' if ($num_exons > 0);

    return ($num_exons, $starts, $ends);
}





# return cds start coord 
sub _cds_start {
    my ($self) = @_;

    # if no coding regions in this RNA, return transcription end
    if (@{$self->CDSs} == 0) {
        return $self->end;

    } else {
        my $cds_start = $self->CDSs->[0]->start;

        foreach my $cds (@{$self->CDSs}) {
            if ($cds->start < $cds_start) {
                $cds_start = $cds->start;
            }
        }
        return $cds_start; 
    }
}





# return cds end coord 
sub _cds_end {
    my ($self) = @_;

    # if no coding regions in this RNA, return transcription end
    if (@{$self->CDSs} == 0) {
        return $self->end;

    } else {
        my $cds_end = $self->CDSs->[0]->end;

        foreach my $cds (@{$self->CDSs}) {
            if ($cds->end > $cds_end) {
                $cds_end = $cds->end;
            }
        }
        return $cds_end; 
    }
}





# add feature: UTR or CDS
sub add_feature {
    my ($self, $feature) = @_;

    # match UTR
    if ($self->package_basename(ref $feature) eq 'UTR') {
        my $UTRs = $self->UTRs;
        push @{$UTRs}, $feature;
        
    # match CDS
    } elsif ($self->package_basename(ref $feature) eq 'CDS') {
        my $CDSs = $self->CDSs;
        push @{$CDSs}, $feature;
    
    } else {
        print "Warning - unrecognized feature type: " . $self->package_basename(ref $feature) . "\n";
    }
}


1;
