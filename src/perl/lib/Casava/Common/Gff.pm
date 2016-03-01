# PROJECT: CASAVA
# MODULE:  $RCSfile$
# AUTHOR:  Lukasz Szajkowski
#
# Copyright (c) 2008 Illumina
# This source file is covered by the "Illumina Public Source License"
# agreement and bound by the terms therein.
#

=pod

=head1 NAME

Casava::Common::Gff.pm - Perl utility library for writing to GFF.

=head1 SYNOPSIS

# include what functions you need... 
use Casava::Common::Gff qw();  

=head1 DESCRIPTION

The library is a simple implementation of GFF2 format. The main advantage of 
the library is that it can be distributed with the Illumina software and 
does not require user to install eny addtional components.

Exports:    
    createGFF(\%);
    addFeature(\%; \@);
    addSNPFeature(\%; $;$;$; $; $; $);
        
# Global variable

=head1 AUTHORSHIP

Copyright (c) 2008 Illumina

This software is covered by the "Illumina Genome Analyzer Software
License Agreement" and the "Illumina Source Code License Agreement",
and certain third party copyright/licenses, and any user of this
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party
copyright/license notices).

Created by Lukasz Szajkowski <lszajkowski@illumina.com>

=head1 SOURCE CODE

The most current release of the Perl source code 
for this module is available through CVS at genie01 a.k.a. 10.44.0.81 
cvs co BullFrog

=cut

package Casava::Common::Gff;

#
# Place functions/variables you want to *export*, ie be visible
# from the caller package into @EXPORT_OK
#
BEGIN {
    use Exporter();
    @ISA    = qw(Exporter);
    @EXPORT = qw($VERSION &createGFF &addSNPFeature &dumpGFFHeader &dumpGFF
      &clearFeatureSet);
    @EXPORT_OK = qw();
}
use warnings FATAL => 'all';
use strict;
use POSIX qw(strftime);

use Carp;
use Casava::Common::Log;

sub createGFF(\%);
sub addGFFFeature(\%; \@);
sub addSNPFeature(\%; $;$;$; $; $; $);
sub dumpGFFHeader(\%; $);
sub dumpGFF(\%; $);
sub clearFeatureSet(\%);

our %featureFields = (
    seqname    => 0,
    source     => 1,
    feature    => 2,
    start      => 3,
    end        => 4,
    score      => 5,
    strand     => 6,
    frame      => 7,
    attributes => 8,
    comments   => 9
);

=pod

=head1 The procedure creates GFF Hash Map.

=over 4

=item createGFF($gffConfRef)

The procedure procedure creates GFF Hash Map.

Parameters:
    $gffConfRef - GFF configurationz

Returns:
    GFF Hash Map Ref
=back

=cut

sub createGFF(\%) {
    croak "ERROR: createGFF\n" unless ( @_ == 1 );
    my ($gffConfRef) = @_;
    my $gffSetRef    = ();
    my @featureSet   = ();
    if ( !defined $gffConfRef->{source} ) {
        $gffConfRef->{source} = 'Illumina_CASAVA';
    }
    if ( !defined $gffConfRef->{title} ) {
        $gffConfRef->{title} = 'Illumina CASAVA SNP file';
    }
    $gffSetRef->{config}     = $gffConfRef;
    $gffSetRef->{featureSet} = \@featureSet;
    if ( !defined $gffSetRef->{config}->{track} ) {
        $gffSetRef->{config}->{track}->{name}        = 'Illumina_SNP';
        $gffSetRef->{config}->{track}->{description} =
          'Illumina SNPs generated by CASAVA';
        $gffSetRef->{config}->{track}->{color} = 'yellow';
    }    # if
    $gffSetRef->{config}->{start} = 99999999999999999;
    $gffSetRef->{config}->{end}   = -1;
    $gffSetRef->{config}->{time}  = strftime '%H:%M:%S %d-%m-%y', localtime;
    $gffSetRef->{config}->{date}  = strftime '%d-%m-%y', localtime;
    return $gffSetRef;
}

=pod

=head1 The procedure prints the GFF header to file.

=over 4

=item dumpGFFHeader($gffConfRef, $file)

The procedure prints the GFF header to file.

Parameters:
    $gffSetRef  - GFF feature set
    $fileOut    - gff file
Returns:
    GFF Hash Map Ref
=back

=cut

sub dumpGFFHeader(\%; $) {
    croak "ERROR: dumpGFFHeader\n" unless ( @_ == 2 );
    my ( $gffSetRef, $fileOut ) = @_;
    print $fileOut "##gff-version 3\n";
    print $fileOut "#$gffSetRef->{config}->{title}\n";
    print $fileOut "#Created $gffSetRef->{config}->{time}\n";
    if (defined $gffSetRef->{config}->{seqname}) {
    print $fileOut "##sequence-region $gffSetRef->{config}->{seqname} "
      . "$gffSetRef->{config}->{start} "
      . "$gffSetRef->{config}->{end}" . "\n";
    print $fileOut
      "browser position $gffSetRef->{config}->{seqname}"
      . ":$gffSetRef->{config}->{start}-$gffSetRef->{config}->{end}\n";
    }
    print $fileOut
      "track name=$gffSetRef->{config}->{track}->{name} description=\""
      . $gffSetRef->{config}->{track}->{description}
      . "\" color=0,0,255\n"    
      ;
}

=pod

=head1 The procedure prints the GFF feature set to file.

=over 4

=item dumpGFFHeader($gffConfRef, $file)

The procedure prints the GFF feature set to file.

Parameters:
    $gffSetRef  - GFF feature set
    $fileOut    - gff file
Returns:
    GFF Hash Map Ref
=back

=cut

sub dumpGFF(\%; $) {
    croak "ERROR: dumpGFFHeader\n" unless ( @_ == 2 );
    my ( $gffSetRef, $fileOut ) = @_;
    foreach my $featureRef ( @{ $gffSetRef->{featureSet} } ) {
        my $featureStr = join( "\t", @{$featureRef} );
        print $fileOut $featureStr, "\n";
    }
}

=pod

=head1 The procedure creates GFF Hash Map.

=over 4

=item addFeature($gffConfRef)

The procedure procedure creates GFF Hash Map.

Parameters:
     $gffSetRef  - GFF feature set

Returns:
    Nothing
=back

=cut

sub addGFFFeature(\%; \@) {
    croak "ERROR: addFeature\n" unless ( @_ == 2 );
    my ( $gffSetRef, $featureRef ) = @_;
    my $start = $featureRef->[ $featureFields{start} ];
    my $end   = $featureRef->[ $featureFields{end} ];
    if ( !defined $gffSetRef->{config}->{seqname} ) {
        $gffSetRef->{config}->{seqname} =
          $featureRef->[ $featureFields{seqname} ];
    }    # if
    else {
        if ( $gffSetRef->{config}->{seqname} ne
            $featureRef->[ $featureFields{seqname} ] )
        {
            errorExit "ERROR: addGFFFeature seqname(chrom) " . "should be unique\n";
        }    # if
    }    # else
    if ( $gffSetRef->{config}->{start} > $start ) {
        $gffSetRef->{config}->{start} = $start;
    }    # if
    if ( $gffSetRef->{config}->{end} < $end ) {
        $gffSetRef->{config}->{end} = $end;
    }    # if
    push @{ $gffSetRef->{featureSet} }, $featureRef;
}

=pod

=head1 The procedure add SNP feature to FeatureSet.

=over 4

=item addFeature($gffConfRef, $chrom, $start, $end, 
    $score, $strand, $attribute)

The procedure add SNP feature to FeatureSet.

Parameters:
    $gffSetRef  - GFF feature set
    $chrom      - chromosome name
    $start      - start position
    $end        - end position
    $score      - snp score
    $strand     - strand
    $attribute  - unused 

Returns:
    Nothing
=back

=cut

sub addSNPFeature(\%;$;$;$;$;$;$) {
    croak "ERROR: addSNPFeature\n" unless ( @_ == 7 );
    my ( $gffSetRef, $chrom, $start, $end, $score, $strand, $attribute ) = @_;
    my @feature = (
        $chrom, $gffSetRef->{config}->{source},
        'SNP', $start, $end, $score, $strand, '.', $attribute
    );
    addGFFFeature( %{$gffSetRef}, @feature );
}

=pod

=head1 The procedure removes all features from FeatureSet.

=over 4

=item clearFeatureSet($gffConfRef)

The procedure removes all features from FeatureSet.

Parameters:
    $gffSetRef  - GFF feature set

Returns:
    Nothing
=back

=cut

sub clearFeatureSet(\%) {
    croak "ERROR: addSNPFeature\n" unless ( @_ == 7 );
    my ($gffSetRef) = @_;
    delete $gffSetRef->{featureSet};
}
1;    # says use was ok
__END__

