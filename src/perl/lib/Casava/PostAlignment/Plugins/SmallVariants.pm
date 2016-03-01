# =head1 LICENSE
#
# Copyright (c) 2007-2010 Illumina, Inc.
#
# This software is covered by the "Illumina Genome Analyzer Software
# License Agreement" and the "Illumina Source Code License Agreement",
# and certain third party copyright/licenses, and any user of this
# source file is bound by the terms therein (see accompanying files
# Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
# Illumina_Source_Code_License_Agreement.pdf and third party
# copyright/license notices).
#
# This file is part of the Consensus Assessment of Sequence And VAriation
# (CASAVA) software package.
#

=head1 SUMMARY

callSmallVariants - uses aligned reads and local de-novo assembly
contigs (from e.g. GROUPER) to locally realign reads and call snps and
small indels under a model which assumes a diploid sample.

=head1 DESCRIPTION

The small variant calling module can be used to call snps and small
indels, in addition to providing genotype calls and counts for every
site in the genome.

All snp and indel calls made by the module are given an associated
quality score which expresses the probability of the snp or indel
being present in the sample. Both models assume that the sample is a
single diploid individual.

The indel caller can detect indels directly from gapped alignments to
find relatively small indels (typically less than 20 bases). It can
also be run following the GROUPER local de-novo contig assembly
module. In this case the GROUPER contig and gapped alignment
information can be combined together to detect much larger indels.

The small variant caller will perform local read realigment around
indels to reduce noise during snp-calling and site count
summarization. The reads which are realigned by the small variant
caller may be optionally written out to separate BAM files for
diagnostic purposes. These realigned reads cannot be automatically
integrated back into the primary BAM files.

After snp and indel calling has completed, the small variants module
will generate summary tables and html reports for snps, indels, and
site coverage for the build.

B<IMPORTANT:> Note that in order for the I<callSmallVariants> target
to work, the I<sort> target must have been completed or be part of the
current workflow, and the I<assembleIndels> target must have been
completed for paired-end DNA workflows, unless this module is run with
the "--variantsSkipContigs" option, which causes the variant caller to
ignore any information provided from the indel assembler.

=head1 BASIC USAGE

The callSmallVariants module is designed to use the results of the
assembleIndels module if available, so a new paired-end build could be
run with the following minimum set of targets:

"--targets sort assembleIndels callSmallVariants"

...if assembleIndels (GROUPER) cannot be run, an alternative workflow is:

"--targets sort callSmallVariants --variantsSkipContigs"

=head1 ADVANCED_USAGE

To have the plugin provide a BAM file containing all reads which have
had their alignments altered during realignment, add the following to
the configuration command-line:

"--variantsWriteRealigned"

These reads will appear in the file "sorted.realigned.bam" in the
chromosome realigned bam directory.

=head1 RESULTS

The plugin provides depth and single position genotype call scores for
every mapped site in the reference genome. These calls can be found in
each bin directory in the gzipped file "sites.txt.gz", e.g.:

F<Project_Dir/Parsed_NN-NN-NN/chr1/0000/sites.txt.gz>

Note that this output can be omitted with the "--variantsNoSitesFiles"
option.

The snps for each reference sequence are aggregated and filtered
according to the --variantsSnpCovCutoff setting and summarized in the
chromosome-level file "snps.txt", e.g.:

F<Project_Dir/Parsed_NN-NN-NN/chr1/snps.txt>

The indels for each reference sequence are aggregated and filtered
according to the --variantsIndelCovCutoff setting and summarized in the
chromosome-level file "indels.txt", e.g.:

F<Project_Dir/Parsed_NN-NN-NN/chr1/indels.txt>

If any snps and indels are removed by the high-depth filter, they can
be found in their corresponding bin directory as, respectively:

F<Project_Dir/Parsed_NN-NN-NN/chr1/snps.removed.txt>
F<Project_Dir/Parsed_NN-NN-NN/chr1/indels.removed.txt>

When the "--variantsWriteRealigned" option is selected, there will
alse be a BAM file written to each reference sequence realigned bam
directory containing only those reads realigned by the variant
caller. The BAM filename is "sorted.realigned.bam", e.g.:

F<Project_Dir/Parsed_NN-NN-NN/chr1/bam/realigned/sorted.realigned.bam>

Statistics for coverage, as well as snp and indel calls for all
reference sequences are found in the 'stats' directory as:

F<Project_Dir/stats/coverage.summary.txt>

F<Project_Dir/stats/snps.summary.txt>

F<Project_Dir/stats/indels.summary.txt>

A summary of the same information is also available on the following
html pages:

F<Project_Dir/html/coverage.html>

F<Project_Dir/html/snps.html>

F<Project_Dir/html/indels.html>

To summarize the snps and indels in the stats and html directories
above, quality thresholds are used to select a subset of snps and
indels for summary reporting. The default thresholds are Q(snp) >= 20
and Q(indel) >= 20.  These values may be changed using the options
"--variantsSummaryMinQsnp" and "--variantsSummaryMinQindel".

=head1 DATA OUTPUT FORMATS

The details of the small variant caller output files are described in
the sections below. All of these files share a similar formating
structure. They are written as text files and composed of one header
segment followed by one data segment.

All lines in the header segment begin with the '#' character. Header
lines beginning with the sequence '#$' contain a key,value pair. The
reserved key 'COLUMNS' has an associated value containing the set of
column labels in the following data segment.

The data segment contains one entry per line, where each line is a set
of tab-delimited columns. Wherever appropriate, columns for sequence
name and position number are included such that the files are tabix
compatible.

=head1 SNP/SITE FILE CONTENTS

Snp and site files contain the same set of default columns. These are:

=over 4

=item 1. seq_name - chromosome/reference sequence name

=item 2. pos - position

=item 3. bcalls_used - basecalls used to make the genotype call for this site

=item 4. bcalls_filt - basecalls mapped to the site but filtered out before genotype calling

=item 5. ref - reference base

=item 6. Q(snp) - A Q-value expressing the probability of the homozygous reference genotype, subject to the expected rate of haplotype difference as expressed by the (Watterson) theta parameter

=item 7. max_gt - The most likely genotype (subject to theta, as above).

=item 8. Q(max_gt) - A Q-value expressing the probability that the genotype is not the most likely genotype above (subject to theta).

=item 9. max_gt|poly_site - The most likely genotype assuming this site is polymorphic with an expected allele frequency of 0.5 (theta is still used to calculate the probability of a third allele -- i.e. the chance of observing two non-reference alleles).

=item 10. Q(max_gt|poly_site) - A Q-value expressing the probability that the genotype is not the most likely genotype above.

=item 11. A_used - 'A' basecalls used

=item 12. C_used - 'C' basecalls used

=item 13. G_used - 'G' basecalls used

=item 14. T_used - 'T' basecalls used

=back

For many analyses where only likely polymorphic sites are being
examined (e.g. only sites in dbSNP), max_gt|poly_site should be used
instead of max_gt.

=head1 INDEL FILE CONTENTS

The indels.txt file follows the general variant caller output file
structure. The data segment of this file consists of 16 tab-delimited
fields. The fields are described in the table below (note that all
information is given with respect to the forward strand of the
reference sequence):

=over 4

=item 1. seq_name - chromosome/reference sequence name

=item 2. pos - start position of the indel (1)

=item 3. type - String summarizing the indel type. One of:

=over 8

=item nI - insertion of length n (e.g. 10I is a 10 base insertion)

=item nD - deletion of length n (e.g. 10D is a 10 base deletion)

=item BP_LEFT - left-side breakpoint

=item BP_RIGHT - right-side breakpoint

=item

=back

=item 4. ref_upstream - segment of the reference sequence 5' of the indel event. For right-side breakpoints this field is set to the value 'N/A'

=item 5. ref/indel - equal length sequences corresponding to the reference and indel cases which span the indel event. The character '-' indicates a gap in the reference or indel sequence.

=item 6. ref_downstream - segment of the reference sequence 3' of the indel event. For left-side breakpoints this field is set to the value 'N/A'.

=item 7. Q(indel) - phred scaled quality score of the indel (2)  By default the variant caller reports all indels with Q(indel) > 0.

=item 8. max_gtype - most probable indel genotype {het,hom,ref} (3)

=item 9. Q(max_gtype) - phred scaled quality score of the most probable indel genotype (2)

=item 10. depth - except for right-side breakpoints, this field reports the depth of the position preceding the left-most indel breakpoint. For right-side breakpoints this is the depth of the position following the breakpoint.

=item 11. alt_reads - number of reads strongly supporting either the reference path or an alternate indel path (4)

=item 12. indel_reads - number of reads strongly supporting the indel path (4)

=item 13. other_reads - number of reads intersecting the indel, but not strongly supporting either the reference or any one indel path (4)

=item 14. repeat_unit - the smallest repeating sequence unit within the inserted or deleted sequence. For breakpoints this field is set to the value 'N/A'.

=item 15. ref_repeat_count - number of times the repeat_unit sequence is contiguously repeated starting from the indel start position in the reference case.

=item 16. indel_repeat_count - number of times the repeat_unit sequence is contiguously repeated starting from the indel start position in the indel case.

=back

Footnotes:

(1) Except for right-side breakpoints, the reported start position of
the indel is the first (left-most) reference position following the
indel breakpoint. For right-side breakpoints the reported position is
the right-most position preceding the breakpoint. Also note that
wherever the same indel could be represented in a range of locations,
the caller attempts to report it in the left-most position possible.


(2) 'Q(indel)' refers to the Phred-scaled probability that this indel
does not exist at the given position. 'Q(max_gtype)' refers to
the Phred-scaled probability that the genotype of the indel is not
that given as 'max_gtype'. Phred-scaled Q-values are derived
from the corresponding probability P by the relationship Q =
-10log10(P). The Q-values given only reflect those error conditions
which can be represented in the indel calling model, which is not
comprehensive.

(3) The indel genotype categories are as follows: 'hom' refers
to a homozygous indel, 'het' refers to a heterozygous indel and
'ref' refers to no indel at this position. Note that these do
refer to true genotypes where indels overlap because the model is not
capable of jointly calling overlapping indels. In the case of
overlapping indels max_gtype refers to the most likely copy-number of
the indel. Note that indel calls where 'ref' is the most-likely
genotype will be reported. These correspond to indels with very low
Q(indel) values.

(4) For a read to strongly support either the reference or the indel
alignment, it must overlap an indel breakpoint by at least 6 bases and
the probability of the read's alignment following either the
reference or the indel path must be at least 0.999 (see Appendix for
more information).



=head1 OPTIONS

=head2 WORKFLOW OPTIONS

=over 4

=item --variantsSkipContigs

By default information from the GROUPER indel finder is used (and
required) in paired-end DNA-Seq analysis. This option disables use of
indel contigs during variant calling, and only uses gapped alignment
to find indels.

=item --variantsNoSitesFiles

Don't write out sites files.

=item --variantsNoReadTrim

By default, the ends of reads can be trimmed if the alignment path
through an indel is ambiguous. This option disables read trimming and
chooses the ungapped sequence alignment for any ambiguous read
segment. Note that this can trigger spurious SNP calls near indels.

=item --variantsWriteRealigned

Write only those reads which have been realigned to bam file:
'sorted.realigned.bam' for each reference sequence.

=back

=head2 READ MAPPING OPTIONS

=over 4

=item --variantsSEMapScoreRescue

Include reads if they have an SE mapping score equal to or above that set
by the "--QVCutoffSingle" option, even if the read pair fails the PE
mapping score threshold.

=item --variantsIncludeSingleton

Include paired-end reads which have unmapped mate reads. Note that
"--variantsSEMapScoreRescue" must also be specified because ELAND
gives singleton reads a PE mapping score of zero.

=item --variantsIncludeAnomalous

Include paired-end reads which have anomalous insert-size or
orientation. Note that "--variantsSEMapScoreRescue" must also be
specified because ELAND gives anomalous reads a PE mapping score of
zero.

=back

=begin comment

# =item --variantsMinPEMapScore=INTEGER

# Paired end reads must have at least this mapping score to be used
# (default is 90).

# =item --variantsMinSEMapScore=INTEGER

# Single ended reads must have at least this mapping score to be used
# (default is 90).

=end comment

=head2 SNP AND INDEL OPTIONS

=over 4

=item --variantsNoCovCutoff

Disables the snp and indel coverage filters detailed below for the
options: --variantsSnpCovCutoff and --variantsIndelCovCutoff. This
setting is recommended for targeted resequencing and RNA-Seq (Note it
is already set by default for RNA-Seq).

=back

=head2 SNP OPTIONS

=over 4

=item --variantsHetBias=FLOAT

Myrax bias factor, must be in range [0.,0.5), where 0. turns bias off
such that the results are equivilent to Hyrax.

=item --variantsSnpTheta=FLOAT

The frequency with which single base differences are expected between
two chromosomes (default is 0.001).

=item --variantsSnpCovCutoff=FLOAT

snps are filtered out of the final output if the used-depth is greater
than this value times the mean chromosomal used-depth. (default 3.0)

The filter may be disabled for targeted resequencing or other
applications by setting this value to -1 (or any negative number).

=item --variantsSnpCovCutoffAll

By default the mean chromosomal depth filter uses "used-depth"
calculated from all known sites in the reference sequence. When this
option is set, the threshold and the filtration use the full depth at
all known sites in the reference sequence.

=item --variantsMDFilterCount=INTEGER

The mismatch density filter removes all basecalls from consideration
during snp-calling where greater then 'variantsMDFilterCount'
mismatches to the reference occur on a read within a window of
1+2*'variantsMDFilterFlank' positions encompassing the current
position. The default value for 'variantsMDFilterCount' is 2 and for
'variantsMDFilterFlank is 20. Set either value to less than 0 to
disable the filter.

=item --variantsMDFilterFlank=INTEGER

(see description of mismatch density filter above)

=item --variantsIndependentErrorModel

This switch turns off all error dependency terms in the snp-calling
model, resulting in a simpler model where each basecall at a
site is treated as an independent observation.

=item --variantsMinQbasecall=INTEGER

The minimum basecall quality used for snp-calling. (default is 0).

=item --variantsSummaryMinQsnp=INTEGER

The snps.txt files contain all positions where Q(snp) > 0, however it
is expected that only a higher Q(snp) subset of these will be used
dependent upon the false positive tolerance of a user's workflow. For
this reason summary statistics about the called snps are created at a
higher "average-application" threshold, which can be set using this
option (default is 20).

=back

=head2 INDEL OPTIONS

=over 4

=item --variantsIndelTheta=FLOAT

The frequency with which indels are expected between two chromosomes
(default is 0.0001).

=item --variantsIndelCovCutoff=FLOAT

Indels are filtered out of the final output if the local sequence
depth is greater than this value times the mean chromosomal depth. The
sequence depth of the indel is approximated by the depth of the site
5' of the indel. (default 3.0)

The filter may be disabled for targeted resequencing or other
applications by setting this value to -1 (or any negative number).

=item --variantsCanIndelMin=INTEGER

Unless an indel is observed in at least this many gapped or GROUPER
reads, the indel cannot become a candidate for realignment and
genotype calling. (default: 3)

=item --variantsCanIndelMinFrac=FLOAT

Unless an indel is observed in at least this fraction of intersecting
reads, the indel cannot become a candidate for realignment and
genotype calling. (default: 0.02)

=item --variantsSmallCanIndelMinFrac=FLOAT

In addition to the above filter for all indels, for indels of size 4
or less, unless the indel is observed in at least this fraction of
intersecting reads, the indel cannot become a candidate for
realignment and genotype calling. (default: 0.1)

=item --variantsIndelErrorRate=FLOAT

Set the indel error rate used in the indel genotype caller to a
constant value of f (0<=f<=1).  The default indel error rate is taken
from an empirical function accounting for homopolymer length and indel
type (i.e. insertion or deletion). This flag overrides the default
behavior with a constant error rate for all indels.

=item --variantsSummaryMinQindel=INTEGER

The indels.txt files contain all positions where Q(indel) > 0, however
it is expected that only a higher Q(indel) subset of these will be
used dependent upon the false positive tolerance of a user's
workflow. For this reason summary statistics about the called snps are
created at a higher "average-application" threshold, which can be set
using this option (default is 20).

=item --variantsMaxIndelSize=INTEGER

Sets the maximum indel size for realignment and indel genotype
calling. Whenever an indel larger than this size is nominated by a
de-novo assembly contig it is handled as two independent
breakpoints. Note that increasing this value should lead to an
approximately linear increase in variant caller memory
consumption. The default value is 300 for paired-end builds and 50 for
single-end builds.

=back

=begin comment

=head2 RD OPTIONS

=over 4

=item --variantsConsensusVCF

Set consensusVCF (or genomeVCF) output mode. If not already defined,
this also sets 'variantsUsedAlleleCountMinQ' to 20 and switches on
'variantsPrintExtraSnpInfo'.

=item --variantsUsedAlleleCountMinQ=INTEGER

If printing used allele counts, only print basecalls with qscore greater than
this (default: no filtration)

=item --variantsPrintExtraSnpInfo

Append extra fields onto snp records used for filtration in gVCF output.
(default: off)

=item --variantsPrintAllGT

Print all polymorphic-site genotype probabilities in the sites and snps files

=item --variantsNoIndels

Don't call indels. Note that this refers to something entirely
different than not using GROUPER contigs (the non-RD option
--variantsSkipContigs above). Indels can be called when grouper isn't
run or used, and GROUPER can be run when indels aren't called. This
should probably not be a external option.

=item --no-variantsPrintUsedAlleleCounts

Don't print A,C,G,T used base counts for each allele in the sites and snps files

=item --variantsSnpMaxFilterFrac=FLOAT

If the fraction of basecalls at a site exceeds this value ([0-1]), then
don't call snps at this location. (default: 1.)

=item --variantsAddBacon

Run the BaCON snp and allele caller using the same realigned basecall
input provided to the default snp-caller (note this includes the
mapping score thresholds). This will create 'sort.count' and 'snp.txt'
files in each bin.

=item --variantsCanIndelMaxDensity=FLOAT

If the number of candidate indels per base exceeds this value for a
read then then the read realignment is curtailed to only allow
isolated indels. (default: 0.15)

=item --variantsSiteErrorDepDefault=FLOAT

The degree to which same-strand observations of an allele that
disagrees with a genotype are considered dependent (default is 0.6)

=item --variantsSiteErrorDepNoMismatch=FLOAT

The degree to which same-strand observations of an allele that
disagrees with a genotype are considered dependent when no neighboring
mismatches exist within the mismatch density filter window (default is
0.35)

=item --variantsMinVexp=FLOAT

Defines a limit on the dependent error penalty such that snp calling
is stable at higher depths (default is 0.25).

=item --variantsVexpIter=INTEGER

RD

=back

=end comment

=cut


package Casava::PostAlignment::Plugins::SmallVariants;

use warnings FATAL => 'all';
use strict;

use File::Spec;
use Getopt::Long;
use Sys::Hostname;

use lib '@CASAVA_FULL_PERL_LIBDIR@';
use Casava::Common::Log;
use Casava::Common::IOLib qw(createDirs executeCmd);
use Casava::PostAlignment::Sequencing::Config
  qw(loadConfiguration %CONF_APP %runsConfig isSpliceJunctionChrom);
use Casava::PostAlignment::Sequencing::SamLib qw(checkSamtoolsBin);

use Casava::TaskManager qw(executeSingleTask checkPointInit
               executeArrayTaskEx executeCheckPoint addTask2checkPoint tryAddTask2CheckPoint);
use Casava::PostAlignment::Plugins::RefSeq qw(getDoneCheckPoint);

#------------------------------------------------------------------------------
#
#=pod
#
#=head2 getTarget()
#
#All plugins are required to implement this method. The framework will execute
#the C<executeOrGenerateWorkflow> if user has requested the target to be added
#to the build.
#
#=over 4
#
#=item *
#
#Parameters:
#
#  No parameters
#
#=item *
#
#Returns:
#
#  name of target for which the plugin can perform workflow configuration
#
#=item *
#
#Exceptions
#
#  None
#
#=back
#
#=cut

sub getTarget { return 'callSmallVariants'; }

#------------------------------------------------------------------------------
#
#=pod
#
#=head2 getOptionsMapping
#
#Returns command-line options understood by plugin. Unless the option is one of
#the standard ones or supported by a plugin, the CASAVA command line will not
#accept it. All options supplied via command line are accessible through
#%CONF_PROJ hash
#
#I<Parameters:>
#
#=over 4
#
#=item *
#
#\%PARAMS
#
#reference to target hash
#
#=back
#
#I<Returns:>
#
#Mappings suitable for Getopt::Long::GetOptions
#
#I<Exceptions:>
#
#=over 4
#
#=item *
#
#None
#
#=back
#
#=cut

sub getOptionsMapping(\%) {
    my ( $PARAMS ) = @_;

    return (
            # workflow:
            'variantsWriteRealigned!' => \$PARAMS->{variantsWriteRealigned},
            'variantsSkipContigs!' => \$PARAMS->{variantsSkipContigs},
            'variantsNoSitesFiles!' => \$PARAMS->{variantsNoSitesFiles},
            'variantsNoReadTrim!' => \$PARAMS->{variantsNoReadTrim},

            # snps and indels
            'variantsNoCovCutoff!' => \$PARAMS->{variantsNoCovCutoff},

            # mapping:
            # 'variantsMinPEMapScore=i' => \$PARAMS->{variantsMinPEMapScore},
            # 'variantsMinSEMapScore=i' => \$PARAMS->{variantsMinSEMapScore},
            'variantsSEMapScoreRescue!' => \$PARAMS->{variantsSEMapScoreRescue},
            'variantsIncludeSingleton!' => \$PARAMS->{variantsIncludeSingleton},
            'variantsIncludeAnomalous!' => \$PARAMS->{variantsIncludeAnomalous},

            # snps:
            'variantsHetBias=f' => \$PARAMS->{variantsHetBias},
            'variantsSnpTheta=f' => \$PARAMS->{variantsSnpTheta},
            'variantsSnpCovCutoff=f' => \$PARAMS->{variantsSnpCovCutoff},
            'variantsSnpCovCutoffAll!' => \$PARAMS->{variantsSnpCovCutoffAll},
            'variantsMDFilterCount=i' => \$PARAMS->{variantsMDFilterCount},
            'variantsMDFilterFlank=i' => \$PARAMS->{variantsMDFilterFlank},
            'variantsSiteErrorDepDefault=f' => \$PARAMS->{variantsSiteErrorDepDefault},
            'variantsSiteErrorDepNoMismatch=f' => \$PARAMS->{variantsSiteErrorDepNoMismatch},
            'variantsMinVexp=f' => \$PARAMS->{variantsMinVexp},
            'variantsIndependentErrorModel!' => \$PARAMS->{variantsIndependentErrorModel},
            'variantsMinQbasecall=i' => \$PARAMS->{variantsMinQbasecall},
            'variantsSummaryMinQsnp=i' =>\$PARAMS->{variantsSummaryMinQsnp},

            # indels:
            'variantsIndelTheta=f' => \$PARAMS->{variantsIndelTheta},
            'variantsIndelCovCutoff=f' => \$PARAMS->{variantsIndelCovCutoff},
            'variantsCanIndelMin=i' => \$PARAMS->{variantsCanIndelMin},
            'variantsCanIndelMinFrac=f' => \$PARAMS->{variantsCanIndelMinFrac},
            'variantsSmallCanIndelMinFrac=f' => \$PARAMS->{variantsSmallCanIndelMinFrac},
            'variantsIndelErrorRate=f' => \$PARAMS->{variantsIndelErrorRate},
            'variantsSummaryMinQindel=i' =>\$PARAMS->{variantsSummaryMinQindel},
            'variantsMaxIndelSize=i' =>\$PARAMS->{variantsMaxIndelSize},

            # rd:
            'variantsConsensusVCF!' => \$PARAMS->{variantsConsensusVCF},
            'variantsUsedAlleleCountMinQ=i' => \$PARAMS->{variantsUsedAlleleCountMinQ},
            'variantsPrintExtraSnpInfo!' => \$PARAMS->{variantsPrintExtraSnpInfo},
            'variantsPrintAllGT!' => \$PARAMS->{variantsPrintAllGT},
            'variantsPrintUsedAlleleCounts!' => \$PARAMS->{variantsPrintUsedAlleleCounts},
            'variantsNoIndels!' => \$PARAMS->{variantsNoIndels},
            'variantsSnpMaxFilterFrac=f' => \$PARAMS->{variantsSnpMaxFilterFrac},
            'variantsAddBacon!' => \$PARAMS->{variantsAddBacon},
            'variantsCanIndelMaxDensity=f' => \$PARAMS->{variantsCanIndelMaxDensity},
            'variantsVexpIter=i' => \$PARAMS->{variantsVexpIter});
}

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#
#=pod
#
#=head2 executeOrGenerateWorkflow
#
#All plugins are required to implement this method. The framework calls it when
#the target returned by L</getTarget> has been requested to be added to the
#CASAVA build.
#
#I<Parameters:>
#
#=over 4
#
#=item *
#
#$prevTarget
#
#if C<start>, the target is the first target in the workflow. Otherwise, plugin
#is expected to specify dependencies on the other target tasks (if needed).
#
#=item *
#
#$projectDir
#
#Folder containing chromosome folders.
#
#=item *
#
#$CONF_PROJ_Ref
#
#HASH MAP Ref to project configuration
#
#=item *
#
#$chromsBinSizesRef
#
#HASH MAP Ref to chromosome sizes
#
#=back
#
#I<Returns:>
#
#Nothing
#
#I<Exceptions:>
#
#=over 4
#
#=item *
#
#Incorrect analysis type - indels are possible only for paired DNA
#
#=back
#
#=cut

sub executeOrGenerateWorkflow($$\%\%)
{
    my ($prevTarget, $projectDir, $CONF_PROJ_Ref, $chromsBinSizesRef) = @_;

    my $readMode            = $CONF_PROJ_Ref->{readMode};
    my $dirCurrentBuild     = $CONF_PROJ_Ref->{dirBuildParsed};
    my $isSkipContigs       = ((defined $CONF_PROJ_Ref->{variantsSkipContigs}) and $CONF_PROJ_Ref->{variantsSkipContigs});
    $isSkipContigs          = ($isSkipContigs or ($readMode ne 'paired'));

    my $applicationType     = $CONF_PROJ_Ref->{applicationType};

    my $libexecDir          = '@CASAVA_FULL_LIBEXECDIR@';
    my $scriptDir           = $libexecDir;
    my $variantsBin         = File::Spec->catfile($scriptDir, 'callSmallVariants.pl');
    my $filterBin           = File::Spec->catfile($scriptDir, 'filterSmallVariants.pl');
    my $mergeRealignedBin   = File::Spec->catfile($scriptDir, 'mergeRealignedBam.pl');
    my $statsBin            = File::Spec->catfile($scriptDir , "generateVariantsStats.pl");
    my $reportsBin          = File::Spec->catfile($scriptDir , "generateVariantsReports.pl");

    my $isWriteRealigned  = ((defined $CONF_PROJ_Ref->{variantsWriteRealigned}) and $CONF_PROJ_Ref->{variantsWriteRealigned});

    my $taskIdRefSeq = Casava::PostAlignment::Plugins::RefSeq::getDoneCheckPoint();

    my $checkPointName = "VARIANTS";
    my $checkPointRef  = checkPointInit($checkPointName , "Small variant caller sequence-level processing" );

    my @chroms = keys(%{$chromsBinSizesRef});

    for my $chrom (@chroms)
    {
        next if ( isSpliceJunctionChrom($chrom, %{$CONF_PROJ_Ref}) );
        my $checkPointChromName = 'N/A';

        {
            $checkPointChromName = "VARIANTS_$chrom";
            my $checkPointChromRef  = checkPointInit($checkPointChromName , "Small variant caller bin-level processing for chrom $chrom");

            my $binCount = $chromsBinSizesRef->{$chrom};
            for ( my $i = 0 ; $i < $binCount ; $i++ ) {
                my $binId = sprintf "%04d", $i;
                my $cmd = "$variantsBin --projectDir=$projectDir --chrom='$chrom' --binId=$binId";
                if ($CONF_PROJ_Ref->{alignSortedBam})
                {
                    $cmd .= " --bamFile='$CONF_PROJ_Ref->{alignSortedBam}'";
                }

                # enumerate all prerequisite tasks for this bin:
                #
                my $checkPointNamePrereq = "PREREQ_VARIANTS_$chrom/$binId";
                my $checkPointPrereq  = checkPointInit( $checkPointNamePrereq, "Variants prereq for $chrom/$binId" );

                addTask2checkPoint( $taskIdRefSeq , %{$checkPointPrereq} );

                # dependent on whole BAM file now:
                my $taskIdSplitSameBinPrereq = "SYNCH_MERGE_$chrom\_DONE";
                tryAddTask2CheckPoint( $taskIdSplitSameBinPrereq, %{$checkPointPrereq} );

                # my $taskIdSplitSameBinPrereq = "build$chrom/$binId";
                # tryAddTask2CheckPoint( $taskIdSplitSameBinPrereq, %{$checkPointPrereq} );
                # if (0 != $i) {
                #     my $prevBinId = sprintf "%04d", ($i - 1);
                #     my $taskIdSplitPrevBinPrereq = "build$chrom/$prevBinId";
                #     tryAddTask2CheckPoint( $taskIdSplitPrevBinPrereq, %{$checkPointPrereq} );
                # }
                if(not $isSkipContigs) {
                    my $checkPointIndelChromName = "GROUPER_$chrom/$binId";
                    my $ret = tryAddTask2CheckPoint( $checkPointIndelChromName, %{$checkPointPrereq} );
                    if( not $ret) {
                        $checkPointIndelChromName = "GROUPER_$chrom";
                        tryAddTask2CheckPoint( $checkPointIndelChromName, %{$checkPointPrereq} );
                    }
                }
                executeCheckPoint( %{$checkPointPrereq} );

                my $taskId = "VARIANTS_$chrom/$binId";
                executeArrayTaskEx( $cmd, $taskId, "Computing small variants for $chrom $binId", $checkPointNamePrereq, %{$checkPointChromRef} );
            }

            executeCheckPoint( %$checkPointChromRef );
        }
        my $cmd = "$filterBin --projectDir=$projectDir --chrom='$chrom'";
        executeArrayTaskEx( $cmd, "VARIANTS_FILTER_$chrom", "Aggregating small variant calls for $chrom", $checkPointChromName , %$checkPointRef );

        if($isWriteRealigned)
        {
            # \TODO a second checkpoint to indicate that all
            # realignments were merged would be slightly better, as
            # this would allow stats/html to proceed without waiting
            # for this:
            my $cmd = "$mergeRealignedBin --projectDir=$projectDir --chrom='$chrom'";
            executeArrayTaskEx( $cmd, "VARIANTS_MERGE_REALIGN_$chrom", "Merging realigned BAM files for $chrom", $checkPointChromName , %$checkPointRef );
        }
    }
    executeCheckPoint( %$checkPointRef );

    ## Generate stats table and html report for variants:
    my $taskId = "VARIANTS_STATS";
    my $statsCmd = "$statsBin --projectDir=$projectDir";
    executeSingleTask( $statsCmd , $taskId, "Small variants stats", $checkPointName );

    # use "REPORT" taskId to prevent reports from being run in
    # parallel (this makes html page registration safer)
    my $reportCmd = "$reportsBin --projectDir=$projectDir";
    executeSingleTask( $reportCmd , "REPORT", "Small variants report", $taskId );
}


#------------------------------------------------------------------------------

1;
__END__

# =pod
#
# =head1 CONFIGURATION AND ENVIRONMENT
#
# Uses CASAVA build configuration data.
#
# =head1 DEPENDENCIES
#
# =over 4
#
# =item Standard perl modules
#
# strict, warnings
#
# =item External perl modules
#
# =item Casava perl modules
#
# =head1 BUGS AND LIMITATIONS
#
# Please report problems to Illumina Technical Support (support@illumina.com)
# 
# Patches are welcome.
#
# =head1 AUTHOR
#
# Chris Saunders
#
# =cut

#------------------------------------------------------------------------------
