#!/usr/bin/perl -w
#########################################################
# Author:	Nancy F. Hansen
# Program:	report_lofreq_novel_variant_haplotypes.pl
# Function:	Read in annotated allele counts and
#               HMM state predictions and report haplotypes
#               where novel variants are observed.
# Date:		May 21, 2018
##########################################################

use strict;
use Getopt::Long;
use Pod::Usage;

use GTB::File qw(Open);

use vars qw($VERSION);

our %Opt;

my $Id = q$Id:$;
$VERSION = sprintf "%.4f", substr(q$Rev:10000$, 4)/10000;

my $Usage = "Usage: report_lofreq_novel_variant_haplotypes.pl <allele count file> <HMM results directory> <founder genotype_file> <unmapped in dm3 positions> <file of 6million Huang et al. vars> <alias file>\n";

process_commandline();

($#ARGV == 5)
    or die "$Usage";

my $ac_file = $ARGV[0];
my $hmm_dir = $ARGV[1];
my $genotype_file = $ARGV[2];
my $unmapped_in_dm3_bed = $ARGV[3];
my $six_million_huang_file = $ARGV[4];
my $alias_file = $ARGV[5];

# create sample/chrom hash with lists of [start, end, state] arrays:
my $rh_sample_founder_states = read_founder_haplotypes($hmm_dir);

# read in founder genotypes to look for missing gt sites with "0.00" founder fractions:
my $rh_founder_genotypes = read_founder_genotypes($genotype_file);

my $rh_unmappableindm3_positions = read_unmappable_positions($unmapped_in_dm3_bed);

# file of 6 million originally discovered DGRP vars:
my $rh_huang_six_million = read_huang_6million($six_million_huang_file);

# read alias file and build a hash of selection populations for each sample:
my $rh_selection_population = get_selection_pops($alias_file);

my $ac_fh = Open($ac_file);

my $ra_sample_ids;
my $rh_sample_names;
while (<$ac_fh>) { # examine observed alleles at each site
    chomp;
    my @fields = split /\t/, $_;
    if ($fields[0] eq 'Chrom') { # determine sample columns from header line
        ($ra_sample_ids, $rh_sample_names) = parse_ac_header($_);
        next;
    }

    my ($chrom, $pos, $ref_allele, $alt_allele, $linefrac, $lofreq_detections,
           $lofreq_score, $lofreq_strandbias) = @fields[0..7];

    @fields = @fields[8..$#fields];

    #next if ($linefrac ne 'NA') && ($linefrac != 0.00); # only interested in novel for now
    if ($Opt{filterqual} && $lofreq_score < $Opt{filterqual}) {
        print "$chrom:$pos:$alt_allele\tLOW_LOFREQ_SCORE\t$lofreq_score\t$linefrac\t$ref_allele\t$alt_allele\tNA\tNA\tNA\tNA\n";
        next;
    } # filter low scoring variants

    # divide sequenced samples into those with >=90%, >=10%, and <10% alternate allele in reads
    my %sample_alt_level = (); # sample name, then "Hom", "Het", or "Low";
    my %sample_founder_state = (); # sample name, then ref to list of founder states
    my %sample_founder_state_type = (); # sample name, then Hom or Het depending on founder state type, or NA if no founder state

    foreach my $sample_id (@{$ra_sample_ids}) {

        # determine HMM state for this sample at this position:
        my @founderstates = grep {$_->[0] <=$pos && $_->[1] >= $pos} @{$rh_sample_founder_states->{$rh_sample_names->{$sample_id}}->{$chrom}};
        my $state = $founderstates[0]->[2] || 'NA';
        my @founders = ($state =~ /RAL\_\d+/g);
        $sample_founder_state{$rh_sample_names->{$sample_id}} = [@founders];
        $sample_founder_state_type{$rh_sample_names->{$sample_id}} = (@founders < 2) ? 'NA' : (($founders[0] eq $founders[1]) ? 'Hom' : 'Het');

        my $total_reads = $fields[$sample_id];
        my $alt_reads = $fields[$sample_id - 1];
        my $homhet = ($alt_reads >= 0.9*$total_reads) ? 'Hom' : ($alt_reads >= 0.1*$total_reads) ? 'Het' : (($alt_reads) ? 'Low' : 'NA');
        $sample_alt_level{$rh_sample_names->{$sample_id}} = $homhet;
    }

    # assess overall explainability of data at this site:
    if (($linefrac ne 'NA') && ($linefrac == 0.000)) { # a genotyped founder site for which no lines showed the alt allele
        my $rh_founder_gts = $rh_founder_genotypes->{"$chrom:$pos:$alt_allele"};
        my @potential_variant_founders = ($rh_founder_gts) ?
                     grep { $rh_founder_gts->{$_} ne '0' } keys %{$rh_founder_gts} : ();
        my @unexplained_fstates = ();
        if (@potential_variant_founders) { # is everything consistent?
            foreach my $sample (keys %sample_alt_level) {
                my $level = $sample_alt_level{$sample};
                my @founder_states = @{$sample_founder_state{$sample}};
                if ($level eq 'Hom') { # should be only potential variant founders
                    foreach my $founder_state (@founder_states) {
                        if (!grep {$_ eq $founder_state} @potential_variant_founders) {
                            push @unexplained_fstates, $founder_state;
                        }
                    }
                }
                elsif ($level eq 'Het') { # should be at least one potential variant founder
                    my $unexplained = 1;
                    foreach my $founder_state (@founder_states) {
                        if (grep {$_ eq $founder_state} @potential_variant_founders) {
                            $unexplained = 0;
                            last;
                        }
                    }
                    if ($unexplained) { # this het state doesn't make sense
                        my $fstate = join '', @founder_states;
                        push @unexplained_fstates, $fstate;
                    }
                }
            }
        }
        my $unexplainedstring = (@unexplained_fstates) ? join ',', @unexplained_fstates : 'None';
        my $foundergtstring = (@potential_variant_founders) ? join ':', @potential_variant_founders : 'None';
        my @hom_samples = grep {$sample_alt_level{$_} eq 'Hom'} keys %sample_alt_level;
        my @het_samples = grep {$sample_alt_level{$_} eq 'Het'} keys %sample_alt_level;
        my @low_samples = grep {$sample_alt_level{$_} eq 'Low'} keys %sample_alt_level;

        my $hom_string = get_founderstate_string(\@hom_samples, \%sample_founder_state);
        my $het_string = get_founderstate_string(\@het_samples, \%sample_founder_state);
        my $low_string = get_founderstate_string(\@low_samples, \%sample_founder_state);

        if ($unexplainedstring eq 'None') {
            print "$chrom:$pos:$alt_allele\tDGRP_UNGENOTYPED_SNP\t$lofreq_score\t$linefrac\t$ref_allele\t$alt_allele\t$foundergtstring\tHOM:$hom_string\tHET:$het_string\tLOW:$low_string\tUNEXPLAINED:$unexplainedstring\n";
        }
        else {
            print "$chrom:$pos:$alt_allele\tUNEXPLAINED\t$lofreq_score\t$linefrac\t$ref_allele\t$alt_allele\t$foundergtstring\tHOM:$hom_string\tHET:$het_string\tLOW:$low_string\tUNEXPLAINED:$unexplainedstring\n";
        }
    }
    elsif ($linefrac eq 'NA') { # novel alt allele -- look for a single associated haplotype indicating a true de novo variant
        if ($rh_huang_six_million->{"$chrom:$pos:$alt_allele"}) {
            print "$chrom:$pos:$alt_allele\tDGRP_FILTERED_SNP\t$lofreq_score\t$linefrac\t$ref_allele\t$alt_allele\tNA\tNA\tNA\tNA\tFILTERED_HUANG_VARIANT\n";
            next;
        }
        if ($rh_unmappableindm3_positions->{$chrom}->{$pos}) {
            print "$chrom:$pos:$alt_allele\tUNMAPPED_IN_DM3\t$lofreq_score\t$linefrac\t$ref_allele\t$alt_allele\tNA\tNA\tNA\tNA\tUNMAPPABLE_IN_DM3\n";
            next;
        }

        # de novo presumably in a region that should have been discovered in Huang et al.
        # Is it always associated with a single founder? Or always with a single founder in a single population?
        my $foundergtstring = 'NA';
        
        my @hom_samples = grep {$sample_alt_level{$_} eq 'Hom'} keys %sample_alt_level;
        my @het_samples = grep {$sample_alt_level{$_} eq 'Het'} keys %sample_alt_level;
        my @low_samples = grep {$sample_alt_level{$_} eq 'Low'} keys %sample_alt_level;
        my %hom_counts = get_founderstate_counts(\@hom_samples, \%sample_founder_state);
        my $hom_string = get_founderstate_string(\@hom_samples, \%sample_founder_state);
        my $het_string = get_founderstate_string(\@het_samples, \%sample_founder_state);
        my $low_string = get_founderstate_string(\@low_samples, \%sample_founder_state);
        
        # Logic: If there are samples homozygous for this variant, they should all have
        # a single homozygous founder state
        
        my @observed_homstates = keys %hom_counts;
        my $associated_haplotype = 'NA';
        if ((@observed_homstates == 1) && $observed_homstates[0] =~ /^(RAL\_\d+)\1$/) { # only one homozygous founder state observed in homozygous samples--could be a true de novo
        
            my ($proposed_haplotype) = ($observed_homstates[0] =~ /^(RAL\_\d+)/) ? $1 : '';
            # check all samples for correct association:
            # Are there any samples with this founder haplotype that don't show the variant allele?
            my $agreeing_samples = 0;
            my @hap_problems = (); # samples that have the proposed founder hap but not enough variant allele
            my @variant_problems = (); # samples that have the variant but not the proposed founder hap
            my %mutated_pops = (); # keep track of which populations carry the founder hap and mutation
            my %unmutated_pops = (); # keep track of which populations carry the founder hap but not the mutation
            foreach my $sample (keys %sample_founder_state) {
                my $alt_level = $sample_alt_level{$sample};
                my @founder_haps = @{$sample_founder_state{$sample}};
                next if (!@founder_haps);
        
                my $no_haps = grep {$_ eq $proposed_haplotype} @founder_haps; # number of times proposed founder hap is predicted in this sample
                my $pop = $rh_selection_population->{$sample};
                if ($no_haps == 2) { # should be a Hom
                    if ($alt_level ne 'Hom') {
                        push @hap_problems, "$sample:$no_haps:$alt_level";
                        $unmutated_pops{$pop}++;
                    }
                    else {
                        $agreeing_samples++;
                        $mutated_pops{$pop}++;
                    }
                }
                elsif ($no_haps == 1) { # should be a Het
                    if ($alt_level ne 'Het') {
                        push @hap_problems, "$sample:$no_haps:$alt_level";
                        $unmutated_pops{$pop}++;
                    }
                    else {
                        $agreeing_samples++;
                        $mutated_pops{$pop}++;
                    }
                }
                elsif ($no_haps == 0) { # should be Low or NA
                    if ($alt_level eq 'Hom' || $alt_level eq 'Het') {
                        push @variant_problems, "$sample:$no_haps:$alt_level";
                        $mutated_pops{$pop}++;
                    }
                    else {
                        $agreeing_samples++;
                        $unmutated_pops{$pop}++;
                    }
                }
            }
            my $no_hap_problems = @hap_problems;
            my $assoc_hap_string;
            if (!@hap_problems && !@variant_problems) {
                $assoc_hap_string = "ASSOCHAP:$proposed_haplotype:$agreeing_samples";
            }
            elsif (@hap_problems) { # could be selected mutation in just one population
                my @varpops = keys %mutated_pops;
                my @unmutpops = keys %unmutated_pops;
                if (@varpops == 1 && @unmutpops) { # consistent with mutation in single selected population
                    my $unmutstring = join '/', @unmutpops;
                    $assoc_hap_string = "ASSOCPOP:$proposed_haplotype:$varpops[0],$unmutstring";
                }
                else {
                    $assoc_hap_string = "ASSOCHAPPROBS:$proposed_haplotype:$agreeing_samples/$no_hap_problems";
                }
            }
            else { # variant problems
                $assoc_hap_string = "ASSOCVARPROBS:$proposed_haplotype:$agreeing_samples/$no_hap_problems";
            }
            my $assoc_type = $assoc_hap_string;
            $assoc_type =~ s/:.*//;
            $assoc_type = ($assoc_type eq 'ASSOCHAP') ? 'DENOVO_SNP' :
                          ($assoc_type eq 'ASSOCPOP') ? 'SELECTED_DENOVO' : 'UNEXPLAINED';
            print "$chrom:$pos:$alt_allele\t$assoc_type\t$lofreq_score\t$linefrac\t$ref_allele\t$alt_allele\t$foundergtstring\tHOM:$hom_string\tHET:$het_string\tLOW:$low_string\t$assoc_hap_string\n";
        }
        else { # homozygous samples have either no founder states or many founder states or one heterozygous or empty state?
            my $obs_homstatestring = join '/', @observed_homstates;
            print "$chrom:$pos:$alt_allele\tUNEXPLAINED\t$lofreq_score\t$linefrac\t$ref_allele\t$alt_allele\tNA\tHOM:$hom_string\tHET:$het_string\tLOW:$low_string\tNOFOUNDERHAP\n";
        }
    }
    else { # DGRP polymorphism known to be present in at least one founder line
        print "$chrom:$pos:$alt_allele\tDGRP_SNP\t$lofreq_score\t$linefrac\n";
    }
}
close $ac_fh;

# parse the command line arguments to determine program options
sub process_commandline {
    
    # Set defaults here
    %Opt = (  filterqual => 0
           );
    GetOptions(\%Opt, qw(
                bwa novo outfile=s filterqual=i help+ version 
               )) || pod2usage(0);
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}); }
    if ($Opt{version}) { die "$0, ", q$Revision: $, "\n"; }

}

sub parse_ac_header {
    my $line = shift;

    my @sample_ids = ();
    my %sample_names = ();
    if ($line =~ /LoFreq_StrandBias\s(.*)$/) {
        my @all_fields = split /\s/, $1;
        my $propstring = ($Opt{'bwa'}) ? 'bwa.TotalReads' : ($Opt{'novo'}) ? 'novo.TotalReads' : 'TotalReads';
        @sample_ids = ($Opt{'bwa'}) ? grep {$all_fields[$_] =~ /\.bwa\.TotalReads/} (0..$#all_fields) :
                      ($Opt{'novo'}) ? grep {$all_fields[$_] =~ /\.novo\.TotalReads/} (0..$#all_fields) : 
                      grep {$all_fields[$_] =~ /TotalReads/} (0..$#all_fields);
        %sample_names = map {$all_fields[$_] =~ s/\.$propstring//; $_ => $all_fields[$_]} @sample_ids;
    }
    return ([@sample_ids], {%sample_names});
}

sub read_founder_haplotypes {
    my $dir = shift;

    opendir DIR, $dir
        or die "Couldn\'t open HMM directory $dir: $!\n";

    my @tdf_files = grep /\.tdf(\.gz){0,1}$/, readdir DIR;
    closedir DIR;

    my %founderstates = ();
    foreach my $tdf_file (@tdf_files) {
        next if ($Opt{bwa} && $tdf_file !~ /bwa/);
        next if ($Opt{novo} && $tdf_file !~ /novo/);

        my $tdf_fh = Open("$dir/$tdf_file");
        my ($this_sample, $state_start, $current_state, $last_pos, $last_chrom);
        while (<$tdf_fh>) {
            next if (/^Sample/);
            if (/^(\S+)\t(\S+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)\t/) {
                my ($sample, $chrom, $pos, $state, $prob, $homhet) = ($1, $2, $3, $4, $5, $6);
                $sample =~ s/\.bwa// if ($Opt{'bwa'});
                $sample =~ s/\.novo// if ($Opt{'novo'});

                if ($current_state && $state ne $current_state) { # state switch--print region
                    #print "$chrom\t$state_start\t$pos\t$sample\t$current_state\n";
                    push @{$founderstates{$sample}->{$chrom}}, [$state_start, $last_pos, $current_state];
                    $current_state = $state;
                    $state_start = $pos;
                }

                $current_state = $state if (!$current_state);
                $state_start = $pos if (!$state_start);
                $last_pos = $pos;
                $last_chrom = $chrom;
                $this_sample = $sample;
            }
        }
        close $tdf_fh;
        if ($current_state) {
            push @{$founderstates{$this_sample}->{$last_chrom}}, [$state_start, $last_pos, $current_state];
        }
    }

    return {%founderstates};
}

sub read_founder_genotypes {
    my $gt_file = shift;
    my $geno_fh = Open("$gt_file");

    my %founder_gts = ();
    my $no_fields;
    my @header_fields = ();
    my %founder_field_numbers = (); #where to find each sample's genotype
    while (<$geno_fh>) {
        chomp;
        my @fields = split /\s/, $_;
        if ($no_fields && $#fields != $no_fields-1) { # inconsistent file!
            die "Line has wrong number of fields (should be $no_fields)\n$_\n";
        }
        elsif (!$no_fields) {
            $no_fields = $#fields + 1;
        }
        # header
        if (!@header_fields) {
            @header_fields = @fields; 
            my %field_numbers = map {$header_fields[$_] => $_} (0..$#header_fields);
            foreach my $fieldname (keys %field_numbers) {
                if ($fieldname =~ /^line\_(\d+)$/) {
                    my $founder_field_number = $field_numbers{$fieldname};
                    $fieldname =~ s/line/RAL/;
                    $founder_field_numbers{$fieldname} = $founder_field_number;
                }
            }
            next;
        }
        else { # data!
            my $chrom = $fields[0];
            my $pos = $fields[1];
            my $allele = $fields[4];
            foreach my $founder (keys %founder_field_numbers) {
                $founder_gts{"$chrom:$pos:$allele"}->{$founder} = $fields[$founder_field_numbers{$founder}];
            }
        }
    }
    close $geno_fh;
    return {%founder_gts};
}

sub read_unmappable_positions {
    my $unmappable_file = shift;
    my $fh = Open("$unmappable_file");

    my %unmappable_positions = ();
    while (<$fh>) {
        chomp;
        my ($chrom, $pmo, $pos, $rest) = split /\s/, $_;
        if ($pos - $pmo == 1) { # SNP site
            $unmappable_positions{$chrom}->{$pos} = 1;
        }
    }
    close $fh;

    return {%unmappable_positions};
}

sub read_huang_6million {
    my $unfiltered_file = shift;
    my $unfiltered_fh = Open("$unfiltered_file");

    my %unfiltered_gts = ();
    my $no_fields;
    my @header_fields = ();
    my %founder_field_numbers = (); #where to find each sample's genotype
    while (<$unfiltered_fh>) {
        next if (/^##/); # header
        chomp;
        my @fields = split /\s/, $_;
        if ($no_fields && $#fields != $no_fields-1) { # inconsistent file!
            die "Line has wrong number of fields (should be $no_fields)\n$_\n";
        }
        elsif (!$no_fields) {
            $no_fields = $#fields + 1;
        }
        # header
        if ($fields[0] eq '#CHROM') {
            @header_fields = @fields; 
            my %field_numbers = map {$header_fields[$_] => $_} (0..$#header_fields);
            foreach my $fieldname (keys %field_numbers) {
                if ($fieldname =~ /^DGRP\-(\d+)$/) {
                    my $founder_field_number = $field_numbers{$fieldname};
                    $fieldname =~ s/DGRP\-/RAL_/;
                    $founder_field_numbers{$fieldname} = $founder_field_number;
                }
            }
            next;
        }
        else { # data!
            my $chrom = $fields[0];
            my $pos = $fields[1];
            my $allele = $fields[4];
            $unfiltered_gts{"$chrom:$pos:$allele"} = 1;
            #foreach my $founder (keys %founder_field_numbers) {
                #$unfiltered_gts{"$chrom:$pos:$allele"}->{$founder} = $fields[$founder_field_numbers{$founder}];
            #}
        }
    }
    close $unfiltered_fh;
    return {%unfiltered_gts};
}

sub get_selection_pops {
    my $alias_file = shift;

    my $alias_fh = Open("$alias_file");

    my %sample_pops = ();
    while (<$alias_fh>) {
        if (/^(\S+),(\S+)/) {
            my ($sample, $alias) = ($1, $2);
            my $pop = ($alias =~ /^SIP\_(\S+)\_/) ? $1 : 'NA';
            $sample_pops{$sample}=$pop;

            if ($pop eq 'NA') {
                print STDERR "Couldn\'t determine selected population for alias $alias\n";
            }
        }
    }
    close $alias_fh;

    return {%sample_pops};

}

sub get_founderstate_counts {
    my $ra_samples = shift;
    my $rh_founder_states = shift;

    my %state_counts = ();
    foreach my $sample (@{$ra_samples}) { 
        my $founder_state = join '', @{$rh_founder_states->{$sample}};
        $state_counts{$founder_state}++;
    }

    return %state_counts;
}

sub get_founderstate_string {
    my $ra_samples = shift;
    my $rh_founder_states = shift;

    my %state_counts = get_founderstate_counts($ra_samples, $rh_founder_states);
    my @state_strings = map {"$_:$state_counts{$_}"} keys %state_counts;
    my $state_string = (@state_strings) ? join ',', @state_strings : 'NA';

    return $state_string;
}

=pod

=head1 NAME

report_lofreq_novel_variant_haplotypes.pl - examine the founder haplotype predictions of "novel" lofreq variants.

=head1 SYNOPSIS

	report_lofreq_novel_variant_haplotypes.pl <file of allele counts> <directory of HMM founder state predictions>

Run report_lofreq_novel_variant_haplotypes.pl -man for a detailed description of options and the output files.

=head1 DESCRIPTION

=head1 OPTIONS

=over 5

=item B<--outfile>

This option specifies the location of the file to be written.

=back

=head1 OUTPUT

=head1 AUTHOR

 Nancy F. Hansen - nhansen@mail.nih.gov

=head1 LEGAL

This software/database is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software/database is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use. 

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S.  Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose. 

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
Genome Technology Branch" as the citation. 

=cut

__END__
