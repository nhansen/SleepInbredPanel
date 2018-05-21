#!/usr/bin/perl -w
#########################################################
# Author:	Nancy F. Hansen
# Program:	gather_annotations.pl
# Function:	Read in the DGRP line genotype file
#               and LoFreq VCFs and write annotations
#               for Susan's counts file.
# Date:		September 15, 2014
##########################################################

use strict;
use Getopt::Long;
use Pod::Usage;

use GTB::File qw(Open);
use GTB::File::SortedIter;

use vars qw($VERSION);

our %Opt;

my $Id = q$Id:$;
$VERSION = sprintf "%.4f", substr(q$Rev:10000$, 4)/10000;

my $Usage = "Usage: gather_annotations.pl --sites < BED file of SNP sites with alternate alleles> --dgrp <tab-delimited DGRP site file> --vcfs <file of vcf filenamess> --outfile <path to output file>\nFor more information, type \"perldoc gather_annotations.pl.";

process_commandline();

($#ARGV < 0)
    or die "$Usage";

my $sites_file = $Opt{'sites'};
my $dgrp_file = $Opt{'dgrp'};
my $vcf_fof = $Opt{'vcfs'};
my $output_file = $Opt{'outfile'};

if (!$dgrp_file || !$vcf_fof || !$output_file) {
    die $Usage;
}

print STDERR "File of sites to include: $sites_file\n";
print STDERR "DGRP file of alleles: $dgrp_file\n";
print STDERR "File of LoFreq VCF files: $vcf_fof\n";
print STDERR "Writing annotations to: $output_file\n";

my $rh_sites = read_site_file($sites_file);
my $rh_dgrp_counts = read_dgrp_file($dgrp_file);
my $rh_vcf_data = read_lofreq_vcfs($vcf_fof, $rh_sites);

print_annotation_file($output_file, $rh_sites, $rh_dgrp_counts, $rh_vcf_data);

# parse the command line arguments to determine program options
sub process_commandline {
    
    # Set defaults here
    %Opt = ( 
           );
    GetOptions(\%Opt, qw(
                sites=s dgrp=s vcfs=s outfile=s help+ version
               )) || pod2usage(0);
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}); }
    if ($Opt{version}) { die "$0, ", q$Revision: $, "\n"; }

}

sub read_site_file {
    my $file = shift;

    my $sites_fh = Open($file);
    my %all_sites = ();
    while (<$sites_fh>) {
        chomp;
        if (/^(\S+)\s(\d+)\s(\d+)/) { # will be all SNPs for now
            my ($chr, $posmo, $pos, $ref_allele, $var_allele) = split /\s/, $_;
            $ref_allele = uc($ref_allele);
            $var_allele = uc($var_allele);
            $all_sites{$chr}->{$pos}->{$ref_allele}->{$var_allele}++;
        }
        elsif (!(/^#/)) {
            print "Skipping unorthodox BED file line:\n$_\n";
        }
    }
    close $sites_fh;
    return {%all_sites};
}

sub read_dgrp_file {
    my $file = shift;

    my $dgrp_fh = Open($file);
    my %all_sites = ();
    while (<$dgrp_fh>) {
        chomp;
        if (/^(\S+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)/) { # will be all SNPs for now
            my ($chr, $pos, $type, $ref, $var, @other_fields) = split /\t/, $_;
            $ref = uc($ref);
            $var = uc($var);
            next if ($type !~ /SNP/);
            if (@other_fields == 14) {
                for (my $i=0; $i<4; $i++) {
                    shift @other_fields;
                }
            }
            if (@other_fields != 10) { # not enough genotypes!
                die "Not enough DGRP genotypes for position $chr:$pos ($type $ref/$var) in $file!\n";
            }
            my $variant_lines = grep {$_ eq '2'} @other_fields;
            my $ref_lines = grep {$_ eq '0'} @other_fields;
            my $frac = 'NA';
            if (!($variant_lines+$ref_lines)) {
                print "No genotypes for site $chr:$pos ($type $ref/$var) in $file!\n";
            }
            else {
                $frac = $variant_lines/($ref_lines + $variant_lines);
                $frac = sprintf("%4.3f", $frac);
            }
            $all_sites{$chr}->{$pos}->{$ref}->{$var} = {'LineFraction' => $frac};
            #print STDERR "Recording $frac for $chr:$pos:$ref:$var\n";
        }
        elsif (!(/^#/)) {
            print "Skipping oddly formatted DGRP file line:\n$_\n";
        }
    }
    close $dgrp_fh;
    return {%all_sites};
}

sub read_lofreq_vcfs {
    my $file_of_files = shift;
    my $rh_sites = shift;

    my $fof_fh = Open($file_of_files);
    my @all_files = ();
    my @iters = ();
    my @cols = qw( CHROM POS ID REF ALT QUAL FILTER INFO );

    while (<$fof_fh>) {
        chomp;
        if (/^(\S+)\s(\S+)/) {
            my ($sample, $vcffile) = ($1, $2);
            if (-r $vcffile) {
                push @all_files, {'sample' => $sample, 'vcf' => $vcffile };
                my $fh = Open($vcffile);
                while (<$fh>) {
                    last if (/^#CHROM/);
                }
                    
                push @iters, GTB::File::SortedIter->new( 
                         fh => $fh,
                         cols => \@cols,
                         skip => '/^#/',
                         key => sub { my $ra = shift->current_array() || return;
                                      my $chr = $ra->[0];
                                      my $pos = $ra->[1];
                                      my $var = $ra->[4];
                                      $var = uc($var);
                                      return sprintf "%s\t%09d\t%s", $chr, $pos, $var},
                         );
            }
            else {
                warn "Skipping file $vcffile--not readable!\n";
            }
        }
    }
    close $fof_fh;

    my %all_sites = ();
    foreach my $chr (sort keys %{$rh_sites}) {
        foreach my $pos (sort {$a <=> $b} keys %{$rh_sites->{$chr}}) {
            foreach my $ref (sort keys %{$rh_sites->{$chr}->{$pos}}) {
                $ref = uc($ref);
                foreach my $var (sort keys %{$rh_sites->{$chr}->{$pos}->{$ref}}) {
                    $var = uc($var);
                    my $pos_string = sprintf "%s\t%09d\t%s", $chr, $pos, $var;
                    my ($max_qual, $max_sb) = (0, 0); 
                    my %samples_with_pass = ();
                    for (my $i=0; $i<=$#iters; $i++) {
                        my $iter_obj = $iters[$i];
                        my $vcffile = $all_files[$i]->{vcf};
                        my $sample = $all_files[$i]->{sample};
                        # attempt to get the file's line for this position:
                        my $file_pos = $iter_obj->current_key();
                        #print "$vcffile:$file_pos compared to $pos_string\n";
                        while (($file_pos) && ($file_pos lt $pos_string)) {
                            #print "$vcffile:$file_pos (looking for $pos_string)\n";
                            $file_pos = $iter_obj->next_key();
                        }
       
                        if (($file_pos) && ($file_pos eq $pos_string)) { # found a record!
                            #print "Found $file_pos in $vcffile!\n";
                            my @fields = $iter_obj->current_array();
                            my $file_chr = $fields[0];
                            my $file_pos = $fields[1];
                            if ($file_chr ne $chr || $file_pos != $pos) {
                                warn "Our iteration has returned the wrong position: $file_chr:$file_pos for $chr:$pos in file $vcffile!\n";
                            }
                            my $qual = $fields[5];
                            if ($qual =~ /^(\-{0,1}\d+)$/ && $qual > $max_qual) {
                                $max_qual = $qual;
                            }
                            my $filter = $fields[6];
                            if ($filter eq 'PASS') {
                                $samples_with_pass{$sample}++;
                            }
                            my $info = $fields[7];
                            my $sb = ($info =~ /SB=(\d+)/) ? $1 : 'NA';
                            if ($sb =~ /^(\-{0,1}\d+)$/ && $sb > $max_sb) {
                                $max_sb = $sb;
                            }
                        }
                    }
                    my $files_with_pass = keys %samples_with_pass;
                    $all_sites{$chr}->{$pos}->{$ref}->{$var} = {'LoFreq_Detections' => $files_with_pass,
                                                        'LoFreq_MaxQual' => $max_qual,
                                                        'LoFreq_StrandBias' => $max_sb };
                    #print "LoFreq\t$chr\t$pos\t$var\t$files_with_pass\t$max_qual\t$max_sb\n";
                }
            }
        }
    }

    return {%all_sites};
}

sub print_annotation_file {
    my $output_file = shift;
    my $rh_sites = shift;
    my $rh_dgrp = shift;
    my $rh_vcf_data = shift;

    my $out_fh = Open($output_file, "w");

    my @fields = qw( LineFraction LoFreq_Detections LoFreq_MaxQual LoFreq_StrandBias );
    my $header_string = join "\t", @fields;
    $header_string = "Chrom\tPos\tRefBase\tVarBase\t$header_string\n";

    print $out_fh $header_string;

    foreach my $chr (keys %{$rh_sites}) {
        foreach my $pos (sort {$a <=> $b} keys %{$rh_sites->{$chr}}) {
            foreach my $ref (sort keys %{$rh_sites->{$chr}->{$pos}}) {
                foreach my $var (sort keys %{$rh_sites->{$chr}->{$pos}->{$ref}}) {
                    my @var_fields = ( $chr, $pos, $ref, $var );
                    my $rh_dgrp_data = ($rh_dgrp && $rh_dgrp->{$chr} && $rh_dgrp->{$chr}->{$pos} &&
                             $rh_dgrp->{$chr}->{$pos}->{$ref} && $rh_dgrp->{$chr}->{$pos}->{$ref}->{$var}) ?
                             $rh_dgrp->{$chr}->{$pos}->{$ref}->{$var} : {};
                    push @var_fields, (defined ($rh_dgrp_data->{'LineFraction'})) ? $rh_dgrp_data->{'LineFraction'} : 'NA';
                    my $rh_lofreq_data = ($rh_vcf_data && $rh_vcf_data->{$chr} && $rh_vcf_data->{$chr}->{$pos} &&
                             $rh_vcf_data->{$chr}->{$pos}->{$ref} && $rh_vcf_data->{$chr}->{$pos}->{$ref}->{$var}) ? 
                             $rh_vcf_data->{$chr}->{$pos}->{$ref}->{$var} : {};
                    push @var_fields, (defined ($rh_lofreq_data->{'LoFreq_Detections'})) ? $rh_lofreq_data->{'LoFreq_Detections'} : 'NA';
                    push @var_fields, (defined ($rh_lofreq_data->{'LoFreq_MaxQual'})) ? $rh_lofreq_data->{'LoFreq_MaxQual'} : 'NA';
                    push @var_fields, (defined ($rh_lofreq_data->{'LoFreq_StrandBias'})) ? $rh_lofreq_data->{'LoFreq_StrandBias'} : 'NA';
                    my $var_string = join "\t", @var_fields;
                    print $out_fh "$var_string\n";
                }
            }
        }
    }
    close $out_fh;
}

=pod

=head1 NAME

gather_annotations.pl - read in a set of SNP sites with variant alleles, then a DGRP file of inbred line genotypes, then the LoFreq VCF files for all samples, and generate an annotation file to be used with "combine_allele_counts.pl".

=head1 SYNOPSIS

	gather_annotations.pl --sites <bed file with site locations> --dgrp <file with DGRP line information for each polymorphic site> --vcfs <file of LoFreq vcf files> --outfile <location to write output>

Run gather_annotations.pl -man for a detailed description of options and the output files.

=head1 DESCRIPTION

=head1 OPTIONS

=over 5

=item B<--sites>

This option specifies the location of a bedfile (assumed for now to be single base sites of SNPs, with an additional field containing the alternate alleles).

=item B<--dgrp>

This option specifies the location of a tab-delimited file with one record for each polymorphic site in the DGRP inbred line set (see "freeze2_informative_polymorphisms_corrected.tdf.txt" in the "project_data" directory).

=item B<--vcfs>

This option specifies the location of a file of file locations for LoFreq VCF files--two fields, sample name, then VCF file path.

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
