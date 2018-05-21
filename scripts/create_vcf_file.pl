#!/usr/bin/perl -w
#########################################################
# Author:	Nancy F. Hansen
# Program:	create_vcf_file.pl
# Function:	Read in annotated allele counts and
#               output them in VCF format
# Date:		May 18, 2018
##########################################################

use strict;
use Getopt::Long;
use Pod::Usage;

use GTB::File qw(Open);
use GTB::File::SortedIter;

use NISC::Sequencing::Date;

use vars qw($VERSION);

our %Opt;

my $Id = q$Id:$;
$VERSION = sprintf "%.4f", substr(q$Rev:10000$, 4)/10000;

my $Usage = "Usage: create_vcf_file.pl <allele count file>\n";

process_commandline();

($#ARGV == 0)
    or die "$Usage";

my $ac_file = $ARGV[0];
my $output_file = $Opt{'outfile'};

if (!$output_file) {
    $output_file = $ac_file;
    $output_file =~ s:.*/::;
    $output_file =~ s:\.gz$::;
    $output_file .= ".vcf.gz";
}

my $out_fh = Open($output_file, "w");
print_vcf_header($out_fh);

my $ac_fh = Open($ac_file);

my @sample_ids = ();
my @sample_names;
while (<$ac_fh>) {
    chomp;
    my @fields = split /\t/, $_;
    if ($fields[0] eq 'Chrom') {
        if (/LoFreq_StrandBias\s(.*)$/) {
            my @all_fields = split /\s/, $1;
            my $propstring = ($Opt{'bwa'}) ? 'bwa.TotalReads' : ($Opt{'novo'}) ? 'novo.TotalReads' : 'TotalReads';
            @sample_ids = ($Opt{'bwa'}) ? grep {$all_fields[$_] =~ /\.bwa\.TotalReads/} (0..$#all_fields) :
                          ($Opt{'novo'}) ? grep {$all_fields[$_] =~ /\.novo\.TotalReads/} (0..$#all_fields) : 
                          grep {$all_fields[$_] =~ /TotalReads/} (0..$#all_fields);
            @sample_names = map {$all_fields[$_] =~ /^(.*)\.$propstring/ ? $1 : 'NA'} @sample_ids;
        }
        print $out_fh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

        foreach my $samplename (@sample_names) {
            print $out_fh "\t$samplename";
        }
        print $out_fh "\n";

        next;
    }

    my $chrom = shift @fields;
    my $pos = shift @fields;
    my $ref_allele = shift @fields;
    my $alt_allele = shift @fields;
    my $linefrac = shift @fields;
    my $lofreq_detections = shift @fields;
    my $lofreq_score = shift @fields;
    my $lofreq_strandbias = shift @fields;

    #my @sample_fields = map {$fields[$_]} @sample_ids;
    #my @detected_ids = grep {my @limits = ($sample_fields[$_] eq 'NA') ? (0, 0) : split /:/, $sample_fields[$_]; $limits[0]>0 } (0..$#sample_ids);
    #my $sample_string = join ':', @detected_ids;
    my $no_samples = @sample_names;

    my @info_fields = ("LOFREQDETECTIONS=$lofreq_detections",
                       "LOFREQSTRANDBIAS=$lofreq_strandbias",
                       "FOUNDERLINEFRACTION=$linefrac");
    my $info_field = (@info_fields) ? join ';', @info_fields : '.';
    my $filter_field = ($lofreq_score >= 1000) ? 'PASS' : 'LowQual';
    print $out_fh "$chrom\t$pos\t.\t$ref_allele\t$alt_allele\t$lofreq_score\t$filter_field\t$info_field\tDP:AC";
    foreach my $sample_id (@sample_ids) {
        my $total_reads = $fields[$sample_id];
        my $alt_reads = $fields[$sample_id - 1];
        print $out_fh "\t$total_reads:$alt_reads";
    }
    print $out_fh "\n";
}
close $ac_fh;

# parse the command line arguments to determine program options
sub process_commandline {
    
    # Set defaults here
    %Opt = ( 
           );
    GetOptions(\%Opt, qw(
                bwa novo outfile=s help+ version 
               )) || pod2usage(0);
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}); }
    if ($Opt{version}) { die "$0, ", q$Revision: $, "\n"; }

}

sub print_vcf_header {
    my $fh = shift;
    print $fh "##fileformat=VCFv4.3\n";
    my $date_obj = NISC::Sequencing::Date->new(-plain_language => 'today');
    my $year = $date_obj->year();
    my $month = $date_obj->month();
    $month =~ s/^(\d)$/0$1/;
    my $day = $date_obj->day();
    $day =~ s/^(\d)$/0$1/;
    print $fh "##fileDate=$year$month$day\n";
    print $fh "##reference=$Opt{refname}\n" if ($Opt{refname});
    print $fh "##INFO=<ID=LOFREQDETECTIONS,Number=1,Type=Integer,Description=\"Number of lines allele was detected in by LoFreq\">\n";
    print $fh "##INFO=<ID=LOFREQSTRANDBIAS,Number=1,Type=Float,Description=\"Maximum strand bias of LoFreq detections\">\n";
    print $fh "##INFO=<ID=FOUNDERLINEFRACTION,Number=1,Type=Float,Description=\"Fraction of genotyped founders with alternate allele\">\n";
    print $fh "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total depth of coverage\">\n";
    print $fh "##FORMAT=<ID=AC,Number=A,Type=Integer,Description=\"Alternate allele count\">\n";

}

=pod

=head1 NAME

create_vcf_file.pl - write a VCF file from a file of allele counts

=head1 SYNOPSIS

	create_vcf_file.pl <file of allele counts> --outfile <location to write output>

Run create_vcf_file.pl -man for a detailed description of options and the output files.

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
