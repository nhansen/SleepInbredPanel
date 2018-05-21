#!/usr/bin/perl -w
#########################################################
# Author:	Nancy F. Hansen
# Program:	combine_allele_counts.pl
# Function:	Read in allele counts (output of 
#               bamcounts and allele_proportions.R)
#               for single samples, and output a large
#               table with counts for all samples.
# Date:		September 10, 2014
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

my $Usage = "Usage: combine_allele_counts.pl --ci_files <file of filenames,header strings> --sites <BED file with sites to include> --annotations <annot_file> --outfile <path to output file>\nFor more information, type \"perldoc combine_allele_counts.pl.";

process_commandline();

($#ARGV < 0)
    or die "$Usage";

my $ci_fof = $Opt{'ci_files'};
my $site_bedfile = $Opt{'sites'};
my $annotations_file = $Opt{'annotations'}; # optional
my $output_file = $Opt{'outfile'};

if (!$ci_fof || !$site_bedfile) {
    die $Usage;
}

if (!$output_file) {
    $output_file = $ci_fof;
    $output_file =~ s:.*/::;
    $output_file .= ".allcounts.out";
}

print STDERR "CI sample file of files: $ci_fof\n";
print STDERR "BED file of sites: $site_bedfile\n";
print STDERR "Writing counts to: $output_file\n";

my $rh_count_sites = read_site_file($site_bedfile);
my ($annot_header, $rh_annotations) = read_annotations_file($annotations_file) if ($annotations_file);
print STDERR "Read site file $site_bedfile for sites!\n";
read_and_write_counts($ci_fof, $rh_count_sites, $annot_header, $rh_annotations, $output_file);

# parse the command line arguments to determine program options
sub process_commandline {
    
    # Set defaults here
    %Opt = ( 
           );
    GetOptions(\%Opt, qw(
                sites=s ci_files=s outfile=s annotations=s includens help+ version 
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
            my ($chr, $posmo, $pos, $dummy) = split /\s/, $_;
            $all_sites{$chr}->{$pos}++;
        }
        elsif (!(/^#/)) {
            print STDERR "Skipping unorthodox BED file line:\n$_\n";
        }
    }
    close $sites_fh;
    return {%all_sites};
}

sub read_annotations_file {
    my $file = shift;

    my $annots_fh = Open($file);
    my %all_annotations = ();
    my $header_string;
    while (<$annots_fh>) {
        chomp;
        if (/^Chrom\tPos\tRefBase\tVarBase\t(.*)$/) { # header!
            $header_string = $1;
        }
        elsif (!$header_string) {
            die qq!Annotations file must be tab-delimited and have a header with first columns "Chrom", "Pos", "RefBase", "VarBase".\n!;
        }
        elsif (/^(\S+)\t(\d+)\t(\S+)\t(\S+)\t(.*)$/) { # will be all SNPs for now
            my ($chr, $pos, $ref, $var, $annot) = ($1, $2, $3, $4, $5);
            $all_annotations{$chr}->{$pos}->{$ref}->{$var} = $annot;
        }
        elsif (!(/^#/)) {
            print STDERR "Skipping misformatted annotation file line:\n$_\n";
        }
    }
    close $annots_fh;
    return ($header_string, {%all_annotations});
}

sub read_and_write_counts {
    my $file_of_files = shift;
    my $rh_sites = shift;
    my $annot_header = shift;
    my $rh_annots = shift;
    my $outfile = shift;

    my $output_fh = Open($outfile, "w");

    my $fof_fh = Open($file_of_files);
    my @all_files = ();
    my @iters = ();

    my @cols = qw( Index Chrom Pos Ref Acount Tcount Gcount Ccount acount tcount gcount ccount TotalReads VarBases VarCounts VarProps VarCIs );
    my @output_header = qw( Chrom Pos RefBase VarBase );
    my $no_annotation_fields = 0;
    if ($annot_header) {
        my @annot_fields = split /\t/, $annot_header;
        push @output_header, @annot_fields;
        $no_annotation_fields = @annot_fields;
    }
    while (<$fof_fh>) {
        chomp;
        if (/^(\S+)\s(\S+)/) {
            my ($samplename, $countsfile) = split /\s/, $_;
            if (-r $countsfile) {
                push @all_files, {'sample' => $samplename, 'counts' => $countsfile };
                push @output_header, ("$samplename.VarCount", "$samplename.TotalReads", "$samplename.PropVar", "$samplename.PropCI");
                my $fh = Open($countsfile);
                my $dummy_line = <$fh>; # unused, because R header is not useful for now
                push @iters, GTB::File::SortedIter->new( 
                         fh => $fh,
                         cols => \@cols,
                         skip => '/^V1/',
                         key => sub { my $ra = shift->current_array() || return;
                                      my $chr = $ra->[1];
                                      my $pos = $ra->[2];
                                      return sprintf "%s\t%09d", $chr, $pos; },
                         );
            }
            else {
                print STDERR "Skipping file $countsfile--not readable!\n";
            }
        }
    }
    close $fof_fh;
    my $header_string = join "\t", @output_header;
    print $output_fh "$header_string\n";

    foreach my $chr (sort keys %{$rh_sites}) {
        foreach my $pos (sort {$a <=> $b} keys %{$rh_sites->{$chr}}) {
            my $pos_string = sprintf "%s\t%09d", $chr, $pos;
            my $rh_sample_data = {};
            my $ref_allele = 'NA'; # store for later
            for (my $i=0; $i<=$#iters; $i++) {
                my $iter_obj = $iters[$i];
                my $sample = $all_files[$i]->{sample};
                my $countsfile = $all_files[$i]->{counts};
                # attempt to get the file's line for this position:
                my $file_pos = $iter_obj->current_key();
                while (($file_pos) && ($file_pos lt $pos_string)) {
                    #print STDERR "$countsfile:$file_pos (looking for $pos_string)\n";
                    $file_pos = $iter_obj->next_key();
                }

                if (($file_pos) && ($file_pos eq $pos_string)) { # found a record!
                    #print STDERR "Found $file_pos in $countsfile!\n";
                    my @fields = $iter_obj->current_array();
                    my $file_chr = $fields[1];
                    my $file_pos = $fields[2];
                    if ($file_chr ne $chr || $file_pos != $pos) {
                        warn "Our iteration has returned the wrong position: $file_chr:$file_pos for $chr:$pos in sample $sample!\n";
                    }
                    if ($ref_allele ne 'NA' && $ref_allele ne uc($fields[3])) {
                        warn "Multiple reference alleles ($ref_allele, $fields[3]) at position $chr:$pos!\n";
                    }
                    $ref_allele = uc($fields[3]); # 4-11 are allele counts for each base
                    my $total_reads = $fields[12];
                    my @variant_bases = split /,/, $fields[13];
                    my @variant_counts = split /,/, $fields[14];
                    my @variant_props = split /;/, $fields[15];
                    my @variant_cis = split /;/, $fields[16];
                    for (my $i=0; $i<=$#variant_bases; $i++) {
                        my $variant_base = $variant_bases[$i];
                        my $variant_count = $variant_counts[$i];
                        my $variant_prop = $variant_props[$i];
                        my $variant_ci = $variant_cis[$i];
                        if (!defined ($variant_count) || !defined($variant_prop) || !defined($variant_ci)) {
                            warn "Misformatted line at $countsfile:line $fields[0]\n";
                            next;
                        }
                        $rh_sample_data->{$variant_base}->{$sample} = {
                                   'total_reads' => $total_reads,
                                   'variant_count' => $variant_count, 'variant_prop' => $variant_prop,
                                   'variant_ci' => $variant_ci };
                    }
                }
            }

            # now write a record for each alternate allele:

            foreach my $allele (sort keys %{$rh_sample_data}) {
                next if ((!$Opt{includens}) && ($allele eq 'N'));
                my @fields = ($chr, $pos, $ref_allele, $allele);
                if ($rh_annots && $rh_annots->{$chr} && $rh_annots->{$chr}->{$pos} &&
                     $rh_annots->{$chr}->{$pos}->{$ref_allele} &&
                     $rh_annots->{$chr}->{$pos}->{$ref_allele}->{$allele}) {
                     my @annots = split /\t/, $rh_annots->{$chr}->{$pos}->{$ref_allele}->{$allele};
                     push @fields, @annots;
                }
                else {
                     for (my $i=0; $i<$no_annotation_fields; $i++) {
                         push @fields, 'NA';
                     }
                }
                    
                for (my $i=0; $i<=$#iters; $i++) {
                    my $sample = $all_files[$i]->{sample};
                    my $data_allele = $allele;
                    my $rh_data = $rh_sample_data->{$allele}->{$sample};
                    if (!$rh_data && $rh_sample_data->{'N'}) { # do we have an entry for "N"?
                        $rh_data = $rh_sample_data->{'N'}->{$sample}; 
                        $data_allele = 'N' if ($rh_data);
                    }
                    if (!$rh_data) { # do we have any entries for this sample?
                        foreach my $diff_allele (sort keys %{$rh_sample_data}) {
                            if ($rh_sample_data->{$diff_allele}->{$sample}) {
                                $rh_data = $rh_sample_data->{$diff_allele}->{$sample};
                                $data_allele = $diff_allele;
                                last;
                            }
                        }
                    }
                    my $total_reads = ($rh_data) ? $rh_data->{'total_reads'} : 0;
                    my $count = ($rh_data && $data_allele eq $allele) ? $rh_data->{'variant_count'} : 0;
                    my $prop = ($rh_data && $data_allele eq $allele) ? $rh_data->{'variant_prop'} : 0;
                    my $ci = ($rh_data && $data_allele eq $allele) ? $rh_data->{'variant_ci'} : 'NA';
                    push @fields, ($count, $total_reads, $prop, $ci);
                }
                my $record_string = join "\t", @fields;
                print $output_fh "$record_string\n";
            }
        }
    }

    return;
}

=pod

=head1 NAME

combine_allele_counts.pl - read in counts files for each sample, then write out a grand table with all counts.

=head1 SYNOPSIS

	combine_allele_counts.pl --sites <bed file with site locations> --ci_files <file of files with header strings> --outfile <location to write output>

Run combine_allele_counts.pl -man for a detailed description of options and the output files.

=head1 DESCRIPTION

=head1 OPTIONS

=over 5

=item B<--sites>

This option specifies the location of a bedfile (assumed for now to be single base sites of SNPs).

=item B<--ci_files>

This option specifies the location of a file with one record for each sample to be combined.  The first column contains the sample "name" string for the column headers, and the second contains the path to the file of counts for that sample.

=item B<--annotations>

This option specifies the location of a file with annotation columns to be added for particular Chrom/Pos/VarBase combinations.  It must be tab-delimited with the same number of fields in each row, and it must contain a header with first three fields "Chrom", "Pos", and "VarBase".  Subsequent columns will be reported in the columns following the first three in the output file.

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
