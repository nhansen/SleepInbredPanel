#!/usr/bin/perl -w

use strict;

my $Usage="report_states.pl <HMM output file>\n";

$#ARGV==0
    or die $Usage;

my $hmm_file = $ARGV[0];

open HMM, $hmm_file
    or die "Couldn\'t open $hmm_file: $!\n";
my $firstline = <HMM>;
chop $firstline;
my @fieldnames = split /,/, $firstline;
while ($fieldnames[$#fieldnames] =~ /^A/ || $fieldnames[$#fieldnames] =~ /^Total/) {
    pop @fieldnames;
}

my @statenames = map {$fieldnames[$_]} (4..$#fieldnames);
my $no_states = @statenames;
#print "$no_states states reported\n";
print "Sample\tChrom\tPos\tMLState\tMLProb\tHomHet\tSwitch\tQuality\n";
my %index_fieldnames = map {$_ => $fieldnames[$_]} (0..$#fieldnames);

my $last_state = '';
my $qual_threshold = 0.99;
while (<HMM>) {
    chomp;
    my @fields = split /,/, $_;
    my $max_prob = 0;
    my $max_state = 'NA';
    for (my $i=4; $i<=$#fieldnames; $i++) {
        if ($fields[$i] > $max_prob) {
            $max_prob = $fields[$i];
            $max_state = $index_fieldnames{$i};
        }
    }
    my $qual = ($max_prob >= $qual_threshold) ? 'H' : 'L';
    my $switch = ($max_state ne $last_state) ? 'D' : 'C';
    my $type = ($max_state =~ /^(.+?)(?=\1)/g) ? 'Hom' : 'Het';
    print "$fields[0]\t$fields[2]\t$fields[3]\t$max_state\t$max_prob\t$type\t$switch\t$qual\n";
    $last_state = $max_state;
}
close HMM;
