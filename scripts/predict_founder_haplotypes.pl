#!/usr/bin/perl -w

use warnings;
use strict;
use Math::Cephes qw (:gammas);
use Math::Cephes qw (:betas);

use Getopt::Long;

$| = 1;

my $Usage = "predict_founder_haplotypes.pl --recomb <file of recombination data> --founder_haps <file of founder haplotypes> --allele_counts <file of allele counts> --nogens <number of generations> [--chroms <comma-delimited list of chroms>] [--rils <comma-delimited list of RILs>]\n";

my ($recomb_file, $chromstring_opt, $rils_opt, $comeron_opt, $no_gens, $verbose_opt);
my $trans_scale = 1.0; # scaling down by 1/10 was used in King et al., 2012
my $founder_gt_file = 'freeze2_informative_polymorphisms_corrected_new.dm6.tdf.txt';
my $allele_file = 'STHInbred_annotated_counts.dm6.txt.gz';

GetOptions("recomb=s" => \$recomb_file, "founder_haps=s"=> \$founder_gt_file,
           "allele_counts=s" => \$allele_file, "chroms=s" => \$chromstring_opt,
           "rils=s"=> \$rils_opt, "comeron" => \$comeron_opt, "nogens=i" => \$no_gens,
           "trans_scale=f" => \$trans_scale, "verbose" => \$verbose_opt);

($#ARGV<0 && $recomb_file && $founder_gt_file && $allele_file && $no_gens)
    or die "$Usage";

my @desired_chromarms = ($chromstring_opt) ? split /,/, $chromstring_opt : (); # if empty, will process all
my @desired_rils = ($rils_opt) ? split /,/, $rils_opt : (); # if empty, will process all

# read in founder haplotype data (hash of chrom, pos, alleles, founder):
my ($ra_founders, $rh_founder_probs) = read_founder_haplotypes($founder_gt_file, \@desired_chromarms);
my $n_founders = scalar(@{$ra_founders});
my @chromosomes = keys %{$rh_founder_probs}; # will only include desired chroms if specified

# observed allele counts in recombinant inbred lines (hash of chrom, pos, RIL sample):
my ($ra_rils, $rh_ril_counts) = read_allele_counts($allele_file, $rh_founder_probs, \@desired_chromarms, \@desired_rils);

# codes for all states (homozygous and heterozygous):
my %statecodes = create_statecodes($ra_founders);
my @statecodes = sort keys(%statecodes);

# initial state probabilities:
my %initiation = create_initialprobs($ra_founders, \%statecodes);

# recombination data (to be used to calculate transition probabilities):
my %recombdat = read_recombination_data($recomb_file);

# Genotyping error rate
my $err = 0.005;
my $alpha = 0.5;
my $beta = 0.5;

my %binom_coef_cache = (); # cache binomial coefficients to save time
my %ps_cache = (); # cache probability sums to save time

foreach my $arm (@chromosomes) { 

    ############################################
    # Make index and position hash for arm
    # Positions are stored as 10 x MB

    my @chromarm_positions = sort {$a <=> $b} keys %{$rh_founder_probs->{$arm}};

    ############################################
    #Loop through RILS

    foreach my $ril_sample (@{$ra_rils}) {

        open(PROBFILE, ">$ril_sample.state_probabilities.$arm.csv");
    
        #print header
        print PROBFILE "rilID", "," ,"index", "," , "chr", "," , "pos", "," ;
        foreach my $state (@statecodes){
            print PROBFILE "$state", ",";
        }
        foreach my $founder (@{$ra_founders})
        {
            print PROBFILE "Total$founder",",";
        }
        print PROBFILE "\n";

        ############################################
        # HMM

        ############################################
        #forward equations (alphas)
        my %alphas = ();

        # For each position ($t) get alpha: count forward
        for (my $t=0; $t<=$#chromarm_positions; ++$t){
            print "$ril_sample, $t\n" if ($verbose_opt);
            my $thispos = $chromarm_positions[$t];
            my ($thisposalleles) = keys %{$rh_founder_probs->{$arm}->{$thispos}};

            # for first position
            if ($t == 0){
                my %emission = emission($t, $rh_founder_probs->{$arm}->{$thispos}->{$thisposalleles},
                                            $rh_ril_counts->{$arm}->{$thispos}->{$ril_sample},
                                            \%binom_coef_cache, \%ps_cache );
                foreach my $state (@statecodes){
                    $alphas{$t}{$state} = elnproduct($initiation{$state}, $emission{$state});
                }
            }else{
                # Get transition probabilities for position and previous position
                my $lastpos = $chromarm_positions[$t-1];
                my %transitions = trans($lastpos,$thispos,$arm);

                # Get emission probabilities for position
                my %emission = emission($thispos, $rh_founder_probs->{$arm}->{$thispos}->{$thisposalleles},
                                            $rh_ril_counts->{$arm}->{$thispos}->{$ril_sample},
                                            \%binom_coef_cache, \%ps_cache );

                foreach my $state1 (@statecodes){
                    my $sum = 'logzero';

                    # Sum over founders
                    foreach my $state2 (@statecodes){

                        # The alpha is already eln from above etc.
                        my $a1_trans = elnproduct($alphas{$t-1}{$state2}, $transitions{$state1}{$state2});
                        $sum = elnsum($sum, $a1_trans);
                    }
                    $alphas{$t}{$state1}=elnproduct($sum, $emission{$state1});
                }
            }
        }

        ############################################
        # Backward equations (betas)
        my %betas = ();

        # For each position ($t) get beta: count backwards
        for (my $t=$#chromarm_positions; $t>=0; --$t){
            
            print "$ril_sample, $t\n" if ($verbose_opt);
            #for last position
            if ($t == $#chromarm_positions){
                foreach my $state (@statecodes){
                    $betas{$t}{$state} = eln(1);
                }
            }else{
                my $lastpos = $chromarm_positions[$t+1];
                my $thispos = $chromarm_positions[$t];
                # Get transition for position and position ahead
                my %transitions = trans($thispos,$lastpos,$arm);
                my ($lastposalleles) = keys %{$rh_founder_probs->{$arm}->{$lastpos}};
                # Get emission for position ahead
                my %emission = emission($lastpos, $rh_founder_probs->{$arm}->{$lastpos}->{$lastposalleles},
                                              $rh_ril_counts->{$arm}->{$lastpos}->{$ril_sample},
                                              \%binom_coef_cache, \%ps_cache);
                foreach my $state1 (@statecodes){
                    my $sum = 'logzero';
                    # Sum over founders
                    foreach my $state2 (@statecodes){
                        my $b_tm1 = elnproduct(elnproduct($betas{$t+1}{$state2}, $transitions{$state1}{$state2}), $emission{$state2});
                        $sum = elnsum($sum, $b_tm1);
                    }
                    $betas{$t}{$state1}=$sum;
                }
            }
        }

        # Hash to store probability of state for each position
        # Get Y for each position
        my %prob_state = ();
        for (my $t=0; $t<=$#chromarm_positions; ++$t){
            print "$ril_sample, $t\n" if ($verbose_opt);
            foreach my $state1 (@statecodes){
                my $sum = 'logzero';
                foreach my $state2 (@statecodes){
                    my $value = elnproduct($alphas{$t}{$state2},$betas{$t}{$state2});
                    $sum = elnsum($sum, $value);
                }
                # Get Y out of log space
                my $thisposition = $chromarm_positions[$t];
                $prob_state{$thisposition}->{$state1} = eexp(elnproduct(elnproduct($alphas{$t}{$state1},$betas{$t}{$state1}), -$sum));
            }
        }
        for (my $t=0; $t<=$#chromarm_positions; ++$t) {
            print "$ril_sample, $t\n" if ($verbose_opt);
            my $thispos = $chromarm_positions[$t];
            print PROBFILE "$ril_sample", "," ,"$t", "," , "$arm", "," , "$thispos", "," ;
            foreach my $state (@statecodes){
                print PROBFILE $prob_state{$thispos}->{$state}, ",";
            }
            foreach my $founder (@{$ra_founders}){
                my $addprob = additive(\%prob_state, $thispos,$founder);
                print PROBFILE "$addprob",",";
            }
            print PROBFILE "\n";
        }
        close(PROBFILE);
    } #close foreach for rils
} #close foreach for arms


############################################
############################################
############################################
# Subroutines

sub read_founder_haplotypes {
    my $founder_file = shift;
    my $ra_desiredchroms = shift;

    my $founderstring = ($founder_file =~ /\.gz$/) ? "gunzip -c $founder_file | " : $founder_file;

    open FOUNDERS, $founderstring
        or die "Couldn\'t open $founderstring for reading $!\n";

    my $header_line = <FOUNDERS>;
    chomp $header_line;
    ($header_line =~ s/^Chrom\tPos\tType\tRef\tAlt\t//)
        or die "Improperly formatted header in $founder_file! First fields must be Chrom, Pos, Type, Ref, and Alt, and file must be tab-delimited.\n";

    # header line no longer has chrom, pos, type, ref, and alt
    my @header_fields = split /\t/, $header_line;

    my $ra_founders = \@header_fields;
    
    my %founderprob = ();
    my $sites_read = 0;
    while (<FOUNDERS>) {
        chomp;
        my @fields = split /\t/, $_;
        my $chrom = shift @fields;
        next if ((@{$ra_desiredchroms}) && !(grep {$_ eq $chrom} @{$ra_desiredchroms}));
        my $pos = shift @fields;
        my $type = shift @fields;
        my $refallele = shift @fields;
        my $altallele = shift @fields;

        if ($#fields != $#{$ra_founders}) { # should be 10 lines
            die "Line with wrong number of fields in $founder_file!\n";
        }
        else {
            # switching to 50/50 probability for alleles when "-" is reported--more accurate
            #my $avg_allele_freq = (grep {$_ eq '-'} @fields) ? calc_avg_prob(@fields) : undef;
            for (my $i=0; $i<=$#fields; $i++) {
                my $founder = $ra_founders->[$i];
                my $foundergt = $fields[$i];
                my $founder_af = ($foundergt eq 2) ? 1.0 : ($foundergt eq 0) ? 0.0 : 0.5;
                if (!(defined($founder_af))) {
                    die "Unable to calculate founder allele frequency for genotype $foundergt: file $founder_file should have 0, 2, or -.\n";
                }
                $founderprob{$chrom}->{$pos}->{"$refallele:$altallele"}->{$founder} = $founder_af;
            }
        }
        $sites_read++;
        if ($sites_read == 100000*int($sites_read/100000)) {
            print "Read $sites_read founder sites\n" if ($verbose_opt);
        }
    }
    close FOUNDERS;

    return ([@header_fields], {%founderprob});
}

sub calc_avg_prob {
    my @founder_gts = @_;

    my ($alt_total, $total_gts) = (0, 0);
    foreach my $foundergt (@founder_gts) {
        if ($foundergt eq '2') { # alt allele
            $alt_total++;
            $total_gts++;
        }
        elsif ($foundergt eq '0') { # ref allele
            $total_gts++;
        }
    }

    my $avg_prob = ($total_gts) ? $alt_total/$total_gts : undef;
    return $avg_prob;
}

sub read_allele_counts {
    my $allele_file = shift;
    my $rh_desiredsites = shift;
    my $ra_desiredchroms = shift;
    my $ra_desiredrils = shift;

    my $allelefilestring = ($allele_file =~ /\.gz$/) ? "gunzip -c $allele_file | " : $allele_file;

    open AC, $allelefilestring
        or die "Couldn\'t open $allelefilestring for reading $!\n";

    my $header_line = <AC>;
    chomp $header_line;
    ($header_line =~ s/^Chrom\tPos\tRefBase\tVarBase\t//)
        or die "Improperly formatted header in $allele_file! First fields must be Chrom, Pos, RefBase, and VarBase, and file must be tab-delimited.\n";

    my @header_fields = split /\t/, $header_line;
    my %field_indices = map {$header_fields[$_] => $_} (0..$#header_fields);
    my @samples = ();
    foreach my $varfield (grep /\.VarCount/, @header_fields) {
        my $samplebase = $varfield;
        $samplebase =~ s/\.VarCount//;
        next if ((@{$ra_desiredrils}) && (!grep {$_ eq $samplebase} @{$ra_desiredrils}));
        my $totalfield = "$samplebase.TotalReads";
        if (!grep {$_ eq $totalfield} @header_fields) {
            die "File $allele_file has a column $varfield but no field for total count $totalfield!\n";
        }
        push @samples, $samplebase;
    }

    my %ril_allele_counts = ();
    my $sites_read = 0;
    while (<AC>) {
        chomp;
        my @fields = split /\t/, $_;
        if ($#fields != $#header_fields + 4) {
            die "Wrong number of fields ($#fields, $#header_fields):\nFile $allele_file\n$_\n";
        }
        my $chromarm = shift @fields;
        my $pos = shift @fields;
        my $refallele = shift @fields;
        my $varallele = shift @fields;
        if (!$rh_desiredsites->{$chromarm}->{$pos}) {
            next;
        }

        my %fieldhash = map {$header_fields[$_] => $fields[$_]} (0..$#fields);
        foreach my $sample (@samples) {
            next if ((@{$ra_desiredrils}) && (!grep {$_ eq $sample} @{$ra_desiredrils}));
            my $total = $fieldhash{"$sample.TotalReads"};
            if ($total ne 'NA') {
                my $varcount = $fieldhash{"$sample.VarCount"};
                $ril_allele_counts{$chromarm}->{$pos}->{$sample} = [$varcount, $total - $varcount, "$refallele:$varallele"];
                print "$chromarm:$pos:$sample\t$varcount\t$total\t$refallele:$varallele\n" if ($verbose_opt);
            }
        }
        $sites_read++;
        if ($sites_read == 100000*int($sites_read/100000)) {
            print "Read $sites_read RIL sites\n" if ($verbose_opt);
        }
    }
    close AC;

    return ([@samples], {%ril_allele_counts});
}

sub read_recombination_data {
    my $recomb_file = shift;

    # RECOMBINATION FREQUENCY INPUTS ARE IN CM/MB/generation
    open (RF, $recomb_file) || die ("Can't open $recomb_file: $!");
    
    # Store arm, position (in units of 10kbp), and RF (in CM/bp/50 generations)
    my %recombdat = ();
    
    ## Populate %recombdat using reads from the table
    ## Using new file generated by script at http://petrov.stanford.edu/cgi-bin/recombination-rates_updateR5.pl
    foreach my $line (<RF>) {
        chomp $line;
        my @fields = split /\s/, $line;
        my $chrom = $fields[0];
        my $pos1 = $fields[1];
        my $pos2 = $fields[2];
        my $mareymapstart = $fields[$#fields-5]; # in CM/Mb/gens
        my $mareymapcenter = $fields[$#fields-4];
        my $mareymapend = $fields[$#fields-3];
        my $comeronmapstart = $fields[$#fields-2];
        my $comeronmapcenter = $fields[$#fields-1];
        my $comeronmapend = $fields[$#fields];
        my $avgpos = int(0.5*($pos1+$pos2));
        my $nochromgens = ($chrom =~ /X/) ? 0.6667*$no_gens : 0.5*$no_gens;
        $recombdat{$chrom}{$avgpos} = ($comeron_opt) ? # now in CM/bp/$no_gens gens
                         1E-6 * $comeronmapcenter * $nochromgens :
                         1E-6 * $mareymapcenter * $nochromgens;
        print "Set recomb rate for $chrom $avgpos to $recombdat{$chrom}{$avgpos} CM/bp/$no_gens generations\n" if ($verbose_opt);
    }
    close RF;

    return %recombdat;
}

sub create_statecodes {
    my $ra_founders = shift;

    # Make states hash with homozygous and heterozygous states and founders array (0,1)
    my %statecodes=();
    for (my $i=0; $i<@{$ra_founders}; ++$i){
        # homozygous states:
        $statecodes{$ra_founders->[$i].$ra_founders->[$i]}=[$ra_founders->[$i], $ra_founders->[$i]];
        for (my $j=0; $j < @{$ra_founders}; ++$j){
            if ($i < $j){
                # heterozygous states:
                $statecodes{$ra_founders->[$i].$ra_founders->[$j]}=[$ra_founders->[$i], $ra_founders->[$j]];
            }
        }
    }
    return %statecodes;    
}

sub create_initialprobs {
    my $ra_founders = shift;
    my $rh_statecodes = shift;

    my $n_founders = @{$ra_founders};
    my %statecodes = %{$rh_statecodes};

    # Set up initiation probability hash
    my %initiation = ();
    
    # Number of heterozygous states
    # This will be N*(N-1)/2, where N is the number of founder states
    my $n_hets=0;
    for (my $i=$n_founders-1; $i > 0; --$i) {
        $n_hets = $n_hets+$i;
    }
    
    # Set parameters needed
    # Frequency of homozygous states
    my $freq_homo = 0.95;
    
    # Frequncy of heterozygous states
    my $freq_hets = 1 - $freq_homo;
    
    # Fill initiation hash with initiation probabilities
    my $total_prob = 0;
    foreach my $state (@statecodes){
        if ($statecodes{$state}->[0] eq $statecodes{$state}->[1]) { # homs have total freq_homo, evenly dist
           $initiation{$state} = eln((1/$n_founders) * $freq_homo);
           $total_prob += (1/$n_founders) * $freq_homo;
        } else {
            $initiation{$state} = eln((1/$n_hets) * $freq_hets); # hets have total freq_hets, evenly dist
           $total_prob += (1/$n_hets) * $freq_hets;
        }
    }
    print "Total initial probabilities: $total_prob\n" if ($verbose_opt);

    return %initiation;
}

sub read_inbredline_names {
    my $rilfile = shift;

    my @rilnames = ();
    open RIL, $rilfile
        or die "Couldn\'t open $rilfile for reading: $!\n";
    while (<RIL>) {
        next if (/^#/); # skip comments

        if (/^\s*(\S+)/) {
            push @rilnames, $1;
        }
    }
    close RIL;

    return @rilnames;
}

############################################
# Get additive model

sub additive {
    my ($rh_probstate, $pos, $founder ) = @_;
    my $addsum=0;
    foreach my $state (@statecodes){
        my $found1 = $statecodes{$state}->[0];
        my $found2 = $statecodes{$state}->[1];
        if ($found1 eq $founder && $found2 eq $founder){
            $addsum=$addsum+$rh_probstate->{$pos}->{$state};
        }
        elsif (($found1 eq $founder) || ($found2 eq $founder)){
            $addsum=$addsum+(0.5*$rh_probstate->{$pos}->{$state});
        }
    }
    return $addsum;
}

############################################
# Transition probabilities

sub trans {
    my ($pos1, $pos2, $arm) = @_;

    # Set parameters for transition probabilities
    my $prob_homo = 0.95;
    my $prob_transhomo = 0.85;

    my %transitions = ();
    my $trans_prob = recomb_prob($pos1,$pos2,$arm); # probability of recombination between pos1 and pos2

    # fill transitions hash with trans probs from all states at pos1 to all states at pos2
    foreach my $state (@statecodes){
        foreach my $state1 (@statecodes){
            if ($state eq $state1){
                # Probability of staying in state -- could not recombine, or recombine and randomly select same state
                if ($statecodes{$state}->[0] eq $statecodes{$state}->[1]){ # homozygous
                    $transitions{$state}{$state1} = eln((1 - $trans_prob) + (($prob_homo*$trans_prob)/$n_founders));
                    print "Probability of staying in state $state: $transitions{$state}{$state1}\n" if ($verbose_opt);
                }else{ # heterozygous -- event could be on either chrom, but must bring back to same founder allele
                    $transitions{$state}{$state1} = eln((1 - $trans_prob) + 2 * (1-
                        $prob_transhomo) * $trans_prob * (1/(2*($n_founders-1))));
                }
            }else{
                # If original is homozygous:
                if ($statecodes{$state}->[0] eq $statecodes{$state}->[1]){
                    # And if new is homozygous
                    if ($statecodes{$state1}->[0] eq $statecodes{$state1}->[1]){
                        # Probability of moving to another homozygous state
                        $transitions{$state}{$state1} = eln($trans_prob*$prob_homo*(1/($n_founders)));
                    }else{
                        # Original homozygous, new is heterozygous with one copy same as original
                        if (($statecodes{$state}->[0] eq $statecodes{$state1}->[0]) or
                             ($statecodes{$state}->[0] eq $statecodes{$state1}->[1])){
                            # Probability of comment above
                            $transitions{$state}{$state1} = eln($trans_prob*(1- $prob_homo)*(1/($n_founders-1)));
                        }else{
                            # All other moves from homozygous state involve 2 events, prob near 0
                            $transitions{$state}{$state1} = eln(0);
                        }
                    }
                }else{
                    # If original heterozygous, new is one of homo states (only btw two represented in het)
                    if (($state1 eq $statecodes{$state}->[0].$statecodes{$state}->[0]) or
                         ($state1 eq $statecodes{$state}->[1].$statecodes{$state}->[1])){
                        $transitions{$state}{$state1} = eln($trans_prob * $prob_transhomo * 0.5);
                    }else{
                        # Origninal heterozygous, new is het with one founder same as original het 
                        if (($statecodes{$state1}->[0] eq $statecodes{$state}->[0]) or
                              ($statecodes{$state1}->[0] eq $statecodes{$state}->[1]) or
                              ($statecodes{$state1}->[1] eq $statecodes{$state}->[0]) or
                              ($statecodes{$state1}->[1] eq $statecodes{$state}->[1])){
                            $transitions{$state}{$state1} = eln($trans_prob * (1-$prob_transhomo) * (1/(2*($n_founders-1))));
                        }else{
                            # All other moves from het involve two events, prob near zero
                            $transitions{$state}{$state1} = eln(0);
                        }
                    }
                }
            }
        }
    }
    return %transitions;
}

sub recomb_prob {
    # positions are chromosome coordinates
    my ($pos1, $pos2, $arm) = @_;

    #get distance betweeen markers in bp
    my $dist = $pos2 - $pos1;

    #get midpoint between markers for lookup (in 10kbp units)
    my $midpt = ($pos1 + $pos2)/2;

    #find positions for table lookup--uses RF TABLE
    my $LB = int($midpt);
    my $UB = $LB + 1;
    #below is reg expression to test if midpoint is an integer
    ##if ($midpt =~ /^\d+$/) {
    if (defined($recombdat{$arm}{$midpt})) {
        #if is at exact integer, just get from table
        my $rf = $recombdat{$arm}{$midpt} * 0.01; # one CM = 0.01 recomb prob

        #and multiply by distance
        my $trans_prob = (1-exp -($rf*$dist*$trans_scale));
        print "Returning defined $trans_prob at $arm:$pos1-$pos2\n" if ($verbose_opt);
        return $trans_prob;
    }else{
        #use nearest value, within reason (max 50000 bases away)
        my $steps = 0;
        while(!defined($recombdat{$arm}{$UB})) { # use nearest
            $UB++;
            $steps++;
            if ($steps >= 50000) {
                last;
            }
        }
        $steps = 0;
        while (!(defined($recombdat{$arm}{$LB}))) {
            $LB--;
            $steps++;
            if ($steps >= 50000) {
                last;
            }
        }
        if (!(defined($recombdat{$arm}{$LB})) || !(defined($recombdat{$arm}{$UB}))) {
            my $trans_prob = 1 - exp(-(0.05*$dist/1000000*$trans_scale));
            print STDERR "Returning default $trans_prob at $arm:$pos1-$pos2\n";
            return $trans_prob;
        }
        my $slope = ($recombdat{$arm}{$UB} - $recombdat{$arm}{$LB})/($UB-$LB);
        my $rf = $recombdat{$arm}{$LB} + $slope * ($midpt - $LB);

        #and multiply by distance
        my $trans_prob = (1-exp -($rf*$dist*$trans_scale));
        print "Returning interpolated $trans_prob at $arm:$pos1-$pos2\n" if ($verbose_opt);
        return $trans_prob;
    }
}

############################################
# Emission probabilities

sub emission {
    #pass index
    my $pos = shift;
    my $rh_position_probs = shift;
    my $rh_position_sample_acs = shift;
    my $rh_cache_binom = shift;
    my $rh_cache_ps = shift;

    my %emission;

    if ($rh_position_sample_acs){
        #get observed counts (k,N)
        my $O_n1 = $rh_position_sample_acs->[0]; #k
        my $O_n2 = $rh_position_sample_acs->[1];
        my $O_N = $O_n1 + $O_n2; #N
        #put in if for 0,0
        if ($O_N>0){
            foreach my $state (@statecodes){
                #get founders for each state (homo or het)
                my $found1 = $statecodes{$state}->[0];
                my $found2 = $statecodes{$state}->[1];
                #get probability of obs "A" for each founder

                my $p1 = $rh_position_probs->{$found1};
                my $p2 = $rh_position_probs->{$found2};

                if (!(defined($p1)) || !(defined($p2))) {
                    print "No position prob for $found1/$found2 ($pos, $O_n1, $O_n2)!\n";
                }

                if ($p1 >0.995) {
                    $p1=0.995;
                }
                if ($p1 <0.005) {
                    $p1=0.005;
                }
                if ($p2 >0.995) {
                    $p2=0.995;
                }
                if ($p2 <0.005) {
                    $p2=0.005;
                } 
                #get probability of obs k "A"'s in N observations given state.
                #sum over possible genotypes of ril
                my $prob_sum;
                if (defined($rh_cache_ps->{$p1}->{$p2}->{$O_N}->{$O_n1})) {
                    $prob_sum = $rh_cache_ps->{$p1}->{$p2}->{$O_N}->{$O_n1};
                }
                else {
                    $prob_sum = eln(($p1 * $p2 * (1-$err)**$O_n1 * $err**($O_N-$O_n1))+
                       ((1-$p1) * $p2 * (Math::Cephes::beta($alpha+$O_n1, $O_N+$beta- $O_n1)/Math::Cephes::beta($alpha, $beta)))+
                       ($p1 * (1-$p2) * (Math::Cephes::beta($alpha+$O_n1, $O_N+$beta- $O_n1)/Math::Cephes::beta($alpha, $beta)))+
                       ((1-$p1) * (1-$p2) * $err**$O_n1 * (1-$err)**($O_N-$O_n1)));
                    $rh_cache_ps->{$p1}->{$p2}->{$O_N}->{$O_n1} = $prob_sum;
                }

                #get binoial coefficient
                my $binom_coef;
                if (defined($rh_cache_binom->{$O_N}->{$O_n1})) {
                    $binom_coef = $rh_cache_binom->{$O_N}->{O_n1};
                }
                else {
                    $binom_coef = Math::Cephes::lgam($O_N + 1) - Math::Cephes::lgam($O_n1+1) - Math::Cephes::lgam($O_N-$O_n1+1);
                    $rh_cache_binom->{$O_N}->{O_n1} = $binom_coef;
                }
                #calculate probability
                if ($prob_sum eq 'logzero'){
                    $emission{$state} = 'logzero';
                }
                else{
                    $emission{$state} = $binom_coef + $prob_sum;
                }
            }#foreach close
        }else{
            foreach my $state (@statecodes){
                $emission{$state} = eln(1);
            }
        }
    }else{
        foreach my $state (@statecodes){
            $emission{$state} = eln(1);
        }
    }#ifelse close

    return %emission;

} #sub close

############################################
#FUNCTIONS FOR LOG SPACE#

# Extended exp
sub eexp{
    #$X is a ln prob
    my ($X) = @_;
    if ($X eq 'logzero'){
        return 0;
    }else{
        return exp $X;
    }
}

# Extended ln
sub eln{
    #$X is raw prob
    my ($X) = @_;
    if ($X eq 0){
        return 'logzero';
    }else{
        return log $X;
    }
}

# Extended ln sum
sub elnsum{
    #$X, $Y are ln probs
    my ($X, $Y) = @_;
    if (($X eq 'logzero') or ($Y eq 'logzero')){
        if ($X eq 'logzero'){
            return $Y;
        }else{
            return $X;
        }
    }else{
        if ($X > $Y){
            return $X + eln(1 + exp($Y-$X));
        }else{
            return $Y + eln(1 + exp($X-$Y));
        }
    }
}

# Extended product
sub elnproduct{
    #$X, $Y are ln probs
    my ($X, $Y) = @_;
    if (($X eq 'logzero') or ($Y eq 'logzero')){
        return 'logzero';
    }else{
        return $X + $Y;
    }
}
