#!/usr/bin/perl
use strict;
use warnings;

use Bio::Seq;
use Bio::AlignIO;
use Bio::Tools::IUPAC;

my $usage = "USAGE:\t./degenerate_code_generator REFERENCE_FILE.fasta Probe_size OUTFILE.txt\n";

my $infile = $ARGV[0] or die $usage . "\n";
my $probesize = $ARGV[1] or die $usage . "\n";
my $stem = $ARGV[2] or die $usage . "\n";

my $alnio = Bio::AlignIO->new(-file => $infile, -format => 'fasta');

open OUTFILE, ">$stem";
select OUTFILE;

my $cs;

while (my $aln = $alnio->next_aln) {
    $cs = $aln->consensus_iupac;
    my $seq = Bio::Seq->new(-seq => $cs, -alphabet => 'dna');
}

$cs =~ s/-//g;
print "<FULL>$cs\n";

my $score;

my $winsize = $probesize;
my $i = 1;
my $cs_len = length($cs);
my %hash;

while ($i < $cs_len) {
        my $window = substr($cs, $i,$winsize);
        $i++;

        foreach (split(//,$window)){
                if(m/A|C|T|G|U/){
                        $score += "5";
                }elsif(m/R|Y|S|W|K|M/){
                        $score += "1";
                }
        }
        $hash{$window} = $score;
        $score = 0;
}

my @keys = sort { $hash{$b} <=> $hash{$a} } keys(%hash);
my @vals = @hash{@keys};

my $c = 0;

#print "Score\tProbe\n";

while ($c < 51){
        print "$vals[$c]\t$keys[$c]\n";
        $c++;
}
