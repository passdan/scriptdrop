#!/usr/bin/perl
use strict;
use warnings;

#####################################
##Dan Pass | 18/01/13              ##
##daniel.antony.pass@googlemail.com##
#####################################

use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::IUPAC;

my $usage = "USAGE:\t./primer_extract INFILE.fasta forward_seq reverse_seq\nALTERNATIVE\t./primer_extract INFILE.fasta leading_offset length_of_extract\nOUTPUT:\tINFILE_trim.fasta and non_matching.fasta";

my $infile = $ARGV[0] or die $usage . "\n";
my $primer_f = $ARGV[1] or die $usage . "\n";
my $primer_r = $ARGV[2] or die $usage . "\n";
my $iupac;

my $in = Bio::SeqIO->new(-file => "$infile", -format => 'fasta');

#####
##Generate forward matching regex 
#####
my $f_flag = 0;
if ($primer_f =~ m/[0-9]+/){
	$f_flag = 1;
}

my $regexp_f;
if ($primer_f =~ m/[A-Z]+/i){
	my $ambiseq_f = Bio::Seq->new(-seq => "$primer_f", -alphabet => 'dna');

	# Create all possible non-degenerate sequences
	$iupac = Bio::Tools::IUPAC->new(-seq => $ambiseq_f);

	$regexp_f = $iupac->regexp($ambiseq_f);
	print "Forward regexp: $regexp_f\n";
}

#####
##Generate reverse matching regex 
#####
my $r_flag = 0;
if ($primer_r =~ m/[0-9]+/){
	$r_flag = 1;
}

my $regexp_r;
if ($primer_r =~ m/[A-Z]+/i){
	my $ambiseq_r = Bio::Seq->new(-seq => "$primer_r", -alphabet => 'dna');

	# Create all possible non-degenerate sequences
	$iupac = Bio::Tools::IUPAC->new(-seq => $ambiseq_r);

	$regexp_r = $iupac->regexp($ambiseq_r);
	print "Reverse regexp: $regexp_r\n";
}

#####
##Trim to matching region
#####
my @missed;
my $all_count = 0;
my $trim_count = 0;
my $fnm_count = 0;
my $rnm_count = 0;

$infile =~ s/\..*//;
open OUTFILE, ">$infile" . "_trim.fasta";
select OUTFILE;
while ( my $seq = $in->next_seq() ){
	my $trim5 = 0;
	my $trim53 = 0;
	$all_count++;

	##Forward trim##
	if ($f_flag == 1){
		my $newseq = $seq->seq;
		$trim5 = substr($newseq, $primer_f);
		#print "5:\t$trim5\n";

	}elsif($f_flag == 0){
		if ( $seq->seq =~ m/^.*$regexp_f(.*)/i ){
			$trim5 = $1;
			#print "5:\t$trim5\n";
		}else{
			push (@missed, $seq->display_id);
			$fnm_count++;
			next;
		}
	}
	##Reverse trim##
	if ($r_flag == 1){
		if(length($trim5) > $primer_r){
			$trim53 = substr($trim5, 1, $primer_r);
		}else{
	                push (@missed, $seq->display_id);
	                $rnm_count++;
        	        next;
		}
	}elsif($r_flag == 0){
		if ($trim5 =~ m/^(.*)$regexp_r.*/i ){
			$trim53 = $1;
			#print "53:\t$trim53\n";
		}else{
			push (@missed, $seq->display_id);
			$rnm_count++;
			next;
		}
	}
	print ">" . $seq->display_id . "\n" . $trim53 . "\n";
	$trim_count++;
}

###
##NON MATCHING SEQUENCES
###
open NONMATCH, ">non_matching.txt";
select NONMATCH;
foreach (@missed){
	print $_ . "\n";
}

select STDOUT;
print "Input sequences: $all_count\n";
print "Trimmed sequences: $trim_count\n";
my $nmc = $fnm_count + $rnm_count;
print "Non-Matching sequences: $nmc  (forward: $fnm_count reverse: $rnm_count)\n";
