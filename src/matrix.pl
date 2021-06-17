#!/usr/bin/perl -w

use strict;
use warnings;

# get command-line arguments, or die with a usage statement
my $usage = "usage: matrix.pl UTR.tab\n";
my $file  = shift or die $usage;

open (READ, "$file") || die "cannot open $file: $!";
while (<READ>){
	if (/^(\S+)\t(\S+)$/){
		my ($seq_id,$seq_nuc) = ($1,$2);
		print $seq_id,"\n\n";

		# Forward matrix (F)
		my @seq = split (//,$seq_nuc);
		print "F\t",join ("\t",@seq),"\n";
		foreach my $x (@seq) {
			print $x;
			foreach my $y (@seq) {
				if ($x eq $y) {
					print "\t5";
				}
				elsif ($x ne $y) {
					print "\t-4";
				}								
			}
			print "\n";
		}
		print "\n";
		undef @seq;
	
		# Reverse complement matrix (R)
		$seq_nuc = reverse $seq_nuc;
		$seq_nuc =~ tr/ATGCatgc/TACGtacg/;
		@seq = split (//,$seq_nuc);
		print "R\t",join ("\t",@seq),"\n";
		foreach my $x (@seq) {
			print $x;
			foreach my $y (@seq) {
				if ($x eq $y) {
					print "\t5";
				}
				elsif ($x ne $y) {
					print "\t-4";
				}								
			}
			print "\n";
		}
	}
	print "\n\n";
}
close (READ);
