#!/usr/bin/env perl

use strict;
use warnings;

use lib "/home/sejj745/miniconda3/envs/atacseq/bin";
use HomerSVGLogo;

my ($motif_file, $output_file) = @ARGV;

die "Usage: motif_to_homer_svg.pl <motif_file> <output_svg>\n"
    unless defined $motif_file && defined $output_file;

open(my $IN, "<", $motif_file)
    or die "Could not open $motif_file: $!\n";

my @matrix;

while (<$IN>) {
    chomp;
    s/\r//g;

    next if /^\s*$/;
    next if /^>/;

    my @line = split /\s+/;

    next unless @line >= 4;

    push @matrix, [
        $line[0] + 0,
        $line[1] + 0,
        $line[2] + 0,
        $line[3] + 0
    ];
}

close($IN);

die "No PWM matrix found in $motif_file\n" unless @matrix;

my $motif = {
    matrix => \@matrix
};

# Same settings used by findKnownMotifs.pl:
# bitFlag = 0 = probability/frequency scaling
# rnaFlag = 0 = DNA (A/C/G/T)
my $svg = HomerSVGLogo::motif2svg($motif, 0, 0);

open(my $OUT, ">", $output_file)
    or die "Could not write $output_file: $!\n";

print $OUT $svg;
close($OUT);

print "Created $output_file\n";
