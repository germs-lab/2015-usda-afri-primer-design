#!/usr/bin/perl -w

# Lookup table of degenerate IUPAC nucleotide codes.
my %deg2nuc = (
    "R" => ["A", "G"],
    "Y" => ["C", "T"],
    "S" => ["G", "C"],
    "W" => ["A", "T"],
    "K" => ["G", "T"],
    "M" => ["A", "C"],
    "B" => ["C", "G", "T"],
    "D" => ["A", "G", "T"],
    "H" => ["A", "C", "T"],
    "V" => ["A", "C", "G"],
    "N" => ["A", "C", "G", "T"],
    "I" => ["A", "C", "G", "T"]
    );

# Recursive function that replaces degenerate nucleotides with all combinations.
sub generate
{
    if ($_[0] =~ /(.*)([RYSWKBDHVNMI])(.*)/) {
        my $head = $1;
        my $tail = $3;
        my @seqs;
        foreach my $nuc (@{$deg2nuc{$2}}) {
            push @seqs, generate($head.$nuc.$tail);
        }
        return @seqs;
    }
    else {
        return $_[0];
    }
}

# Demo: print all sequences generated from ANCRG.
print join("\n", generate("GGACTACHVGGGTWTCTAAT")), "\n";
