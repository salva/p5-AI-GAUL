#!/usr/bin/perl

use strict;
use warnings;

use 5.010;

use blib;
use Data::Dumper;
use Algorithm::GAUL;

$| = 1;

sub evaluate {
    shift;
    # print ".";
    my $c = 0;
    $c += unpack("b*", $_) =~ tr/1// for @_;
    ($c & 1 ? undef : $c)
}

sub adapt {
    my $pop = shift;
    # print "o";
    for (@_) {
        # say ">", unpack("b*", $_), "<";
        vec($_, rand($pop->len_chromo), 1) = 1
    }
    '';
}

sub stop {
    my ($pop, $evolution) = @_;
    my $fitness = $pop->fitness_from_rank(0);
    my $len = $pop->len_chromo * $pop->num_chromo;
    print "evolution $evolution, fitness $fitness/$len\n";

    $fitness >= $len;
}

my $p = Algorithm::GAUL->new(num_chromo => 5,
                             population_size => 100,
                             len_chromo => 30,
                             evaluate => \&evaluate,
                             adapt => \&adapt,
                             stop => \&stop,
                             scheme => "lamarck_children" );

say "evolution...";

for (1) {
    $p->evolution(2000);
    my @chromo = $p->chromosomes_from_rank(0);
    my $fit = $p->fitness_from_rank(0);
    print "$_ " for (map unpack("b*", $_), @chromo);
    print "[$fit]\n"
}

say "bye";
