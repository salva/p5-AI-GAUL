#!/usr/bin/perl

use 5.010;

use strict;
use warnings;

use Algorithm::GAUL;
use Math::Vector::Real;
use Math::Vector::Real::Random;
use Math::Vector::Real::kdTree;
use Bit::Util qw(bu_count bu_first);
use GD;

my $width = 1024;


$| = 1;

my $dist = shift @ARGV || 1;
my $pop_size = shift @ARGV || 100;
my @p;

while (<>) {
    if (my ($v) = /\{\s*([^\s,]+(?:\s*,\s*[^\s,]+)*)\s*\}/) {
        push @p, V(split /\s*,\s*/, $v);
    }
}

#my $t = Math::Vector::Real::kdTree->new(@p);
#@p = @p[$t->ordered_by_proximity];

my @cover;
my $t = Math::Vector::Real::kdTree->new(@p);
for my $p (@p) {
    push @cover, [$t->find_in_ball($p, $dist)];
}

for my $ix (0..$#p) {
    my $n = 0;
    my $p = $p[$ix];
    for my $q (@p) {
        $n++ if $p->dist2($q) < $dist * $dist
    }
    my $c = @{$cover[$ix]};
    if ($c != $n) {
        warn "find in ball failed for ix: $ix, c: $c, n: $n";
    }
}

say "cover generated";

my $missing = '';
vec($missing, $_, 1) = 1 for (0..$#p);

sub weight {
    for my $s (@{$_[1]}) {
        my $m = $missing;
        my $ix = -1;
        while (defined($ix = bu_first($s, $ix + 1))) {
            vec($m, $_, 1) = 0 for @{$cover[$ix]};
        }
        return @p - (bu_count($m) + bu_count($s));
    }
}

sub stop { 0 }

sub repair {
    # say "repairing chromosomes";
    for my $s (@{$_[1]}) {
        my $m = $missing;
        my $ix = -1;
        while (defined($ix = bu_first($s, $ix + 1))) {
            vec($m, $_, 1) = 0 for @{$cover[$ix]};
        }
        $ix = -1;
        while (defined($ix = bu_first($m, $ix + 1))) {
            my $r = $cover[$ix][rand scalar @{$cover[$ix]}];
            vec($s, $r, 1) = 1;
            vec($m, $_, 1) = 0 for @{$cover[$r]};
        }
    }
}

sub mutate {
    # say "mutate";
    for my $s (@{$_[1]}) {
        for (0..9) {
            my $ix = int rand scalar @p;
            if (vec($s, $ix, 1)) {
                vec($s, $ix, 1) = 0;
                my $c = $cover[$ix];
                vec($s, $c->[rand @$c],  1) = 1;
            }
        }
    }
}

my ($b0, $b1) = Math::Vector::Real->box(@p);
my $db = $b1 - $b0;

sub crossover {
    # say "crossing...";
    my $z = $b0 + $db->random_in_box;
    my $dir = Math::Vector::Real->random_versor(scalar @$b0);

    my $s = \$_[3][0];
    my $d = \$_[4][0];

    my $p = $_[1][0];
    my $m = $_[2][0];

    for my $ix (0..$#p) {
        if ($dir * ($p[$ix] - $z) >= 0) {
            vec($$s, $ix, 1) = vec($p, $ix, 1);
            vec($$d, $ix, 1) = vec($m, $ix, 1);
        }
        else {
            vec($$s, $ix, 1) = vec($m, $ix, 1);
            vec($$d, $ix, 1) = vec($p, $ix, 1);
        }
    }
}

sub seed {
    for my $s (@{$_[1]}) {
        $s = "\x00" x length $s;
        my $m = $missing;
        my $ix = -1;
        while (defined ($ix = bu_first($m, $ix + 1))) {
            my $c = $cover[$ix][rand scalar @{$cover[$ix]}];
            vec($s, $c, 1) = 1;
            vec($s, $_, 1) = 0 for @{$cover[$c]};
        }
    }
}

my $gaul = Algorithm::GAUL->new(len_chromo => scalar(@p),
                                population_size => $pop_size,
                                select_two => 'random',
                                selec_one => 'random',
                                mutation_ratio => 0.05,
                                seed => \&seed,
                                evaluate => \&weight,
                                adapt => \&repair,
                                stop => \&stop,
                                mutate => \&mutate,
                                elitism => 'parents_survive',
                                scheme => 'lamarck_children',
                                crossover => \&crossover);


sub scl {
    my $p = shift;
    @{0.125 * $width * ($p + [4,4])};
}

my $scl_dist = $dist * 0.125 * 2 * $width;

say "scl_dist: $scl_dist";

my $iteration = 0;
while (1) {
    $iteration++;
    $gaul->evolution(1);
    my $fit = $gaul->fitness_by_rank(0);
    say "iteration: $iteration, fit: $fit";

    my $s = $gaul->entity_by_rank(0)->[0];

    my $len = length($s) * 8;
    my $p = scalar @p;

    my $im = GD::Image->new($width, $width, 0);
    my $white = $im->colorAllocate(255, 255, 255);
    $im->filledRectangle(0, 0, $width, $width, $white);

    my $black = $im->colorAllocate(0, 0, 0);
    my $gray = $im->colorAllocate(200, 200, 200);
    my $green = $im->colorAllocate(0, 150, 0);
    my $red = $im->colorAllocate(255, 0, 0);
    my $yellow = $im->colorAllocate(255, 255, 0);

    my $m = $missing;
    my $ix = -1;

    while (defined($ix = bu_first($s, $ix + 1))) {
        vec($m, $_, 1) = 0 for @{$cover[$ix]};
        $im->filledEllipse(scl($p[$ix]), $scl_dist, $scl_dist, $yellow);
    }

    for $ix (0..$#p) {
        $im->filledEllipse(scl($p[$ix]), 6, 6, (vec($m, $ix, 1) ? $green : $gray));
    }

    $ix = -1;
    while (defined ($ix = bu_first($s, $ix + 1))) {
        $im->filledEllipse(scl($p[$ix]), 6, 6, $red);
    }

    my $mis = bu_count($m);
    my $sol = bu_count($s) + $mis;

    $im->string(gdSmallFont, 4, 4, "iteration: $iteration, sol: $sol, missing: $mis, fit: $fit, len: $len, p: $p", $black);

    my $name = sprintf("cover-%05d.png", $iteration);
    open my $fh, '>', $name;
    print $fh $im->png;
    close $fh;
};

