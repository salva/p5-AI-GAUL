#!/usr/bin/perl

use 5.010;

use strict;
use warnings;

use Algorithm::GAUL;
use Math::Vector::Real;
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

my $t = Math::Vector::Real::kdTree->new(@p);
@p = @p[$t->ordered_by_proximity];

my @cover;
$t = Math::Vector::Real::kdTree->new(@p);
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

sub weight {
    say "weighting";
    -bu_count($_[1]);
}

my $missing = '';
vec($missing, $_, 1) = 1 for (0..$#p);

sub weight1 {
    for my $s ($_[1]) {
        my $m = $missing;
        my $ix = -1;
        while (defined($ix = bu_first($s, $ix + 1))) {
            vec($m, $_, 1) = 0 for @{$cover[$ix]};
        }
        return @p - (bu_count($m) + bu_count($s));
    }
}

sub stop {
    say "stop?";
    return 0;
}

sub repair {
    # say "repairing chromosome";
    for my $s ($_[1]) {
        my $m = $missing;
        my $ix = -1;
        while (defined($ix = bu_first($s, $ix + 1))) {
            vec($m, $_, 1) = 0 for @{$cover[$ix]};
        }
        $ix = -1;
        while (defined($ix = bu_first($m, $ix + 1))) {
            vec($s, $ix, 1) = 1;
            vec($m, $_, 1) = 0 for @{$cover[$ix]};
        }
    }
}

sub mutate {
    # say "mutate";
    for (0..9) {
        my $ix = int rand scalar @p;
        if (vec($_[1], $ix, 1)) {
            vec($_[1], $ix, 1) = 0;
            my $c = $cover[$ix];
            vec($_[1], $c->[rand @$c],  1) = 1;
        }
    }
}

my $gaul = Algorithm::GAUL->new(len_chromo => scalar(@p),
                                population_size => $pop_size,
                                select_two => 'linearrank',
                                mutation_ratio => 0.05,
                                seed => 'random',
                                evaluate => \&weight1,
                                adapt => \&repair,
                                stop => \&stop,
                                mutate => \&mutate,
                                elitism => 'parents_survive',
                                scheme => 'lamarck_children',
                                crossover => 'doublepoints');


sub scl {
    my $p = shift;
    @{0.25 * $width * ($p + [2,2])};
}

my $iteration = 0;
while (1) {
    $iteration++;
    say "iteration $iteration";
    $gaul->evolution(1);
    say "iterating...";
    say "iteration: $iteration, fit: ", $gaul->fitness_from_rank(0);

    my $s = $gaul->chromosomes_from_rank(0);

    my $len = length($s) * 8;
    my $p = scalar @p;

    my $im = GD::Image->new($width, $width, 0);
    my $white = $im->colorAllocate(255, 255, 255);
    $im->filledRectangle(0, 0, $width, $width, $white);

    my $black = $im->colorAllocate(0, 0, 0);
    my $gray = $im->colorAllocate(200, 200, 200);
    my $green = $im->colorAllocate(0, 150, 0);
    my $red = $im->colorAllocate(255, 0, 0);

    my $m = $missing;
    my $ix = -1;

    my $s2 = 0;
    while (defined($ix = bu_first($s, $ix + 1))) {
        $s2++;
        vec($m, $_, 1) = 0 for @{$cover[$ix]};
    }

    my $s3 = 0;
    for my $ix (0..$#p) {
        my $color = (vec($s, $ix, 1) ? ($s3++, $red) :
                     vec($m, $ix, 1) ? $green   :
                                       $gray  );
        $im->filledEllipse(scl($p[$ix]), 6, 6, $color);
    }

    my $mis = bu_count($m);
    my $sol = bu_count($s) + $mis;

    $im->string(gdSmallFont, 4, 4, "iteration: $iteration, sol: $sol, missing: $mis, s2: $s2, s3: $s3, len: $len, p: $p", $black);

    my $name = sprintf("cover-%05d.png", $iteration);
    open my $fh, '>', $name;
    print $fh $im->png;
    close $fh;
};

