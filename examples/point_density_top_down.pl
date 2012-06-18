#!/usr/bin/perl

use 5.010;

use strict;
use warnings;

use Math::Vector::Real;
use Math::Vector::Real::Random;
use Math::Vector::Real::kdTree;
use Bit::Util qw(bu_count bu_first);
use GD;

$| = 1;

my $dist = shift @ARGV || 1;
my $pop_size = shift @ARGV || 100;
my @p;

while (<>) {
    if (my ($v) = /\{\s*([^\s,]+(?:\s*,\s*[^\s,]+)*)\s*\}/) {
        push @p, V(split /\s*,\s*/, $v);
    }
}

my $width = 1024;
my $scl_dist = $dist * 0.125 * 2 * $width;
say "scl_dist: $scl_dist";
sub scl {
    my $p = shift;
    @{0.125 * $width * ($p + [4,4])};
}



my $t = Math::Vector::Real::kdTree->new(@p);
my @cover = map [$t->find_in_ball($_, $dist)], @p;
say "cover generated";

my @n = map scalar(@$_), @cover;
my $s = '';
vec($s, $_, 1) = 1 for (0..$#p);

my $iteration = 0;
OUT: while (1) {
    $iteration++;
    my $ix = bu_first($s);
    my $other = $t->find_nearest_neighbor_internal($ix);
    my $min_d2 = $p[$ix]->dist2($p[$other]);
    my $min_ix = $ix;
    my $min_other = $other;
    while (defined($ix = bu_first($s, $ix + 1))) {
        my $other = $t->find_nearest_neighbor_internal($ix);
        my $d2 = $p[$ix]->dist2($p[$other]);
        if ($d2 < $min_d2) {
            $min_d2 = $d2;
            $min_ix = $ix;
            $min_other = $other;
        }
    }
    draw($s, $min_ix, $min_other) unless $iteration % 100;
    my $victim = ($min_ix, $min_other)[rand 2];
    say "victim: $victim, ix: $min_ix, other: $min_other, $p[$min_ix], $p[$min_other]";
    vec($s, $victim, 1) == 1 or say "$victim is not in selected set";
    for (@{$cover[$victim]}) {
        last OUT if --$n[$_] <= 0;
    }
    $t->hide($victim);
    # use Data::Dumper;
    # print Dumper $t->{tree};
    vec($s, $victim, 1) = 0;
}
draw($s);

sub draw {
    my $s = shift;
    say "iteration: $iteration";

    my $im = GD::Image->new($width, $width, 0);
    my $white = $im->colorAllocate(255, 255, 255);
    $im->filledRectangle(0, 0, $width, $width, $white);

    my $black = $im->colorAllocate(0, 0, 0);
    my $gray = $im->colorAllocate(200, 200, 200);
    my $green = $im->colorAllocate(0, 150, 0);
    my $red = $im->colorAllocate(255, 0, 0);
    my $yellow = $im->colorAllocate(255, 255, 0);

    my $ix = -1;
    while (defined($ix = bu_first($s, $ix + 1))) {
        $im->filledEllipse(scl($p[$ix]), $scl_dist, $scl_dist, $yellow);
    }
    for (@p) {
        $im->filledEllipse(scl($_), 6, 6, $gray);
    }
    $ix = -1;
    while (defined ($ix = bu_first($s, $ix + 1))) {
        my $other = $t->find_nearest_neighbor_internal($ix);
        $im->line(scl($p[$ix]), scl($p[$other]), $black);
        $im->filledEllipse(scl($p[$ix]), 6, 6, (grep($_ == $ix, @_) ? $green : $red));
    }
    my $sol = bu_count($s);
    $im->string(gdSmallFont, 4, 4, "top_down, iteration: $iteration, sol: $sol", $black);

    my $name = sprintf("top_down-%05d.png", $iteration);
    open my $fh, '>', $name;
    print $fh $im->png;
    close $fh;
};
