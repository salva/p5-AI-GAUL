#!/usr/bin/perl

use 5.010;

use strict;
use warnings;

use Algorithm::GAUL;
use Math::Vector::Real;
use Math::Vector::Real::Random;
use Math::Vector::Real::kdTree;
use GD;
use Devel::Peek;
use Sort::Key qw(nkeysort);

$| = 1;
use Getopt::Long;

my $width = 1000;
my $height = 800;
my $cookies = 10;
my $pop_size = 2000;

GetOptions("width|w=i", \&width,
           "height|h=i", \&height,
           "cookies|n=i", \&cookies,
           "pop_size|s=i", \&pop_size) or die "Error in command line arguments";

my $b0 = V(0, 0);
my $b1 = V($height, $width);

sub chromo2vs {
    my @p = unpack 'd*' => $_[0];
    my $len = @p / 2;
    map {
        my $ix = $_ * 2;
        [@p[$ix, $ix+1]];
    } 0 .. $len - 1;
}

sub vs2chromo { pack 'd*', map @$_, @_ }

sub weight {
    my @v = chromo2vs($_[1][0]);
    my $t = Math::Vector::Real::kdTree->new(@v);
    my $d = $t->find_two_nearest_vectors;
    for my $v (@v) {
        $d = $v->[0] if $v->[0] < $d;
        $d = $v->[1] if $v->[1] < $d;
        $d = $width - $v->[0] if $width - $v->[0] < $d;
        $d = $height - $v->[1] if $height - $v->[1] < $d;
    }
    ($d > 0 ? $d : 0);
}

sub crossover {
    my ($o0, $o1) = map V((0, $width)[rand 2], (0, $height)[rand 2]), 0, 1;
    my @v0 = nkeysort { $o0->dist2($_) } chromo2vs($_[1][0]);
    my @v1 = nkeysort { $o1->dist2($_) } chromo2vs($_[2][0]);
    my $cut = int(rand $cookies);
    my @c0 = (splice(@v0, 0, $cut), splice(@v1, 0, $cookies - $cut));
    my @c1 = (@v0, @v1);
    $_[3][0] = vs2chromo(@c0);
    $_[4][0] = vs2chromo(@c1);
}

#sub mutate {
#    print "mutate!\n";
#    Dump $_[1][0];
#}
#sub crossover {
#    print "crossover!\n";
#    Dump $_[1][0];
#    Dump $_[2][0];
#}


#sub repair {}

sub seed {
    $_[1][0] = pack "d*", map @{Math::Vector::Real::random_in_box(V($width, $height))}, 1..$cookies;
}

my $gaul = Algorithm::GAUL->new(len_chromo => 2 * $cookies,
                                type_chromo => 'double',
                                population_size => $pop_size,
                                select_two => 'aggressive',
                                selec_one => 'aggressive',
                                mutation_ratio => 0.5,
                                seed => \&seed,
                                evaluate => \&weight,
                                mutate => 'singlepoint_drift',
                                elitism => 'one_parent_survives',
                                scheme => 'darwin',
                                crossover => \&crossover,
                               );

say "iterating...";

my $iteration = 0;
while (1) {
    $iteration++;
    $gaul->evolution(1);
    my $fit = $gaul->fitness_by_rank(0);
    say "iteration: $iteration, fit: $fit";

    my $im = GD::Image->new($width, $height, 0);
    my $white = $im->colorAllocate(255, 255, 255);
    $im->filledRectangle(0, 0, $width, $height, $white);

    my $black = $im->colorAllocate(0, 0, 0);
    my $gray = $im->colorAllocate(200, 200, 200);
    my $green = $im->colorAllocate(0, 150, 0);
    my $red = $im->colorAllocate(255, 0, 0);
    my $yellow = $im->colorAllocate(255, 255, 0);

    my @v = chromo2vs($gaul->entity_by_rank(0)->[0]);
    my $d = $gaul->fitness_by_rank(0);

    for my $v (@v) {
        #print "ellipse $v->[0], $v->[1] => $d\n";
        $im->filledEllipse($v->[0], $v->[1], $d, $d, $green);
    }
    #$im->filledEllipse(scl($p[$ix]), $scl_dist, $scl_dist, $yellow);

    my $next = join(' ', map { $gaul->fitness_by_rank($_) } 1..25);
    $im->string(gdSmallFont, 4, 4, "iteration: $iteration, d: $d, next: $next", $black);

    my $name = sprintf("cover-%05d.png", $iteration);
    open my $fh, '>', $name;
    print $fh $im->png;
    close $fh;
};

