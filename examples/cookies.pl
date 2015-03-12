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

GetOptions("width|w=i",    \$width,
           "height|h=i",   \$height,
           "cookies|n=i",  \$cookies,
           "pop_size|s=i", \$pop_size) or die "Error in command line arguments";

my $box = V($height, $width);

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

sub mutate {
    my $lin = length $_[1][0];
    my @v = chromo2vs($_[1][0]);
    if (rand > 0.5) {
        my $t = Math::Vector::Real::kdTree->new(@v);
        my ($ix0, $ix1, $d) = $t->find_two_nearest_vectors;
        my $ix = ($ix0, $ix1)[rand 2];
        $d *= 0.5;
        for my $i (0 .. $#v) {
            my $v = $v[$i];
            for my $ds ($v->[0], $v->[1], $width - $v->[0], $height - $v->[1]) {
                if ($ds < $d) {
                    $d = $ds;
                    $ix = $i;
                }
            }
        }
        $d = 0 if $d < 0;
        $v[$ix] = Math::Vector::Real->random_in_box($box);
    }
    else {
        for (0 .. rand(3) + 1) {
            my $ix = int rand @v;
            my $n =  Math::Vector::Real->random_in_box($box);
            $v[$ix] = $n;
            printf "replacing vector at %d: %s (%s)\n", $ix, V(@{$v[$ix]}), $n;
        }
    }
    my $out = vs2chromo(@v);
    my $lout = length $out;
    if ($lin != $lout) {
        my @out = chromo2vs($out);
        for (0..$#v) {
            my $v0 = V(@{$v[$_]});
            my $v1 = V(@{$out[$_]});
            if ("$v0" ne "$v1") {
                print "vectors at $_ differ: $v0, $v1, $box\n";
            }
        }
    }
    $_[1][0] = $out;
}

sub min_dists {
    my $t = Math::Vector::Real::kdTree->new(@_);
    my @near = $t->find_nearest_vector_all_internal;
    my @d;
    if (@near == @_ and not grep !defined, @near) {
        my $worst_d2 = 'inf' + 0;
        my $sum_sqrt_d = 0;
        for my $i (0..$#_) {
            my $v = $_[$i];
            my $d = 0.5 * Math::Vector::Real::dist($v, $_[$near[$i]]);
            $d = $v->[0] if $v->[0] < $d;
            $d = $v->[1] if $v->[1] < $d;
            $d = $width - $v->[0] if $width - $v->[0] < $d;
            $d = $height - $v->[1] if $height - $v->[1] < $d;
            $d = 0 if $d < 0;
            $worst_d2 = $d if $d < $worst_d2;
            push @d, $d;
        }
        return ($worst_d2, @d);
    }
    else {
        warn "find_nearest_vector_all_internal failed: [@near]\n".$t->dump_to_string;
        return ((0) x @near);
    }
}

sub weight_2 {
    my @v = chromo2vs($_[1][0]);
    my ($worst_d, @d) = min_dists(@v);
    my $sum_d = 0;
    $sum_d += $_ for @d;
    return $worst_d + $sum_d / @v;
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
    $_[1][0] = pack "d*", map @{Math::Vector::Real::random_in_box($box)}, 1..$cookies;
}

my $gaul = Algorithm::GAUL->new(len_chromo => 2 * $cookies,
                                type_chromo => 'double',
                                population_size => $pop_size,
                                select_two => 'roulette',
                                selec_one => 'roulette',
                                mutation_ratio => 0.5,
                                seed => \&seed,
                                evaluate => \&weight_2,
                                mutate => \&mutate, # singlepoint_drift',
                                elitism => 'parents_survive',
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
    my ($worst_d, @d) = min_dists(@v);

    for my $i (0 .. $#v) {
        my $v = $v[$i];
        my $d = $d[$i];
        $im->filledEllipse($v->[0], $v->[1], $d * 2, $d * 2, $yellow);
    }
    for my $i (0 .. $#v) {
        my $v = $v[$i];
        $im->filledEllipse($v->[0], $v->[1], $worst_d * 2, $worst_d * 2, $green);
    }
    
    #$im->filledEllipse(scl($p[$ix]), $scl_dist, $scl_dist, $yellow);

    my $next = join(' ', map { $gaul->fitness_by_rank($_) } 1..25);
    $im->string(gdSmallFont, 4, 4, "iteration: $iteration, worst_d: $worst_d, next: $next", $black);

    my $name = sprintf("cover-%05d.png", $iteration);
    open my $fh, '>', $name;
    print $fh $im->png;
    close $fh;
};

