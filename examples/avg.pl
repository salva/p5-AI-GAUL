#!/usr/bin/perl

use strict;
use warnings;

use 5.010;

use Bit::Grep qw(bg_count_and_sum);
use Bit::Util qw(bu_first bu_last);
use List::Util qw(sum);
use AI::GAUL;

my $fn = shift // die "test data file name missing";

open my $fh, '<', $fn or die "unable to open $fn";
my @size = map { int $_ } grep /^\d/, <$fh>;
close $fh;

$#size = 999;
@size = sort { $a <=> $b } @size;

my $size = scalar @size;

my ($bottom, $top) = (4000, 4100);

my $first = 0;
my $last = $#size;

my @cp = @size;
my $sum = sum @size;
while (@cp) {
    my $avg = $sum / @cp;
    if ($avg < $bottom) {
        # say "deleting $cp[0]";
        $sum -= shift @cp;
        $first++;
    }
    elsif ($avg > $top) {
        # say "deleting $cp[-1]";
        $sum -= pop @cp;
        $last--;
    }
    else {
        last;
    }
}

my $inside = grep { $_ >= $bottom and $_ <= $top } @size;

if (@cp) {
    say "naive solution: ", scalar @cp, " inside: $inside";
}
else {
    say "no naive solution";
}

my $seed = '';
vec($seed, $_, 1) = ($_ >= $first and $_ <= $last) for 0..$#size;

sub evaluate {
    my ($cnt, $sum) = bg_count_and_sum($_[1], @size);
    if ($cnt) {
        my $avg = $sum / $cnt;
        $avg >= $bottom and $avg <= $top and return $cnt;
        # warn "average $avg out of limits";
        return undef;
    }
    return 0;
}

sub adapt {
    my ($cnt, $sum) = bg_count_and_sum($_[1], @size);
    if ($cnt) {
        my $avg = ($sum / $cnt);
        if ($avg < $bottom or $avg > $top) {
            my $first = bu_first($_[1]);
            my $last = bu_last($_[1]);

            while (1) {
                my $avg = ($sum / $cnt);
                if ($avg < $bottom) {
                    vec($_[1], $first, 1) = 0;
                    $sum -= $size[$first];
                    $first = bu_first($_[1], $first+1);
                    return if $first < 0;
                }
                elsif ($avg > $top) {
                    vec($_[1], $last, 1) = 0;
                    $sum -= $size[$last];
                    $last = bu_last($_[1], $last - 1);
                    return if $last < 0;
                }
                else {
                    return;
                }
                $cnt--;
            }
        }
    }
}

my $pop = AI::GAUL->new(len_chromo => $size,
                        population_size => 2000,
                        seed => "random",
                        mutate_prob => 0.01,
                        mutate => "multipoint",
                        evaluate => \&evaluate,
                        adapt => \&adapt,
                        scheme => "lamarck_children",
                        # seed => sub { $_[1] = $seed }
                       );

while (1) {
    $pop->evolution(1);
    say $pop->fitness_from_rank(0), "/$size";
}
