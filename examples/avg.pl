#!/usr/bin/perl

use strict;
use warnings;

use 5.010;

use Bit::Grep qw(bg_count_sum_and_sum2 bg_count_and_sum);
use Bit::Util qw(bu_first bu_last);
use List::Util qw(sum);
use Algorithm::GAUL;

my $fn = shift // die "test data file name missing";

my $n = shift // 1000;
my $bottom = shift // 1000;
my $top = shift // 1001;

print "n: $n, bottom: $bottom, top: $top\n";

open my $fh, '<', $fn or die "unable to open $fn";
my @size = map { int $_ } grep /^\d/, <$fh>;
close $fh;

$#size = $n - 1;
@size = sort { $a <=> $b } @size;

my $size = scalar @size;

my $naive = '';
for (0..$#size) { vec($naive, $_, 1) = 1 };
adapt(undef, $naive);

my ($cnt, $sum) = bg_count_and_sum($naive, @size);

my $inside = grep { $_ >= $bottom and $_ <= $top } @size;

if ($cnt) {
    say "naive solution cnt: $cnt, inside: $inside";
}
else {
    say "no naive solution";
}

my $last_pow = 0;

my $iteration = 1;
my $pow = 1;

sub evaluate {
    say "evaluating...";
    my ($cnt, $sum, $sum2) = bg_count_sum_and_sum2($_[1], @size);
    if ($cnt) {
        my $avg = $sum / $cnt;
        if ($avg >= $bottom and $avg <= $top) {
            my $dev = ($sum2 - $avg * $avg) / $cnt;
            return $cnt * $dev ** $pow;
        }
        return undef;
    }
    return 0;
}

sub stop {
    $iteration++;
    $pow = 0.5 * $iteration ** -0.9;
    undef;
}

sub find_le_ix {
    my $n = shift;
    my $l = -1;
    my $r = $#size;
    while ($l < $r) {
        my $p = int(($l + $r + 1) / 2);
        if ($size[$p] > $n) {
            $r = $p - 1;
        }
        else {
            $l = $p;
        }
    }
    $l;
}

sub find_ge_ix {
    my $n = shift;
    my $l = 0;
    my $r = @size;
    while ($l < $r) {
        my $p = int(($l + $r) / 2);
        if ($size[$p] < $n) {
            $l = $p + 1;
        }
        else {
            $r = $p;
        }
    }
    $l;
}

sub adapt {
    say "adapting...";
    my ($cnt, $sum) = bg_count_and_sum($_[1], @size);
    while ($cnt) {
        my $avg = ($sum / $cnt);
        if ($avg < $bottom) {
            my $ix = find_le_ix($sum - $bottom * ($cnt - 1));
            my $ix1 = bu_last($_[1], $ix) // bu_first($_[1], $ix) // die "internal error";
            vec($_[1], $ix1, 1) = 0;
            $sum -= $size[$ix1];
            $cnt--;
            # say "element at $ix1 ($size[$ix1]) deleted avg: $avg, cnt: $cnt, dir: <==";

        }
        elsif ($avg > $top) {
            my $ix = find_ge_ix($sum - $top * ($cnt - 1));
            # my $last = bu_last($_[1]);
            my $ix1 = bu_first($_[1], $ix) // bu_last($_[1], $ix) // die "internal error";
            vec($_[1], $ix1, 1) = 0;
            $sum -= $size[$ix1];
            $cnt--;
            # say "element at $ix1 ($size[$ix1]) deleted, last: $last, avg: $avg, cnt: $cnt, dir: ==>", ($last != $ix1 ? '*': '');
        }
        else {
            # say "range reached, avg: $avg, cnt: $cnt";
            return
        }
    }
}

my $pop = Algorithm::GAUL->new(len_chromo => $size,
                               population_size => 300,
                               select_two => "linearrank",
                               mutate_ratio => 0.2,
                               mutate => "multipoint",
                               crossover => "allele_mixing",
                               evaluate => \&evaluate,
                               adapt => \&adapt,
                               stop => \&stop,
                               scheme => "lamarck_children",
                               elitism => "rescore_parents",
                               seed => "zero",
                               #seed => sub { $_[1] = $naive }
                              );

while (1) {
    $pop->evolution(1);
    my $chromo = $pop->chromosomes_from_rank(0);
    my ($cnt, $sum, $sum2) = bg_count_sum_and_sum2($chromo, @size);
    $cnt ||= 1;
    my $avg = $sum / $cnt;
    my $dev = ($sum2 - $avg * $avg) / $cnt;
    say "$iteration $cnt/$size fit: ", $pop->fitness_from_rank(0), " avg: ", $avg, " dev1/2: ", sqrt($dev), " pow: ", $pow, " fact: ", $dev ** $pow;
    
}
