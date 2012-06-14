#!/usr/bin/perl

use strict;
use warnings;

use Math::Vector::Real;
use Math::Vector::Real::Random;

my ($d, $n) = @ARGV;
print "# dim: $d, n: $n\n";

for (1..$n) {
    print Math::Vector::Real->random_normal($d), "\n"
}

