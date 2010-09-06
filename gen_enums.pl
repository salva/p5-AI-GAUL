#!/usr/bin/perl

use 5.010;

use strict;
use warnings;

use Data::Dumper;

my %enum;

my %cb = map { $_ => [] } qw(ga_select_one ga_select_two
                             ga_replace ga_generation_hook
                             ga_evaluate ga_adapt
                             ga_seed_bitstring ga_seed_char ga_seed_boolean
                             ga_seed_integer ga_seed_double
                             ga_mutate_bitstring ga_mutate_char ga_mutate_boolean
                             ga_mutate_integer ga_mutate_double
                             ga_crossover_bitstring ga_crossover_char ga_crossover_boolean
                             ga_crossover_integer ga_crossover_double
                           );

push @{$cb{$_}}, "_perl" for qw(ga_evaluate);

my $base = '';
my $table = '';

while (<>) {
    given ($_) {
        when (/typedef\s+enum\s+(ga_(\w+)_type_t)/) {
            $table = $1;
            $base = $2;
        }
        when (/typedef\s+enum\s+((\w+)_t)/) {
            $table = $1;
            $base = $2;
        }
        when (/}\s+ga_${base}_type/) {
            $base = '';
        }
        when (/(\w+)\s*=\s*(\d+)/) {
            if (length $base) {
                my ($name, $value) = (lc $1, $1);
                if ($name =~ /^ga_${base}_(\w+)/) {
                    $enum{$table}{$1} = $value;
                }
                else {
                    warn "enum entry $name does not match $base";
                }
            }
        }
        when (/FUNCPROTO(?:\s+\w+)+\s+(ga_\w+)/) {
            my $func = $1;
            for my $k (keys %cb) {
                if ($func =~ /^${k}_(\w+)/) {
                    push @{$cb{$k}}, $1;
                }
                my $k1 = $k;
                if ($k1 =~ s/_char$/_printable/) {
                    if ($func =~ /^${k1}_(\w+)/) {
                        push @{$cb{$k}}, "printable_$1";
                    }
                }
            }
        }
    }
}

open my $enum_out, '>', 'enum-tables.h' or die $!;

for my $enum_name (sort keys %enum) {
    my $enum = $enum{$enum_name};
    print $enum_out "static struct enum_table table_for_${enum_name}[] = {\n";
    for my $name (sort keys %$enum) {
        print $enum_out qq(    { "$name", $enum->{$name} },\n);
    }
    print $enum_out "    { NULL, 0 },\n};\n\n";
}

open my $cb_out, '>', 'cb-tables.h' or die $!;

for my $cb_name (sort keys %cb) {
    my @cb = sort @{$cb{$cb_name}};
    print $cb_out "static struct cb_table table_for_${cb_name}[] = {\n";
    for my $name (@cb) {
        my $cb = "${cb_name}_$name";
        $cb =~ s/_char_printable_/_printable_/;
        print $cb_out qq(    { "$name", (cb_C)&$cb },\n);
    }
    print $cb_out "    { NULL, 0 },\n};\n\n";
}

