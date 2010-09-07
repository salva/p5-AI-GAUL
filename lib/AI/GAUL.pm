package AI::GAUL;

our $VERSION = '0.01';

use strict;
use warnings;

require XSLoader;
XSLoader::load('AI::GAUL', $VERSION);

sub new {
    my ($class, %opts) = @_;
    $class->genesis_any(\%opts);
}

1;
__END__

=head1 NAME

AI::GAUL - Perl wrapper for the Genetic Algorithm Utility Library

=head1 SYNOPSIS

  use AI::GAUL;

  sub evaluate {
    my ($pop, $chromo) = @_;
    # count bits in chromosome
    (unpack "b*", $chromo) =~ tr/1//;
  }

  my $pop = AI::GAUL->new(len_chromo => 100,
                          population_size => 1000,
                          evaluate => \&evaluate);

  $pop->evolution(100);

  my $best    = $pop->chromosomes_by_rank(0);
  my $fitness = $pop->fitness_by_rank(0);

  my $best_bin = unpack "b*", $best;

  say "best: $best_bin, fitness: $fitness";

=head1 DESCRIPTION

AI::GAUL is a perl wrapper for the Genetic Algorithm Utility Library
(GAUL).

 WARNING:

   This is an early release with limited functionality and probably
   tons of bugs.

   The API should not be considered stable either.

=head2 Why another GA in Perl?

The aim behind AI::GAUL was to have an efficient and easy to use GA
usable from Perl.

Using the GAUL library as the base, has also resulted in a flexible
and feature rich solution.

=head2 Chromosome representation

Internally the base unit for the GA implemented by GAUL is not the
chromosome but the "entity" that is a data structure holding a given
number of chromosomes (num_chromo).

Every chromosome, has a given number of alleles (len_chromo). Alleles
can be of type bitstring, boolean, char, integer or double (boolean
alleles are integers accepting only values 0 or 1, probably this
type is useless in Perl as the bitstring type has the same exact
behavior and uses far less memory).

Using entities instead of chromosomes as the base for the GA allows
some fine grainer control over the operators used on the evolution,
but for most practical cases it can be just ignored using entities
holding just one chromosome (this is the default).

Entities are never visible at the Perl level. Instead, the set of
chromosomes in the entity is passed directly to perl. For instance,
when callback functions operating over the entities are called, for
every chromosome in the entity a perl scalar is passed to the
callback in its argument list.

The chromosomes are always passed to Perl as packed scalars
independently of the type of the alleles. Perl C<unpack> builtin (or
C<vec> for bitstrings) can be used to extract the data.

=head1 API

=over 4

=item $pop = AI::GAUL->new(%opts)

Creates a new object of type AI::GAUL (AKA the population object)

The options accepted are:

=over 4

=item population_size => $n

This argument is mandatory.

=item len_chromo => $n

Number of alleles per chromosome

This argument is mandatory.

=item type_chromo => $type

Allele types. Valid values are C<bitstring>, C<boolean>, C<integer>, C<double> and C<char>.

Default is C<bitstring>

=item num_chromo => $n

Number of chromosomes, for most cases, the default, 1, is what you want.

=item evalutate => \&cb

Evaluation callback. It is called as follows:

  evaluate($pop, @chromosomes)

and must return the fitness of the solution or undef if it is invalid

=item adapt => \&cb

Adapt callback, used for lamarckian evolution. Called as follows:

  adapt($pop, @chromosomes)

The callback should modify the chromosomes in @_ inplace

=item stop => \&cb

Callback that allows to stop the evolution when some condition is meet. Called as follows:

  stop($pop, $iteration)

Where C<$iteration> is the iteration number.

If the callback returns a true value, the evolution is stopped.

=item scheme => $scheme

Valid values are C<darwin>, C<baldwin_all>, C<baldwin_children>,
C<baldwin_parents>, C<lamarck_all>, C<lamarck_children> and
C<lamark_parents>. Default is C<darwin>.

=item elitism => $elitism

Accepted values are C<one_parent_survives>, C<parents_die>,
C<parents_survive>, C<rescore_parents>, C<null>. Default is ???.

=item crossover_prob => $prob

=item mutation_prob => $prob

=item migration_prob => $prob

=item allele_migration_prob => $prob

=item allele_min => $min

=item allele_max => $max

minimum and maximum value for integer or double alleles

=item seed => $cb

callback to initialize the chromosomes.

Accepted values are C<random>, C<zero>.

For alleles of type C<char> the value C<printable_random> is also
valid.

For alleles of type C<double> the value C<random_unit_gaussian> is
also valid.

Currently, custom seed callbacks are not accepted (tell me if you need
them ;-)

=item crossover => $cb

Accepted values are C<allele_mixing>, C<doublepoints>, C<mixing> and
C<singlepoints>.

For alleles of type C<double> and C<integer>, the value C<mean> is
also valid.

Currently, custom crossover callbacks are not supported (tell me if
you need them ;-)

=item mutate => $cb

Accepted values are C<multipoint> and C<singlepoint>.

For types C<char>, C<integer> and C<double>, the values
C<singlepoint_drift> and S<multipoint_drift> are also valid.

For type C<char>, the values C<printable_allpoint>,
C<printable_multipoint>, C<printable_singlepoint_drift> and
C<printable_singlepoint_randomize> are also valid.

Currently, custom mutate callbacks are not supported (again, tell me
if you need them ;-)

=back

=item $pop->evolution($iterations)

Run the evolution algorithm for the given number of iterations.

If a C<stop> callback is also used, it will be able to stop the
algorithm before.

=item $pop->chromosomes_from_rank($rank)

Returns the chromosome/chromosomes that is/are ranked at the given position.

=item $pop->fitness_from_rank($rank)

Returns the fitness of the chromosome/chromosomes that is/are ranked
at the given position.

=back

=head1 SEE ALSO

For a mature and flexible implementation of genetic algorithms in
Perl, see L<Algorithm::Evolutionary>.

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by Salvador Fandiño E<lt>sfandino@yahoo.comE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.1 or,
at your option, any later version of Perl 5 you may have available.


=cut
