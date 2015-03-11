package Algorithm::GAUL;

our $VERSION = '0.01';

use strict;
use warnings;

require XSLoader;
XSLoader::load('Algorithm::GAUL', $VERSION);

sub new {
    my ($class, %opts) = @_;
    @_ = ($class, \%opts);
    goto &genesis_any;
}

1;
__END__

=head1 NAME

Algorithm::GAUL - Perl wrapper for the Genetic Algorithm Utility Library

=head1 SYNOPSIS

  use Algorithm::GAUL;

  sub evaluate {
    my ($pop, $chromo) = @_;
    # count bits in chromosome
    (unpack "b*", $chromo) =~ tr/1//;
  }

  my $pop = Algorithm::GAUL->new(len_chromo => 100,
                          population_size => 1000,
                          evaluate => \&evaluate);

  $pop->evolution(100);

  my $best    = $pop->chromosomes_by_rank(0);
  my $fitness = $pop->fitness_by_rank(0);

  my $best_bin = unpack "b*", $best;

  say "best: $best_bin, fitness: $fitness";

=head1 DESCRIPTION

Algorithm::GAUL is a perl wrapper for the Genetic Algorithm Utility Library
(GAUL).

 WARNING:

   This is an early release with limited functionality and probably
   tons of bugs.

   The API should not be considered stable either.

=head2 Why another GA in Perl?

I wrote Algorithm::GAUL because I wanted an efficient implementation
of GA for Perl.

Using the GAUL library as the base, has also resulted in a easy to
use, flexible and feature rich solution.



=head2 Chromosome representation

Internally the base unit for the GA implemented by GAUL is not the
chromosome but the "entity" that is a data structure holding several
chromosomes (num_chromo). Though, in most cases, you will be using
entities containing just one chromosome.

Entities are never visible at the Perl level. Instead, the set of
chromosomes in the entity is passed directly to perl. For instance,
when callback functions operating over the entities are called, for
every chromosome in the entity a perl scalar is passed as an argument
to the callback.

Every chromosome, has a given number of alleles (len_chromo) that can
be of type bitstring, boolean, char, integer or double (boolean
alleles are integers accepting only values 0 or 1, probably this type
is useless in Perl as the bitstring type has the same exact behavior
and uses far less memory).

The chromosomes are always passed to Perl as packed scalars
independently of the type of the alleles. Perl C<unpack> builtin (or
C<vec> for bitstrings) can be used to extract the data.

=head1 API

=over 4

=item $pop = Algorithm::GAUL->new(%opts)

Creates a new object of type Algorithm::GAUL (AKA the population object)

The options accepted are:

=over 4

=item population_size => $n

The total number of individuals in the population.

This argument is mandatory.

=item len_chromo => $n

Number of alleles per chromosome

This argument is mandatory.

=item type_chromo => $type

Allele type.

Valid values are C<bitstring>, C<boolean>, C<integer>, C<double> and C<char>. Default is C<bitstring>

=item num_chromo => $n

Number of chromosomes per individual, for most cases, the default, 1, is probably what you want.

=item evalutate => \&cb

Evaluation callback. It is called as follows:

  evaluate($pop, @chromosomes)

and must return the fitness of the solution or undef if it is invalid

This argument is mandatory.

=item adapt => \&cb

Adapt callback, used for lamarckian evolution. Called as follows:

  adapt($pop, @chromosomes)

The callback should modify the chromosomes passed into @_ inplace

=item stop => \&cb

Callback that allows to stop the evolution when some condition is
meet. Called as follows:

  stop($pop, $iteration)

Where C<$iteration> is the iteration number.

If the callback returns a true value, the evolution is stopped.

This callback can also be used to perform some operation between
iterations of the algorithm. For instance printing some information to
the user or modifying parameters used in other callbacks.

=item scheme => $scheme

Evolutionary scheme indicates the form of adaptation that should be
used. C<darwin> relates to a standard GA, while the others relate to
so-called hybrid genetic algorithms.

Valid values are C<darwin>, C<baldwin_all>, C<baldwin_children>,
C<baldwin_parents>, C<lamarck_all>, C<lamarck_children> and
C<lamark_parents>. Default is C<darwin>.

This excerpt from the GAUL documentation explains the different
schemes:

  Adaptation

  During the lifetime of an individual in the natural world, its
  fitness may increase, either through "physical adaptation"
  (e.g. bigger muscles resulting from exercise) or through "learnt
  adaptations" (e.g. learning how to find plentiful food supplies). In
  GAUL, there is no distinction between the two types of adaptation,
  it simply provides a mechanism with which the programmer may
  reevaluate solutions. (Of course, the programmer may consider the
  two mechanisms for adaptation seperately if they wish.)

  Additionally, in our computational optimisation procedures we are
  not resticted to following the exact procedure from nature, whatever
  that may be. This provides us with the flexibility to produce,
  so-called, hybrid-GA algorithms.

  "The Baldwin Effect"

  If learning helps survival or procreation, then the individuals
  which are best able to learn will have the most
  offspring. Therefore, the genes responsible for learning will be
  more frequent in future generations.

  "The Lamarkian Hypothesis"

  The (largely discredited, but well-known) Lamarckian hypothesis
  states that traits aquired during the lifetime of an organism may be
  genetically transmitted to its offspring. Although the required
  mechanism for "reverse transcription" of aquired traits into the
  individual's genes does not exist in nature, computationally it may
  form the basis of very powerful optimisation procedures.


=item elitism => $elitism

The elitism paramter indicates how parents are passed into subsequent
generations.

Accepted values are as follows:

=over

=item parents_survive

All parents that rank sufficiently highly will pass to the next
generation.

=item one_parent_survives

Single fittest parent will pass to the next generation if it ranks
sufficiently well.

=item parents_die

No parents pass to next generation, regardless of their fitness.

=item rescore_parents

All parents are re-evalutated, and those that subsequently rank
sufficiently highly will pass to the next generation.

=back

Note that some elitism options make little sense when coupled to
certian evolutionary schemes. For example, using Lamarckian or
Baldwinian evolution and applying that to the parent individuals is
pointless when none of those parents will ever pass into subsequent
generations.

=item crossover_ratio => $prob

=item mutation_ratio => $prob

=item migration_ratio => $prob

=item allele_migration_ratio => $prob

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

=item crossover => $cb

Accepted values are C<allele_mixing>, C<doublepoints>, C<mixing> and
C<singlepoints>.

For alleles of type C<double> and C<integer>, the value C<mean> is
also valid.

=item mutate => $cb

Accepted values are C<multipoint> and C<singlepoint>.

For types C<char>, C<integer> and C<double>, the values
C<singlepoint_drift> and S<multipoint_drift> are also valid.

For type C<char>, the values C<printable_allpoint>,
C<printable_multipoint>, C<printable_singlepoint_drift> and
C<printable_singlepoint_randomize> are also valid.

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

The documentation of the GAUL library L<http://gaul.sourceforge.net/>.

For a mature and flexible implementation of genetic algorithms in
Perl, see L<Algorithm::Evolutionary>.

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by Salvador Fandiño E<lt>sfandino@yahoo.comE<gt>

Parts of this module documentation has been copied or adapted from the
GAUL docs, Copyright (C) 2000-2009 by Stewart Adcock.

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.1 or,
at your option, any later version of Perl 5 you may have available.


=cut
