use 5.010001;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile( NAME              => 'AI::GAUL',
               VERSION_FROM      => 'lib/AI/GAUL.pm', # finds $VERSION
               PREREQ_PM         => {}, # e.g., Module::Name => 1.1
               AUTHOR         => 'Salvador Fandiño <salva@>',
               LIBS              => ['-L/usr/local/gaul/lib -lgaul -lgaul_util'], # e.g., '-lm'
               DEFINE            => '', # e.g., '-DHAVE_SOMETHING'
               INC               => '-I. -I/usr/local/gaul/include', # e.g., '-I. -I/usr/include/other'
               OPTIMIZE => '-g -O0',
             );
