use 5.010001;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile( NAME              => 'Algorithm::GAUL',
               VERSION_FROM      => 'lib/Algorithm/GAUL.pm', # finds $VERSION
               PREREQ_PM         => {},
               AUTHOR         => 'Salvador Fandiño <sfandino@yahoo.com>',
               LIBS              => ['-L/usr/local/gaul/lib -lgaul -lgaul_util'], # e.g., '-lm'
               INC               => '-I. -I/usr/local/gaul/include', # e.g., '-I. -I/usr/include/other'
               # OPTIMIZE => '-g -O0',
             );
