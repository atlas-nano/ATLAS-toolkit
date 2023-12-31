use 5.008008;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'PDBsub',
    VERSION_FROM      => 'lib/PDBsub.pm', # finds $VERSION
    PREREQ_PM         => {}, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/PDBsub.pm', # retrieve abstract from module
       AUTHOR         => 'Tohru Terada <tterada@iu.a.u-tokyo.ac.jp>') : ()),
    LIBS              => ['-L/home/tpascal/Software/LAPACK -llapack -lrefblas -lm'], # e.g., '-lm'
    LD                => 'gfortran',
    DEFINE            => '', # e.g., '-DHAVE_SOMETHING'
    INC               => '-I.', # e.g., '-I. -I/usr/include/other'
	# Un-comment this if you add C files to link with later:
    # OBJECT            => '$(O_FILES)', # link all the C files too
);
