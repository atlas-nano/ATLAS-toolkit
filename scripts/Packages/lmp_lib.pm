package lmp_box;
use Carp;
local $Carp::CarpLevel = 1; # short messages
my ($cpack, $cfile) = caller();

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $data = {
		'LO' => [$_[0], $_[1], $_[2]],
		'HI' => [$_[3], $_[4], $_[5]],
		'XY' => $_[6],
		'YZ' => $_[7],
		'XZ' => $_[8],
		'P'  => [$_[9], $_[10], $_[11]],
		'C'  => $_[12],
	};
	my $self = sub : lvalue {
		my $field = shift;
		croak "Not a valid field '$field'"
			unless exists $data->{$field};
		if (@_) { $data->{$field} = shift }
		return $data->{$field};
	};
	bless($self,$class);
	return $self;
}
#generate method names
for my $field (qw(lo hi xy yz xz p c)) {
	no strict "refs"; #for access to the symbol table
	*$field = sub : lvalue {
		my $self = shift;
		return $self->(uc $field, @_);
	};
}
1;

package lmp_atom;
use Carp;
local $Carp::CarpLevel = 1; # short messages

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $data = {
		'ID'    => 0,
		'T'     => 0,
		'X'     => [0,0,0],
		'V'     => [0,0,0],
		'IMAGE' => 0,
	};
	my $self = sub {
		my $field = shift;
		croak "Not a valid field '$field' in lmp_atom"
			unless exists $data->{$field};
		$data->{$field} = shift if (@_);
		return $data->{$field};
	};
	bless $self, $class;
	return $self;
}
#generate method names
for my $field (qw(id t x v image)) {
	no strict "refs"; #for access to the symbol table
	*$field = sub {
		my $self = shift;
		return $self->(uc $field, @_);
	};
}
1;

package myself;
use Carp;
local $Carp::CarpLevel = 1; # short messages
our @ISA = "lmp_box";

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $data = {
		'PTR'        => undef,
		'IS_OPEN'    => 0,
		'BOX'        => new lmp_box(),
		'NATOMS'     => 0,
		'NMOLECULES' => 0,
		'NTYPES'     => 0,
	};
	my $self = sub : lvalue {
		my $field = shift;
		croak "Not a valid field '$field'"
			unless exists $data->{$field};
		if (@_) { $data->{$field} = shift }
		return $data->{$field};
	};
	bless($self,$class);
	return $self;
}
#generate method names
for my $field (qw(ptr is_open box natoms nmolecules ntypes)) {
	no strict "refs"; #for access to the symbol table
	*$field = sub : lvalue {
		my $self = shift;
		return $self->(uc $field, @_);
	};
}
1;

package lmp_libc;
bootstrap lmp_lib; #load shared library file
1;

package lmp_lib;
require v5.6;

use base qw(Exporter);
use base qw(DynaLoader);
use strict;
use Encode qw(decode encode);

use vars qw(@EXPORT %OWNER %ITERATORS %BLESSEDMEMBERS);

@EXPORT = qw();
our @ISA = "myself";

# ---------- BASE METHODS -------------

package lmp_lib;

my %atom_ref = (
	'mass'     => {
					'n'    => 'ntypes',
					'type' => '1',
				  },
	'id'       => {
					'n'    => 'natoms',
					'type' => '0',
				  },
	'type'     => {
					'n'    => 'ntypes',
					'type' => '0',
				  },
	'mask'     => {
					'n'    => 'natoms',
					'type' => '0',
				  },
	'image'    => {
					'n'    => 'natoms',
					'type' => '0',
				  },
	'x'        => {
					'n'    => 'natoms',
					'm'    => 3,
					'type' => '1',
				  },
	'v'        => {
					'n'    => 'natoms',
					'm'    => 3,
					'type' => '1',
				  },
	'f'        => {
					'n'    => 'natoms',
					'm'    => 3,
					'type' => '1',
				  },
	'molecule' => {
					'n'    => 'nmolecules',
					'type' => '0',
				  },
	'q'        => {
					'n'    => 'natoms',
					'type' => '1',
				  },
	'mu'       => {
					'n'    => 'natoms',
					'm'    => 4,
					'type' => '1',
				  },
	'angmom'   => {
					'n'    => 'natoms',
					'm'    => 3,
					'type' => '1',
				  },
	'torque'   => {
					'n'    => 'natoms',
					'm'    => 3,
					'type' => '1',
				  },
	'radius'   => {
					'n'    => 'natoms',
					'type' => '1',
				  },
	'rmass'    => {
					'n'    => 'natoms',
					'type' => '1',
				  },
	'ellipsoid'=> {
					'n'    => 'natoms',
					'type' => '0',
				  },
	'line'     => {
					'n'    => 'natoms',
					'type' => '0',
				  },
	'tri'      => {
					'n'    => 'natoms',
					'type' => '0',
				  },
);

sub init {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $self = new myself;
	my (@cargs, $nargs, $i, $tmpstr, @args);

	@args = @_;
	$OWNER{$self} = 1;
	bless $self, $class;

	my $__open_no_mpi = sub {
		if(@args) {
			$tmpstr = join(" ",@args);
			@cargs = split /\s+/,$tmpstr;
			unshift @cargs, "dummyvar";
		}
		$nargs = scalar(@cargs);
		$self->ptr = lmp_libc::lammps_open_no_mpi($nargs,\@cargs);
	};
	$self->$__open_no_mpi;
	$self->natoms = $self->extract_global('natoms',0);
	$self->ntypes = $self->extract_global('ntypes',0);
	$self->nmolecules = $self->extract_global('nmolecules',0);
	$self->is_open = 1;
	return $self;
}

sub DESTROY {
	return unless $_[0]->isa('HASH');
	my $self = tied(%{$_[0]});
	delete $ITERATORS{$self};
	if (exists $OWNER{$self}) {
		$self->close();
		delete $OWNER{$self};
	 }
}

# ------- FUNCTION WRAPPERS --------
sub file {
	my $self = shift;
	my $fname = shift;

	return undef if(!$self->is_open);
	lmp_libc::lammps_file($self->ptr,$fname);
	$self->natoms = $self->extract_global('natoms',0);
	$self->ntypes = $self->extract_global('ntypes',0);
	$self->nmolecules = $self->extract_global('nmolecules',0);
}

sub close {
	my $self = shift;
	
	return undef if(!$self->is_open);
	lmp_libc::lammps_close($self->ptr);
	$self->ptr = undef;
	$self->is_open = 0;
}

sub free {
	my $self = shift;
	my $ptr = shift;

	return undef if(!$self->is_open);
	lmp_libc::lammps_free($ptr);
}

sub version {
	my $self = shift;
	return undef if(!$self->is_open);
	return lmp_libc::lammps_version($self->ptr);
}

sub run {
	my $self = shift;
	return undef if(!$self->is_open);
	my $nrun = do { my $arg = shift; defined($arg) ? $arg =~ /^\d+$/ ? $arg : 0 : 0};

	$self->command("run $nrun");
}

sub command {
	my $self = shift;
	my $cmd = shift;
	return undef if(!$self->is_open);
	return if (!defined($cmd));

	$cmd = encode('UTF-8', $cmd, Encode::FB_CROAK);
	lmp_libc::lammps_command($self->ptr,$cmd);
}

sub commands_list {
	my $self = shift;
	my @clist = map {encode('UTF-8', $_, Encode::FB_CROAK) } split /\s+/,shift;
	my $nlist = scalar(@clist);
	return if($nlist == -1);
	return undef if(!$self->is_open);
	lmp_libc::lammps_commands_list($self->ptr,$nlist,\@clist);
}

sub commands_string {
	my $self = shift;
	my $cmdstr = shift;
	return undef if(!$self->is_open);
	return if (!defined($cmdstr));
	$cmdstr = encode('UTF-8', $cmdstr, Encode::FB_CROAK);
	lmp_libc::lammps_commands_string($self->ptr,$cmdstr);
}

sub get_natoms {
	my $self = shift;
	return undef if(!$self->is_open);

	return lmp_libc::lammps_get_natoms($self->ptr);
}

sub get_thermo {
	my $self = shift;
	my $tstr = shift;
	return undef if(!$self->is_open);
	return if (!defined($tstr));
	$tstr = encode('UTF-8', $tstr, Encode::FB_CROAK);
	return lmp_libc::lammps_get_thermo($self->ptr,$tstr);
}

sub set_variable {
	my $self = shift;
	return undef if(!$self->is_open);

	my $name = shift;
	return undef if(!defined($name));
	$name = encode('UTF-8', $name, Encode::FB_CROAK);

	my $val = shift;
	return undef if(!defined($val));
	$val = encode('UTF-8', $val, Encode::FB_CROAK);

	my $rval = lmp_libc::lammps_set_variable($self->ptr,$name,$val);
	die "ERROR: Cannot set variable $name $val!\n" if($rval == -1);
}


sub extract_global {
	my $self = shift;
	my $class = ref($self) || $self;
	return undef if(!$self->is_open);
	my $str = shift;
	die "ERROR: Need to specify global variable!\n" if (! defined($str));
	my $rtype = shift;
	die "ERROR: Need to specify a type (0=int,1=double)!\n" if (!defined($rtype));

	my $ret;
	my $ptr = lmp_libc::lammps_extract_global($self->ptr,$str);
	if ($rtype == 0) { #integer
		$ret = lmp_libc::c_void_p_int($ptr);
	} elsif ($rtype == 1) { # double
		$ret = lmp_libc::c_void_p_dbl($ptr);
	}
	undef $ptr;
	return $ret;
}

sub extract_atom {
	my $self = shift;
	my($class,$str,$ptr,$rptr,$ret);
	my($n,$m,$r,$t);

	$class = ref($self) || $self;
	return undef if(!$self->is_open);

	$str = shift;
	die "ERROR: Need to specify global variable!\n" if (! defined($str));

	return undef if (!exists($atom_ref{$str}));
	$n = $self->extract_global($atom_ref{$str}{n},0);
	$n++;
	$m = 1;
	if(exists($atom_ref{$str}{m})) {
		$m = $atom_ref{$str}{m};
		$m = $self->extract_global($m,0) if($m !~ /^\d+$/);
	}
	$r = $atom_ref{$str}{type};
	$t = $n * $m;

	$ptr = lmp_libc::lammps_extract_atom($self->ptr,$str);
	return undef if (! defined($ptr));

	if($m == 1) {
		$ret = lmp_libc::c_void_p_array($ptr,$t,$r);
	} else {
		$rptr = lmp_libc::c_void_p_2d_array($ptr,$n,$m,$r);
		$ret = $self->create_2d_array($rptr,$n,$m);
	}
	undef $ptr;
	undef $rptr;
	return $ret;
}

sub extract_setting {
	my $self = shift;
	return undef if(!$self->is_open);

	my $str = shift;
	die "ERROR: Need to specify global variable!\n" if (! defined($str));

	return lmp_libc::lammps_extract_settings($self->ptr,$str);

}

sub extract_box {
	my $self = shift;
	return undef if(!$self->is_open);

	my $ptr = lmp_libc::_perl_wrap_extract_box($self->ptr);
	$self->box = new lmp_box(@{ $ptr });
	undef $ptr;
}

sub reset_box {
	my $self = shift;
	return undef if(!$self->is_open);

	my $box = new lmp_box(@_);
	lmp_libc::lammps_reset_box($self->ptr,$box->lo,$box->hi,$box->xy,$box->yz,$box->xz);
	undef $box;
}

sub extract_compute {
	my $self = shift;
	return undef if(!$self->is_open);

	my $name = shift;
	die "ERROR: Need to specify compute name!\n" if (! defined($name));
	$name = encode('UTF-8', $name, Encode::FB_CROAK);

	my $type = shift;
	die "ERROR: Need to specify compute type! (0 = scalar, 1=vector, 2 = array)\n!\n" if (! defined($type));
	die "ERROR: invalid type! Expected 0 = scalar, 1=vector, 2 = array!" if($type < 0 || $type > 2);

	my $style = lmp_libc::lammps_compute_get_style($self->ptr,$name);
	return undef if($style < 0);

	my $ptr = lmp_libc::lammps_extract_compute($self->ptr,$name,$style,$type);
	return undef if (!defined($ptr));

	my ($n,$m) = lmp_libc::lammps_compute_get_dim($self->ptr,$name);
	$n++;
	$m++;
	my ($rptr, $ret);

	if($type == 0) {
		return undef if($style > 0);
		$rptr = lmp_libc::c_void_p_array($ptr,$n,1);
		undef $ptr;
		$ret = $rptr->[0];
	} elsif($type == 1) {
		$ret = lmp_libc::c_void_p_array($ptr,$n,1);
	} elsif($type == 2) {
		$rptr = lmp_libc::c_void_p_2d_array($ptr,$n,$m,1);
		undef $ptr;
		$ret = $self->create_2d_array($rptr,$n,$m);
	}
	$self->free($ptr);
	undef $rptr;
	return $ret;

}

sub extract_fix {
	my $self = shift;
	return undef if(!$self->is_open);

	my $name = shift;
	die "ERROR: Need to specify fix name!\n" if (! defined($name));
	$name = encode('UTF-8', $name, Encode::FB_CROAK);

	my $type = shift;
	die "ERROR: Need to specify fix type! (0 = scalar, 1=vector, 2 = array)\n!\n" if (! defined($type));
	die "ERROR: invalid type! Expected 0 = scalar, 1=vector, 2 = array!" if($type < 0 || $type > 2);

	my $style = lmp_libc::lammps_fix_get_style($self->ptr,$name);
	return undef if($style < 0);

	my $ptr = lmp_libc::lammps_extract_fix($self->ptr,$name,$style,$type,1,1);
	return undef if (!defined($ptr));

	my ($n,$m) = lmp_libc::lammps_fix_get_dim($self->ptr,$name);
	$n++;
	$m++;
	my ($rptr, $ret);

	if($type == 0) {
		return undef if($style > 0);
		$rptr = lmp_libc::c_void_p_array($ptr,$n,1);
		undef $ptr;
		$ret = $rptr->[0];
	} elsif($type == 1) {
		$ret = lmp_libc::c_void_p_array($ptr,$n,1);
	} elsif($type == 2) {
		$rptr = lmp_libc::c_void_p_2d_array($ptr,$n,$m,1);
		undef $ptr;
		$ret = $self->create_2d_array($rptr,$n,$m);
	}
	$self->free($ptr);
	undef $rptr;
	return $ret;
}

sub extract_variable {
	my $self = shift;
	return undef if(!$self->is_open);

	my $name = shift;
	die "ERROR: Need to specify variable name!\n" if (! defined($name));
	$name = encode('UTF-8', $name, Encode::FB_CROAK);

	my $type = shift;
	die "ERROR: Need to specify variable type! (0 = global, 1=local array)\n!\n" if (! defined($type));
	die "ERROR: invalid type! Expected 0 = global, 1=local array!" if($type < 0 || $type > 1);

	my $group = shift;

	my ($ptr, $rptr, $ret);

	if($type == 0) {
		$group = "all";
		$ptr = lmp_libc::lammps_extract_variable($self->ptr,$name,$group);
		return undef if (!defined($ptr));
		$rptr = lmp_libc::c_void_p_array($ptr,1,1);
		$ret = $rptr->[0];
	} elsif($type == 1) {
		$group = "all" if (!defined($group));
		$ptr = lmp_libc::lammps_extract_variable($self->ptr,$name,$group);
		$self->command("variable gc equal count($group)");
		my $n = $self->extract_variable("gc",0);
		$n++;
		$ret = lmp_libc::c_void_p_array($ptr,$n,1);
	}
	$self->free($ptr);
	undef $rptr;
	return $ret;
}

sub gather_atoms {
	my $self = shift;
	return undef if(!$self->is_open);

	my $name = shift;
	die "ERROR: Need to specify gather_atoms name!\n" if (! defined($name));
	$name = encode('UTF-8', $name, Encode::FB_CROAK);
	return undef if (!exists($atom_ref{$name}));

	my $type = $atom_ref{$name}{type};
	my $n = $self->natoms;
	my $m = $atom_ref{$name}{m};
	$m = 1 if (! defined($m));
	my $t = $n*$m;
	my $rptr = lmp_libc::_perl_wrap_gather_atoms($self->ptr,$name,$type,$m,$t);
	my $ret;
	if ($m == 1) {
		@{ $ret  } = @{ $rptr };
	} else  {
		$ret = $self->create_2d_array($rptr,$n,$m);
	}
	undef $rptr;
	return $ret;
}

sub scatter_atoms {
	my $self = shift;
	return undef if(!$self->is_open);

	my $name = shift;
	die "ERROR: Need to specify scatter_atoms name!\n" if (! defined($name));
	$name = encode('UTF-8', $name, Encode::FB_CROAK);
	return undef if (!exists($atom_ref{$name}));

	my $nval = shift;
	die "ERROR: Need to specify new value in scatter_atoms!\n" if (! defined($nval));

	my $type = $atom_ref{$name}{type};
	my $n = $self->natoms;
	my $m = $atom_ref{$name}{m};
	$m = 1 if (! defined($m));
	my $t = $n*$m;
	my $ref = $nval;
	$ref = $self->create_1d_array($nval) if ($m > 1);
	lmp_libc::_perl_wrap_lammps_scatter_atoms($self->ptr,$name,$type,$m,$t,$ref);
}

sub create_atom {
	my $self = shift;
	return undef if(!$self->is_open);

	my ($i, $j);
	my @atoms = @_;
	my $n = $#atoms;
	die "ERROR: Need at least 1 atom!\n"
		if($n == -1);


	for $i (0 .. $n) {
		die "ERROR: argument ($i+1) is not an atom object!\n"
			if(ref($atoms[$i]) ne "lmp_atom");
	}
	my $natoms = $self->natoms;
	my ($ids,$types,$x,$v,$images);
	for $i (@atoms) {
		$natoms++;
		push @{ $ids }, $natoms;
		push @{ $images }, $i->image;
		push @{ $types }, $i->t;
		for $j (0 .. 2) { 
			push @{ $x }, $i->x->[$j]; 
			push @{ $v }, $i->v->[$j]; 
		}
	}
	lmp_libc::lammps_create_atoms($self->ptr, $n+1, $ids, $types, $x, $v, $images, 1);
	$self->natoms += ($n+1);
	($ids,$types,$x,$v,$images) = ((),(),(),(),());
	
}

sub create_2d_array {
	my $self = shift;
	my ($arr,$n,$m)  = @_;
	my ($i,$j,$c,$ret);

	$c = 0;

	for $i (1 .. $n) {
		for $j (1 .. $m) {
			$ret->[$i][$j] = $arr->[$c];
			$c++;
		}
	}
	undef $arr;
	return $ret;
};

sub create_1d_array {
	my $self = shift;
	my $arr = shift;

	my ($i,$j,$ret);

	for $i (1 .. $#{ $arr }) {
		for $j (1 .. $#{ $arr->[$i] }) { 
			push @{ $ret }, $arr->[$i][$j];
		}
	}

	return $ret;
}
1;

# ------- FUNCTION WRAPPERS --------

package lmp_lib;

#====DONE====
#*c_void_p_int = *lmp_libc::c_void_p_int;
#*c_void_p_dbl = *lmp_libc::c_void_p_dbl;
#*c_void_p_int_array = *lmp_libc::c_void_p_int_array;
#*c_void_p_dbl_array = *lmp_libc::c_void_p_dbl_array;
#*c_void_p_2d_int_array = *lmp_libc::c_void_p_2d_int_array;
#*c_void_p_2d_dbl_array = *lmp_libc::c_void_p_2d_dbl_array;
#*lammps_open_no_mpi = *lmp_libc::lammps_open_no_mpi;
#*lammps_close = *lmp_libc::lammps_close;
#*lammps_version = *lmp_libc::lammps_version;
#*lammps_file = *lmp_libc::lammps_file;
#*lammps_command = *lmp_libc::lammps_command;
#*lammps_commands_list = *lmp_libc::lammps_commands_list;
#*lammps_commands_string = *lmp_libc::lammps_commands_string;
#*lammps_free = *lmp_libc::lammps_free;
#*lammps_get_thermo = *lmp_libc::lammps_get_thermo;
#*lammps_get_natoms = *lmp_libc::lammps_get_natoms;
#*lammps_extract_setting = *lmp_libc::lammps_extract_setting;
#*lammps_extract_global = *lmp_libc::lammps_extract_global;
#*lammps_extract_box = *lmp_libc::lammps_extract_box;
#*lammps_extract_atom = *lmp_libc::lammps_extract_atom;
#*lammps_extract_compute = *lmp_libc::lammps_extract_compute;
#*lammps_reset_box = *lmp_libc::lammps_reset_box;
#*lammps_set_variable = *lmp_libc::lammps_set_variable;
#*lammps_extract_variable = *lmp_libc::lammps_extract_variable;
#*lammps_gather_atoms = *lmp_libc::lammps_gather_atoms;
#*lammps_scatter_atoms = *lmp_libc::lammps_scatter_atoms;
#*lammps_create_atoms = *lmp_libc::lammps_create_atoms;
#*lammps_extract_fix = *lmp_libc::lammps_extract_fix;
#
#====TODO=====
#*lammps_open = *lmp_libc::lammps_open;

# ------- VARIABLE STUBS --------


package lmp_lib;

1;
