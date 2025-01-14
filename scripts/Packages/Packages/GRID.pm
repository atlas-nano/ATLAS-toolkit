package GRID;
use Carp;
local $Carp::CarpLevel = 1; # short messages
my ($cpack, $cfile) = caller();

sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my ($data);

	if (scalar(@_) == 2) { #cubic grid
		$data->{nx} = $data->{ny} = $data->{nz} = $_[0];
		$data->{sx} = $data->{sy} = $data->{sz} = $_[1];
	} elsif (scalar(@_) == 6) { #orthorhombic
		($data->{nx},$data->{ny},$data->{nz},$data->{sx},$data->{sy},$data->{sz}) = @_;
	} elsif (scalar(@_) == 4) { #orthorhombic with same sized cells
		($data->{nx},$data->{ny},$data->{nz},$data->{sx}) = @_;
		$data->{sy} = $data->{sz} = $data->{sx};
	} elsif (scalar(@_) > 0) {
		die "ERROR: Invalid number of arguments used to initialize grid!\n";
	}

	$data->{is_ready} = 0;
	$data->{grid} = ();
	my $self = sub : lvalue {
		my $field = shift;
		croak "Not a valid field '$field'"
			unless exists $data->{$field};
		if (@_) { $data->{$field} = shift }
		return $data->{$field};
	};

	bless($self,$class);
	$self->initGrid() if (scalar(@_) > 0); 
	return $self;
}
#generate method names
for my $field (qw(nx ny nz sx sy sz grid is_ready)) {
	no strict "refs"; #for access to the symbol table
	*$field = sub : lvalue {
		my $self = shift;
		return $self->(lc $field, @_);
	};
}

#######################################################
#PRIVATE SUBROUTINES
#######################################################
my $find_cell_atoms = sub : lvalue {
	my $self = shift;
	return () if (! $self->is_ready);
	return () if (scalar(@_) < 3);

	my ($ix,$iy,$iz) = @_;

	return $self->grid->{$ix}{$iy}{$iz}{ATOMS};
};


#######################################################
#PUBLIC SUBROUTINES
#######################################################

sub molOverlap {
	my $self = shift;
	my $mol = shift;
	my $offset = shift;

	return () if (! $self->is_ready);
	return if (! defined($mol));
	$offset->{XCOORD} = $offset->{YCOORD} = $offset->{ZCOORD} = 0.0 if(! defined($offset));
	my $isOverlap = 0;
	my ($i, $j);

	for $i (keys %{ $mol }) {
		$isOverlap = $self->atomOverlap($mol->{$i}, $offset);
		last if ($isOverlap);
	}

	return $isOverlap;
}

sub atomOverlap {
	my $self = shift;
	my $atom = shift;
	my $offset = shift;

	return () if (! $self->is_ready);
	return if (! defined($atom));
	$offset->{XCOORD} = $offset->{YCOORD} = $offset->{ZCOORD} = 0.0 if(! defined($offset));
	my ($i, $idx, $gs, $tmp);

	for $i ("X", "Y", "Z") {
		$tmp = "s" . lc $i; # grid dimension size
		$gs = $self->$tmp;
		$tmp = "n" . lc $i; # grid dimension max index
		$gi = $self->$tmp;
		$idx->{$i} = int(($atom->{"${i}COORD"} + $offset->{"${i}COORD"})/$gs)+1;
		while($idx->{$i} > $gi) { $idx->{$i} -= $gi; }
		while($idx->{$i} < 1)   { $idx->{$i} += $gi; }
	}
	$tmp = $self->$find_cell_atoms($idx->{X},$idx->{Y},$idx->{Z});
	my $isOverlap = 0;
	$isOverlap = 1 if (keys %{ $tmp });
	return $isOverlap;
}

sub initGrid {
	my $self = shift;

	my ($i, $j, $k, $l, $m, $n, $c, $idx, $t);
	my ($nx, $ny, $nz) = ($self->nx, $self->ny, $self->nz);
	my $grid = \%{ $self->grid };

	for $i (1 .. $nx) {
		for $j (1 .. $ny) {
			for $k (1 .. $nz) {
				$grid->{$i}{$j}{$k}{ATOMS} = ();
				$c = 1;
				for $l (-1 .. 1) {
					$idx->{X} = $l + $i;
					if($idx->{X} > $nx) { $idx->{X} -= $nx; }
					if($idx->{X} < 1)   { $idx->{X} += $nx; }
					for $m (-1 .. 1) {
						$idx->{Y} = $m + $j;
						if($idx->{Y} > $ny) { $idx->{Y} -= $ny; }
						if($idx->{Y} < 1)   { $idx->{Y} += $ny; }
						for $n (-1 .. 1) {
							next if ($l == $m and $m == $n and $n == 0);
							$idx->{Z} = $n + $k;
							if($idx->{Z} > $nz) { $idx->{Z} -= $nz; }
							if($idx->{Z} < 1)   { $idx->{Z} += $nz; }
							$grid->{$i}{$j}{$k}{NEIGH}{$c++} = $grid->{$idx->{X}}{$idx->{Y}}{$idx->{Z}};
						}
					}
				}
			}
		}
	}
	$self->is_ready = 1;
}

sub store {
	my $self = shift;
	my $atom = shift;

	return () if (! $self->is_ready);
	return () if (!defined($atom));
	my ($i, $idx, $c, $gs, $gn, $tmp, $grid);

	$grid = \%{ $self->grid }; #pointer to grid
	for $i ("X", "Y", "Z") {
		$tmp = "s" . lc $i; #grid size
		$gs = $self->$tmp;
		$tmp = "n" . lc $i; #grid max index
		$gn = $self->$tmp;
		$idx->{$i} = int($atom->{"${i}COORD"}/$gs)+1;
		while($idx->{$i} > $gn) { $idx->{$i} -= $gn };
		while($idx->{$i} < 1)   { $idx->{$i} += $gn };
	}
	$c = 1;
	while(exists($grid->{$idx->{X}}{$idx->{Y}}{$idx->{Z}}{ATOMS}{$c})) { $c++;}
	$grid->{$idx->{X}}{$idx->{Y}}{$idx->{Z}}{ATOMS}{$c} = $atom;
}

sub find_voids {
	my $self = shift;
	return () if (! $self->is_ready);

	my ($i, $j, $k, $voids, $aList, $empty, $elist, $nvoid);
	my ($nx, $ny, $nz) = ($self->nx, $self->ny, $self->nz);

	for $i (1 .. $nx) {
		for $j (1 .. $ny) {
			for $k (1 .. $nz) {
				$aList = $self->$find_cell_atoms($i,$j,$k);
				next if (keys %{ $aList });
				$empty->{$i}{$j}{$k} = 1;
			}
		}
	}
	for $i (keys %{ $empty }) {
		for $j (keys %{ $empty->{$i} }) {
			for $k (keys %{ $empty->{$i}{$j} }) {
				next if ($empty->{$i}{$j}{$k} == 0);
				$nvoid = 0;
				$nvoid = $self->find_empty_neighbors($empty, $i, $j, $k, $nvoid);
				print "";
			}
		}
	}

	return $voids;
}

sub find_empty_neighbors {
	my $self = shift;
	return () if (! $self->is_ready);
	return () if (scalar(@_) < 5);

	my ($empty, $ix, $iy, $iz, $count) = @_;
	my ($i, $j, $k, $l, $m, $n, $nlist);
	my ($nx, $ny, $nz) = ($self->nx, $self->ny, $self->nz);

	for $l (-1 .. 1) {
		$i = $l + $ix;
		if ($i > $nx) { $i -= $nx };
		if ($i < 1)   { $i += $nx };
		for $m (-1 .. 1) {
			$j = $m + $iy;
			if ($j > $ny) { $j -= $ny };
			if ($j < 1)   { $j += $ny };
			for $n (-1 .. 1) {
				$k = $n + $iz;
				if ($k > $nz) { $k -= $nz };
				if ($k < 1)   { $k += $nz };
				if(exists($empty->{$i}) and exists($empty->{$i}{$j}) and exists($empty->{$i}{$j}{$k}) and $empty->{$i}{$j}{$k} == 1) {
					$rec = { i => $i, j => $j, k => $k }; 
					push @{ $nlist }, $rec;
				}
			}
		}
	}
	return $count if (!$nlist);
	$empty->{$ix}{$iy}{$iz} = 0;
	for $i (@{ $nlist }) {
		$count++;
		$self->find_empty_neighbors($empty, $i->{i}, $i->{j}, $i->{k}, $count);
	}
}

1;
