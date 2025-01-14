package MolData;

$DB::deep = 50000; # or more if necessary
use strict;
use Carp;
use Atom;
use Molecule;
use LinkedList;
use File::Basename qw(basename); 
no warnings 'recursion';
use constant SMALL => -9999999999999;
use constant LARGE => 9999999999999;

require Exporter;

our (@ISA, @EXPORT, $VERSION, @EXPORT_OK);
local $Carp::Carplevel = 1;
my ($cpack, $cfile) = caller();
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
$VERSION = "1.00";

# Initialize the atom object
sub new {
	my $invocant = shift;
	my $class = ref($invocant) || $invocant;
	my $data = {
		atomlist  => LinkedList->new(), # this is the underlying linked list data self
		atoms	 => {}, # this is the accessor (hash) interface to the above linked list
		header	=> {},
		mollist   => {}, # underlying data self
		molecule  => {}, # accessor interface
		cell	  => { valid => 0, },
		count	 => {},
		shash	 => {}, #hash holding all field values for every atom
		vdwbox	=>	{
					"a" =>  {
							"max" => SMALL,
							"min" => LARGE,
						},
						"b" =>  {
							"max" => SMALL,
							"min" => LARGE,
						},
						"c" =>  {
							"max" => SMALL,
							"min" => LARGE,
						},
					},
	};

	my $self = sub {
		my $field = lc shift;

		# Access checks
		croak "No valid field '$field' in object"
		   unless (exists($data->{$field}) or ($field =~ /exists/));
		if ($field =~ /exists/) {
			return 0 if (! exists($data->{$field}));
			return 1;
		} elsif (@_) { 
			my $val = shift;
			return $data->{atoms}->{$val} if ($field eq "atoms");
			return $data->{mollist}->{$val} if ($field eq "molecule");
			$data->{$field} = $val;
		}
		return $data->{$field};
	};
	bless($self, $class);
	return $self;
}

for my $field (qw(atomlist header molecule cell atoms count mollist vdwbox shash)) {
	no strict "refs";
	*$field = sub {
		my $self = shift;
		return $self->(lc $field, @_);
	 };
}

#File format readers
my $readBGF = sub {
	my $self = shift;
	my $bgfFile = shift;
	my $verbose = shift;
	my ($atmPatern, $atmCounter, $molCounter, $isAtoms, $newMol); 
	my ($in_data, $boxPatern, $atomC, $type, $atom, $map); 
	my (@fields, $i, @tmp, $bondPatern, $bondAtom, $bondCounter, $bondData);

	$map = { "x" => "a", "y" => "b", "z" => "c", };
	$verbose = 1 if (! defined($verbose));
	$atmPatern = '^ATOM\s*(\d+)\s+(.{5})\s+(\S+)\s\w?\s*(\d+)\s*' .
		'(\-?\d+\.\d{5})\s*(\-?\d+\.\d{5})\s*(\-?\d+\.\d{5})\s+(\S+)\s+' .
		'(\d+)\s+(\d+)\s+(\-?\d+\.\d+)\s*(\d+\s+\d+\s+\d+\.\d+)?';
	$bondPatern = '^(CONECT|ORDER |DISP. )(.{6})(.*)';
	$boxPatern = '^(CRYSTX|CELLS)\s+(\-?\d+\.?\d*\s+.+)';

	$atmCounter = $bondCounter = $molCounter = $isAtoms = 0;
	@fields = ("index", "atmname", "resname", "resid", "x", "y",
			   "z", "fftype", "numbonds", "lonepairs", "charge");
	open BGFFILE, $bgfFile or die "ERROR: Cannot open $bgfFile: $!\n";
	while (<BGFFILE>) {
		chomp;
		$in_data = $_;
		$in_data =~ s/^HETATM/ATOM  /;
		if ($in_data =~ /$atmPatern/) {
			$isAtoms = 1;
			$atmCounter++;
			$atom = Atom->new();
			$atom->("code", $atmCounter);
			for $i (1 .. 11) {
				$atom->($fields[$i-1],eval('$' . $i));		
				$self->shash->{$fields[$i-1]}{eval('$' . $i)}{$atmCounter} = $atom;
			}
			for ("x", "y", "z") {
				$self->vdwbox->{$map->{$_}}{max} = $atom->$_ if ($atom->$_ > $self->vdwbox->{$map->{$_}}{max});
				$self->vdwbox->{$map->{$_}}{min} = $atom->$_ if ($atom->$_ < $self->vdwbox->{$map->{$_}}{min});
			}
			if ($12) {
				@tmp = split /\s+/, $12;
				$atom->("occupancy", $tmp[0]);
				$atom->("resonance", $tmp[1]);
				$atom->("radii", $tmp[2]);
				$self->shash->{"radii"}{$tmp[2]}{$atmCounter} = $atom;
			}
			# now get optional chain information
			$atom->("chain", "A");
			$atom->("chain", $1) if ($in_data =~ /^ATOM\s+\d+\s+\S+\s+\S+\s+(\w)\s+/);

			# finally get label information
			$atom->("label", "ATOM");
			$atom->("label", "HETATM") if ($_ =~ /^HETATM/);

			# add the atom mass
			$atom->("mass", 1); # placeholder for now, will have to look up element and add mass

			#insert curr atom into atoms linked list
			$self->insertAtom($atom);
			#$self->atoms->{$atmCounter} = $atom;
		} elsif ($in_data =~ /$bondPatern/) {
			($type, $atomC, $bondData) = (lc($1), $2, $3);
			$atomC =~ s/\s//g;
			$atom = $self->atoms->{$atomC};
			next if (! $atom);
			@tmp = ();
			while ($bondData =~ /(.{6})/g) {
				push @tmp, $1;
			}
			if ($type =~ /conect/) { # create a regular bond
				$i = 0;
				while ($i <= $#tmp) {
					$tmp[$i] =~ s/\s//g;
					$bondAtom = $self->atoms->{$tmp[$i]};
					if (! $bondAtom) {
						splice @tmp, $i, 1;
						next;
					}
					$bondCounter++;
					#$atom->bondlist->create($bond);
					#$bond->bondlist->create($atom);
					$atom->createBond($bondAtom, 0);
					$i++;
				}
			} elsif ($type =~ /order|disp(x|y|z)/i) { # store the cell display index/order info
				$type =~ s/\s//g;
				$bondAtom = $atom->bondlist->first;
				#$tmp[0] =~ s/\s//g;
				#$atom->bondlist->store($bond, $type, shift @tmp);
				for $i (@tmp) {
					$i =~ s/\s//g;
					$atom->bondlist->store($bondAtom, $type, $i);
					$bondAtom = $bondAtom->next;
				}
			}
		} elsif ($_ =~ /$boxPatern/i) {
			$self->cell->{"valid"} = 1;
			if (uc($1) eq  "CRYSTX") {
				print "Using Box Information from file..." if ($verbose);
				($self->cell->{a}, $self->cell->{b}, $self->cell->{c},
				 $self->cell->{alpha}, $self->cell->{beta}, $self->cell->{gamma}) = split /\s+/, $2;
			} else {
				($self->cell->{image}{a}{p}, $self->cell->{image}{a}{n}, $self->cell->{image}{b}{p},
				$self->cell->{image}{b}{n}, $self->cell->{image}{c}{p}, $self->cell->{image}{c}{n}) = split /\s+/, $2;
			}
		} elsif (! $isAtoms and $in_data =~ /^(\S+)\s(.+)/) {
			$self->header->{lc($1)} .= "\n" if ($self->header->{lc($1)});
			$self->header->{lc($1)} .= $2;
		}
	}
	close BGFFILE;
	die "ERROR: No valid data found in $bgfFile: $!\n" if (! $atmCounter);

	$self->count->{"atoms"} = $atmCounter;
	$self->count->{"bonds"} = $bondCounter;
	$self->header->{"footer"} = "END\n";
};

my $readPDB = sub {
	my $self = shift;
	my $pdbFile = shift;
	my $verbose = shift;
	my ($atmPatern, $atmCounter, $molCounter, $isAtoms, $newMol); 
	my ($in_data, $boxPatern, $molid, $atomC, $type, $atom, $map); 
	my (@fields, $i, @tmp, $bondPatern, $bondAtom, $bondData, $bondCounter, $fftype, $q);

	$map = { "x" => "a", "y" => "b", "z" => "c", };
	$verbose = 1 if (! defined($verbose));
	$atmPatern = '^ATOM\s*(\d+)\s+(.{4})\s+(\S+)\s\w?\s*(\d+)\s*' .
		'(\-?\d+\.\d{3})\s*(\-?\d+\.\d{3})\s*(\-?\d+\.\d{3})\s*(\-?\d+\.\d{3})';
	$bondPatern = '^(CONECT)(.{5})(.+)';
	$boxPatern = '^(CRYST1)\s+(\-?\d+\.?\d*\s+.+)';

	$atmCounter = $bondCounter = $molCounter = $isAtoms = 0;
	@fields = ("index", "atmname", "resname", "resid", "x", "y", "z");
	open PDBFILE, $pdbFile or die "ERROR: Cannot open $pdbFile: $!\n";
	while (<PDBFILE>) {
		chomp;
		$in_data = $_;
		$in_data =~ s/^HETATM/ATOM  /;
		if ($in_data =~ /$atmPatern/) {
			$q = 0;
			$q = $8 if (defined($8));
			$fftype = $2;
			$isAtoms = 1;
			$atmCounter++;
			$atom = Atom->new();
			$atom->("code", $atmCounter);
			for $i (1 .. 7) {
				$atom->($fields[$i-1],eval('$' . $i));		
				$self->shash->{$fields[$i-1]}{eval('$' . $i)}{$atmCounter} = $atom;
			}
			for ("x", "y", "z") {
				$self->vdwbox->{$map->{$_}}{max} = $atom->$_ if ($atom->$_ > $self->vdwbox->{$map->{$_}}{max});
				$self->vdwbox->{$map->{$_}}{min} = $atom->$_ if ($atom->$_ < $self->vdwbox->{$map->{$_}}{min});
			}
			# now get optional chain information
			$atom->("chain", "A");
			#$atom->("chain", $1) if ($in_data =~ /^ATOM\s+\d+\s+\S+\s+\S+\s+(\w)\s+/);

			# finally get label information
			$atom->("label", "ATOM");
			$atom->("label", "HETATM") if ($_ =~ /^HETATM/);

			# add the atom mass
			$atom->("mass", 1); # placeholder for now, will have to look up element and add mass

			#set the fftype to be the atmname
			$fftype =~ s/\s//g;
			$atom->("fftype", $fftype);

			#lonepairs
			$atom->("lonepairs", 0);

			#charge
			$atom->("charge", $q);

			#insert curr atom into atoms linked list
			$self->insertAtom($atom);
			#$self->atoms->{$atmCounter} = $atom;
		} elsif ($in_data =~ /$bondPatern/) {
			($type, $atomC, $bondData) = (lc($1), $2, $3);
			$atomC =~ s/\s//g;
			$atom = $self->atoms->{$atomC};
			next if (! $atom);
			@tmp = ();
			while ($bondData =~ /(.{6})/g) {
				push @tmp, $1+0;
			}
			$i = 0;
			undef($molid);
			$molid = $atom->molid if (defined($atom->molid));
			while ($i <= $#tmp) {
				$tmp[$i] =~ s/\s//g;
				$bondAtom = $self->atoms->{$tmp[$i]};
				if (! $bondAtom) {
					splice @tmp, $i, 1;
					next;
				}
				$molid = $bondAtom->molid if (! defined($molid) and defined($bondAtom->molid));
				$bondCounter++;
				$atom->createBond($bondAtom);
				#$atom->bondlist->create($bond);
				#$atom->("numbonds",$atom->bondlist->count);
				#$bond->bondlist->create($atom);
				#$bond->("numbonds",$bond->bondlist->count);
				$i++;
			}
			
		} elsif ($_ =~ /$boxPatern/i) {
			$self->cell->{"valid"} = 1;
			print "Using Box Information from file..." if ($verbose);
			($self->cell->{a}, $self->cell->{b}, $self->cell->{c},
			$self->cell->{alpha}, $self->cell->{beta}, $self->cell->{gamma}) = split /\s+/, $2;
		} elsif (! $isAtoms and $in_data =~ /^(\S+)\s(.+)/) {
			$self->header->{lc($1)} .= "\n" if ($self->header->{lc($1)});
			$self->header->{lc($1)} .= $2;
		}
	}
	close PDBFILE;
	die "ERROR: No valid data found in $pdbFile: $!\n" if (! $atmCounter);

	$self->count->{"atoms"} = $atmCounter;
	$self->count->{"bonds"} = $bondCounter;
	$self->header->{"footer"} = "END\n";
};

my $readMOL2 = sub {
	my $self = shift;
	my $mol2File = shift;
	my $verbose = shift;
	my ($atmPatern, $atmCounter, $molCounter, $isAtoms, $newMol, $map); 
	my ($in_data, $molid, $atomC, $order, $atom, $isBonds, $hC, $bondC); 
	my (@fields, $bondPatern, $bondAtom, $bondData, $bondCounter, $isHeader, $i, $isFooter);

	$map = { "x" => "a", "y" => "b", "z" => "c", };
	$verbose = 1 if (! defined($verbose));
	$atmPatern = '^\s+(\d+)\s+(\S+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)' .
		'\s+(\-?\d+\.\d+)\s+(\S+)\s+(\d+)\s+(\S*)\s*(\-?\d+\.\d+)';
	$bondPatern = '^\s+\d+\s+(\d+)\s+(\d+)\s+(\w+)';

	$atmCounter = $bondCounter = $molCounter = $hC = 0; 
	$isBonds = $isAtoms = $isHeader = $isFooter = 0;
	@fields = ("index", "atmname", "x", "y", "z", 
			   "fftype","resid", "resname", "charge");
	open MOL2FILE, $mol2File or die "ERROR: Cannot open $mol2File: $!\n";
	while (<MOL2FILE>) {
		chomp;
		$in_data = $_;
		if ($_ =~ /@<TRIPOS>MOLECULE/) {
			$isHeader = 1;
			$isAtoms = $isBonds = $isFooter = 0;
		} elsif ($_ =~ /^@<TRIPOS>ATOM/) {
			$isAtoms = 1;
			$isBonds = $isHeader = $isFooter = 0;
		} elsif ($_ =~ /^@<TRIPOS>BOND/) {
			$isBonds = 1;
			$isAtoms = $isHeader = $isFooter = 0;
		} elsif ($_ =~ /^@<TRIPOS>SUBself/) {
			$isFooter = 1;
			 $isBonds = $isAtoms = $isHeader = 0;
			$self->header->{"footer"} = "$_\n";
		} elsif ($in_data =~ /$atmPatern/ and $isAtoms) {
			$atmCounter++;
			$atom = Atom->new();
			$atom->("code", $atmCounter);
			for $i (1 .. 9) {
				$atom->($fields[$i-1],eval('$' . $i));
				$self->shash->{$fields[$i-1]}{eval('$' . $i)}{$atmCounter} = $atom;
			}
			#insert curr atom into atoms linked list
			$atom->("numbonds",0);
			$atom->("chain","X");
			$atom->("label", "ATOM");
			$atom->("lonepairs", 0);
			$atom->("resname", substr($atom->resname,0,3)) if (length($atom->resname) > 3);
			$self->insertAtom($atom);
			#$self->atoms->{$atmCounter} = $atom;
			for ("x", "y", "z") {
				$self->vdwbox->{$map->{$_}}{max} = $atom->$_ if ($atom->$_ > $self->vdwbox->{$map->{$_}}{max});
				$self->vdwbox->{$map->{$_}}{min} = $atom->$_ if ($atom->$_ < $self->vdwbox->{$map->{$_}}{min});
			}

		} elsif ($in_data =~ /$bondPatern/ and $isBonds) {
			($atomC, $bondC, $order) = ($1, $2, $3);
			$atom = $self->atoms->{$atomC};
			next if (! $atom);
			# create a regular bond
			undef($molid);
			$bondAtom = $self->atoms->{$bondC};
			next if (! $bondAtom);
			$bondCounter++;
			$atom->("numbonds",$atom->numbonds +1);
			$bondAtom->("numbonds",$bondAtom->numbonds +1);
			$atom->createBond($bondAtom);
			#$atom->bondlist->create($bond);
			#$bond->bondlist->create($atom);
			#now store the type information
			$atom->bondlist->store($atom->bondlist->last, "order", $order);
			$bondAtom->bondlist->store($bondAtom->bondlist->last, "order", $order);
		} elsif ($isHeader) {
			$hC++;
			if ($hC == 1) {
				$self->header->{"descrp"} .= $_;
			} elsif ($hC == 2) {
				$self->header->{"mol2count"} .= $_;
			} else {
				$self->header->{"remark"} .= "$_\n" if ($_ =~ /\S+/);
			}
		} elsif ($isFooter) {
			$self->header->{"footer"} .= "$_\n";
		}
	}
	close MOL2FILE;
	die "ERROR: No valid data found in $mol2File: $!\n" if (! $atmCounter);
	
	chomp($self->header->{"remark"});
	$self->count->{"atoms"} = $atmCounter;
	$self->count->{"bonds"} = $bondCounter;
	$self->count->{"molecule"} = $molCounter;
};

my $readXYZ = sub {
	my $self = shift;
	my $xyzFile = shift;
	my $verbose = shift;
	my ($atmPatern, $atmCounter, $molCounter, $isAtoms, $newMol, $totAtoms); 
	my ($in_data, $molid, $atomC, $type, $atom, $map, $lineCounter); 
	my (@fields, $i, @tmp, $bondPatern, $bondAtom, $bondData, $bondCounter, $fftype, $q);

	$map = { "x" => "a", "y" => "b", "z" => "c", };
	$verbose = 1 if (! defined($verbose));
	$atmPatern = '^\s*(\S+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)';
	$bondPatern = '^(CONECT)(.{5})(.+)';

	$atmCounter = $bondCounter = $molCounter = $isAtoms = $lineCounter = $totAtoms = 0;
	@fields = ("atmname", "x", "y", "z");
	open XYZFILE, $xyzFile or die "ERROR: Cannot open $xyzFile: $!\n";
	while (<XYZFILE>) {
		chomp;
		$in_data = $_;
		$lineCounter++;
		if($lineCounter == 1) {
			#first line has the number of atoms
			if ($_ =~ /^\s*(\d+)/) {
				$totAtoms = $1;
			}
		} elsif ( $lineCounter == 2) {
			$self->header->{"remark"} = "$_\n";
		} elsif ($in_data =~ /$atmPatern/) {
			$fftype = $1; #fftype is atomname
			$isAtoms = 1;
			$atmCounter++;
			$atom = Atom->new();
			$atom->("code", $atmCounter);
			for $i (1 .. 4) {
				$atom->($fields[$i-1],eval('$' . $i));		
				$self->shash->{$fields[$i-1]}{eval('$' . $i)}{$atmCounter} = $atom;
			}
			for ("x", "y", "z") {
				$self->vdwbox->{$map->{$_}}{max} = $atom->$_ if ($atom->$_ > $self->vdwbox->{$map->{$_}}{max});
				$self->vdwbox->{$map->{$_}}{min} = $atom->$_ if ($atom->$_ < $self->vdwbox->{$map->{$_}}{min});
			}
			# now set missing fields
			#chain
			$atom->("chain", "A");

			#label
			$atom->("label", "ATOM");

			#residue
			$atom->("resname", "RES");

			#resnum
			$atom->("resnum", 444);
			
			#atom mass
			$atom->("mass", 1); # placeholder for now, will have to look up element and add mass

			#set the fftype to be the atmname
			$fftype =~ s/\s//g;
			$atom->("fftype", $fftype);

			#lonepairs
			$atom->("lonepairs", 0);

			#charge
			$atom->("charge", 0);

			#insert curr atom into atoms linked list
			$self->insertAtom($atom);
			#$self->atoms->{$atmCounter} = $atom;
		} elsif ($in_data =~ /$bondPatern/) {
			($type, $atomC, $bondData) = (lc($1), $2, $3);
			$atomC =~ s/\s//g;
			$atom = $self->atoms->{$atomC};
			next if (! $atom);
			@tmp = ();
			while ($bondData =~ /(.{6})/g) {
				push @tmp, $1+0;
			}
			$i = 0;
			undef($molid);
			$molid = $atom->molid if (defined($atom->molid));
			while ($i <= $#tmp) {
				$tmp[$i] =~ s/\s//g;
				$bondAtom = $self->atoms->{$tmp[$i]};
				if (! $bondAtom) {
					splice @tmp, $i, 1;
					next;
				}
				$molid = $bondAtom->molid if (! defined($molid) and defined($bondAtom->molid));
				$bondCounter++;
				$atom->createBond($bondAtom);
				#$atom->bondlist->create($bond);
				#$atom->("numbonds",$atom->bondlist->count);
				#$bond->bondlist->create($atom);
				#$bond->("numbonds",$bond->bondlist->count);
				$i++;
			}			
		} 
	}
	close XYZFILE;
	die "ERROR: No valid data found in $xyzFile: $!\n" if (! $atmCounter);
	die "ERROR: Incorrect XYZ file format! (#atoms in line 1: '$totAtoms' != #atoms read: '$atmCounter'\n";

	$self->count->{"atoms"} = $atmCounter;
	$self->count->{"bonds"} = $bondCounter;
	$self->header->{"footer"} = "END\n";
};

my $readCharmmPSF = sub {
	die "Not implemented yet...\n";
};

my $readAmberPRMTOP = sub {
	die "Not implemented yet...\n";
};

my $readGromacsGro = sub {
	my $self = shift;
	my $groFile = shift;
	my $verbose = shift;

	my ($map, $atmPatern, $atmCounter, $isAtoms, $totAtoms);
	my (@fields, $in_data, $fftype, $atom, $i);

	$map = { "x" => "a", "y" => "b", "z" => "c", };
	$verbose = 1 if (! defined($verbose));
	$atmPatern = '^\s*(\d+)(\S+)\s+(\S+)\s+(\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)';

	$atmCounter = 0;
	@fields = ("resid", "resname", "atmname", "index", "x", "y", "z");
	open GROFILE, $groFile or die "ERROR: Cannot open $groFile: $!\n";
	while (<GROFILE>) {
		chomp;
		$in_data = $_;
		if ($in_data =~ /$atmPatern/) {
			$fftype = $3; #fftype is atmname - placeholder for now
			$atmCounter++;
			$atom = Atom->new();
			$atom->("code", $atmCounter);
			for $i (1 .. 7) {
				$atom->($fields[$i-1],eval('$' . $i));		
				$self->shash->{$fields[$i-1]}{eval('$' . $i)}{$atmCounter} = $atom;
			}
			for ("x", "y", "z") {
				$atom->($_, $atom->$_ * 10); #convert from nm to Angstroms
				$self->vdwbox->{$map->{$_}}{max} = $atom->$_ if ($atom->$_ > $self->vdwbox->{$map->{$_}}{max});
				$self->vdwbox->{$map->{$_}}{min} = $atom->$_ if ($atom->$_ < $self->vdwbox->{$map->{$_}}{min});
			}
			# now set missing fields
			#chain
			$atom->("chain", "A");

			#label
			$atom->("label", "ATOM");

			#atom mass
			$atom->("mass", 1); # placeholder for now, will have to look up element and add mass

			#set the fftype to be the atmname
			$fftype =~ s/\s//g;
			$atom->("fftype", $fftype);

			#lonepairs
			$atom->("lonepairs", 0);

			#charge
			$atom->("charge", 0);

			#insert curr atom into atoms linked list
			$self->insertAtom($atom);
			#$self->atoms->{$atmCounter} = $atom;
		} 
	}
	close XYZFILE;
	die "ERROR: No valid data found in $groFile: $!\n" if (! $atmCounter);

	$self->count->{"atoms"} = $atmCounter;
};

my $readGromacsTop = sub {
	my $self = shift;
	my $topFile = shift;
	my $verbose = shift;

	my ($groFile, $in_data, $mode, $i, $j);
	my ($atomi, $atomj, $indxFromName, $bondCounter);

	$verbose = 1 if (! defined($verbose));
	$bondCounter = $mode = 0;

	open GTOPFILE, $topFile or die "ERROR: Cannot open $topFile: $!\n";
	while (<GTOPFILE>) {
		chomp;
		$in_data = $_;
		if ($in_data =~ /^\{\s+atomType\s+mass/i) { #masses
			$mode = 1; $i = 1;
		} elsif ($in_data =~ /^\{\s+atomName\s+atomType\s+charge/i) { #atom types and charge
			$mode = 2; $i = 1;
		} elsif ($in_data =~ /^\{\s+Bonds:/i) { #bonds
			$mode = 3;
		} elsif ($mode==1 and $in_data =~ /^MASS\s+(\S+)\s+(\d+\.\d+)/) {
			$atomi = $self->atoms->{$i};
			$atomi->("atmname", $1);
			$atomi->("mass", $2);
			$indxFromName->{$1} = $i; $i++;
		} elsif ($mode==2 and $in_data =~ /^ATOM\s+(\S+)\s+TYPE=\s+(\S+)\s+CHARGE=\s+(\-?\d+\.\d+)/) {
			$atomi = $self->atoms->{$i};
			$atomi->("atmname", $1);
			$atomi->("fftype", $2);
			$atomi->("charge", $3);
			$indxFromName->{$1} = $i; $i++;
		} elsif($mode==3 and $in_data =~ /^BOND\s+(\S+)\s+(\S+)/) {
			$i = $indxFromName->{$1}; $j = $indxFromName->{$2};
			$atomi = $self->atoms->{$i};
			$atomj = $self->atoms->{$j};
			$bondCounter++;
			$atomi->createBond($atomj);
			#$atomi->bondlist->create($atomj);
			#$atomj->bondlist->create($atomi);
			#$atomi->("numbonds",$atomi->bondlist->count);
			#$atomj->("numbonds",$atomj->bondlist->count);
		}
	}
	close GTOPFILE;
	$self->count->{"bonds"} = $bondCounter;
	$self->header->{"footer"} = "END\n";	
};

my $readGromacs = sub {
	my $self = shift;
	my $groFile = shift;
	my $verbose = shift;
	my ($topFile) = shift;

	if (! defined($topFile)) {
		$topFile = $groFile;
		$topFile =~ s/\.\w+$/\.top/;
		print "Searcing for $groFile..." if ($verbose);
		if (! $self->testFile($groFile,0, $verbose)) {
			$topFile = $groFile;
			$topFile =~ s/\.w+$/\.xplor/;
			$self->testFile($groFile);
		}
	} else {
		$self->testFile($groFile);
	}
	print "Reading topology info from $groFile..." if($verbose);
	$self->$readGromacsGro($groFile,$verbose);
	print "Reading coordinate infor from $topFile..." if ($verbose);
	$self->$readGromacsTop($topFile, $verbose);
};

#some helper functions
my $createHeader = sub {
	my $self = shift;
	my $type = shift;
	my $npairs = shift;

	my ($outData, $i, @tmp);

	if ($type eq "bgf") {
		for $i ("biogrf", "descrp", "remark", "forcefield", "period", "axes", "sgname") {
			next if (! $self->header->{$i});
			@tmp = split /\n/, $self->header->{$i};
			for (@tmp) {
				$outData .= uc($i) . " $_\n";
			}
		}
		if ($self->cell->{valid}) {
			$outData .= sprintf("CRYSTX %11.5f%11.5f%11.5f%11.5f%11.5f%11.5f\n", 
						$self->cell->{a}, $self->cell->{b}, $self->cell->{c}, 
						$self->cell->{alpha}, $self->cell->{beta}, $self->cell->{gamma});
			$outData .= sprintf("CELLS  %5d%5d%5d%5d%5d%5d\n", $self->cell->{image}{a}{p},
						$self->cell->{image}{a}{n}, $self->cell->{image}{b}{p}, $self->cell->{image}{b}{n},
						$self->cell->{image}{c}{p}, $self->cell->{image}{c}{n}) if(defined($self->cell->{image}));
		}
		$outData .= "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,f10.5)\n";
	} elsif ($type eq "mol2") {
		$outData = '@<TRIPOS>MOLECULE' . "\n";
		$self->header->{mol2count} = sprintf("%-5d%-5d%-5d", $self->count->{atoms}, $npairs, $self->count->{molecule});
		for $i ("descrp", "mol2count", "remark") {
			next if (! $self->header->{$i});
			@tmp = split /\n/, $self->header->{$i};
			for (@tmp) {
				$outData .= "$_\n";
			}
		}
		$outData .= "\n\@<TRIPOS>ATOM\n";
	} elsif ($type eq "pdb") {
		$outData = "HEADER	PDB file\n";
		$outData .= "TITLE	 " . $self->header->{"descrp"} if ($self->header->{"descrp"});
	} elsif ($type eq "xyz") {
		$outData = $self->count->{atoms} . "\nCreated on " . localtime() . "\n";
	}

	return $outData;
};

my $getAtomFileFormat = sub {
	my $self = shift;
	my $type = shift;
	my (@fields, @fmt);

	if ($type eq "bgf") {
		@fmt = ('%-6s ','%5d ','%-5s ','%4s ','%1s',' %5d','%10.5f','%10.5f','%10.5f ','%-5s','%3d','%2d ','%9.6f');
		@fields = ("label", "index", "atmname", "resname", "chain", "resid", "x", "y", "z", "fftype", "numbonds", "lonepairs", "charge");
		#$self->atoms(1)->("chain","") if (!$self->atoms(1)->chain);
	} elsif ($type eq "mol2") {
		@fmt = ('%7d ','%-4s ','%13.4f','%10.4f','%10.4f ','%-4s','%7d','%4s','%16.5f');
		@fields = ("index", "atmname", "x", "y", "z", "fftype", "resid", "resname", "charge");
	} elsif ($type eq "pdb") {
		@fmt = ('%-6s', '%5d ', '%4s', '%3s ', '%1s', '%4d	', '%8.3f', '%8.3f', '%8.3f');
		@fields = ("label", "index", "atmname", "resname", "chain", "resid", "x", "y", "z");
	} elsif ($type eq "xyz") {
		@fmt = ('%-6s', '%8.3f', '%8.3f', '%8.3f');
		@fields = ("label", "x", "y", "z");
	}

	return (\@fields, \@fmt);
};

my $writeBondsData = sub {
	my $self = shift;
	my $type = shift;
	my ($OUTDATA, $count, $currAtom, $bondAtom, $bondCount, %outData, $i, $bondData, $tot);

	$tot = 0;
	$count = 1;
	$currAtom = $self->atomlist->first;
	$OUTDATA = "";
	if ($type eq "bgf") {
		$OUTDATA  = "FORMAT CONECT (a6,12i6)\n";
		while ($count <= $self->count->{"atoms"}) {
			%outData = ();
			$outData{conect} = "";
			$bondCount = 1;
			$bondAtom = $currAtom->data->bondlist->first;
			while ($bondCount <= $currAtom->data->bondlist->count) {
				next if (! $bondAtom);
				$outData{conect} .= sprintf("%6d", $bondAtom->data->index);
				for $i ("order", "dispx", "dispy", "dispz") {
					if (defined($bondAtom->$i)) {
						$bondData = $bondAtom->$i;
						if ($bondData =~ /ar|am|[a-z]+/) {
							$bondData = 2;
						}
						$outData{$i} .= sprintf("%6d", $bondData);
					}
				}
				$bondCount++;
				$bondAtom = $bondAtom->next;
			}
			for $i ("conect", "order", "dispx", "dispy", "dispz") {
				next if (! exists($outData{$i}) or ($i ne "conect" and $outData{$i} !~ /[1-9]/));
				$outData{$i} =~ s/(\s+0)*$//;
				$OUTDATA .= sprintf("%-6s%6d$outData{$i}\n",uc($i), $currAtom->data->index);
			}
			$count++;
			$currAtom = $currAtom->next;
		}
	} elsif ($type eq "pdb") {
		while ($count <= $self->count->{"atoms"}) {
			%outData = ();
			$bondCount = 1;
			$bondAtom = $currAtom->data->bondlist->first;
			while ($bondCount <= $currAtom->data->bondlist->count) {
				next if (! $bondAtom);
				$outData{conect} .= sprintf("%5d", $bondAtom->data->index);
				$bondCount++;		
				$bondAtom = $bondAtom->next;
			}
			$OUTDATA .= sprintf("%-6s%5d$outData{conect}\n","CONECT", $currAtom->data->index) if(keys %outData);
			$count++;
			$currAtom = $currAtom->next;
		}
	} elsif ($type eq "mol2") {
		$OUTDATA .= '@<TRIPOS>BOND' . "\n";
		while ($count <= $self->count->{"atoms"}) {
			$bondCount = 1;
			$bondAtom = $currAtom->data->bondlist->first;
			while ($bondCount <= $currAtom->data->bondlist->count) {
				next if (! $bondAtom);
				if ($currAtom->data->index > $bondAtom->data->index) {
					$bondCount++;		
					$bondAtom = $bondAtom->next;
					next;
				}
					
				$tot++;
				$bondData = 1;
				$bondData = $bondAtom->order if ($bondAtom->order);
				$OUTDATA .= sprintf("%6d %4d %4d %-4s\n", $tot, $currAtom->data->index, $bondAtom->data->index, $bondData);
				$bondCount++;		
				$bondAtom = $bondAtom->next;
			}
			$count++;
			$currAtom = $currAtom->next;
		}
	} else {
		printf "";
	}
	return ($OUTDATA, $tot);
};

my $getContainer = sub {
	my $self = shift;
	my $atom = shift;
	my ($i, $curr);
	return if (! $atom);
	$curr = $self->atomlist->first;
	return $curr if ($curr->data == $atom);
	$curr = $curr->next;
	while ($curr != $self->atomlist->first) {
		return $curr if ($curr->data == $atom);
		$curr = $curr->next;
	}
	return;
};

my $findPointCoord = sub {
	my $self = shift;
	my $point = shift;
	my $com = shift;
	my ($ptCoord, $map, $dim, $box);
	
	$map = ({ "a" => "x", "b" => "y", "c" => "z" });
	$point = "origin" if (! defined($point));
	if (! ref($point) or ref($point) eq "SCALAR") { # if it's a string
		if ($point =~ /^\s*(origin|box.?origin|cell.?origin)/) { #move to the origin
			$ptCoord = ({ "x" => 0, "y" => 0, "z" => 0 });
		} elsif ($point =~ /^\s*(center|box.?center|cell.?center)/) {
			if ($self->cell->{valid}) {
				for $dim ("a", "b", "c") {
					$ptCoord->{ $map->{$dim} } = ($self->cell->{$dim}/2);
				}
			} else {
				for $dim ("a", "b", "c") {
					$ptCoord->{ $map->{$dim} } = ($self->box->{$dim}->{max} - $self->box->{$dim}->{min})/2;
				}
			}
		} elsif ($point =~ /^\s*(com|com.?center)/) { # place the center of mass of the molecules in the box/cell center
			if ($self->cell->{valid}) {
				for $dim ("a", "b", "c") {
					$ptCoord->{ $map->{$dim} } = ($self->cell->{$dim}/2) - $com->{ $map->{$dim} };
				}
			} else {
				for $dim ("a", "b", "c") {
					$ptCoord->{ $map->{$dim} } = (($self->box->{$dim}->{max} - $self->box->{$dim}->{min})/2) - $com->{$dim};
				}
			}
		}
	} elsif (ref($point) eq "HASH") {
		$ptCoord = $point;
	} else {
		return;
	}
	return $ptCoord;
};

sub find {
# invoked as self->find("field", "value", "operator") where operator is:
# <,=,>,!=,<=,>= for integers
# eq,neq,gt,lt for strings
	my $self = shift;
	my $field = lc shift;
	my $val = shift;
	my ($hashResults, @tmp, $isNum, $i, $sResults, $operator, $sStr);

	return () if (! defined($field) or ! defined($val));
	return () if (! exists($self->shash->{$field}));

	$operator = shift if (@_);

	$field =~ s/\s//g;
	$val =~ s/\s//g;

	$hashResults = $self->shash->{$field};
	@tmp = keys %{ $hashResults };
	$isNum = 1;
	for $i (@tmp) {
		if ($i =~ /[A-Za-z]/) {
			$isNum = 0;
			last;
		}
	}

	if ($isNum) {
		$operator = "==" if (! $operator);
		return () if ('< = == > != <= >=' !~ /$operator/);
		return () if ($val !~ /^\-?\d+\.?\d*e?\-?\d*/);
	} else {
		$operator = "eq" if (! $operator);
		return () if ('q neq gt lt' !~ /$operator/i);
		return () if ($val !~ /[A-Za-z]/);
	}

	$operator = "==" if ($operator eq "=");
	$operator = lc $operator;
	for $i (@tmp) {
		$sStr = "$i $operator $val" if ($isNum);
		$sStr = "\"$i\" $operator \"$val\"" if (! $isNum);
		if (eval($sStr)) {
			$sResults->{$_} = $hashResults->{$i}{$_} for (keys %{ $hashResults->{$i} });
		}
	}

	return $sResults;
}

sub extractAtoms {
#will return a hash of atoms as well as their bonds only to atoms also in the list
	my $self = shift;
	my $atomsList = shift;
	my ($returnList, $i, $atom, $bond, $j);

	return () if (! $atomsList or ! keys %{ $atomsList });

	for $i (keys %{ $atomsList }) {
		$atom = $atomsList->{$i};
		$bond = $atom->bondlist->first;
		$j = 1;
		while ($j <= $atom->bondlist->count) {
			if (! exists($atomsList->{$bond->data->index})) {
				$bond = $bond->next;
				$atom->bondlist->remove($bond->prev->data);
			} else {
				$bond = $bond->next;
				$j++;
			}
		}
		$returnList->{$i} = $atom;
	}

	return $returnList;
}

sub deleteRes {
	my $self = shift;
	my $resID = shift;
	my $storeBonds = shift;
	my ($resAtoms, $bondAtm, $bondAtmResID); 
	my ($atom, $count, $rec, $store, $i);

	return if (! defined($resID));
	$resAtoms = $self->shash->{"resid"}{$resID};
	return if (! $resAtoms );
	$store = 1;
	$store = 0 if (! defined($storeBonds));
	$count = 0;

	for $atom (values %{ $resAtoms }) {
		if ($atom->bondlist->count > 0) {
			$bondAtm = $atom->bondlist->first;
			for $i (1 .. $atom->bondlist->count) {
				$bondAtmResID = $bondAtm->data->resid;
				if ($store and $bondAtmResID != $resID) {
					$count++;
					$rec = { 
							atom => $atom->atmname, 
							bond => $bondAtm->data 
						   };
					$storeBonds->{$count} = $rec;
				}
				$bondAtm = $bondAtm->next;
			}
		}
		$self->deleteAtom($atom);
	}
	delete $self->shash->{"resid"}{$resID};
}

sub insertRes {
	return () if (scalar(@_) < 4);
	my $self = shift;
	my $resID = shift;
	my $newAtmList = shift;
	my $bondOpts = shift;
	my ($atom, $i, $insertPoint, $nHash, $atmName, @tmp, $bondAtom, $atomCount);
	
	# first find the insertion point = first atom in res+1 or set to first atom 
	# (so that will insert as the last atom: first->prev = last)
	$insertPoint = $self->atomlist->first;
	if (exists($self->shash->{"resid"}{($resID + 1)})) {
		@tmp = keys %{ $self->shash->{"resid"}{($resID + 1)} };
		@tmp = sort { $a<=>$b } @tmp; #sort it numerically
		$insertPoint = $self->$getContainer($self->atoms(shift @tmp));
	}

	$atomCount = $self->count->{"atoms"};
	@tmp = sort { $a<=>$b} keys %{ $newAtmList };
	for $i (@tmp) {
		$atomCount++;
		$atom = $newAtmList->{$i};
		$atom->("resid", $resID);
		$atom->("code", $atomCount);
		$self->insertAtom($atom, $insertPoint);
		$atmName = $atom->atmname;
		$nHash->{ $atmName } = $atom;
	}

	for $i (values %{ $bondOpts }) {
		$atom = $nHash->{ $i->{atom} };
		$bondAtom = $i->{bond};
		$atom->createBond($bondAtom);
		#$atom->bondlist->create($bond);
		#$bond->bondlist->create($atom);
	}

	$self->count->{"atoms"} = $atomCount;
}

sub findContainer {
	my $self = shift;
	my $atom = shift;
	
	return () if (! $atom);
	return $self->$getContainer($atom);
}

sub createAtom {
	my $self = shift;
	my ($newAtom);

	$newAtom = Atom->new();

	return $newAtom;
}

sub cloneAtom {
	my $self = shift;
	my $currAtom = shift;
	my ($newAtom, $fields, $i, $val);
	return () if (! $currAtom);

   $newAtom = Atom->new();
   $fields = $currAtom->fields;

   for $i (split(/\s/, $fields)) {
		next if ($i =~ /bond|fields/);
		$val = $currAtom->$i;
		$newAtom->($i, $val);
   }

   $currAtom->bondlist->clone($newAtom->bondlist);
   $newAtom->bondlist->("count", $currAtom->bondlist->count);
   return $newAtom;
}

sub insertAtom {
	my $self = shift;
	my $atom = shift;
	my $insertPoint = shift;
	my ($fields, $i);
	return () if (! $atom);

	$self->atomlist->insert($atom) if (! defined($insertPoint));
	$self->atomlist->insert($atom, $insertPoint) if (defined($insertPoint));
	
	@{ $fields } = split /\s/, $atom->fields;
	for $i (@{ $fields }) {
		next if ($i =~ /bond/ or ! exists($self->shash->{$i}));
		$self->shash->{$i}{ $atom->$i }{ $atom->code } = $atom;
	}
	if (! defined($insertPoint)) {
		$self->count->{atoms}++;
		$self->atoms->{ $self->count->{atoms} } = $atom;
	}
}

sub read {
	my $self = shift;
	my $nameStr = shift;
	my $type = lc shift;

	my $opts;
	@{ $opts } = split /\s+/,$nameStr;
	my $name = shift @{ $opts };

	if ($type eq "bgf") {
		$self->$readBGF($name, @_, @{ $opts });
	} elsif ($type eq "mol2") {
		$self->$readMOL2($name, @_, @{ $opts });
	} elsif ($type eq "pdb" or $type eq "pqr") {
		$self->$readPDB($name, @_, @{ $opts });
	} elsif ($type eq "xyz") {
		$self->$readXYZ($name, @_, @{ $opts });
	} elsif ($type eq "psf") {
		$self->$readCharmmPSF($name, @_, @{ $opts });
	} elsif ($type eq "prmtop") {
		$self->$readAmberPRMTOP($name, @_, @{ $opts });
	} elsif ($type eq "gro") {
		$self->$readGromacs($name, @_, @{ $opts });
	} else {
		die "ERROR: Filetype '$type' cannot be read\n";
	}
	&getMols($self);
}

sub write {
	my $self = shift;
	my $name = shift;
	my $type = lc shift;
	my ($count, $fmt, $currAtom, $outLine, $fields, $i, $cField, $currBond, $val, $bdata, $npairs, $adata);

	open OUTDATA, "> $name" or die "Cannot write to $name: $!\n";
	# bond data

	($fields, $fmt) = $self->$getAtomFileFormat($type);

	$count = 1;
	$currAtom = $self->atomlist->first;
	$adata = "";
	while ($count <= $self->count->{"atoms"}) {
		$outLine = "";
		$currAtom->data->("index", $count);
		$currAtom->data->("numbonds", $currAtom->data->bondlist->count);
		if($currAtom->data->bondlist->count > 10) {
			$currAtom->data->("numbonds", 9);
		}	
		for $i (0 .. $#{ $fields }) {
			$cField = $fields->[$i];
			$val = $currAtom->data->$cField;
			$outLine .= sprintf($fmt->[$i], $val);
		}
		$adata .= "$outLine\n";
		$currAtom = $currAtom->next;
		$count++;
	}

	($bdata,$npairs) = $self->$writeBondsData($type);
	print OUTDATA $self->$createHeader($type, $npairs);
	print OUTDATA $adata;
	print OUTDATA $bdata;
	print OUTDATA $self->header->{"footer"} if ($self->header->{"footer"});
	close OUTDATA;
}

sub testFile {
	my $self = shift;

	die "Expected file in testFile!\n" if (scalar(@_) == 0);
	my $name = shift;
	my $terminate = shift;
	my $verbose = shift;

	$terminate = 1 if(! defined($terminate));
	$verbose = 1 if (! defined($verbose));
	if ($terminate) {
		die "Cannot access file $name: $!\n"
			if (! -e $name or ! -r $name or ! -T $name);
	} else {
		print "WARNING: Cannot access file $name: $! ..." if ($verbose);
	}
}

sub deleteAtom {
	my $self = shift;
	my $currAtom = shift;
	my ($i, $aVal);
	return () if (! $currAtom);

	my $currBond = $currAtom->bondlist->last;
	my $atmCounter = $self->count->{"atoms"};
	my $bondCounter = $self->count->{"bonds"};
	my $atmCode = $currAtom->code;
	my $fields = $currAtom->fields;
	chop $fields;

	while ($currBond and $bondCounter>0) {
		$currBond->data->bondlist->remove($currAtom);
		$currAtom->bondlist->remove($currBond->data);
		$currBond = $currAtom->bondlist->last;
		$bondCounter--;
	}
	delete($self->atoms->{ $currAtom->code });
	$self->atomlist->delete($currAtom);
	$self->count->{"bonds"} = $bondCounter;

	while ($fields =~ /(\S+)/g) {
		$aVal = $currAtom->$1;
		delete $self->shash->{$1}{$aVal}{$atmCode} if (exists($self->shash->{$1}) and defined($aVal));
	}

	$currAtom = undef;
	$atmCounter--;
	$self->count->{"atoms"} = $atmCounter;
	$self->count->{"bonds"} = $bondCounter;
}

sub updateAtomIndex {
	my $self = shift;
	
	my $count = 1;
	my $currAtom = $self->atomlist->first;
	while ($count <= $self->count->{"atoms"}) {
		$currAtom->data->("index", $count);
		$currAtom->data->("numbonds", $currAtom->data->bondlist->count);
		$currAtom = $currAtom->next;
		$count++;
	}
}
	
sub deleteMol {
	my $self = shift;
	my $currMol = shift;
	my ($index, $currAtom);
	return () if (! $currMol);

	$index = $currMol->id;
	while (@{ $currMol->atom }) {
		$currAtom = shift @{ $currMol->atom };
		$self->deleteAtom($currAtom);
	}
	delete($self->mollist->{$index});
	$currMol = undef;
	$self->count->{molecule}--;
}

sub getFileName {
	my $self = shift;
	my $fileName = shift;

	my $saveName = basename($fileName);
	$saveName =~ /(\S+)\.(\S+)$/;
	$saveName = "${1}_mod.${2}";

	return $saveName;
}

sub getFileType {
	my $self = shift;
	my $fileName = shift;

	#first strip leading spaces
	$fileName =~ /^\s*(\S+)/; 
	my ($fileType) = $1; 
	#now get extension
	$fileType =~ /\.(\w+)$/;
	$fileType = lc $1;

	return $fileType;
}

sub CoM {
	my $self = shift;
	my $mols = shift;
	my ($mass, $com, $i, $tmp);

	$mols = $self->mollist if (!defined($mols) or (! ref($mols) and $mols =~ /all/i));
	return if (! $mols);
	$mass = 0;
	for $i (keys %{ $mols }) {
		$tmp = $mols->{$i}->CoM;
		for ("x", "y", "z") {
			$com->{$_} += ($tmp->{$_} * $mols->{$i}->mass);
		}
		$mass += $mols->{$i}->mass;
	}
	for ("x", "y", "z") {
		$com->{$_} /= $mass;
	}

	return $com;
}

sub moveMol {
	my $self = shift;
	my $mols = shift;
	my $point = shift;
	my ($com, $atom, $i, $offset);

	return if (! defined($mols));
	$mols = $self->mollist if (! ref($mols) and $mols =~ /all/i);
	$com = $self->CoM($mols);
	$point = $self->$findPointCoord($point, $com);
	return if (! $point);

	for $i ("x", "y", "z") {
		next if (! exists($point->{$i}));
		$offset = "+ " . $point->{$i} if ($point->{$i} >= 0);
		$offset = "- " . (-1*$point->{$i}) if ($point->{$i} < 0);
		for (keys %{ $mols }) {
			$mols->{$_}->trans($i, $offset);
		}
	}
}

sub imageMol {
	my $self = shift;
	my $mol = shift;

	return if (! $mol or ! $self->cell->valid);
	$mol = $self->mollist if (! defined($mol) or $mol =~ /all/i);

}

sub unwrapMol {
	my $self = shift;
	my $mol = shift;

	return if (! $mol or ! $self->cell->valid);
	$mol = $self->mollist if (! defined($mol) or $mol =~ /all/i);

}

sub stressCell {
# either compresses (cell length < 1) or expands (> 1) cell along specified axis(s)
	my $self = shift;
	my $axesStr = shift;
	my $factor = shift;
	my ($com, $i, $mol, $dim, $offset, $axes, $tmp, $map);

	return if (! defined($factor) or ($factor !~ /^\d+/ and $factor !~ /^\.\d+/));
	$map = { "x" => "a", "y" => "b", "z" => "c", };
	while ($axesStr =~ /(x|y|z)/gi) {
		$axes->{$1} = 1;
	}
	for $i (1 .. $self->count->{molecule}) {
		$mol = $self->molecule($i);
		$com = $mol->CoM;
		$tmp = ();
		for $dim (keys %{ $axes }) {
			$offset->{$dim} = $com->{$dim} * ($factor - 1);
		}
		$tmp->{$mol->id} = $mol;
		$self->moveMol($tmp, $offset);
	}
	for $dim (keys %{ $axes }) {
		$self->cell->{ $map->{$dim} } *= $factor 
			if(exists($map->{$dim}) and exists($self->cell->{ $map->{$dim} }));
	}
}

sub getExtrema {
	my $self = shift;
	my $which = shift;
	my ($i, $map, $atom, $extrema);

	return if ($which !~ /(min|max)/i);
	$which = lc $1;
	return if (! exists($self->vdwbox->{a}{$which}));
	$map = ({ "x" => "a", "y" => "b", "z" => "c" });
	for $i (1 .. $self->count->{atoms}) {
		$atom = $self->atoms($i);
		for ("x", "y", "z") {
			$self->vdwbox->{$map->{$_}}{max} = $atom->$_ if ($atom->$_ > $self->vdwbox->{$map->{$_}}{max});
			$self->vdwbox->{$map->{$_}}{min} = $atom->$_ if ($atom->$_ < $self->vdwbox->{$map->{$_}}{min});
		}
	}
	for $i ("a", "b", "c") {
		$extrema->{$i} = $self->vdwbox->{$i}{$which};
	}

	return $extrema;
}

sub dist {
	my $self = shift;
	my $atom1 = shift;
	my $atom2 = shift;
	my $dist = 0;
	return () if (! defined($atom1) or ! defined($atom2));
	$dist = sqrt(($atom1->x - $atom2->x)**2 + ($atom1->y - $atom2->y)**2 + ($atom1->z - $atom2->z)**2);

	return $dist;
}

sub addMass {
	my $self = shift;
	my $atom = shift;
	my $mass = shift;
	return () if (! defined($atom) or ! defined($mass));

	$atom->("mass",$mass);
}

sub ShuffleList {
# fisher_yates_shuffle
	my $array = $_[0];
	my $i;
	for ($i = @$array; --$i; ) {
		my $j = int rand ($i+1);
		next if $i == $j;
		@$array[$i,$j] = @$array[$j,$i];
	}
}

sub findAtoms {
	my $self = shift;
	my $atomSel = shift;
	my $tot = shift;
	my $molOpt = shift;
	my $randOpt = shift;

	my ($atomID, $atom, $tmpSel, $fields, $found);
	my ($i, $val, $count, $atomList);

	return () if (! $self or ! $atomSel);

	@{ $fields } = sort { $a cmp $b } keys %{ $self->shash };
	unshift @{ $fields }, "index";
	$tot = $self->count->{atoms} if (! $tot or $tot == -1);
	$count = 0;
	$found = ();
	@{ $atomList } = sort {$a<=>$b} keys %{ $self->shash->{index} };
	&ShuffleList($atomList) if ($randOpt);
	for $atomID (@{ $atomList }) {
		$atom = $self->atoms($atomID);
		next if (! defined($atom));
		$tmpSel = $atomSel;
		for $i (@{ $fields }) {
			$val = $atom->$i;
			$tmpSel =~ s/$i/$val/ig;
		}
		if (eval($tmpSel)) {
			if(! $molOpt) {
				$found->{$atom->index} = $atom;
				$count++;
				last if ($count == $tot);
			} else {
				next if (exists($found->{$atom->molid}));
				$found->{$atom->molid} = $self->molecule($atom->molid);
				$count++;
				last if ($count == $tot);
			}
		}
	}

	return $found;
}

sub rotate {
	my $self = shift;
	my $rotMatrix = shift;
	my ($i, $atom, $j, $currPos, $newPos, $dims, $cdim);

	return 0 if (! $rotMatrix);

	$dims = (["x", "y", "z"]);
	for $i (1 .. $self->count->{atoms}) {
		$atom = $self->atoms($i);
		for $j (0 .. 2) {
			$cdim = $dims->[$j];
			$currPos = $atom->$cdim;
			$newPos = ($atom->x * $rotMatrix->[$j][0] + $atom->y * $rotMatrix->[$j][1] + $atom->z * $rotMatrix->[$j][2]);
			$atom->move($cdim, $newPos);
		}
	}

}

sub getMols {
	my $self = shift;
	my ($i, $USED, $molCounter, $currAtom, $newMol);

	$molCounter = 0;
	$USED->{bcounter} = 0;
	for $i (1 .. $self->count->{"atoms"}) {
		if (! exists($USED->{$i})) {
			$molCounter++;
			$newMol = Molecule->create();
			$newMol->("id", $molCounter);
			$currAtom = $self->atoms($i);
			&addMolAtoms($self, $currAtom, $newMol, $molCounter, \%{ $USED });
			$self->mollist->{$molCounter} = $newMol;
		}
	}
	$self->count->{molecule} = $molCounter;
}

sub addMolAtoms {
	my $self = shift;
	my $catom = shift;
	my $mol = shift;
	my $molCounter = shift;
	my $used = shift;
	my ($i, $batom, $bdata);
	
	$mol->addAtom($catom);
	$catom->("molid", $molCounter);
	$used->{$catom->index} = 1;

	$batom = $catom->bondlist->first;
	for $i (1 .. $catom->bondlist->count) {
		if(exists($used->{$batom->data->index})) {
			$batom = $batom->next;
			next;
		}
		$used->{"bcounter"}++;
		$bdata = $batom->data;
		&addMolAtoms($self,$bdata, $mol, $molCounter, \%{ $used });
		$batom = $batom->next;
	}
};
1;
