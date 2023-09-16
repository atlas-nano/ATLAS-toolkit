#!/usr/bin/perl -w

# version 01-11-2002

# given a master file of multiple ligands in mol2 format and a template
# receptor bgf file, create a hybrid bgf file for each ligand, and run
# mpsim ONE_EF on it.

die "Usage: clean_bgf.pl ligands.mol2 ligands.bgf\n" unless @ARGV;

# master table of Sybyl to Dreiding atom types

# jaw - changed the following types ( 1-1-05 )
# C.cat/C.CAT -> C_R  (not C_3)
# N.am/N.AM -> N_R  (not N_2)
# N.pl3/N.PL3 -> N_R  (not N_1)

%to_bgf = ( "Du"    => "B_3",     # not X.  Why? It was H_ in previous versions
            "As"    => "As",
            "Ca"    => "Ca",
            "Fe"    => "Fe",
            "Na"    => "Na",
            "Cl"    => "Cl",
            "Mg"    => "Mg",
            "C.3"   => "C_3",
            "C.2"   => "C_2",
            "C.2ar" => "C_R",
            "C.2AR" => "C_R",
            "C.ar"  => "C_R",
            "C.AR"  => "C_R",
            "C.cat" => "C_R",
            "C.CAT" => "C_R",
            "C.1"   => "C_1",
            "N.4"   => "N_3",
            "N.3"   => "N_3",
            "N.2"   => "N_2",
            "N.am"  => "N_R",
            "N.AM"  => "N_R",
            "N.2ar" => "N_R",
            "N.2AR" => "N_R",
            "N.ar"  => "N_R",
            "N.AR"  => "N_R",
            "N.3ar" => "N_R",
            "N.3AR" => "N_R",
            "N.1"   => "N_1",
            "N.pl3" => "N_R",
            "N.PL3" => "N_R",
            "O.3"   => "O_3",
            "O.2"   => "O_2",
            "O.2ar" => "O_R",
            "O.2AR" => "O_R",
            "O.ar"  => "O_R",
            "O.AR"  => "O_R",
            "O.1"   => "O_1",
            "O"     => "O_3",
            "O.co2" => "O_2",
            "O.CO2" => "O_2",
            "S.3"   => "S_3",
            "P.3"   => "P_3",
            "P.3p"  => "P_4",  # phosphate
            "H"     => "H_",
            "Hhb"   => "H___A",
            "D"     => "H_",
            "Dhb"   => "H___A",
            "F"     => "F_",
            "Br"    => "Br",
            "I"     => "I_",
            "Li"    => "Li",
            "Al"    => "Al3",
            "Si"    => "Si3",
            "ZN"    => "Zn",   # check capitalization
            "Zn"    => "Zn",
            "S.2"   => "S_3",
            "S.o"   => "S_3",
            "S.O"   => "S_3",
            "S.o2"  => "S_3",
            "S.O2"  => "S_3",
            "As" => "As3",
          );

$mol_file = shift @ARGV;
$bgf_file = shift @ARGV;

open(MOL,"< $mol_file") || die "File $mol_file not found\n";
open(BGF,"< $bgf_file") || die "File $bgf_file not found\n";


LINE: for (;;) {
  last if eof(MOL);
  
  # go to top of next block
  do {
    $_ = <MOL>;
  } until (/^\@\<TRIPOS\>(\w+)/);

  $keyword = $1;

  if  ($keyword eq "MOLECULE")     { $keyword = &read_mol;    }
  if  ($keyword eq "ATOM")         { $keyword = &read_atoms;  }
  if  ($keyword eq "BOND")         { $keyword = &read_bonds;  }
  if  ($keyword eq "SUBSTRUCTURE") { $keyword = &read_substr; }

# read in bgf file, replace atomtypes and spit out new file

  do {
    $_ = <BGF>;
    print;
  } until (/^FORMAT ATOM/ || eof(BGF)) ;

# check to see if bgf file has FORMAT ATOM CARD
# to avoid getting stuck because of incomplete files

  if (/^FORMAT ATOM/) {
     ($ff_off, $ff_len) = &parse_format_atom($_);
  
#    die "off = $ff_off   len=$ff_len";

    $atom = 0;
    while (<BGF>) {
      last if (! (/^ATOM/ || /^HETATM/));
      $atom++;
#      print "atom = $atom   type = $type[$atom]  bgf=$to_bgf{$type[$atom]}\n";
      substr($_,$ff_off,$ff_len) = sprintf "%-${ff_len}s", $to_bgf{$type[$atom]};
      print;
    }

  } #end of if for FORMAT ATOM

  if (/^FORMAT\s+CONECT/) {print;}   # FORMAT CONECT line

  do {
    $_ = <BGF>;
    print;
  } until (/^END/ || eof(BGF));

#  die "end of first structure";

}

sub read_mol {
  do {
    $_ = <MOL>;
  } until /^\@<TRIPOS>(\w+)/;
  return $1;
}

sub read_substr {
  do {
    $_ = <MOL>;
  } until (/^\@<TRIPOS>(\w+)/ || eof(MOL));
  return $1;
}

sub read_atoms {

  $atom=0;
  while (<MOL>) {
    ++$atom;
    if (/^\@<TRIPOS>(\w+)/) { return $1; }
    @fields = split;
    if ((scalar @fields) != 9) 
      { die "No of fields is not 9: ",scalar @fields; }
    $type[$atom] = $fields[5];
    $atmlabel[$atom] = $fields[1];
#    print "READ:  atom=$atom  type=$type[$atom]\n";
    if (! exists $to_bgf{$type[$atom]}) { die "$type[$atom] is not defined!"; }
  }
}


sub read_bonds {

  $n_P_o = 0; $i_P = 0;
  $i_Du = 0; $n_Du = 0;

  # change types on the fly, including hydrogens

  while (<MOL>) {
    if (/^\@<TRIPOS>(\w+)/) { return $1; }
    (undef,$i,$j,$bond) = split;

#    print "parsing i=$i  j=$j  bond=$bond\n";
    # case 0: Du
    if ($type[$j] eq "Du" ) {
       $i_Du = $j; $n_Du++;
       # Atom is a sulphur
       if ($atmlabel[$j] =~ /^S/) {
          $type[$i_Du] = "S.3";
       }
       # Atom is a P (check for P4 is done later)
       if ($atmlabel[$j] =~ /^P/) {
          $type[$i_Du] = "P.3";
       }
       # Atom is a carbon
       if ($atmlabel[$j] =~ /^C/) {
          if ($n_Du <= 3) {
             $type[$i_Du] = "C." . $n_Du;
          }
          else {
             $type[$i_Du] = "C.3";
          }
       }
     #print "$i_Du type=$type[$j] atomlabel=$atmlabel[$j], n_Du=$n_Du\n"; #debug
    }

    if ($type[$i] eq "Du" ) {
       $i_Du = $i; $n_Du++;
       # Atom is a sulphur
       if ($atmlabel[$i] =~ /^S/) {
          $type[$i_Du] = "S.3";
       }
       # Atom is a P (check for P4 is done later)
       if ($atmlabel[$i] =~ /^P/) {
          $type[$i_Du] = "P.3";
       }
       # Atom is a carbon
       if ($atmlabel[$i] =~ /^C/) {
          if ($n_Du <= 3) {
             $type[$i_Du] = "C." . $n_Du;
          }
          else {
             $type[$i_Du] = "C.3";
          }
       }
     #print "$i_Du type=$type[$i_Du] atomlabel=$atmlabel[$i], n_Du=$n_Du\n"; #debug
    }


    # case 1:  $i is hydrogen
    if ($type[$i] eq "H" or $type[$i] eq "D") {
      if ($type[$j] =~ /^[ONS]/) { $type[$i] .= "hb"; }
    }

    # case 2:  $j is hydrogen
    if ($type[$j] eq "H" or $type[$j] eq "D") {
      if ($type[$i] =~ /^[ONS]/) { $type[$j] .= "hb"; }
    }

    # case 3: Phosphate
    if ($type[$j] eq "P.3" and ($type[$i] eq "O.2" or $type[$i] eq "O.3" or 
        $type[$i] eq "O.co2")) {
       $i_P = $j; $n_P_o++;
       if ($n_P_o == 4) {$type[$i_P] .= "p";}
    }
    if ($type[$i] eq "P.3" and ($type[$j] eq "O.2" or $type[$j] eq "O.3" or
        $type[$j] eq "O.co2")) {
       $i_P = $i; $n_P_o++;
       if ($n_P_o == 4) {$type[$i_P] .= "p";}
    }

    # case 4:  bond is aromatic
    if ($bond eq "ar") {
      if (( $type[$i] eq "C.2" or $type[$i] eq "N.2" or
            $type[$i] eq "N.3" or $type[$i] eq "O.2" ) && 
          (! ($type[$i] =~ /ar$/ or $type[$i] =~ /^[HD]/))) {
        $type[$i] .= "ar";
      }
      if (( $type[$j] eq "C.2" or $type[$j] eq "N.2" or
            $type[$j] eq "N.3" or $type[$j] eq "O.2" ) && 
          (! ($type[$j] =~ /ar$/ or $type[$j] =~ /^[HD]/))) {
        $type[$j] .= "ar";
      }
    }

  }
}

sub parse_format_atom {
  my $line = shift @_;

  my ($field, @fields);

  my ($ff_off, $ff_len) = (0,0);

  $line =~ s/^.*\(//;
  $line =~ s/\).*$//;

  # FF type is first field after the floats

  @fields = split ',', $line;

  # find xyz fields
  FIELDS: for (;;) {
    $field = shift @fields;
    if ($field =~ /f/) {
      ($mul,$len) = ($field =~ /(\d*)f(\d*)/);
      $ff_off += $mul*$len;
      for (;;) {
        $field = shift @fields || die "out of fields";
#        print "field = $field\n";
        if ($field =~ /a(\d+)/) {
          $ff_len = $1;
          last FIELDS;
        }
        if ($field =~ /(\d*)x/) {
          $off = $1 || 1;
          $ff_off += $off;
        };
      }

    }
    else {
      my ($off) = ($field =~ /(\d+)/);
      ($ff_off) += $off;
    }
  }
  return ($ff_off,$ff_len);
}

