package PDB;
use strict;
use Carp;
use FindBin qw($Bin);
use lib "$FindBin::Bin/Packages";
use lib "${Bin}/Packages/PDBsub/blib/arch";
use lib "${Bin}/Packages/PDBsub/blib/lib";
use PDB;
use Storable qw(dclone);
use PDBsub;
use Math::Trig;
use Math::Complex;

my %aa1 = (
  "ALA","A",
  "CYS","C",
  "CYX","C",
  "CYM","C",
  "ASP","D",
  "ASH","D",
  "GLU","E",
  "GLH","E",
  "PHE","F",
  "GLY","G",
  "HIS","H",
  "HID","H",
  "HIE","H",
  "HIP","H",
  "ILE","I",
  "LYS","K",
  "LEU","L",
  "MET","M",
  "MSE","M",
  "ASN","N",
  "PRO","P",
  "GLN","Q",
  "ARG","R",
  "SER","S",
  "THR","T",
  "VAL","V",
  "TRP","W",
  "TYR","Y"
);

my %hbond = (
  ALA_N   =>  1,
  ALA_O   => -1,

  CYS_N   =>  1,
  CYS_O   => -1,
  CYS_SG  =>  2,
  CYX_N   =>  1,
  CYX_O   => -1,
  CYX_SG  => -1,
  CYM_N   =>  1,
  CYM_O   => -1,
  CYM_SG  => -1,

  ASP_N   =>  1,
  ASP_O   => -1,
  ASP_OD1 => -1,
  ASP_OD2 => -1,
  ASH_N   =>  1,
  ASH_O   => -1,
  ASH_OD1 =>  2,
  ASH_OD2 =>  2,

  GLU_N   =>  1,
  GLU_O   => -1,
  GLU_OE1 => -1,
  GLU_OE2 => -1,
  GLH_N   =>  1,
  GLH_O   => -1,
  GLH_OE1 =>  2,
  GLH_OE2 =>  2,

  PHE_N   =>  1,
  PHE_O   => -1,

  GLY_N   =>  1,
  GLY_O   => -1,

  HIS_N   =>  1,
  HIS_O   => -1,
  HIS_ND1 =>  2,
  HIS_NE2 =>  2,
  HID_N   =>  1,
  HID_O   => -1,
  HID_ND1 =>  1,
  HID_NE2 => -1,
  HIE_N   =>  1,
  HIE_O   => -1,
  HIE_ND1 => -1,
  HIE_NE2 =>  1,
  HIP_N   =>  1,
  HIP_O   => -1,
  HIP_ND1 =>  1,
  HIP_NE2 =>  1,

  ILE_N   =>  1,
  ILE_O   => -1,

  LYS_N   =>  1,
  LYS_O   => -1,
  LYS_NZ  =>  1,

  LEU_N   =>  1,
  LEU_O   => -1,

  MET_N   =>  1,
  MET_O   => -1,
  MSE_N   =>  1,
  MSE_O   => -1,

  ASN_N   =>  1,
  ASN_O   => -1,
  ASN_OD1 => -1,
  ASN_ND2 =>  1,

  PRO_O   => -1,

  GLN_N   =>  1,
  GLN_O   => -1,
  GLN_OE1 => -1,
  GLN_NE2 =>  1,

  ARG_N   =>  1,
  ARG_O   => -1,
  ARG_NE  =>  1,
  ARG_NH1 =>  1,
  ARG_NH2 =>  1,

  SER_N   =>  1,
  SER_O   => -1,
  SER_OG  =>  2,

  THR_N   =>  1,
  THR_O   => -1,
  THR_OG1 =>  2,

  VAL_N   =>  1,
  VAL_O   => -1,

  TRP_N   =>  1,
  TRP_O   => -1,
  TRP_NE1 =>  1,

  TYR_N   =>  1,
  TYR_O   => -1,
  TYR_OZ  =>  2,

  HOH_O   =>  2,
  WAT_O   =>  2
);

use subs qw(fname natom seqres recname atnum atnam alt resnam resnam1 insert chain resnum x y z q b segname element ter);

sub new {
  my $invocant=shift;
  my $class=ref($invocant) || $invocant;
  my $self={"fname"   => "",
            "natom"   => 0,
            "seqres"  => {},
            "recname" => [],
            "atnum"   => [],
            "atnam"   => [],
            "alt"     => [],
            "resnam"  => [],
            "resnam1" => [],
            "chain"   => [],
            "resnum"  => [],
            "insert"  => [],
            "x"       => [],
            "y"       => [],
            "z"       => [],
            "q"       => [],
            "b"       => [],
            "segname" => [],
            "element" => [],
            "ter"     => [],
            @_};
  bless($self,$class);
  return $self;
}

sub init {
  my $self=shift;
  $self->{natom}   = 0;
  $self->{seqres}  = {};
  $self->{recname} = [];
  $self->{atnum}   = [];
  $self->{atnam}   = [];
  $self->{alt}     = [];
  $self->{resnam}  = [];
  $self->{resnam1} = [];
  $self->{chain}   = [];
  $self->{resnum}  = [];
  $self->{insert}  = [];
  $self->{x}       = [];
  $self->{y}       = [];
  $self->{z}       = [];
  $self->{q}       = [];
  $self->{b}       = [];
  $self->{segname} = [];
  $self->{ter}     = [];
}

sub copy {
  my $self=shift;
  my $b;
  my ($i,$j,$k);
  my %args=@_;
  my @list;

  if(defined($args{selection})) {
    @list=@{$args{selection}};
  } else {
    @list=(0 .. $self->{natom}-1);
  }
  $b=PDB->new();
  $b->{natom}=@list;
  $b->{fname}=$self->{fname};
  foreach $k (keys(%{$self->{seqres}})) {
    $b->{seqres}->{$k}=$self->{seqres}->{$k};
  }

  $i=0;
  foreach $j (@list) {
    $b->{recname}->[$i]=$self->{recname}->[$j];
    $b->{atnum}->[$i]  =$self->{atnum}->[$j]  ;
    $b->{atnam}->[$i]  =$self->{atnam}->[$j]  ;
    $b->{alt}->[$i]    =$self->{alt}->[$j]    ;
    $b->{resnam}->[$i] =$self->{resnam}->[$j] ;
    $b->{resnam1}->[$i]=$self->{resnam1}->[$j];
    $b->{chain}->[$i]  =$self->{chain}->[$j]  ;
    $b->{resnum}->[$i] =$self->{resnum}->[$j] ;
    $b->{insert}->[$i] =$self->{insert}->[$j] ;
    $b->{x}->[$i]      =$self->{x}->[$j]      ;
    $b->{y}->[$i]      =$self->{y}->[$j]      ;
    $b->{z}->[$i]      =$self->{z}->[$j]      ;
    $b->{q}->[$i]      =$self->{q}->[$j]      ;
    $b->{b}->[$i]      =$self->{b}->[$j]      ;
    $b->{segname}->[$i]=$self->{segname}->[$j];
    $b->{element}->[$i]=$self->{element}->[$j];
    $b->{ter}->[$i]    =$self->{ter}->[$j]    ;
    $i++;
  }
  return $b;
}

sub append {
# fname and seqres are not modified.
  my $self=shift;
  my $b=shift;

  $self->{natom}+=$b->{natom};

  $self->{recname}=[@{$self->{recname}},@{$b->{recname}}];
  $self->{atnum}  =[@{$self->{atnum}}  ,@{$b->{atnum}}  ];
  $self->{atnam}  =[@{$self->{atnam}}  ,@{$b->{atnam}}  ];
  $self->{alt}    =[@{$self->{alt}}    ,@{$b->{alt}}    ];
  $self->{resnam} =[@{$self->{resnam}} ,@{$b->{resnam}} ];
  $self->{resnam1}=[@{$self->{resnam1}},@{$b->{resnam1}}];
  $self->{chain}  =[@{$self->{chain}}  ,@{$b->{chain}}  ];
  $self->{resnum} =[@{$self->{resnum}} ,@{$b->{resnum}} ];
  $self->{insert} =[@{$self->{insert}} ,@{$b->{insert}} ];
  $self->{x}      =[@{$self->{x}}      ,@{$b->{x}}      ];
  $self->{y}      =[@{$self->{y}}      ,@{$b->{y}}      ];
  $self->{z}      =[@{$self->{z}}      ,@{$b->{z}}      ];
  $self->{q}      =[@{$self->{q}}      ,@{$b->{q}}      ];
  $self->{b}      =[@{$self->{b}}      ,@{$b->{b}}      ];
  $self->{segname}=[@{$self->{segname}},@{$b->{segname}}];
  $self->{element}=[@{$self->{element}},@{$b->{element}}];
  $self->{ter}    =[@{$self->{ter}}    ,@{$b->{ter}}    ];
}

sub read {
  my $self=shift;
  my %args=@_;
  my $unwrap=0;

#delete all data
  $self->init();

#start reading
  if(defined($args{fname})) {
    $self->{fname}=$args{fname};
  }
  if(defined($args{unwrap})) {
    $unwrap=1;
  }
  my $i=0;
  my %seqres;
  my ($s,$c);
  my (@recname,@atnum,@atnam,@alt,@resnam,@resnam1,@chain,@resnum,@insert,@x,@y,@z,@q,@b,@segname,@element,@ter);
  my $fh;
  my $big;
  if(defined($args{fname})) {
    open($fh,$self->{fname});
  } else {
    $fh=$args{fh};
  }
  if (defined($args{big})) {
    $big=$args{big};
  } else {
    $big=0;
  }
  while(<$fh>) {
    chomp;
    last if(defined($args{firstmodel}) and /^ENDMDL/);
    if(/^ATOM/ or /^HETATM/) {
      if(/^ATOM/ && $big ) {
        $recname[$i] = substr($_, 0,4);
        $atnum[$i]   = substr($_, 4,7);
      } else {
        $recname[$i] = substr($_, 0,6);
        $atnum[$i]   = substr($_, 6,5);
      }
      $atnum[$i]   =~s/^\s+//g;
      $atnum[$i]   =~s/\s+$//g;
      if($unwrap) {
        $atnam[$i] = substr($_,13,3) . substr($_,12,1);
        $atnam[$i] =~s/\s+//g;
      } else {
        $atnam[$i]   = substr($_,12,4);
      }
      $atnam[$i]   =~s/^\s+//g;
      $atnam[$i]   =~s/\s+$//g;
      $alt[$i]     = substr($_,16,1);
      $resnam[$i]  = substr($_,17,3);
#     $resnam[$i]  =~s/^\s+//g;
#     $resnam[$i]  =~s/\s+$//g;
      if(defined($aa1{$resnam[$i]})) {
        $resnam1[$i]=$aa1{$resnam[$i]};
      } else {
        $resnam1[$i]="X";
      }
      $chain[$i]   = substr($_,21,1);
      $resnum[$i]  = substr($_,22,4);
      $resnum[$i]  =~s/^\s+//g;
      $resnum[$i]  =~s/\s+$//g;
      $insert[$i]  = substr($_,26,1);
      $x[$i]       = substr($_,30,8);
      $x[$i]       =~s/^\s+//g;
      $x[$i]       =~s/\s+$//g;
      $y[$i]       = substr($_,38,8);
      $y[$i]       =~s/^\s+//g;
      $y[$i]       =~s/\s+$//g;
      $z[$i]       = substr($_,46,8);
      $z[$i]       =~s/^\s+//g;
      $z[$i]       =~s/\s+$//g;
      if(length($_) >= 60) {
        $q[$i]       = substr($_,54,6);
        $q[$i]       =~s/^\s+//g;
        $q[$i]       =~s/\s+$//g;
      } else {
        $q[$i]       = 1.0;
      }
      if(length($_) >= 66) {
        $b[$i]       = substr($_,60,6);
        $b[$i]       =~s/^\s+//g;
        $b[$i]       =~s/\s+$//g;
      } else {
        $b[$i]       = 0.0;
      }
      if(length($_) >=76) {
        $segname[$i] = substr($_,72,4);
        $segname[$i] =~s/^\s+//g;
        $segname[$i] =~s/\s+$//g;
      } else {
        $segname[$i] = "";
      }
      if(length($_) >=78) {
        $element[$i] = substr($_,76,2);
        $element[$i] =~s/^\s+//g;
        $element[$i] =~s/\s+$//g;
      } else {
        $element[$i] = "";
      }
      $ter[$i]       = 0;
      ++$i;
    } elsif(/^TER/) {
      $ter[$i-1]     = 1;
    } elsif(/^SEQRES/) {
      $c=substr($_,11,1);
      $s=substr($_,18);
      $s =~ s/\s+$//g;
      $seqres{$c}.=$s;
    }
  }
  if(defined($args{fname})) {
    close($fh);
  }
  $self->{natom}=$i;
  $self->{seqres}={%seqres};
  $self->{recname}=[@recname];
  $self->{atnum}=[@atnum];
  $self->{atnam}=[@atnam];
  $self->{alt}=[@alt];
  $self->{resnam}=[@resnam];
  $self->{resnam1}=[@resnam1];
  $self->{chain}=[@chain];
  $self->{resnum}=[@resnum];
  $self->{insert}=[@insert];
  $self->{x}=[@x];
  $self->{y}=[@y];
  $self->{z}=[@z];
  $self->{q}=[@q];
  $self->{b}=[@b];
  $self->{segname}=[@segname];
  $self->{element}=[@element];
  $self->{ter}=[@ter];
}

sub write {
  my $self = shift;
  my %args=@_;
  my ($i,$len,$fh,$no_close,$atnam,$use_ter,$print_ter);
  my @list;
  my $prev_chain="";
  my $prev_rec="";
  my $prev_resnum="";
  my ($form,$big);

  if(defined($args{fname})) {
    $self->{fname}=$args{fname};
  }
  if(defined($args{selection})) {
    @list=@{$args{selection}};
  } else {
    @list=(0 .. $self->{natom}-1);
  }
  if(defined($args{fh})) {
    $fh=$args{fh};
    $no_close=1;
  } else {
    open($fh,">".$self->{fname});
    $no_close=0;
  }
  if(defined($args{use_ter})) {
    $use_ter=$args{use_ter};
  } else {
    $use_ter=0;
  }
  if(defined($args{big})) {
    $big=$args{big};
  } else {
    $big=0;
  }
  foreach $i (@list) {
    $print_ter=0;
    $len=length($self->{atnam}->[$i]);
    if($len < 4 && $self->{atnam}->[$i] =~ /^\D/) {
      $atnam=sprintf(" %-3s",$self->{atnam}->[$i]);
    } elsif($len == 4 && defined($args{wrap})) {
      $self->{atnam}->[$i] =~ /^(...)(.)/;
      $atnam=sprintf("%s%s",$2,$1);
    } else {
      $atnam=sprintf("%-4s",$self->{atnam}->[$i]);
    }
    unless($prev_chain) {
      $prev_chain=$self->{chain}->[$i];
    }
    if($prev_chain ne $self->{chain}->[$i]) {
      $print_ter=1;
      $prev_chain=$self->{chain}->[$i];
    }
    unless($prev_rec) {
      $prev_rec=$self->{recname}->[$i];
    }
    if($prev_rec ne $self->{recname}->[$i]) {
      $print_ter=2;
      $prev_rec=$self->{recname}->[$i];
    }
    unless($prev_resnum) {
      $prev_resnum=$self->{resnum}->[$i];
    }
    if($prev_resnum != $self->{resnum}->[$i] &&
       $self->{recname}->[$i] eq "HETATM") {
      $print_ter=3;
      $prev_resnum=$self->{resnum}->[$i];
    }
    if($print_ter && !$use_ter) {
#     printf($fh "TER %d\n",$print_ter);
      printf($fh "TER\n");
    }
    if($self->{recname}->[$i] eq "ATOM" && $big) {
      $form="%-4s%7d %-4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s\n";
    } else {
      $form="%-6s%5d %-4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s\n";
    }
    printf($fh $form,
           $self->{recname}->[$i],$self->{atnum}->[$i],
           $atnam,$self->{alt}->[$i],$self->{resnam}->[$i],
           $self->{chain}->[$i],$self->{resnum}->[$i],$self->{insert}->[$i],
           $self->{x}->[$i],$self->{y}->[$i],$self->{z}->[$i],
           $self->{q}->[$i],$self->{b}->[$i],
           $self->{segname}->[$i],$self->{element}->[$i]);
    if($use_ter && $self->{ter}->[$i] == 1) {
      printf($fh "TER\n");
    }
  }
  if($no_close == 0) {
    printf($fh "END\n");
    close($fh);
  }
}

sub selection_atnam {
  my $self=shift;
  my %a;
  my @list;
  my $i;

  foreach $i (@_) {
    $a{$i}=1;
  }
  for($i=0;$i<$self->{natom};$i++) {
    if($a{$self->{atnam}->[$i]}) {
      push(@list,$i);
    }
  }
  return @list;
}

sub selection_resnum {
  my $self=shift;
  my %a;
  my @list;
  my $i;

  foreach $i (@_) {
    $a{$i}=1;
  }
  for($i=0;$i<$self->{natom};$i++) {
    if($a{$self->{resnum}->[$i]}) {
      push(@list,$i);
    }
  }
  return @list;
}

sub selection_alt {
  my $self=shift;
  my %a;
  my @list;
  my $i;

  foreach $i (@_) {
    $a{$i}=1;
  }
  for($i=0;$i<$self->{natom};$i++) {
    if($a{$self->{alt}->[$i]}) {
      push(@list,$i);
    }
  }
  return @list;
}

sub selection_chain {
  my $self=shift;
  my %a;
  my @list;
  my $i;

  foreach $i (@_) {
    $a{$i}=1;
  }
  for($i=0;$i<$self->{natom};$i++) {
    if($a{$self->{chain}->[$i]}) {
      push(@list,$i);
    }
  }
  return @list;
}

sub selection_resnam {
  my $self=shift;
  my %a;
  my @list;
  my $i;

  foreach $i (@_) {
    $a{$i}=1;
  }
  for($i=0;$i<$self->{natom};$i++) {
    if($a{$self->{resnam}->[$i]}) {
      push(@list,$i);
    }
  }
  return @list;
}

sub selection_sidechain {
  my $self=shift;
  my @list;
  my ($i,$atnam);

  for($i=0;$i<$self->{natom};$i++) {
    $atnam=$self->{atnam}->[$i];
    if($atnam ne "N"   and $atnam ne "CA"  and $atnam ne "C"  and
       $atnam ne "O"   and $atnam ne "H"   and $atnam ne "H1" and
       $atnam ne "H2"  and $atnam ne "H3"  and $atnam ne "HA" and
       $atnam ne "HA1" and $atnam ne "HA2" and $atnam ne "HA3"and
       $atnam ne "1HA" and $atnam ne "2HA" and $atnam ne "3HA"and
       $atnam ne "1H"  and $atnam ne "2H"  and $atnam ne "3H") {
      push(@list,$i);
    }
  }
  return @list;
}

sub selection_heavy {
  my $self=shift;
  my @list;
  my ($i,$atnam);

  for($i=0;$i<$self->{natom};$i++) {
    $atnam=$self->{atnam}->[$i];
    if(!($atnam =~ /^H/) and !($atnam =~ /^\d/)) {
      push(@list,$i);
    }
  }
  return @list;
}

sub selection_element {
  my $self=shift;
  my %a;
  my @list;
  my $i;

  foreach $i (@_) {
    $a{$i}=1;
  }
  for($i=0;$i<$self->{natom};$i++) {
    if($a{$self->{element}->[$i]}) {
      push(@list,$i);
    }
  }
  return @list;
}

sub selection_all {
  my $self=shift;
  my @list;
  my $i;

  for($i=0;$i<$self->{natom};$i++) {
    push(@list,$i);
  }
  return @list;
}

sub selection_intersection {
  my $self=shift;
  my @list;
  my ($i,$l);
  my %a;

  foreach $l (@_) {
    if(!ref($l)) {
      croak "Error in selection_intersection().";
    }
    foreach $i (@{$l}) {
      $a{$i}++;
    }
  }
  foreach $i (sort {$a <=> $b} (keys(%a))) {
    if($a{$i} == @_) {
      push(@list,$i);
    }
  }
  return @list;
}

sub selection_union {
  my $self=shift;
  my (@list,$i,$l,%a);

  foreach $l (@_) {
    if(!ref($l)) {
      croak "Error in selection_union().";
    }
    foreach $i (@{$l}) {
      $a{$i}++;
    }
  }
  foreach $i (sort {$a <=> $b} (keys(%a))) {
    push(@list,$i);
  }
  return @list;
}

sub select {
  my $self=shift;
  my %sel=@_;
  my (@list1,@list2,@list3,@list4,@list5,@list6);
  my ($i,$atnam);

  if(defined($sel{"atnam"})) {
    if(ref($sel{"atnam"})) {
      @list1=$self->selection_atnam(@{$sel{"atnam"}});
    } elsif($sel{"atnam"} eq "sidechain") {
      @list1=$self->selection_sidechain();
    } elsif($sel{"atnam"} eq "heavy") {
      @list1=$self->selection_heavy();
    } else {
      croak "selection atnam => ". $sel{"atnam"} ." not supported";
    }
  } else {
    @list1=$self->selection_all();
  }
  if(defined($sel{"resnum"}) and ref($sel{"resnum"})) {
    @list2=$self->selection_resnum(@{$sel{"resnum"}});
  } else {
    @list2=$self->selection_all();
  }
  if(defined($sel{"resnam"}) and ref($sel{"resnam"})) {
    @list3=$self->selection_resnam(@{$sel{"resnam"}});
  } else {
    @list3=$self->selection_all();
  }
  if(defined($sel{"chain"}) and ref($sel{"chain"})) {
    @list4=$self->selection_chain(@{$sel{"chain"}});
  } else {
    @list4=$self->selection_all();
  }
  if(defined($sel{"alt"}) and ref($sel{"alt"})) {
    @list5=$self->selection_alt(@{$sel{"alt"}});
  } else {
    @list5=$self->selection_all();
  }
  if(defined($sel{"element"}) and ref($sel{"element"})) {
    @list6=$self->selection_element(@{$sel{"element"}});
  } else {
    @list6=$self->selection_all();
  }
  return $self->selection_intersection([@list1],[@list2],[@list3],[@list4],[@list5],[@list6]);
}

sub zone {
  my $self=shift;
  my %args=@_;
  my (@center,@list);
  my ($i,$cut,$cut2,$d2);

  if(defined($args{"center_id"})) {
    $i=$args{"center_id"};
    @center=($self->{x}->[$i],$self->{y}->[$i],$self->{z}->[$i]);
  } elsif(defined($args{"center_xyz"})) {
    @center=@{$args{"center_xyz"}};
  } else {
    croak "Specify center_id or center_xyz";
  }
  if(defined($args{"cut"})) {
    $cut=$args{"cut"};
    $cut2=$cut**2;
  } else {
    croak "Specify cut";
  }
  for($i=0;$i<$self->{natom};$i++) {
    $d2=($self->{x}->[$i]-$center[0])**2
       +($self->{y}->[$i]-$center[1])**2
       +($self->{z}->[$i]-$center[2])**2;
    if($d2 < $cut2) {
      push(@list,$i);
    }
  }
  return @list;
}

sub distance {
  my $self=shift;
  my ($xi,$xj,$yi,$yj,$zi,$zj,$d);
  my ($i,$j)=@_;

  $xi=$self->{x}->[$i];
  $yi=$self->{y}->[$i];
  $zi=$self->{z}->[$i];
  $xj=$self->{x}->[$j];
  $yj=$self->{y}->[$j];
  $zj=$self->{z}->[$j];
  $d=sqrt(($xi-$xj)**2+($yi-$yj)**2+($zi-$zj)**2);
  return $d;
}
  
sub angle {
  my $self=shift;
  my ($xi,$xj,$yi,$yj,$zi,$zj,$xk,$yk,$zk,$dij,$dkj,$p);
  my ($i,$j,$k)=@_;

  $xi=$self->{x}->[$i];
  $yi=$self->{y}->[$i];
  $zi=$self->{z}->[$i];
  $xj=$self->{x}->[$j];
  $yj=$self->{y}->[$j];
  $zj=$self->{z}->[$j];
  $xk=$self->{x}->[$k];
  $yk=$self->{y}->[$k];
  $zk=$self->{z}->[$k];
  $dij=sqrt(($xi-$xj)**2+($yi-$yj)**2+($zi-$zj)**2);
  $dkj=sqrt(($xk-$xj)**2+($yk-$yj)**2+($zk-$zj)**2);
  $p=($xi-$xj)*($xk-$xj)+($yi-$yj)*($yk-$yj)+($zi-$zj)*($zk-$zj);
  $p/=($dij*$dkj);
  return Math::Trig::acos($p)/3.14159265358979*180.0;
}

sub dihedral { 
  my $self=shift;
  my ($xi,$xj,$yi,$yj,$zi,$zj,$xk,$yk,$zk,$xl,$yl,$zl,$p,$phi);
  my ($ax,$ay,$az,$bx,$by,$bz,$cx,$cy,$cz,$dx,$dy,$dz,$ex,$ey,$ez,$fx,$fy,$fz);
  my ($i,$j,$k,$l)=@_;
    
  $xi=$self->{x}->[$i];
  $yi=$self->{y}->[$i];
  $zi=$self->{z}->[$i];
  $xj=$self->{x}->[$j];
  $yj=$self->{y}->[$j];
  $zj=$self->{z}->[$j];
  $xk=$self->{x}->[$k];
  $yk=$self->{y}->[$k];
  $zk=$self->{z}->[$k];
  $xl=$self->{x}->[$l];
  $yl=$self->{y}->[$l];
  $zl=$self->{z}->[$l];
  $ax=$xi-$xj;
  $ay=$yi-$yj;
  $az=$zi-$zj;
  $bx=$xk-$xj;
  $by=$yk-$yj;
  $bz=$zk-$zj;
  $cx=$xk-$xl;
  $cy=$yk-$yl;
  $cz=$zk-$zl;
  $dx=$ay*$bz-$az*$by;
  $dy=$az*$bx-$ax*$bz;
  $dz=$ax*$by-$ay*$bx;
  $ex=$by*$cz-$bz*$cy;
  $ey=$bz*$cx-$bx*$cz;
  $ez=$bx*$cy-$by*$cx;
  $p=$dx*$ex+$dy*$ey+$dz*$ez;
  $p/=sqrt(($dx*$dx+$dy*$dy+$dz*$dz)*($ex*$ex+$ey*$ey+$ez*$ez));
  $phi=Math::Trig::acos($p)/3.14159265358979*180.0;
  $fx=$dy*$ez-$dz*$ey;
  $fy=$dz*$ex-$dx*$ez;
  $fz=$dx*$ey-$dy*$ex;
  if($fx*$bx+$fy*$by+$fz*$bz > 0.0) {
    return $phi;
  } else {
    return -1.0*$phi;
  }
}

sub set_dihedral {
  my $self=shift;
  my %args=@_;
  my ($xi,$xj,$yi,$yj,$zi,$zj,$xk,$yk,$zk);
  my ($ax,$ay,$az,$bx,$by,$bz,$cx,$cy,$cz);
  my ($i,$j,$k,$l,$dih,$dih0,$theta,$ok3,$ok4,$n,$norm,$m,$p);
  my ($xm,$ym,$zm,$px,$py,$pz,$qx,$qy,$qz);
  my @list;


  if(defined($args{dih_atoms})) {
    ($i,$j,$k,$l)=@{$args{dih_atoms}};
  } else {
    croak "Specify dihedral angle atoms with dih_atoms option";
  }
  if(defined($args{angle})) {
    $dih=$args{angle};
  } else {
    croak "Specify angle with angle option";
  }
  if(defined($args{moving_atoms})) {
    @list=@{$args{moving_atoms}};
  } else {
    croak "Specify moving atoms with moving_atoms option";
  }
  $ok4=0;
  foreach $n (@list) {
    $ok4=1 if($n == $l);
  }
  if($ok4 == 0) {
    push(@list,$l);
  }

  $dih0=$self->dihedral($i,$j,$k,$l);
  printf(STDERR "Initial dihedral angle: %f\n",$dih0);

  $xi=$self->{x}->[$i];
  $yi=$self->{y}->[$i];
  $zi=$self->{z}->[$i];
  $xj=$self->{x}->[$j];
  $yj=$self->{y}->[$j];
  $zj=$self->{z}->[$j];
  $xk=$self->{x}->[$k];
  $yk=$self->{y}->[$k];
  $zk=$self->{z}->[$k];
  $ax=$xk-$xj;
  $ay=$yk-$yj;
  $az=$zk-$zj;
  $norm=sqrt($ax*$ax+$ay*$ay+$az*$az);
  $ax/=$norm;
  $ay/=$norm;
  $az/=$norm;
  $bx=$xi-$xj;
  $by=$yi-$yj;
  $bz=$zi-$zj;
  $p=$bx*$ax+$by*$ay+$bz*$az;
  $bx-=$p*$ax;
  $by-=$p*$ay;
  $bz-=$p*$az;
  $norm=sqrt($bx*$bx+$by*$by+$bz*$bz);
  $bx/=$norm;
  $by/=$norm;
  $bz/=$norm;
  $cx=$ay*$bz-$az*$by;
  $cy=$az*$bx-$ax*$bz;
  $cz=$ax*$by-$ay*$bx;
  $norm=sqrt($cx*$cx+$cy*$cy+$cz*$cz);
  $cx/=$norm;
  $cy/=$norm;
  $cz/=$norm;

  $theta=($dih-$dih0)*3.14159265358979/180.0;
  foreach $m (@list) {
    $xm=$self->{x}->[$m]-$xk;
    $ym=$self->{y}->[$m]-$yk;
    $zm=$self->{z}->[$m]-$zk;
    $px=$xm*$ax+$ym*$ay+$zm*$az;
    $py=$xm*$bx+$ym*$by+$zm*$bz;
    $pz=$xm*$cx+$ym*$cy+$zm*$cz;
    $qx=$px;
    $qy=cos($theta)*$py-sin($theta)*$pz;
    $qz=sin($theta)*$py+cos($theta)*$pz;
    $self->{x}->[$m]=$xk+$qx*$ax+$qy*$bx+$qz*$cx;
    $self->{y}->[$m]=$yk+$qx*$ay+$qy*$by+$qz*$cy;
    $self->{z}->[$m]=$zk+$qx*$az+$qy*$bz+$qz*$cz;
  }
  $dih0=$self->dihedral($i,$j,$k,$l);
  printf(STDERR "Final dihedral angle: %f\n",$dih0);
}

sub getseq {
  my $self=shift;
  my %args=@_;
  my ($old_resnum,$old_chain,$seq,$first,$chain,$resnum,$i,@seq3,$c);

  $c="";
  if(defined($args{chain})) {
    $c=$args{chain};
  }
  $first=1;
  for($i=0;$i<$self->{natom};$i++) {
    if($self->{recname}->[$i] eq "HETATM" && $self->{resnam}->[$i] ne "MSE") {
      next;
    }
    $chain=$self->{chain}->[$i];
    $resnum=$self->{resnum}->[$i];
    if($first and (!$c or $c eq $chain)) {
      push(@seq3,$self->{resnam}->[$i]);
      if(defined($aa1{$self->{resnam}->[$i]})) {
        $seq=$aa1{$self->{resnam}->[$i]};
      } else {
        printf(STDERR "getseq(): %s is not defined.\n",$self->{resnam}->[$i]);
      }
      $old_resnum=$resnum;
      $old_chain=$chain;
      $first=0;
    }
    if($c and $c eq $chain and $old_resnum != $resnum) {
      push(@seq3,$self->{resnam}->[$i]);
      if(defined($aa1{$self->{resnam}->[$i]})) {
        $seq.=$aa1{$self->{resnam}->[$i]};
      } else {
        printf(STDERR "getseq(): %s is not defined.\n",$self->{resnam}->[$i]);
      }
      $old_resnum=$resnum;
      $old_chain=$chain;
    }
    if(!$c and ($old_resnum != $resnum || $old_chain ne $chain)) {
      push(@seq3,$self->{resnam}->[$i]);
      if(defined($aa1{$self->{resnam}->[$i]})) {
        $seq.=$aa1{$self->{resnam}->[$i]};
      } else {
        printf(STDERR "getseq(): %s is not defined.\n",$self->{resnam}->[$i]);
      }
      $old_resnum=$resnum;
      $old_chain=$chain;
    }
  }
  if(defined($args{three})) {
    return @seq3;
  } else {
    return $seq;
  }
}

sub rmsd {
  my $a=shift;
  my %args=@_;
  my ($anum,$rmsd,$ii,$jj,$i,$temp);
  my ($b,@alist,@blist);

  if(!defined($args{comp_pdb})) {
    croak "Specify comp_pdb.";
  } else {
    $b=$args{comp_pdb};
  }
  if(!defined($args{ref_atoms})) {
    croak "Specify ref_atoms.";
  } else {
    $temp=$args{ref_atoms};
    @alist=@{$temp};
  }
  if(!defined($args{comp_atoms})) {
    @blist=@alist;
  } else {
    $temp=$args{comp_atoms};
    @blist=@{$temp};
  }

  if(@alist != @blist) {
    croak "Number of atoms differs.";
  }
  $rmsd=0.0;
  for($i=0;$i<@alist;$i++) {
    $ii=$alist[$i];
    $jj=$blist[$i];
    $rmsd+=($a->{x}->[$ii] - $b->{x}->[$jj])**2;
    $rmsd+=($a->{y}->[$ii] - $b->{y}->[$jj])**2;
    $rmsd+=($a->{z}->[$ii] - $b->{z}->[$jj])**2;
  }
  $rmsd/=@alist;
  return sqrt($rmsd);
}


sub fit {
  my $a=shift;
  my %args=@_;
  my ($list_num,$rmsd);
  my (@acrd,@bcrd);
  my ($b,$alist,$blist,$i,$j,$k,$index,$det,$sign,$x,$y,$z);
  my (@cma,@cmb,@R,@V,@e,@L,@U);

  if(!defined($args{moving_pdb})) {
    croak "Specify moving_pdb.";
  } else {
    $b=$args{moving_pdb};
  }
  if(!defined($args{ref_atoms})) {
    croak "Specify ref_atoms.";
  } else {
    $alist=$args{ref_atoms};
  }
  if(!defined($args{moving_atoms})) {
    $blist=$alist;
  } else {
    $blist=$args{moving_atoms};
  }

  $list_num=@{$alist};
  if($list_num != @{$blist}) {
    croak "Number of atom differ.";
  }

# Translate center of mass
  $i=0;
  @cma=(0.0,0.0,0.0);
  for($i=0;$i<$list_num;$i++) {
    $j=$alist->[$i];
    $cma[0]+=$acrd[$i*3  ]=$a->{x}->[$j];
    $cma[1]+=$acrd[$i*3+1]=$a->{y}->[$j];
    $cma[2]+=$acrd[$i*3+2]=$a->{z}->[$j];
  }
  for($i=0;$i<3;$i++) {
    $cma[$i]/=$list_num;
  }
  for($i=0;$i<$list_num;$i++) {
    $acrd[$i*3  ]-=$cma[0];
    $acrd[$i*3+1]-=$cma[1];
    $acrd[$i*3+2]-=$cma[2];
  }
  $i=0;
  @cmb=(0.0,0.0,0.0);
  for($i=0;$i<$list_num;$i++) {
    $j=$blist->[$i];
    $cmb[0]+=$bcrd[$i*3  ]=$b->{x}->[$j];
    $cmb[1]+=$bcrd[$i*3+1]=$b->{y}->[$j];
    $cmb[2]+=$bcrd[$i*3+2]=$b->{z}->[$j];
  }
  for($i=0;$i<3;$i++) {
    $cmb[$i]/=$list_num;
  }
  for($i=0;$i<$list_num;$i++) {
    $bcrd[$i*3  ]-=$cmb[0];
    $bcrd[$i*3+1]-=$cmb[1];
    $bcrd[$i*3+2]-=$cmb[2];
  }

# Calculate R matrix
  for($j=0;$j<3;$j++) {
    for($k=0;$k<3;$k++) {
      $index=$k*3+$j;
      $R[$index]=0.0;
      for($i=0;$i<$list_num;$i++) {
        $R[$index]+=$acrd[$i*3+$j]*$bcrd[$i*3+$k];
      }
    }
  }

# Right-handed or left-handed?
  $det=$R[0]*$R[4]*$R[8]+$R[3]*$R[7]*$R[2]+$R[6]*$R[5]*$R[1]
      -$R[6]*$R[4]*$R[2]-$R[7]*$R[5]*$R[0]-$R[8]*$R[1]*$R[3];
  if($det < 0.0) {
    $sign=-1.0;
  } else {
    $sign=1.0;
  }

# Perform singular value decomposition
  @e=(0.0,0.0,0.0);
  @V=(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  PDBsub::do_svd3(\@R,\@e,\@V);
  $e[0]=1.0/$e[0];
  $e[1]=1.0/$e[1];
  $e[2]=$sign/$e[2];

# Calculate inverse of lambda matrix
# V is transpose of P
  for($j=0;$j<3;$j++) {
    for($k=0;$k<3;$k++) {
      $index=$k*3+$j;
      $L[$index]=0.0;
      for($i=0;$i<3;$i++) {
        $L[$index]+=$e[$i]*$V[$j*3+$i]*$V[$k*3+$i];
      }
    }
  }

# Calculate rotation matrix
  for($j=0;$j<3;$j++) {
    for($k=0;$k<3;$k++) {
      $index=$k*3+$j;
      $U[$index]=0.0;
      for($i=0;$i<3;$i++) {
        $U[$index]+=$R[$i*3+$j]*$L[$k*3+$i];
      }
    }
  }

# Transform
  for($i=0;$i<$b->{natom};$i++) {
    $x=$b->{x}->[$i]-$cmb[0];
    $y=$b->{y}->[$i]-$cmb[1];
    $z=$b->{z}->[$i]-$cmb[2];
    $b->{x}->[$i]=$U[0]*$x+$U[3]*$y+$U[6]*$z+$cma[0];
    $b->{y}->[$i]=$U[1]*$x+$U[4]*$y+$U[7]*$z+$cma[1];
    $b->{z}->[$i]=$U[2]*$x+$U[5]*$y+$U[8]*$z+$cma[2];
  }

# Calculate RMSD
  $rmsd=0.0;
  for($i=0;$i<$list_num;$i++) {
    $j=$alist->[$i];
    $k=$blist->[$i];
    $rmsd+=($a->{x}->[$j]-$b->{x}->[$k])**2
          +($a->{y}->[$j]-$b->{y}->[$k])**2
          +($a->{z}->[$j]-$b->{z}->[$k])**2;
  }
  $rmsd/=$list_num;
  return sqrt($rmsd);
}

sub is_donor {
  my $self=shift;
  my ($i)=@_;
  my $name;

  $name=$self->{resnam}->[$i] . "_" . $self->{atnam}->[$i];
  if(!defined($hbond{$name})) {
    return 0;
  } elsif($hbond{$name} == 2 || $hbond{$name} == 1) {
    return 1;
  }
}

sub is_acceptor {
  my $self=shift;
  my ($i)=@_;
  my $name;

  $name=$self->{resnam}->[$i] . "_" . $self->{atnam}->[$i];
  if(!defined($hbond{$name})) {
    return 0;
  } elsif($hbond{$name} == 2 || $hbond{$name} == -1) {
    return 1;
  }
}

sub getchains {
  my $self=shift;
  my %chains=();
  my $i;

  for($i=0;$i<$self->{natom};$i++) {
    $chains{$self->{chain}->[$i]}++;
  }
  return keys(%chains);
}

sub getseqres {
  my $self=shift;
  my %args=@_;
  my @seq3;
  my ($c,$seq);
  my (@missing,@seq_atom);
  my ($i,$j);

  unless(defined($args{chain})) {
    croak "Specify chain.";
  }
  $c=$args{chain};
  $seq=$self->{seqres}->{$c};
  $seq=~s/^\s+//g;
  $seq=~s/\s+$//g;
  @seq3=split(/\s+/,$seq);
  @seq_atom=$self->getseq(chain => $c, three => 1);
  $i=0;
  $j=0;
  for($i=0;$i<@seq3;$i++) {
    if($j < @seq_atom && $seq3[$i] eq $seq_atom[$j]) {
      $missing[$i]=0;
      $j++;
    } else {
      $missing[$i]=1;
    }
  }
      
  $seq="";
  for($i=0;$i<@seq3;$i++) {
    if(defined($args{missing}) && $missing[$i]) {
      $seq.=$args{missing};
    } else {
      $seq.=$aa1{$seq3[$i]};
    }
  }
  if(defined($args{three})) {
    return @seq3;
  } else {
    return $seq;
  }
}

sub AUTOLOAD {
  my $self=shift;
  if(!ref($self)) {
    croak "$self not an object";
  }
  my $name=our $AUTOLOAD;
  if($name =~ /::DESTROY$/) {
    return;
  }
  $name =~ s/^PDB:://g;
  if(!defined($self->{$name})) {
    croak "Can't access $name field in $self";
  }
  if(@_) {
    return $self->{$name}=shift;
  } else {
    return $self->{$name};
  }
}

1;
