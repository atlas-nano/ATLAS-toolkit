# NOTE: Derived from blib/lib/Math/ematica.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Math::ematica;

#line 460 "blib/lib/Math/ematica.pm (autosplit into blib/lib/auto/Math/ematica/dispatch.al)"
sub dispatch {
  my $link = shift;
  
  $link->NewPacket;
  while (my $packet = $link->NextPacket) {
    if ($packet == RETURNPKT) {
      return $link->read_packet;
    } elsif ($packet == MESSAGEPKT) {
      return $link->read_packet;
    } elsif ($packet == TEXTPKT) {
      return $link->read_packet;
    } elsif ($packet == CALLPKT) {
      $link->do_callback;
    } elsif ($packet == DISPLAYPKT) {
      $link->GetNext() == MLTKSTR  or die "Expected DISPLAYPKT to start with 'MLTKSTR'";
      $link->GetByteString() eq '' or die "Expected DISPLAYPKT to start with empty string";
      # $link->GetNext() == MLTKFUNC or die "Expected DISPLAYPKT to contain 'MLTKFUNC'";
      my $result = '';
      while ($link->GetNext() == MLTKFUNC) {
        my ($name, $nargs) = $link->GetFunction();
        $$name eq "DisplayPacket"  or
          $$name eq "DisplayEndPacket" or die "Expected 'DisplayPacket' symbol in DISPLAYPKT, not '$$name'";
        $nargs == 1                  or die "Expected 'DisplayPacket'to habe one argument only";
        $result .= $link->GetByteString();
        return $result if $$name eq "DisplayEndPacket";
      }
    } elsif ($packet == INPUTNAMEPKT) {
      next;
    } else {
      warn "Ignoring Unkown packet: $packet\n";
      return;
    }
  }
  $link->NewPacket;
}

# end of Math::ematica::dispatch
1;
