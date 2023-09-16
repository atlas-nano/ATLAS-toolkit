# NOTE: Derived from blib/lib/Math/ematica.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Math::ematica;

#line 496 "blib/lib/Math/ematica.pm (autosplit into blib/lib/auto/Math/ematica/main.al)"
sub main {
  my $link = shift;

  $link->PutSymbol('End');
  $link->Flush;
  delete $link->{passive};
  $link->NewPacket;
  while (my $packet = $link->NextPacket) {
    if ($packet == CALLPKT) {
      $link->do_callback;
    } else {
      warn "Ignoring Unkown packet: $packet\n";
    }
  }
}

# end of Math::ematica::main
1;
