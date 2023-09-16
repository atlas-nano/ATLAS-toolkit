# NOTE: Derived from blib/lib/Math/ematica.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Math::ematica;

#line 361 "blib/lib/Math/ematica.pm (autosplit into blib/lib/auto/Math/ematica/send_packet.al)"
sub send_packet {
  my $link = shift;

  $link->_send_packet(@_);
  $link->EndPacket;
  $link->Flush;
}

# end of Math::ematica::send_packet
1;
