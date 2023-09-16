# NOTE: Derived from blib/lib/Math/ematica.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Math::ematica;

#line 369 "blib/lib/Math/ematica.pm (autosplit into blib/lib/auto/Math/ematica/_send_packet.al)"
# The following is a cludge. The goal ist to make the Perl syntax the
# same as the mathematica syntax execpt that the opening '[' are move
# one token left:
# Mathematica Perl
# Sin[x]      [Sin, x]
# Pi          Pi
# Blank[]     [Blank]

sub _send_packet {
  my $link  = shift;

  while (@_) {
    my $elem = shift;
    if (ref $elem eq 'ARRAY') {
      $link->_send_call(@$elem);
    } else {
      $link->PutToken($elem);   # PutSymbol in doubt
    }
  }
}

# end of Math::ematica::_send_packet
1;
