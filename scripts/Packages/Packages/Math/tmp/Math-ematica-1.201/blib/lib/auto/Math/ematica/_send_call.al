# NOTE: Derived from blib/lib/Math/ematica.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Math::ematica;

#line 390 "blib/lib/Math/ematica.pm (autosplit into blib/lib/auto/Math/ematica/_send_call.al)"
sub _send_call {
  my ($link, $head, @tail)  = @_;

  if (ref $head eq 'ARRAY') {
    $link->_send_call(@$head);
  } else {
    $link->PutToken($head, scalar @tail) # PutFunction in doubt
  }
  while (@tail) {
    my $elem = shift @tail;
    if (ref $elem eq 'ARRAY') {
      $link->_send_call(@$elem);
    } else {
      $link->PutToken($elem); # PutSymbol in doubt
    }
  }
}

# end of Math::ematica::_send_call
1;
