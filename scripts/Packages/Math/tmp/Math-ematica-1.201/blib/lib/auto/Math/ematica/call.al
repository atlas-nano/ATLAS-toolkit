# NOTE: Derived from blib/lib/Math/ematica.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Math::ematica;

#line 450 "blib/lib/Math/ematica.pm (autosplit into blib/lib/auto/Math/ematica/call.al)"
sub call {
  my $link = shift;
  my $fname = shift;

  # first argument may be symbol name instead of a symbol
  $fname = symbol $fname unless ref $fname;
  $link->send_packet($fname, @_);
  $link->dispatch unless $link->{passive};
}

# end of Math::ematica::call
1;
