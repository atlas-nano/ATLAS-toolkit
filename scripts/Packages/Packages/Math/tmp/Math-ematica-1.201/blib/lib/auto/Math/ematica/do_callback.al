# NOTE: Derived from blib/lib/Math/ematica.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Math::ematica;

#line 434 "blib/lib/Math/ematica.pm (autosplit into blib/lib/auto/Math/ematica/do_callback.al)"
sub do_callback {
  my $link = shift;

  my $func_no = $link->read_packet;
  my $args = $link->read_packet;
  if ($FTABLE[$func_no]) {
    my @result;
    if (ref $args eq 'ARRAY') {
      @result = $FTABLE[$func_no]->(@$args);
    } else {
      @result = $FTABLE[$func_no]->();
    }
    $link->send_packet(@result);
  }
}

# end of Math::ematica::do_callback
1;
