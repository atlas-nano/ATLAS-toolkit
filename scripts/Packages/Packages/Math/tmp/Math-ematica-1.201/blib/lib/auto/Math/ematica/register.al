# NOTE: Derived from blib/lib/Math/ematica.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Math::ematica;

#line 408 "blib/lib/Math/ematica.pm (autosplit into blib/lib/auto/Math/ematica/register.al)"
sub register {
  my ($link, $name, $code, @args) = @_;
  my $fno = @FTABLE;
  my @parm;
  my $var = 'aaaa';

  push @FTABLE, $code;
  $FTABLE[$fno] = $code;
  my @list = (symbol 'List');

  for my $type (@args) {
    if (defined $type) {
      push @parm, [symbol 'Pattern', symbol $var, [ symbol 'Blank', symbol $type ]];
    } else {
      push @parm, [symbol 'Pattern', symbol $var, [ symbol 'Blank']];
    }
    push @list, symbol $var;
    $var++;
  }
  $link->call([symbol 'SetDelayed',
               [symbol $name, @parm],
               [symbol 'ExternalCall',
                [symbol 'LinkObject', "ParentLink", 1, 1],
                [symbol 'CallPacket', $fno, \@list]]]);
}

# end of Math::ematica::register
1;
