# NOTE: Derived from blib/lib/Math/ematica.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Math::ematica;

#line 512 "blib/lib/Math/ematica.pm (autosplit into blib/lib/auto/Math/ematica/install.al)"
sub install {
  my ($link, $name, $nargs, $alias) = @_;
  my $package   = caller;
  $alias ||= $name;

  # This is very bad style. We steal the C pointer from $link since
  # DESTROY will not be called for it unless the function we generate
  # is undefined. Perl would die horribly when Perl_destruct would
  # encouter a blessed reference in the padlist of the function.
  # So *never* call the installed function after dropping $link!!!

  my $ptr = $link->{'mlink'};
  my $func = sub {
    my $link = bless {mlink => $ptr}, 'Math::ematica'; # this is the *nono*!
    my $result;
    
    if (defined $nargs) {
      die "${package}::$alias must be called with $nargs arguments\n"
        if $nargs != @_;
        $result    = $link->call([symbol($name), @_]);
    } else {
      die "${package}::$alias must be called with $nargs arguments\n"
        if @_;
        $result    = $link->call(symbol($name));
    }
    $link->{mlink} = 0;         # make DESTROY less harmfull
    $result;
  };

  no strict 'refs';
  *{"${package}::$alias"} = $func;
}

1;
# end of Math::ematica::install
