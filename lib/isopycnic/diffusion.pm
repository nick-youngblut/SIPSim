package isopycnic::diffusion;

use 5.006;
use strict;
use warnings;

=head1 NAME

isopycnic::diffusion - scripts for running isopycnic::diffusion.pl

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

=head1 EXPORT


=head1 SUBROUTINES/METHODS

=cut

use base 'Exporter';
our @EXPORT_OK = '';
use Carp qw/ confess carp croak /;
use Data::Dumper;
use List::Util qw/ sum /;
use Math::Random qw/
		     random_normal
		   /;


=head2 diffusion_browian

Simulating browian motion diffusion on a fragment.


=head3 IN

BD -- string with BD value
argv_diff -- '-diffusion' arg of getopt::euclid ARGV. 'n' and 'stdev'

=head3 OUT

=cut

push @EXPORT_OK, 'diffusion_brownian';
sub diffusion_brownian{
  my $BD = shift;
  defined $BD or confess "Provide BD value\n";
  my $argv_diff = shift or confess "Provide hashref with 'n' and 'stdev' keys\n";

  # input check
  map{ confess "ERROR: provide $_\n" 
	 unless defined $argv_diff->{$_} }
    qw/ n stdev /;

#  my $BD_orig = $BD;
    
  # getting all moves & summing 
  $BD += sum( random_normal($argv_diff->{n}, 0, $argv_diff->{stdev}) )
    unless $argv_diff->{n} < 1;

#  print $BD_orig - $BD, "\n";

  # floor
  $BD = 0 if $BD < 0;

  # return
  return $BD;
}





=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-isopycnic::diffusion at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=isopycnic::diffusion>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc isopycnic::diffusion


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=isopycnic::diffusion>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/isopycnic::diffusion>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/isopycnic::diffusion>

=item * Search CPAN

L<http://search.cpan.org/dist/isopycnic::diffusion/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of isopycnic::diffusion
