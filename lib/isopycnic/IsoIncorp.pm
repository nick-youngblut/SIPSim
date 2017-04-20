package isopycnic::IsoIncorp;

use 5.006;
use strict;
use warnings;

=head1 NAME

isopycnic::IsoIncorp - scripts for running isopycnic::IsoIncorp.pl

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
use FindBin qw /$Bin/;
use IPC::Cmd qw/can_run/;
use List::Util qw/min max/;

use isopycnic::t;
our $argv_err;


=head2 call_rTraitContW

Calling R script to simulate trait evolution.

=head3 IN

-file :  tree file name
-format :  tree format
-weight :  weight param ($)

=head3 OUT

{taxon => trait_value}

=cut

push @EXPORT_OK, 'call_rTraitContW';

sub call_rTraitContW{
  my %h = @_;
  my $tree = exists $h{-tree} ? $h{-tree} : confess "Provide -tree";
  my $format = exists $h{-format} ? $h{-format} : confess "Provide -format";
  my $weight = exists $h{-weight} ? $h{-weight} : confess "Provide -weight";

  my $exe;
  if( can_run( "$Bin/rTraitContW.r" ) ){
    $exe = "$Bin/rTraitContW.r";
  }
  elsif( can_run( 'rTraitContW.r') ){
    $exe = 'rTraitContW.r';
  }
  else{
    confess "ERROR: cannot find rTraitContW.r\n";
  }

  my $cmd = "Rscript $exe -t $tree -f $format -w $weight | ";

  open PIPE, $cmd or confess $!;

  my %trait;
  while(<PIPE>){
    chomp;
    next if /^\s/;  # skipping header

    my @l = split / +(?=[0-9-])/, $_, 2;
    $trait{$l[0]} = $l[1];
  }
  close PIPE or confess $!;

  #print Dumper %trait; exit;
  return \%trait;
}


=head2 setIncorpDist

Setting % isotope incorporoation 
min-max & percent of taxa

=head3 IN

$trait_r :   {taxon => trait_value}
$argv_r :  {-min => min, -max => max, -percent => percent, ...}

=head3 OUT

$trait_r editted 

=cut

push @EXPORT_OK, 'setIncorpDist';

sub setIncorpDist{
  my $trait_r = shift || confess "Provide \$trait_r\n";
  my $argv_r = shift || confess "Provide \%ARGV\n";

  # IO check
  my $percent = exists $argv_r->{-percent} ?
    $argv_r->{-percent} : confess "Provide '-percent'";
  exists $argv_r->{-range} or confess "Provide -range min-max";
  confess "ERROR: percent must be between 0 and 100\n"
    unless $argv_r->{-percent} >= 0 and $argv_r->{-percent} <= 100;
  map {confess "ERROR: $_ must be between 0 and 100\n"
         unless $argv_r->{-range}{$_} >= 0
           and $argv_r->{-range}{$_} <= 100 } qw/min max/;


  # precentile (how many taxa should have incorporation?)
  ## ranking taxa by percent incorp
  my %ranks; # rank as percent of total taxa
  my $rank = 0;
  foreach my $taxon (sort{ $trait_r->{$b} <=> $trait_r->{$a} }
                     keys %$trait_r){
    $rank++;
    $ranks{$taxon} = $rank / (scalar keys %$trait_r);
  }


  ## dropping all taxa < $rank to incorp of 0
#  foreach my $taxon (keys %ranks){
#    $trait_r->{$taxon} = 0 if $ranks{$taxon} * 100 > $percent;
#  }
  # parsing taxa by incorp vs no incorp
  my %parsed;
  foreach my $taxon (keys %ranks){
    if($ranks{$taxon} * 100 > $percent){
      $parsed{noIncorp}{$taxon} = $trait_r->{$taxon};
    }
    else{
      $parsed{Incorp}{$taxon} = $trait_r->{$taxon};
    }
  }


  ## sanity check
  my $cnt = 0;
#  map{ $cnt++ if $parsed{Incorp}->{$_} > 0 } keys %$trait_r;
#  my $perc = $cnt / (scalar keys %$trait_r) * 100;
  my $perc_taxa_incorp = (scalar keys %{$parsed{Incorp}}) / (scalar keys %$trait_r) * 100;
  confess "ERROR; max percent taxa with incorporation ($percent)" .
    " exceeded ($perc_taxa_incorp)\n" if $perc_taxa_incorp > $percent;

  ## changing to percentiles of value for incorp taxa
  my %traits_perc;
  if( exists $parsed{Incorp} ){
    my $max_val = max values %{$parsed{Incorp}};
    map{ $_ = $_ / $max_val * 100 } values %{$parsed{Incorp}}
      unless $max_val <= 0;

    # rescale to min-max range
    sub rescale{
      # function(x) = (( max_scale - min_scale )( x - min_value)) / (max_value - min_value) + min_scale
      # return 0 if x = 0  # no incorporation
      # x = value in range of values
      # min_val & max_val = range of values
      # min_scale & max_scale = desired range of values
      my ($x, $min_val, $max_val, $min_scale, $max_scale) = @_;
      map{ defined $_ or confess "Value missing for rescale" } @_[0..4];
      if($x == 0){ return $x; }
      else{
	return (( $max_scale - $min_scale )*( $x - $min_val)) /
	  ($max_val - $min_val) + $min_scale;
      }
    }
    
    ## rescaling 
    my $min_val = min values %{$parsed{Incorp}};
    $max_val = max values %{$parsed{Incorp}};
    foreach my $k ( keys %{$parsed{Incorp}} ){
      $traits_perc{$k} = rescale( $parsed{Incorp}{$k}, $min_val, $max_val,
			       $argv_r->{-range}{min},
			       $argv_r->{-range}{max});
    }
  }

  # adding no-incorp taxa
  map{ $traits_perc{$_} = 0 } keys %{$parsed{noIncorp}};



#print Dumper %traits_perc; exit;
return \%traits_perc;
}


=head2 writeIncorpFile

Write incorportion file.

=head3 IN

$trait_r -- trait hashref
$argv_r -- getopt::Euclid argv

=head3 OUT

=cut

push @EXPORT_OK, 'writeIncorpFile';

sub writeIncorpFile{
  my $trait_r = shift || confess "Provide triats\n";
  my $argv_r = shift || confess $argv_err;

  # input check
  map{ exists $argv_r->{$_} || confess "Provide arg: '$_'\n" }
    qw/ -replicates /;

  # header
  my @header;
  ## control
  unless ( exists $argv_r->{-no_control} ){
    for my $i (1..$argv_r->{-replicates}){
      push @header, "Control_rep$i";   # 0% incorp
    }
  }
  ## incorp
  for my $i (1..$argv_r->{-replicates}){
    push @header, "Treatment_rep$i";
  }
  print join("\t", "#Taxon", @header), "\n";

  # body
  foreach my $taxon (keys %$trait_r){
    my @row;

    # control
    unless ( exists $argv_r->{-no_control} ){
      for my $i (1..$argv_r->{-replicates}){
	push @row, 0;   # 0% incorp
      }
    }

    # incorp
    for my $i (1..$argv_r->{-replicates}){
      push @row, $trait_r->{$taxon};
    }

    # writting row
    print join("\t", $taxon, @row), "\n";
  }

}


=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-isopycnic::IsoIncorp at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=isopycnic::IsoIncorp>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc isopycnic::IsoIncorp


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=isopycnic::IsoIncorp>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/isopycnic::IsoIncorp>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/isopycnic::IsoIncorp>

=item * Search CPAN

L<http://search.cpan.org/dist/isopycnic::IsoIncorp/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of isopycnic::IsoIncorp
