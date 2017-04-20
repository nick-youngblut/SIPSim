package isopycnic::PCR;

use 5.006;
use strict;
use warnings;

=head1 NAME

isopycnic::PCR - scripts for running isopycnic::PCR.pl

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
use Set::IntervalTree;
use Bio::Tools::AmpliconSearch;

=head2 makeAmpliconRegionIndex


=head3 IN

=head3 OUT

=cut

push @EXPORT_OK, 'makeAmpliconRegionIndex';

sub makeAmpliconRegionIndex{
  my %h = @_;

  # IO check
  map{ exists $h{$_} || confess "Provide $_\n" }
    qw/ -argv -genome_db /;

  # status
  print STDERR "Making amplicon region index...\n"
    unless $h{-argv}{'--quiet'};

  # Initialize search for amplicons
  my $amp_search = Bio::Tools::AmpliconSearch->new(
						   -primer_file => $h{-argv}{-forward_reverse},
						  );

  # all genome IDs
  my @genomes = $h{-genome_db}->ids;

  # foreach genome scaf/chromo
  my %amp_region_index;
  foreach my $genome_scaf_id (@genomes){
    ## stripping off \|SCAFID=\| from genome ID (to associate with af and if files)
    #(my $genome_id = $genome_scaf_id) =~ s/\|SCAFID=\d+\|//;

    # bio::seq object of genome scaf/chromo
    my $genomeO = $h{-genome_db}->get_Seq_by_id($genome_scaf_id);

    # amplicon search
    $amp_search->template($genomeO);

    while (my $amp = $amp_search->next_amplicon) {
      push @{$amp_region_index{$genome_scaf_id}}, [$amp->start, $amp->end];
    }    
  }

  #print Dumper %amp_region_index; exit;
  return \%amp_region_index;
}


=head2 fragsWithAmp

Screening out all fragments that to do not contain
the amplicon of interest.

Loading all fragment start-stops into an interval tree.
Then screening all amplicon regions

=head3 IN

fragStarts_r -- $@ of starts
fragEnds_r -- $@ of ends
kwargs

=head3 OUT

=cut

push @EXPORT_OK, 'fragsWithAmp';
sub fragsWithAmp{
  my $fragStarts_r = shift or confess "Provide frag starts\n";
  my $fragEnds_r = shift or confess "Provide frag ends\n";
  my %h = @_;

  # kwarg check
  map{ exists $h{$_} or confess "ERROR: Provide $_\n" }
    qw/ -input -genome_scaf_id /;
  my $amp_region_index = exists $h{-input}{-amp_region_index} ?
    $h{-input}{-amp_region_index} : 
      confess "ERROR: Provide -amp_region_index\n";
  
  # loading interval tree
  my $itree = Set::IntervalTree->new();
  for my $i (0..$#$fragStarts_r){
    $itree->insert([$fragStarts_r->[$i], $fragEnds_r->[$i]],
		  $fragStarts_r->[$i], $fragEnds_r->[$i]);
  }

  # fetching fragments that amplify
  my (@fragStarts, @fragEnds);
  if(exists $amp_region_index->{ $h{-genome_scaf_id} }){
    foreach my $amp_region ( @{$amp_region_index->{ $h{-genome_scaf_id} }} ){
      my $ret = $itree->fetch($amp_region->[0], $amp_region->[1]);
      
      foreach my $SE (@$ret){
	push @fragStarts, $SE->[0];
	push @fragEnds, $SE->[1];
      }
    }
  }
  else{  # for genome, no amplification with primers
    printf STDERR "WARNING: the provided primers did not target any regions in '%s'. No fragments kept!\n",
      $h{-genome_scaf_id};
  }

  # status
  unless( $h{-argv}{'--quiet'} ){
    printf STDERR "Amplicon region filtering for genome: '%s':\n", $h{-genome_scaf_id};
    printf STDERR "     Number of fragments pre-filter:  %i\n", scalar @$fragStarts_r;
    printf STDERR "     Number of fragments post-filter: %i\n", scalar @fragStarts;
  }


#  print Dumper @fragStarts[0..9], @fragEnds[0..9]; exit;
  return \@fragStarts, \@fragEnds;
}



=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS





=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc isopycnic::PCR


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=isopycnic::PCR>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/isopycnic::PCR>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/isopycnic::PCR>

=item * Search CPAN

L<http://search.cpan.org/dist/isopycnic::PCR/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of isopycnic::PCR
