package isopycnic::t;

use 5.006;
use strict;
use warnings;
use POSIX;

=head1 NAME

isopycnic::t - sanity checks

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
use Bio::SeqIO;


#--- error handling variables ---#
our @EXPORT = ();
# GetOpt::Euclid
our $argv_err = 'Provide GetOpt::Euclid \%ARGV hashref\n';
our $genome_db_err = 'Provide a Bio::DB::Fasta object\n';
push @EXPORT, qw/$argv_err $genome_db_err/;



=head2 check_length_model

Checking to make sure the correct length model
was provided and that it has the correct
number of args or defaults provided

=cut

push @EXPORT_OK, 'check_length_model';

sub check_length_model{
  my $argv_r = shift || confess $argv_err;


  exists $argv_r->{-length_model} || 
    confess "-length_model arg not found\n";
  $argv_r->{-length_model} =~ s/ .+//;

  # index of defaults
  my %index = ( 
	       uniform => [],
	       normal => [12500, 2000],
	       'skewed-normal' => [10],
	       exponential => [12500],
	       poisson => [10],
	       f => [90, 10]
	      );

  my $lm = $argv_r->{-length_model};

  # setting defaults if not found
  my %new_lm;
  my $model = defined $lm->[0] ? 
    shift @$lm : confess "No length model provided\n";
  exists $index{ $model } || confess "Length model '$model' not supported\n";
  for my $i ( 0..$#{$index{$model}} ){    
    $new_lm{ $model }->[$i] = defined $lm->[$i] ?
      $lm->[$i] : $index{ $model }->[$i];
  }
  
  $argv_r->{-length_model} = \%new_lm;
  $argv_r->{-lm} = \%new_lm;

  #print Dumper $argv_r->{-length_model}; exit;
  return 1;
}


=head2 fragEncompassRead

Is the fragment encompassing the read?

=head3 IN

=head3 OUT

=cut

push @EXPORT_OK, 'fragEncompassRead';

sub fragEncompassRead{
  my $frag_start = shift || confess;
  my $frag_seq = shift || confess;
  my $genome_len = shift || confess;
  my $read_start = shift || confess;
  my $read_end = shift || confess;
  my $genome = shift || confess;

  # sanity check: fragment is encompassing read
  my ($frag_start_chk, $frag_end_chk) = ($frag_start, $frag_start + length($frag_seq));
  $frag_start_chk = 0 if $frag_start_chk < 0;
  $frag_end_chk = $genome_len if $frag_end_chk > $genome_len;
  if($frag_start_chk > $read_start or $frag_end_chk < $read_end){ 
    # fragment must encompass the amplicon
    warn "WARNING for genome $genome: fragment does not encompass the read (frag_start=$frag_start_chk, frag_end=$frag_end_chk, read_start=$read_start, read_end=$read_end). Skipping\n";
    return 0;
  }
  else{ return 1; } # success
}


=head2 check_label_match

Checking that tree labels match rank-file taxon
labels.

=head3 IN

$tree :  TreeIO object
$rank :  Rank file object

=head3 OUT

=cut

push @EXPORT_OK, 'check_label_match';

sub check_label_match{
  my $tree = shift || confess "Provide TreeIO object";
  my $ranks = shift || confess "Provide rank file object";

  # intersection of leaf & rank-taxonID
  my @leaf_ids = map{ $_->id } $tree->get_leaf_nodes;

  my (%union, %isect);
  foreach my $e (@leaf_ids, keys $ranks){
    $union{$e}++ && $isect{$e}++; 
  }

  die "ERROR: the leaf IDs do not match total match taxa in *-rank.txt file.
\tAre there commas or quotes in the leaf lables?
\tCommas and quotes are not allowed.\n"
    unless scalar keys %isect == scalar keys %$ranks;
}


=head2 correct_fasta

Making sure sequence lines for each entry are the same length.

=head IN

genome file

=head OUT

=cut

push @EXPORT_OK, 'correct_fasta';

sub correct_fasta{
  my $genome_file = shift || confess "Provide a fasta as $%";

  (my $genome_out = $genome_file) =~ s/\.[^.]+$/_fix.fasta/;

  my $seqin = Bio::SeqIO->new(-file => $genome_file, -format => 'fasta');
  my $seqout = Bio::SeqIO->new(-file => ">$genome_out", -format => 'fasta');

  while(my $seq = $seqin->next_seq){
    $seqout->write_seq($seq);
  }

  print STDERR "Corrected genome file written: '$genome_out'\n";
  return $genome_out;
}


=head2 check_abs_abund

Checking that the abs abund of each sample (af table)
matches the total genome copies.

=head3 IN

$af_r -- af table
$argv_r -- getopt::euclid

=cut

push @EXPORT_OK, 'check_abs_abund';

sub check_abs_abund{
  my $af_r = shift or confess "Provide af_r\n";
  my $argv_r = shift or confess "Provide argv_r\n";

  # input check
  exists $argv_r->{-total_copies} or confess "Provide '-total_copies'\n";

  # totaling by sample
  my %sample_sum;
  foreach my $genome_id (keys %$af_r){
    foreach my $sample (keys %{$af_r->{$genome_id}}){
      $sample_sum{$sample} += $af_r->{$genome_id}{$sample};
    }
  }


  # checking that grinder was approx the required total genome copies
  ## 95% accuracy  
  my $cutoff = 0.95;
  map{ $sample_sum{$_} / $argv_r->{-total_copies} >= $cutoff or 
	 die "\nGrinder ERROR: genome copies in sample '$_'" .
	   " differed by >" . sprintf("%.2f", (1-$cutoff) * 100) .
	   "% of -total_copies (sample copies: " .
	     ceil($sample_sum{$_}) . "; -total_copies: " . 
	       ceil($argv_r->{-total_copies}) .
		 "). Retry run.\n" }  keys %sample_sum;

}


=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-isopycnic::t at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=isopycnic::t>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc isopycnic::t


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=isopycnic::t>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/isopycnic::t>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/isopycnic::t>

=item * Search CPAN

L<http://search.cpan.org/dist/isopycnic::t/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of isopycnic::t
