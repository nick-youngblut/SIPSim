package isopycnic::seq;

use 5.006;
use strict;
use warnings;

=head1 NAME

isopycnic::seq - scripts for running isopycnic::seq.pl

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
use Text::ParseWords;


=head2 parse_desc

parsing sequence description of reads from grinder

=head3 input

* description from Bio::DB::Fasta seq object <string>
* shotgun read? [false] <bool>
* read info \%

=head3 output

editted read_info

=cut

sub parse_desc{
  my ($desc, $shot_bool, $read_info_r) = @_;

  $desc =~ s/^>//;
  my %vals = quotewords("=| ", 0, $desc);

  # checking read descriptions
  my @chk = ('reference');
  if(! $shot_bool){ push @chk, 'amplicon'; }
  else{ push @chk, 'position'; }
  map{ die "ERROR: cannot find '$_' in read!\n"
         unless exists $vals{$_} } @chk;

  # strand
  my $strand;
  if($vals{$chk[1]} =~ /complement/){
    $strand = -1;
  }
  else{
    $strand = 1;
  }
  # postion
  $vals{$chk[1]} =~ /(\d+)\.\.(\d+)/;

  # loading hash
  $read_info_r->{ref} = $vals{reference};
  $read_info_r->{read_start} = $1;
  $read_info_r->{read_end} = $2;
  $read_info_r->{strand} = $strand;
}


=head2 get_read_info

Getting info on a read from a Bio::DB::Fasta object

=head3 IN

=head3 OUT

read_info editted

=cut

push @EXPORT_OK, 'get_read_info';

sub get_read_info{
  my $UID = shift || confess "Provide a read UID";
  my $in_r = shift || confess "Provide input object";
  my $argv_r = shift || confess "Provide GetOpt::Euclid \%ARGV";
  my $read_info_r = shift || confess "Provide \%read_info";

  # getting read seq object
  my $seqo = $in_r->{read_db}->get_Seq_by_id($UID);
  
  ## getting read info
  parse_desc($seqo->desc, $argv_r->{shotgun}, $read_info_r);
  $read_info_r->{read_len} = abs($read_info_r->{read_end} - 
				 $read_info_r->{read_start}) + 1;
  $read_info_r->{read_seq} = $seqo->seq;


  # genome_name
  # unless(exists $read_info_r->{genome_name}){
  #   print Dumper $in_r->{read_db}->header($UID); exit;

  #   my @l = grep(/description=/, quotewords('\s+', 0, $in_r->{read_db}->header($UID)) );
  #   croak("ERROR: no 'description=' found for read '$UID'\n")
  #     unless @l;
  #   map{ $read_info_r->{genome_name} = $1 if /description="*(.+)"*/ } $l[0];
  # }

}

=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-isopycnic::seq at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=isopycnic::seq>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc isopycnic::seq


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=isopycnic::seq>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/isopycnic::seq>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/isopycnic::seq>

=item * Search CPAN

L<http://search.cpan.org/dist/isopycnic::seq/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of isopycnic::seq
