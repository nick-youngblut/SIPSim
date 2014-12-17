package isopycnic::GC;

use 5.006;
use strict;
use warnings;

=head1 NAME

isopycnic::GC - scripts for running isopycnic::GC.pl

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

use isopycnic::t;
our $genome_db_err;
our $argv_err;


=head2 sliding_window_GC

Sliding window GC analysis across
all genomes. min-max GC of all fragments
are recorded.

=head3 IN

$genome_db -- bio::db::fasta object
$argv_r -- GetOpt::Euclid \%ARGV

=head3 OUT

=cut

push @EXPORT_OK, 'sliding_window_GC';

sub sliding_window_GC{
  my $genome_db = shift || confess $genome_db_err;
  my $argv_r = shift || confess $argv_err;

  # TODO: add back gap fraction

  # window & jump size based on min frag size
  my $min_frag_len = exists $argv_r->{-range}{min_length} ?
    $argv_r->{-range}{min_length} : 
      confess "Cannot find -range => min_length";
  my $window = $min_frag_len / 2;
  my $jump = int($window / 2);

  # staus
  print STDERR "Calculating min-max GC of any potential genome fragment\n";

  # sliding window
  my %GCrange = (min => undef, max => undef);
  foreach my $id ($genome_db->ids){
    print STDERR "genome: $id\n";
   
    my $genome_seq = $genome_db->get_Seq_by_id($id);
    my $genome_len = $genome_seq->length;

    for (my $i=1; $i<=$genome_len; $i+=$jump){     
      my ($GC, $frag_len)  = calc_GC( $genome_seq->subseq( $i, $i+$window-1), 75 );      

      $GCrange{min} = $GC unless defined $GCrange{min};
      $GCrange{max} = $GC unless defined $GCrange{max};      
      $GCrange{min} = $GC if $GC < $GCrange{min};
      $GCrange{max} = $GC if $GC > $GCrange{max};
    }
  }

  print Dumper %GCrange;
}




=head2 get_frag_GC

Fragmenting genomes 

Creating fragments for each read & calculating GC on each.
Return (gather): a hash that replicates input hash along with frag_GC, frag_len, & frag_start

=cut

push @EXPORT_OK, 'get_frag_GC';

sub get_frag_GC{
  my ($genome, $in_r, $argv_r) = @_;
   
  # IO check
  map{ confess "ERROR: required arg '$_' is missing\n" 
	 unless exists $in_r->{$_} } qw/genome_db read_db 
					read_ids incorp/;
  map{ confess "ERROR: required arg '$_' is missing\n"
	 unless exists $argv_r->{$_} } qw/ -gap_fraction 
					   -isotope/;

  # isotope incorporation for the taxon
  my $perc_incorp = exists $in_r->{incorp}{ $genome } ?
    $in_r->{incorp}{ $genome } : 
      confess "Cannot find genome '$genome' in incorp file\n";

  # read_IDs for that genome
  my @spec_ids = grep /$genome\__/, @{$in_r->{read_ids}};
  croak "ERROR: no read IDs found for genome '$genome'\n" 
    unless @spec_ids;

  # genome_info
  my $genome_len = $in_r->{genome_db}->length($genome);
  die "ERROR: cannot find '$genome' in genomes database\n"
    unless defined $genome_len;

  
  # processing each read
  my %ret;  # \%% of info
  my %read_info;
  foreach my $UID (@spec_ids){

    # getting read info  (loaded into \%read_info)
    get_read_info( $UID, $in_r, $argv_r, \%read_info );

    # getting fragment size
    my $frag_size = getFragSize( $in_r, $argv_r, $genome);

    # sanity check: frag_size > read_len
    if($read_info{read_len} > $frag_size){
      printf STDERR "WARNING for %s -> %s: read_len:% > fragment_length:%i. Skipping!\n",
	$genome, $UID, $read_info{read_len}, $frag_size;
      next;
    }
    

    # determine fragment start-end based on read start-end
    my ($frag_start, $frag_end) = getFragStartEnd($frag_size, $argv_r, \%read_info);

    # sanity check: fragment start within range of read?
    printf STDERR  "WARNING for genome '%s': frag_start is too far from read!\n", $genome
	  if $frag_start - ($read_info{read_end} + $argv_r->{-primer_buffer}) > $frag_size;
    
    # wrapping fragment (if circular genome)
    my $pos_r = wrapFrag($frag_start, $frag_size, $genome_len);

    # getting sequence
    my $frag_seq = "";
    foreach (keys %$pos_r){
      $frag_seq .= $in_r->{genome_db}->seq($genome, 
					   ${$pos_r->{$_}}[0], 
					   ${$pos_r->{$_}}[1]);
    }

    # sanity check: fragment is encompassing read
    fragEncompassRead($frag_start, $frag_seq, $genome_len,
		      $read_info{read_start},
		      $read_info{read_end}, $genome);
			
    
    # calculating fragment GC (finally!)
    my ($frag_GC, $frag_length) = calc_GC( $frag_seq, $argv_r->{-gap_fraction} );

    # calculating fragment BD
    ## no isotope 
    my($frag_BD, $frag_BD_iso) = calc_BD($frag_GC, $perc_incorp, $argv_r->{-isotope});

    # calculating read GC 
    confess "Cannot find read_info{read_seq}" unless exists $read_info{read_seq};
    my ($read_GC, $read_length) = calc_GC( $read_info{read_seq}, $argv_r->{-gap_fraction} );

    # calculating read BD
    ## no isotope 
    my($read_BD, $read_BD_iso) = calc_BD($read_GC, $perc_incorp, $argv_r->{-isotope});


    # loading hash of GC and BD info
    # UID => cat => value
    $ret{$UID} = { 
		  frag_GC => $frag_GC,
		  frag_BD => $frag_BD,
		  frag_BD_iso => $frag_BD_iso,
		  read_GC => $read_GC,
		  read_BD => $read_BD,
		  read_BD_iso => $read_BD_iso,
		  read_start => $read_info{read_start},
		  read_len => $read_info{read_len},
		  frag_start => $frag_start,
		  frag_len => $frag_size
		 }; 	
  }

  return \%ret;
}


=head2 calc_GC

Calculating GC content for a string (DNA)

=head3 IN

1) dna_string 
2) fraction of string length that can be gaps ([-.]+); otherwise 'NA'
 (value range: 0-1)

=head3 OUT

$GC_content
$len :  length of sequence (no gaps)

=cut

push @EXPORT_OK, 'calc_GC';

sub calc_GC {
# calculting GC of a sequence
  my $seq = shift || confess "ERROR: no sequence provided\n";
  my $gap_frac = shift;
  $gap_frac = 0.05 unless defined $gap_frac;

  # IOcheck
  confess "gap_frac range: 0-1\n" unless $gap_frac >= 0 and $gap_frac <= 1;

  # length
  my $raw_len = length $seq;

  # removing gaps; checking to see if too many gaps
  $seq =~ s/[.-]+//g;
  my $noGap_len = length $seq;
  my $val = ($raw_len - $noGap_len) / $raw_len;
  return ('NA',$noGap_len) if ($raw_len - $noGap_len) / $raw_len > $gap_frac;

  # uppercase
  $seq =~ tr/a-z/A-Z/;

  # scoring table (included ambiguous nucleotides)
  my %score = (G => 1,
               C => 1,
               R => 0.5,
               Y => 0.5,
               S => 1,
               K => 0.5,
               M => 0.5,
               B => 0.66,
               D => 0.33,
               H => 0.33,
               V => 0.66,
               N => 0.5
               );

  # GC
  my $q = join("", "[", keys %score, "]");
  $q = qr/$q/;
  my $GC_sum = 0;
  $GC_sum += $score{$1} while $seq =~ /($q)/g;
  my $GC_content = $GC_sum / $noGap_len * 100;

  return $GC_content;
}

push @EXPORT_OK, 'calc_GC_fast';
sub calc_GC_fast {
# calculting GC of a sequence
  my $seq = shift || confess "ERROR: no sequence provided\n";
#  my $gap_frac = shift;
#  $gap_frac = 0.05 unless defined $gap_frac;

  # IOcheck
#  confess "gap_frac range: 0-1\n" unless $gap_frac >= 0 and $gap_frac <= 1;

  # removing gaps; checking to see if too many gaps
  $seq =~ s/[.-]+//g;

  # seq length
  my $seq_len = length $seq;

  # uppercase
  $seq =~ tr/a-z/A-Z/;

  # scoring table (included ambiguous nucleotides)
  # my %score = (G => 1,
  #              C => 1,
  #              S => 1,
  #              B => 0.66,
  #              V => 0.66,
  #              R => 0.5,
  #              Y => 0.5,
  #              K => 0.5,
  #              M => 0.5,
  #              N => 0.5,
  #              D => 0.33,
  #              H => 0.33
  #              );
  
  # summing each value category
  my $GCS = 0;
  $GCS++ while $seq =~ /[GCS]/g;
  my $BV = 0;
  $BV++ while $seq =~ /[BV]/g;
  my $RYKMN = 0;
  $RYKMN++ while $seq =~ /[RYKMN]/g;
  my $DH = 0;
  $DH++ while $seq =~ /[DH]/g;
 

  # GC
  my $GC_sum = $GCS * 1 + $BV * 0.66 + $RYKMN * 0.5 + $DH * 0.33;
 
  return $GC_sum / $seq_len * 100;
}



=head2 calc_BD

calculating fragment buoyant density
based on [GC]

=head3 IN

fragGC -- fragment GC value (percentage GC)
perc_incorp -- percent isotope incorporation. {taxon=>percent}

=head3 OUT

=cut

sub calc_BD{
  my $frag_GC = shift or confess;
  my $perc_incorp = shift or confess;
  my $isotope_r = shift or confess;


  # IO check
  return ('NA','NA') unless $frag_GC =~ /^[\d.]+$/;  # 'NA' unless fragGC real number

  # Buoyant density
  my $frag_BD = ($frag_GC * 0.098 / 100) + 1.660;

  # isotope incorporation
  ## isotope index (BD shift for 100% incorporation
  my %isoIndex = ( '13C' => 0.036,
		   '15N' => 0.016 );
  ## check isotopes exists
  map{ die "Isotope '$_' no supported\n"
	 unless exists $isoIndex{$_} } @$isotope_r;

  ## adding BD to frag_BD
  my $frag_BD_iso = $frag_BD;
  map{ $frag_BD_iso += $isoIndex{$_} * ($perc_incorp / 100) }
    @$isotope_r;

  return $frag_BD, $frag_BD_iso;
}



=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-isopycnic::GC at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=isopycnic::GC>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc isopycnic::GC


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=isopycnic::GC>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/isopycnic::GC>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/isopycnic::GC>

=item * Search CPAN

L<http://search.cpan.org/dist/isopycnic::GC/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of isopycnic::GC
