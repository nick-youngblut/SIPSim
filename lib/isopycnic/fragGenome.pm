package isopycnic::fragGenome;

use 5.006;
use strict;
use warnings;

=head1 NAME

isopycnic::fragGenome - scripts for running isopycnic::fragGenome.pl

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
use List::Util qw/ min max sum/;
use List::MoreUtils qw/minmax/;
use Set::IntervalTree;
use Regexp::Common;
use MCE::Map;
use Math::Gauss ':all';
use Math::Random::SkewNormal qw/generate_sn/;
use Math::Random qw/
random_uniform_integer
random_uniform
random_normal
random_exponential
random_poisson
random_f/;

use isopycnic::GC qw/calc_GC calc_GC_fast/;
use isopycnic::IO qw/writeFrags/;
use isopycnic::PCR qw/fragsWithAmp/;
use isopycnic::fraction qw/countByFrac/;
use isopycnic::diffusion qw/diffusion_brownian/;


=head2 fragGenomes

fragment each genome to achive the absolute
copy number (coverage) needed for that genome.
For each fragment:
  Calculate its GC + BD
  Write it to the correct fraction file.

=head3 IN

argv -- GetOpt::Euclid ARGV
fractions -- fraction filehandles
input:
   sample -- sample ID (treatment/control)
   if -- parsed incorp file
   input -- contains af, if, header, genome_db

=head3 OUT

=cut

push @EXPORT_OK, 'fragGenomes';
sub fragGenomes{
  my %h = @_;

  # IO check
  map{ exists $h{$_} || confess "Provide $_\n" }
    qw/ -sample -argv -frac_ranges -input /;


  # genome ids & status
  my @genomes = $h{-input}{-genome_db}->ids;


  # foreach genome scaf/chromo
  for my $i (0..$#genomes){
    # genome ID
    my $genome_scaf_id = $genomes[$i];
    ## stripping off \|SCAFID=\| from genome ID (to associate with if file)
    (my $genome_id = $genome_scaf_id) =~ s/\|SCAFID=\d+\|//;

    # check to make sure genome_id in abundanc file & incorp file
    exists $h{-input}{-if}{$genome_id} ||
      confess "Could not find '$genome_id' in incorp file\n";

    # getting genome scaf/chromo length
    my $genome_len = $h{-input}{-genome_db}->length($genome_scaf_id);

    ## calculating coverage needed
    ## making the fragments required
    ## calculating GC + BD of each fragment
    makeFrags(
	      -genome_scaf_id => $genome_scaf_id, # with |SCAFID=\d+|
	      -genome_id => $genome_id,
	      -genome_len => $genome_len,
	      %h
	     );
  }

  return 1; # end of main function for fragGenome
}



=head2 makeFrags

    ## calculating coverage needed
    ## making the fragments required
   ## calculating GC + BD
    ## writing to file

=head3 IN

-genome_scaf_id
-genome_id -- genome id as string
%h -- see fragGenomes

=head3 OUT

[seq : sequence,
 start : frag_start, 
 end : frag_end,
 size : frag_size,
 GC : frag_GC,
 BD_noIso : BD with no isotope incorp
 BD_iso : BD with isotope incorp]

=cut

sub makeFrags{
  my %h = @_;
 
  ## IOcheck
  map{confess "Provide $_\n" unless exists $h{$_} }
    qw/ -genome_scaf_id -genome_id -argv 
	-input -sample -genome_len/;
 
  my $genome_scaf_id = $h{-genome_scaf_id};
  my $sample = $h{-sample};

  # getting all frag sizes needed for the coverage; returning hashref
  my $fragSizes_r = getFragSizes( %h ); 

  # getting a genome start-end positions for each fragment  
  my ($fragStarts_r, $fragEnds_r, $ampsPerFrag_r) = @_;
  if( $h{-argv}{-forward_reverse} ){
    ## fragments must contain amplicon region
    ($fragStarts_r, $fragEnds_r) = getAmpFragStartEnd( $fragSizes_r, %h );    
    
    ## return if no frag starts/ends (no amplicons found)
    return 0 unless defined $fragStarts_r and defined $fragEnds_r;

    ## determine number of complete amplicon templates per fragment
    $ampsPerFrag_r = getAmpTempPerFrag( $fragStarts_r, $fragEnds_r, %h );
  }
  else{
    ## fragments from anywhere in the genome
    ($fragStarts_r, $fragEnds_r) = getFragStartEnd( $fragSizes_r, %h );
  }


  # extracting each fragment 
  ## calc GC & BD
  ## binning by fraction
  my $allFragInfo_r = getGenomeFrags( 
				     $fragStarts_r, 
				     $fragEnds_r, 
				     $fragSizes_r,
				     %h 
				    );
 
  # writing
  foreach my $frag (@$allFragInfo_r){
    print join("\t", $sample, $genome_scaf_id, @$frag), "\n";
  }
 
  # counting frags by the gradient fractions they fall into
#  countByFrac( $allFragInfo_r, $ampsPerFrag_r, %h );

  return 1;
}


=head2 getFragSizes

Simulating fragment size encompassing 
the focal read. The size is drawn 
from a user-selected distribution.

=head3 IN

$coverage -- amout of coverage needed for genome (string)
%h -- see makeFrags

=head3 OUT

$@ -- fragsizes 
$ -- index to number of value in fragsizes actually needed

=cut

sub getFragSizes{
#  my $coverage = shift;
  my %h = @_;

  # IO check
  my $argv_r = exists $h{-argv} ? 
    $h{-argv} : confess "No -argv\n";  
  map{ exists $argv_r->{$_} or confess "Provide arg: $_" }
    qw/ -length_model -range -total_fragments/;
  map{ exists $argv_r->{-range}{$_} or confess "Provide -range -> $_\n" }
    qw/ min_length max_length /;


  # status
  print STDERR "Getting fragment sizes...\n" if $h{-argv}{'--debug'} > 1;

  # setting variables
  my $frag_min = $argv_r->{-range}{min_length};
  my $frag_max = $argv_r->{-range}{max_length};


  # frag sizes
  my $frag_sizes_r =  fragSizeFromDist( 
				       -N_frags => $h{-argv}{-total_fragments},  
				       -frag_min => $frag_min, 
				       -frag_max => $frag_max,
				       %h
				      );

  rescaleSizes($frag_sizes_r, $frag_min, $frag_max);
 

  # return
#  print Dumper $frag_sizes_r; exit;
  return $frag_sizes_r;
}


=head2 rescaleSizes

Rescaling frag sizes to defined min-max

=cut

sub rescaleSizes{
  my $frag_sizes_r = shift or confess "Provide frag sizes\n";
  my $frag_min = shift or confess "Provide frag_min\n";
  my $frag_max = shift or confess "Provide frag_max\n";
  
  # min-max value of 'raw' frag sizes (pulled from distribution)
  my ($min_val,$max_val) = minmax @$frag_sizes_r;

  # prevent 'divide by zero' error
  if($max_val - $min_val == 0){  
    $frag_sizes_r->[0] += 1;
    ($min_val,$max_val) = minmax @$frag_sizes_r;
  }

  # rescale each and 
  @$frag_sizes_r = map{
    int( 
	(( $frag_max - $frag_min )*( $_ - $min_val ))/
	($max_val - $min_val) + $frag_min
      )
  } @$frag_sizes_r;

}


=head2 fragSizeFromDist

subroutine for getFragSizes

=head3 IN

hash of args:
-N_frags -- number of frag sizes needed
-frag_min -- min frag length
-frag_max -- max frag length
%h -- args passed from getFragSizes

=cut

sub fragSizeFromDist{
  my %h = @_;

  # input check
  map{ exists $h{$_} or confess "Provdie arg: $_\n" }
    qw/ -N_frags -frag_min -frag_max -argv -input/;
  map{ exists $h{-argv}{$_} or confess "Provide -argv -> $_\n" }
    qw/ -nw -cs -lm /;

  my $lm_r = $h{-argv}{-lm};
  


  # fragment size
  my @frag_sizes;
  if( exists $lm_r->{uniform} ){	
    @frag_sizes = random_uniform_integer($h{-N_frags}, 
					$h{-frag_min}, 
					$h{-frag_max});
  }
  elsif( exists $lm_r->{normal} ){
    @frag_sizes = random_normal($h{-N_frags}, 
			       @{$lm_r->{normal}} );
  }
  elsif( exists $lm_r->{exponential} ){
    @frag_sizes = random_exponential($h{-N_frags}, 
				    @{$lm_r->{exponential}} );
    }
  elsif( exists $lm_r->{f} ){
    @frag_sizes = random_f( $h{-N_frags}, 
			   @{$lm_r->{f}} );  # scaling by 5; values >5 will be excluded
  }
  elsif( exists $lm_r->{poisson} ){
    @frag_sizes = random_poisson( $h{-N_frags}, 
				   @{$lm_r->{poisson}} );
    }
  elsif( exists $lm_r->{'skewed-normal'} ){
    my ($skewness) = @{$lm_r->{'skewed-normal'}};
    @frag_sizes = map{ generate_sn( -$skewness ) } 1..$h{-N_frags};
  }
  else{ # model not supported
    foreach my $lm (keys %$lm_r){
      confess "ERROR: do not recognize length model: $lm \n"; 
    }
  }


  #print Dumper 'frag_sizes', @frag_sizes; exit;
  return \@frag_sizes;
}


=head2 getFragStartEnd

Determining genome start-end position for each fragment.
Start chosen from a uniform distribution.
End = start + length -1

WARNING: end can be > genome_length

=head3 IN


=head3 OUT

=cut

sub getFragStartEnd{
  my $frag_sizes_r = shift or confess "Provide frag_sizes_r\n";
  my %h = @_;

  #Input check
  my $genome_len = exists $h{-genome_len} ?
    $h{-genome_len} : confess "Provide -genome_len\n";
  map{ exists $h{-argv}{$_} or confess "Provide -argv -> $_\n" }
    qw/ -nw -cs -lm /;
  

  # start for each fragment
  my @starts = random_uniform_integer( scalar @$frag_sizes_r,
				       1, $genome_len );

  # getting fragment ends
  my @ends;
  for my $i (0..$#$frag_sizes_r){
    push @ends, $starts[$i] + $frag_sizes_r->[$i] -1;
  }

#  print Dumper @starts[0..9], @ends[0..9]; exit;
  return \@starts, \@ends;
}


=head2 getAmpFragStartEnd

Determining genome start-end position for each fragment.
Start chosen from a uniform distribution.
End = start + length -1

WARNING: end can be > genome_length

=head3 IN


=head3 OUT

=cut

sub getAmpFragStartEnd{
  my $fragSizes_r = shift or confess "Provide fragSizes_r\n";
  my %h = @_;

  # primer buffer
  my $primer_buffer = 50;  # bp

  #Input check
  map{ exists $h{$_} or confess "Provide $_\n" }
    qw/ -genome_len -genome_scaf_id/;
  map{ exists $h{-input}{$_} or confess "Provide input -> $_\n"}
    qw/ -amp_region_index /;
  map{ exists $h{-argv}{$_} or confess "Provide -argv -> $_\n" }
    qw/ -nw -cs -lm /;
  my $amp_region_index = $h{-input}{-amp_region_index};
  my $genome_scaf_id = $h{-genome_scaf_id};

  # no amplicons for scaffold?
  unless( exists $amp_region_index->{ $genome_scaf_id } ){
    warn "WARNING: Cannot find $genome_scaf_id in amplicon region index.\n";
    warn "         This was likely caused by no amplicons generated from provided primers.\n";
    return undef, undef;
  }


  # making random index for selecting amplicons
  my @rand_index = random_uniform_integer(
					  scalar @$fragSizes_r,
					  0, 
					  $#{$amp_region_index->{$genome_scaf_id}},					  
					 );
  
  # foreach fragment, select start-end 
  my (@fragStarts, @fragEnds);
  for my $i (0..$#rand_index){
    my $frag_size = defined $fragSizes_r->[$i] ?
      $fragSizes_r->[$i] : confess "Internal error\n";

    # start-end of the amplicon
    my $SE_r = defined $amp_region_index->{$h{-genome_scaf_id}}->[$rand_index[$i]] ?
      $amp_region_index->{$h{-genome_scaf_id}}->[$rand_index[$i]] :
	confess "Undefined index: ".$rand_index[$i]."\n";
    scalar @$SE_r == 2 or confess "Internal error\n";

    # amplicon center postion
    my $amp_len = $SE_r->[1] - $SE_r->[0] + 1;
    my $amp_center = $SE_r->[0] + $amp_len / 2;

    ## fragment_start = amp_center - (frag_size * x);
    ## x = random draw from unifrom distribution (range: 0-1)
    my $x = random_uniform();

    ### frag_start_floor = (primer_buffer + 0.5*amp_len) / frag_size
    my $floor  = ($primer_buffer + 0.5 * $amp_len) / $frag_size;
    $x = $floor if $x < $floor;

    ### frag_start_ceiling = 1 - floor
    my $ceiling = 1 - $floor;                               # ceiling = 1 - floor
    $x = $ceiling if $x > $ceiling;
    my $frag_start = int( $amp_center - ($frag_size * $x) );

    # frag end
    my $frag_end = $frag_start + $frag_size -1;
    
    # sanity check: fragment encompasses amplicon
    ($frag_start <= $SE_r->[0] and $frag_end >= $SE_r->[1]) or
      confess "ERROR: fragment does not fully encompass the amplicon template\n";
    

    # loading fragment start-end
    push @fragStarts, $frag_start;
    push @fragEnds, $frag_end;   
  }

  #print Dumper @fragStarts[0..9], @fragEnds[0..9]; exit;
  return \@fragStarts, \@fragEnds;
}

=head2 getAmpTempPerFrag

Getting the number of amplicon templates per fragment.

=head3 IN

fragStarts_r -- $@ of starts
fragEnds_r -- $@ of ends
kwargs

=head3 OUT

=cut

sub getAmpTempPerFrag{
  my $fragStarts_r = shift or confess "Provide fragStarts\n";
  my $fragEnds_r = shift or confess "Provide fragEnds\n";
  my %h = @_;

  # input check
  my $amp_region_index = $h{-input}{-amp_region_index};
  my $genome_scaf_id = $h{-genome_scaf_id};
  exists $amp_region_index->{ $genome_scaf_id } or
    confess "ERROR: Cannot find $genome_scaf_id in amplicon region index\n";

  # counting amps in fragment for each frag start-stop
  my @amp_cnts;
  for my $i (0..$#$fragStarts_r){
    $amp_cnts[$i] = 0 unless defined $amp_cnts[$i];
    # each amplicon position in the genome
    foreach my $amp_SE_r ( @{$amp_region_index->{$genome_scaf_id}} ){
      $amp_cnts[$i]++ if 
	$fragStarts_r->[$i] <= $amp_SE_r->[0] and
	  $fragEnds_r->[$i] >= $amp_SE_r->[1];	

    }
  }

 # print Dumper @amp_cnts; exit;
  return \@amp_cnts;
}


=head2 processGenomeFrags

Extracting fragment sequences from the genome.
Foreach fragment:
 * calculating GC
 * calculating BD

Fragment will be binned by fraction.
If only amplicon fragments:
 * fraction abundance = frags_in_frac / total_frags_made = relative amplicon template copy number
ElseIf all genome fragments:
 * fraction abundance = length_frags_in_frac / total_length_frags = relative genome copy number


IF frag_end > genome_len:
    need to wrap frag (end of frag at start of linear genome)

=head3 IN

=head3 OUT

=cut

sub processGenomeFrags{
  my $fragStarts_r = shift or confess "Provide fragStarts_r\n";
  my $fragEnds_r = shift or confess "Provide fragEnds_r\n";
  my $fragSizes_r = shift or confess "Provide fragSizes_r\n";
  my %h = @_;


  # IOcheck
  map{ exists $h{$_} or confess "Provide arg: $_\n" }
    qw/ -genome_scaf_id -genome_id -genome_len -frac_ranges/;
  exists $h{-input}{-genome_db} or confess "Provide -input -> -genome_db\n";
  my $genome_len = exists $h{-genome_len} ?
    $h{-genome_len} : confess "ERROR: provide -genome_len\n";

  ## MCE/fork variables
  map{ exists $h{-argv}{$_} or confess "Provide -argv -> $_\n" }
    qw/ -nw -cs /;

  
  # gettig bio::db::fasta object of genomes
  my $seqo = $h{-input}{-genome_db}->get_Seq_by_id( $h{-genome_scaf_id} );

  # getting fragment sequences
  my @seqs;
  foreach my $i (0..$#$fragStarts_r){
    # getting genome fragment
    ## frag may 'wrap' around linear genome. Need to account for this
    my $frag_seq = "";
    if( $fragEnds_r->[$i] > $genome_len ){  # need to wrap  
      # adding end of genome portion
      $frag_seq = $seqo->subseq( $fragStarts_r->[$i], $genome_len );
      
      # adding beginning of genome protion
      ## wrapped fragment end
      my $wrap_end = $fragEnds_r->[$i] - $genome_len;  # bp from genome start
      $frag_seq .= $seqo->subseq( 1, $wrap_end);
    }
    else{  # no need to wrap
      $frag_seq = $seqo->subseq( $fragStarts_r->[$i], $fragEnds_r->[$i] );
    }
  #  push @seqs, $frag_seq;
  }    

}


=head2 getGenomeFrags

Extracting fragment sequences from the genome.
Adding to $frags_r.

IF frag_end > genome_len:
    need to wrap frag (end of frag at start of linear genome)

=head3 IN


=head3 OUT

=cut

sub getGenomeFrags{
  my $fragStarts_r = shift or confess "Provide fragStarts_r\n";
  my $fragEnds_r = shift or confess "Provide fragEnds_r\n";
  my $fragSizes_r = shift or confess "Provide fragSizes_r\n";
#  my $ampsPerFrag_r = shift;
  my %h = @_;


  # IOcheck
  map{ exists $h{$_} or confess "Provide arg: $_\n" }
    qw/ -genome_scaf_id -genome_id -genome_len/;
  exists $h{-input}{-genome_db} or confess "Provide -input -> -genome_db\n";
  ## MCE/fork variables
  map{ exists $h{-argv}{$_} or confess "Provide -argv -> $_\n" }
    qw/ -nw -cs /;

  
  # gettig bio::db::fasta object of genomes
  my $seqo = $h{-input}{-genome_db}->get_Seq_by_id( $h{-genome_scaf_id} );


  # extracting sequences and calculating GC / BD  
  ## forking setup
  my $pm2 = Parallel::ForkManager->new( $h{-argv}{-nw} );
  my @allFragInfo;
  $pm2->run_on_finish(
  		     sub{
  		        my ($pid, $exit_code, $ident, 
  			    $exit_signal, $core_dump, $frag_info_r) = @_;
			
			# compiling info on all fragments
			push @allFragInfo, @$frag_info_r;
  		       }
  		     );

  # parsing fragments into chunks 
  my $chunk_size = $h{-argv}{-cs};
  my @chunks;
  for (my $i=0; $i<=$#$fragStarts_r; $i += $chunk_size){
    my $chunk_end = $i + $chunk_size - 1;
    $chunk_end = $#$fragStarts_r if $chunk_end > $#$fragStarts_r;
    push @chunks, [$i , $chunk_end];
  }
  

  ## iterating through chunks
  for my $chunk_r (@chunks){
    # forking by chunk
    $pm2->start and next;

    my @frag_info;
    for my $i ($chunk_r->[0]..$chunk_r->[1]){
      # sanity check
      confess "Internal error at line $.\n"
  	unless defined $fragEnds_r->[$i];
      
      # getting genome fragment
      ## frag may 'wrap' around linear genome. Need to account for this
      my $frag_seq = "";
      if( $fragEnds_r->[$i] > $h{-genome_len} ){  # need to wrap  
  	# adding end of genome portion
  	$frag_seq = $seqo->subseq( $fragStarts_r->[$i], $h{-genome_len} );
	
  	# adding beginning of genome protion
  	## wrapped fragment end
  	my $wrap_end = $fragEnds_r->[$i] - $h{-genome_len};  # bp from genome start
       $frag_seq .= $seqo->subseq( 1, $wrap_end);
      }
      else{  # no need to wrap
  	$frag_seq = $seqo->subseq( $fragStarts_r->[$i], $fragEnds_r->[$i] );
      }    
      
      
      # calculating GC
      my $GC = calc_GC_fast( $frag_seq );
      
      # calculating BD
      my $BD = calc_BD( $GC, 
  			$h{-genome_id}, 
  			$h{-input}{-if}, 
  			$h{-sample},
  			$h{-argv}{-isotope});

      # simulating diffusion
      #$BD = diffusion_brownian($BD, $h{-argv}{-diffusion});

      # loading values
      push @frag_info, [ $GC, $BD, $fragSizes_r->[$i] ];
    }
    
    # end fork; return values to parent    
    $pm2->finish( 0, \@frag_info );
  }
 $pm2->wait_all_children;


  # returning fragments info: BD & fragment length
  return \@allFragInfo; 

}


=head2 calc_BD

Calculating buoyant density for a fragment.
Adding isotope shift if provided. 

=head3 IN

=head3 OUT

=cut

sub calc_BD{
  my $GC = shift or confess "Provide GC content\n";
  my $genome_id = shift or confess "Provide genome_id\n";
  my $if_r = shift or confess "Provide if_r\n";
  my $sample = shift or confess "Provide sample\n";
  my $isotopes = shift or confess "Provide isotops\n";
  my %h = @_;

  #IO check
  ## checking GC
  return 'NA' if $GC !~ /$RE{num}{real}/; 


  ## using other sample unless sampleID exists in incorp file (if incorp file has 1 sample)
  if( ! exists $if_r->{$genome_id}{$sample} ){    
    foreach my $randSample ( keys %{$if_r->{$genome_id}} ){   # using selecting sample ID randomly from incorp file
      print STDERR "WARNING: could not find: '" .
	$if_r->{$genome_id}{$sample} . "' in incorp file. Using $randSample\n";
      $sample = $randSample;
      last;
    }
  }
  ## percent incorp IO check
  my $perc_incorp = exists $if_r->{$genome_id}{$sample} ?
    $if_r->{$genome_id}{$sample} :
      confess "ERROR: Cannot find genome:$genome_id, sample:$sample in -if file\n";
  
  
  # making istope => BD-shift index
  my %iso_index = ( '13C' => 0.036,
		    '15N' => 0.016 ); 
  

  # calculating BD
  ## no incorp
  my $BD  = ($GC / 100) * 0.098 + 1.66;  # 0% isotope incorp
   
  ## adding isotop incorp (possilbe mutlple isotopes)
  foreach my $iso ( @$isotopes ){
    my $max_shift = exists $iso_index{ $iso } ?
      $iso_index{$iso} : confess "Isotope: $iso not supported\n";
    
    ## max shift for isotope * percent incorp
    $BD += $max_shift * ($if_r->{$genome_id}{$sample} / 100);
  }

  #print Dumper $BD; exit;
  return $BD;
}








=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-isopycnic::fragGenome at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=isopycnic::fragGenome>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc isopycnic::fragGenome


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=isopycnic::fragGenome>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/isopycnic::fragGenome>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/isopycnic::fragGenome>

=item * Search CPAN

L<http://search.cpan.org/dist/isopycnic::fragGenome/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of isopycnic::fragGenome
