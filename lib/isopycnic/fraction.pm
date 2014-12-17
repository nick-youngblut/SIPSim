package isopycnic::fraction;

use 5.006;
use strict;
use warnings;

=head1 NAME

isopycnic::fraction - scripts for running isopycnic::fraction.pl

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
use Math::Random qw/random_normal/;
use File::Spec;
use File::Path qw/rmtree/;
use IO::Compress::Gzip qw/gzip $GzipError/;

use isopycnic::GC qw/sliding_window_GC/;
use isopycnic::t;
our $argv_err;


=head2 makeFractions

Workflow for making fraction file handles

=head3 IN

=head3 OUT

=cut

push @EXPORT_OK, 'makeFractions';

sub makeFractions{
  my $argv_r = shift || confess $argv_err;
  my %h = @_;
  
  # input check
  my $genome_db = exists $h{-input}{-genome_db} ?
    $h{-input}{-genome_db} : confess "Provide a Bio::DB::Fasta object\n";
  my $if_r = exists $h{-input}{-if} ? 
    $h{-input}{-if} : confess "Provide incorp hashref\n";
  

  # min-max BD determined by theoretical min-max GC + iso incorp
  my %BD_range = (
		  iso_min => 1.66 + 0,                    # 0% incorp, 0% GC
		  iso_max => 0.098 + 1.66 + 0.036 + 0.016 # 100% incorp of 13C & 15N, 100% GC
		 );
  

  # making BD fractions
  ## getting sampleIDs
  my $samples_r = exists $h{-input}{-header} ? 
    $h{-input}{-header} : confess "Provide -input => header\n";
  
  ## iterating over sample 
  my %fracs;
  foreach my $sample (values %$samples_r){
    $fracs{$sample} = make_BD_fracs(\%BD_range, $argv_r);
#    make_frac_ranges($sample,  $fracs{$sample}, $argv_r, %h);
  }
  
  # return faction file handles (& BD value ranges)
#  print Dumper %fracs;
  return \%fracs;
}



=head2 make_BD_fracs

Making BD fractions, which
are ranges of BD values
resulting from isopynic gradient
fractionation.

Making fractions for both control 
& isotope assimilation.

=head3 IN

$BD_range -- Min & max BD for control & assimilation. {
$argv_r -- GetOpt::Euclid \%ARGV

=head3 OUT

{iso|no_iso}{'fracs'} => [frac_start, frac_end]

=cut

push @EXPORT_OK, 'make_BD_fracs';

sub make_BD_fracs{
  my $BD_range = shift or confess "Provide BD_range as \%.\n";
  my $argv_r = shift or confess "Provde GetOpt::Euclid \%ARGV\n";
  

  # IO check
  ## setting defaults
  unless (exists $argv_r->{-fraction}){
    $argv_r->{-fraction} = { frac_mean => 0.004,
			     frac_stdev => 0.0015,
			     frac_min => 0.001,
			     frac_max => 0.006 };
  }
  map{ confess "Cannot find '$_' in args\n"
	 unless exists $argv_r->{-fraction}{$_} } qw/frac_mean
						     frac_stdev
						     frac_min
						     frac_max/;
  map{ confess "Cannot find '$_' in args\n"
	 unless exists $BD_range->{$_} } qw/iso_min iso_max/;

  sub cat_fracs{  # sub for making all fractions spanning BD range
    my $BD_min = shift or die $!;
    my $BD_max = shift or die $!;
    my $argv_r = shift or die $!;

    my @fracs;
    while(1){
      # get fraction size; BD_min will be upated with last BD range
      my $frac_size = get_frac($argv_r);
      push @fracs, {min => $BD_min, max => $BD_min + $frac_size};  # frac_start-end

      $BD_min += $frac_size;
      if($BD_min >= $BD_max){ # ending loop; should have enough fractions to span gradient
	$fracs[$#fracs]->{max} += 0.0001;   # making max a bit larger so that all fragments added to a fraction
	return \@fracs; 
      }
    }
  }

  sub get_frac{  # sub for finding fraction size; using a normal distribution
    my $argv_r = shift or confess "Provide \$argv_r\n";

    my $timeout = 0;
    my $timeout_max = 9999;
    while(1){
      my $frac_size = random_normal(1, $argv_r->{-fraction}{frac_mean},
				    $argv_r->{-fraction}{frac_stdev});

      return $frac_size if $frac_size >= $argv_r->{-fraction}{frac_min} and
	$frac_size <= $argv_r->{-fraction}{frac_max};
      
      $timeout++;
      my $err = join(" ", "ERROR: could not make fraction >=",
		     $argv_r->{-fraction}{frac_min},
		     " & <= ", $argv_r->{-fraction}{frac_max},
		     " in $timeout_max tries");
      confess "ERROR: $err\n" if $timeout >= $timeout_max;
    }
  }
  

  # actually getting fractions
  ## isotope
  my $fracs_r = cat_fracs($BD_range->{iso_min}, $BD_range->{iso_max}, $argv_r);
  
  #print Dumper $fracs_r; exit;
  return $fracs_r;  
}



=head2 make_frac_files

Making a file handle for each fraction.
Using -output flag as file prefix.

=head3 IN

$fracs_r -- {iso|no_iso}=> [{min|max => value}]
$argv_r -- GetOpt::Euclid

=head3 OUT

edited fracs_r

=cut 

push @EXPORT_OK, 'make_frac_files';

sub make_frac_files{
  my $sample = shift or confess "Provide a sample\n";
  my $fracs_r = shift or confess "Provide a \$fracs_r\n";
  my $argv_r = shift or confess "Provide GetOpt::Euclid \%ARGV\n";
  my %h = @_;
  
  # IO check
  confess "Cannot find '-output' arg\n"
    unless exists $argv_r->{-output};

  # making sample subdirectory
  my $subDir = File::Spec->catdir($argv_r->{-output}, $sample);
  mkdir $subDir or croak $! unless -d $subDir;  
    
  # making file handles    
  ## iterating over fractions
  foreach my $frac ( @$fracs_r ){  #$frac = \%
    # file name: contains BD min-max for fraction
    my $file_name = join(".",
			 join("-", 
			      sprintf("%.4f", $frac->{min}),
			      sprintf("%.4f", $frac->{max})
			     ),
			 "fasta");
    $file_name = File::Spec->catfile($subDir, $file_name);  # adding directory
      
    open $frac->{fh}, ">$file_name" or confess $!; 
    $frac->{file} = $file_name;
  }

  #  print Dumper $fracs_r; exit;
  return 1;
}


=head2 countByFrac

Fragment will be binned by fraction.
If only amplicon fragments:
 * fraction abundance = frags_in_frac / total_frags_made = relative amplicon template copy number
ElseIf all genome fragments:
 * fraction abundance = length_frags_in_frac / total_length_frags = relative genome copy number

=head3 IN

$allFragInfo_r -- [ [BD, frag_len] ]
kwargs

=head3 OUT

STDOUT

=cut

push @EXPORT_OK, 'countByFrac';
sub countByFrac{
  my $allFragInfo_r = shift or confess "Provide allFragInfo\n";
  my $ampsPerFrag_r = shift;
  my %h = @_;

  # input check
  map{ exists $h{$_} or confess "Provide $_\n" }
    qw/ -genome_scaf_id -genome_id -genome_len -frac_ranges -argv -sample/;
  map{ exists $h{-argv}{$_} or confess "Provide -argv->$_\n"}
    qw/ -total_fragments /;
  map{ exists $h{-input}{$_} or confess "Provide -input->$_\n"}
    qw/ -amp_region_index/;

  # summing totals for starting community
  my ($totalAmpTemp, $totalFragLen) = (0,0);
  if(exists $h{-input}{-amp_region_index} ){
    $#$allFragInfo_r == $#$ampsPerFrag_r or confess "Internal Error\n";
    for my $i (0..$#$allFragInfo_r){
      $totalAmpTemp += $ampsPerFrag_r->[$i];
    }
  }
  else{
    # summing total fragment length if needed
    map{ $totalFragLen += $_->[1] } @$allFragInfo_r;
  }

  # by gradient fraction
  foreach my $frac_range_r ( sort{$a->{min} <=> $b->{min}} @{$h{-frac_ranges}} ){

    ## just amplicon fragments
    if( exists $h{-input}{-amp_region_index} ){

      ## summing by total amplicon templates in fraction
      my $totalInFrac = 0;
      for my $i (0..$#$allFragInfo_r){
	if( $allFragInfo_r->[$i][0] >= $frac_range_r->{min}        # is BD in fraction?
	    and  $allFragInfo_r->[$i][0] < $frac_range_r->{max} ){
	  $totalInFrac += $ampsPerFrag_r->[$i];
	}
      }

      # skipping unless taxon in fraction
      next unless $totalInFrac > 0;

      # writing line
      print join("\t", 
		 $h{-genome_scaf_id},
#		 $h{-genome_id},
		 $h{-sample}, 
		 sprintf("%.3f", $frac_range_r->{min}),
		 sprintf("%.3f", $frac_range_r->{max}),
		 $totalInFrac, 
		 $totalAmpTemp,
		 $totalInFrac / $totalAmpTemp
		 ), "\n";
		 
    }
    ## fragments from across the genome
    else{  
      my $totalLen = 0;
      map{ $totalLen += $_->[1] 
	     if $_->[0] >= $frac_range_r->{min} 
	       and $_->[0] < $frac_range_r->{max} } @$allFragInfo_r;
      
      # next unless taxon in fraction
      next unless $totalLen > 0;

      # writing line
      print join("\t", 
		 $h{-genome_id},
		 $h{-sample}, 
		 sprintf("%.3f", $frac_range_r->{min}),
		 sprintf("%.3f", $frac_range_r->{max}),		 
		 $totalLen, 
		 $totalFragLen,
		 $totalLen / $totalFragLen
		 ), "\n";      
    }
  }
}



=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-isopycnic::fraction at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=isopycnic::fraction>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc isopycnic::fraction


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=isopycnic::fraction>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/isopycnic::fraction>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/isopycnic::fraction>

=item * Search CPAN

L<http://search.cpan.org/dist/isopycnic::fraction/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of isopycnic::fraction
