package isopycnic::IO;

use 5.006;
use strict;
use warnings;

=head1 NAME

isopycnic::IO - scripts for running isopycnic::IO.pl

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
use Bio::TreeIO;
use Regexp::Common;
use IPC::Cmd qw/can_run run/;
use Parallel::ForkManager;
use Hash::MultiKey;

use isopycnic::t;
our $argv_err;


=head2 read_tree

Reading in tree as a TreeIO object.
Just returning 1st tree (if multiple provided).

=head3 IN

hash of args for Bio::TreeIO->new

=head3 OUT

TreeIO object

=cut

push @EXPORT_OK, 'read_tree';

sub read_tree{
  my %h = @_;

  my $treeio = Bio::TreeIO->new( %h );

  while( my $tree = $treeio->next_tree ){
    return $tree;
  }
}



=head2 read_af_if

Reading csv (either abundance_file or incorp_file)

=head3 IN

-fh --  file handle
-file --  file name
-header -- return header? (1=True)

=head3 OUT

{taxon}{sample}{abund}

=cut

push @EXPORT_OK, 'read_af_if';

sub read_af_if{
  my %h = @_;
  confess "ERROR: provide -fh or -file\n"
    unless exists $h{-file} or exists $h{-fh};

  my $fh;
  exists $h{-fh} ? $fh = $h{-fh} :
    open $fh, $h{-file} or confess $!;

  my %tbl;
  my %header;
  while(<$fh>){
    chomp;
    next if /^\s*$/;
    s/[\t ]+$//;
    my @l = split / *\t */;
    croak "ERROR: line $. does not contain >=2 columns!\n"
      unless scalar @l >= 2;

    # header
    if( $. == 1 ){
      if(/^#/ ){  # header provided
	for my $i (1..$#l){
	  $header{$i} = $l[$i];
	}
      }
      else{  # no header assumed
	for my $i (1..$#l){
	  $header{$i} = "Sample$i";  # making header IDs
	}	
	for my $i (1..$#l){  # 1 line of body
	  $tbl{$l[0]}{$header{$i}} = $l[$i];
	}	
      }
    }
    else{  # body
      for my $i (1..$#l){
	$tbl{$l[0]}{$header{$i}} = $l[$i];
      }
    }
  }
  close $fh;

  # return
  if( $h{-header} ){
    return \%tbl, \%header;
  }
  else{
    return \%tbl;
  }
}


=head2 af_from_factory

Making af DS from grinder factory object

=head3 IN

factory object

=head3 OUT

{taxon}{sample}{abund}

=cut

push @EXPORT_OK, 'af_from_factory';
sub af_from_factory{
  my $factory = shift || confess "Provide Grinder \$factory object\n";
  my %opts = @_;
  my $header_r = $opts{-header} if defined $opts{-header};  # optional header

  my %af;
  my $sample = 0;   # sample count
  while( my $struct = $factory->next_lib ){
    $sample++;
    for my $i (0..$#{$struct->{ids}}){  # each taxon ID
      my $taxon_id = $struct->{ids}->[$i];

      # editing the taxon_id (grinder adds extra to each name)
      $taxon_id =~ s/\/.+//;      
      
      # determining sample ID
      my $sampleID;
      if( defined $header_r and exists $header_r->{$sample} ){
	$sampleID = $header_r->{$sample};
      }
      else{
	print STDERR "WARNING: cannot find header value number $sample." .
	  " Using Sample$sample\n"
	    if exists $opts{-verbose} and $opts{-verbose} > 0;
	$sampleID = "Sample$sample";
      }

      # loading hash
      $af{ $taxon_id }{ $sampleID } = $struct->{abs}->[$i];
    }
  }

#  print Dumper \%af; exit;
  return \%af;
}


=head2 writeFrags

Writing fragments to fraction files.

=head3 IN

kwargs:
-frag_info -- [frag_sequence, BD]
-genome_id -- $genome_ID
-frac_fh -- [ file : filename, fh : filehandle, min : BD_min, max : BD_max ]

=head3 OUT

=cut

push @EXPORT_OK, 'writeFrags';

sub writeFrags{
  my %h = @_;

  # IOcheck
  map{ exists $h{$_} or confess "Provide arg: $_\n" }
    qw/ -frag_info -chunk_start -genome_scaf_id -frac_fh /;
  my $frag_info_r = $h{-frag_info};

  # writing fragments to file
  
  foreach my $i ( 0..$#$frag_info_r ){
    # fragID based on chunk_start
    my $frag_cnt_id = $i + $h{-chunk_start};

    # next unless real number for BD values (check for no fragment in fraction)
    next unless $frag_info_r->[$i][1] =~ /$RE{num}{real}/;
    
    # placing the fragment in a the correct fraction file
    foreach my $frac_r ( @{$h{-frac_fh}} ){

      # if within range, get sequence, write
      if( $frag_info_r->[$i][1] >= $frac_r->{min} and
	  $frag_info_r->[$i][1] < $frac_r->{max} ){
	
   	## writing to file
   	my $fh = $frac_r->{fh};	
   	print $fh join("\n", 
		       join("", ">", $h{-genome_scaf_id}, "__$frag_cnt_id"),
		       $frag_info_r->[$i][0]), "\n";
      }
    }    
  }
}


# extracting sequences
sub extract_seq{
  my $ret_r = shift or confess "Provide ret_r\n";
  my $start = shift or confess "Provide the sequence start\n";
  my $end = shift or confess "Provide the sequence end\n";
  my %h = @_;


  # input check
  exists $ret_r->{-genome_scaf_id} or confess "Provide -genome_scaf_id\n";
  my $genome_len = exists $ret_r->{-genome_len} ?
    $ret_r->{-genome_len} : confess "Provide -genome_len\n";
  exists $h{-input}{genome_db} or confess "Provide -input -> genome_db\n";
  my $seqo = $h{-input}{genome_db}->get_Seq_by_id( $ret_r->{-genome_scaf_id} );

  my $frag_seq = "";
  if( $end > $genome_len ){  # need to wrap
    # adding end of genome portion
    $frag_seq = $seqo->subseq( $start, $genome_len );
    
    # adding beginning of genome protion
    ## wrapped fragment end
    my $wrap_end = $end - $genome_len;  # bp from genome start
    $frag_seq .= $seqo->subseq( 1, $wrap_end);
  }
  else{
    $frag_seq = $seqo->subseq( $start, $end );
  }
  
  return $frag_seq;
}



=head2 closeFracFh

closing all fraction file handles

=head3 IN

=head3 OUT

=cut

push @EXPORT_OK, 'closeFracFh';

sub closeFracFh{
  my $fracs_r = shift or confess "Provide a \$fracs_r\n";
  my $argv_r = shift or confess "Provide GetOpt::Euclid \%ARGV\n";
  my %h = @_;

  # IO check
  confess "Cannot find '-output' arg\n"
    unless exists $argv_r->{-output};


  # closing all filehandles
  foreach my $frac ( @$fracs_r  ){
    close $frac->{fh} or die $!;
  }
}


=head2 gzipFiles

Using gzip to compress all files

=head3 IN

fracs_r : {sample : [ 'file' : value ]}

=head3 OUT

=cut

push @EXPORT_OK, 'gzipFiles';

sub gzipFiles{
  my $fracs_r = shift or confess "Provide a \$fracs_r\n";
  my $argv_r = shift or confess $argv_err;
  my %h = @_;

  # status
  print STDERR "Compressing all fasta files with gzip...\n"
    unless $argv_r->{'--quiet'};

  # checking for gzip
  can_run('gzip') || confess "Cannot gzip fasta files! gzip is not in \$PATH\n";

  # forking
  exists $argv_r->{'-gzip'} or confess "ERROR: provide '-gzip'\n";
  my $pm = Parallel::ForkManager->new( $argv_r->{'-gzip'} );

  # closing all filehandles
  foreach my $sample ( keys %$fracs_r ){
    foreach my $frac ( @{$fracs_r->{$sample}} ){
      $pm->start and next;
      
      -e $frac->{file} 
	or confess "Cannot find: " . $frac->{file}, "\n";
      
      my $cmd =  join(" ", 'gzip', $frac->{file});
      
      my $verbose = $argv_r->{'--quiet'} ? 0 : 1;  # inverse of --quiet
      my( $success, $error_message, $full_buf,
	  $stdout_buf, $stderr_buf ) =
	    run( command => $cmd,
		 verbose => $verbose );
      
      $pm->finish;
    }
  }
  $pm->wait_all_children;
}




=head2 write_af_ranks

Writing out *-ranks file for community abundance.
Taking af file as hashref and writing as
Grinder-formatted *-ranks.txt file.

Ranks file: # rank  seq_id  rel_abund_perc  abs_abund

=head3 IN

$af_r :  hashref {taxon : sample : value}
$argv_r :  getopt::euclide \%ARGV 

=head3 OUT

Ranks by sort order (sorting by relative abundance)

=cut

push @EXPORT_OK, 'write_af_ranks';

sub write_af_ranks{
  my $af_r = shift || confess "Provide af as hashref\n";
  my $argv_r = shift || confesss $argv_err;

#  print Dumper $argv_r; exit;


  # output file
  my ($outfh, $outfile);
  if( exists $argv_r->{-output_prefix} ){
    $outfile = $argv_r->{-output_prefix} . "-ranks.txt";
    open $outfh, ">$outfile" or die $!;
  }
  else{
    $outfh = \*STDOUT;
  }
 
  # header 
  print $outfh join("\t", qw/sample rank seq_id rel_abund_perc abs_abund/), "\n";

  # body
  ## getting samples
  my %samples;
  my $rank = 0; 
  foreach my $taxon (keys %$af_r){
    foreach my $sample ( keys %{$af_r->{$taxon}} ){
      $samples{ $sample } = 1;
    }
  }


  ## body by sample
  foreach my $sample (keys %samples){
    # getting taxa in sample
    my @taxa;
    map{ push @taxa, $_ if exists $af_r->{$_}{$sample} }keys %$af_r;
    
    # sorting by rank
    $rank = 0;
    foreach my $taxon ( sort{$af_r->{$b}{$sample} <=> $af_r->{$a}{$sample}} 
			@taxa ){
      $rank += 1;
      print $outfh join("\t", 
			$sample, 
			$rank, 
			$taxon, 
			$af_r->{$taxon}{$sample} * 100, 
			$af_r->{$taxon}{$sample} * $argv_r->{-total_copies}
		       ), "\n";
    }
  }
  
  close $outfh;

  # status
  print STDERR "  *-ranks.txt file with pre-isopycnic gradient taxon abundances written to: $outfile\n"
    if ! $argv_r->{'--quiet'} and defined $outfile;

}


=head2 readFracFile

Reading in tab-delimited file provided by isopycnic.pl

=head3 IN

infile -- $ infile name

=head3 OUT

[genome_id : sample_id : [ cat : value] ]

=cut

push @EXPORT_OK, 'readFracFile';
sub readFracFile{
  my $infile = shift or confess "Provide infile\n";

  open IN, $infile or confess $!;
  my %fracTbl;
  my %header; 
  while(<IN>){
    chomp;
    next if /^\s*$/;
    
    my @l = split /\t/;
    scalar @l >= 7 or 
      confess "ERROR: line $. of '$infile' has < 7 columns!\n";
    
    # fraction table
    my $BD_range = join("-", $l[2], $l[3]);
    $fracTbl{$l[0]}{$l[1]}{ $BD_range } = { 
					   N_temp_by_frac => $l[4],
					   N_temp_total => $l[5],
					   rel_temp => $l[6] 
					  };
    

    # header (samples)
    $header{ $l[1] } = scalar(keys %header) + 1
      unless exists $header{$l[1]};
  }
  # flipping key-values for header
  my %new_header;
  map{ $new_header{ $header{$_} } = $_ } keys %header;

#  print Dumper %fracTbl; exit;
  return \%fracTbl, \%new_header;
}


=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-isopycnic::IO at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=isopycnic::IO>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc isopycnic::IO


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=isopycnic::IO>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/isopycnic::IO>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/isopycnic::IO>

=item * Search CPAN

L<http://search.cpan.org/dist/isopycnic::IO/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of isopycnic::IO
