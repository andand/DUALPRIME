#!/usr/bin/perl -w

##### FindSharedRegions.pl #####

#2 mm: en exakt 6mer & minst en exakt 4mer.         666666111122223333   hexapos: 0  quadrapos: 0,4,8 (position i återstående 12mer)
#                                                   111122666666223333            6             0,4,8
#                                                   111122223333666666            12            0,4,8
#
#                                                   1111,2222,3333 = 4merer
#                                                           666666 = 6merer
#I alla gener i genom 1: Flytta ett fönster a 18 baser över genen = 18 mer. Spara 18 meren på 9 sätt i hashen
# $mer{position på 6mer}{6mer sekvens} {position på 4mer i förhållande till 6mer}{4mer sekvens}

=head1 NAME

FindSharedRegions.pl - Finds 18mers shared between genes of two genomes. Allows two mismatches within the 18 bases. Output file includes all shared 18mer sequences and positions of these in each genome.

=head1 SYNOPSIS

FindSharedRegions.pl -fasta1 genome1.fasta -fasta2 genome2.fasta

=head1 DESCRIPTION

The program searches each 18mer of each gene sequence in the first genome against each 18mer of each gene sequence of the second genome. If an 18mer occurs in both genomes it will be typed on a row in the output file. Two mismatches are allowed in the 18mers. Except the two exact 18mer sequences, the genes in which they occur, the positions in the genes, and the orientation of the 18mer in the second genome (forward or reverse complementary) is printed. The output file will be named genome1_genome2_regions.tmp

=head1 OPTIONS

(these are required)

-fasta1 <genome1.fasta> - gene sequences of genome 1 in fasta format

-fasta2 <genome2.fasta> - gene sequences of genome 2 in fasta format

=head1 AUTHOR - Anders Andersson

Anders Andersson E<lt>ason@biotech.kth.seE<gt>

=cut

use Getopt::Long;

$infile1 = undef;
$infile2 = undef;
$help = undef;
&GetOptions('fasta1=s' => \$infile1, 'fasta2=s' => \$infile2, 'h!' => \$help);
if ($help or !$infile1 or !$infile2) {
  system ('perldoc', $0);
  exit;
}
if (substr($infile1,-6,6) ne ".fasta") {
  die ("$infile1 is no fasta file");
} else {
  $prefix1 = $infile1;
  substr($prefix1, -6,6) = "";
}
if (substr($infile2,-6,6) ne ".fasta") {
  die ("$infile2 is no fasta file");
} else {
  $prefix2 = $infile2;
  substr($prefix2, -6,6) = "";
}
$organism1 = "$prefix1"."_tab.tmp";
$organism2 = "$prefix2"."_tab.tmp";

$outfile = $prefix1."_".$prefix2."_regions.tmp";

print"--- Finding shared regions ---\n";

&convert_file($infile1,$organism1);
&convert_file($infile2,$organism2);
$start = 0;
$end = 601;
$genesleft = defined;
open(OUT, ">$outfile");
while (defined $genesleft) {
  &build_db($organism1);
  &search_db($organism2);
  &rev_comp_search_db($organism2);
  $start = $start + 600;
  $end = $end + 600;
  %mer = ();
  @hexa_pos = ();
  @quadra_pos = ();
}
close(OUT);
unlink($organism1);
unlink($organism2);

sub convert_file {
  local($infile) = $_[0];
  local($outfile) = $_[1];
  open (INFILE, $infile) || die ("could not open $ARGV[1] $!");
  open (OUT, ">$outfile") || die ("could not open $outfile");
  $line = 0;
  while (<INFILE>) {
    $line++;
    chomp $_;
    if (/^>/) {
      if ($line > 1) {
	print OUT $name."\t".$seq."\n";
      }
      $name = $_;
      $seq ="";
    } elsif(/^\w/) {
      $seq = $seq.$_;
    } else {
      print"$_\n";
      die ("something strange with $infile, line $line \n"); 
    }
  }
  print OUT $name."\t".$seq."\n";
  close (OUT);
}

sub build_db {
  $infile_db = $_[0];
  open (ORFFILE1, $infile_db) || die ("could not open!");
  @hexa_pos = (0,6,12);
  @quadra_pos = (0,4,8);
  $nr=1;
  while (<ORFFILE1>) {
    if ($nr > $start) {
      if ($nr < $end) {
	#print"$nr\n";
	$genesleft = undef;
	chomp $_;
	@fields = split(/\t/);
	$geneseq = $fields[1];
	$l=length($geneseq);
	for ($i=0; $i<$l-17; $i++) {
	  $window = substr($geneseq, $i, 18);
	  foreach $hexapos (@hexa_pos) { 
	    $hexamer = substr($window, $hexapos, 6);
	    $decamer = substr($window, 0, $hexapos).substr($window, $hexapos+6, 12); #sista substr out of range?
	    foreach $quadrapos (@quadra_pos) {
	      $quadramer = substr($decamer, $quadrapos, 4);
	      $octamer = substr($decamer, 0, $quadrapos).substr($decamer, $quadrapos+4, 12);
	      push(@{$mer{$hexapos}{$hexamer}{$quadrapos}{$quadramer}},$octamer."\t$nr\t$i");
	    }
	  }
	}
      } else { $genesleft = defined; }
    }
    $nr++;
  }
  close (ORFFILE1);
}

sub search_db {
  $infile_query = $_[0];
  @hexa_pos = (0,6,12);
  @quadra_pos = (0,4,8);
  local($nr);
  $nr=0;
  open (ORFFILE2, "$infile_query") || die ("could not open !");
  while (<ORFFILE2>) {
    $nr++;
    #print"\t$nr\n";
    chomp $_;
    @fields = split(/\t/);
    $geneseq = $fields[1];
    $l=length($geneseq);
    for ($i=0; $i<$l-17; $i++) {
      $window = substr($geneseq, $i, 18);
      %allready = ();
      foreach $hexapos (@hexa_pos) {
	$hexamer = substr($window, $hexapos, 6);
	$decamer = substr($window, 0, $hexapos).substr($window, $hexapos+6, 12);
	foreach $quadrapos (@quadra_pos) {
	  $quadramer = substr($decamer, $quadrapos, 4);
	  $octamer = substr($decamer, 0, $quadrapos).substr($decamer, $quadrapos+4, 12);
	  if (defined $mer{$hexapos}{$hexamer}{$quadrapos}{$quadramer}) {
	    foreach $otherocta (@{$mer{$hexapos}{$hexamer}{$quadrapos}{$quadramer}}) {
	      @fields = split(/\t/,$otherocta);
	      if (defined $allready{"$fields[1]\t$fields[2]"}) {
	      } else {
		$allready{"$fields[1]\t$fields[2]"} = 1;#if the window allready hit this position and gene
		$match = 0;
		for ($pos=0; $pos<8; $pos++) {
		  if (substr($octamer,$pos,1) eq substr($otherocta,$pos,1)) {
		    $match++;
		  }
		}
		if ($match > 5) {
		  $slide = substr($otherocta, $quadrapos, (8-$quadrapos));
		  $rest = substr($otherocta, 0, $quadrapos);
		  $otherdeca = $rest.$quadramer.$slide;
		  $slide = substr($otherdeca, $hexapos, (12-$hexapos));
		  $rest = substr($otherdeca, 0, $hexapos);
		  $otherwindow = $rest.$hexamer.$slide;
		  print OUT "$infile_db\t$fields[1]\t$fields[2]\t$otherwindow\tfor\t$infile_query\t$nr\t$i\t$window\n";
		}
	      }
	    }
	  }
	}
      }
    }	
  }
  close(ORFFILE2);
}

sub rev_comp_search_db {
  $infile_query = $_[0];
  @hexa_pos = (0,6,12);
  @quadra_pos = (0,4,8);
  local($nr);
  $nr=1;
  open (ORFFILE2, "$infile_query") || die ("could not open !");
  while (<ORFFILE2>) {
    #print"\tR$nr\n";
    chomp $_;
    @fields = split(/\t/);
    $geneseq = $fields[1];
    $l=length($geneseq);
    $geneseq = &make_revcomp;
    for ($i=0; $i<$l-17; $i++) {
      $window = substr($geneseq, $i, 18);
      %allready = ();
      foreach $hexapos (@hexa_pos) {
	$hexamer = substr($window, $hexapos, 6);
	$decamer = substr($window, 0, $hexapos).substr($window, $hexapos+6, 12);
	foreach $quadrapos (@quadra_pos) {
	  $quadramer = substr($decamer, $quadrapos, 4);
	  $octamer = substr($decamer, 0, $quadrapos).substr($decamer, $quadrapos+4, 12);
	  if (defined $mer{$hexapos}{$hexamer}{$quadrapos}{$quadramer}) {
	    foreach $otherocta (@{$mer{$hexapos}{$hexamer}{$quadrapos}{$quadramer}}) {
	      @fields = split(/\t/,$otherocta);
	      if (defined $allready{"$fields[1]\t$fields[2]"}) {
	      } else {
		$allready{"$fields[1]\t$fields[2]"} = 1;
		$match = 0;
		for ($pos=0; $pos<8; $pos++) {
		  if (substr($octamer,$pos,1) eq substr($otherocta,$pos,1)) {
		    $match++;
		  }
		}
		if ($match > 5) {
		  $slide = substr($otherocta, $quadrapos, (8-$quadrapos));
		  $rest = substr($otherocta, 0, $quadrapos);
		  $otherdeca = $rest.$quadramer.$slide;
		  $slide = substr($otherdeca, $hexapos, (12-$hexapos));
		  $rest = substr($otherdeca, 0, $hexapos);
		  $otherwindow = $rest.$hexamer.$slide;
		  print OUT "$infile_db\t$fields[1]\t$fields[2]\t$otherwindow\trev\t$infile_query\t$nr\t$i\t$window\n";
		}
	      }
	    }
	  }
	}
      }
    }
    ++$nr;	
  }
  close(ORFFILE2);
}

sub make_revcomp {
  local($base);
  local($i);
  $revcomp = "";
  for ($i=0; $i<$l; $i++) {
    $base = substr($geneseq, $i, 1);
    if ($base eq "A") {
      $revcomp = "T".$revcomp;
    } elsif ($base eq "T") {
      $revcomp = "A".$revcomp;
    } elsif ($base eq "C") {
      $revcomp = "G".$revcomp;
    } elsif ($base eq "G") {
      $revcomp = "C".$revcomp;
    } elsif ($base eq "N") {
      $revcomp = "N".$revcomp;
    } else {
      print"\nStrange base in gene sequence: \n$geneseq\n";
      die;
    }
  }
  $revcomp;
}
