#!/usr/bin/perl -w

##### FindSelfSimilarity.pl #####

=head1 NAME

FindSelfSimilarity.pl - Finds regions with similarity to other regions within the same genome/genomes

=head1 SYNOPSIS

FindSelfSimilarity.pl -tab1 genome1.tab -tab2 genome2.tab -l regionlength -s percentage similarity

=head1 DESCRIPTION

The program builds a BLAST database for each of the one or two genomes given in the tabfile(s). Each gene will be BLAST searched aginst each other gene in the same genome, and regions exceeding the lengths (number of nucleotides) and similarity (percentage identical nucleotides) given in the input will be printed in the output file (used by subsequent scripts). If no values are given on these the default values 30 nucleotides and 66% similarity will be used.

=head1 OPTIONS

(these are required)

-tab1 <genome1.tab> - gene sequences of genome 1 in tab format (output of MakeTabFiles.pl)

(optional)

-tab2 <genome2.tab> - gene sequences of genome 2 in tab format (output of MakeTabFiles.pl)

-l   - minimum value for region length is 30

-s   - minimum and maximum value for similarity is 66 and 99, respectively

=head1 AUTHOR - Anders Andersson

Anders Andersson E<lt>ason@biotech.kth.seE<gt>

=cut

use Getopt::Long;

$infile1 = undef;
$infile2 = undef;
$help = undef;
$maxhitlength = 30;
$maxpercentage = 66;
&GetOptions('fasta1=s' => \$infile1, 'fasta2=s' => \$infile2, 'h!' => \$help, 'l=i' => \$maxhitlength, 's=i' => \$maxpercentage);
if ($help or !$infile1) {
  system ('perldoc', $0);
  exit;
}
if (($maxhitlength < 30) or ($maxpercentage < 66) or ($maxpercentage > 99)) {
  system ('perldoc', $0);
  exit;
}

if (substr($infile1,-6,6) ne ".fasta") {
  die ("$infile1 is no fasta file");
} else {
  $prefix1 = $infile1;
  substr($prefix1, -6,6) = "";
  $organism1 = "$prefix1"."_tab.tmp";
}

if ($infile2) {
  if (substr($infile2,-6,6) ne ".fasta") {
    die ("$infile2 is no fasta file");
  } else {
    $prefix2 = $infile2;
    substr($prefix2, -6,6) = "";
    $organism2 = "$prefix2"."_tab.tmp";
    $outfile = $prefix1."_".$prefix2."_selfsim.tmp";
    $infile_mers = $prefix1."_".$prefix2."_regions.tmp";
    $ok = open (MERFILE, $infile_mers);
    if ($ok) {
      close (MERFILE);
    } else {
      $outfile = $prefix2."_".$prefix1."_selfsim.tmp";
      $infile_mers = $prefix2."_".$prefix1."_regions.tmp";
      $ok = open (MERFILE, $infile_mers);
      if ($ok) {
	close (MERFILE);
      } else {
	print"Couldn't find files with shared primers!\n";
	die;
      }
    }
  }
} else {
  $outfile = $prefix1."_selfsim.tmp";
}
$maxpercentage = $maxpercentage/100;

$blastinput = "blastinput.tmp";
$blastoutput = "blastoutput.tmp";

#--------------- Main Program -----------------------
print"--- Finding self similarity ---\n";
&convert_file($infile1,$organism1);
$infile_ORFs = $organism1;
&make_blastdb($prefix1);
&get_selfhits($prefix1);
if ($organism2) {
  &convert_file($infile2,$organism2);
  $infile_ORFs = $organism2;
  &make_blastdb($prefix2);
  &get_selfhits($prefix2);
  unlink($organism2,$infile2.".nhr",$infile2.".nin",$infile2.".nsq");
}
&print_selfhits;
unlink($organism1,$infile1.".nhr",$infile1.".nin",$infile1.".nsq");
unlink($blastinput,$blastoutput);
#----------------------------------------------------
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

sub make_blastdb {
  $fastafile = $_[0].".fasta";
  system ("formatdb -t $infile_ORFs -i $fastafile -p F");
}

sub get_selfhits {
  open (ORFFILE, $infile_ORFs) || die ("could not open ! $!");
  $genome = $infile_ORFs;
  $nr = 0;
  $fastafile = $_[0].".fasta";
  while (<ORFFILE>) {
    $nr++;
    #if ($nr > 50) { last; }
    chomp $_;
    @fields = split(/\t/);
    $realname = $fields[0];
    $gene = $nr;
    $geneseq = $fields[1];
    $realn{$genome}{$gene} = $realname;
    open (OUT, ">$blastinput");
    print OUT ">$gene\n$geneseq\n";
    close(OUT);
    system("blastall -p blastn -d  $fastafile -i $blastinput -q -2 -e 100 -b 100000 -o $blastoutput"); #Varför -q 2 ? Bättre på att plocka upp homologier.
    &extract_crossmatches($realname);
    #print"gene: $nr\tselfhits: $hitindex\n"; 
    #} else {
    #last;
    #}
    #}
  }
  close(ORFFILE);
}

sub print_selfhits {
  open(OUT, ">$outfile");
  print OUT "maxhitlength:\t$maxhitlength\tmaxpercentage:\t$maxpercentage\n";
  foreach $genome (keys %realn) {
    foreach $gene (keys %{$realn{$genome}}) {
      print OUT "$genome\t$gene\t$realn{$genome}{$gene}";
      if (defined $hitstarts{$genome}{$gene}) {
	$l = (@{$hitstarts{$genome}{$gene}});
	for ($i=0; $i<$l; $i++) {
	  print OUT "\t$hitstarts{$genome}{$gene}[$i]_$hitends{$genome}{$gene}[$i]";
	}
      } 
      print OUT "\n";
    }
  }
  close(OUT);
}

sub extract_crossmatches {
  local($l);
  if (length($_[0]) < 20) {
    $l = length($_[0]);
  } else {
    $l = 20;
  }
  $query = substr($_[0],0,$l);
  $hitindex = 0;
  $crossmatch = 0;
  $potentialcrossmatch = 0;
  open (BLASTOUT, "$blastoutput") || die ("could not open !");
  while (<BLASTOUT>) {
    if (/^>\w/) {
      $crossmatch = 0;
      $pothit = substr($_,0,$l);
      if ($pothit ne $query) {
	$potentialcrossmatch = 1;
      } else {
	$potentialcrossmatch = 0;
      }
    }
    if (/Identities = /) {
      if ($potentialcrossmatch == 1) {
	($ratio) = $_ =~ /((\d)+\/(\d)+)/;
	@fields = split(/\//,$ratio);
	$hitlength = $fields[1];
	if ($hitlength > $maxhitlength) {
	  $matches = $fields[0];
	  if (($matches / $hitlength) > $maxpercentage) {
	    $crossmatch = 1;
	    $hitread = -1;
	  }
	}
      }
    }
    if (/^Query:/) {
      if ($crossmatch == 1) {
	if ($hitread == -1) {
	  $hitread = 1;
	  ($temp) = $_ =~ /(^Query: (\d)+)/;
	  @fields = split(/\s+/,$temp);
	  $hitstart = $fields[1];
	  $hitstarts{$genome}{$gene}[$hitindex] = $hitstart-1;
	  $hitends{$genome}{$gene}[$hitindex] = $hitstart+$hitlength - 2;	
	  $hitindex++;
	}
      }
    }
  }
  close (BLASTOUT);
}
