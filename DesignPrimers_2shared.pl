#!/usr/bin/perl -w

##### DesignPrimers_2shared.pl #####

=head1 NAME

DesignPrimers_2shared.pl - Designs primer pairs for each gene within regions shared by a gene of the other genome

=head1 SYNOPSIS

DesignPrimers_2shared.pl -fasta1 genome1.fasta -fasta2 genome2.fasta -param parameterfile

=head1 DESCRIPTION

For as many genes as possible in each of the two genomes given in the fastafiles, the program will design primers within regions shared by a gene of the other genome. Shared regions are localised by using the output file from FindSharedRegions.pl. The output file will be named genome1_genome2_2sharedprimers.tmp

=head1 OPTIONS

(these are required)

-fasta1 <genome1.fasta> - gene sequences of genome 1 in fasta format

-fasta2 <genome2.fasta> - gene sequences of genome 2 in fasta format

(optional)

-param <parameterfile> - file with input parameters for Primer3. Easiest to change in the included parameterfile (parameter.txt). More parameters can be added in the end of the list according to the Primer3 v 0.9 readme (http://frodo.wi.mit.edu/primer3/primer3_code.html). Importantant not to change format of any parameter. Newline after last parameter. No empty lines at end. Without param default parameters (same as in in parameters.txt) will be used.


=head1 AUTHOR - Anders Andersson

Anders Andersson E<lt>ason@biotech.kth.seE<gt>

=cut

use Getopt::Long;

$infile1 = undef;
$infile2 = undef;
$help = undef;
$selfsim = undef;
$parameterfile = undef;
&GetOptions('fasta1=s' => \$infile1, 'fasta2=s' => \$infile2, 'h!' => \$help, 'selfsim=s' => \$selfsim, 'param=s' => \$parameterfile);
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
$output_file = $prefix1."_".$prefix2."_2allprimers.tmp";
$infile_mers = $prefix1."_".$prefix2."_regions.tmp";
$ok = undef;
$ok = open (MERFILE, $infile_mers);
if ($ok) {
  close (MERFILE);
} else {
  $output_file = $prefix2."_".$prefix1."_2allprimers.tmp";
  $infile_mers = $prefix2."_".$prefix1."_regions.tmp";
  $ok = open (MERFILE, $infile_mers);
  if ($ok) {
    close (MERFILE);
  } else {
    print"Couldn't find shared regions file!\n";
    die;
  }
}
if ($parameterfile) {
  local($ok) = undef;
  $ok = open (PARAMFILE, $parameterfile);
  if ($ok) {
    close (PARAMFILE);
  } else {
    print"Couldn't find primer parameter file: $parameterfile!\n\n";
    die;
  }
}
$inputp3 = "inputp3.tmp";
$outputp3 = "outputp3.tmp";

#--------------- Main Program -----------------------
print"--- Designing shared primer pairs ---\n";
&convert_file($infile1,$organism1);
&convert_file($infile2,$organism2);
open(OUT, ">$output_file") || die;
if ($parameterfile) {
  read_primer_parameters($parameterfile);
}
&get_gene_lengths($organism1);
&get_gene_lengths($organism2);
&find_multipels;
&map_mers;
&merge_mers($organism1);
&merge_mers($organism2);
close(OUT);
unlink($organism1,$organism2);
unlink($inputp3,$outputp3);
#----------------------------------------------------
#$infile_db\tgene_no\tposition\t$otherwindow\tdirection\t$infile_query\t$nr\t$i\t$window\n;
#   0            1       2          3            4              5        6    7     8

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

sub read_primer_parameters {
  $infile_primerparam = $_[0];
  $parameters = "";
  open (PARAMFILE, $infile_primerparam) || die ("couldn't open the file $infile_primerparam");
  while(<PARAMFILE>) {
    $row = $_;
    $parameters = $parameters.$row;
  }
  close (PARAMFILE);
}

sub get_gene_lengths {
  $infile_ORFs = $_[0];
  open (ORFFILE, $infile_ORFs) || die ("could not open ! $!");
  $genenr=0;
  while (<ORFFILE>) {
    $genenr++;
    chomp $_;
    @fields = split(/\t/);
    $geneseq = $fields[1];
    $genel = length($geneseq);
    $gene_length{$infile_ORFs}{$genenr} = $genel;
  }
}

sub find_multipels {
  #print"finding multipels\n";
  open (MERFILE, $infile_mers) || die ("could not open ! $!");
  while(<MERFILE>) {
    chomp $_;
    @fields = split(/\t/);
    if (defined $allready{$fields[0]}{$fields[1]}{$fields[2]}) {
      $multipel{$fields[0]}{$fields[1]}{$fields[2]} = 1;
    } else {
      $allready{$fields[0]}{$fields[1]}{$fields[2]} = 1;
    }
    if ($fields[4] eq "for") {
      if (defined $allready{$fields[5]}{$fields[6]}{$fields[7]}) {
	$multipel{$fields[5]}{$fields[6]}{$fields[7]} = 1;
      } else {
	$allready{$fields[5]}{$fields[6]}{$fields[7]} = 1;
      }
    } else {
      $position = ($gene_length{$fields[5]}{$fields[6]} - $fields[7]);
      if (defined $allready{$fields[5]}{$fields[6]}{$position}) {
	$multipel{$fields[5]}{$fields[6]}{$position} = 1;
      } else {
	$allready{$fields[5]}{$fields[6]}{$position} = 1;
      }
    }
  }
  %allready = ();
  close (MERFILE);
}

#$infile_db\tgene_no\tposition\t$otherwindow\tdirection\t$infile_query\t$nr\t$i\t$window\n;
#   0         1         2          3             4              5        6    7     8

sub map_mers {
  #print "mapping mers\n";
  open (MERFILE, $infile_mers) || die ("could not open ! $!");
  while(<MERFILE>) {
    chomp $_;
    @fields = split(/\t/);
    #map search mer onto genome1
    $genome = $fields[0];
    $genenr = $fields[1];
    $position = $fields[2];
    $genepartnernr = "$fields[4]\t$fields[6]";
    if (!defined $multipel{$genome}{$genenr}{$position}) {
      for ($i=0; $i<18; $i++) {
	if (substr($fields[3], $i, 1) ne substr($fields[8],$i,1)) {
	  $mm = $position + $i;
	  $mmpos{$genome}{$genenr}{$genepartnernr}{$mm} = 1;
	  #obs specialfall: inga mm!
	}
      }
      push (@{$hompos{$genome}{$genenr}{$genepartnernr}}, $position);
    }
    #map hit mer onto genome2
    $genome = $fields[5];
    $genenr = $fields[6];
    $genepartnernr = "$fields[4]\t$fields[1]";
    if ($fields[4] eq "rev") {
      $position = ($gene_length{$fields[5]}{$fields[6]} - $fields[7]);
      if (!defined $multipel{$genome}{$genenr}{$position}) {
	for ($i=0; $i<18; $i++) {
	  if (substr($fields[3], $i, 1) ne substr($fields[8],$i,1)) {
	    $mm = $position - $i;
	    $mmpos{$genome}{$genenr}{$genepartnernr}{$mm} = 1;
	  }
	}
	push (@{$hompos{$genome}{$genenr}{$genepartnernr}}, $position);
      }
    } else {
      $position = $fields[7];
      if (!defined $multipel{$genome}{$genenr}{$position}) {
	for ($i=0; $i<18; $i++) {
	  if (substr($fields[3], $i, 1) ne substr($fields[8],$i,1)) {
	    $mm = $position + $i;
	    $mmpos{$genome}{$genenr}{$genepartnernr}{$mm} = 1;
	  }
	}
	push (@{$hompos{$genome}{$genenr}{$genepartnernr}}, $position);
      }
    }
  }
  %multipel = ();
  close(MERFILE);
}

sub merge_mers {
  #print "mnm\n";
  $infile_ORFs = $_[0];
  open (ORFFILE, $infile_ORFs) || die ("could not open ! $!");
  $genome = $infile_ORFs;
  $genenr = 0;
  $nrregions = 0;
  while (<ORFFILE>) {
    $genenr++; 
    #print"$nr\n";
    chomp $_;
    @fields = split(/\t/);
    #$realname = $fields[0];
    $geneseq = $fields[1];
    $genel = length($geneseq);
    #@mergedpos = (); Nätverks varianten
    #@mergedmm = (); Nätverks varianten
    #$nrregions = 0; Nätverks varianten
    if (defined $hompos{$genome}{$genenr}) {
      foreach $genepartnernr (keys %{$hompos{$genome}{$genenr}}) {
	@mergedpos = (); #Ej i nätverks varianten
	@mergedmm = (); #Ej i nätverks varianten
	$nrregions = 0; #Ej i nätverks varianten
	@tempmerged = ();
	@{$sortedpos{$genome}{$genenr}{$genepartnernr}} = sort by_number (@{$hompos{$genome}{$genenr}{$genepartnernr}});
	$l = @{$sortedpos{$genome}{$genenr}{$genepartnernr}};
	$tempmerged[0] = $a = $b = $sortedpos{$genome}{$genenr}{$genepartnernr}[0];
	$j = 0;
	for ($i=1; $i<$l; $i++) {
	  $b = $sortedpos{$genome}{$genenr}{$genepartnernr}[$i];
	  if (($b - $a) > 18) {
	    $tempmerged[$j] = $tempmerged[$j]."\t$a";
	    $j++;
	    $tempmerged[$j] = $b;
	  }
	  $a = $b;
	}
	if ($a == $b) {
	  $tempmerged[$j] = $tempmerged[$j]."\t$b";
	}
	$l = (@tempmerged);
	for ($i = 0; $i < $l; $i++) {
	  @fields = split(/\t/,$tempmerged[$i]);
	  $mm = 0;
	  if (defined $mmpos{$genome}{$genenr}{$genepartnernr}) {
	    foreach $key (keys %{$mmpos{$genome}{$genenr}{$genepartnernr}}) {
	      if (($fields[0] <= $key) && (($fields[1]+18) > $key)) {
		$mm++;
		$mergedmm[$i+$nrregions]{$key} = 1;
	      }
	    }
	  }
	}
	$nrregions = $nrregions + $l;
	@mergedpos = (@mergedpos, @tempmerged);
	#print"$infile_ORFs\t$genenr\t$genepartnernr\t$nrregions\n"; For evaluation
	if ($nrregions > 1) { #Ej i nätverks varianten
	  &combineNrun; #Ej i nätverks varianten
	} #Ej i nätverks varianten
      }
      #&combineNrun; Här ska den va i nätverksvarianten
    }
  }
  close(ORFFILE);
}

sub combineNrun {
  open (P3INPUT, ">$inputp3") || die;
  $nrregions = @mergedpos;
  $nrrecords = 0;
  for ($i=0; $i<$nrregions; $i++) {
    $mutseq1 = $geneseq;
    @fields = split(/\t/,$mergedpos[$i]);
    $start1 = $fields[0];
    $end1 = $fields[1] + 18;
    foreach $mm (keys %{$mergedmm[$i]}) {
      substr($mutseq1, $mm, 1) = "N";
    }
    for ($j=0; $j<$i; $j++) {
      @fields = split(/\t/,$mergedpos[$j]);
      $start2 = $fields[0];
      $end2 = $fields[1] + 18;
      $go = 0;
      if (($end2-$start1 >= 0)) {
	$firststart = $start1;
	$firstend = $end1;
	$laststart = $start2;
	$lastend = $end2;
	$go = 1;
      } elsif (($end1-$start2 >= 0)) {
	$firststart =$start2;
	$firstend = $end2;
	$laststart = $start1;
	$lastend = $end1;
	$go = 1;
      }
      if ($go == 1) {
	$mutseq2 = $mutseq1;
	$nrrecords++;
	if(($nrrecords % 10000) == 9999) {
	  close (P3INPUT);
	  system ("primer3<$inputp3>$outputp3");
	  &extract_primers;
	  open (P3INPUT, ">$inputp3") || die;
	}
	$excluded = "0,$firststart ".$firstend.",".($laststart-$firstend)." ".($lastend).",".($genel-($lastend));
	foreach $mm (keys %{$mergedmm[$j]}) {
	  substr($mutseq2, $mm, 1) = "N";
	}
	&print_input;
      }
    }
  }
  #A record without excluded regions: (if the gene sequence is identical to another gene)
  $excluded = "";
  &print_input;
  close (P3INPUT);
  system ("primer3<$inputp3>$outputp3");
  &extract_primers;
}

sub by_number {
  if ($a < $b) {
    -1;
  } elsif ($a == $b) {
    0;
  } elsif ($a > $b) {
    1;
  }
}

sub print_input {
  print P3INPUT "PRIMER_SEQUENCE_ID=$genome\t$genenr\n";
  print P3INPUT "SEQUENCE=$mutseq2\n";
  print P3INPUT "EXCLUDED_REGION=$excluded\n";
  if ($parameterfile) {
    print P3INPUT "$parameters";
  } else {
    print P3INPUT "PRIMER_OPT_SIZE=19\n";
    print P3INPUT "PRIMER_MIN_SIZE=18\n";
    print P3INPUT "PRIMER_MAX_SIZE=22\n";
    print P3INPUT "PRIMER_MIN_TM=48.0\n";
    print P3INPUT "PRIMER_OPT_TM=52.0\n";
    print P3INPUT "PRIMER_MAX_TM=65.0\n";
    print P3INPUT "PRIMER_MAX_DIFF_TM=6.0\n";
    #print P3INPUT "PRIMER_SELF_ANY=8\n";#default=8
    print P3INPUT "PRIMER_MAX_POLY_X=4\n";
    print P3INPUT "PRIMER_PRODUCT_SIZE_RANGE=100-800\n";
    #print P3INPUT "PRIMER_WT_TM_GT=0\n";
    #print P3INPUT "PRIMER_WT_TM_LT=0\n";
    #print P3INPUT "PRIMER_WT_SIZE_GT=0\n";
    #print P3INPUT "PRIMER_WT_SIZE_LT=0\n";
    #print P3INPUT "PRIMER_WT_GC_PERCENT_GT=0\n";
    #print P3INPUT "PRIMER_WT_GC_PERCENT_LT=0\n";
    #print P3INPUT "PRIMER_WT_NUM_NS=100\n";
    #print P3INPUT "PRIMER_PAIR_WT_DIFF_TM=1\n";
    #print P3INPUT "PRIMER_PAIR_WT_COMPL_END=1\n";
    #print P3INPUT "PRIMER_PAIR_WT_COMPL_ANY=1\n";
  }
  print P3INPUT "PRIMER_NUM_NS_ACCEPTED=2\n";
  print P3INPUT "PRIMER_NUM_RETURN=50\n";
  print P3INPUT "=\n";
}

sub extract_primers {
  open (P3OUTPUT, "$outputp3") || die ("could not open $!");
  local($n);
  while (<P3OUTPUT>) {
    chomp $_;
    if (/PRIMER_SEQUENCE_ID=/) {
      @fields = split(/=/);
      $name = $fields[1];
    }
    if (/PRIMER_LEFT(_(\d)+)?_SEQUENCE=/) {
      @fields = split(/=/);
      $seq1 = $fields[1];
      $n = 0;
      if (substr($seq1,-1,1) eq "N") {
	$n = 1;
      }
    }
    if (/PRIMER_RIGHT(_(\d)+)?_SEQUENCE=/) {
      @fields = split(/=/);
      $seq2 = $fields[1];
      if (substr($seq2,-1,1) eq "N") {
	$n = 1;
      }
    }
    if (/PRIMER_LEFT(_(\d)+)?=/) {
      @fields = split(/=/);
      $pos1 = $fields[1];
    }
    if (/PRIMER_RIGHT(_(\d)+)?=/) {
      @fields = split(/=/);
      $pos2 = $fields[1];
      if ($n == 0) {
	print OUT "$seq1\t$seq2\t$name\t$pos1\t$pos2\n";
      }
    }
  }
  close (P3OUTPUT);
}

sub make_revcomp {
  local($seq) = $_[0];
  local($l) = length($seq);
  local($base);
  local($i);
  $revcomp = "";
  for ($i=0; $i<$l; $i++) {
    $base = substr($seq, $i, 1);
    if ($base eq "A") {
      $revcomp = "T".$revcomp;
    }
    if ($base eq "T") {
      $revcomp = "A".$revcomp;
    }
    if ($base eq "C") {
      $revcomp = "G".$revcomp;
    }
    if ($base eq "G") {
      $revcomp = "C".$revcomp;
    }
    if ($base eq "N") {
      $revcomp = "N".$revcomp;
    }
  }
  $revcomp;
}
