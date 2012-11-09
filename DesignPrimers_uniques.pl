#!/usr/bin/perl -w

##### DesignPrimers_uniques.pl #####

=head1 NAME

DesignPrimers_uniques - Designs unique primer pairs for one or two genomes

=head1 SYNOPSIS

DesignPrimers_2shared -fasta1 genome1.fasta -fasta2 genome2.fasta -selfsim -param parameterfile

=head1 DESCRIPTION

This script can be used in two ways; as a part of the DualPrime strategy for primer design or for designing primers for a single genome.

DualPrime: Allready primed genes will be determined by reading the file genome1_genome2_2shared.tab and genome1_genome2_1shared.tab (required). Among the remaing genes: for as many genes as possible in each of the two genomes given in the fastafiles, the program will design primers. The user can choose to avoid primer pairs generating intra-genome unspecific amplicons. For this option to work, the program FindSelfSimilarity.pl must have been run in beforehand. The output file will be named genome1_genome2_uniques.tab

Single genome: Prime pairs will be designed for each gene in the fasta file. The user can choose to avoid primer pairs generating intra-genome unspecific amplicons. For this option to work, the program FindSelfSimilarity.pl must have been run in beforehand. The output file will be named genome1_uniques.tab

=head1 OPTIONS

(these are required)

-fasta1 <genome1.fasta> - gene sequences of genome 1 in fasta format

(optional)

-fasta2 <genome2.fasta> - gene sequences of genome 2 in fasta format

-selfsim            - The program avoids primer pairs generating unspecific amplicons. Requires that a file genome1_genome2_selfsim.tmp (or genome1_selfsim.tmp) exists in the current directory

-param <parameterfile> - file with input parameters for Primer3. Easiest to change in the included parameterfile (parameter.txt). More parameters can be added in the end of the list according to the Primer3 v 0.9 readme (http://frodo.wi.mit.edu/primer3/primer3_code.html). Importantant not to change format of any parameter. Newline after last parameter. No empty lines at end. Without param default parameters (same as in in parameters.txt) will be used.

=head1 AUTHOR - Anders Andersson

Anders Andersson E<lt>ason@biotech.kth.seE<gt>

=cut

use Getopt::Long;

$infile1 = undef;
$infile2 = undef;
$help = undef;
$filterselfsim = undef;
$parameterfile = undef;
&GetOptions('fasta1=s' => \$infile1, 'fasta2=s' => \$infile2, 'h!' => \$help, 'selfsim!' => \$filterselfsim, 'param=s' => \$parameterfile);
if ($help or !$infile1) {
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
    $output_file = $prefix1."_".$prefix2."_uniques.tab";
    $infile_primedgenes1 = $prefix1."_".$prefix2."_2shared.tab";
    $infile_primedgenes2 = $prefix1."_".$prefix2."_1shared.tab";
    $infile_selfhits = $prefix1."_".$prefix2."_selfsim.tmp";
    $ok = open (INFILE, $infile_primedgenes2);
    if ($ok) {
      close (INFILE);
    } else {
      $output_file = $prefix2."_".$prefix1."_uniques.tab";
      $infile_primedgenes1 = $prefix2."_".$prefix1."_2shared.tab";
      $infile_primedgenes2 = $prefix2."_".$prefix1."_1shared.tab";
      $infile_selfhits = $prefix2."_".$prefix1."_selfsim.tmp";
      $ok = open (INFILE, $infile_primedgenes2);
      if ($ok) {
	close (INFILE);
      } else {
	print"Couldn't find one or both of files with shared primers!\n";
	die;
      }
    }
  }
} else {
  $output_file = $prefix1."_uniques.tab";
  $infile_selfhits = $prefix1."_selfsim.tmp";
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
print"--- Designing unique primers ---\n";
&convert_file($infile1,$organism1);
open(OUTFILE, ">$output_file") || die;
if ($parameterfile) {
  read_primer_parameters($parameterfile);
}
if ($organism2) {
  &convert_file($infile2,$organism2);
  &read_allready_primed($infile_primedgenes1);
  &read_allready_primed($infile_primedgenes2);
}
if ($filterselfsim) {
  &read_selfblasthits;
}
&print_header;
&run_uniques($organism1,$prefix1);
if ($organism2) {
  &run_uniques($organism2,$prefix2);
  unlink($organism2);
}
close(OUTFILE);
unlink($organism1);
unlink($inputp3,$outputp3);
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

sub read_allready_primed {
  $infile_primedgenes = $_[0];
  local($nr) = 0;
  open (PRIMEDFILE, $infile_primedgenes) || die ("couldn't open file!");
  while(<PRIMEDFILE>) {
    $nr++;
    if ($nr == 1) { next; }
    chomp $_;
    @fields = split(/\t/);
    local($genome1) = $fields[0]."_tab.tmp";
    local($genome2) = $fields[1]."_tab.tmp";
    local($genename1) = $fields[2];
    local($genename2) = $fields[3];
    $allready_primed{$genome1}{$genename1} = 1;
    $allready_primed{$genome2}{$genename2} = 1;
  }
}

sub read_selfblasthits {
  open (HITSFILE, $infile_selfhits) || die ("Could not open self hit file. $!");
  local($nr) = 0;
  local(%hits);
  while (<HITSFILE>) {
    $nr++;
    %hits = ();
    chomp $_;
    if ($nr == 1) {
      @fields =  split(/\t/);
      $maxOLhit = $fields[1];
      next;
    }
    @fields = split(/\t/);
    $genome = $fields[0];
    #$genenr =  $fields[1];
    $genename = $fields[2];
    $l = @fields;
    for ($i=3; $i<$l; $i++) {
      $hits{$fields[$i]} = 1;
    }
    $i=0;
    foreach $hit (keys %hits) {
      @fields = split(/_/,$hit);
      $hitstarts{$genome}{$genename}[$i] = $fields[0]; 
      $hitends{$genome}{$genename}[$i] = $fields[1];
      $i++;
    }
  }
  close(HITSFILE);
}

sub print_header {
  print OUTFILE"Organism\t";
  print OUTFILE"Gene\t";
  print OUTFILE"Primer1\t";
  print OUTFILE"Primer2\t";
  print OUTFILE"Primer1 Length\t";
  print OUTFILE"Primer2 Length\t";
  print OUTFILE"Primer1 tm\t";
  print OUTFILE"Primer2 tm\t";
  print OUTFILE"Average tm\t";
  print OUTFILE"Start Pos (in gene)\t";
  print OUTFILE"PCR Length\n";
}

sub run_uniques {
  $nr_madeit = 0;
  $nr_failed = 0;
  $genome = $infile_ORFs = $_[0];
  $prefix = $_[1];
  open (ORFFILE, $infile_ORFs) || die ("couldn't open file!");
  $genenr=0;
  while(<ORFFILE>) {
    $genenr++;
    chomp $_;
    @fields = split(/\t/);
    $realname{$genome}{$genenr} = $genename = $fields[0];
    $geneseq{$genome}{$genenr} = $geneseq = $fields[1];
    if (!defined $allready_primed{$infile_ORFs}{$genename}) {
      open (P3INPUT, ">$inputp3") || die;
      &print_input;
      close (P3INPUT);
      system ("primer3<$inputp3>$outputp3");
      @primers = ();
      &extract_primers;
      %hit = ();
      while(@primers) {
	$_ = @primers;
	$these_primers = shift @primers;
	($primer1,$primer2,$pos1,$pos2) = split(/:/,$these_primers);
	$primerpair = "$primer1\t$primer2";
	@fields = split(/,/,$pos1);
	$start1 = $fields[0];
	$l1 = $fields[1];
	@fields = split(/,/,$pos2);
	$start2 = $fields[0];
        $l2 = $fields[1];
	$_ = $start2 - $start1;
	$failed = 0;
	if ($filterselfsim) {
	  if (&filter_selfblast1($start1, $start2) == 1) {
	    $failed = 1;
	  }
	}
	if ($failed == 0) {
	  $nr_madeit++; $madeit{$genome}{$genenr} = 1; #print"made it: $nr_madeit\n";
	  &get_primer_info;
	  $avtm = ($primer_tm{$primer1}{$genome} + $primer_tm{$primer2}{$genome})/2;
	  print OUTFILE"$prefix\t";
	  print OUTFILE"$genename\t";
	  print OUTFILE"$primer1\t";
	  print OUTFILE"$primer2\t";
	  print OUTFILE"$primer_l{$primer1}{$genome}\t";
	  print OUTFILE"$primer_l{$primer2}{$genome}\t";
	  print OUTFILE"$primer_tm{$primer1}{$genome}\t";
	  print OUTFILE"$primer_tm{$primer2}{$genome}\t";
	  print OUTFILE"$avtm\t";
	  print OUTFILE"$primerpair_PCRstart{$primerpair}{$genome}\t";
	  print OUTFILE"$primerpair_PCRlength{$primerpair}{$genome}\n";
	  last;
	}
      }
      if (!defined $madeit{$genome}{$genenr}) {
	$nr_failed++; #print "failed: $nr_failed\n";
      }
    }
  }
}

sub filter_selfblast1 {
  local($start) = $_[0];
  local($end) = $_[1];
  if (defined $hitstarts{$genome}{$genename}) {
    $nrselfhits = @{$hitstarts{$genome}{$genename}};
    $discard = 0;
    for ($i=0; $i<$nrselfhits; $i++) {
      $hitstart = $hitstarts{$genome}{$genename}[$i];
      $hitend = $hitends{$genome}{$genename}[$i];
      if ($hitstart < $start) {
	$hitstart = $start;
      }
      if ($hitend > $end) {
	$hitend = $end;
      }
      if (($hitend - $start) > $maxOLhit) {
	$discard = 1;
	last;
      }
    }
    if ($discard == 0) {
      return 0;
    } else {
      return 1;
    }
  }
}

sub get_primer_info {
  $_ = $geneseq{$genome}{$genenr};
  $exakt1 = substr($geneseq{$genome}{$genenr}, $start1 , $l1);
  $gc1 = 0;
  for ($i=0; $i<$l1; $i++) {
    if (substr($exakt1, $i, 1) eq "G" || substr($exakt1, $i, 1) eq "C")  {
      $gc1++;
    }
  }
  $tm1 = 4*$gc1 + 2*($l1 - $gc1);
  $exakt2 = substr($geneseq{$genome}{$genenr}, ($start2-$l2+1), $l2);
  $_ = $geneseq{$genome}{$genenr};
  #OBS REVCOMP to exact sequence
  $gc2=0;
  for ($i=0; $i<$l2; $i++) {
    if (substr($exakt2, $i, 1) eq "G" || substr($exakt2, $i, 1) eq "C")  {
      $gc2++;
    }
  }
  $tm2 = 4*$gc2 + 2*($l2 - $gc2);
  $primer_l{$primer1}{$genome} = $l1;
  $primer_l{$primer2}{$genome} = $l2;
  $primer_exact{$primer1}{$genome} = $exakt1; 
  $primer_exact{$primer2}{$genome} = &make_revcomp($exakt2);
  $primer_tm{$primer1}{$genome} = $tm1;
  $primer_tm{$primer2}{$genome} = $tm2;
  $primerpair_PCRlength{"$primer1\t$primer2"}{$genome} = $start2 - $start1;
  $primerpair_PCRlength{"$primer2\t$primer1"}{$genome} = $start2 - $start1;
  $primer_genename{$primer1}{$genome} = $realname{$genome}{$genenr};
  $primer_genename{$primer2}{$genome} = $realname{$genome}{$genenr};
  $primerpair_PCRstart{"$primer1\t$primer2"}{$genome} = $start1;
  $primerpair_PCRstart{"$primer2\t$primer1"}{$genome} = $start1;
}

sub print_input {
  print P3INPUT "PRIMER_SEQUENCE_ID=$genome\t$genenr\n";
  print P3INPUT "SEQUENCE=$geneseq\n";
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
  print P3INPUT "PRIMER_NUM_RETURN=200\n";
  print P3INPUT "PRIMER_NUM_NS_ACCEPTED=2\n";
  print P3INPUT "=\n";
}

sub extract_primers {
  open (P3OUTPUT, $outputp3) || die ("could not open $!");
  local($n);
  local($these_primers);
  local($temp);
  local($pcrl);
  local($otherpcrl);
  local($i);
  while (<P3OUTPUT>) {
    chomp $_;
    #if (/PRIMER_SEQUENCE_ID=/) {
    #  @fields = split(/=/);
    #  $name = $fields[1];
    #}
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
	$these_primers = "$seq1:$seq2:$pos1:$pos2";
	@fields = split(/,/,$pos1);
	local($start) = $fields[0];
	@fields = split(/,/,$pos2);
	local($end) = $fields[0];
	$pcrl = $pcr_length{$these_primers} = $end - $start;
	push(@primers, $these_primers);
	$i = @primers;
	while ($i > 0) {
	  $otherpcrl = $pcr_length{$primers[$i-1]};
	  if ($pcrl > $otherpcrl) {
	    $temp = $primers[$i-1];
	    $primers[$i-1] = $primers[$i];
	    $primers[$i] = $temp;
	  }
	  $i--;
	}
      }
    }
  }
  close (P3OUTPUT);
}


sub make_revcomp {
  local($seq) = $_[0];
  local($base);
  local($i);
  local($l) = length($seq);
  $revcomp = "";
  for ($i=0; $i<$l; $i++) {
    $base = substr($seq, $i, 1);
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
      die ("strange base in revcomp!");
    }
  }
  return $revcomp;
}
