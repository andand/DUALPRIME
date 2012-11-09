#!/usr/bin/perl -w

##### DualPrime.pl #####

=head1 NAME

DualPrime.pl - Designs shared primers for the genes of two genomes

=head1 SYNOPSIS

DualPrime.pl -fasta1 genome1.fasta -fasta2 genome2.fasta -selfsim -l regionlength -s percentage similarity -param parameterfile

=head1 DESCRIPTION

For the genes of the two input fasta files, DualPrime will design primer pairs belonging to either of three categories: Both primers shared between a pair of genes from both genomes, one primer shared between a pair of genes from both genomes, both primers unique. The output files for each primer category are called: genome1_genome2_2shared.tab, genome1_genome2_1shared.tab and genome1_genome2_uniques.tab. During the process several temporary files (*.temp) will be created (necessary for the process) that will automatically be erased when the process is finished.

If the option -selfsim is used, output primers are filtered such that non gene-specific amplicons are avoided: Amplicons with regions similar to other genes in the same genome, exceeding the length (number of nucleotides) and similarity (percentage identical nucleotides) given in the input (-l and -s) will be filtered out. If no values are given on -l and -s, the default values 30 nucleotides and 66% similarity will be used.

If the option -param is used, input parameters to Primer3 (min and max primer GC-content, length, tm etc.) can be specified in the input parameter file. This is easiest accomplished by changing in the provided file parameters.txt, and saving this in a new name. It is importantant not to change format of any parameters. More parameters can be added according to the Primer3 v 0.9 readme (http://frodo.wi.mit.edu/primer3/primer3_code.html). Newline after last parameter. No empty lines at end. Without -param, default parameters (same as in in parameters.txt) will be used.

=head1 OPTIONS

(these are required)

-fasta1 <genome1.fasta> - gene sequences of genome 1 in fasta format

-fasta2 <genome2.fasta> - gene sequences of genome 2 in fasta format

(these are optional)

-selfsim   - The program avoids primer pairs generating unspecific amplicons.

-l         - minimum value for region length is 30 (default 30)

-s         - minimum and maximum value for similarity is 66 and 99, respectively (default 66)

-param <parameterfile> - file with input parameters for Primer3 (see above)

=head1 AUTHOR - Anders Andersson

Anders Andersson E<lt>ason@biotech.kth.seE<gt>

=cut

use Getopt::Long;

$infile1 = undef;
$infile2 = undef;
$help = undef;
$filterselfsim = undef;
$maxhitlength = undef;
$maxpercentage = undef;
$parameterfile = undef;
&GetOptions('fasta1=s' => \$infile1, 'fasta2=s' => \$infile2, 'h!' => \$help, 'selfsim!' => \$filterselfsim, 'l=i' => \$maxhitlength, 's=i' => \$maxpercentage, 'param=s' => \$parameterfile);
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
if ($maxhitlength and !$filterselfsim) {
  system ('perldoc', $0);
  exit;
}
if ($maxpercentage and !$filterselfsim) {
  system ('perldoc', $0);
  exit;
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
if (!$maxhitlength) {
  $maxhitlength = 30;
}
if (!$maxpercentage) {
  $maxpercentage = 66;
}
if (($maxhitlength < 30) or ($maxpercentage < 66) or ($maxpercentage > 99)) {
  system ('perldoc', $0);
  exit;
}

print"\n### Running the DualPrime scripts ###\n";
if (defined $filterselfsim) {
  print"1 of 7 ";
  $status = system ("perl FindSharedRegions.pl -fasta1 $infile1 -fasta2 $infile2"); &crash if ($status != 0);
  print"2 of 7 ";
  $status = system ("perl FindSelfSimilarity.pl -fasta1 $infile1 -fasta2 $infile2 -l $maxhitlength -s $maxpercentage"); &crash if ($status != 0);
  print"3 of 7 ";
  if (defined $parameterfile) {
    $status = system ("perl DesignPrimers_2shared.pl -fasta1 $infile1 -fasta2 $infile2 -param $parameterfile"); &crash if ($status != 0);
  } else {
    $status = system ("perl DesignPrimers_2shared.pl -fasta1 $infile1 -fasta2 $infile2"); &crash if ($status != 0);
  }
  print"4 of 7 ";
  $status = system ("perl GetGenePairs_2shared.pl -fasta1 $infile1 -fasta2 $infile2 -selfsim"); &crash if ($status != 0);
  print"5 of 7 ";
  if (defined $parameterfile) {
    $status = system ("perl DesignPrimers_1shared.pl -fasta1 $infile1 -fasta2 $infile2 -param $parameterfile"); &crash if ($status != 0);
  } else {
    $status = system ("perl DesignPrimers_1shared.pl -fasta1 $infile1 -fasta2 $infile2"); &crash if ($status != 0);
  }
  print"6 of 7 ";
  $status = system ("perl GetGenePairs_1shared.pl -fasta1 $infile1 -fasta2 $infile2 -selfsim"); &crash if ($status != 0);
  print"7 of 7 ";
  if (defined $parameterfile) {
    $status = system ("perl DesignPrimers_uniques.pl -fasta1 $infile1 -fasta2 $infile2 -selfsim -param $parameterfile"); &crash if ($status != 0);
  } else {
    $status = system ("perl DesignPrimers_uniques.pl -fasta1 $infile1 -fasta2 $infile2 -selfsim"); &crash if ($status != 0);
  }
  unlink($prefix1."_".$prefix2."_regions.tmp",$prefix1."_".$prefix2."_2allprimers.tmp",$prefix1."_".$prefix2."_1allprimers.tmp",$prefix1."_".$prefix2."_selfsim.tmp");
  print"\n### Run Completed ###\n";
} else {
  print"1 of 6 ";
  $status = system ("perl FindSharedRegions.pl -fasta1 $infile1 -fasta2 $infile2"); &crash if ($status != 0);
  print"2 of 6 ";
  if (defined $parameterfile) {
    $status = system ("perl DesignPrimers_2shared.pl -fasta1 $infile1 -fasta2 $infile2 -param $parameterfile"); &crash if ($status != 0);
  } else {
    $status = system ("perl DesignPrimers_2shared.pl -fasta1 $infile1 -fasta2 $infile2"); &crash if ($status != 0);
   }
  print"3 of 6 ";
  $status = system ("perl GetGenePairs_2shared.pl -fasta1 $infile1 -fasta2 $infile2"); &crash if ($status != 0);
  print"4 of 6 ";
  if (defined $parameterfile) {
    $status = system ("perl DesignPrimers_1shared.pl -fasta1 $infile1 -fasta2 $infile2 -param $parameterfile"); &crash if ($status != 0);
  } else {
    $status = system ("perl DesignPrimers_1shared.pl -fasta1 $infile1 -fasta2 $infile2"); &crash if ($status != 0);
  }
  print"5 of 6 ";
  $status = system ("perl GetGenePairs_1shared.pl -fasta1 $infile1 -fasta2 $infile2"); &crash if ($status != 0);
  print"6 of 6 ";
  if (defined $parameterfile) {
    $status = system ("perl DesignPrimers_uniques.pl -fasta1 $infile1 -fasta2 $infile2 -param $parameterfile"); &crash if ($status != 0);
  } else {
    $status = system ("perl DesignPrimers_uniques.pl -fasta1 $infile1 -fasta2 $infile2"); &crash if ($status != 0);
  }
  unlink($prefix1."_".$prefix2."_regions.tmp",$prefix1."_".$prefix2."_2allprimers.tmp",$prefix1."_".$prefix2."_1allprimers.tmp");
  print"\n### Run Completed ###\n";
}

sub crash {
  unlink($prefix1."_tab.tmp",$prefix2."_tab.tmp",$prefix1."_".$prefix2."_regions.tmp",$prefix1."_".$prefix2."_2allprimers.tmp",$prefix1."_".$prefix2."_1allprimers.tmp",$prefix1."_".$prefix2."_selfsim.tmp");
  die ("\nRun aborted. Primer design not completed. Error occured");
}







