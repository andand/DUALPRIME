#!/usr/bin/perl -w

##### GetGenePairs_2shared.pl #####

=head1 NAME

GetGenePairs_2shared.pl - Finds a large set of gene pairs sharing primer pairs

=head1 SYNOPSIS

GetGenePairs_2shared.pl -fasta1 genome1.fasta -fasta2 genome2.fasta -selfsim

=head1 DESCRIPTION

Based on the output from DesignPrimers_2shared.pl, as many non redundant dual-genome gene pairs as possible will be selected, sharing primer pairs. One primer pair per gene pair will be printed to the output file, which will include primer sequences and information. The user can choose to avoid primer pairs generationg intra-genomic unspecific amplicons. For this option to work, the program FindSelfSimilarity.pl must have been run in beforehand. The output file will be named genome1_genome2_2shared.tab

=head1 OPTIONS

(these are required)

-fasta1 <genome1.fasta> - gene sequences of genome 1 in fasta format

-fasta2 <genome2.fasta> - gene sequences of genome 2 in fasta format

(optional)

-selfsim            - The program avoids primer pairs generating unspecific amplicons. Requires that a file genome1_genome2_selfsim.tmp (or genome2_genome1_selfsim.tmp) exists in the current directory, i.e. that the script FindSelfSimilarity.pl has been run.

=head1 AUTHOR - Anders Andersson

Anders Andersson E<lt>ason@biotech.kth.seE<gt>

=cut

use Getopt::Long;

$infile1 = undef;
$infile2 = undef;
$help = undef;
$filterselfsim = undef;
&GetOptions('fasta1=s' => \$infile1, 'fasta2=s' => \$infile2, 'h!' => \$help, 'selfsim!' => \$filterselfsim);
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
$outfile = $prefix1."_".$prefix2."_2shared.tab";
$infile_primers = $prefix1."_".$prefix2."_2allprimers.tmp";
$infile_selfhits = $prefix1."_".$prefix2."_selfsim.tmp";
$ok = undef;
$ok = open (PRIMERFILE, $infile_primers);
if ($ok) {
  close (PRIMERFILE);
} else {
  $outfile = $prefix2."_".$prefix1."_2shared.tab";
  $infile_primers = $prefix2."_".$prefix1."_2allprimers.tmp";
  $infile_selfhits = $prefix2."_".$prefix1."_selfsim.tmp";
  $ok = open (PRIMERFILE, $infile_primers);
  if ($ok) {
    close (PRIMERFILE);
  } else {
    print"Couldn't find primer file!\n";
    die;
  }
}
srand;

#--------------- Main Program -----------------------
print"--- Selecting gene pairs sharing primer pairs ---\n";
&convert_file($infile1,$organism1);
&convert_file($infile2,$organism2);
if ($filterselfsim) {
  &read_selfblasthits;
}
&check_unspec_primers;
&read_primerfile;
&find_all_partners;
&select_maximum_partners;
&read_ORFseqs($organism1);
&read_ORFseqs($organism2);
&get_primer_info_selected_pairs;
&print_primers;
unlink($organism1,$organism2);
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

sub read_selfblasthits {
  open (HITSFILE, $infile_selfhits) || die ("could not open ! $!");
  local($nr) = 0;
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
    $genenr =  $fields[1];
    #$realname = $fields[2];
    $l = @fields;
    for ($i=3; $i<$l; $i++) {
      $hits{$fields[$i]} = 1;
    }
    $i=0;
    foreach $hit (keys %hits) {
      @fields = split(/_/,$hit);
      $hitstarts{$genome}{$genenr}[$i] = $fields[0]; 
      $hitends{$genome}{$genenr}[$i] = $fields[1];
      $i++;
    }
  }
  close(HITSFILE);
}

sub check_unspec_primers {
  %org_primer_gene = ();
  %unspec_primer = ();
  #print"check_unspec_primers\n";
  open(PRIMERFILE, $infile_primers) || die ("couldn't open $infile_primers");
  while(<PRIMERFILE>) {
    $row = $_;
    chomp $row;
    @fields = split(/\t/,$row);
    $primer1 = $fields[0];
    $primer2 = $fields[1];
    $genome = $fields[2];
    $gene = $fields[3];
    if (defined $org_primer_gene{$genome}{$primer1}) {
      if ($org_primer_gene{$genome}{$primer1} ne $gene) {
	$unspec_primer{$genome}{$primer1} = 1;
      }
    } else {
      $org_primer_gene{$genome}{$primer1} = $gene;
    }
    if (defined $org_primer_gene{$genome}{$primer2}) {
      if ($org_primer_gene{$genome}{$primer2} ne $gene) {
	$unspec_primer{$genome}{$primer2} = 1;
      }
    } else {
      $org_primer_gene{$genome}{$primer2} = $gene;
    }
  }
}

sub read_primerfile {
  %org1_primerpair_gene = ();
  %org2_primerpair_gene = ();
  %org1_gene_nrpartner = ();
  %org2_gene =();
  #print"read_primerfile\n";
  open(PRIMERFILE, $infile_primers) || die ("couldn't open $infile_primers");
  while(<PRIMERFILE>) {
    $row = $_;
    chomp $row;
    $discard = 0;
    if ($filterselfsim) {
      @fields = split(/\t/,$row);
      $genome = $fields[2];
      $gene = $fields[3];
      $first = $fields[4];
      $last = $fields[5];
      @fields = split(/,/,$first);
      $firststart = $fields[0];
      @fields = split(/,/,$last);
      $lastend = $fields[0];
      if (defined $hitstarts{$genome}{$gene}) {
	$nrselfhits = @{$hitstarts{$genome}{$gene}};
	for ($i=0; $i<$nrselfhits; $i++) {
	  $hitstart = $hitstarts{$genome}{$gene}[$i];
	  $hitend = $hitends{$genome}{$gene}[$i];
	  if ($hitstart < $firststart) {
	    $hitstart = $firststart;
	  }
	  if ($hitend > $lastend) {
	    $hitend = $lastend;
	  }
	  if (($hitend - $hitstart) > $maxOLhit) {
	    $discard = 1;
	    last;
	  }
	}
      }
    }
    if ($discard == 0) {
      @fields = split(/\t/,$row);
      if (defined $unspec_primer{$fields[2]}{$fields[0]}) {
	next;
      }
      if (defined $unspec_primer{$fields[2]}{$fields[1]}) {
	next;
      }
      $gene = $fields[3];
      if ($fields[2] eq $organism1) {
	$org1_gene_nrpartner{$gene} = 0;
	$org1_primerpair_gene{"$fields[0]\t$fields[1]"} = $gene;
      } else {
	$org2_gene{$gene} = 1;
	$org2_primerpair_gene{"$fields[0]\t$fields[1]"} = $gene;
      }
    }
  }
  $_ = (keys %org2_gene);
  #print"$_ nr of org2 genes primed in total\n";
  $_ = (keys %org1_gene_nrpartner);
  #print"$_ nr of org1 genes primed in total\n";
}

sub find_all_partners {
  %partners = ();
  %flip_partners = ();
  $flip = undef;
  #print"find_all_partners\n";
  foreach $primerpair (keys %org2_primerpair_gene) {
    if (defined $org1_primerpair_gene{$primerpair}) {
      if (defined $partners{ $org2_primerpair_gene{$primerpair} }{ $org1_primerpair_gene{$primerpair} }) {
	$partners{ $org2_primerpair_gene{$primerpair} }{ $org1_primerpair_gene{$primerpair} }++;
      } else {
	$partners{ $org2_primerpair_gene{$primerpair} }{ $org1_primerpair_gene{$primerpair} } = 1;
      }
    }
  }
  foreach $primerpair (keys %org1_primerpair_gene) {
    if (defined $org2_primerpair_gene{$primerpair}) {
      if (defined $flip_partners{ $org1_primerpair_gene{$primerpair} }{ $org2_primerpair_gene{$primerpair} }) {
	$flip_partners{ $org1_primerpair_gene{$primerpair} }{ $org2_primerpair_gene{$primerpair} }++;
      } else {
	$flip_partners{ $org1_primerpair_gene{$primerpair} }{ $org2_primerpair_gene{$primerpair} } = 1;
      }
    }
  }
  $_ = (keys %partners);
  #print"$_ partners..................\n";
  $_ = (keys %flip_partners);
  #print"$_ flip_partners.............\n";
  if ((keys %partners) > (keys %flip_partners)) {
    #print"Flipping\n";
    $flip = defined;
  }
}

sub get_primer_info_selected_pairs {
  #print"get_primer_info_selected_gene_pairs\n";
  open(PRIMERFILE, $infile_primers) || die ("couldn't open file!");
  while(<PRIMERFILE>) {
    chomp $_;
    @fields = split(/\t/);
    if (defined $unspec_primer{$fields[2]}{$fields[0]}) {
      next;
    }
    if (defined $unspec_primer{$fields[2]}{$fields[1]}) {
      next;
    }
    $primer1 = $fields[0];
    $primer2 = $fields[1];
    $genome = $fields[2];
    $genenr = $fields[3];
    $pos1 = $fields[4];
    $pos2 = $fields[5];
    @fields = split(/,/,$pos1);
    $start1 = $fields[0];
    $l1 = $fields[1];
    @fields = split(/,/,$pos2);
    $start2 = $fields[0];
    $l2 = $fields[1];
    if (defined $org2_primerpair_gene{"$primer1\t$primer2"}) {
      $org2gene = $org2_primerpair_gene{"$primer1\t$primer2"};
      if (defined $org1_primerpair_gene{"$primer1\t$primer2"}) {
	$org1gene = $org1_primerpair_gene{"$primer1\t$primer2"};
	if (defined $org2_gene_partner{$org2gene}) {
	  if ($org1gene eq $org2_gene_partner{$org2gene}) {
	    &get_primer_info;
	    push (@{$primerpair{$org2gene}}, "$primer1\t$primer2");
	  }
	}
      }
    }
  }
}

sub print_primers {
  open(OUTFILE, ">$outfile") || die ("couldn't open");
  ($genome1, $genome2) = ($organism1, $organism2);
  print OUTFILE"Organism1\t";
  print OUTFILE"Organism2\t";
  print OUTFILE"Org1 Gene\t";
  print OUTFILE"Org2 Gene\t";
  print OUTFILE"Primer1 (common)\t";
  print OUTFILE"Primer2 (common)\t";
  print OUTFILE"Org1 Primer1 Exact\t";
  print OUTFILE"Org1 Primer2 Exact\t";
  print OUTFILE"Org2 Primer1 Exact\t";
  print OUTFILE"Org1 Primer2 Exact\t";
  print OUTFILE"Primer1 Length\t";
  print OUTFILE"Primer2 Length\t";
  print OUTFILE"Org1 Primer1 tm\t";
  print OUTFILE"Org1 Primer2 tm\t";
  print OUTFILE"Org2 Primer1 tm\t";
  print OUTFILE"Org2 Primer2 tm\t";
  print OUTFILE"Average tm\t";
  print OUTFILE"Org1 Start Pos (in gene)\t";
  print OUTFILE"Org2 Start Pos (in gene)\t";
  print OUTFILE"Org1 PCR Length\t";
  print OUTFILE"Org2 PCR Length\n";
  foreach $org2gene (keys %org2_gene_partner) {
    $longest_pcrl = 0;
    while (@{$primerpair{$org2gene}} > 0) {
      $primerpair = shift(@{$primerpair{$org2gene}});
      $pcrl1 = $primerpair_PCRlength{$primerpair}{$genome1};
      $pcrl2 = $primerpair_PCRlength{$primerpair}{$genome2};
      if ($pcrl1 < $pcrl2) {
	$pcrl = $pcrl1;
      } else {
	$pcrl = $pcrl2;
      }
	if ($pcrl > $longest_pcrl) {
	  $longest_pcrl = $pcrl;
	  $longest_primerpair = $primerpair;
	}
    }
    @fields  = split(/\t/,$longest_primerpair);
    $primer1 = $fields[0];
    $primer2 = $fields[1];
    $avavtm = ($primer_tm{$primer1}{$genome1} + $primer_tm{$primer2}{$genome1} + $primer_tm{$primer1}{$genome2} + $primer_tm{$primer2}{$genome2})/4;
    $exact1_1 = $primer_exact{$primer1}{$genome1};
    $exact1_2 = $primer_exact{$primer1}{$genome2};
    $exact2_1 = $primer_exact{$primer2}{$genome1};
    $exact2_2 = $primer_exact{$primer2}{$genome2};
    $wobble1 = &make_wobbleprimer($exact1_1, $exact1_2);
    $wobble2 = &make_wobbleprimer($exact2_1, $exact2_2);
    print OUTFILE"$prefix1\t";
    print OUTFILE"$prefix2\t";
    print OUTFILE"$primer_genename{$primer1}{$genome1}\t";
    print OUTFILE"$primer_genename{$primer1}{$genome2}\t";
    print OUTFILE"$wobble1\t";#gemensam1
    print OUTFILE"$wobble2\t";#gemensam2
    print OUTFILE"$primer_exact{$primer1}{$genome1}\t";
    print OUTFILE"$primer_exact{$primer2}{$genome1}\t";
    print OUTFILE"$primer_exact{$primer1}{$genome2}\t";
    print OUTFILE"$primer_exact{$primer2}{$genome2}\t";
    print OUTFILE"$primer_l{$primer1}{$genome1}\t";
    print OUTFILE"$primer_l{$primer2}{$genome1}\t";
    print OUTFILE"$primer_tm{$primer1}{$genome1}\t";
    print OUTFILE"$primer_tm{$primer2}{$genome1}\t";
    print OUTFILE"$primer_tm{$primer1}{$genome2}\t";
    print OUTFILE"$primer_tm{$primer2}{$genome2}\t";
    print OUTFILE"$avavtm\t";
    print OUTFILE"$primerpair_PCRstart{$primerpair}{$genome1}\t";
    print OUTFILE"$primerpair_PCRstart{$primerpair}{$genome2}\t";
    print OUTFILE"$primerpair_PCRlength{$primerpair}{$genome1}\t";
    print OUTFILE"$primerpair_PCRlength{$primerpair}{$genome2}\n";
  }
  close(OUTFILE);
}

sub get_primer_info {
  $exakt1 = substr($geneseq{$genome}{$genenr}, $start1 , $l1);
  $gc1 = 0;
  for ($i=0; $i<$l1; $i++) {
    if (substr($exakt1, $i, 1) eq "G" || substr($exakt1, $i, 1) eq "C")  {
      $gc1++;
    }
  }
  $tm1 = 4*$gc1 + 2*($l1 - $gc1);
  $exakt2 = substr($geneseq{$genome}{$genenr}, ($start2-$l2+1), $l2);
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

sub read_ORFseqs {
  $genome = $_[0];
  $infile_ORFs = $_[0];
  open (ORFFILE, $infile_ORFs) || die ("could not open ! $!");
  $genenr = 0;
  while (<ORFFILE>) {
    $genenr++;
    chomp $_;
    @fields = split(/\t/);
    $seq = $fields[1];
    $geneseq{$genome}{$genenr} = $seq;
    $realname{$genome}{$genenr} = $fields[0];
  }
}


sub make_wobbleprimer {
  $nrmm=0;
  local($base1);
  local($base2);
  local($seq1) = $_[0]; 
  local($seq2) = $_[1];
  if (length($seq1) != length($seq2)) {
    die ("different primers to degenerate!");
  }
  local($l) = length($seq1); 
  $modseq = "";
  for ($i = 0; $i < $l; $i++) {
    $base1 = substr($seq1, $i, 1);
    $base2 = substr($seq2, $i, 1);
    if ($base1 ne $base2) {
      $nrmm++;
      if ($base1 eq "A") {
	if ($base2 eq "T") {
	  $modseq = $modseq."W";
	} elsif ($base2 eq "C") {
	  $modseq = $modseq."M";
	} elsif ($base2 eq "G") {
	  $modseq = $modseq."R";
	} else {
	  die ("unrecogniced base!");
	}
      } elsif ($base1 eq "T") {
	if ($base2 eq "A") {
	  $modseq = $modseq."W";
	} elsif ($base2 eq "C") {
	  $modseq = $modseq."Y";
	} elsif ($base2 eq "G") {
	  $modseq = $modseq."K";
	} else {
	  die ("unrecogniced base!");
	}
      } elsif ($base1 eq "C") {
	if ($base2 eq "A") {
	  $modseq = $modseq."M";
	} elsif ($base2 eq "T") {
	  $modseq = $modseq."Y";
	} elsif ($base2 eq "G") {
	  $modseq = $modseq."S";
	} else {
	  die ("unrecogniced base!");
	}
      } elsif ($base1 eq "G") {
	if ($base2 eq "A") {
	  $modseq = $modseq."R";
	} elsif ($base2 eq "T") {
	  $modseq = $modseq."K";
	} elsif ($base2 eq "C") {
	  $modseq = $modseq."S";
	} else {
	  die ("unrecogniced base!");
	}
      }
    } else {
      $modseq = $modseq.$base1;
    }
  }
  if ($nrmm > 2) {
    die ("more than 2 mismatches!");
  }
  return $modseq;
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


sub select_maximum_partners {
  local(%genes);
  local(%pairs);
  local(%path);
  local(%dist);
  %org2_gene_partner = ();
  foreach $org2gene (keys %partners) {
    $pairs{"source"}{$org2gene."\t".$organism2} = 1;
    $pairs{$org2gene."\t".$organism2}{"source"} = -1;
    $genes{$org2gene."\t".$organism2} = 1;
    foreach $org1gene (keys %{$partners{$org2gene}}) {
      $pairs{$org2gene."\t".$organism2}{$org1gene."\t".$organism1} = 1;
      $pairs{$org1gene."\t".$organism1}{$org2gene."\t".$organism2} = -1;
      $pairs{$org1gene."\t".$organism1}{"sink"} = 1;
      $pairs{"sink"}{$org1gene."\t".$organism1} = -1;
      $genes{$org1gene."\t".$organism1} = 1;
    }
  }
  $genes{"source"} = 1;
  $genes{"sink"} = 1;
  $nr_genes = (keys %genes);
  $path_found = 1;
  while ($path_found) {
    $path_found = undef;
    foreach $gene (keys %genes) {
      $dist{$gene} = undef;
    }
    $dist{"source"} = 1;
    %path = ();
    for ($i = 1; $i < ($nr_genes + 1); $i++) {
      foreach $gene (keys %genes) {
	if ($dist{$gene}) {
	  if ($dist{$gene} == $i) {
	    if ($gene eq "sink") {
	      #update residual graph
	      $us_gene = "sink";
	      while ($us_gene ne "source") {
		$ds_gene = $us_gene;
		$us_gene = $path{$ds_gene};
		$pairs{$us_gene}{$ds_gene} = -1 * $pairs{$us_gene}{$ds_gene};
		$pairs{$ds_gene}{$us_gene} = -1 * $pairs{$ds_gene}{$us_gene};
		$selected_sometime{$us_gene}{$ds_gene} = $pairs{$us_gene}{$ds_gene};
		$selected_sometime{$ds_gene}{$us_gene} = $pairs{$ds_gene}{$us_gene};
	      }
	      $path_found = 1;
	      goto PATHFOUND;
	    } else {
	      foreach $partner (keys %{$pairs{$gene}}) {
		if (!$dist{$partner}) {
		  if ($pairs{$gene}{$partner} == 1) {
		    $dist{$partner} = $i + 1;
		    $path{$partner} = $gene;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  PATHFOUND:
  }
  #transfer matchings from graf to partner genes
  foreach $gene (keys %selected_sometime) {
    if (($gene ne "source") && ($gene ne "sink")) {
      @gene_fields = split(/\t/,$gene);
      if ($gene_fields[1] eq $organism2) {
	foreach $partner (keys %{$selected_sometime{$gene}}) {
	  if ($selected_sometime{$gene}{$partner} == -1) {
	    if (defined $org2_gene_partner{$gene}) {
	      die ("\norg2gene with multiple partner genes!\n");
	    } else {
	      @partner_fields = split(/\t/,$partner);
	      $org2_gene_partner{$gene_fields[0]} = $partner_fields[0];
	    }
	  }
	}
      }
    }
  }
  $_ = keys (%org2_gene_partner);
}
