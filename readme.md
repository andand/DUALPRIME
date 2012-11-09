DualPrime

Anders Andersson (2004)

CITATION
____________
Dual-genome primer design for construction of DNA microarrays.
Andersson A, Bernander R, Nilsson P.
Bioinformatics. 2005 Feb 1;21(3):325-32. Epub 2004 Aug 27.


INTRODUCTION
------------
DualPrime selects PCR primers for construction of gene specific tag (GST) microarrays. It tries to minimize the numbers of primers needed to prime all genes of two genomes (or other two lists of genes).

For as many dual-genome pairs of genes as possible, common primer pairs will be selected. Among the remaining genes, as many dual-genome pairs as possible sharing single primers will be selected. Finally, unique primers will be used for the remaining genes. 

To assure high gene-specificity of the resulting probes (amplicons), only those that do not have high sequence similarity to another gene within the same genome will be selected. If no specific amplicon is found for a gene it will remain unprimed.

The primer design is performed by Primer3 (http://www-genome.wi.mit.edu/genome_software/other/primer3.html) and the specificity test is conducted using the software BLAST (http://www.ncbi.nlm.nih.gov/BLAST/). Both of these need to be properly installed on the workstation.

DualPrime can be run on Linux/Unix platforms.

Input is two files with gene sequences in fasta format. Output is three tabulator separated text files, one for each category of primers, wich can be viewed with a text editor or with a spread sheet editor.

Requirements regarding primer gc contents, lengths etc. can be provided in a primer parameter file.


INSTALLATION
------------
DualPrime requires that Perl 5.8 or later (URL: http://www.perl.org), Primer3 (URL: http://www-genome.wi.mit.edu/genome_software/other/primer3.html) and BLAST (URL: http://www.ncbi.nlm.nih.gov/BLAST/) are installed, and that search paths are properly set to find these programs (see installation instructions for the individual programs).

Download the compressed file DualPrime.tgz. Extract the files by $ tar -zxvf DualPrime.tgz. A folder named DualPrime with a set of perl scripts (.pl), two .txt files and two .fasta files is created.  

If Perl is not installed in the directory /usr/bin/perl, the first row of each script (.pl files) need to be changed to the proper directory.

BuchneraSg.fasta and BuchneraBp.fasta are test input files, suitable for a test run of the program. parameters.txt is an example of a primer parameter file. All other files (with pl extensions) are required for running the program.


RUNNING THE PROGRAM
-------------------
The simplest way to run the program is to run the script DualPrime.pl. This script excecutes the other scripts automatically, rendering three output files with primer sequences and information on the primers.
 
Copy the two fasta files with gene sequences (here called Genome1.fasta and Genome2.fasta) to the DualPrime folder.

Change to the DualPrime directory. In the terminal window, write:

   perl DualPrime.pl -fasta1 Genome1.fasta -fasta2 Genome2.fasta
   
The option -selfsim can be used to avoid primer pairs generating amplicons with sequence similarities to other genes than the intended. The default cut-off values are (more than) 66% identity over (more than) 30 nucleotides (without gaps). However, by setting the parameters -l and/or -s the tolerance can be increased (but not decreased below the default values). For instance, by typing:
	 
      perl DualPrime.pl -fasta1 Genome1.fasta -fasta2 Genome2.fasta -selfsim -l 50 -s 99

amplicons with more than 99% identity over more than 50 nucleotides will be avoided.

The option -param <parameterfile> can be used to specify parameters regarding primer quality, as gc-content, length etc. Required parameters are given in the input parameterfile. It's easiest to change in the provided parameterfile parameters.txt, and giving this as parameterfile. More parameters can be added according to the Primer3 v 0.9 readme (http://frodo.wi.mit.edu/primer3/primer3_code.html). No extra line after last parameter.



To get help, type:

   perl DualPrime.pl? or perl DualPrime.pl h


Output
------
Three output files will be created in the DualPrime folder during the process:

    Genome1_Genome2_2shared.tab (both primers in a primer pair are shared between a pair of genes)
    Genome1_Genome2_1shared.tab (one primer in a primer pair is shared between a pair of genes) 
    Genome1_Genome2_uniques.tab (non of the primers are shared)

Genome1_Genome2_2shared.tab:	On each row is sequences and information for a primer pair targeting a gene in each genome. Information includes gene names, wobble primer sequences (including a maximum of two mixed (IUB code) positions), genome specific exact sequences, lengths, melting temperatures (calculated with the 2 for A/T + 4 for G/C role), in-gene start positions of amplicons, and lengths of amplicons. The first row includes descriptive headers.

Genome1_Genome2_1shared.tab:	On each row is sequences and information for three primers (two primer pairs with one common primer) targeting a gene in each genome. Information includes gene names, wobble primer sequence (including a maximum of two mixed (IUB code) positions), genome specific exact sequences, lengths, melting temperatures, in-gene start positions of amplicons, and lengths of amplicons. The first row includes descriptive headers.

Genome1_Genome2_uniques.tab:	On each row is sequences and information for a primer pair targeting a gene in on of the genomes. Inforamtion includes gene name, primer sequences, lengths, melting temperatures, in-gene start positions of amplicon, and lengths of amplicon. The first row includes descriptive headers.	


Output on screen will be:

       ### Running the DualPrime scripts ###
       1 of 6 --- Finding shared regions ---	
       2 of 6 --- Designing shared primer pairs ---
       3 of 6 --- Selecting gene pairs sharing primer pairs --- 
       4 of 6 --- Designing shared primers ---
       5 of 6 --- Selecting gene pairs sharing one primer ---
       6 of 6 --- Designing unique primers ---
       ### Run Completed ###

or, if -selfsim is used:

    ### Running the DualPrime scripts ###
    1 of 7 --- Finding shared regions ---
    2 of 7 --- Finding self similarity ---
    3 of 7 --- Designing shared primer pairs ---
    4 of 7 --- Selecting gene pairs sharing primer pairs --- 
    5 of 7 --- Designing shared primers ---
    6 of 7 --- Selecting gene pairs sharing one primer ---
    7 of 7 --- Designing unique primers ---
    ### Run Completed ###


During a run, several *.tmp files will be created, necessary for the process. These will be deleted before completion of the run. If the run is aborted they will, however, remain, and may be deleted manually.


To test DualPrime
-----------------
Use the (relatively small) test input files BuchneraSg.fasta and BuchneraBp.fasta:

    perl DualPrime.pl -fasta1 Genome1.fasta -fasta2 Genome2.fasta -selfsim

Using these input files, the process takes approximatly 1h to complete on with a Pentium 4 processor, where the first script, which finds shared regions, is the most time consuming.


Running the scripts individually
--------------------------------
In som cases it is advantagous to run the scripts one by one, for instance if one wants to test how different values on -s and -l influence the number of primed genes. It is then inconvenient to re-run the whole process from the start. By running the scripts individually, the *.tmp files (with information on shared regions etc.), used by subsequent scripts in the flow-scheme, will not be erased, and each step can be run several times using different parameters. This is the order to run the scripts when the -selfsim option is to be used:

FindSharedRegions.pl	        (Finds sequence similarities between the two genomes)       
DesignPrimers_2shared.pl	(Designs primer pairs, both primers within shared regions)    
FindSelfSimilarity.pl		(Finds sequence similarities within the genomes. Cut-off can be specified, or default is used)
GetGenePairs_2shared.pl		(Selects gene pairs sharing primer pairs. Option to avoid unsecific amplicons)  
DesignPrimers_1shared.pl	(Designs primer pairs, one primer within shared regions)  
GetGenePairs_1shared.pl		(Selects gene pairs sharing single primers. Option to avoid unsecific amplicons)  
DesignPrimers_uniques.pl	(Designs and selects unique primer pairs. Option to avoid unsecific amplicons)

To change cut-off values on self similarity, FindSelfSimilarity.pl is re-run with new parameters, followed by running subsequent scripts.


Inputs required (optional in parenthesis) for the different scripts:

FindSharedRegions.pl	        -fasta1 Genome1.fasta -fasta2 Genome2.fasta       
DesignPrimers_2shared.pl	-fasta1 Genome1.fasta -fasta2 Genome2.fasta
FindSelfSimilarity.pl		-fasta1 Genome1.fasta -fasta2 Genome2.fasta (-l number -s number)*
GetGenePairs_2shared.pl		-fasta1 Genome1.fasta -fasta2 Genome2.fasta (-selfsim))**		
DesignPrimers_1shared.pl	-fasta1 Genome1.fasta -fasta2 Genome2.fasta***
GetGenePairs_1shared.pl		-fasta1 Genome1.fasta -fasta2 Genome2.fasta (-selfsim)**	  
DesignPrimers_uniques.pl	-fasta1 Genome1.fasta -fasta2 Genome2.fasta**** (-selfsim)**

*Default: -l 30 -s 66
**Requires that FindSelfSimilarity.pl has been run before
***Requires that GetGenePairs_2shared.pl has been run before
****Requires that GetGenePairs_2shared.pl and GetGenePairs_1shared.pl has been run before

(*.tmp files are loaded automatically by the scripts. As these are not erased automatically, this may be done manually after a completed primer design procedure.) 


To get help, type:
   
   perl Scriptname.pl? or perl Scriptname.pl h
 


Single genome primer design
---------------------------
It is also possible to design primers for a single genome, with or without the -selfsim option. 

With -selfsim:

FindSelfSimilarity.pl		-fasta1 Genome1.fasta (-l number -s number)
DesignPrimers_uniques.pl	-fasta1 Genome1.fasta (-selfsim)*
*Requires that FindSelfSimilarity.pl has been run before

Without -selfsim:

DesignPrimers_uniques.pl	-fasta1 Genome1.fasta

Output will be: Genome1_uniques.tab
 
