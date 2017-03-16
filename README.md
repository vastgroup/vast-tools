VAST-TOOLS
==========

Table of Contents:

- [Summary](#summary)
- [Requirements](#requirements)
- [Installation](#installation)
	- [VAST-TOOLS](#vast-tools-1)
	- [VASTDB](#vastdb)
- [Usage](#usage)
	- [Help](#help)
	- [Quick Usage](#quick-usage)
	- [Alignment](#alignment)
	- [Merging Outputs](#merging-outputs)
	- [Combining Results](#combining-results)
	- [Comparing PSIs Between Samples](#comparing-psis-between-samples)
	- [Differential Splicing Analysis](#differential-splicing-analysis)
	- [Plotting](#plotting)
- [Combine output format](#combine-output-format)
- [Issues](#issues)
- [Contributions](#contributions)
- [Citation](#citation)
- [References](#references)
	
Summary
-------
Vertebrate Alternative Splicing and Transcription Tools (VAST-TOOLS) is a toolset for profiling and comparing alternative splicing events in RNA-Seq data. 

Requirements
------------

VAST-TOOLS requires the following software:
 * bowtie 1.0.0 (Langmead et al., 2009), http://bowtie-bio.sourceforge.net/index.shtml
 * R 3.1 or higher, with the following packages installed (see Installation Section):
   * optparse
   * RColorBrewer
   * reshape2
   * ggplot2 >= v2.0
   * MASS
   * devtools
   * [psiplot](https://github.com/kcha/psiplot)
 * Perl 5.10.1 or higher
 * GNU coreutils `sort` (all versions)
 
VAST-TOOLS also requires species-specific library files (collectively known as
VASTDB), which must be downloaded separately. VASTDB consists of pre-made Bowtie
indices, annotation, and template files. See below for installation instructions.

Installation
------------
There are two components that need to be installed: VAST-TOOLS software package
and VASTDB library files.

### VAST-TOOLS

The VAST-TOOLS software package can be downloaded, unpacked and run as is. It
does not require any additional installation steps. See [Releases](https://github.com/vastgroup/vast-tools/releases)
for the latest release or to get the latest development version, simply clone this repo:

~~~~
> git clone https://github.com/vastgroup/vast-tools.git
~~~~

### VASTDB

VASTDB must be downloaded separately and can be saved in the VAST-TOOLS 
directory, or in an external location. If the latter, the path of VASTDB 
must be supplied to `vast-tools` via `--dbDir` or alternatively, a symbolic 
link can be created in the root of VAST-TOOLS directory. By default, 
VAST-TOOLS looks for VASTDB inside its own directory 
(e.g. `~/bin/vast-tools/VASTDB`).

**Automatic DB Installation:**

You can automatically install the database files (and any pre-requisite R
packages) using:
~~~~
> ./install.R
~~~~
Follow the command prompt to install automatically, and that should be it!

Of course you should add the vast-tools directory (e.g. ~/bin/vast-tools) to your PATH:
~~~~
$ export PATH=~/bin/vast-tools:$PATH
$ echo 'export PATH=~/bin/vast-tools:$PATH' >> ~/.bashrc
~~~~
**Manual DB Installation:**

For manual, install human (hsa), mouse (mmu), chicken (gga), or all of them to any location by:

Human (hg19) - 6.2G [vastdb.hsa.22.06.16.tar.gz](http://vastdb.crg.eu/libs/vastdb.hsa.22.06.16.tar.gz):
~~~~
> wget http://vastdb.crg.eu/libs/vastdb.hsa.22.06.16.tar.gz
> tar xzvf vastdb.hsa.22.06.16.tar.gz
~~~~
Mouse (mm9) - 5.6G [vastdb.mmu.22.06.16.tar.gz](http://vastdb.crg.eu/libs/vastdb.mmu.22.06.16.tar.gz):
~~~~
> wget http://vastdb.crg.eu/libs/vastdb.mmu.22.06.16.tar.gz
> tar xzvf vastdb.mmu.22.06.16.tar.gz
~~~~
Chicken (galGal3) - 1.4G [vastdb.gga.13.11.15.tar.gz](http://vastdb.crg.eu/libs/vastdb.gga.13.11.15.tar.gz):
~~~~
> wget http://vastdb.crg.eu/libs/vastdb.gga.13.11.15.tar.gz
> tar xzvf vastdb.gga.13.11.15.tar.gz
~~~~
Planarian (v31) - 942M [vastdb.sme.31.11.15.tar.gz](http://vastdb.crg.eu/libs/vastdb.sme.31.11.15.tar.gz):
~~~~
> wget http://vastdb.crg.eu/libs/vastdb.sme.31.11.15.tar.gz
> tar xzvf vastdb.sme.31.11.15.tar.gz
~~~~

If manually installed to central location, link the database files to vast-tools directory using:
~~~~
> ln -s <path to VASTDB> VASTDB
~~~~

If you did a manual install, you can test and see if you have everything that is 
necessary to run all of vast-tools, OR vast-tools will try and install R packages 
on the fly when necessary.

~~~~
> ./install.R
~~~~

### R packages

If you did not use `install.R` as described above, you can manually install the
required R packages using the following commands:

~~~~
> R -e 'install.packages(c("optparse", "RColorBrewer", "reshape2", "ggplot2", "devtools"))'
> R -e 'devtools::install_github("kcha/psiplot")'
~~~~


Usage
-----

VAST-TOOLS contains four main sub-commands: align, combine, diff, and plot. The
following provides a quick introduction on how to run the tool:

### Help

Command usage can be retrieved through the -h (--help) flag to any sub-command:

~~~~
> vast-tools [sub-command] -h
~~~~

### Quick Usage

NOTE: Unless specified, all options are default, for example the output
directory is assumed to be 'vast_out', the database to be
``<path>/vast-tools/VASTDB``, and the species ``Hsa``. To change these use
the ``--output``, ``-dbDir`` and ``-sp`` flags!

VAST-TOOLS can be run as simply as:

~~~~
> vast-tools align tissueA_rep1.fq.gz
> vast-tools align tissueA_rep2.fq.gz
> vast-tools align tissueB_rep1.fq.gz
> vast-tools align tissueB_rep2.fq.gz

> vast-tools combine

> vast-tools compare -a tissueA_rep1,tissueA_rep2 -b tissueB_rep1,tissueB_rep2
OR
> vast-tools diff -a tissueA_rep1,tissueA_rep2 -b tissueB_rep1,tissueB_rep2 > INCLUSION-FILTERED.tab
> vast-tools plot INCLUSION-FILTERED.tab
~~~~

You can speed up vast-tools significantly by allowing it to use multiple cores
by running it on a cluster.  The ``-c`` flag can be passed to both ``align`` and
``diff``.

~~~~
> vast-tools align tissueA_rep1.fq.gz -c 8
~~~~ 
AND
~~~~
> vast-tools diff -a tissueA_rep1,tissueA_rep2 -b tissueB_rep1,tissueB_rep2 -c 8 > INCLUSION-FILTERED.tab
~~~~

### Alignment

In this step, to increase the fraction of mapping junction reads within each RNA-Seq sample, each read is first split into 50-nucleotide (nt) read groups, using by default a sliding window of 25 nt (``--stepSize`` option). For example, a 100-nt read would produce 3 overlapping reads (from positons 1-50, 26-75, and 51-100). In addition, both read mates from the paired-end sequencing are pooled, if available. For quantification, only one random count per read group (i.e. all sub-reads coming from the same original read) is considered to avoid multiple counting of the same original sequenced molecule. VAST-TOOLS ``align`` can also be used with pre-trimmed reads (``--pretrimmed`` option), but *only* if reads have been trimmed by VAST-TOOLS. (Read headings need a special format so that they are properly recognized by VAST-TOOLS subscripts, and these are generated during the trimming process). Also, it is highly recommended that special characters ('-', '.', etc.) are not part of the fastq file names as this may cause unforeseen problems; use '_' instead. (The use of '-' is reserved for providing the read length (legacy) or specify the reads have been genome substracted; see below).

Next, these 50-nt split reads are aligned against a reference genome to obtain
unmapped reads, and these are then aligned to predefined splice junction libraries. Unmapped reads are saved
in the output directory as ``<sample>-50-e.fa.gz``, where ``sample`` is the sample
name. The input reads can be compressed (via gzip) or uncompressed.

Currently, VAST-TOOLS supports three species, human (Hsa), mouse (Mmu), chicken (Gga), and planarian (Sme). By
default, the ``-sp`` option is ``Hsa``.

To enable gene expression analysis, use either the option ``--expr`` (PSI/PSU/PIRs pluscRPKM calculations [corrected-for-mappability Reads per Kbp and Million mapped reads; see Labbé *et al*, 2012 for details]) or ``--exprONLY`` (cRPKMs only). cRPKMs are obtained by mapping only the first 50 nucleotides of each read, or only the first 50 nucleotides of the forward read if paired-end reads are provided.

In addition to a file named *.cRPKM containing cRPKMs for each gene, a file name *.3bias will be created. This file contains information to estimate 3′ sequencing biases in the RNA-seq sample. Each of *.3bias file contains two rows:
- For all mRNAs with >200 mapped reads throughout at least 2500 nt, it provides the percentage of reads within the last five 500-nt windows, starting from the 3′-most window. The last column corresponds to the number of probed genes.
- For all mRNAs with >500 mapped reads throughout at least 5000 nt, it provides the percentage of reads within the last five 1000-nt windows, starting from the 3′-most window. The last column corresponds to the number of probed genes.

For example, to perform alignment with expression and 3′bias analysis on mouse data:

~~~~
> vast-tools align mouse_tissue.fq.gz -sp Mmu --expr
~~~~

If this alignment step needs to be repeated, the initial genome alignment step
can be skipped by supplying the ``<sample>-50-e.fa.gz`` file as input. VAST-TOOLS
will recognize the \"-e.fa\" suffix and start at the splice junction alignment
step. Gene expression and intron retention analyses *cannot* be run from this stage (you must start
from the raw reads).

~~~~
> vast-tools align mouse_tissue-50-e.fa.gz -sp Mmu
~~~~

Although you can specify two fastq files to vast-tools in a 'paired-end' format,
the program treats both mates independently because of trimming, but will not double
count the any trim or mate pair more than once  (see above). Reads must be given to the program
such that `vast-tools align fwd_mate_1.fq.gz rev_mate_2.fq.gz` refers to two fastq
files of identical line number where Read1 from file_1 is mated to Read1 from file_2. NOTE: if reads are downloaded from SRA as sra files, use ``fastq-dump --split-file ./sample.sra`` to generate separate fastq files for each paired-end (plus a third file with unmatched mates).

Finally, from `vast-tools v1.0.0` it is possible to calculate intron retention using two slightly different approaches with the ``--IR_version 1/2`` option. `Version 1` is the one used in Braunschweig *et al* 2014. `Version 2` incorporates alternative exon-exon junction for skipping reads to obtain a more representative Percent Intron Retention at the transcript level.

### Merging Outputs

``vast-tools merge`` can be used to pull the ``align`` outputs from various samples together into a new set of output files. This can be used to merge replicates when the read coverage for the independent replicates is not deep enough for a proper AS analysis. Unlike gene expression analysis, where 20-30 million reads is usually enough, differential AS analyses require deep sequencing coverage. Otherwise, only PSIs for AS events in highly expressed genes will be confidentely estimated, biasing the results. Since VAST-TOOLS uses only junction reads for quantification, we recommend at least 70 million reads per sample, and ideally >150 million reads. In our experience, the fraction of AS events whose PSI can be confidently estimated with VAST-TOOLS scales linearly with both read depth and length (specially if paired end) up to 200-300 million reads (depending on library complexity, read length, etc.). ``vast-tools merge`` is OPTIONAL; it is not needed in order to run ``vast-tools combine``.

``vast-tools merge`` needs a configuration file  that is provided through the ``--groups`` option. This file must have two TAB-separated columns like this:

~~~~
Sub_sample_A1	Sample_A
Sub_sample_A2	Sample_A
Sub_sample_A2	Sample_A
Sub_sample_B1	Sample_B
Sub_sample_B2	Sample_B
...
~~~~

The subsample names must be identical to the first part of the file names of the RNAseq data before the first occurrence of a dot in the file name; usually everything until the file ending .fq.gz .
For paired-end RNAseq data, the subsample name is defined by the first RNAseq file used during the align step.
``vast-tools merge`` will then merge the outputs of subsamples A1, A2 and A3 into a new Sample_A set of outputs that can then be used for combine.

~~~~
> vast-tools merge --groups config_file
~~~~

To merge cRPKM files from gene expression analysis, the option ``--expr`` needs to be provided. If so, ``vast-tools merge`` needs the species key (``--sp``) to upload the file with effective mappable positions per gene to recalculate cRPKMs. It is also possible to merge only cRPKM files by activating the option flag ``--exprONLY``.

Finally, the subsample files can be moved to a subfolder (`output_folder/PARTS`) by using the option ``--move_to_PARTS``. 


### Combining Results 

``vast-tools combine`` will join all of the files sent to the same output
directory found in <output_dir>/to_combine/, to form one final table in the main
<output_dir> folder.  This is the file you give to ``compare`` or ``diff`` in the case that you
intend to compare multiple samples.  This output file contains a value for the percent of sequence inclusion (PSI/PSU/PIR) and a qual column for each sample. Details on the output format are provided below. At least two samples must be combined. In addition, the version of intron retention used in ``align`` needs can be specified using ``--IR_version``. 
From release v1.0.0-beta.3, it is possible to get the output in mm10 and hg38 coordinates. vast-tools will still work on mm9 and hg19, respectively, but all output coordinates are then lifted to the newer assemblies. This need to be provided using the ``-a`` option. (Output in mm10 and hg38 need an extra file in VASTDB. If you have downloaded a VASTDB version older than vastdb.*.22.06.16, you will need to download the following patch: [PATCH_mm10-hg38.tar.gz](http://vastdb.crg.eu/libs/PATCH_mm10-hg38.tar.gz)

~~~~
> vast-tools combine -o outputdir -sp [Hsa|Mmu|Gga] --IR_version [1|2]
~~~~

### Comparing PSIs Between Samples

``vast-tools compare`` identifies differentially spliced AS events between two groups (A and B) based mainly on the difference in their average inclusion levels (i.e. ΔPSI = average_PSI_B - average_PSI_A). 

``vast-tools compare`` takes any number of replicates for each group, provided either as sample names or column numbers (0-based) in the INCLUSION table. Then, it first filters out those AS events that do not have enough read coverage in *ALL* replicates (see [Combine Output Format](#combine-output-format) below for information on coverage thresholds). For intron retention, it is possible to filter out introns also based on a binomial test to assess appropriate balance of reads at the two exon-intron junctions by using ``--p_IR`` (see Braunschweig *et al* 2014 for details). 

For valid AS events, ``vast-tools compare`` then requires that the absolute value of ΔPSI is higher than a threshold provided as ``--min_dPSI``. In addition, it requires that the PSI distribution of the two groups do not overlap. This can be modified with the ``--min_range`` option, to provide higher or lower stringency. For example: 

~~~~
> vast-tools compare INCLUSION_TABLE.tab -a sampleA1,sampleA2 -b sampleB1,sampleB2 --min_dPSI 25 --min_range 5
~~~~
 
will identify those AS events with an absolute ΔPSI between A and B higher than 25 and a minimum difference between the ranges of 5 (i.e. the maximum PSI value of group A and the minimum value of group B, if ΔPSI is positive, and the minimum PSI value of group A and the maximum value of group B, if ΔPSI is negative). ``--min_range`` can also take negative values to allow the PSI distributions of groups A and B to overlap to a defined extent (e.g. if ``--min_range -100``, events will be called as differentially spliced as long as |ΔPSI| > N when ``--min_dPSI N``).

By default, the comparison is not paired. A paired analysis can be done using the option ``--paired``. If provided, each A replicate is compared with the corresponding one in group B. In this case, ``--min_dPSI`` is the minimum threshold for the average between each pair's ΔPSI, and ``--min_range`` is the minimum value *any* individual pair's ΔPSI can be (i.e. all individual ΔPSI's absolute values have to be larger than the value provided as ``--min_range`` and they all have to have the same sign [positive or negative]). 

Summary statistics are printed for a tally of AS events by type that show higher sequence inclusion in group A or B, the total number of events in ``vast-tools`` that have been profiled (i.e. that pass any applied filter), and the subset of the latter that are alternatively spliced in at least one of the compared samples (10<PSI<90).

Differentially spliced AS events are printed out to an output file. The name of this file is generated by default based on the parameters used, but it can also be provided using ``--outFile``. The file is created in the directory of the input file. The selected events are then plotted automatically using ``vast-tools plot`` (see below), using a config file generated based on the type and order of the replicates provided as input. By default, all samples in the INCLUSION file are plotted, but it is possible to plot only the compared samples providing the ``only_samples`` option. The plotting step can be skipped by activating the ``--no_plot`` flag.

Finally, ``vast-tools compare`` can also produce list of gene IDs for the selected events to run Gene Ontology (GO) analyses using ``--GO``. In particular, it generates four list: (i) differentially spliced cassette exons and microexons, (ii) introns with higher retention in B (IR_UP), (iii) introns with higher retention in A (IR_DOWN), and (iv) backgroup set, for all multiexonic genes that meet similar read coverage criteria. (The latter is crucial to avoid GO enrichment of highly expressed genes in the specific cell or tissue type of study). To generate the list of gene IDs, VAST-TOOLS needs to access VASTDB and thus needs the species key provided with ``-sp``. Alternatively, a custom list of gene IDs for each AS event (`event1\tgene_id1`) can be provided using ``--GO_file``, or gene symbols from the first colum of the INCLUSION table can be used instead activating the ``--use_names`` flag.

### Differential Splicing Analysis

*IMPORTANT NOTE*: `diff` is still an experimental part of this package and is currently under development and testing. Please use at your own knowledge and risk.

Bayesian inference followed by differential analysis of posterior distributions with respect to
PSI/PSU/PIR.  With replicate data, joint posterior distributions for a sample are estimated from 
emperical posterior distributions of the replicates using maximum-likelihood (MLE) fitting.

~~~~
> vast-tools diff -a sampA_r1,sampA_r2,sampA_r3 -b sampB_r1,sampB_r2 -o outputdir > outputdir/diff_output.tab
~~~~
Note: Sample names do not have to follow any specific convention as long as they remain valid ASCII words and any number of replicate sample names may be given 1:N.


*Statistics Options*

Probably the most important extra options to consider are ``-r PROB (--prob)``,
``-m MINDIFF (--minDiff)``, ``-e MINREADS (--minReads)``, and `-S MINSAMPLES (--minSamples)` 
These represent the stringency criterion for filtering of visual output and textual 
data sent to file.  `-S` is the minimum number of samples for each set `-a` and `-b`
that have to have at least `-e` reads each to be considered in the downstream
statistical comparison.

The ``-r`` flag represents the
minimal probability of acceptance that is required to consider a comparison to
be 'believable'.  By default this is 0.95, but it can be altered depending on
stringency requirements.  

The ``-m`` flag represents the minimum value of difference (`MV`, see example below) 
between psi1 and psi2 that you will accept, such that we are are sure with at least 
probability ``-r`` that there is a difference of at least ``-m``.  `-m` does not 
currently alter the output sent to STDOUT, but does filter what is plotted to PDF
and printed to file.

The ``-e`` flag specifies the minimum number of reads for a sample/event to be
compared.  In cases where the prior distribution has been methodically calculated
and/or is believable beyond an uninformative prior (like the uniform default),
this may not be necessary, however it is still highly recommended.  The default 
value for ``-e`` is 10, though this could easily be higher.

Additionally, ``diff`` allows you to alter the parameters of the conjugate beta
prior distribution, which is set as a uniform beta with shape parameters
``--alpha`` and ``--beta`` as 1 and 1 respectively.
Beta shape parameters greater than one alter this probability distribution, and
may be more or less applicable to certain uses, see: [beta
distribution](http://www.wolframalpha.com/input/?i=beta+distribution) 
NOTE: When considering differential analysis of event types like intron retention
it may be more appropriate to use a custom prior model that is able to more accurately
reflect the lower expectation of inclusion levels.

In the case that you have paired samples, where NormalA is dependent on
PerturbationA, it is appropriate to use the ``--paired=TRUE`` flag.  For
example when considering NormalA and NormalB, to compare to PerturbationA and
PerturbationB, the probability that P( joint_psi1 - joint_psi2 > ``-m`` ) is
calculated such that NormalA is only compared to PerturbationA, and then NormalB
is compared to PerturbationB.  No MLE fitting is used in this case.

In all multireplicate cases where `--paired=FALSE`, the posterior distributions
of the individual replicates are used to estimate a 'best fit joint posterior' distribution
over psi for each sample.
 
*Performance Options*

The ``-s`` flag can be used to specify the ``-s SIZE`` of the emperical
posterior distribution to sample, lower numbers decrease accuracy but increase
performance.

The ``diff`` command is also able to run in parallel.., specify the number of
cores to use with ``-c INT``
Obviously more cores will increase the speed of ``diff``, though it may increase
the RAM usage as well..

Using the ``-n`` flag to specify the number of lines to read/process at a time,
will set a max threshold to the RAM used by parallel processing with the ``-c``
flag.  A lower number means that ``diff`` will use significantly less memory,
however by decreasing ``-n`` you have increased the number of times that the
``mclapply`` function must calculate the parallel processing overhead.  The
default is 100, which works well.

*Output Format*

The text output of diff looks like:

  GENE	| EVENT		| SampleA	| SampleB	|    E[dPsi]	| MV[abs(dPsi)]_at_0.95
:-------|:------------- |:------------- |:------------- |:------------- | ---------------------------
 BOD1L	| HsaEX0008312	| 0.124353	| 0.700205	| -0.575851	| 0.3		   	
 KARS	| HsaEX0032865	| 0.172134	| 0.460027	| -0.287892	| 0.22              
 NISCH	| HsaEX0043017	| 0.247743	| 0.500657	| -0.252915	| 0.09             
 ALAS1	| HsaEX0003568	| 0.293333	| 0.537553	| -0.244220	| 0.09		   
 VPS13D	| HsaEX0070518	| 0.984622	| 0.657589	| 0.327033	| 0.1		   
 TRO	| HsaEX0067335	| 0.757929	| 0.474551	| 0.283378	| 0.04		   
 USP33	| HsaEX0069762	| 0.337669	| 0.845228	| -0.507560	| 0.14		   
 BCORL1	| HsaEX0007940	| 0.213452	| 0.500425	| -0.286973	| 0.05		   

Where for example the first event HsaEX0008312 in the BOD1L gene has multireplicate point estimate
for SampleA of 0.12 and 0.7 for SampleB.  While this gives an expected value for the difference of
Psi (deltaPsi) between SampleA and SampleB of -0.57, the minimum value (`MV`) for |dPsi| at 0.95 is
0.3, meaning that there is a 0.95 probability that |deltaPsi| is greater than 0.3. Use this value 
to filter for events that are statistically likely to have at least a minimal difference of some 
magnitude that you deem to be biologically relevant.  The `-m` argument (default 0.1) provides a
lower bound for events that will be plotted to PDF and printed to file based on `MV`.  As a cutoff,
the default is meant to provide a reasonably stringent baseline, however you could relax this if you
would rather view more events that may have significant but modest changes.


![Diff](https://raw.githubusercontent.com/vastgroup/vast-tools/master/R/sample_data/DiffExample.png "Example") 

The output plot above shows in the left panel the two joint posterior distributions over psi, and the
point estimates for each replicate plotted as points below the histograms.  In the right panel:
the y-axis represents the probability of delta psi being greater than some magnitude value of x (shown on the x-axis).
The red line indicates the maximal value of x where P(deltaPsi > x) > `-r`, or the 0.95 default.  

### Plotting

VAST-TOOLS comes with a plotting script written in R. As of version 0.2.0, the
core functionality of this script has been ported to the R package
[psiplot](https://github.com/kcha/psiplot). Advanced users who would like to
generate plots interactively in R are encouraged to use this package directly.

The input format follows the same format from the ``combine`` step described
above. The output is a pdf of scatterplots (one per AS event) of PSI values.
To execute from VAST-TOOLS, use the subcommand ``plot``:

~~~~
> vast-tools plot outputdir/significant_events.tab
~~~~

It is recommended to filter the input file to a subset of events of interest
before plotting, such as those obtained from ``diff``. Otherwise, the
resulting pdf file will be very large. 

Plot customizations such as coloring and ordering of samples can be applied
using a configuration file. For more details on this advanced usage, see the
help description: ``vast-tools plot -h`` or see the
[README](https://github.com/kcha/psiplot#the-config-file-way) of the psiplot
package. An example of sample input data and
configuration file template can be found under ``R/sample_data``:

~~~~
> vast-tools plot -c R/sample_data/sample_psi_data.config R/sample_data/sample_psi_data.tab
~~~~

![Diff](https://raw.githubusercontent.com/vastgroup/vast-tools/master/R/sample_data/PsiExample.png "Example") 

Combine output format
---------------------
The output of ``combine`` is a tab-separated table with an entry (row) for each predefined alternative splicing event. For each event, there are six columns with basic information about it, and then a pair of columns for each sample from ``align`` that is combined. 

 * **Column 1**: Official gene symbol.
 * **Column 2**: VAST-DB event ID. Formed by: 
  * Species identifier: Hsa (Human), Mmu (Mouse), or Gga (Chicken);
  * Type of alternative splicing event: alternative exon skipping (EX), retained intron (INT), alternative splice site donor choice (ALTD), or alternative splice site acceptor choice (ALTA). In the case of ALTD/ALTA, each splice site within the event is indicated (from exonic internal to external) over the total number of alternative splice sites in the event (e.g. HsaALTA0000011-1/2).
  * Numerical identifier.
 * **Column 3**: Genomic coordinate of the alternative sequence.
 * **Column 4**: Length of the alternative sequence. In ALTD/ALTA events, the first splice site within each event has a length of 0 nt, by definition.
 * **Column 5**: Full set of genomic coordinates of the alternative splicing event. 
  * For EX: *chromosome:C1donor,Aexon,C2acceptor*. Where C1donor is the reference upstream exon's donor, C2acceptor the reference downstream exon's acceptor, and A the alternative exon. Strand is "+" if C1donor < C2acceptor. If multiple acceptor/donors exist in any of the exons, they are shown separated by "+". 
  * For ALTD: *chromosome:Aexon,C2acceptor*. Multiple donors of the event are separated by "+".
  * For ALTA: *chromosome:C1donor,Aexon*. Multiple acceptors of the event are separated by "+".
  * For INT: *chromosome:C1exon=C2exon:strand*.
 * **Column 6**: Type of event.
  * S, C1, C2, C3: exon skipping (EX) events quantified by the *splice site-based* or *transcript-based* modules, with increasing degrees of complexity (based on *Score 5* for a wide panel of RNA-seq samples; see below and Irimia *et al.* 2014 for further information).
  * MIC: exon skipping (EX) events quantified by the microexon pipeline.
  * IR-S: intron retention event with no other annotated overlapping alternative splicing event and/or alternative first/last exons.
  * IR-C: intron retention event with other annotated overlapping alternative splicing event(s) and/or alternative first/last exons (similar to Type C introns in Braunschweig *et al*, 2014).  
  * Alt3: ALTA events.
  * Alt5: ALTD events.

Then, for each combined sample, a pair of columns: 
 * **Column 7**: Estimated percent of sequence inclusion (PSI/PSU/PIR). PSI: percent spliced in (for EX). PSU: percent splice site usage (for ALTD and ALTA). PIR: percent intron retention (for INT).
 * **Column 8**: Quality scores, and number of corrected inclusion and exclusion reads (qual@inc,exc).
  * *Score 1*: Read coverage, based on actual reads (as used in Irimia *et al*, Cell 2014):
    - For EX: OK/LOW/VLOW: (i) ≥20/15/10 actual reads (i.e. before mappability correction) mapping to all exclusion splice junctions, OR (ii) ≥20/15/10 actual reads mapping to one of the two groups of inclusion splice junctions (upstream or downstream the alternative exon), and ≥15/10/5 to the other group of inclusion splice junctions.
    - For EX (microexon module): OK/LOW/VLOW: (i) ≥20/15/10 actual reads mapping to the sum of exclusion splice junctions, OR (ii) ≥20/15/10 actual reads mapping to the sum of inclusion splice junctions.
    - For INT: OK/LOW/VLOW: (i) ≥20/15/10 actual reads mapping to the sum of skipping splice junctions, OR (ii) ≥20/15/10 actual reads mapping to one of the two inclusion exon-intron junctions (the 5' or 3' of the intron), and ≥15/10/5 to the other inclusion splice junctions.
    - For ALTD and ALTA: OK/LOW/VLOW: (i) ≥40/20/10 actual reads mapping to the sum of all splice junctions involved in the specific event.
    - For any type of event: SOK: same thresholds as OK, but a total number of reads ≥100.
    - For any type of event: N: does not meet the minimum threshold (VLOW).

  * *Score 2*: Read coverage, based on corrected reads (similar values as per *Score 1*).
  * *Score 3*: Read coverage, based on uncorrected reads mapping only to the reference C1A, AC2 or C1C2 splice junctions (similar values as per *Score 1*). Always NA for intron retention events.
  * *Score 4*: Imbalance of reads mapping to inclusion splice junctions (only for exon skipping events quantified by the *splice site-based* or *transcript-based* modules; For intron retention events, numbers of reads mapping to the upstream exon-intron junction, downstream intron-exon junction, and exon-exon junction in the format A=B=C)
    - OK: the ratio between the total number of reads supporting inclusion for splice junctions upstream and downstream the alternative exon is < 2.
    - B1: the ratio between the total number of reads supporting inclusion for splice junctions upstream and downstream the alternative exon is > 2 but < 5.
    - B2: the ratio between the total number of reads supporting inclusion for splice junctions upstream and downstream the alternative exon is > 5.
    - Bl/Bn: low/no read coverage for splice junctions supporting inclusion.
  * *Score 5*: Complexity of the event (only for exon skipping events quantified by the *splice site-based* or *transcript-based* modules); For intron retention events, p-value of a binomial test of balance between reads mapping to the upstream and downstream exon-intron junctions, modified by reads mapping to a 200-bp window in the centre of the intron (see Braunschweig et al., 2014).
    - S: percent of complex reads (i.e. those inclusion- and exclusion-supporting reads that do not map to the reference C1A, AC2 or C1C2 splice junctions) is < 5%.
    - C1: percent of complex reads is > 5% but < 20%.
    - C2: percent of complex reads is > 20% but < 50%.
    - C3: percent of complex reads is > 50%.
    - NA: low coverage event.
  * inc,exc: total number of reads, corrected for mappability, supporting inclusion and exclusion.

Issues
------
Please report all bugs and issues using the GitHub [issue tracker]
(https://github.com/vastgroup/vast-tools/issues).

Contributions
-------------
* Manuel Irimia
* Nuno Barbosa-Morais 
* Ulrich Braunschweig 
* Sandy Pan
* Kevin Ha
* Tim Sterne-Weiler

Citation
--------

* VAST-TOOLS:

Irimia, M., Weatheritt, R.J., Ellis, J., Parikshak, N.N., Gonatopoulos-Pournatzis, T., Babor, M., Quesnel-Vallières, M., Tapial, J., Raj, B., O’Hanlon, D., Barrios-Rodiles, M., Sternberg, M.J.E., Cordes, S.P., Roth, F.P., Wrana, J.L., Geschwind, D.H., Blencowe, B.B. (2014). A highly conserved program of neuronal microexons is misregulated in autistic brains. *Cell*, 59:1511-23.

* Intron retention analysis:

Braunschweig, U., Barbosa-Morais, N.L., Pan, Q., Nachman, E., Alipahani, B., Gonatopoulos-Pournatzis, T., Frey, B., Irimia, M., Blencowe, B.J. (2014). Widespread intron retention in mammals functionally tunes transcriptomes. *Genome Research*, 24:1774-86

* Chicken database:

Gueroussov, S., Gonatopoulos-Pournatzis, T., Irimia, M., Raj, B., Lin, Z.Y., Gingras, A.C., Blencowe, B.J. (2015). An alternative splicing event amplifies evolutionary differences between vertebrates. *Science*, 349:868-73

References
----------

* cRPKMs:

Labbé, R.M., Irimia, M., Currie, K.W., Lin, A., Zhu, S.J., Brown, D.D., Ross, E.J., Voisin, V., Bader, G.D., Blencowe, B.J., Pearson, B.J. (2012). A comparative transcriptomic analysis reveals conserved features of stem cell pluripotency in planarians and mammals. *Stem Cells*, 30:1734-45.

* *bowtie*:

Langmead, B., Trapnell, C., Pop, M., Salzberg, S.L. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. *Genome Biology*, 10:R25.
