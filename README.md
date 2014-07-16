VAST-TOOLS
==========

Summary
-------
VAST-TOOLS is a toolset for profiling alternative splicing events in RNA-Seq data. 

Requirements
------------

VAST-TOOLS requires the following software:
 * bowtie 1.0.0 (Langmead et al., 2009), http://bowtie-bio.sourceforge.net/index.shtml
 * R 3.0.1 or higher, with the following packages installed (see Installation Section):
  * optparse
  * RColorBrewer
  * reshape2
  * ggplot2
 * Perl 5.10.1 or higher
 * GNU coreutils `sort` version v8.1 or greater (`sort --version` to check)
 
VAST-TOOLS also requires species-specific library files (collectively known as
VASTDB), which must be downloaded separately. VASTDB consists of pre-made Bowtie
indices, annotation, and template files. See below for installation instructions.

Installation
------------
There are two components that need to be installed: VAST-TOOLS software package
and VASTDB library files.

### VAST-TOOLS

The VAST-TOOLS software package can be downloaded, unpacked and run as is. It
does not require any additional installation steps.

### VASTDB

VASTDB must be downloaded separately and can be saved in the VAST-TOOLS 
directory, or in an external location. If the latter, the path of VASTDB 
must be supplied to `vast-tools` via `--dbDir` or alternatively, a symbolic 
link can be created in the root of VAST-TOOLS directory. By default, 
VAST-TOOLS looks for VASTDB inside its own directory 
(e.g. `~/bin/vast-tools/VASTDB`).

~~~~
> git clone https://github.com/vastgroup/vast-tools.git
> cd vast-tools
~~~~
**Automatic DB Installation:**

Then you can automatically install the database files using:
~~~~
> ./install.R
~~~~
Follow the command prompt to install automatically, and that should be it!

Of course you should add the vast-tools directory to your PATH:
~~~~
$ export PATH=~/bin/vast-tools:$PATH
$ echo 'export PATH=~/bin/vast-tools:$PATH' >> ~/.bashrc
~~~~
**Manual DB Installation:**

For manual, install human (hsa), or mouse (mmu), or both to any location by:

Human (hg19) - 6.0G [vastdb.hsa.7.3.14.tar.gz](http://vastdb.crg.eu/libs/vastdb.hsa.7.3.14.tar.gz):
~~~~
> wget http://vastdb.crg.eu/libs/vastdb.hsa.7.3.14.tar.gz
> tar xzvf vastdb.hsa.7.3.14.tar.gz
~~~~
Mouse (mm9) - 7.9G [vastdb.mmu.7.3.14.tar.gz](http://vastdb.crg.eu/libs/vastdb.mmu.7.3.14.tar.gz):
~~~~
> wget http://vastdb.crg.eu/libs/vastdb.mmu.7.3.14.tar.gz
> tar xzvf vastdb.mmu.7.3.14.tar.gz
~~~~

If manually installed to central location, link the database files to vast-tools.0.0.1
directory using:
~~~~
> ln -s <path to VASTDB> VASTDB
~~~~

If you did a manual install, you can test and see if you have everything that is 
necessary to run all of vast-tools, OR vast-tools will try and install R packages 
on the fly when necessary.

~~~~
> ./install.R
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

VAST-TOOLS can be run as simply as:
NOTE: Unless specified, all options are default, for example the output
directory is assumed to be 'vast_out', the database to be
``<path>/vast-tools/VASTDB``, and the species ``Hsa``. To change these use
the ``--output``, ``-dbDir`` and ``-sp`` flags!

~~~~
> vast-tools align tissueA-rep1.fq.gz
> vast-tools align tissueA-rep2.fq.gz
> vast-tools align tissueB-rep1.fq.gz
> vast-tools align tissueB-rep2.fq.gz

> vast-tools combine

> vast-tools diff -a tissueA-rep1@tissueA-rep2 -b tissueB-rep1@tissueB-rep2 > INCLUSION-FILTERED.tab

> vast-tools plot INCLUSION-FILTERED.tab
~~~~

You can speed up vast-tools significantly by allowing it to use multiple cores
by running it on a cluster.  The ``-c`` flag can be passed to both ``align`` and
``diff``.

~~~~
> vast-tools align tissueA-rep1.fq.gz -c 8
~~~~ 
AND
~~~~
> vast-tools diff -a tissueA-rep1@tissueA-rep2 -b tissueB-rep1@tissueB-rep2 -c 8 > INCLUSION-FILTERED.tab
~~~~

### Alignment

In this step, reads are first aligned against a reference genome to obtain
unmapped reads, and these are then aligned to predefined splice junction libraries. Unmapped reads are saved
in the output directory as ``<sample>-<length>-e.fq``, where ``sample`` is the sample
name and ``length`` is the trimmed read length (e.g. 50). The input reads can be
compressed (via gzip) or uncompressed.

Currently, VAST-TOOLS supports two species, human (Hsa) and mouse (Mmu). By
default, the ``-sp`` option is ``Hsa``.

To enable gene expression analysis, use either the option ``-expr`` (PSI/PSU/PIRs plus
cRPKM calculations [corrected-for-mappability Reads per Kbp and Million mapped reads; see Labbé *et al*, 2012 for details]) or ``-exprONLY`` (cRPKMs only). For example, to perform
alignment with expression analysis on mouse data:

~~~~
> vast-tools align mouse_tissue.fq.gz -sp Mmu -expr
~~~~

If this alignment step needs to be repeated, the initial genome alignment step
can be skipped by supplying the ``<sample>-<length>-e.fq`` file as input. VAST-TOOLS
will recognize the \"-e.fq\" suffix and start at the splice junction alignment
step. Gene expression analysis *cannot* be run from this stage (you must start
from the raw reads).

~~~~
> vast-tools align mouse_tissue-e.fq.gz -sp Mmu
~~~~

Although you can specify two fastq files to vast-tools in a 'paired-end' format,
the program treats both mates independently because of trimming such that
`> vast-tools align fwd.fq.gz rev.fq.gz` is the same as running the program with
`> vast-tools align <(cat fwd.fq.gz rev.fq.gz)` catenated reads.

### Combining Results 

``vast-tools combine`` will join all of the files sent to the same output
directory found in <output_dir>/to_combine/, to form one final table in the main
<output_dir> folder.  This is the file you give to ``diff`` in the case that you
intend to compare multiple samples.  This output file contains a value for the percent of sequence inclusion (PSI/PSU/PIR) and a qual column for each sample. Details on the output format are provided below.

~~~~
> vast-tools combine -o outputdir -sp [Hsa|Mmu]
~~~~

### Differential Splicing Analysis

Bayesian inference followed by differential analysis of joint emperical posterior distributions with respect to
PSI/PSU/PIR.  

Diff Specific Inquiries: Tim Sterne-Weiler [email](mailto:tim.sterne.weiler@utoronto.ca) - [web](http://www.utoronto.ca/intron/sterne-weiler.html)


~~~~
> vast-tools diff -a sampleA_rep1@sampleA_rep2 -b sampleB_rep1@sampleB_rep2 -o outputdir > outputdir/significant_events.tab
~~~~

*Statistics Options*

Probably the most important extra options to consider are ``-r PROB (--prob)``,
``-m MINDIFF (--minDiff)`` and ``-e MINREADS (--minReads)`` These represent the stringency criterion for
visual output and filtering of input to STDOUT.

The ``-r`` flag represents the
minimal probability of acceptance that is required to consider a comparison to
be 'believable'.  By default this is 0.95, but it can be altered depending on
stringency requirements.  

The ``-m`` flag represents the minimum difference between psi1 and psi2 that you
will accept, such that we are are sure with at least probability ``-r`` that
there is a difference of at least ``-m``.

The ``-e`` flag specifies the minimum number of reads for a sample/event to be
compared.  In cases where the prior distribution has been methodically calculated
and/or is believable beyond an uninformative prior (like the uniform default),
this may not be necessary, however it is still highly recommended.  The default 
value for ``-e`` is 5, though this could easily be higher.

Additionally, ``diff`` allows you to alter the parameters of the conjugate beta
prior distribution, which is set as a uniform beta with shape parameters
``--alpha`` and ``--beta`` as 1 and 1 respectively.
Beta shape parameters greater than one alter this probability distribution, and
may be more or less applicable to certain uses, see: [beta
distribution](http://mathworld.wolfram.com/BetaDistribution.html) 

In the case that you have paired samples, where NormalA is dependent on
PerturbationA, it is appropriate to use the ``--paired=TRUE`` flag.  When
calculating the joint emperical posterior distribution,
for example from NormalA and NormalB, to compare to PerturbationA and
PerturbationB, the probability that P( joint_psi1 - joint_psi2 > ``-m`` ) is not
resampled such that NormalA is only compared to PerturbationA, and then NormalB
is compared to PerturbationB.  

 
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

### Plotting

VAST-TOOLS comes with a plotting script written in R.
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
help description: ``vast-tools plot -h``. An example of sample input data and
configuration file template can be found under ``R/sample_data``:

~~~~
> vast-tools plot -c R/sample_data/sample_psi_data.config R/sample_data/sample_psi_data.tab
~~~~

Combine output format
---------------------
The output of ``combine`` is a tab-separated table with an entry (row) for each predefined alternative splicing event. For each event, there are six columns with basic information about it, and then a pair of columns for each sample from ``align`` that is combined. 

 * **Column 1**: Official gene symbol.
 * **Column 2**: VAST-DB event ID. Formed by: 
  * Species identifier: Hsa (Human) or Mmu (Mouse);
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
  * S, C1, C2, C3: exon skipping (EX) events quantified by the *a posteriori* or *a priori* modules, with increasing degrees of complexity (based on *Score 5* for a wide panel of RNA-seq samples; see below and Irimia *et al.* 2014 for further information).
  * MIC: exon skipping (EX) events quantified by the microexon pipeline.
  * IR-S: intron retention event with no other annotated overlapping alternative splicing event and/or alternative first/last exons.
  * IR-C: intron retention event with other annotated overlapping alternative splicing event(s) and/or alternative first/last exons (similar to Type C introns in Braunschweig *et al*, 2014).  
  * Alt3: ALTA events.
  * Alt5: ALTD events.

Then, for each combined sample, a pair of columns: 
 * **Column 7**: Estimated percent of sequence inclusion (PSI/PSU/PIR). PSI: percent spliced in (for EX). PSU: percent splice site usage (for ALTD and ALTA). PIR: percent intron retention (for INT).
 * **Column 8**: Quality scores, and number of corrected inclusion and exclusion reads (qual@inc,exc).
  * *Score 1*: Read coverage, based on actual reads:
    - For EX: OK/LOW/VLOW: (i) ≥20/15/10 actual reads (i.e. before mappability correction) mapping to all exclusion splice junctions, OR (ii) ≥20/15/10 actual reads mapping to one of the two groups of inclusion splice junctions (upstream or downstream the alternative exon), and ≥15/10/5 to the other group of inclusion splice junctions.
    - For EX (microexon module): OK/LOW/VLOW: (i) ≥20/15/10 actual reads mapping to the sum of exclusion splice junctions, OR (ii) ≥20/15/10 actual reads mapping to the sum of inclusion splice junctions.
    - For INT: OK/LOW/VLOW: (i) ≥20/15/10 actual reads mapping to the sum of skipping splice junctions, OR (ii) ≥20/15/10 actual reads mapping to one of the two inclusion exon-intron junctions (the 5' or 3' of the intron), and ≥15/10/5 to the other inclusion splice junctions.
    - For ALTD and ALTA: OK/LOW/VLOW: (i) ≥40/20/10 actual reads mapping to the sum of all splice junctions involved in the specific event.
    - For any type of event: SOK: same thresholds as OK, but a total number of reads ≥100.
    - For any type of event: N: does not meet the minimum threshold (VLOW).

  * *Score 2*: Read coverage, based on corrected reads (similar values as per *Score 1*).
  * *Score 3*: Read coverage, based on uncorrected reads mapping only to the reference C1A, AC2 or C1C2 splice junctions (similar values as per *Score 1*).
  * *Score 4*: Imbalance of reads mapping to inclusion splice junctions (only for exon skipping events quantified by the *a posteriori* or *a priori* modules; for intron retention, see Braunschweig *et al*, 2014).
    - OK: the ratio between the total number of reads supporting inclusion for splice junctions upstream and downstream the alternative exon is < 2.
    - B1: the ratio between the total number of reads supporting inclusion for splice junctions upstream and downstream the alternative exon is > 2 but < 5.
    - B2: the ratio between the total number of reads supporting inclusion for splice junctions upstream and downstream the alternative exon is > 5.
    - Bl/Bn: low/no read coverage for splice junctions supporting inclusion.
  * *Score 5*: Complexity of the event (only for exon skipping events quantified by the *a posteriori* or *a priori* modules).
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

Citations
---------
Coming soon...

References
----------

Labbé, R.M., Irimia, M., Currie, K.W., Lin, A., Zhu, S.J., Brown, D.D., Ross, E.J., Voisin, V., Bader, G.D., Blencowe, B.J., Pearson, B.J., 2012. A comparative transcriptomic analysis reveals conserved features of stem cell pluripotency in planarians and mammals. Stem Cells. 30 (8):1734-45.

Langmead, B., Trapnell, C., Pop, M., Salzberg, S.L., 2009. Ultrafast and
memory-efficient alignment of short DNA sequences to the human genome. Genome
Biol. 10:R25.
