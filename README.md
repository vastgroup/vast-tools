VAST-TOOLS
==========

Summary
-------
VAST-TOOLS is a toolset for profiling alternative splicing events in RNA-Seq data. 

Requirements
------------

VAST-TOOLS requires the following software:
 * bowtie 1.0.0 (Langmead et al., 2009), http://bowtie-bio.sourceforge.net/index.shtml
 * R 3.0.1 or higher, with the following packages installed:
  * optparse
  * RColorBrewer
  * reshape2
  * ggplot2
 * Perl 5.10.1 or higher

VAST-TOOLS also requires species-specific library files (collectively known as
VASTDB), which must be downloaded separately. VASTDB consists of pre-made Bowtie
indices, annotation, and template files.

Installation
------------
There are two components that need to be installed: VAST-TOOLS software package
and VASTDB library files.

### VAST-TOOLS
The VAST-TOOLS software package can be downloaded, unpacked and run as is. It
does not require any additional installation steps.

### VASTDB
(VASTDB download instructions coming soon)

VASTDB must be downloaded separately and can be saved in the VAST-TOOLS 
directory, or in an external location. If the latter, the path of VASTDB must be
supplied to `vast-tools` via ``--dbDir`` or alternatively, a symbolic link can be created in the
root of VAST-TOOLS directory. By default, VAST-TOOLS looks for VASTDB
inside its own directory (e.g. `~/bin/vast-tools-0.0.1/VASTDB`).

~~~~
> cd ~/bin/vast-tools-0.0.1/
> ln -s <path to VASTDB> VASTDB
~~~~

You can test and see if you have everything installed that is necessary to run all 
of vast-tools, OR vast-tools will try and install R packages on the fly when necessary.

~~~~
> ./install.packages.R
~~~~

And of course you should add the vast-tools-0.0.1 directory to your PATH:
IN bash
~~~~
$ export PATH=~/bin/vast-tools-0.0.1:$PATH
$ echo 'export PATH=~/bin/vast-tools-0.0.1:$PATH' >> ~/.bashrc
~~~~

OR IN csh
~~~~
% setenv PATH ~/bin/vast-tools-0.0.1/:$PATH
% echo 'setenv PATH ~/bin/vast-tools-0.0.1/:$PATH' >> ~/.cshrc
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
NOTE: Unless specified, all options are default, for example the output directory 
is assumed to be 'vast_out', the database to be <path>/vast-tools-0.0.1/VASTDB, 
and the species Hsa.. to change these use the ``--output``, ``-dbDir`` and ``-sp`` flags!

~~~~
> vast-tools align tissueA-rep1.fq.gz
> vast-tools align tissueA-rep2.fq.gz
> vast-tools align tissueB-rep1.fq.gz
> vast-tools align tissueB-rep2.fq.gz

> vast-tools combine

> vast-tools diff -a tissueA-rep1@tissueA-rep2 -b tissueB-rep1@tissueB-rep2 > INCLUSION-FILTERED.tab

> vast-tools plot INCLUSION-FILTERED.tab
~~~~

You can speed up vast-tools significantly by allowing it to use multiple cores by running it on a cluster.
The ``-c`` flag can be passed to both 'align' and 'diff'.

~~~~
> vast-tools align tissueA-rep1.fq.gz -c 8
~~~~
AND
~~~~
> vast-tools diff -a tissueA-rep1@tissueA-rep2 -b tissueB-rep1@tissueB-rep2 --filter=T -c 8 > INCLUSION-FILTERED.tab
~~~~

### Alignment

In this step, reads are first aligned against a reference genome to obtain
unmapped reads, followed by splice junction libraries. Unmapped reads are saved
in the output directory as <sample>-<length>-e.fq, where <sample> is the sample
name and <length> is the trimmed read length (e.g. 50). The input reads can be
compressed (via gzip) or uncompressed.

To enable gene expression analysis, use either the option -expr (PSIs plus cRPKM
calculations) or -exprONLY (cRPKMs only).

If this alignment step needs to be repeated, the initial genome alignment step
can be skipped by supplying the <sample>-<length>-e.fq file as input. VAST-TOOLS
will recognize the \"-e.fq\" suffix and start at the splice junction alignment
step.

### Combining Results 

``vast-tools combine`` will join all of the files sent to the same output
directory found in <output_dir>/to_combine/, to form one final table in the main
<output_dir> folder.  This is the file you give to ``diff`` in the case that you
intend to compare multiple samples.  This output file contains a psi value and a
qual column for each sample.

~~~~
> vast-tools combine -o outputdir
~~~~

### Differential Splicing Analysis

Bayesian inference to estimate the posterior distribution over Psi followed by
differential analysis of joint emperical posterior distributions with respect to
Psi.  

Author Inquiries: Tim Sterne-Weiler [tim dot sterne dot weiler at
utoronto.ca](mailto:tim.sterne.weiler@utoronto.ca)
[web](www.utoronto.ca/intron/sterne-weiler.html)

~~~~
> vast-tools diff -a sampleA_rep1@sampleA_rep2 -b sampleB_rep1@sampleB_rep2 -o outputdir > significant_events.tab
~~~~

*Statistics Options*

Probably the most important extra options to consider are ``-r PROB (--prob)``
and ``-m MINDIFF (--minDiff)`` These represent the stringency criterion for
visual output and filtering of input to STDOUT The ``-r`` flag represents the
minimal probability of acceptance that is required to consider a comparison to
be 'believable'.  By default this is 0.9, but it can be increased depending on
stringency requirements.  
The ``-m`` flag represents the minimum difference between psi1 and psi2 that you
will accept, such that we are are sure with at least probability ``-r`` that
there is a difference of at least ``-m``

Additionally, 'diff' allows you to alter the parameters of the conjugate beta
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

The 'diff' command is also able to run in parallel.., specify the number of
cores to use with ``-c INT``
Obviously more cores will increase the speed of 'diff', though it may increase
the RAM usage as well..

Using the ``-n`` flag to specify the number of lines to read/process at a time,
will set a max threshold to the RAM used by parallel processing with the ``-c``
flag.  A lower number means that 'diff' will use significantly less memory,
however by decreasing ``-n`` you have increased the number of times that the
``mclapply`` function must calculate the parallel processing overhead.  The
default is 100, which works well.

### Plotting

VAST-TOOLS comes with a plotting script written in R.
The input format follows the same format from the "combine" step described
above. The output is a pdf of scatterplots (one per AS event) of PSI values.
To execute from VAST-TOOLS, use the subcommand "plot":

~~~~
> vast-tools plot significant_events.tab
~~~~

It is recommended to filter the input file to a subset of events of interest
before plotting, such as those obtained from "vast-tools diff". Otherwise, the
resulting pdf file will be very large.

For advanced usage (e.g. plot customizations), see the help description:
``vast-tools plot -h``.

Questions 
---------
Send questions to Manuel Irimia: [mirimia@gmail.com](mailto:mirimia@gmail.com).

Please send the smallest example necessary order to reproduce the error. 

Issues
------
Please report all bugs and issues using the issue tracker on GitHub (url
coming soon!)

Authors
-------
Manuel Irimia
Nuno Barbosa-Morais
Ulrich Braunschweig
Sandy Pan
Kevin Ha
Tim Sterne-Weiler

Citation
--------
Coming soon...

References
----------

Langmead, B., Trapnell, C., Pop, M., Salzberg, S.L., 2009. Ultrafast and
memory-efficient alignment of short DNA sequences to the human genome. Genome
Biol. 10, R25.
