# Twisst: Topology weightying by iterative sampling of sub-trees

Topology weighting is a means to quantify relationships between taxa that are not necessarily monophyletic. It's a simple, descriptive method, designed for exploring how relationship vary across the genome using population genomic data.

The relationship among a given set of taxa can be defined by a number of possible topologies. For example, for four taxa labelled A, B, C and D, there are three possible (unrooted) bifurcating topologies:


```
  A-\    /-C
     |--|   
  B-/    \-D


  A-\    /-B
     |--|
  C-/    \-D


  A-\    /-B
     |--|
  D-/    \-C
```
Given a tree with any number of tips (or leaves), each belonging to a particular taxon, the weighting of each taxon topology is defined as the fraction of all unique sub-trees, in which each taxon is represented by a single tip, that match that topology. Topology weighting therefore reduces the complexity of the full tree to a number of values, each giving the proportionate contribution of a particular taxon tree to the full tree. 

This code implements the method *Twisst* (topology weighting by iterative sampling of sub-trees), which does what it says: it computes the weightings by iteratively sampling sub-trees from the full tree and checking their topology. This can be slow if there are many tips (e.g. 4 taxa with ten tips each gives 10 000 unique subtrees to consider. But there are some shortcuts to speed things up.

Firstly, monophyletic groups of samples from the same taxon can be collapsed and weighted appropriately. This dramatically reduces the number of unique subtrees to consider. However, in some cases there are still too many to compute the weighting in reasonable time. In these cases it is best to randomly sample a limited number of subtrees and take the estimates weightings. The errors are binomially distributed, so it's possible to also estimate confidence intervals.

---

### Code

The main script, `twisst.py` implements the topology weighting.

It requires Pythion 2.7 and the libraries [`ete3`](http://etetoolkit.org/download/) and `numpy` (tested on version 1.8).

A typical command looks like this:

```bash
python twisst.py -t input.trees.gz -w output.weights.csv.gz -o topologies.trees -g A 1,2,3,4,5 -g B 6,7,8,9,10 -g C 11,12,13,14,15 -g D 16,17,18,19,20 --method complete
```

You can get a full list ot command options with `python twisst.py -h`.

The script `run_twisst_parallel.py` allows parallelisation using python `multiprocessing`. It requires `twisst.py` to be present in the same directory or in your python path. The command line is the same, except that a number of threads must be specifies with the `-T` flag. The parallel version offers considerable speedups if the trees are large and complex, but little improvement (or even a slowdown) for small, simple trees that can be analysed very rapidly. 

---

### Input

The main input is a tree file containing one or more trees in newick format. Most variants of newick are accepted - see the [ETE documentation](http://etetoolkit.org/docs/latest/reference/index.html) for details.

Multiple trees should be listed on separate lines.

All trees must contain all specified individual names as tip labels (see below.)

---

### Output

There are two outputs:

1. The topologies file is specified with the `-o` flag. This file contains the possible taxon topologies in newick format.

2. The weights file is specified with the `-w` flag. This is a comma separated file, with one column for each topology, giving its weighting. The weightings are given in absolute counts (rather than proportions). This allows for estimation of confidence intervals downstream, if desired.

---

### Specifying taxa and individuals 

Taxa (groups) must be specified in the command line, using the `-g` flag. This flag must be present at least four times (with three groups there is only one possible unrooted topology).

The name of the group must be given after the flag, (optionally) followed by the names of the individuals (tip labels) that it contains, separated by commas (e.g: `-g A 1,2,3,4,5`).

Alternatively (or additionally), a tab-delimited file of tip labels and their corresponding groups can be provided with the `--groupsFile` flag. This should have tip labels in the first column and group names in the second column. The group names given in the groups file must match those in the command line.

For example, the example command above could alternatively be specified as 

```bash
python twisst.py -t input.trees.gz -w output.weights.csv.gz -o topologies.trees -g A g B -g C -g D --method complete --groupsFile groups.tsv
```

Where groups.tsv is a text file containing the following:

```
1	A
2	A
3	A
4	A
5	A
6	B
7	B
8	B
9	B
10	B
11	C
12	C
13	C
14	C
15	C
16	D
17	D
18	D
19	D
20	D
```
---

### Weighting method

There are three options for the weighting method, specified with the `--method` flag.

The recommended method is `complete`. This will calculate the exact weightings by considering all subtrees. It performs a smart shortcut by collapsing monophyletic nodes, and can be quite fast if the taxa are well resolved and/or the tree is small. However, if the taxa are unresolved and the tree is large, this method may be too slow to be useful.
If there are too many combinations to test (default > 100000), the script will abort and try one of the two approximate methods described below. You can control this behaviour using the `--abortCutoff` and `--backupMethod` flags. Note, `complete` is currently not the default method, and must be specified explicitly.

The default method is `fixed`. This estimates the weightings by randomly sampling a fixed number of subtrees. Sampling is done with replacement, so that the errors will fit a simple binomial distribution. The default number of iterations is 400, which gives fairly narrow error margins, but this can be modified using the `--iterations` flag. One issue with this method is that the confidence in the weightings will be much higher for some trees than others, because of the nature of the binomial distribution. For example, if there are three possible topologies that occur at frequencies 400, 0 and 0, respectively, we can be highly confident of their weightings. But if they occur at frequencies of 120, 130 and 150, we have less certainty of their exact weightings. In other words, to get a certain level of confidence, some trees need less sampling than others. For this reason, the threshold sampling option is recommended over this one.

The third option is `threshold`. This is similar to the sampling method above, except that the sampling is repeated until a certain dynamic threshold is reached. After *n* iterations, each topology must have been observed **fewer** than *k* times OR **more** than *n-k* times, with combinations of *n* and *k* being specified by the user. This allows specification of a unique threshold dependng on the number of iterations (*n*) and the number of observations (*k*). For example, you can set the thresholds such that the 95% binomial confidence interval around each weighting is less than 5%. Thresholds must be specified in a file, using the flag `--thresholdTable`. The file gives iterations (*n*) in the first column and the threshold number of observations (*k*) in the second column. If tyhe particular *n* has no acceptible *k* (often the case for low *n*), then *k* of -1 can be specified. Two example threshold tables are provided. These give thresholds to ensure that the 95% CI (Calculated using the "Wilson" method) is less than 5% or 10%.

### Pipeline to generate the input trees file

There are various options for producing trees for windows across the genome. If you have whole genome sequence data, it is recommended to infer trees for narrow genomic intervals. 50 SNPs proved a useful window size in various simulations. If the window is too large, you may be averaging over regions of distinct ancestry, which can eliminate subtle quantitative variation in taxon relationships. However, if the interval is too small, you may have insufficient signal to infer a good tree.

My approach is to subset the alignment into windows and infer trees for each separately, using either maximum likelihgood or neighbour joining methods. [PhyML](http://www.atgc-montpellier.fr/phyml/) can do both. [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/) has an option to do the sliding windows automatically. I don't recommend [Saguaro](http://saguarogw.sourceforge.net/) as it tends to be biased towards the most abundant topologies.

#### My pipeline from BAM to trees

* Starting with bam files, I genotype with [GATK](https://software.broadinstitute.org/gatk/), using the `HaplotypeCaller` and `GenotypeGVCFs` tools. Here are example commands:

```bash
#HaplotypeCaller
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -nct 16 -R reference.fa -I input.bam -o output.g.vcf --emitRefConfidence GVCF --output_mode EMIT_ALL_CONFIDENT_SITES
#GenotypeGVCFs
java -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -nt 16 -R reference.fa -V output.g.vcf --includeNonVariantSites -o output.vcf
```

* I filter vcfs using [Bcftools](https://samtools.github.io/bcftools/). This allows you to remove indels and invariant sites. I also like to convert uncertain genotypes to missing (`./.`) before phasing. Here is an example command that will convert all genotypes with < 5x depth of coverage and GQ < 30 to missing.

```bash
bcftools filter -e 'FORMAT/DP < 5 | FORMAT/GQ < 30' --set-GTs . input.vcf.gz -O u | bcftools view -U -i 'TYPE=="snp" & MAC >= 2' -O z > output.vcf.gz
``` 

* For diploids, I infer phase using [Beagle 4](https://faculty.washington.edu/browning/beagle/beagle.html), although I plan to start using [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html). Here is an example Beagle command:

```bash
java -Xmx12g -jar beagle.jar gt=input.vcf.gz out=output.vcf.gz impute=true nthreads=20 window=10000 overlap=1000 gprobs=false
```

* My script to generate the trees takes a simple genotype format as input, which gives the scaffold, position and genotype for each sample. My code to generate this file from a vcf is in my repo [genomics_general](https://github.com/simonhmartin/genomics_general). Here is an example command:

```bash
python parseVCF.py -i input.vcf.gz --skipIndel --minQual 30 --gtf flag=DP min=5 | gzip > output.geno.gz
```

* To get neighbour joining trees for snp windows, I have a script that runs [Phyml](http://www.atgc-montpellier.fr/phyml/) for windows, using parallelisation, and outputs a single trees file. You can get the script from my [genomics_general](https://github.com/simonhmartin/genomics_general) repo. Here is an example command:
```bash
python phyml_sliding_windows.py -T 10 -g input.phased.geno.gz --prefix output.phyml_bionj.w50 -w 50 --windType sites --model GTR --genoFormat phased
```
**NOTE: if you use phased diploid genotypes in `phyml_sliding_windows.py`, the output trees will include two tips for each sample, with suffixes "_A" and "_B". You will need to ensure that the groupd defined for `Twisst` match these names.**













