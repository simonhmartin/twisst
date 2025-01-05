# Twisst: Topology weightying by iterative sampling of sub-trees

#### NOTE: Twisst2 is ready for testing. It provides significant improvements. It has it's own repository [here](https://github.com/simonhmartin/twisst2).

Topology weighting is a means to quantify relationships between taxa that are not necessarily monophyletic. It's a simple, descriptive method, designed for exploring how relationship vary across the genome using population genomic data.

The relationship among a given set of taxa can be defined by a number of possible topologies. For example, for four taxa labelled A, B, C and D, there are three possible (unrooted) bifurcating topologies:


```
  A-\    /-C         A-\    /-B         A-\    /-B
     |--|               |--|               |--|
  B-/    \-D         C-/    \-D         D-/    \-C

```
Given a tree with any number of tips (or leaves), each belonging to a particular taxon, the weighting of each taxon topology is defined as the fraction of all unique sub-trees, in which each taxon is represented by a single tip, that match that topology. Topology weighting therefore reduces the complexity of the full tree to a number of values, each giving the proportionate contribution of a particular taxon tree to the full tree. 

This code implements the method *Twisst* (topology weighting by iterative sampling of sub-trees), which does what it says: it computes the weightings by iteratively sampling sub-trees from the full tree and checking their topology. This can be slow if there are many tips (e.g. 4 taxa with ten tips each gives 10 000 unique subtrees to consider. But there are some shortcuts to speed things up - see [Weighting Method](#weighting-method) below.

---
### Papers
[Martin and Van Belleghem 2017](http://doi.org/10.1534/genetics.116.194720) is where we present the *Twisst* method in full and test it on simulated and real data. Please cite this when using the software.

[Van Belleghem et al. 2017](http://doi.org/10.1038/s41559-016-0052) is where we first applied *Twisst* to identify narrow shared blocks in regulatory reagions.

---

### Code

The main script, `twisst.py` implements the topology weighting.

It requires [`ete3`](http://etetoolkit.org/download/) and `numpy` (tested on version 1.8).

A typical command looks like this:

```bash
python twisst.py -t input.trees.gz -w output.weights.csv.gz -g A 1,2,3,4,5 -g B 6,7,8,9,10 -g C 11,12,13,14,15 -g D 16,17,18,19,20
```

You can get a full list ot command options with `python twisst.py -h`.

**NOTE** Due to ongoing improvements, the script `run_twisst_parallel.py` is currently not working. However, the speed improvements to the core script make it usable even for quite large data sets.

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

### Weighting Method

There are three options for the weighting method, specified with the `--method` flag.

The default method is `complete`. This will calculate the exact weightings by considering all subtrees. It performs a smart shortcut by collapsing monophyletic nodes, and can be quite fast if the taxa are well resolved and/or the tree is small. However, if the taxa are unresolved and the tree is large, this method may be too slow to be useful.
If there are too many combinations to test (default > 1,000,000), the script will abort and revert to the `fixed` method described below. You can control this behaviour using the `--abortCutoff` option.

`--method fixed` will cause the script to estimates the weightings by randomly sampling a fixed number of subtrees. Sampling is done with replacement, so that the errors will fit a simple binomial distribution. The default number of iterations is 10,000, which gives very narrow error margins, but this can be modified using the `--iterations` flag.

The original version of `twisst` had an option to keep sampling until some level of confidence was reached. However, with speed improvements, sampling of thousands of combinations is achievable and the error margins become tiny and essentially irellevant. There is a supplementary figure in [Martin and Van Belleghem 2017](http://doi.org/10.1534/genetics.116.194720) demonstrating this.

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

* For diploids, you can infer phase using [Beagle 4](https://faculty.washington.edu/browning/beagle/beagle.html) or [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html). Here is an example Beagle command:

```bash
java -Xmx12g -jar beagle.jar gt=input.vcf.gz out=output.vcf.gz impute=true nthreads=20 window=10000 overlap=1000 gprobs=false
```

* My script to generate the trees takes a simple genotype format as input, which gives the scaffold, position and genotype for each sample. My code to generate this file from a vcf is in my repo [genomics_general](https://github.com/simonhmartin/genomics_general). Here is an example command:

```bash
python parseVCF.py -i input.vcf.gz --skipIndels --minQual 30 --gtf flag=DP min=5 | gzip > output.geno.gz
```

* To get neighbour joining trees for snp windows, I have a script that runs [Phyml](http://www.atgc-montpellier.fr/phyml/) for windows, using parallelisation, and outputs a single trees file. You can get the script from my [genomics_general](https://github.com/simonhmartin/genomics_general) repo. Here is an example command:
```bash
python phyml_sliding_windows.py -T 10 -g input.phased.geno.gz --prefix output.phyml_bionj.w50 -w 50 --windType sites --model GTR --optimise n
```
**NOTE: if you use phased diploid genotypes in `phyml_sliding_windows.py`, the output trees will include two tips for each sample, with suffixes "_A" and "_B". You will need to ensure that the groupd defined for `Twisst` match these names.**













