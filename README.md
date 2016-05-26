#Twisst: Topology weightying by iterative sampling of subtrees

Topology weighting is a means to quantify relationships between taxa that are not necessarily monophyletic. It's a simple, descriptive method, designed for exploring relationships in population genomic data.

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

###Code

The main script, `twisst.py` implements the topology weighting.

It required Pythion 2.7 and the libraries `ete3` and `numpy` (tested on version 1.8).

A typical command looks like this:

```bash
python twisst.py -t input.trees.gz -w output.weights.csv.gz -o topologies.trees -g A 1,2,3,4,5 -g B 6,7,8,9,10 -g C 11,12,13,14,15 -g D 16,17,18,19,20 --method complete
```

You can get a full list ot command options with `python twisst.py -h`.

The script `run_twisst_parallel.py` allows parallelisation using python `multiprocessing`. The command line is the same, except that a number of threads must be specifies with the `-T` flag. The parallel version offers considerable speedups if the trees are large and complex, but little improvement for small, simple trees that can be analysed very rapidly. 



###Input

The main input is a tree file containing one or more trees in newick format. Most variants of newick are accepted - see the `ete3` documentation for details.

Multiple trees should be listed on separate lines.

All trees must contain all specified tip labels (see below.)

###Output

There are two outputs:

1. The topologies file is specified with the `-o` flag. This file contains the possible taxon topologies in newick format.

2. The weights file is specified with the `-w` flag. This is a comma separated file, with one column for each topology, giving its weighting. The weightings are given in absolute counts (rather than proportions). This allows for estimation of confidence intervals downstream, if desired.

###Specifying taxa and individuals 

Taxa (groups) must be specified in the command line, using the `-g` flag. This flag must be present at least four times (otherwise there's only one possible unrooted topology).

The name of the group must be given after the flag, and optionally the names of the individuals (tip labels) that it contains, separated by commas (e.g: `-g A 1,2,3,4,5`).


Alternatively (or additionally), a tab-delimited file of tip labels and their corresponding groups can be provided with the `--groupsFile` flag. This should have tip labels in the first column and group names in the second column. The group names given in the groups file must match those in the command line.

For example, the command above could alternatively be specified as 

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

###Weighting method

There are three options for the weighting method, specified with the `--method` flag.

The recommended method is `complete`. This will calculate the exact weightings by considering all subtrees. It performs a smart shortcut by collapsing monophyletic nodes, and can be quite fast if the taxa are well resolved. However, if the taxa are unresolved and the tree is large, this method may be too slow to be useful. One of the other two approximate methods is recommended. For this reason, `complete` is not the default method, and must be specified explicitly.

The default method is `fixed`. This estimates the weightings by randomly sampling a fixed number of subtrees. Sampling is with replacement, so that the errors will fit a simple binomial distribution. The default number of iterations is 400, but this can be modified using the `--iterations` flag. One issue with this method is that the confidence in the weightings will be much higher for some trees than others, because of the nature of the binomial distribution. For example, if there are three possible topologies that occur at frequencies 400, 0 and 0, respectively, we can be highly confident of their weightings. But if they occur at frequencies of 120, 130 and 150, we have less certainty of their exact weightings.

The third option is `threshold`. This is similar to the sampling method above, except that the sampling is repeated until a certain dynamic threshold is reached. In other words, after a *n* iterations, each topology must have been observed **fewer** than *k* times OR **more** than *n-k* times. There is therefore a separate threshold dependng on the number of iterations and the number of observations. This allows you to be more thourough when the tree is more complex and the weightings less certain. For example, you can set the thresholds such that the 95% binomial confidence interval around each weighting is less than 5%. Thresholds must be specified in a file, using the flag `--thresholdTable`. Two example threshold tables are provided. These give thresholds to ensure that the 95% CI is less than 5% or 10%.




