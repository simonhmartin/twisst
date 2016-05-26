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
Given a tree with any number of tips (or leaves), each belonging to a particular taxon,the weighting of each taxon topology is defined as the fraction of all unique sub-trees, in which each taxon is represented by a single tip, that match that topology. Topology weighting therefore reduces the complexity of the full tree to a number of values, each giving the proportionate contribution of a particular taxon tree to the full tree. 


This code implements the method *Twisst* (topology weighting by iterative sampling of sub-trees), which does what it says: it computes the weightings by iteratively sampling sub-trees from the full tree and checking their topology. This can be slow if there are many tips (e.g. 4 taxa with ten tips each gives 10 000 unique subtrees to consider. But there are some shortcuts to speed things up.

Firstly, monophyletic groups of samples from the same taxon can be collapsed an weighted appropriately. This dramatically reduced the number of unique subtrees to consider. However, in some cases there are still too many to compute the weighting in reasonable time. In these cases it's best to randomly sample a limited number of subtrees and take the estimates weightings. The errors are binomially distributed, so it's possible to estimate confidence intervals.

The main script, `twisst.py` runs in Pythion 2.7 and requires the library `ETE3`.

###Running

A typical command looks like this:

```bash
python twisst.py -t input.trees.gz -w output.weights.weights.csv.gz -o topologies.trees -g A 1,2,3,4,5 -B 6,7,8,9,10 -g C 11,12,13,14,15 -D 16,17,18,19,20 --method complete
```

You can get a full list ot command options with:
```bash
python twisst.py -h
``` 



####Input

The main input for twisst is one or more trees in newick format. Most variants of newick are accepted - see the `ETE3` documentation for details.

Multiple trees should be listed on separate lines.

All trees must contain all specified tip labels (see below.)

####Specifying taxa and individuals 

Taxa (groups) must be specified in the command line, using the `-g` flag. The name of the group must be given, and optionally the names of the individuals (tip labels) that it contains, separated by commas. Alternatively, a tab-delimited file of tip labels and their corresponding groups can be provided with the `--groupsFile` flag.


