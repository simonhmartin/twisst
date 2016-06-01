
import argparse
import itertools
import sys
import gzip
import operator
import ete3
import random
import numpy as np

np.seterr(divide='ignore', invalid="ignore")
##############################################################################################################################

def prod(iterable): return reduce(operator.mul, iterable, 1)

def sample(things, n = None, replace = False):
    if n == None: n = len(things)
    if replace == False: return random.sample(things,n)
    else: return [random.choice(things) for i in range(n)]

def randomComboGen(lists):
    while True: yield tuple(random.choice(l) for l in lists)



#new version that does not do simplification, but has a few improvements
def weightTree(tree, taxa, taxonNames, nIts, topos=None, getDists = False):
    nTaxa = len(taxonNames)
    taxonDict = {}
    for x in range(nTaxa):
        for y in taxa[x]:
            taxonDict[y] = taxonNames[x]
    leaves = tree.get_leaves()
    leafNames = [leaf.name for leaf in leaves]
    #we make a generator object for all combos
    #if there are more combos than Its, we need to make random samples
    nCombos = prod([len(t) for t in taxa])
    if nIts >= nCombos:
        comboGenerator = itertools.product(*taxa)
        nIts = nCombos
    else: comboGenerator = randomComboGen(taxa)
    if not topos: topos=allTrees(taxonNames, [])
    #unique id for each topo
    topoIDs = [t.get_topology_id() for t in topos]
    counts = np.zeros(len(topos))
    if getDists: dists = np.zeros([nTaxa, nTaxa, len(topos)])
    for iteration in xrange(nIts):
        combo = comboGenerator.next()
        pruned = tree.copy("newick")
        pruned.prune(combo, preserve_branch_length = getDists)
        for leaf in pruned.iter_leaves(): leaf.name = taxonDict[leaf.name]
        pruned.unroot()
        #get pairwise dists if necessary
        if getDists: 
            currentDists = np.zeros([nTaxa,nTaxa])
            for pair in itertools.combinations(range(nTaxa), 2):
                currentDists[pair[0],pair[1]] = currentDists[pair[1],pair[0]] = pruned.get_distance(taxonNames[pair[0]], taxonNames[pair[1]])
        #find topology match
        prunedID = pruned.get_topology_id()
        x = topoIDs.index(prunedID)
        counts[x] += 1
        if getDists: dists[:,:,x] += currentDists
    if getDists: meanDists = dists/counts
    else: meanDists = np.NaN
    return {"topos":topos,"weights":counts,"dists":meanDists}




#new version that does not do simplification, but has a few improvements
def weightTreeThreshold(tree, taxa, taxonNames, thresholdDict, topos=None, getDists = False):
    nTaxa = len(taxonNames)
    taxonDict = {}
    for x in range(nTaxa):
        for y in taxa[x]:
            taxonDict[y] = taxonNames[x]
    leaves = tree.get_leaves()
    leafNames = [leaf.name for leaf in leaves]
    #make random combinations for sampling
    comboGenerator = randomComboGen(taxa)
    if not topos: allTrees(taxonNames, [])
    #unique id for each topo
    topoIDs = [t.get_topology_id() for t in topos]
    counts = np.zeros(len(topos))
    if getDists: dists = np.zeros([nTaxa, nTaxa, len(topos)])
    total=0
    while True:
        combo = comboGenerator.next()
        pruned = tree.copy("newick")
        pruned.prune(combo, preserve_branch_length = getDists)
        for leaf in pruned.iter_leaves(): leaf.name = taxonDict[leaf.name]
        pruned.unroot()
        #get pairwise dists if necessary
        if getDists: 
            currentDists = np.zeros([nTaxa,nTaxa])
            for pair in itertools.combinations(range(nTaxa), 2):
                currentDists[pair[0],pair[1]] = currentDists[pair[1],pair[0]] = pruned.get_distance(taxonNames[pair[0]], taxonNames[pair[1]])
        #find topology match
        prunedID = pruned.get_topology_id()
        x = topoIDs.index(prunedID)
        counts[x] += 1
        if getDists: dists[:,:,x] += currentDists
        total+=1
        #now check is we're under the threshold
        minCounts = np.minimum(counts,total-counts)
        if total not in thresholdDict or np.all(minCounts<=thresholdDict[total]): break
    if getDists: meanDists = dists/counts
    else: meanDists = np.NaN
    return {"topos":topos,"weights":counts,"dists":meanDists}




#a test version - this is designed to show how the values approach the truth with increasing iterations
def weightTreeEachIter(tree, taxa, taxonNames, nIts, topos):
    nTaxa = len(taxonNames)
    taxonDict = {}
    for x in range(nTaxa):
        for y in taxa[x]:
            taxonDict[y] = taxonNames[x]
    leaves = tree.get_leaves()
    leafNames = [leaf.name for leaf in leaves]
    #we make a generator object for all combos
    #if there are more combos than Its, we need to make random samples
    nCombos = prod([len(t) for t in taxa])
    if nIts >= nCombos:
        comboGenerator = itertools.product(*taxa)
        nIts = nCombos
    else: comboGenerator = randomComboGen(taxa)
    #unique id for each topo
    topoIDs = [t.get_topology_id() for t in topos]
    topoIDsSet = set(topoIDs)
    counts = np.zeros([len(topos), nIts])
    for i in xrange(nIts):
        combo = comboGenerator.next()
        pruned = tree.copy("newick")
        pruned.prune(combo, preserve_branch_length = False)
        for leaf in pruned.iter_leaves(): leaf.name = taxonDict[leaf.name]
        pruned.unroot()
        #check for topology match
        prunedID = pruned.get_topology_id()
        x = topoIDs.index(prunedID)
        counts[x,i:] += 1
    total = np.sum(counts,0)
    weights =  counts/total
    return {"topos":topos,"counts":counts,"weights":weights}



''' below is a bunch of code for a potentially much faster method.
It first simplifies the tree, and weights each of the collapsed nodes
accordingly, and then does the weighting as before, but taking node weights into account.
It is super fast, and works well when trees are small, such that all combinations
can be tested, but when only a subset of combos are tested it is unreliable.
It seems to increase the variance dramatically, so in order to get a good result,
you need to test far more combos, so the efficiency improvement goes away. I couldn't
figure out if this makes sense, so I just dropped it for now.'''


def simplifyClade(tree, node):
    #will Collapse a node, but change its branch length to the mean of all decendents
    #get lengths to all decendents
    leaves = node.get_leaves()
    leafDists = [node.get_distance(leaf) for leaf in leaves]
    meanLeafDist = sum(leafDists)/len(leafDists)
    #now remove the children from this node
    for child in node.get_children():
        node.remove_child(child)
    #rename and add length
    node.name = leaves[0].name
    node.dist += meanLeafDist
    node.add_feature("weight", len(leaves))



def simplifyTree(tree,taxonDict):
    simpTree = tree.copy("newick")
    #remove all leaves that are not in the taxonDict
    simpTree.prune(taxonDict.keys(), preserve_branch_length = True)
    #we will check all nodes for monophyly
    for node in simpTree.traverse("levelorder"):
        #print "Checking node", node.name
        #check if node is in tree (not sure if this is necessary)
        if node not in simpTree: # probably not a necessary line if traversing from top down
            continue
        leafNames = node.get_leaf_names()
        #print len(leafNames), "leaves:"
        #print leafNames
        if len(leafNames) >= 2:
            sameTaxon = True
            #get taxon of first node, and then ckeck others
            taxon = taxonDict[leafNames[0]]
            for leafName in leafNames[1:]:
                if taxonDict[leafName] != taxon:
                    sameTaxon = False
                    break
            #if all same taxon, get mean branch length and collapse
            if sameTaxon:
                #print "simplifying node."
                simplifyClade(simpTree, node)
                #print "new name:", node.name, "weight:", node.weight
        elif node.is_leaf():
            node.add_feature("weight", 1)
    return simpTree


#new version that does the tree simplification
def weightTreeSimp(tree, taxa, taxonNames, topos = None):
    taxonDict = {}
    for x in range(len(taxa)):
        for y in taxa[x]:
            taxonDict[y] = taxonNames[x]
    
    #simplify the tree
    simpTree = simplifyTree(tree, taxonDict)
    leaves = simpTree.get_leaves()
    leafNames = [leaf.name for leaf in leaves]
    leafWeights = dict(zip(leafNames, [leaf.weight for leaf in leaves]))
    #leafWeights = dict(zip(leafNames, [1]*len(leaves)))
    simpTaxa = [[t for t in taxon if t in leafNames] for taxon in taxa]
    #we make a generator object for all combos
    nCombos = prod([len(t) for t in simpTaxa])
    comboGenerator = itertools.product(*simpTaxa)

    if not topos: allTrees(taxonNames, [])
    topoIDs = [t.get_topology_id() for t in topos]
    
    counts = np.zeros(len(topos))
    for combo in comboGenerator:
        comboWeight = prod(leafWeights[leafName] for leafName in combo)
        pruned = simpTree.copy("newick")
        pruned.prune(combo)
        #generify pruned tree
        for leaf in pruned.iter_leaves(): leaf.name = taxonDict[leaf.name]
        pruned.unroot()
        #check for topology match
        prunedID = pruned.get_topology_id()
        x = topoIDs.index(prunedID)
        counts[x] += comboWeight
    
    return {"topos":topos,"weights":counts}


def listToNwk(t):
    t = str(t)
    t = t.replace("[","(")
    t = t.replace("]",")")
    t = t.replace("'","")
    t += ";"
    return(t)

def allTrees(branches, trees = []):
    assert 4 <= len(branches) <= 8, "Please specify between 4 and 8 unique taxon names."
    #print "current tree is:", branches 
    for x in range(len(branches)-1):
        for y in range(x+1,len(branches)):
            #print "Joining branch", x, branches[x], "with branch", y, branches[y]
            new_branches = list(branches)
            new_branches[x] = [new_branches[x],new_branches.pop(y)]
            #print "New tree is:", new_branches
            if len(new_branches) == 3:
                #print "Tree has three branches, so appending to trees."
                #now check that the tree doesn't match a topology already in trees, and if not add it
                tree = ete3.Tree(listToNwk(new_branches))
                if len(trees) == 0 or min([tree.robinson_foulds(t, unrooted_trees = True)[0] for t in trees]) > 0:
                    trees.append(tree)
            else:
                #print "Tree still unresolved, so re-calling function."
                trees = allTrees(new_branches, trees)
    return(trees)


#################################################################################################################################


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--treeFile", help="File containing tree(s) to analyse", action = "store", required = True)
    parser.add_argument("-o", "--topoFile", help="Output file of all topologies", action = "store", required = True)
    parser.add_argument("-w", "--weightsFile", help="Output file of all weights", action = "store", required = True)
    parser.add_argument("-D", "--distsFile", help="Output file of mean pairwise dists", action = "store", required = False)
    parser.add_argument("--method", help="Tree sampling method", choices=["fixed", "threshold", "complete"], action = "store", default = "fixed")
    parser.add_argument("--iterations", help="Number of iterations for fixed partial sampling", type=int, action = "store", default = 400)
    parser.add_argument("--thresholdTable", help="Lookup_table_for_sampling_thresholds", action = "store")
    parser.add_argument("-g", "--group", help="Group name and individual names (separated by commas)", action='append', nargs="+", required = True, metavar=("name","[inds]"))
    parser.add_argument("--groupsFile", help="Optional file of sample names and groups", action = "store", required = False)
    parser.add_argument("--verbose", help="Verbose output", action="store_true")


    args = parser.parse_args()
    #args = parser.parse_args("-n 5 -t test.trees -o test.topos.txt -w test.weights.B.csv -g A a,b,c -g B d,e,f -g C g,h,i -g D j,k,l".split())

    treeFileName = args.treeFile

    topoFileName = args.topoFile

    weightsFileName = args.weightsFile    

    distsFileName = args.distsFile
    if distsFileName: getDists = True
    else: getDists = False

    method = args.method

    verbose = args.verbose

    #################################################################################################################################
    #parse taxa
    assert len(args.group) >= 4, "Please specify at least four groups."
    taxonNames = []
    taxa = []
    for g in args.group:
        taxonNames.append(g[0])
        if len(g) > 1: taxa.append(g[1].split(","))
        else: taxa.append([])
    
    if args.groupsFile:
        with open(args.groupsFile, "r") as gf: groupDict = dict([ln.split() for ln in gf.readlines()])
        for sample in groupDict.keys():
            try: taxa[taxonNames.index(groupDict[sample])].append(sample)
            except: pass
    
    assert min([len(t) for t in taxa]) >= 1, "Please specify at least one sample name per group."

    #get all topologies
    topos = allTrees(taxonNames, [])
    
    for topo in topos: print >> sys.stderr, topo


    #make a rooted set of topos, just for printing - this doesn't affect the analysis
    #toposRooted = [topo.copy("newick") for topo in topos]
    #for topo in toposRooted: topo.set_outgroup(taxonNames[-1])

    topoFile = open(topoFileName, "w")

    topoFile.write("\n".join([t.write(format = 9) for t in topos]) + "\n")

    topoFile.close()
    
    #################################################################################################################################
    
    # check method

    if method == "fixed":
        nIts = args.iterations
        if nIts >= prod([len(t) for t in taxa]):
            print >> sys.stderr, "Warning: number of iterations is equal or greater than possible combinations.\n"
            nIts = prod([len(t) for t in taxa])
            print >> sys.stderr, "This could be very slow. Use method 'complete' for fast(er) exhaustive sampling."
    
    elif method == "threshold":
        assert args.thresholdTable, "A threshold table must be provided using argument --thresholdTable."
        thresholdTableFileName = args.thresholdTable
        with open(thresholdTableFileName) as ttf:
            thresholdDict = dict([(int(tries),int(threshold)) for line in ttf.readlines() for tries,threshold in (line.split(),)])

    #################################################################################################################################
    ### file for weights

    if weightsFileName[-3:] == ".gz": weightsFile = gzip.open(weightsFileName, "w")
    else: weightsFile = open(weightsFileName, "w")

    weightsFile.write(",".join(["topo" + str(x) for x in range(len(topos))]) + "\n")

    ### file for lengths

    if getDists:
        if distsFileName[-3:] == ".gz": distsFile = gzip.open(distsFileName, "w")
        else: distsFile = open(distsFileName, "w")
        for x in range(len(topos)):
            distsFile.write("topo" + str(x) + "_".join([pair for pair in itertools.combinations(taxonNames,2)]))
        distsFile.write("\n")

    ################################################################################################################################

    #open tree file

    if treeFileName[-3:] == ".gz": treeFile = gzip.open(treeFileName, "r")
    else: treeFile = open(treeFileName, "r")

    line = treeFile.readline()
    
    ################################################################################################################################

    n = 0

    while len(line) >= 1:
        tree = ete3.Tree(line.rstrip())
        if method == "fixed":
            weightsData = weightTree(tree=tree, taxa=taxa, taxonNames=taxonNames, nIts=nIts, topos=topos, getDists=getDists)
        elif method == "threshold":
            weightsData = weightTreeThreshold(tree=tree, taxa=taxa, taxonNames=taxonNames, thresholdDict=thresholdDict, topos=topos, getDists=getDists)
        elif method == "complete":
            weightsData = weightTreeSimp(tree=tree, taxa=taxa, taxonNames=taxonNames, topos=topos)
        
        weightsFile.write(",".join([str(round(x,3)) for x in weightsData["weights"]]) + "\n")
        if getDists:
            for x in range(len(topos)):
                distsFile.write(",".join([str(weightsData["dists"][pair[0],pair[1]]) for pair in itertools.combinations(range(nTaxa, 2))]) + ",")
            distsFile.write("\n")
        print n
        n += 1
        line = treeFile.readline()

    treeFile.close()
    weightsFile.close()

    sys.exit()




#############################################################################################################################################



