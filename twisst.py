#!/usr/bin/env python

import argparse
import itertools
import sys
import gzip
import ete3
import random
import numpy as np

np.seterr(divide='ignore', invalid="ignore")

verbose = False
##############################################################################################################################

def sample(things, n = None, replace = False):
    if n == None: n = len(things)
    if replace == False: return random.sample(things,n)
    else: return [random.choice(things) for i in range(n)]

def randomComboGen(lists):
    while True: yield tuple(random.choice(l) for l in lists)

def readTree(newick_tree):
    try:
        if newick_tree[0] == "[": return ete3.Tree(newick_tree[newick_tree.index("]")+1:])
        else: return ete3.Tree(newick_tree)
    except:
        return None

def asciiTrees(trees, nColumns = 5):
    treeLines = [tree.get_ascii().split("\n") for tree in trees]
    maxLines = max(map(len,treeLines))
    for tl in treeLines:
        #add lines if needed
        tl += [""]*(maxLines - len(tl))
        #add spaces to each line to make all even
        lineLengths = map(len,tl)
        maxLen = max(lineLengths)
        for i in range(len(tl)): tl[i] += "".join([" "]*(maxLen-len(tl[i])))
    #now join lines that will be on the same row and print
    treeLinesChunked = [treeLines[x:(x+nColumns)] for x in range(0,len(trees),nColumns)]
    zippedLinesChunked = [zip(*chunk) for chunk in treeLinesChunked]
    return "\n\n".join(["\n".join(["    ".join(l) for l in chunk]) for chunk in zippedLinesChunked])


def getPrunedCopy(tree, leavesToKeep, preserve_branch_length):
    pruned = tree.copy("newick")
    ##prune function was too slow for big trees
    ## speeding up by first deleting all oter leaves
    for leaf in pruned.iter_leaves():
        if leaf.name not in leavesToKeep: leaf.delete(preserve_branch_length=preserve_branch_length)
    #and then prune to fix the root (not sure why this is necessary, but it is)
    #but at least it's faster than pruning the full tree
    pruned.prune(leavesToKeep, preserve_branch_length = preserve_branch_length)
    return pruned

def isRooted(tree): return len(tree.get_children()) == 2

def treeToEdges(tree):
    edges=[]
    for n in tree.traverse(strategy="postorder"):
        if n.up is not None: edges.append((n.name, n.up.name, n.dist,))
    return edges

def addNodeNames(tree):
    for n in tree.traverse():
        if n.name is "": n.name = "_".join(sorted(n.get_leaf_names()))


#new version that does not do simplification, but has a few improvements
def weightTree(tree, taxa, taxonNames, taxonDict=None, nIts=None, topos=None, getLengths = False, getDists = False):
    nTaxa = len(taxonNames)
    if not taxonDict:
        taxonDict = {}
        for x in range(nTaxa):
            for y in taxa[x]: taxonDict[y] = taxonNames[x]
    #we make a generator object for all combos
    #if there are more combos than Its, we need to make random samples
    nCombos = np.prod([len(t) for t in taxa])
    if nIts is None: nIts = nCombos
    if nIts >= nCombos:
        comboGenerator = itertools.product(*taxa)
        nIts = nCombos
    else: comboGenerator = randomComboGen(taxa)
    if not topos: topos=allTopos(taxonNames, [])
    nTopos = len(topos)
    for t in topos: t.set_outgroup(taxonNames[-1])
    #if getting branch lengths get parents, children and make an array for the lengths
    if getLengths:
        #nodes need names to allow branch comparisons
        for t in topos: addNodeNames(t)
        parentsAndChildren = [[(n.up.name,n.name,) for n in t.traverse() if n.up is not None] for t in topos]
        parents = [[pc[0] for pc in parentsAndChildren[x]] for x in range(nTopos)]
        children = [[pc[1] for pc in parentsAndChildren[x]] for x in range(nTopos)]
        lengths = [np.zeros(len(children[x])) for x in range(nTopos)]
    else: parents = children = None
    if getDists: dists = np.zeros([nTaxa, nTaxa, nTopos])
    #unique id for each topo
    topoIDs = [t.get_topology_id() for t in topos]
    counts = np.zeros(nTopos, dtype=int)
    for iteration in xrange(nIts):
        combo = comboGenerator.next()
        pruned = getPrunedCopy(tree, combo, preserve_branch_length = getDists)
        #rename leaves
        for leaf in pruned.iter_leaves(): leaf.name = taxonDict[leaf.name]
        pruned.set_outgroup(taxonNames[-1])
        #find topology match
        #pruned.unroot()
        prunedID = pruned.get_topology_id()
        x = topoIDs.index(prunedID)
        counts[x] += 1
        #get pairwise dists if necessary
        if getDists:
            currentDists = np.zeros([nTaxa,nTaxa])
            for pair in itertools.combinations(range(nTaxa), 2):
                currentDists[pair[0],pair[1]] = currentDists[pair[1],pair[0]] = pruned.get_distance(taxonNames[pair[0]], taxonNames[pair[1]])
            dists[:,:,x] += currentDists
        #if getting average lengths, add branch lengths to average topos 
        if getLengths:
            addNodeNames(pruned)
            lengthDict = dict([(n.name,n.dist,) for n in pruned.traverse()])
            lengths[x] += [lengthDict[name] for name in children[x]]
    meanDists = dists/counts if getDists else np.NaN
    meanLengths = [lengths[x]/counts[x] for x in range(nTopos)] if getLengths else np.NaN
    return {"topos":topos,"weights":counts,"dists":meanDists,"parents":parents,"children":children,"lengths":meanLengths}


#new version that does not do simplification, but has a few improvements
def weightTreeThreshold(tree, taxa, taxonNames, thresholdDict, taxonDict=None, topos=None, getLengths = False, getDists = False):
    nTaxa = len(taxonNames)
    if not taxonDict:
        taxonDict = {}
        for x in range(nTaxa):
            for y in taxa[x]: taxonDict[y] = taxonNames[x]
    #make random combinations for sampling
    comboGenerator = randomComboGen(taxa)
    if not topos: topos=allTopos(taxonNames, [])
    for t in topos: t.set_outgroup(taxonNames[-1])
    if getLengths:
        #nodes need names to allow branch comparisons
        for t in topos: addNodeNames(t)
    nTopos = len(topos)
    #unique id for each topo
    topoIDs = [t.get_topology_id() for t in topos]
    counts = np.zeros(nTopos, dtype=int)
    #if getting branch lengths get parents, children and make an array for the lengths
    if getLengths:
        parentsAndChildren = [[(n.up.name,n.name,) for n in t.traverse() if n.up is not None] for t in topos]
        parents = [[pc[0] for pc in parentsAndChildren[x]] for x in range(nTopos)]
        children = [[pc[1] for pc in parentsAndChildren[x]] for x in range(nTopos)]
        lengths = [np.zeros(len(children[x])) for x in range(nTopos)]
    else: parents = children = None
    if getDists: dists = np.zeros([nTaxa, nTaxa, nTopos])
    total=0
    while True:
        combo = comboGenerator.next()
        pruned = getPrunedCopy(tree, combo, preserve_branch_length = getDists)
        for leaf in pruned.iter_leaves(): leaf.name = taxonDict[leaf.name]
        #get pairwise dists if necessary
        pruned.set_outgroup(taxonNames[-1])
        if getDists: 
            currentDists = np.zeros([nTaxa,nTaxa])
            for pair in itertools.combinations(range(nTaxa), 2):
                currentDists[pair[0],pair[1]] = currentDists[pair[1],pair[0]] = pruned.get_distance(taxonNames[pair[0]], taxonNames[pair[1]])
        #find topology match
        prunedID = pruned.get_topology_id()
        x = topoIDs.index(prunedID)
        counts[x] += 1
        if getDists: dists[:,:,x] += currentDists
        #if getting average lengths, add branch lengths to average topos 
        if getLengths:
            addNodeNames(pruned)
            lengthDict = dict([(n.name,n.dist,) for n in pruned.traverse()])
            lengths[x] += [lengthDict[name] for name in children[x]]
        total+=1
        #now check is we're under the threshold
        minCounts = np.minimum(counts,total-counts)
        if total not in thresholdDict or np.all(minCounts<=thresholdDict[total]): break
    meanDists = dists/counts if getDists else np.NaN
    meanLengths = [lengths[x]/counts[x] for x in range(nTopos)] if getLengths else np.NaN
    return {"topos":topos,"weights":counts,"dists":meanDists,"parents":parents,"children":children,"lengths":meanLengths}




#a test version - this is designed to show how the values approach the truth with increasing iterations
def weightTreeEachIter(tree, taxa, taxonNames, nIts, topos):
    nTaxa = len(taxonNames)
    taxonDict = {}
    for x in range(nTaxa):
        for y in taxa[x]:
            taxonDict[y] = taxonNames[x]
    #we make a generator object for all combos
    #if there are more combos than Its, we need to make random samples
    nCombos = np.prod([len(t) for t in taxa])
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
        pruned = getPrunedCopy(tree, combo, preserve_branch_length = False)
        pruned.unroot()
        for leaf in pruned.iter_leaves(): leaf.name = taxonDict[leaf.name]
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
you need to test far more combos, so the efficiency improvement goes away.
So this is only used for "complete" weighting'''


def simplifyClade(node,preserve_branch_length=True):
    #will Collapse a node, but change its branch length to the mean of all decendents
    #get lengths to all decendents
    leaves = node.get_leaves()
    if preserve_branch_length:
        leafDists = [1.*node.get_distance(leaf) for leaf in leaves]
        meanLeafDist = sum(leafDists)/len(leafDists)
    #now remove the children from this node
    for child in node.get_children(): node.remove_child(child)
    #rename and add length
    node.name = leaves[0].name
    if preserve_branch_length: node.dist += meanLeafDist
    node.add_feature("weight", len(leaves))


def simplifyClade(node,preserve_branch_length=True):
    #will Collapse a node, but change its branch length to the mean of all decendents
    #get lengths to all decendents
    leaves = node.get_leaves()
    if preserve_branch_length:
        leafDists = [1.*node.get_distance(leaf) for leaf in leaves]
        meanLeafDist = sum(leafDists)/len(leafDists)
    #now remove the children from this node
    for child in node.get_children(): node.remove_child(child)
    #rename and add length
    node.name = leaves[0].name
    if preserve_branch_length: node.dist += meanLeafDist
    node.add_feature("weight", len(leaves))


'''If a leaf's nearest outgroup is a leaf from the same taxon, the inner leaf can be removed.
Example in the tree (C1,(B1,(B2,A1))); B2 can be removed, because A1's closest relative will stlill be from taxon B.
We remove B2, and correct branch lengths by moving the base downward (ie increase their  base), moving A1 upward,
and giving B1 the average length of B1 and B2.'''
def removeInnerLeaf(base, outerLeaf, innerNode, innerLeaf,preserve_branch_length=True):
    #correct lengths
    relativeWeight = 1.*innerLeaf.weight/(innerLeaf.weight+outerLeaf.weight)
    base.dist += innerNode.dist * relativeWeight
    innerNode.dist = innerNode.dist * (1-relativeWeight)
    outerLeaf.dist = innerLeaf.dist*relativeWeight + outerLeaf.dist*(1-relativeWeight)
    #transfer weight
    outerLeaf.weight += innerLeaf.weight
    #delete
    innerLeaf.delete(preserve_branch_length=preserve_branch_length)

'''If a pair of nodes both have two children, and those children are in the same taxa, then the nodes are twins.
We can remove one (along with it's children) after transfering the weights and averaging the lengths'''

'''I COULD NOT GET THIS STEP TO GIVE CORRECT DISTS - SO I GAVE UP'''

def simplifyTwins(nodeAB, nodeA, nodeB, leafA0, leafA1, leafB0, leafB1):
    
    WA0=leafA0.weight
    WA1=leafA1.weight
    WB0=leafB0.weight
    WB1=leafB1.weight
    
    leafB0.delete()
    leafB1.delete()
    nodeB.delete()
    nodeAB.delete()

    leafA0.weight = WA0+WB0
    leafA1.weight = WA1+WB1    

'''An EXPERIMENTAL even more drastic simplicifcation can be achieved by simplifying any clade that has
only two taxa represented. Function takes the node and two sets of decendent leaves.

'''

#def simplifyCladeTo2Leaves(node, leavesA, leavesB):
    ##collapse clade to two nodes, weighted accordingly
    ##get mean dist for each leaf, and mean dist among leaves
    #DA = np.mean([1.*node.get_distance(leaf) for leaf in leavesA])
    #DB = np.mean([1.*node.get_distance(leaf) for leaf in leavesB])
    #DAB = np.mean([1.*leafA.get_distance(leafB) for leafA,leafB in itertools.product(leavesA,leavesB)])
    
    ##compute new distances for the new nodes
    #newAdist = (DAB - DB + DA)/2
    #newBdist = DB - (DA - newAdist)
    #newNdist = DA - newAdist
    
    ##now remove the children from this node
    #for child in node.get_children(): node.remove_child(child)
    
    ##add the first leaf in each set as a child (arbitrary)
    #node.add_child(leavesA[0])
    #node.add_child(leavesB[0])
    
    ##set new distances
    #node.dist = newNdist
    #leavesA[0].dist = newAdist
    #leavesB[0].dist = newBdist
    
    ##add weights
    #leavesA[0].add_feature("weight", len(leavesA))
    #leavesB[0].add_feature("weight", len(leavesB))


##experimental version of simplifyTree Jan 2018, uses new fucnton to simplify clades with two taxa
#def simplifyTree(tree,taxonDict,preserve_branch_length=True):
    #simpTree = tree.copy("newick")
    ##remove all leaves that are not in the taxonDict
    #simpTree.prune(taxonDict.keys(), preserve_branch_length = True)
    ##we will check all nodes for monophyly
    #for node in simpTree.traverse("levelorder"):
        ##print "Checking node", node.name
        ##check if node is in tree (not sure if this is necessary)
        #if node not in simpTree: continue # probably not a necessary line if traversing from top down
        #leaves = node.get_leaves()
        ##print len(leafNames), "leaves:"
        ##print leafNames
        #if len(leaves) >= 2:
            #leafNames = [l.name for l in leaves]
            #leafTaxa = [taxonDict[leafName] for leafName in leafNames]
            #taxaPresent = list(set(leafTaxa))
            #if len(taxaPresent) == 1:
                ##print "simplifying 1-taxon node."
                #simplifyClade(node)
                ##print "new name:", node.name, "weight:", node.weight
            #elif len(taxaPresent) == 2:
                #leavesA = [leaves[x] for x in range(len(leaves)) if leafTaxa[x] == taxaPresent[0]]
                #leavesB = [leaves[x] for x in range(len(leaves)) if leafTaxa[x] == taxaPresent[1]]
                ##print "simplifying 2-taxon node."
                #simplifyCladeTo2Leaves(node, leavesA, leavesB)
        #elif node.is_leaf():
            #node.add_feature("weight", 1)
    
    #return simpTree


def simplifyTree(tree,taxonDict,preserve_branch_length=True):
    simpTree = tree.copy("newick")
    #remove all leaves that are not in the taxonDict
    simpTree.prune(taxonDict.keys(), preserve_branch_length = True)
    #we will check all nodes for monophyly
    for node in simpTree.traverse("levelorder"):
        #print "Checking node", node.name
        #check if node is in tree (not sure if this is necessary)
        if node not in simpTree: continue # probably not a necessary line if traversing from top down
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
                simplifyClade(node)
                #print "new name:", node.name, "weight:", node.weight
        elif node.is_leaf():
            node.add_feature("weight", 1)
    #now we can do a second and third step.
    #Working up from tips, we remove any redundant leaves and transfer their weight and correct the lengths
    #We aslo check whether each node has twin pairs below it.
    #if so, we remove one and pass weights to the other, and average lengths.
    #this all gets run over and over until no new changes are made
    changed = True
    while changed:
        changed = False
        for node in simpTree.traverse("postorder"):
            cdn0 = node.get_children()
            if len(cdn0) == 2:
                if cdn0[0].is_leaf():
                    #print "testing", cdn0[0].name
                    cdn1 = cdn0[1].get_children()
                    if len(cdn1) == 2:
                        if cdn1[0].is_leaf() and taxonDict[cdn0[0].name] == taxonDict[cdn1[0].name]:
                            #print "removing", cdn1[0].name
                            removeInnerLeaf(base=node, outerLeaf=cdn0[0], innerNode=cdn0[1], innerLeaf=cdn1[0])
                            changed = True
                        elif cdn1[1].is_leaf() and taxonDict[cdn0[0].name] == taxonDict[cdn1[1].name]:
                            #print "removing", cdn1[1].name
                            removeInnerLeaf(base=node, outerLeaf=cdn0[0], innerNode=cdn0[1], innerLeaf=cdn1[1])
                            changed = True
                elif cdn0[1].is_leaf():
                    #print "testing", cdn0[1].name
                    cdn1 = cdn0[0].get_children()
                    if len(cdn1) == 2:
                        if cdn1[0].is_leaf() and taxonDict[cdn0[1].name] == taxonDict[cdn1[0].name]:
                            #print "removing", cdn1[0].name
                            removeInnerLeaf(base=node, outerLeaf=cdn0[1], innerNode=cdn0[0], innerLeaf=cdn1[0])
                            changed = True
                        elif cdn1[1].is_leaf() and taxonDict[cdn0[1].name] == taxonDict[cdn1[1].name]:
                            #print "removing", cdn1[1].name
                            removeInnerLeaf(base=node, outerLeaf=cdn0[1], innerNode=cdn0[0], innerLeaf=cdn1[1])
                            changed = True
                #if we get here the pair are both not leaves, but we can check if they're twins
                #This code only simplifies pair twins, but in theory all twins could be simplified'''
                #DOES NOT PRESERVE BRANCH LENGTH
                elif len(node.get_leaves()) == 4 and not preserve_branch_length:
                    #print "Checking for twins"
                    lvsA,lvsB = [child.get_leaves() for child in cdn0]
                    if len(lvsA) == len(lvsB) == 2:
                        taxaA = [taxonDict[lf.name] for lf in lvsA]
                        taxaB = [taxonDict[lf.name] for lf in lvsB]
                        if taxaA == taxaB:
                            #print "removing", lvsB[0].name, "and", lvsB[1].name 
                            simplifyTwins(nodeAB=node,nodeA=cdn0[0],nodeB=cdn0[1],leafA0=lvsA[0],leafA1=lvsA[1],leafB0=lvsB[0],leafB1=lvsB[1])
                            changed = True
                        elif taxaA[0] == taxaB[1] and taxaA[1] == taxaB[0]:
                            #print "removing", lvsB[1].name, "and", lvsB[0].name 
                            simplifyTwins(nodeAB=node,nodeA=cdn0[0],nodeB=cdn0[1],leafA0=lvsA[0],leafA1=lvsA[1],leafB0=lvsB[1],leafB1=lvsB[0])
                            changed = True
    
    return simpTree


##a fake version of simpTree for testing, that does nothing but add weight to leaves
#def simplifyTree(tree,taxonDict):
    #simpTree = tree.copy("newick")
    ##remove all leaves that are not in the taxonDict
    #for node in simpTree.traverse("levelorder"):
        #if node.is_leaf():
            #node.add_feature("weight", 1)
    
    #return simpTree


#new version that does the tree simplification
def weightTreeSimp(tree, taxa, taxonNames, taxonDict=None, topos = None,
                   getLengths=False, getDists = False, abortCutoff = 100000, simpTreeWeightsDict = None):
    #assert isRooted(tree), "Tree must be rooted"
    nTaxa = len(taxonNames)
    if not topos: topos = allTopos(taxonNames, [])
    for t in topos: t.set_outgroup(taxonNames[-1])
    if getLengths:
        #nodes need names to allow branch comparisons
        for t in topos: addNodeNames(t)
    nTopos = len(topos)
    topoIDs = [t.get_topology_id() for t in topos]
    nCombosTotal = np.prod([len(t) for t in taxa])
    if not taxonDict:
        taxonDict = {}
        for x in range(nTaxa):
            for y in taxa[x]: taxonDict[y] = taxonNames[x]
    #simplify the tree - if getting dists preserve branch length (slower)
    simpTree = simplifyTree(tree, taxonDict, preserve_branch_length=getDists)
    if verbose: print >> sys.stderr, simpTree
    leaves = simpTree.get_leaves()
    if verbose: print >> sys.stderr, "Simplified tree has", len(leaves), "leaves."
    leafNames = [leaf.name for leaf in leaves]
    leafWeights = dict(zip(leafNames, [leaf.weight for leaf in leaves]))
    #leafWeights = dict(zip(leafNames, [1]*len(leaves)))
    simpTaxa = [[t for t in taxon if t in leafNames] for taxon in taxa]
    #if possible, first check if there are already weights available in the dict
    if simpTreeWeightsDict != None and not getDists and not getLengths:
        simpTreeGeneric = simpTree.copy("newick")
        for leaf in simpTreeGeneric.iter_leaves(): leaf.name = taxonDict[leaf.name]
        if verbose: print >> sys.stderr, simpTreeGeneric
        simpTreeID = simpTreeGeneric.get_topology_id()
        if simpTreeID in simpTreeWeightsDict:
            if verbose: print >> sys.stderr, "using recorded counts for", simpTreeID
            counts = simpTreeWeightsDict[simpTreeID]
            return {"topos":topos,"weights":counts,"dists":np.NaN,"parents":None,"children":None,"lengths":np.NaN}
   #we make a generator object for all combos
    nCombos = np.prod([len(t) for t in simpTaxa])
    if nCombos > abortCutoff: return None
    if verbose: print >> sys.stderr, "There are", nCombos, "combinations to test."
    comboGenerator = itertools.product(*simpTaxa)
    counts = np.zeros(nTopos, dtype=int)
    #if getting branch lengths get parents, children and make an array for the lengths
    if getLengths:
        parentsAndChildren = [[(n.up.name,n.name,) for n in t.traverse() if n.up is not None] for t in topos]
        parents = [[pc[0] for pc in parentsAndChildren[x]] for x in range(nTopos)]
        children = [[pc[1] for pc in parentsAndChildren[x]] for x in range(nTopos)]
        lengths = [np.zeros(len(children[x])) for x in range(nTopos)]
    else: parents = children = None
    if getDists: dists = np.zeros([nTaxa, nTaxa, nTopos])
    for combo in comboGenerator:
        comboWeight = np.prod([leafWeights[leafName] for leafName in combo])
        pruned = getPrunedCopy(simpTree, combo, preserve_branch_length = getDists)
        #generify pruned tree
        for leaf in pruned.iter_leaves(): leaf.name = taxonDict[leaf.name]
        pruned.set_outgroup(taxonNames[-1])

        #check for topology match
        prunedID = pruned.get_topology_id()
        x = topoIDs.index(prunedID)
        counts[x] += comboWeight
                
        #if getting dists
        if getDists:
            currentDists = np.zeros([nTaxa,nTaxa])
            for pair in itertools.combinations(range(nTaxa), 2):
                currentDists[pair[0],pair[1]] = currentDists[pair[1],pair[0]] = pruned.get_distance(taxonNames[pair[0]], taxonNames[pair[1]])
            dists[:,:,x] += currentDists*comboWeight
        
        #if getting average lengths, add branch lengths to average topos 
        if getLengths:
            addNodeNames(pruned)
            lengthDict = dict([(n.name,n.dist,) for n in pruned.traverse()])
            lengths[x] += [lengthDict[name]*comboWeight for name in children[x]]
    
    #add to dict if present
    if simpTreeWeightsDict != None and not getDists and not getLengths:
        simpTreeGeneric = simpTree.copy("newick")
        for leaf in simpTreeGeneric.iter_leaves(): leaf.name = taxonDict[leaf.name]
        simpTreeID = simpTreeGeneric.get_topology_id()
        simpTreeWeightsDict[simpTreeID] = counts
        if verbose: print >> sys.stderr, simpTreeID, "recorded" 
    
    meanDists = dists/counts if getDists else np.NaN
    meanLengths = [lengths[x]/counts[x] for x in range(nTopos)] if getLengths else np.NaN
    
    return {"topos":topos,"weights":counts,"dists":meanDists,"parents":parents,"children":children,"lengths":meanLengths}




def listToNwk(t):
    t = str(t)
    t = t.replace("[","(")
    t = t.replace("]",")")
    t = t.replace("'","")
    t += ";"
    return(t)

def allTopos(branches, _topos=None, _topo_IDs=None):
    if _topos is None or _topo_IDs is None:
        _topos = []
        _topo_IDs = set([])
    assert 4 <= len(branches) <= 8, "Please specify between 4 and 8 unique taxon names."
    #print "topos contains", len(_topos), "topologies."
    #print "current tree is:", branches 
    for x in range(len(branches)-1):
        for y in range(x+1,len(branches)):
            #print "Joining branch", x, branches[x], "with branch", y, branches[y]
            new_branches = list(branches)
            new_branches[x] = [new_branches[x],new_branches.pop(y)]
            #print "New tree is:", new_branches
            if len(new_branches) == 3:
                #print "Tree has three branches, so appending to topos."
                #now check that the topo doesn't match a topology already in trees, and if not add it
                t = ete3.Tree(listToNwk(new_branches))
                ID = t.get_topology_id()
                if ID not in _topo_IDs:
                    _topos.append(t)
                    _topo_IDs.add(ID)
            else:
                #print "Tree still unresolved, so re-calling function."
                _topos = allTopos(new_branches, _topos, _topo_IDs)
    #print _topo_IDs
    #print [t.write(format=9) for t in _topos]
    return(_topos)


#################################################################################################################################


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--treeFile", help="File containing tree(s) to analyse", action = "store")
    parser.add_argument("-w", "--weightsFile", help="Output file of all weights", action = "store")
    parser.add_argument("-D", "--distsFile", help="Output file of mean pairwise dists", action = "store", required = False)
    parser.add_argument("-L", "--lengthsFile", help="Output file of average branch lengths", action = "store", required = False)
    parser.add_argument("--inputTopos", help="Input file for user-defined topologies (optional)", action = "store", required = False)
    parser.add_argument("--outputTopos", help="Output file for topologies used", action = "store", required = False)
    parser.add_argument("--outgroup", help="Outgroup for rooting - only affects speed", action = "store")
    parser.add_argument("--method", help="Tree sampling method", choices=["fixed", "threshold", "complete"], action = "store", default = "fixed")
    parser.add_argument("--backupMethod", help="Backup method if aborting complete", choices=["fixed", "threshold"], action = "store", default = "fixed")
    parser.add_argument("--iterations", help="Number of iterations for fixed partial sampling", type=int, action = "store", default = 400)
    parser.add_argument("--abortCutoff", help="# tips in simplified tree to abort 'complete' weighting", type=int, action = "store", default = 100000)
    parser.add_argument("--thresholdTable", help="Lookup_table_for_sampling_thresholds", action = "store")
    parser.add_argument("--noRecord", help="Do not use prior results wind doing complete", action = "store_true")
    parser.add_argument("-g", "--group", help="Group name and individual names (separated by commas)", action='append', nargs="+", required = True, metavar=("name","[inds]"))
    parser.add_argument("--groupsFile", help="Optional file of sample names and groups", action = "store", required = False)
    parser.add_argument("--verbose", help="Verbose output", action="store_true")
    parser.add_argument("--silent", help="No stderr output", action="store_true")


    args = parser.parse_args()
    #args = parser.parse_args("-t examples/msms_4of10_l1Mb_r10k_sweep.seq_gen.SNP.w50sites.phyml_bionj.trees.gz  -g A 1,2,3,4,5,6,7,8,9,10 -g B 11,12,13,14,15,16,17,18,19,20 -g C 21,22,23,24,25,26,27,28,29,30 -g D 31,32,33,34,35,36,37,38,39,40 --method complete > test.weights".split())

    getDists  = args.distsFile is not None
    getLengths = args.lengthsFile is not None

    method = args.method
    backupMethod = args.backupMethod
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
    
    nTaxa=len(taxa)
    
    assert min([len(t) for t in taxa]) >= 1, "Please specify at least one sample name per group."
    
    names = [t for taxon in taxa for t in taxon]
    namesSet = set(names)
    assert len(names) == len(namesSet), "Each sample should only be in one group."
    
    taxonDict = {}
    for x in range(nTaxa):
        for y in taxa[x]: taxonDict[y] = taxonNames[x]
    
    #get all topologies
    if args.inputTopos:
        with open(args.inputTopos, "r") as tf: topos = [ete3.Tree(ln) for ln in tf.readlines()]
    else: topos = allTopos(taxonNames, [])
    
    nTopos = len(topos)
    
    for topo in topos: topo.set_outgroup(taxonNames[-1])

    sys.stderr.write(asciiTrees(topos,5) + "\n")

    if args.outputTopos:
        with open(args.outputTopos, "w") as tf:
            tf.write("\n".join([t.write(format = 9) for t in topos]) + "\n")
    
    #################################################################################################################################
    
    # check method

    if method == "fixed" or (method == "complete" and backupMethod == "fixed"):
        nIts = args.iterations
        if nIts >= np.prod([len(t) for t in taxa]):
            print >> sys.stderr, "Warning: number of iterations is equal or greater than possible combinations.\n"
            nIts = np.prod([len(t) for t in taxa])
            print >> sys.stderr, "This could be very slow. Use method 'complete' for fast(er) exhaustive sampling."
    
    elif method == "threshold" or (method == "complete" and backupMethod == "threshold"):
        assert args.thresholdTable, "A threshold table must be provided using argument --thresholdTable."
        thresholdTableFileName = args.thresholdTable
        with open(thresholdTableFileName) as ttf:
            thresholdDict = dict([(int(tries),int(threshold)) for line in ttf.readlines() for tries,threshold in (line.split(),)])
    
    
    if method == "complete" and not getDists and not getLengths and not args.noRecord:
        simpTreeWeightsDict = {}
    else: simpTreeWeightsDict = None
    
    #################################################################################################################################
    ### file for weights

    if args.weightsFile:
       weightsFile = gzip.open(args.weightsFile, "w") if args.weightsFile.endswith(".gz") else open(args.weightsFile, "w")
    else: weightsFile = sys.stdout
    
    for x in range(nTopos): weightsFile.write("#topo" + str(x+1) + " " + topos[x].write(format = 9) + "\n") 

    weightsFile.write("\t".join(["topo" + str(x+1) for x in range(nTopos)]) + "\n")

    ### file for lengths

    if getDists:
        if args.distsFile[-3:] == ".gz": distsFile = gzip.open(args.distsFile, "w")
        else: distsFile = open(args.distsFile, "w")
        for x in range(nTopos):
            distsFile.write("\t".join(["topo" + str(x+1) + "_" + "_".join(pair) for pair in itertools.combinations(taxonNames,2)]) + "\t")
        distsFile.write("\n")

    if getLengths:
        lengthsFile = gzip.open(args.lengthsFile, "w") if args.lengthsFile[-3:] == ".gz" else open(args.lengthsFile, "w")
        #name internal nodes
        for t in topos: addNodeNames(t)
        parentsAndChildren = [[(n.up.name,n.name,) for n in t.traverse() if n.up is not None] for t in topos]
        children = [[pc[1] for pc in parentsAndChildren[x]] for x in range(nTopos)]
        for x in range(nTopos): lengthsFile.write("#" + "topo" + str(x+1) + "\t" + " ".join(["--".join(pair) for pair in parentsAndChildren[x]]) + "\n")
        #lengthsFile.write("\t".join(["\t".join(["topo" + str(x+1) + "_" + nodeName for nodeName in children[x]]) for x in range(nTopos)]) + "\n")
        
    
    ################################################################################################################################

    #open tree file
    
    if args.treeFile:
        treeFile = gzip.open(args.treeFile, "r") if args.treeFile.endswith(".gz") else open(args.treeFile, "r")
    else: treeFile = sys.stdin

    line = treeFile.readline()
    
    ################################################################################################################################

    nTrees = 0

    while len(line) >= 1:

        tree = readTree(line)
        
        if tree:
            #remove unneccesary leaves (speeds up downstream steps)
            leafNamesSet = set([leaf.name for leaf in tree.get_leaves()])
            if namesSet != leafNamesSet:
                assert namesSet.issubset(leafNamesSet), "Named samples not present in tree:" + " ".join(list(namesSet.difference(leafNamesSet)))
                tree = getPrunedCopy(tree, leavesToKeep=names, preserve_branch_length=True)
            
            #root if necessary
            if args.outgroup: tree.set_outgroup(args.outgroup)
            
            weightsData = None
            
            if method == "complete":
                weightsData = weightTreeSimp(tree=tree, taxa=taxa, taxonNames=taxonNames, taxonDict=taxonDict,
                                            topos=topos, getLengths=getLengths, getDists=getDists, abortCutoff=args.abortCutoff,
                                            simpTreeWeightsDict=simpTreeWeightsDict)
            
            if method == "fixed" or (method == "complete" and backupMethod == "fixed" and weightsData == None):
                weightsData = weightTree(tree=tree, taxa=taxa, taxonNames=taxonNames, taxonDict=taxonDict,
                                         nIts=nIts, topos=topos, getLengths=getLengths, getDists=getDists)
            
            if method == "threshold" or (method == "complete" and backupMethod == "threshold" and weightsData == None):
                weightsData = weightTreeThreshold(tree=tree, taxa=taxa, taxonNames=taxonNames, thresholdDict=thresholdDict, taxonDict=taxonDict,
                                                  topos=topos, getLengths=getLengths, getDists=getDists)
                    
            weightsLine = "\t".join([str(x) for x in weightsData["weights"]])
        
            if getDists:
                distsByTopo = []
                for x in range(nTopos):
                    distsByTopo.append("\t".join([str(round(weightsData["dists"][pair[0],pair[1],x], 4)) for pair in itertools.combinations(range(nTaxa), 2)]))
                distsLine = "\t".join(distsByTopo)    
            if getLengths:
                lengthsLine = "\t".join(["\t".join([str(round(l,4)) for l in weightsData["lengths"][x]]) for x in range(nTopos)])
        else:
            weightsLine = "\t".join(["nan"]*nTopos)
            if getDists: distsLine = "\t".join(["nan"]*nTopos*len(list(itertools.combinations(range(nTaxa), 2))))
            if getLengths: lengthsLine = "\t".join(["nan"]*sum(map(len,children)))
        
        weightsFile.write(weightsLine + "\n")
        if getDists: distsFile.write(distsLine + "\n")
        if getLengths: lengthsFile.write(lengthsLine + "\n")
        
        if not args.silent: sys.stderr.write(str(nTrees)+"\n")
        nTrees += 1
        line = treeFile.readline()

    treeFile.close()
    weightsFile.close()
    if getDists: distsFile.close()
    if getLengths: lengthsFile.close()

    sys.exit()




#############################################################################################################################################



