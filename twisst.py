#!/usr/bin/env python

import argparse
import itertools
import sys
import gzip
import ete3
import random
import numpy as np
from collections import defaultdict
from collections import deque

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
    ## speeding up by first deleting all other leaves
    for leaf in pruned.iter_leaves():
        if leaf.name not in leavesToKeep: leaf.delete(preserve_branch_length=preserve_branch_length)
    #and then prune to fix the root (not sure why this is necessary, but it is)
    #but at least it's faster than pruning the full tree
    pruned.prune(leavesToKeep, preserve_branch_length = preserve_branch_length)
    return pruned

class NodeChain(deque):
    def __init__(self, nodeList, dists=None):
        super(NodeChain, self).__init__(nodeList)
        if dists is None: self.dists = None
        else:
            assert len(dists) == len(self)-1, "incorrect number of iternode distances"
            self.dists = deque(dists)
        self._set_ = None
    
    def addNode(self, name, dist=0):
        self.append(name)
        if self.dists is not None: self.dists.append(dist)
    
    def addNodeLeft(self, name, dist=0):
        self.appendleft(name)
        if self.dists is not None: self.dists.appendleft(dist)
    
    def addNodeChain(self, chainToAdd, joinDist=0):
        self.extend(chainToAdd)
        if self.dists is not None:
            assert chainToAdd.dists is not None, "Cannot add a chain without distances to one with distances"
            self.dists.append(joinDist)
            self.dists.extend(chainToAdd.dists)
    
    def addNodeChainLeft(self, chainToAdd, joinDist=0):
        self.extendleft(chainToAdd)
        if self.dists is not None:
            assert chainToAdd.dists is not None, "Cannot add a chain without distances to one with distances"
            self.dists.appendleft(joinDist)
            self.dists.extendleft(chainToAdd.dists)
    
    def chopLeft(self):
        self.popleft()
        if self.dists is not None: self.dists.popleft()
    
    def chop(self):
        self.pop()
        if self.dists is not None: self.dists.pop()
    
    def fuseLeft(self, chainToFuse):
        new = NodeChain(self, self.dists)
        assert new[0] == chainToFuse[0], "No common nodes"
        i = 1
        while new[1] == chainToFuse[i]:
            new.chopLeft()
            i += 1
        m = len(chainToFuse)
        while i < m:
            new.addNodeLeft(chainToFuse[i], chainToFuse.dists[i-1] if self.dists is not None else None)
            i += 1
        return new
    
    def simplifyToEnds(self, newDist=None):
        if self.dists is not None:
            if not newDist: newDist = sum(self.dists)
            self.dists.clear()
        leftNode = self.popleft()
        rightNode = self.pop()
        self.clear()
        self.append(leftNode)
        self.append(rightNode)
        if self.dists is not None:
            self.dists.append(newDist)
    
    def setSet(self):
        self._set_ = set(self)


##simpler version that only collapses monophyletic clades
#def getChainsToLeaves(node, collapseDict = None):
    #children = node.get_children()
    #if children == []:
        #node.add_feature("weight", 1)
        #return [NodeChain(node)]
    #chains = list(itertools.chain(*[getChainsToLeaves(child, collapseDict) for child in children]))
    #if (collapseDict and sum([len(chain) for chain in chains]) == len(chains) and
        #len(set([collapseDict[chain[0].name] for chain in chains])) == 1):
        ##all chains are a leaf from same group, so we collapse
        #newWeight = sum([chain[0].weight for chain in chains])
        #newDist = node.dist + sum([chain[0].dist * chain[0].weight * 1. for chain in chains]) / newWeight
        #chains[0][0].dist = newDist
        #chains[0][0].weight = newWeight
        #chains = [chains[0]]
    #else:
        #for chain in chains:
            #chain.addNodeLeft(node, dist=chain[0].dist)
    
    #return chains


def getChainsToLeaves(node, collapseDict = None, preserveDists = False):
    children = node.get_children()
    if children == []:
        #if it has no children is is a child, so just record a weight for the node and return is as a new 1-node chain
        chain = NodeChain([node], dists = [] if preserveDists else None)
        setattr(chain, "weight", 1)
        return [chain]
    #otherwise get chains for all children
    childrenChains = [getChainsToLeaves(child, collapseDict, preserveDists) for child in children]
    #now we have the chains from all children, we need to add the current node
    for childChains in childrenChains:
        for chain in childChains: chain.addNodeLeft(node, dist=chain[0].dist if preserveDists else None)
    
    #if collapsing, check groups for each node
    if collapseDict:
        nodeGroupsAll = np.array([collapseDict[chain[-1].name] for childChains in childrenChains for chain in childChains])
        nodeGroups = list(set(nodeGroupsAll))
        nGroups = len(nodeGroups)
        
        if (nGroups == 1 and len(nodeGroupsAll) > 1):
            #all leaves are from same group, so collapse to one chain
            #we can also preserve distances when doing this type of collapsing
            #first list all chains
            chains = [chain for childChains in childrenChains for chain in childChains]
            newWeight = sum([chain.weight for chain in chains])
            if preserveDists:
                newDist = sum([sum(chain.dists) * chain.weight * 1. for chain in chains]) / newWeight
                chains[0].simplifyToEnds(newDist = newDist)
            else:
                chains[0].simplifyToEnds()
            chains[0].weight = newWeight
            chains = [chains[0]]
        
        elif (nGroups == 2 and len(nodeGroupsAll) > 2 and preserveDists==False):
            #all chains end in a leaf from one of two groups, so we can simplify.
            #first list all chains
            chains = [chain for childChains in childrenChains for chain in childChains]
            #Start by getting index of each chain for each group
            indices = [(nodeGroupsAll == group).nonzero()[0] for group in nodeGroups]
            #the new weight for each chain we keep will be the total node weight of all from each group 
            newWeights = [sum([chains[i].weight for i in idx]) for idx in indices]
            #now reduce to just a chain for each group 
            chains = [chains[idx[0]] for idx in indices]
            for j in range(nGroups):
                chains[j].simplifyToEnds()
                chains[j].weight = newWeights[j]
        
        #if we couldn't simply collapse completely, we might still be able to merge down a side branch
        #Side branches are child chains ending in a single leaf
        #If there is a lower level child branch that is itself a side branch, we can merge to it
        elif (preserveDists == False and len(childrenChains) == 2 and
            ((len(childrenChains[0]) == 1 and len(childrenChains[1]) > 1) or
            (len(childrenChains[1]) == 1 and len(childrenChains[0]) > 1))):
            chains,sideChain = (childrenChains[1],childrenChains[0][0]) if len(childrenChains[0]) == 1 else (childrenChains[0],childrenChains[1][0])
            #now check if any main chain is suitable (should be length 3, and the only one that is such. and have correct group
            targets = (np.array([len(chain) for chain in chains]) == 3).nonzero()[0]
            if len(targets) == 1 and collapseDict[chains[targets[0]][-1].name] == collapseDict[sideChain[-1].name]:
                #we have found a suitable chain to merge to
                targetChain = chains[targets[0]]
                newWeight = targetChain.weight + sideChain.weight
                targetChain.simplifyToEnds()
                targetChain.weight = newWeight
            else:
                #if we didn't find a suitable match, just add side chain
                chains.append(sideChain)
        else:
            #if there was no side chain, just list all chains
            chains = [chain for childChains in childrenChains for chain in childChains]
    #otherwise we are not collapsing, so just list all chains
    else:
        #chains = list(itertools.chain(*[getChainsToLeaves(child, collapseDict) for child in children]))
        chains = [chain for childChains in childrenChains for chain in childChains]
    #now we have the chains from all children, we need to add the current node
    
    return chains

#version for tree sequence tree format from msprime and tsinfer
def getChainsToLeaves_ts(tree, node=None, collapseDict = None):
    if node is None: node = tree.root
    children = tree.children(node)
    if children == ():
        #if it has no children is is a child
        #if it's in the collapseDict or there is not collapseDict
        #just record a weight for the node and return is as a new 1-node chain
        if collapseDict is None or node in collapseDict:
            chain = NodeChain([node])
            setattr(chain, "weight", 1)
            return [chain]
        else:
            return []
    #otherwise get chains for all children
    childrenChains = [getChainsToLeaves_ts(tree, child, collapseDict) for child in children]
    #now we have the chains from all children, we need to add the current node
    for childChains in childrenChains:
        for chain in childChains: chain.addNodeLeft(node)
    
    #if collapsing, check groups for each node
    if collapseDict:
        nodeGroupsAll = np.array([collapseDict[chain[-1]] for childChains in childrenChains for chain in childChains])
        nodeGroups = list(set(nodeGroupsAll))
        nGroups = len(nodeGroups)
        
        if (nGroups == 1 and len(nodeGroupsAll) > 1) or (nGroups == 2 and len(nodeGroupsAll) > 2):
            #all chains end in a leaf from one or two groups, so we can simplify.
            #first list all chains
            chains = [chain for childChains in childrenChains for chain in childChains]
            #Start by getting index of each chain for each group
            indices = [(nodeGroupsAll == group).nonzero()[0] for group in nodeGroups]
            #the new weight for each chain we keep will be the total node weight of all from each group 
            newWeights = [sum([chains[i].weight for i in idx]) for idx in indices]
            #now reduce to just a chain for each group 
            chains = [chains[idx[0]] for idx in indices]
            for j in range(nGroups):
                chains[j].simplifyToEnds()
                chains[j].weight = newWeights[j]
        
        #if we couldn't simply collapse completely, we might still be able to merge down a side branch
        #Side branches are child chains ending in a single leaf
        #If there is a lower level child branch that is itself a side branch, we can merge to it
        elif (len(childrenChains) == 2 and
            ((len(childrenChains[0]) == 1 and len(childrenChains[1]) > 1) or
            (len(childrenChains[1]) == 1 and len(childrenChains[0]) > 1))):
            chains,sideChain = (childrenChains[1],childrenChains[0][0]) if len(childrenChains[0]) == 1 else (childrenChains[0],childrenChains[1][0])
            #now check if any main chain is suitable (should be length 3, and the only one that is such. and have correct group
            targets = (np.array([len(chain) for chain in chains]) == 3).nonzero()[0]
            if len(targets) == 1 and collapseDict[chains[targets[0]][-1]] == collapseDict[sideChain[-1]]:
                #we have found a suitable internal chain to merge to
                targetChain = chains[targets[0]]
                newWeight = targetChain.weight + sideChain.weight
                targetChain.simplifyToEnds()
                targetChain.weight = newWeight
            else:
                #if we didn't find a suitable match, just add side chain
                chains.append(sideChain)
        else:
            #if there was no side chain, just list all chains
            chains = [chain for childChains in childrenChains for chain in childChains]
    #otherwise we are not collapsing, so just list all chains
    else:
        chains = [chain for childChains in childrenChains for chain in childChains]
    #now we have the chains from all children, we need to add the current node
    
    return chains


def makeRootLeafChainDict(tree, collapseDict = None, preserveDists=False, treeFormat = "ete3"):
    if treeFormat == "ts":
        chains = getChainsToLeaves_ts(tree, collapseDict=collapseDict)
        return dict([(chain[-1],chain) for chain in chains])
    else:
        chains = getChainsToLeaves(tree, collapseDict=collapseDict, preserveDists=preserveDists)
        return dict([(chain[-1].name,chain) for chain in chains])

def makeLeafLeafChainDict(rootLeafChainDict, pairs):
    leafLeafChainDict = defaultdict(defaultdict)
    
    for pair in pairs:
        #get the leaf leaf chain by removing the unshared ancestors and joining root leaf chains end to end
        leafLeafChainDict[pair[0]][pair[1]] = rootLeafChainDict[pair[0]].fuseLeft(rootLeafChainDict[pair[1]])
    
    return leafLeafChainDict


def checkDisjointChains(leafLeafChains, pairsOfPairs, samples=None):
    if not samples:
        return [leafLeafChains[pairs[0][0]][pairs[0][1]]._set_.isdisjoint(leafLeafChains[pairs[1][0]][pairs[1][1]]._set_) for pairs in pairsOfPairs]
    else:
        return [leafLeafChains[samples[pairs[0][0]]][samples[pairs[0][1]]]._set_.isdisjoint(leafLeafChains[samples[pairs[1][0]]][samples[pairs[1][1]]]._set_) for pairs in pairsOfPairs]


def pairsDisjoint(pair1,pair2):
    return pair1[0] != pair2[0] and pair1[0] != pair2[1] and pair1[1] != pair2[0] and pair1[1] != pair2[1]

def makeTopoDict(taxonNames, topos=None, outgroup = None):
    output = {}
    output["topos"] = allTopos(taxonNames, []) if topos is None else topos
    if outgroup:
        for topo in output["topos"]: topo.set_outgroup(outgroup)
    output["n"] = len(output["topos"])
    pairs = list(itertools.combinations(taxonNames,2))
    pairsNumeric = list(itertools.combinations(range(len(taxonNames)),2))
    output["pairsOfPairs"] = [y for y in itertools.combinations(pairs,2) if pairsDisjoint(y[0],y[1])]
    output["pairsOfPairsNumeric"] = [y for y in itertools.combinations(pairsNumeric,2) if pairsDisjoint(y[0],y[1])]
    output["chainsDisjoint"] = []
    for tree in output["topos"]:
        rootLeafChains = makeRootLeafChainDict(tree)
        leafLeafChains = makeLeafLeafChainDict(rootLeafChains, pairs)
        for pair in pairs: leafLeafChains[pair[0]][pair[1]].setSet()
        output["chainsDisjoint"].append(checkDisjointChains(leafLeafChains, output["pairsOfPairs"]))
    return output

def makeGroupDict(groups, names=None):
    groupDict = {}
    for x in range(len(groups)):
        for y in groups[x]: groupDict[y] = x if not names else names[x]
    return groupDict

#Main weighting function that uses "chains" to check topologies and simplifies while generating chains
def weightTree(tree, taxa, taxonDict=None, pairs=None, topoDict=None, nIts=None,
                     getDists=False, simplify=True, abortCutoff=None, treeFormat="ete3", verbose=True,
                     taxonNames=None, outgroup=None):
    
    nTaxa = len(taxa)
    
    if not taxonDict: taxonDict = makeGroupDict(taxa)
    
    if pairs is None:
        pairs = [pair for taxPair in itertools.combinations(taxa,2) for pair in itertools.product(*taxPair)]
    
    rootLeafChains = makeRootLeafChainDict(tree, collapseDict=taxonDict if simplify else None, preserveDists=getDists, treeFormat=treeFormat)
    leavesRetained = rootLeafChains.keys()
    leavesRetainedSet = set(leavesRetained)
    leafWeights = dict([(ind, rootLeafChains[ind].weight) for ind in leavesRetained])
    _pairs = [pair for pair in pairs if pair[0] in leavesRetainedSet and pair[1] in leavesRetainedSet]
    leafLeafChains = makeLeafLeafChainDict(rootLeafChains, pairs=_pairs)
    #make a set for each chain so that 
    for pair in _pairs: leafLeafChains[pair[0]][pair[1]].setSet()
    
    if topoDict is None:
        if taxonNames is None: taxonNames = [str(x) for x in range(len(taxa))]
        topoDict = makeTopoDict(taxonNames, outgroup=outgroup)
    
    _taxa = [[ind for ind in taxon if ind in leavesRetainedSet] for taxon in taxa]
    
    if getDists:
        assert taxonNames is not None, "taxonNames required for recording pairwise distances"
        dists = np.zeros([nTaxa, nTaxa, topoDict["n"]])
    
    #we make a generator object for all combos
    nCombos = np.prod([len(t) for t in _taxa])
    #if not speciified assume all combinations must be considered
    if nIts is None: nIts = nCombos
    #if doing all combos, we use an exhaustive combo generator
    if nIts >= nCombos:
        if verbose: sys.stderr.write("Complete weighting for {} combinations\n".format(nCombos))
        #unless there are too many combos, in which case we abort
        if abortCutoff and nCombos > abortCutoff:
            if verbose: sys.stderr.write("Aborting\n")
            return None
        comboGenerator = itertools.product(*_taxa)
    #if we are doing a subset, then use a random combo generator, but make sure simplify was false
    else:
        #sys.stderr.write("Approximate weighting with {} combinations\n".format(nIts))
        assert not simplify, "Tree simplification should be turned off when considering only a subset of combinations."
        comboGenerator = randomComboGen(_taxa)
    
    #initialise counts array
    counts = [0]*topoDict["n"]
    i=0
    for combo in comboGenerator:
        chainsDisjoint = checkDisjointChains(leafLeafChains, topoDict["pairsOfPairsNumeric"], samples=combo)
        x = topoDict["chainsDisjoint"].index(chainsDisjoint)
        comboWeight = np.prod([leafWeights[ind] for ind in combo]) 
        counts[x] += comboWeight
        
        #get pairwise dists if necessary
        if getDists:
            comboPairs = [(combo[pair[0]], combo[pair[1]],) for pairs in topoDict["pairsOfPairsNumeric"] for pair in pairs]
            currentDists = np.zeros([nTaxa,nTaxa])
            for comboPair in comboPairs:
                taxPair = (taxonNames.index(taxonDict[comboPair[0]]), taxonNames.index(taxonDict[comboPair[1]]))
                currentDists[taxPair[0],taxPair[1]] = currentDists[taxPair[1],taxPair[0]] = sum(leafLeafChains[comboPair[0]][comboPair[1]].dists)
            dists[:,:,x] += currentDists*comboWeight
        
        i += 1
        if i == nIts: break
    
    meanDists = dists/counts if getDists else np.NaN
    return {"topos":topoDict["topos"], "weights":counts, "dists":meanDists}


def weightTrees(trees, taxa=None, taxonDict=None, pairs=None, topoDict=None, nIts=None,
                     getDists=False, simplify=True, abortCutoff=None, treeFormat="ete3", verbose=True,
                     taxonNames=None, outgroup=None):
    
    if taxa is None:
        assert(treeFormat=="ts"), "Taxa must be specified as a list of lists."
        if taxonNames is None: taxonNames = [str(pop.id) for pop in trees.populations()]
        taxa = [[s for s in trees.samples() if str(trees.get_population(s)) == t] for t in taxonNames]
        
    if topoDict is None:
        if taxonNames is None: taxonNames = [str(x) for x in range(len(taxa))]
        topoDict = makeTopoDict(taxonNames, outgroup=outgroup)
    
    if not taxonDict: taxonDict = makeGroupDict(taxa, names=taxonNames)
    
    if pairs is None:
        pairs = [pair for taxPair in itertools.combinations(taxa,2) for pair in itertools.product(*taxPair)]
    
    _trees_ = trees.trees() if treeFormat=="ts" else trees
    
    allTreeData = [weightTree(tree, taxa, taxonDict=taxonDict, pairs=pairs, topoDict=topoDict, nIts=nIts, getDists=getDists, simplify=simplify, abortCutoff=abortCutoff, treeFormat=treeFormat, verbose=verbose) for tree in _trees_]
    
    output = {"topos":allTreeData[0]["topos"]}
    output["dists"] = np.array([x["dists"] for x in allTreeData])
    output["weights"] = np.array([x["weights"] for x in allTreeData])
    output["weights_norm"] = np.apply_along_axis(lambda x: x/x.sum(), 1, output["weights"])
    
    return output


def summary(weightsData):
    if "weights_norm" not in weightsData:
        weights = np.apply_along_axis(lambda x: x/x.sum(), 1, weightsData["weights"])
    else:
        weights =weightsData["weights_norm"]
    meanWeights = weights.mean(axis=0)
    for i in range(len(meanWeights)):
        print("Topo", i+1)
        print(weightsData["topos"][i].get_ascii())
        print(round(meanWeights[i],3))
        print("\n\n")

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
    #print("topos contains", len(_topos), "topologies.")
    #print("current tree is:", branches)
    for x in range(len(branches)-1):
        for y in range(x+1,len(branches)):
            #print("Joining branch", x, branches[x], "with branch", y, branches[y])
            new_branches = list(branches)
            new_branches[x] = [new_branches[x],new_branches.pop(y)]
            #print("New tree is:", new_branches)
            if len(new_branches) == 3:
                #print("Tree has three branches, so appending to topos.")
                #now check that the topo doesn't match a topology already in trees, and if not add it
                t = ete3.Tree(listToNwk(new_branches))
                ID = t.get_topology_id()
                if ID not in _topo_IDs:
                    _topos.append(t)
                    _topo_IDs.add(ID)
            else:
                #print("Tree still unresolved, so re-calling function.")
                _topos = allTopos(new_branches, _topos, _topo_IDs)
    #print(_topo_IDs)
    return(_topos)

def writeWeights(filename, weightsData):
    nTopos = len(weightsData["topos"])
    with gzip.open(filename, "wt") if filename.endswith(".gz") else open(filename, "wt") as weightsFile:
        #write topologies
        for x in range(nTopos): weightsFile.write("#topo" + str(x+1) + " " + weightsData["topos"][x].write(format = 9) + "\n") 
        #write headers
        weightsFile.write("\t".join(["topo" + str(x+1) for x in range(nTopos)]) + "\n")
        #write weights
        weightsFile.write("\n".join(["\t".join(row) for row in weightsData["weights"].astype(str)]) + "\n")


def writeTsWindowData(filename, ts):
    with open("filename", "wt") as dataFile:
        dataFile.write("chrom\tstart\tend\n")
        dataFile.write("\n".join(["\t".join(["chr1", str(tree.interval[0]), str(tree.interval[1])]) for tree in ts.trees()]) + "\n")

#################################################################################################################################


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--treeFile", help="File containing tree(s) to analyse", action = "store")
    parser.add_argument("-w", "--weightsFile", help="Output file of all weights", action = "store")
    parser.add_argument("-D", "--distsFile", help="Output file of mean pairwise dists", action = "store", required = False)
    parser.add_argument("--inputTopos", help="Input file for user-defined topologies (optional)", action = "store", required = False)
    parser.add_argument("--outputTopos", help="Output file for topologies used", action = "store", required = False)
    parser.add_argument("--outgroup", help="Outgroup for rooting - only affects speed", action = "store")
    parser.add_argument("--method", help="Tree sampling method", choices=["fixed", "complete"], action = "store", default = "complete")
    parser.add_argument("--iterations", help="Number of iterations for fixed partial sampling", type=int, action = "store", default = 10000)
    parser.add_argument("--abortCutoff", help="# tips in simplified tree to abort 'complete' weighting", type=int, action = "store", default = 100000)
    parser.add_argument("-g", "--group", help="Group name and individual names (separated by commas)", action='append', nargs="+", required = True, metavar=("name","[inds]"))
    parser.add_argument("--groupsFile", help="Optional file of sample names and groups", action = "store", required = False)
    parser.add_argument("--verbose", help="Verbose output", action="store_true")
    parser.add_argument("--skip_simplify", help="", action="store_true")
    parser.add_argument("--silent", help="No stderr output", action="store_true")
    
    args = parser.parse_args()
    #args = parser.parse_args("-t examples/msms_4of10_l1Mb_r10k_sweep.seq_gen.SNP.w50sites.phyml_bionj.trees.gz  -g A 1,2,3,4,5,6,7,8,9,10 -g B 11,12,13,14,15,16,17,18,19,20 -g C 21,22,23,24,25,26,27,28,29,30 -g D 31,32,33,34,35,36,37,38,39,40".split())
    
    getDists  = args.distsFile is not None
    
    method = args.method
    
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
        with open(args.groupsFile, "rt") as gf: groupDict = dict([ln.split() for ln in gf.readlines()])
        for sample in groupDict.keys():
            try: taxa[taxonNames.index(groupDict[sample])].append(sample)
            except: pass
    
    nTaxa=len(taxa)
    
    assert min([len(t) for t in taxa]) >= 1, "Please specify at least one sample name per group."
    
    names = [t for taxon in taxa for t in taxon]
    namesSet = set(names)
    assert len(names) == len(namesSet), "Each sample should only be in one group."
    
    taxonDict = makeGroupDict(taxa, taxonNames)
    
    #get all topologies
    if args.inputTopos:
        with open(args.inputTopos, "rt") as tf: topos = [ete3.Tree(ln) for ln in tf.readlines()]
    else: topos = None
    
    topoDict = makeTopoDict(taxonNames, topos, args.outgroup if args.outgroup else None)
    
    nTopos = topoDict["n"]
    
    if not args.silent: sys.stderr.write(asciiTrees(topoDict["topos"],5) + "\n")
    
    if args.outputTopos:
        with open(args.outputTopos, "wt") as tf:
            tf.write("\n".join([t.write(format = 9) for t in topoDict["topos"]]) + "\n")
    
    pairs = [pair for taxPair in itertools.combinations(taxa,2) for pair in itertools.product(*taxPair)]
    
    #################################################################################################################################
    ### file for weights
    
    if args.weightsFile:
       weightsFile = gzip.open(args.weightsFile, "wt") if args.weightsFile.endswith(".gz") else open(args.weightsFile, "wt")
    else: weightsFile = sys.stdout
    
    for x in range(nTopos): weightsFile.write("#topo" + str(x+1) + " " + topoDict["topos"][x].write(format = 9) + "\n") 
    
    weightsFile.write("\t".join(["topo" + str(x+1) for x in range(nTopos)]) + "\n")
    
    ### file for lengths
    
    if getDists:
        if args.distsFile[-3:] == ".gz": distsFile = gzip.open(args.distsFile, "wt")
        else: distsFile = open(args.distsFile, "wt")
        for x in range(nTopos):
            distsFile.write("\t".join(["topo" + str(x+1) + "_" + "_".join(pair) for pair in itertools.combinations(taxonNames,2)]) + "\t")
        distsFile.write("\n")
    
    ################################################################################################################################
    
    #open tree file
    
    if args.treeFile:
        treeFile = gzip.open(args.treeFile, "rt") if args.treeFile.endswith(".gz") else open(args.treeFile, "rt")
    else: treeFile = sys.stdin
    
    ################################################################################################################################
    
    nTrees = 0
    
    for line in treeFile:
        
        tree = readTree(line)
        
        if tree:
            #remove unneccesary leaves (speeds up downstream steps)
            leafNamesSet = set([leaf.name for leaf in tree.get_leaves()])
            if namesSet != leafNamesSet:
                assert namesSet.issubset(leafNamesSet), "Named samples not present in tree:" + " ".join(list(namesSet.difference(leafNamesSet)))
                tree = getPrunedCopy(tree, leavesToKeep=names, preserve_branch_length=True)
            
            #root tree (this only helps speed up analysis, but does not change results)
            if args.outgroup: tree.set_outgroup(taxa[taxonNames.index(args.outgroup)][-1])
            
            weightsData = None
            
            if method == "complete":
                
                weightsData = weightTree(tree=tree, taxa=taxa, taxonDict=taxonDict,
                                                pairs=pairs, topoDict=topoDict, getDists=getDists,
                                                simplify=not args.skip_simplify, abortCutoff=args.abortCutoff, verbose=args.verbose, taxonNames=taxonNames)
            
            if method == "fixed" or weightsData == None:
                weightsData = weightTree(tree=tree, taxa=taxa, taxonDict=taxonDict,
                                                pairs=pairs, topoDict=topoDict, nIts=args.iterations, getDists=getDists,
                                                simplify=False, verbose=args.verbose, taxonNames=taxonNames)
            
            weightsLine = "\t".join([str(x) for x in weightsData["weights"]])
        
            if getDists:
                distsByTopo = []
                for x in range(nTopos):
                    distsByTopo.append("\t".join([str(round(weightsData["dists"][pair[0],pair[1],x], 4)) for pair in itertools.combinations(range(nTaxa), 2)]))
                distsLine = "\t".join(distsByTopo)
        
        else:
            if not args.silent: sys.stderr.write("Warning - failed to read tree.\n")
            weightsLine = "\t".join(["nan"]*nTopos)
            if getDists: distsLine = "\t".join(["nan"]*nTopos*len(list(itertools.combinations(range(nTaxa), 2))))
        
        weightsFile.write(weightsLine + "\n")
        if getDists: distsFile.write(distsLine + "\n")
        
        nTrees += 1
        
        if not args.silent:
            print(".", end="", file=sys.stderr, flush=True)
            if nTrees % 100 == 0: sys.stderr.write(str(nTrees)+"\n")
    
    treeFile.close()
    weightsFile.close()
    if getDists: distsFile.close()
    
    if not args.silent: sys.stderr.write(str(nTrees)+"\nDone.\n")
    
    sys.exit()


#############################################################################################################################################
