
import matplotlib.pyplot as plt

import ete3

import numpy as np

np.seterr(divide='ignore', invalid='ignore')

import wquantiles

import itertools,argparse,gzip,sys


def getMidPos_method2(pos1, pos2, depth1, depth2, depth3): return((1.*pos2*(depth2-depth1) + pos1*(depth3-depth1))/(depth3-2*depth1+depth2))

def getNodePos(node, method = 1):
    children = node.get_children()
    if method == 1:
        c1,c2 = children
        assert len(children)==2, "Position method one only works for bifurcating nodes."
        return((1.*c2.pos*(c1.depth-node.depth)+c1.pos*(c2.depth-node.depth))/(c2.depth-2*node.depth+c1.depth))
    elif method == 2:
        return(np.mean([c.pos for c in node.get_children()]))
    else: raise "Position method can only be 1 or 2."


def drawTree(tree, leafPos = None, depthDict = None, depthRangeDict=None, extendTips=False, rootIsZero=False,
             show=True, posMethod=1, col="black",linewidth=2,alpha=1,direction="down",taxColDict=None):
    tree = tree.copy("newick")
    
    #get node depths
    for node in tree.traverse():
        if depthDict is None: node.add_feature("depth", node.get_distance(tree))
        else: node.add_feature("depth", depthDict[node.name])
    
    #extend tips to align them if needed
    if extendTips:
        maxDP = max([l.depth for l in tree.iter_leaves()])
        for l in tree.iter_leaves(): l.depth = maxDP
    
    #adjust depths so that they are aligned at zero
    if not rootIsZero:
        maxDP = max([l.depth for l in tree.iter_leaves()])
        for node in tree.traverse():
            node.depth -= maxDP
    
    #set leaf positions
    if leafPos is None:
        leafNames = [l.name for l in tree.get_leaves()]
        leafPos = dict(zip(leafNames,range(len(leafNames))))
    
    for leaf in tree.iter_leaves(): leaf.add_feature("pos", leafPos[leaf.name])
    
    #set positions for all other nodes relative to their children
    for node in tree.traverse(strategy="postorder"):
        if not node.is_leaf():
            node.add_feature("pos", getNodePos(node,method=posMethod))
        
    #draw
    for node in tree.traverse():
        if direction is "down":
            plt.setp(plt.gca(),xticks=[])
            for child in node.get_children():
                plt.plot([node.pos,child.pos],[node.depth,child.depth],color=col,linewidth=linewidth,alpha=alpha, solid_capstyle="round")
                if depthRangeDict:
                    plt.plot([node.pos]*2,depthRangeDict[node.name],color=col,linewidth=1,alpha=alpha, solid_capstyle="round")
                    plt.plot([node.pos-.1,node.pos+.1],[node.depth]*2,color=col,linewidth=1,alpha=alpha, solid_capstyle="round")
            
            if node.is_leaf(): plt.text(node.pos, node.depth - 0.1, node.name,
                                        horizontalalignment='center', verticalalignment='center',
                                        color=taxColDict[node.name] if taxColDict else "black")
        
        else:
            plt.setp(plt.gca(),yticks=[])
            for child in node.get_children():
                plt.plot([node.depth,child.depth],[node.pos,child.pos],color=col,linewidth=linewidth,alpha=alpha, solid_capstyle="round")
                if depthRangeDict:
                    plt.plot(depthRangeDict[node.name],[node.pos]*2,color=col,linewidth=1,alpha=alpha, solid_capstyle="round")
                    plt.plot([node.depth]*2,[node.pos-.1,node.pos+.1],color=col,linewidth=1,alpha=alpha, solid_capstyle="round")
            
            if node.is_leaf(): plt.text(node.depth - 0.1, node.pos, node.name,
                                        horizontalalignment='left', verticalalignment='center',
                                        color=taxColDict[node.name] if taxColDict else "black")
    
    if show: plt.show()


def subset(things,subLen):
    starts = range(0,len(things),subLen)
    ends = [start+subLen for start in starts]
    return [things[starts[i]:ends[i]] for i in range(len(starts))]

def asColumn(a):
    return a[:,np.newaxis]

def addNodeNames(tree):
    for n in tree.traverse():
        if n.name is not None: n.name = "_".join(sorted(n.get_leaf_names()))

allIsNaN = lambda x: np.all(np.isnan(x))

def normLength(tree, outgroup):
    l = 1.*tree.get_distance(outgroup)
    for node in tree.traverse(): node.dist /= l


def addNodeNames(tree):
    for n in tree.traverse():
        if n.name is "": n.name = "_".join(sorted(n.get_leaf_names()))

def treeToParentChildTable(tree):
        return [(n.up.name,n.name,n.dist) for n in tree.traverse() if n.up is not None]


def getLeafPairs(node):
    assert len(node.children) == 2, "Node {} does not have two children.".format(node.name)
    return itertools.product(node.children[0].get_leaf_names(), node.children[1].get_leaf_names())

def make2DarrayFrom1DupperTriangle(upperTriangle1D,N,includesDiagnol=False):
    a = np.zeros([N,N])
    n=N if includesDiagnol else N-1
    indices = list(np.triu_indices(n))
    if not includesDiagnol: indices[1]+=1
    a[indices] = a[indices[::-1]] = upperTriangle1D
    return a

########################## plot topos with average branch lengths

parser = argparse.ArgumentParser()
parser.add_argument("-w", "--weightsFiles", help="Input weights file(s) from Twisst", action = "store", nargs= "+", required = True)
parser.add_argument("-d", "--distsFiles", help="Input dists file(s) from Twisst", action = "store", nargs = "+", required = True)

parser.add_argument("-f", "--figFile", help="File for output figure", action = "store", required = True)

parser.add_argument("--figFormat", help="Format of figFile", action = "store", default="pdf")

parser.add_argument("--figSize", help="Size of figFile", action = "store", nargs=2, type=float, default=(10,10,))

parser.add_argument("--posMethod", help="Node positioning method", choices=(1,2,), type=int, action = "store", default = 2)

parser.add_argument("--quantiles", help="Add quantiles for each node in tree", type=float, nargs = 2, action = "store")

parser.add_argument("--plotTaxa", help="Prune tree to include on the specifed taxa", nargs = "+", action = "store")

parser.add_argument("--taxOrder", choices = ("levelorder", "preorder", "postorder", "predefined"), action = "store", default="levelorder",
                    help="How to determine order of taxa in plots")

parser.add_argument("--lineWidth", help="Width for tree lines", type=float, action = "store", default= 4)

parser.add_argument("--scaleLinesByWeights", help="Scale tree lines in figure by weights", action = "store_true")

parser.add_argument("--orderByWeights", help="Order tree plots by weights", action = "store_true")

parser.add_argument("--cols", help="Topology colours", nargs = "+", action = "store")

parser.add_argument("--taxCols", help="Taxon name colours", nargs = "+", action = "store")

parser.add_argument("--alpha", help="Topology alpha", type=float, action = "store", default=1.)

parser.add_argument("--layout", help="Rows and columns to plot", nargs=2, type=int, action = "store")

parser.add_argument("--tight", help="Pading for tight edges", nargs=2, type=float, action = "store")


args = parser.parse_args()

sys.stderr.write("\nReading distances file...")
dists = np.vstack([np.loadtxt(f, skiprows=1) for f in args.distsFiles])

#get topologies
sys.stderr.write("\nGetting topologies...")
topos = []
with gzip.open(args.weightsFiles[0], "r") as wf:
    while True:
        try: topos.append(ete3.Tree(wf.readline().split()[-1]))
        except: break

nTopos = len(topos)

for t in topos: addNodeNames(t)

nTaxa = len(topos[0].get_leaves())

sys.stderr.write("\nThere are {} topologies and {} taxa".format(nTopos,nTaxa))

#make a separate set of topologies for plotting
plotTopos = [t.copy("newick") for t in topos]

if args.plotTaxa:
    for t in plotTopos: t.prune(args.plotTaxa)

if args.taxOrder == "predefined":
    assert args.plotTaxa is not None, "Predefined taxa order must be given using --plotTaxa."
    taxOrder = [args.plotTaxa]*nTopos
else:
    taxOrder = [[node.name for node in topo.traverse(strategy=args.taxOrder) if node.is_leaf()] for topo in plotTopos]

if args.taxCols:
    try: taxColDict = dict(zip(args.plotTaxa,args.taxCols))
    except: raise ValueError("To plot coloured taxon labels, you must specify names of taxa using --plotTaxa")
else:
    taxColDict = None

if args.layout: nRow,nCol = args.layout
elif nTopos == 3: nRow,nCol = (1,3,)
elif nTopos == 15: nRow,nCol = (3,5,)
elif nTopos == 105: nRow,nCol = (3,5,)
else: raise ValueError("Please specify number of rows and columns in plot using --layout")


#get pair names in the dists file. The order is essential here.
#We use the order of the first N headers, but assume that the rest follow the same pattern.
#for example, if the taxaare called A, B C and D, the headers should be:
#Topo1_A_B  Topo1_A_C  Topo1_A_D  Topo1_B_C  Topo1_B_D  Topo1_C_D  Topo2_A_B  Topo2_A_C ... ect
with gzip.open(args.distsFiles[0], "r") as df: pairs = df.readline().split()

pairs = [pair.split("_")[1:] for pair in pairs][:nTopos]
taxonNames = pairs[0] + [pair[1] for pair in pairs[1:nTaxa]]


#get columns for dists separated by topology
topo_column_indices = subset(range(dists.shape[1]), dists.shape[1]/nTopos)


#split topologies into a third dimension
dists = np.dstack([dists[:,i] for i in topo_column_indices])
dists = np.swapaxes(dists,1,2)


#set all missing dists to zero. Necessary for averaging, and doesnt impact results, because weighting for these is zero
dists = np.nan_to_num(dists)

############# weights

#read weights and convert to proportions
sys.stderr.write("\nReading weights file...")
weights = np.vstack([np.loadtxt(f, skiprows=nTopos+1) for f in args.weightsFiles])

assert weights.shape[0] == dists.shape[0]
assert weights.shape[1] == dists.shape[1]

#convert to proportions
rowSums = np.apply_along_axis(np.sum, 1, weights)
weights = weights / np.reshape(rowSums,[len(rowSums),1])

#convert any nan to zero
weights = np.nan_to_num(weights)

#get means
meanWeights = np.apply_along_axis(np.mean, 0, weights)


# now we need to get the average distance between leaves for each node in each topo
# the first step here is to get the two sets of leaves that descend from each node

nodes_all = [list(tree.traverse()) for tree in topos]

nodeNames = [[n.name for n in nodes] for nodes in nodes_all]

nodeLeafPairs = [[zip(*[(taxonNames.index(x),taxonNames.index(y),) for x,y in getLeafPairs(node)]) if not node.is_leaf() else None for node in nodes] for nodes in nodes_all]

#make a detpths array that gives the depth of each node for each topo at each window
depths = np.zeros([dists.shape[0], dists.shape[1], len(nodeNames[0])])

#now we go line by line, topology by topology and retrieve the depth
# as the avergae pairwise distance between all leaf pairs for each node
# unless the node is a leaf, in which case depth is zero.
sys.stderr.write("\nComputing depth for each node for each topology for each line in input...")
for x in range(depths.shape[0]):
    for y in range(nTopos):
        distMat = make2DarrayFrom1DupperTriangle(dists[x,y,:], nTaxa)
        for z in range(len(nodeNames[y])):
            depths[x,y,z] = distMat[nodeLeafPairs[y][z]].mean() if not nodes_all[y][z].is_leaf() else 0.0


#scale depths by dividing by the depth of the root
#the first node in each topo is the root, as the traversal goes to the root first
depths = depths / np.repeat(depths[:,:,0,np.newaxis], depths.shape[2], axis=2)
#anyehere we have nan is where the root depth was zero. This happens where we had missing data. So we can set all these tree depths to zero.
depths = np.nan_to_num(depths)


depths_average = np.average(depths, axis = 0, weights=np.repeat(weights[:,:,np.newaxis], depths.shape[2], axis=2))

depths_median = [[wquantiles.median(depths[:,j,k], weights=weights[:,j]) for k in range(depths.shape[2])] for j in range(depths.shape[1])]

if args.quantiles:
    depths_qL = [[wquantiles.quantile(depths[:,j,k], weights[:,j], args.quantiles[0]) for k in range(depths.shape[2])] for j in range(depths.shape[1])]
    depths_qU = [[wquantiles.quantile(depths[:,j,k], weights[:,j], args.quantiles[1]) for k in range(depths.shape[2])] for j in range(depths.shape[1])]


#cols = np.array([
    #"#2BCE48", #Green
    #"#005C31", #Forest
    #"#94FFB5", #Jade
    #"#9DCC00", #Lime
    #"#426600", #Quagmire
    #"#00998F", #Turquoise
    #"#5EF1F2", #Sky
    #"#0075DC", #Blue
    #"#003380", #Navy
    #"#740AFF", #Violet
    #"#FF5005", #Zinnia
    #"#F0A3FF", #Amethyst
    #"#FFA405", #Orpiment
    #"#FF0010", #Red
    #"#C20088"]) #Mallow

cols = args.cols if args.cols else ["#000000"]*nTopos

lineWidths = np.array([args.lineWidth]*nTopos, dtype=float)

if args.scaleLinesByWeights: lineWidths *= (meanWeights/meanWeights.max())

plotOrder = np.argsort(meanWeights)[::-1] if args.orderByWeights else range(nTopos)

sys.stderr.write("\nMaking plot.")

plt.figure(figsize=args.figSize, frameon=False)

for i in range(len(plotOrder)):
    plt.subplot(nRow,nCol,i+1)
    x = plotOrder[i]
    for y in np.arange(0,1.1,0.1): plt.plot([0,len(taxOrder[x])+1],[y,y],color="#CCCCCC")
    drawTree(plotTopos[x], leafPos = dict(zip(taxOrder[x],np.arange(1,len(taxOrder[x])+1))),
             depthDict = dict(zip(nodeNames[x],depths_median[x])),
             depthRangeDict = dict(zip(nodeNames[x],zip(depths_qL[x],depths_qU[x]))) if args.quantiles else None,
             show=False, alpha = args.alpha, posMethod = args.posMethod,
             linewidth = lineWidths[x], col=cols[x], taxColDict=taxColDict)
    axes = plt.gca()
    axes.set_ylim([-0.1,1.1])
    axes.set_xlim([0.5,len(taxOrder[x])+.5])
    axes.spines["top"].set_visible(False)
    axes.spines["right"].set_visible(False)
    axes.spines["bottom"].set_visible(False)
    axes.spines["left"].set_visible(False)
    plt.text(1,.95,"T"+str(x+1),color=cols[x],
             horizontalalignment='left', verticalalignment='center', size=12)

#plt.show()

if args.tight: plt.tight_layout(h_pad=args.tight[0], w_pad = args.tight[1])

plt.savefig(args.figFile, format=args.figFormat, figsize=args.figSize, frameon=False)
plt.close()

sys.stderr.write("\nDone.\n")
