
import argparse
import itertools
import sys
import operator
import ete2, random, numpy
from ete2 import Tree


parser = argparse.ArgumentParser()
parser.add_argument("-n", "--nIts", help="Number of iterations for sub-sampling tree", type=int, action = "store", required = True, metavar="interger")
parser.add_argument("-t", "--treeFile", help="File containing tree(s) to analyse", required = True)
parser.add_argument("-o", "--topoFile", help="Output file of all topologies", required = True)
parser.add_argument("-w", "--weightsFile", help="Output file of all weights", required = True)
parser.add_argument("-g", "--group", help="Group name and individual names (separated by commas)", action='append', nargs=2, required = True, metavar=("name","inds"))


args = parser.parse_args()
#args = parser.parse_args("-n 5 -t test.trees -g A a,b,c -g B d,e,f -g C g,h,i -g D j,k,l".split())

nIts = args.nIts

groups = args.group

if len(groups) < 4:
  print >> sys.stderr, "Please specify at least four groups."
  sys.exit()

treeFileName = args.treeFile

topoFileName = args.topoFile

weightsFileName = args.weightsFile

##############################################################################################################################

def prod(iterable):
  return reduce(operator.mul, iterable, 1)

def sample(things, n = None, replace = False):
  if n == None:
    n = len(things)
  if replace == False:
    return random.sample(things,n)
  else:
    return [random.choice(things) for i in range(n)]

def weightTree(tree, taxa, taxonNames, nIts, topos = None, unrooted = True):
  if nIts > prod([len(t) for t in taxa]):
    combos = itertools.product(*taxa)
    nCombos = prod([len(t) for t in taxa])
  else:
    combos = [tuple(sample(x,1)[0] for x in taxa) for n in range(nIts)]
    nCombos = len(combos)
  if not topos:
    topos = []
  counts = [0]*len(topos)
  for combo in combos:
    pruned = tree.copy()
    pruned.prune(combo)
    #generify pruned tree
    for leaf in pruned.iter_leaves():
      for x in range(len(taxa)):
        if leaf.name in taxa[x]:
          leaf.name = taxonNames[x]
    #check for matches with each topology
    matchFound = False
    for x in range(len(topos)):
      if pruned.robinson_foulds(topos[x], unrooted_trees = unrooted)[0] == 0:
        counts[x] += 1
        matchFound = True
        break
    if not matchFound:
      topos.append(pruned)
      counts.append(1)
  weights = [1.0*c/nCombos for c in counts]
  return zip(topos,counts,weights)

def listToNwk(t):
  t = str(t)
  t = t.replace("[","(")
  t = t.replace("]",")")
  t = t.replace("'","")
  t += ";"
  return(t)

def allTrees(branches, trees = []):
  try:
    assert 4 <= len(branches) <= 8
  except:
    print "Please specify between 4 and 8 unique taxon names."
    return None
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
        tree = Tree(listToNwk(new_branches))
        if len(trees) == 0 or min([tree.robinson_foulds(t, unrooted_trees = True)[0] for t in trees]) > 0:
          trees.append(tree)
      else:
        #print "Tree still unresolved, so re-calling function."
        trees = allTrees(new_branches, trees)
  return(trees)



#parse taxa

taxonNames = [g[0] for g in groups]

taxa = [g[1].split(",") for g in groups]

#get all topologies
topos = allTrees(taxonNames, [])

#now make the last taxon the outgroup for all - just for standardisation
for topo in topos:
  topo.set_outgroup(taxonNames[-1])

topoFile = open(topoFileName, "w")

topoFile.write("\n".join([t.write(format = 9) for t in topos]) + "\n")

topoFile.close()

# check iterations

if nIts > prod([len(t) for t in taxa]):
  print >> sys.stderr, "Warning: number of iterations,", nIts, ", greater than number of possible combinations,", prod([len(t) for t in taxa]), ".\n"


#taxonNames = ["west", "east","guiana","cyd","tim","silv"]


#taxa = [["ros.MK523_A","ros.MK523_B","ros.MK524_A","ros.MK524_B","ros.MK525_A","ros.MK525_B","ros.MK589_A","ros.MK589_B","ros.MK675_A","ros.MK675_B","ros.MK676_A","ros.MK676_B","ros.MK682_A","ros.MK682_B","ros.MK683_A","ros.MK683_B","ros.MK687_A","ros.MK687_B","ros.MK689_A","ros.MK689_B","ros.CJ531_A","ros.CJ531_B","ros.CJ533_A","ros.CJ533_B","ros.CJ546_A","ros.CJ546_B","ros.CJ2071_A","ros.CJ2071_B","melP.CJ18038_A","melP.CJ18038_B","MelP.CJ18097_A","MelP.CJ18097_B","vul.CS519_A","vul.CS519_B","vul.CS10_A","vul.CS10_B","vul.CS11_A","vul.CS11_B","cyth.CJ2856_A","cyth.CJ2856_B"],["moc.CS228_A","moc.CS228_B","moc.CS16_A","moc.CS16_B","moc.CS17_A","moc.CS17_B","ple.CJ9156_A","ple.CJ9156_B","mapl.CJ16042_A","mapl.CJ16042_B","mal.CJ17162_A","mal.CJ17162_B","mal.CS21_A","mal.CS21_B","mal.CS22_A","mal.CS22_B","mal.CS24_A","mal.CS24_B","ecu.CJ9117_A","ecu.CJ9117_B","ama.JM216_A","ama.JM216_B","ama.JM160_A","ama.JM160_B","ama.JM293_A","ama.JM293_B","ama.JM48_A","ama.JM48_B","agl.JM108_A","agl.JM108_B","agl.JM122_A","agl.JM122_B","agl.JM569_A","agl.JM569_B","agl.JM572_A","agl.JM572_B","aman.CS2228_A","aman.CS2228_B","aman.CS2221_A","aman.CS2221_B"],["melG.CJ9315_A","melG.CJ9315_B","melG.CJ9316_A","melG.CJ9316_B","melG.CJ9317_A","melG.CJ9317_B","melG.CJ13435_A","melG.CJ13435_B","thel.CJ13566_A","thel.CJ13566_B"],["cyd.CJ553_A","cyd.CJ553_B","cyd.CJ560_A","cyd.CJ560_B","cyd.CJ564_A","cyd.CJ564_B","cyd.CJ565_A","cyd.CJ565_B"],["tim.JM313_A","tim.JM313_B","tim.JM57_A","tim.JM57_B","tim.JM84_A","tim.JM84_B","tim.JM86_A","tim.JM86_B"],["hec.JM273_A","hec.JM273_B","eth.JM67_A","eth.JM67_B","par.JM371_A","par.JM371_B","ser.JM202_A","ser.JM202_B"]]


### file for weights

weightsFile = open(weightsFileName, "w")

weightsFile.write(",".join(["topo" + str(x) for x in range(len(topos))]) + "\n")


#open tree file
treeFile = open(treeFileName, "r")

line = treeFile.readline()

n = 0

while len(line) > 1:
  tree = Tree(line, format = 9)
  treeWeights = weightTree(tree, taxa, taxonNames, nIts, topos)
  weightsFile.write(",".join([str(x[2]) for x in treeWeights]) + "\n")
  print n
  n += 1
  line = treeFile.readline()
  #line = treeFile.readline() # add second one if there are line breaks between trees

treeFile.close()
weightsFile.close()


