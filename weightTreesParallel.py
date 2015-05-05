
import argparse
import itertools
import sys
import gzip
import operator
import ete2, random, numpy
from ete2 import Tree

from multiprocessing import Process, Queue
from multiprocessing.queues import SimpleQueue
from threading import Thread
from time import sleep


parser = argparse.ArgumentParser()
parser.add_argument("-n", "--nIts", help="Number of iterations for sub-sampling tree", type=int, action = "store", required = True, metavar="iterations")
parser.add_argument("-t", "--treeFile", help="File containing tree(s) to analyse", required = True)
parser.add_argument("-o", "--topoFile", help="Output file of all topologies", required = True)
parser.add_argument("-w", "--weightsFile", help="Output file of all weights", required = True)
parser.add_argument("-g", "--group", help="Group name and individual names (separated by commas)", action='append', nargs=2, required = True, metavar=("name","inds"))
parser.add_argument("-T", "--Threads", help="Number of worker threads for parallel processing", type=int, default=1, required = False, metavar="threads")
parser.add_argument("--verbose", help="Verbose output", action="store_true")


args = parser.parse_args()
#args = parser.parse_args("-n 5 -t test.trees -o test.topos.txt -w test.weights.B.csv -g A a,b,c -g B d,e,f -g C g,h,i -g D j,k,l".split())

nIts = args.nIts

groups = args.group

if len(groups) < 4:
  print >> sys.stderr, "Please specify at least four groups."
  sys.exit()

treeFileName = args.treeFile

topoFileName = args.topoFile

weightsFileName = args.weightsFile

threads = args.Threads

verbose = args.verbose

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

####################################################################################################################################


'''A function that reads from the line queue, calls some other function and writes to the results queue
This function needs to be tailored to the particular analysis funcion(s) you're using. This is the function that will run on each of the N cores.'''
def weightTree_wrapper(lineQueue, resultQueue, taxa, taxonNames, nIts, topos = None, unrooted = True):
  while True:
    lineNumber,line = lineQueue.get()
    if line[-1] == ";":
      tree = Tree(line, format = 9)
      treeWeights = weightTree(tree, taxa, taxonNames, nIts, topos, unrooted)
      result = ",".join([str(x[2]) for x in treeWeights])
    else:
      result = ",".join(["NA"]*len(topos))
    resultQueue.put((lineNumber, result, True))


'''a function that watches the result queue and sorts results. This should be a generic funcion regardless of the result, as long as the first object is the line number, and this increases consecutively.'''
def sorter(resultQueue, writeQueue, verbose):
  global resultsReceived
  sortBuffer = {}
  expect = 0
  while True:
    lineNumber,result,good = resultQueue.get()
    resultsReceived += 1
    if verbose:
      print >> sys.stderr, "Sorter received result", lineNumber
    if lineNumber == expect:
      writeQueue.put((lineNumber,result,good))
      if verbose:
        print >> sys.stderr, "Result", lineNumber, "sent to writer"
      expect +=1
      #now check buffer for further results
      while True:
        try:
          result,good = sortBuffer.pop(str(expect))
          writeQueue.put((expect,result,good))
          if verbose:
            print >> sys.stderr, "Result", expect, "sent to writer"
          expect +=1
        except:
          break
    else:
      #otherwise this line is ahead of us, so add to buffer dictionary
      sortBuffer[str(lineNumber)] = (result,good)

'''a writer function that writes the sorted result. This is also generic'''
def writer(writeQueue, out):
  global resultsWritten
  global resultsHandled
  while True:
    lineNumber,result,good = writeQueue.get()
    if verbose:
      print >> sys.stderr, "Writer received result", lineNumber
      if good:
        print >> sys.stderr, "Writing good result."
      else:
        print >> sys.stderr, "Omitting bad result."
    if good:
      out.write(result + "\n")
      resultsWritten += 1
    resultsHandled += 1


'''loop that checks line stats'''
def checkStats():
  while True:
    sleep(10)
    print >> sys.stderr, linesQueued, "trees queued", resultsReceived, "results received", resultsWritten, "results written."

############################################################################################################################################


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





### file for weights

if weightsFileName[-3:] == ".gz":
  weightsFile = gzip.open(weightsFileName, "w")
else:
  weightsFile = open(weightsFileName, "w")

weightsFile.write(",".join(["topo" + str(x) for x in range(len(topos))]) + "\n")



############################################################################################################################################

#counting stat that will let keep track of how far we are
linesQueued = 0
resultsReceived = 0
resultsWritten = 0
resultsHandled = 0

'''Create queues to hold the data one will hold the line info to be passed to the analysis'''
lineQueue = SimpleQueue()
#one will hold the results (in the order they come)
resultQueue = SimpleQueue()
#one will hold the sorted results to be written
writeQueue = SimpleQueue()


'''start worker Processes for analysis. The comand should be tailored for the analysis wrapper function
of course these will only start doing anything after we put data into the line queue
the function we call is actually a wrapper for another function.(s) This one reads from the line queue, passes to some analysis function(s), gets the results and sends to the result queue'''
for x in range(threads):
  worker = Process(target=weightTree_wrapper, args = (lineQueue, resultQueue, taxa, taxonNames, nIts, topos,))
  worker.daemon = True
  worker.start()


'''thread for sorting results'''
worker = Thread(target=sorter, args=(resultQueue,writeQueue,verbose,))
worker.daemon = True
worker.start()

'''start thread for writing the results'''
worker = Thread(target=writer, args=(writeQueue, weightsFile,))
worker.daemon = True
worker.start()


'''start background Thread that will run a loop to check run statistics and print
We use thread, because I think this is necessary for a process that watches global variables like linesTested'''
worker = Thread(target=checkStats)
worker.daemon = True
worker.start()



############################################################################################################################################

#open tree file

if treeFileName[-3:] == ".gz":
  treeFile = gzip.open(treeFileName, "r")
else:
  treeFile = open(treeFileName, "r")

line = treeFile.readline().rstrip()


while len(line) > 1:
  lineQueue.put((linesQueued,line))
  linesQueued += 1
  line = treeFile.readline().rstrip()
  #line = treeFile.readline() # add second one if there are line breaks between trees

############################################################################################################################################


### wait for queues to empty

#print >> sys.stderr, "\nWaiting for line queue to empty..."
#while lineQueue.empty() == False:
  #sleep(1)

#print >> sys.stderr, "\nWaiting for result queue to empty..."
#while resultQueue.empty() == False:
  #sleep(1)

#print >> sys.stderr, "\nWaiting for write queue to empty...\n"
#while writeQueue.empty() == False:
  #sleep(1)

print >> sys.stderr, "\nWriting final results...\n"
while resultsHandled < linesQueued:
  sleep(1)

sleep(5)

treeFile.close()
weightsFile.close()


print >> sys.stderr, str(linesQueued), "lines were read.\n"
print >> sys.stderr, str(resultsWritten), "results were written.\n"

print "\nDone."

sys.exit()





