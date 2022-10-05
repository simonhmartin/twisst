import numpy as np
import msprime
import demes
import demesdraw
import twisst
import gzip
#import matplotlib.pyplot as plt

########################### Introduction ##################################

# This script will simulate genealogies under a defined population history using msprime
# and then analyse these using Twisst to determine expected topology weights under
# a known scenario. msprime is a stochastic simulator, so if we want to know the expected
# proportions under neutrality, we need to ensure that we simulate enough genealogies.
# This can be achieved by simulating one long chromosome (which is slow)
# or many short blocks of genome (which is faster).
# We send the output weights to a zipped text file.

########################### demographic model ##################################

# Step 1 is to define the demographic history. The model can be defined using msprime
# notation directly, but we instead use the Demes pachage, which as the advantage that
# we can (1) import any model we want from a text file, and (2) generate a nice figure
# to ensure we have specified the model correctly.

#importing from a txt file would look like this:
#demes_graph = demes.load("model_file.yaml")

# Alternatively, we can define a function that builds the model using demes.Builder.
# This has the advantage that we can edit parameters of the model within this script.

# In this example, the species tree has shape (((V,W),(X,Y)),Z).
# We allow all populations and all ancestral populations to have different sizes
# (although all defaults are set to teh same size). We also allow the split times to vary.
# Finally, we have added parameters for migration between two of the populations (in both directions).
# But the migration rates are set to 0 by default.
# You can specify any model you like, or edit this one.

def four_pop_model(t_VW=2e4, t_XY=4e4, t_VWXY=1e5, t_VWXYZ=5e5,    #split times
                    N_V=1e5, N_W=1e5, N_X=1e5, N_Y=1e5, N_Z=1e5,   #population sizes
                    N_VW=1e5, N_XY=1e5, N_VWXY=1e5, N_VWXYZ=1e5,   #ancestral population sizes
                    m_X_to_W=0, m_W_to_X=0):                       #migration rates
    
    #use the deme builder to set up the demographic history
    b = demes.Builder(time_units="generations")
    b.add_deme("VWXYZanc",                      epochs=[{"start_size":N_VWXYZ, "end_time":t_VWXYZ}])
    b.add_deme("VWXYanc",ancestors=["VWXYZanc"], epochs=[{"start_size":N_VWXY, "end_time":t_VWXY}])
    b.add_deme("VWanc", ancestors=["VWXYanc"], epochs=[{"start_size":N_VW,  "end_time":t_VW}])
    b.add_deme("XYanc", ancestors=["VWXYanc"], epochs=[{"start_size":N_XY,  "end_time":t_XY}])
    b.add_deme("V",     ancestors=["VWanc"],  epochs=[{"start_size":N_V}])
    b.add_deme("W",     ancestors=["VWanc"],  epochs=[{"start_size":N_W}])
    b.add_deme("X",     ancestors=["XYanc"],  epochs=[{"start_size":N_X}])
    b.add_deme("Y",     ancestors=["XYanc"],  epochs=[{"start_size":N_Y}])
    b.add_deme("Z",     ancestors=["VWXYZanc"], epochs=[{"start_size":N_Z}])
    
    #add migration to the model, but only if the rate specified is > 0
    if m_W_to_X > 0: b.add_migration(source="W", dest="X", rate=m_W_to_X)
    if m_X_to_W > 0: b.add_migration(source="X", dest="W", rate=m_X_to_W)
    
    graph = b.resolve()
    
    return(graph)


######################## Simulate a single chromosome  ##########################
#This section will simulate a single chromosome with recombination.

#First define the model. Here we have not added any parameters because we will use the defaults for all parameters
demes_graph = four_pop_model()

# Export a figure to check the model looks right
fig=demesdraw.tubes(demes_graph).get_figure()
fig.savefig("model_graph.pdf", bbox_inches='tight')


ts = msprime.sim_ancestry(samples={"V":4, "W":4, "X":4, "Y":4, "Z":4},
                          demography=msprime.Demography.from_demes(demes_graph),
                          sequence_length = 1e5,
                          recombination_rate = 5e-8, ploidy=2)

results = twisst.weightTrees(ts, treeFormat="ts",
                   taxonNames=["4","5","6","7","8"],
                               outgroup="8", verbose=False)

#write to output files
with gzip.open("sim_1chrom_weights.tsv.gz", "wt") as weights_file:
    twisst.writeWeights(weights_file, results)


############################ simulate blocks #########################

n_blocks=500

# run simulations to produce tree sequence objects
ts_blocks =  [msprime.sim_ancestry(samples={"V":4, "W":4, "X":4, "Y":4, "Z":4},
                                   demography=msprime.Demography.from_demes(demes_graph),
                                   sequence_length = 1e4,
                                   recombination_rate = 1e-8, ploidy=2) for i in range(n_blocks)]


#get topology weights from Twisst
results_blocks = [None]*n_blocks

for i in range(n_blocks):
    print("block", i)
    results_blocks[i] = twisst.weightTrees(ts_blocks[i], treeFormat="ts",
                                           taxonNames=["4","5","6","7","8"],
                                           outgroup="8", verbose=False)

#view quick summary
#twisst.summary(results_blocks[0])

#write to output files
with gzip.open("sim_500blocks_weights.tsv.gz", "wt") as weights_file:
    for i in range(n_blocks):
        twisst.writeWeights(weights_file, results_blocks[i], include_topologies=True if i==0 else False, include_header=True if i==0 else False)
