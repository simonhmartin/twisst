
################################# overview #####################################

# The main data produced by Twisst is a weights file which has columns for each
# topology and their number of observations of that topology within each
# genealogy. Weights files produced by Twisst also contain initial comment lines
# speficying the topologies.

# The other data file that may be of interest is the window data. That is, the
# chromosome/scaffold and start and end positions for each of the regions or
# windows represented in the weights file.

# Both of the above files can be read into R, manipulated and plotted however
# you like, but I have written some functions to make these tasks easier.
# These functions are provided in the script plot_twisst.R

################### load helpful plotting functions #############################

source("plot_twisst.R")

############################## input files ######################################

# It is possible to import one or more weights files at a time.
# Here we just import 1

#weights file with a column for each topology
weights_file <- "examples/msms_4of10_l1Mb_r10k_sweep.seq_gen.SNP.w50sites.phyml_bionj.weights.tsv.gz"


# It is not necessary to import window data files, but if you do there should be one for
# each weights file

#coordinates file for each window
window_data_file <- "examples/msms_4of10_l1Mb_r10k_sweep.seq_gen.SNP.w50sites.phyml_bionj.data.tsv.gz"


################################# import data ##################################

# The function import.twisst reads the weights, window data  files into a list object
# If there are multiple weights files, or a single file with different chromosomes/scaffolds/contigs
# in the window data file, these will be separated when importing.

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)


############################## combined plots ##################################
# there are a functions available to plot both the weightings and the topologies

#a summary plot shows all the topologies and a bar plot of their relative weightings
plot.twisst.summary(twisst_data, lwd=3, cex=0.7)


#or plot ALL the data across the chromosome(s)
#Note, this is not recommended if there are large numbers of windows.
# instead, it is recommended to first smooth the weghtings and plot the smoothed values
plot.twisst(twisst_data, mode=1, show_topos=TRUE)


# make smooth weightings and plot those across chromosomes
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 20000, spacing = 1000)
plot.twisst(twisst_data_smooth, mode=2) #mode 2 overlays polygons, mode 3 would stack them


##################### individual plots: raw weights ############################

#plot raw data in "stepped" style, with polygons stacked.
#specify stepped style by providing a matrix of starts and ends for positions
par(mfrow = c(1,1), mar = c(4,4,1,1))
plot.weights(weights_dataframe=twisst_data$weights[[1]], positions=twisst_data$window_data[[1]][,c("start","end")],
             line_cols=topo_cols, fill_cols=topo_cols, stacked=TRUE)

#plot raw data in stepped style, with polygons unstacked (stacked =FLASE)
#use semi-transparent colours for fill
plot.weights(weights_dataframe=twisst_data$weights[[1]], positions=twisst_data$window_data[[1]][,c("start","end")],
             line_cols=topo_cols, fill_cols=paste0(topo_cols,80), stacked=FALSE)


#################### individual plots: smoothed weights ########################

#plot smoothed data with polygons stacked
plot.weights(weights_dataframe=twisst_data_smooth$weights[[1]], positions=twisst_data_smooth$pos[[1]],
             line_cols=topo_cols, fill_cols=topo_cols, stacked=TRUE)

#plot smoothed data with polygons unstacked
plot.weights(weights_dataframe=twisst_data_smooth$weights[[1]], positions=twisst_data_smooth$pos[[1]],
             line_cols=topo_cols, fill_cols=paste0(topo_cols,80), stacked=FALSE)



########################### plot topologies using Ape ##########################
#unrooted trees
for (i in 1:length(twisst_data$topos)) twisst_data$topos[[i]] <- ladderize(unroot(twisst_data$topos[[i]]))

par(mfrow = c(1,length(twisst_data$topos)), mar = c(1,1,2,1), xpd=NA)
for (n in 1:length(twisst_data$topos)){
  plot.phylo(twisst_data$topos[[n]], type = "unrooted", edge.color=topo_cols[n], edge.width=5, rotate.tree = 90, cex = 1, adj = .5, label.offset=.2)
  mtext(side=3,text=paste0("topo",n))
  }


#rooted topologies
for (i in 1:length(twisst_data$topos)) twisst_data$topos[[i]] <- root(twisst_data$topos[[i]], "D", resolve.root = T)

par(mfrow = c(1,length(twisst_data$topos)), mar = c(1,1,2,1), xpd=NA)
for (n in 1:length(twisst_data$topos)){
  plot.phylo(twisst_data$topos[[n]], type = "clad", edge.color=topo_cols[n], edge.width=5, label.offset=.1, cex = 1)
  mtext(side=3,text=paste0("topo",n))
  }

################## subset to only the most abundant topologies #################

#get list of the top two most abundant topologies
top2_topos <- order(twisst_data$weights_overall_mean, decreasing=T)[1:2]

#subset twisst object for these
twisst_data_top2topos <- subset.twisst.by.topos(twisst_data, top2_topos)
#this can then be used in all the same plotting functions above.

###################### subset to only specific regions #########################

#regions to keep (more than one can be specified)
regions <- c("contig0")

#subset twisst object for these
twisst_data_contig0 <- subset.twisst.by.regions(twisst_data, regions)
#this can then be used in all the same plotting functions above.
