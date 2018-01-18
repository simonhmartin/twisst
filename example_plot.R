
source("plot_twisst.R")

########### input data ################

#weights file with a column for each topology
weights_file <- "examples/msms_4of10_l1Mb_r10k.seq_gen.SNP.w50sites.phyml_bionj.weights.tsv"

#coordinates file for each window
window_data_file <- "examples/msms_4of10_l1Mb_r10k.seq_gen.SNP.w50sites.phyml_bionj.data.tsv"


########## read data ##################
weights = read.table(weights_file, header = T)
#normalise rows so weights sum to 1
weights <- weights / apply(weights, 1, sum)
#retrieve the names of the topologies
topoNames = names(weights)

window_data = read.table(window_data_file, header = T)

#exclude any rows where data is missing
good_rows = which(is.na(apply(weights,1,sum)) == F)
weights <- weights[good_rows,]
window_data = window_data[good_rows,]


########### choose colours for plot ########
#some nice contrasted colours. These are from https://en.wikipedia.org/wiki/Help:Distinguishable_colors
cols = c(
"#F0A3FF", #Amethyst
"#0075DC", #Blue
"#993F00", #Caramel
"#4C005C", #Damson
"#191919", #Ebony
"#005C31", #Forest
"#2BCE48", #Green
"#FFCC99", #Honeydew
"#808080", #Iron
"#94FFB5", #Jade
"#8F7C00", #Khaki
"#9DCC00", #Lime
"#C20088", #Mallow
"#003380", #Navy
"#FFA405", #Orpiment
"#FFA8BB", #Pink
"#426600", #Quagmire
"#FF0010", #Red
"#5EF1F2", #Sky
"#00998F", #Turquoise
"#E0FF66", #Uranium
"#740AFF", #Violet
"#990000", #Wine
"#FFFF80", #Xanthin
"#FFFF00", #Yellow
"#FF5005" #Zinnia
)

#semi-transparent version of each colour
trans_cols = paste0(cols, "25")

######### plot raw data #######

#plot raw data in "stepped" style, with polygons stacked.
#specify stepped style by providing a matrix of starts and ends for positions
pdf(file = paste0("examples/example.raw.stepped.stacked.pdf"), width = 10, height = 4)
par(mar = c(4,4,1,1))
plot.weights(weights_dataframe=weights, positions=cbind(window_data$start,window_data$end),
             line_cols=cols, fill_cols=cols, xlim =c(1,1000000),stacked=TRUE)
dev.off()

#plot raw data in stepped style, with polygons unstacked (stacked =FLASE)
#use semi-transparent colours for fill
pdf(file = paste0("examples/example.raw.stepped.unstacked.pdf"), width = 10, height = 4)
par(mar = c(4,4,1,1))
plot.weights(weights_dataframe=weights, positions=cbind(window_data$start,window_data$end),
             line_cols=cols, fill_cols=trans_cols, xlim =c(1,1000000),stacked=FALSE)
dev.off()

#plot raw data against window midpoints, with polygons unstacked.
#specify midpoint style by providing a single vector of midpoints for positions
pdf(file = paste0("examples/example.raw.midpoints.unstacked.pdf"), width = 10, height = 4)
par(mar = c(4,4,1,1))
plot.weights(weights_dataframe=weights, positions=window_data$mid,
             line_cols=cols, fill_cols=trans_cols, xlim =c(1,1000000),stacked=FALSE)
dev.off()

################ plot smoothed data ######
#use loess to smooth weights.
weights_smooth = smooth.weights(window_positions=window_data$mid, weights_dataframe = weights,
                                 span = 0.01, window_sites=window_data$sites)

#plot smoothed data with polygons stacked
pdf(file = paste0("examples/example.smooth.stacked.pdf"), width = 10, height = 4)
par(mar = c(4,4,1,1))
plot.weights(weights_dataframe=weights_smooth, positions=window_data$mid,
             line_cols=cols, fill_cols=cols, xlim =c(1,1000000),stacked=TRUE)
dev.off()

#plot smoothed data with polygons unstacked
pdf(file = paste0("examples/example.smooth.unstacked.pdf"), width = 10, height = 4)
par(mar = c(4,4,1,1))
plot.weights(weights_dataframe=weights_smooth, positions=window_data$mid,
             line_cols=cols, fill_cols=trans_cols, xlim =c(1,1000000),stacked=FALSE)
dev.off()



################ plot topologies #########

library(ape)

topos = read.tree(file="examples/4.topos")

#unrooted trees

pdf(file = paste0("examples/example.topos.unrooted.pdf"), width = 5, height = 2)
par(mfrow = c(1,3), mar = c(1,1,2,1), xpd=NA)
for (n in 1:length(topos)){
  plot.phylo(topos[[n]], type = "unrooted", edge.color=cols[n], edge.width=5, cex = 1, rotate.tree=90, adj = .5, label.offset=.2)
  mtext(side=3,text=paste0("topo",n))
  }
dev.off()


#root trees
for (i in 1:length(topos)) topos[[i]] <- root(topos[[i]], "D", resolve.root = T)

pdf(file = paste0("examples/example.topos.rooted.pdf"), width = 5, height = 2)
par(mfrow = c(1,3), mar = c(1,1,2,1), xpd=NA)
for (n in 1:length(topos)){
  plot.phylo(topos[[n]], type = "clad", edge.color=cols[n], edge.width=5, label.offset=.1, cex = 1)
  mtext(side=3,text=paste0("topo",n))
  }
dev.off()

