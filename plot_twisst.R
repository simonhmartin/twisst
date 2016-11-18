
simple.loess.predict <- function(x, y, span, weights = NULL, max = NULL, min = NULL){
    y.loess <- loess(y ~ x, span = span, weights = weights)
    y.predict <- predict(y.loess,x)
    if (is.null(min) == FALSE) {y.predict = ifelse(y.predict > min, y.predict, min)}
    if (is.null(max) == FALSE) {y.predict = ifelse(y.predict < max, y.predict, max)}
    y.predict
    }

smooth_df <- function(x, df, span, col.names=NULL, weights=NULL, min=NULL, max=NULL){
    smoothed <- df
    if (is.null(col.names)){col.names=colnames(df)}
    for (col.name in col.names){
        print(paste("smoothing",col.name))
        smoothed[,col.name] <- simple.loess.predict(x,df[,col.name],span = span, max = max, min = min, weights = weights)
        }
    smoothed
    }


stack <- function(mat){
    upper <- t(apply(mat, 1, cumsum))
    lower <- upper - mat
    list(upper=upper,lower=lower)
    }

interleave <- function(x1,x2){
    output <- vector(length= length(x1) + length(x2))
    output[seq(1,length(output),2)] <- x1
    output[seq(2,length(output),2)] <- x2
    output
    }


plot_weights <- function(weights_dataframe,positions=NULL,line_cols=NULL,fill_cols=NULL,xlim=NULL,ylim=c(0,1),stacked=FALSE,
                                        ylab="Weights", xlab = "Position", main="",xaxt=NULL,yaxt=NULL,bty="n"){
    #get x axis
    x = positions
    #if a two-column matrix is given - plot step-like weights with start and end of each window    
    if (is.matrix(x)==TRUE) {
        x = interleave(positions[,1],positions[,2])
        yreps=2
        }
    else {
        if (is.null(x)==FALSE) x = positions
        else x = 1:nrow(weights_dataframe)
        yreps=1
        }
    
    #set x limits
    if(is.null(xlim)) xlim = c(min(x), max(x))
    
    plot(0, pch = "", xlim = xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main,xaxt=xaxt,yaxt=yaxt,bty=bty)
    
    if (stacked == TRUE){
        stacked <- stack(weights_dataframe)
        for (n in 1:ncol(weights_dataframe)){
            y_upper = rep(stacked[["upper"]][,n],each=yreps)
            y_lower = rep(stacked[["lower"]][,n],each = yreps)
            polygon(c(x,rev(x)),c(y_upper, rev(y_lower)), col = fill_cols[n], border=NA)
            }
        }
    else{
        for (n in 1:ncol(weights_dataframe)){
            y = rep(weights_dataframe[,n],each=yreps)
            polygon(c(x,rev(x)),c(y, rep(0,length(y))), col=fill_cols[n], border=NA)
            lines(x,y, type = "l", col = line_cols[n])
            }
        }
    }

options(scipen = 5)

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
#some nice contrasted colours.
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
plot_weights(weights_dataframe=weights, positions=cbind(window_data$start,window_data$end),
             line_cols=cols, fill_cols=cols, xlim =c(1,100000),stacked=TRUE)
dev.off()

#plot raw data in stepped style, with polygons unstacked (stacked =FLASE)
#use semi-transparent colours for fill
pdf(file = paste0("examples/example.raw.stepped.unstacked.pdf"), width = 10, height = 4)
par(mar = c(4,4,1,1))
plot_weights(weights_dataframe=weights, positions=cbind(window_data$start,window_data$end),
             line_cols=cols, fill_cols=trans_cols, xlim =c(1,100000),stacked=FALSE)
dev.off()

#plot raw data against window midpoints, with polygons unstacked.
#specify midpoint style by providing a single vector of midpoints for positions
pdf(file = paste0("examples/example.raw.midpoints.unstacked.pdf"), width = 10, height = 4)
par(mar = c(4,4,1,1))
plot_weights(weights_dataframe=weights, positions=window_data$mid,
             line_cols=cols, fill_cols=trans_cols, xlim =c(1,100000),stacked=FALSE)
dev.off()

################ plot smoothed data ######
#use loess to smooth weights.
span = 0.01
weights_smooth <- smooth_df(x=window_data$mid,weights,col.names=topoNames,span=span, min=0, max=1,weights=window_data$sites)
#rescale to sum to 1
weights_smooth <- weights_smooth / apply(weights_smooth, 1, sum)



#plot smoothed data with polygons stacked
pdf(file = paste0("examples/example.smooth.stacked.pdf"), width = 10, height = 4)
par(mar = c(4,4,1,1))
plot_weights(weights_dataframe=weights_smooth, positions=window_data$mid,
             line_cols=cols, fill_cols=cols, xlim =c(1,100000),stacked=TRUE)
dev.off()

#plot smoothed data with polygons unstacked
pdf(file = paste0("examples/example.smooth.unstacked.pdf"), width = 10, height = 4)
par(mar = c(4,4,1,1))
plot_weights(weights_dataframe=weights_smooth, positions=window_data$mid,
             line_cols=cols, fill_cols=trans_cols, xlim =c(1,100000),stacked=FALSE)
dev.off()


