simple.loess.predict <- function(x, y, span, new_x=NULL, weights = NULL, max = NULL, min = NULL, family=NULL){
    y.loess <- loess(y ~ x, span = span, weights = weights, family=family)
    if (is.null(new_x)) {y.predict <- predict(y.loess,x)}
    else {y.predict <- predict(y.loess,new_x)}
    if (is.null(min) == FALSE) {y.predict = ifelse(y.predict > min, y.predict, min)}
    if (is.null(max) == FALSE) {y.predict = ifelse(y.predict < max, y.predict, max)}
    y.predict
    }

smooth.df <- function(x, df, span, new_x = NULL, col.names=NULL, weights=NULL, min=NULL, max=NULL, family=NULL){
    if (is.null(new_x)) {smoothed <- df}
    else smoothed = df[1:length(new_x),]
    if (is.null(col.names)){col.names=colnames(df)}
    for (col.name in col.names){
        print(paste("smoothing",col.name))
        smoothed[,col.name] <- simple.loess.predict(x,df[,col.name],span = span, new_x = new_x, max = max, min = min, weights = weights, family=family)
        }
    smoothed
    }

smooth.weights <- function(window_positions, weights_dataframe, span, new_positions=NULL, window_sites=NULL){
    weights_smooth <- smooth.df(x=window_positions,df=weights_dataframe,
                                span=span, new_x=new_positions, min=0, max=1, weights=window_sites)

    #return rescaled to sum to 1
    weights_smooth / apply(weights_smooth, 1, sum)
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


sum_df_columns <- function(df, columns_list){
    new_df <- df[,0]
    for (x in 1:length(columns_list)){
        if (length(columns_list[[x]]) > 1) new_df[,x] <- apply(df[,columns_list[[x]]], 1, sum, na.rm=T)
        else new_df[,x] <- df[,columns_list[[x]]]
        if (is.null(names(columns_list)[x]) == FALSE) names(new_df)[x] <- names(columns_list)[x]
        }
    new_df
    }


plot.weights <- function(weights_dataframe,positions=NULL,line_cols=NULL,fill_cols=NULL,density=NULL,lwd=1,xlim=NULL,ylim=c(0,1),stacked=FALSE,
                                        ylab="Weighting", xlab = "Position", main="",xaxt=NULL,yaxt=NULL,bty="n", add=FALSE){
    #get x axis
    x = positions
    #if a two-column matrix is given - plot step-like weights with start and end of each window    
    if (dim(as.matrix(x))[2]==2) {
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
    
    #if not adding to an old plot, make a new plot
    if (add==FALSE) plot(0, pch = "", xlim = xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main,xaxt=xaxt,yaxt=yaxt,bty=bty)
    
    if (stacked == TRUE){
        y_stacked <- stack(weights_dataframe)
        for (n in 1:ncol(weights_dataframe)){
            y_upper = rep(y_stacked[["upper"]][,n],each=yreps)
            y_lower = rep(y_stacked[["lower"]][,n],each = yreps)
            polygon(c(x,rev(x)),c(y_upper, rev(y_lower)), col = fill_cols[n], density=density[n], border=NA)
            }
        }
    else{
        for (n in 1:ncol(weights_dataframe)){
            y = rep(weights_dataframe[,n],each=yreps)
            polygon(c(x,rev(x)),c(y, rep(0,length(y))), col=fill_cols[n], border=NA,density=density[n])
            lines(x,y, type = "l", col = line_cols[n],lwd=lwd)
            }
        }
    }

options(scipen = 7)

#Heres a set of 15 colourful colours from https://en.wikipedia.org/wiki/Help:Distinguishable_colors
topo_cols <- c(
"#0075DC", #Blue
"#2BCE48", #Green
"#FFA405", #Orpiment
"#5EF1F2", #Sky
"#FF5005", #Zinnia
"#005C31", #Forest
"#00998F", #Turquoise
"#FF0010", #Red
"#9DCC00", #Lime
"#003380", #Navy
"#F0A3FF", #Amethyst
"#740AFF", #Violet
"#426600", #Quagmire
"#C20088", #Mallow
"#94FFB5") #Jade


########### Below are some more object-oriented tools for working with standard twisst output files

library(ape)

library(tools)

#a function that imports weights 
import.twisst <- function(weights_files, window_data_files, cleanup=TRUE){
    l = list()

    l$n_datasets <- length(weights_files)
    
    l$weights_raw <- lapply(weights_files, read.table ,header=TRUE)
        
    l$weights <- sapply(l$weights_raw, function(raw) raw/apply(raw, 1, sum), simplify=F)
    
    #get window_data if present
    l$window_data <- lapply(window_data_files, read.table ,header=TRUE)
    
    if (cleanup==TRUE){
        for (i in 1:l$n_datasets){
            #remove rows containing NA values
            good_rows = which(is.na(apply(l$weights[[i]],1,sum)) == F)
            l$weights[[i]] <- l$weights[[i]][good_rows,]
            l$window_data[[i]] = l$window_data[[i]][good_rows,]
            }
        }
    
    l$pos = sapply(l$weights, function(w) 1:nrow(w))
    
        #attempt to retrieve topologies
    if ("package:ape" %in% search()){
        n_topos = ncol(l$weights[[1]])
        if (file_ext(weights_files[[1]]) == ".gz") cat="cat" else cat="zcat"
        topos_text = system(paste(cat, weights_files[[1]], "| head -n", n_topos), intern = T)
        l$topos <- read.tree(text = topos_text)
        }
    else l$topos = NULL
    
    l
    }


smooth.twisst <- function(twisst_object, span=0.05) {
    l=list()
    
    l$topos <- twisst_object$topos
    
    l$n_datasets <- twisst_object$n_datasets
    
    l$weights <- list()
    
    l$pos <- list()
    
    for (i in 1:l$n_datasets){
        if (is.null(twisst_object$window_data[[i]]$mid) == TRUE) {
            mid <- (twisst_object$window_data[[i]]$start + twisst_object$window_data[[i]]$end)/2
            }
        else mid = twisst_object$window_data[[i]]$mid
        
        l$pos[[i]] <- seq(mid[1], tail(mid,1), tail(mid,1)*span*.1)
        
        l$weights[[i]] <- smooth.weights(mid, twisst_object$weights[[i]], new_x <- l$pos[[i]], span = span,
                                         window_sites=twisst_object$window_data$sites[[i]])
        }
    l
    }


plot.twisst <- function(twisst_object, show_topos=TRUE, ncol_topos=NULL, show_weights=TRUE, datasets=NULL, ncol_weights=1,
                        cols=topo_cols,xlim=NULL,stacked=TRUE, rel_height=3, tree_type="clad"){
    
    if (is.null(datasets)==TRUE) datasets <- 1:twisst_object$n_datasets
    
    n_topos <- length(twisst_object$topos)
    
    if (is.null(ncol_topos)) ncol_topos <- n_topos
    
    #if we have too few topologies to fill the spaces in the plot, we can pad in the remainder
    topos_pad <- (n_topos * ncol_weights) %% (ncol_topos*ncol_weights) 
    
    topos_layout_matrix <- matrix(c(rep(1:n_topos, each=ncol_weights), rep(0, topos_pad)),
                                  ncol=ncol_topos*ncol_weights, byrow=T)
    
    #if we have too few datasets to fill the spaces in the plot, we pad in the remainder
    data_pad <- (length(datasets)*ncol_topos) %% (ncol_topos*ncol_weights)
    
    weights_layout_matrix <- matrix(c(rep(n_topos+(1:length(datasets)), each=ncol_topos),rep(0,data_pad)),
                                    ncol=ncol_topos*ncol_weights, byrow=T)
    
    layout(rbind(topos_layout_matrix, weights_layout_matrix),
           height=c(rep(1, nrow(topos_layout_matrix)), rep(rel_height, nrow(weights_layout_matrix))))
    
    par(mar=c(1,1,1,1))
    
    for (i in 1:n_topos){
        plot.phylo(twisst_object$topos[[i]], type = tree_type, edge.color=cols[i],
                   edge.width=5, label.offset=.4, cex=1)
        mtext(side=3,text=paste0("topo",i), cex=0.75)
        }
    
    par(mar=c(4,4,2,2))
    
    for (j in datasets){
        if (is.null(twisst_object$window_data[[j]])) positions <- twisst_object$pos[[j]]
        else positions <- twisst_object$window_data[[j]][,c("start","end")]
        plot.weights(twisst_object$weights[[j]], positions, fill_cols = cols, line_cols=NA,lwd=0, stacked=T)
        }
    }


#function for plotting tree that uses ape to get node positions
draw.tree <- function(phy, x, y, x_scale=1, y_scale=1, method=1, direction="right",
                      edge_col="black", label_col="black", add_labels=TRUE, add_symbols=FALSE,
                      label_offset = 1, symbol_offset=0, symbol_col="black",symbol_bg="NA",pch=19){

    n_tips = length(phy$tip.label)

    if (direction=="right") {
        node_x = (node.depth(phy, method=method) - 1) * x_scale * -1
        node_y = node.height(phy) * y_scale
        label_x = node_x[1:n_tips] + label_offset
        label_y = node_y[1:n_tips]
        adj_x = 0
        adj_y = .5
        symbol_x = node_x[1:n_tips] + symbol_offset
        symbol_y = node_y[1:n_tips]
        }
    if (direction=="down") {
        node_y = (node.depth(phy, method=method) - 1) * y_scale * 1
        node_x = node.height(phy) * x_scale
        label_x = node_x[1:n_tips]
        label_y = node_y[1:n_tips] - label_offset
        adj_x = .5
        adj_y = 1
        symbol_x = node_x[1:n_tips]
        symbol_y = node_y[1:n_tips] - symbol_offset
        }
    
    #draw edges
    segments(x + node_x[phy$edge[,1]], y + node_y[phy$edge[,1]],
             x + node_x[phy$edge[,2]], y + node_y[phy$edge[,2]], col=edge_col)
    
    if (add_labels=="TRUE") text(x + label_x, y + label_y, col = label_col, labels=phy$tip.label, adj=c(adj_x,adj_y))
    if (add_symbols=="TRUE") points(x + symbol_x, y + symbol_y, pch = pch, col=symbol_col, bg=symbol_bg)

    }



