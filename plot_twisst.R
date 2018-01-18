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


plot.weights <- function(weights_dataframe,positions=NULL,line_cols=NULL,fill_cols=NULL,density=NULL,lwd=1,xlim=NULL,ylim=c(0,1),stacked=FALSE,
                                        ylab="Weights", xlab = "Position", main="",xaxt=NULL,yaxt=NULL,bty="n", add=FALSE){
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
