library(phangorn)
library(parallel)



subtree <- function(tree, tips){
  all_tips <- tree$tip.label
  for (tip in all_tips){
    if (!(tip %in% tips)){
      tree <- drop.tip(tree, tip)
      }
    }
  tree
  }



weight_tree <- function(tree, taxa, reps, topos = NULL){
  tree <- unroot(tree)
  if (is.null(topos) == TRUE){
    topos <- allTrees(length(taxa), tip.label=names(taxa))
    }
  
  #if there are not many combos we can do them exaustively, otherwise we subsample
  if (reps >= prod(sapply(taxa,length))){
    combs <- expand.grid(taxa, stringsAsFactors = FALSE)
    combs <- lapply(1:nrow(combs), function(x){combs[x,]})
    }else{
    combs <- lapply(1:reps, function(x){sapply(taxa, sample, 1)})
    }
    
  counts = numeric(length = length(topos))
  for (comb in combs){
    #sample one from each taxon and prune
    s <- subtree(tree, comb)
    s <- unroot(s)
    #rename labels for topology comparisson
    s$tip.label[sapply(1:length(taxa), function(x){which(s$tip.label %in% taxa[[x]])})] <- names(taxa)
    #find matching topo
    match <- which(sapply(topos, function(x){RF.dist(s,x)}) == 0)
    counts[match] <- counts[match] + 1
    }
  output <- list(topos=topos, counts=counts, weights=counts/length(combs))
  output
  }


simple.loess.predict <- function(x, y, span, weights = NULL, max = NULL, min = NULL){
  y.loess <- loess(y ~ x, span = span, weights = weights)
  y.predict <- predict(y.loess,x)
  if (is.null(min) == FALSE) {y.predict = ifelse(y.predict > min, y.predict, min)}
  if (is.null(max) == FALSE) {y.predict = ifelse(y.predict < max, y.predict, max)}
  y.predict
  }


####################


treedata <- read.delim("test.raxml.w1m100s1.tsv", sep = "\t", header = TRUE, as.is = TRUE)

treedata <- treedata[is.na(treedata$tree) == FALSE,]

trees <- lapply(treedata$tree, function(x){read.tree(text = x)})

###

taxa <-list(west=c("ros.MK523_A","ros.MK523_B","ros.MK524_A","ros.MK524_B","ros.MK525_A","ros.MK525_B","ros.MK589_A","ros.MK589_B","ros.MK675_A","ros.MK675_B","ros.MK676_A","ros.MK676_B","ros.MK682_A","ros.MK682_B","ros.MK683_A","ros.MK683_B","ros.MK687_A","ros.MK687_B","ros.MK689_A","ros.MK689_B","ros.CJ531_A","ros.CJ531_B","ros.CJ533_A","ros.CJ533_B","ros.CJ546_A","ros.CJ546_B","ros.CJ2071_A","ros.CJ2071_B","melP.CJ18038_A","melP.CJ18038_B","MelP.CJ18097_A","MelP.CJ18097_B","vul.CS519_A","vul.CS519_B","vul.CS10_A","vul.CS10_B","vul.CS11_A","vul.CS11_B","cyth.CJ2856_A","cyth.CJ2856_B"),east=c("moc.CS228_A","moc.CS228_B","moc.CS16_A","moc.CS16_B","moc.CS17_A","moc.CS17_B","ple.CJ9156_A","ple.CJ9156_B","mapl.CJ16042_A","mapl.CJ16042_B","mal.CJ17162_A","mal.CJ17162_B","mal.CS21_A","mal.CS21_B","mal.CS22_A","mal.CS22_B","mal.CS24_A","mal.CS24_B","ecu.CJ9117_A","ecu.CJ9117_B","ama.JM216_A","ama.JM216_B","ama.JM160_A","ama.JM160_B","ama.JM293_A","ama.JM293_B","ama.JM48_A","ama.JM48_B","agl.JM108_A","agl.JM108_B","agl.JM122_A","agl.JM122_B","agl.JM569_A","agl.JM569_B","agl.JM572_A","agl.JM572_B","aman.CS2228_A","aman.CS2228_B","aman.CS2221_A","aman.CS2221_B"),guiana=c("melG.CJ9315_A","melG.CJ9315_B","melG.CJ9316_A","melG.CJ9316_B","melG.CJ9317_A","melG.CJ9317_B","melG.CJ13435_A","melG.CJ13435_B","thel.CJ13566_A","thel.CJ13566_B"),cyd=c("cyd.CJ553_A","cyd.CJ553_B","cyd.CJ560_A","cyd.CJ560_B","cyd.CJ564_A","cyd.CJ564_B","cyd.CJ565_A","cyd.CJ565_B"),tim=c("tim.JM313_A","tim.JM313_B","tim.JM57_A","tim.JM57_B","tim.JM84_A","tim.JM84_B","tim.JM86_A","tim.JM86_B"),silv=c("hec.JM273_A","hec.JM273_B","eth.JM67_A","eth.JM67_B","par.JM371_A","par.JM371_B","ser.JM202_A","ser.JM202_B"))


topos <- allTrees(length(taxa), tip.label=names(taxa))

reps <- 500
weights <- do.call(rbind, mclapply(trees, function(t){weight_tree(t,taxa,reps)$weights}, mc.cores = 6))
overall <- apply(weights, 2, mean)

weights <- weights[,order(overall, decreasing = TRUE)]

topos <- topos[,order(overall, decreasing = TRUE)]






cols <- sample(colors(), 10)

plot(0, pch = "",xlim = c(0,length(trees)), ylim = c(0,1), ylab = "", xlab = "", bty = "n") 

x=1:length(trees)

for (n in 1:10){
  y = weights[,n]
  
  points(x,y,col = cols[n])
  
  y.predict <- simple.loess.predict(x,y,span = 0.1, max = 1, min = 0)

  lines(x,y.predict, lwd = 1.5, col = cols[n])
  }




