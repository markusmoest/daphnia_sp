#Ancestral state reconstruction for sperm length in Daphnia
rm(list=c(ls()))
library("phytools")
library("ape")
library("geiger")

# define and set path to working directory
PATH = 'SET_YOUR_WORKING_DIRECTORY_WITH_INPUT_FILES_HERE'
setwd(PATH)
OUTPUT = 'SET_YOUR_OUTPUT_DIRECTORY_HERE'

# define tree
TREE = 'DAPH_ALL_FINAL_PURE.treefile'
# define sperm length file
LENGTH = 'Measure_sperm_evolution_mean.csv'

#Read unrooted tree file and inspect tree
sperm.raw <- read.tree(TREE)
is.rooted(sperm.raw)
plot(sperm.raw)
sperm.raw$tip.label

# Root tree with Scapholeberis, remove root from plot and rotate nodes
sperm.rooted <- root(sperm.raw, c("Scapholeberis"), resolve.root=T)
plot(sperm.rooted)
sperm.tree <- drop.tip(sperm.rooted,"Scapholeberis") 
plot(sperm.tree)
sperm.tree <- rotateNodes(sperm.tree,c("20","21","22","24","25","27","29","30","31","32"))
plotTree(sperm.tree,type="phylogram", direction="rightwards",ftype="i", fsize=0.5, node.numbers=T)

# create dated ultrametric tree with multiple node calibration using Cornetti et al. 2019 age bounds for the mtDNA - fossil calibration - Late Jurassic dated tree (most conservative). 

node <- c(
  getMRCA(sperm.tree, tip = c("Dpulex","Dmagna") ),
  getMRCA(sperm.tree, tip = c("Dlongispina","Dpulex") ),
  getMRCA(sperm.tree, tip = c("Dmagna","Dsinensis") ),
  getMRCA(sperm.tree, tip = c("Dsimilis","Dsinensis") )
) 

age.min <- c(
  145,
  101.1,
  71.9,
  45.6
)

age.max <- c(
  145.5,
  108.3,
  76.5,
  50.1
)

soft.bounds <- c(
  FALSE,
  FALSE,
  FALSE,
  FALSE
)

calib <- data.frame(node, age.min, age.max, soft.bounds) 

#Compare ultrametric calibrated trees using different methods and lambdas from 0.1-1
comp_methods <- data.frame(method=character(), lambda=numeric(), logLik=numeric(), PHIIC=numeric(), convergence=logical())

k=1
for(lbd in seq(0,1,0.1)){
## relaxed
sperm.time.tree.rel <- chronos(sperm.tree, lambda = lbd, model = "relaxed", calibration = calib, control = chronos.control())
comp_methods[lbd*10+k,1] <- "sperm.time.tree.rel"
comp_methods[lbd*10+k,2] <- lbd
comp_methods[lbd*10+k,3] <- attributes(sperm.time.tree.rel)$PHIIC$logLik
comp_methods[lbd*10+k,4] <- attributes(sperm.time.tree.rel)$PHIIC$PHIIC
comp_methods[lbd*10+k,5] <- attributes(sperm.time.tree.rel)$convergence

## correlated
sperm.time.tree.cor <- chronos(sperm.tree, lambda = lbd, model = "correlated", calibration = calib, control = chronos.control())
comp_methods[lbd*10+1+k,1] <- "sperm.time.tree.cor"
comp_methods[lbd*10+1+k,2] <- lbd
comp_methods[lbd*10+1+k,3] <- attributes(sperm.time.tree.cor)$PHIIC$logLik
comp_methods[lbd*10+1+k,4] <- attributes(sperm.time.tree.cor)$PHIIC$PHIIC
comp_methods[lbd*10+1+k,5] <- attributes(sperm.time.tree.cor)$convergence
k=k+1
}
## discrete
sperm.time.tree.disc <- chronos(sperm.tree, lambda = 1, model = "discrete", calibration = calib, control = chronos.control())
comp_methods[lbd*10+1+k,1] <- "sperm.time.tree.disc"
comp_methods[lbd*10+1+k,2] <- lbd
comp_methods[lbd*10+1+k,3] <- attributes(sperm.time.tree.disc)$PHIIC$logLik
comp_methods[lbd*10+1+k,4] <- attributes(sperm.time.tree.disc)$PHIIC$PHIIC
comp_methods[lbd*10+1+k,5] <- attributes(sperm.time.tree.disc)$convergence

## discrete with strict clock
sperm.time.tree.disc_strict <- chronos(sperm.tree, lambda = 1, model = "discrete", calibration = calib, control = chronos.control(nb.rate.cat=1))
comp_methods[lbd*10+2+k,1] <- "sperm.time.tree.disc_strict"
comp_methods[lbd*10+2+k,2] <- lbd
comp_methods[lbd*10+2+k,3] <- attributes(sperm.time.tree.disc_strict)$PHIIC$logLik
comp_methods[lbd*10+2+k,4] <- attributes(sperm.time.tree.disc_strict)$PHIIC$PHIIC
comp_methods[lbd*10+2+k,5] <- attributes(sperm.time.tree.disc_strict)$convergence

sperm.time.tree.rel_l0 <- chronos(sperm.tree, lambda = 0, model = "relaxed", calibration = calib, control = chronos.control())
sperm.time.tree.rel_l1 <- chronos(sperm.tree, lambda = 1, model = "relaxed", calibration = calib, control = chronos.control())
sperm.time.tree.cor_l0 <- chronos(sperm.tree, lambda = 0, model = "correlated", calibration = calib, control = chronos.control())
sperm.time.tree.cor_l1 <- chronos(sperm.tree, lambda = 1, model = "correlated", calibration = calib, control = chronos.control())
sperm.time.tree.disc_l0 <- chronos(sperm.tree, lambda = 0, model = "discrete", calibration = calib, control = chronos.control())
sperm.time.tree.disc_l1 <- chronos(sperm.tree, lambda = 1, model = "discrete", calibration = calib, control = chronos.control())
sperm.time.tree.disc_strict_l0 <- chronos(sperm.tree, lambda = 0, model = "discrete", calibration = calib, control = chronos.control(nb.rate.cat=1))
sperm.time.tree.disc_strict_l1 <- chronos(sperm.tree, lambda = 1, model = "discrete", calibration = calib, control = chronos.control(nb.rate.cat=1))

########
#read trait file
sperm.length <- read.csv(LENGTH, header=T, row.names=1)

#create named vector for phytools
sl<-setNames(sperm.length$mean,rownames(sperm.length))

## add Ceriodaphnia sperm length value from Wingstrand
sl["Ceriodaphnia"] = 4.45  #mean sperm length for all Ceriodaphnia = 5.09; Ceriodaphnia reticulata mean sperm length = 4.45
## add Simocephalus sperm length value from Wingstrand
sl["Simocephalus"] = 2.775 #mean sperm length for all Simocephalus = 2.775

# Choose a tree:
#tree <- sperm.tree 
#tree <- sperm.time.tree.rel_l0 
tree <- sperm.time.tree.disc_strict_l1
plot(tree)

# assess phylogenetic signal
phylosig(tree,sl,test=TRUE)
phylosig(tree,sl,method="lambda",test=TRUE)

# fit and compare different models of trait evolution and compare; adjust bounds if warnings appear
fitBM<-fitContinuous(tree,sl)
fitOU<-fitContinuous(tree,sl,model="OU", bounds=list(alpha = c(min = exp(-500), max = exp(2))))
fitEB<-fitContinuous(tree,sl,model="EB", bounds=list(a = c(min = -2.637555e-02, max = -1.000000e-09)))

fitBM
fitOU
fitEB
aic.vals<-setNames(c(fitBM$opt$aicc,fitOU$opt$aicc,fitEB$opt$aicc), c("BM","OU","EB"))
aic.vals
aic.w(aic.vals) # choose the model with the highest aic.w

#run reconstruction with anc.ML and choose a model; best-fit often "OU" but does not alwways converge, the simpler "BM" also performs well, results are very similar
fit_sperm<-anc.ML(tree, sl, model = "BM", maxit=10000, vars=TRUE, CI=TRUE)
#fit_sperm<-fastAnc(sperm.4,sl.4, vars=TRUE, CI=TRUE)
print(fit_sperm,printlen=100)
fit_sperm
str(fit_sperm)
#check convergence
fit_sperm$convergence

#Plot the actual states from the fitted model above - I think that may be better, adjust the lims!! 
fit_plot <- contMap(tree,sl,method="user",anc.states=fit_sperm$ace, plot=T) # check that lims is set to min/max CI, otherwise the heatmap lacks contrast

setwd(OUTPUT)
pdf(file = "sperm.time.tree.disc_strict_l1_BM.pdf",height = 8.25, width = 5.875)
plot(fit_plot,type="phylogram",legend=0.7*max(nodeHeights(tree)),
     sig=2,fsize=c(0.8,0.9))
dev.off()

write.table(file = "sperm.time.tree.disc_strict_l1_BM_CI95.txt", x = fit_sperm$CI95)
str(fit_plot)

pdf(file = "sperm.time.tree.disc_strict_l1_BM_nodes.pdf",width = 8.25, height = 5.875)
plotTree(tree,type="phylogram", direction="rightwards",ftype="i", fsize=0.5, node.numbers=T)
dev.off()
print(fit_sperm,printlen=100)
