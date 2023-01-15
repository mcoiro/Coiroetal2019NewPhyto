## Load libraries
library(mapplots)
library(geiger)
library(ape)
library(phytools)
library(paleotree)
library(strap)
library(phangorn)
library(dispRity) 

## Upload the tree
tree <- read.nexus("tree.nex")

## Read the table with the pollen data
dat <- read.table("data.csv", header=TRUE, sep=";")

## Make a vector of the aperture type data
x <- as.factor(dat$aperture_type)
names(x) <- dat$species
x <- x[!is.na(x)]

## Trim the tree to only include the taxa with data avaiable and rearrange taxon names in the data vector
pruned.tree = drop.tip(tree, setdiff(tree$tip.label, names(x)))
x = x[match(pruned.tree$tip.label, names(x))]


## Simulate 200 stochastic character maps
chartrees <- make.simmap(pruned.tree, x,
                         model="ARD", nsim = 500)
## Output summary information
res_simmap <- describe.simmap(chartrees, plot = FALSE)


## Save the ancestral state reconstruction from the stochastic mapping
nodeace <- res_simmap$ace
## Remove Gymnosperms from the tree and reconstruction file, to obtain an LTT containing only angiosperm pollen
gymno = c("Cycas","Ginkgo_biloba","Metasequoia","Pinus","Gnetum_spp.","Welwitschia")
gymno.node = findMRCA(pruned.tree, tips= gymno)
gymno.desc = getDescendants(pruned.tree, gymno.node)
angio =c("Amborella_trichopoda","Lemna")
angio.node = findMRCA(pruned.tree, tips= angio)
pruned.tree = extract.clade(pruned.tree, angio.node)
## Remove nodes from the nodefile and reorganize the node names to correspond to the pruned tree
nodeace1 = nodeace[setdiff(rownames(nodeace),c(gymno.node,gymno.desc)),]
nodeace1 = nodeace1[-1,]
rownames(nodeace1) = c(608:1213)

## Make function for plotting pie charts
	plotNodePie = function(xvals,yvals,pie,pie.colors,radius=radius,...){
	plotPie = function(temp,pie.colors){
	add.pie(x=temp["xvals"],y=temp["yvals"],z=temp[-c(1:2)],radius=radius,col=pie.colors,labels=NA,init.angle=0)
	}
	plot(xvals,yvals,type="n",...)
	temp.matrix = cbind(xvals,yvals,pie)
	apply(temp.matrix,MARGIN=1,plotPie,pie.colors=pie.colors)
	return("done")
	}

## Add the reconstruction at each node for the piecharts
	for.piecharts <- nodeace1

## Calculate ages of all nodes in the tree.
	time.vector <- max(branching.times(pruned.tree))-branching.times(pruned.tree)

## Estimate lineage diversities of all nodes in the tree.
diversity.vector1 <- estDiversity(pruned.tree,x=x,method="simulation", nsim=500,model="ARD")


## Make a vector of colors for the areas.
piecolors<- c("yellow","magenta" ,"green3","black","cyan")


## Plot it.
	plotNodePie(xvals= time.vector,yvals= diversity.vector1,pie=for.piecharts,pie.colors=piecolors,xlab="Time From Root",ylab="Lineage Diversity Estimate",radius=1.5, xlim=c(0,140), ylim=c(0,130))
	