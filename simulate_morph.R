## Script to simulate morph in phylogeny.
## Polly, 2004 and Adams 2014 series of papers does not use the coordinates to simulate shape in a phylogeny.
## This is a simulation framework using multivariate normal distributions. Each step change of the BM model is draw from independent multivariate normal realizations. Covariance among structures can be achieved by using a vcv matrix.
## Code is build upon packages 'geiger' and 'geomorph'.

## Load packages
library(geomorph)
library(geiger)
library(MASS)
library(phytools)

## Get tps coordinates for ancestral state:
## tps data provided by Michael Burns (Michael.Burns at oregonstate dot edu)
#digitize2d(filelist = "./Alestes.jpeg", nlandmarks = 24, scale = 10, tpsfile = "Alestes.tps")
anc <- readland.tps(file = "data/Alestes.tps")

## Simulate the tree:
phy <- sim.bdtree(stop = "taxa", n = 10)
phy <- rescale(phy, model = "depth", 1)
plot(phy, direction = "upward"); axisPhylo(side = 2)

## Function to simulate shape:
source("functions/sim_morph_functions.R")

## Note that there is no covariance in this simulation.
## The vcv matrix is a diagonal matrix.
vcv <- diag(5.0, nrow = nrow(anc))
## Make the simulation:
res <- sim.geo.char(phy, par = vcv, tps = anc[,,1], model = "BM")

## Plotting results:
jpeg("phylomorphspace.jpeg", width = 960, quality = 90)
par(mfrow = c(1,2))
plotGMPhyloMorphoSpace(phy, res, ancStates = FALSE) ## Does this use phyloPCA?
plot(phy, direction = "upward"); axisPhylo(side = 2)
dev.off()

## Plot difference in shape from root for one species:
dev.new()
plotRefToTarget(anc[,,1], res[,,1], method = "vector")

## Plot all the fishes in respect to the root:
nm <- paste("Evol_Alestes_",1:10,".jpeg",sep="")
for(i in 1:dim(res)[3]){
    jpeg(nm[i], quality = 90)
    plotRefToTarget(anc[,,1], res[,,i], method = "vector")
    dev.off()
}
