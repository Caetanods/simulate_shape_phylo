## Script to simulate morph in phylogeny.
## Polly, 2004 and Adams 2014 series of papers does not use the coordinates to simulate shape in a phylogeny.
## This is a simulation framework using multivariate normal distributions. Each step change of the BM model is draw from independent multivariate normal realizations. Covariance among structures can be achieved by using a vcv matrix.
## Code is biuld upon packages 'geiger' and 'geomorph'.

## Load packages>
library(geomorph)
library(geiger)
library(MASS)

## Get tps coordinates for ancestral state.

#digitize2d(filelist = "./Alestes.jpeg", nlandmarks = 24, scale = 10, tpsfile = "Alestes.tps")
anc <- readland.tps(file = "data/Alestes.tps")

## Simulate the tree:
phy <- sim.bdtree(stop = "taxa", n = 10)
phy <- rescale(phy, model = "depth", 1)
plot(phy, direction = "upward"); axisPhylo(side = 2)

## Function to simulate shape:
source("functions/sim_morph_functions.R")
rr <- nrow(anc)
mm <- diag(0.2, nrow = rr)
res <- sim.geo.char(phy, par = mm, tps = anc[,,1], model = "BM")
plotAllSpecimens(anc)









