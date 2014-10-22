## Script to simulate morph in phylogeny.
## Polly, 2004 and Adams 2014 series of papers does not use the coordinates to simulate form in a phylogeny.
## Polly's method simulate in PC space, but use the direction of selection and the G matrix to make generations to generations simulations. This is very good for ecological time, but can be complicate for phylogenies.
## Adams' method is only used to understand the rates of evolution. He does not simulate in PC space, but do a indepdendent points simulation to check the rates of the shape evolution under BM models.

## What I want to do here is to do a 2D evolution of continuous characters.
## In a first test, there will be no morphological integration.
## Each coordinate will be treated as a independent point.
## Will start with a known shape and see where the phylogeny goes without morphological integration.

library(geomorph)
library(geiger)
library(MASS)

## Get ancestral state.

#digitize2d(filelist = "./Alestes.jpeg", nlandmarks = 24, scale = 10, tpsfile = "Alestes.tps")
#dev.copy2pdf()
anc <- readland.tps(file = "data/Alestes.tps") ## Load the tps data from file.
#plotAllSpecimens(anc) ## Simple landmark plot.

## Simulate the tree:
phy <- sim.bdtree(stop = "taxa", n = 100)

## For this version of the simulation I just need to get everything under a BM model.
mm <- anc
col1 <- mm[,1,1]
col2 <- mm[,2,1]

## Starting with the modified version of geiger 'sim.char':
source("functions/sim_morph_functions.R")
rr <- 
mm <- diag(2.0, nrow = )
res <- sim.geo.char(phy, par = mm, tps = anc, model = "BM")

anc[,,1]
