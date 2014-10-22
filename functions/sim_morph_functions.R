sim.geo.char <- function (phy, par, tps, model = c("BM", "speciational"), nsim = 1) {
	## Version of 'sim.char' from geiger. Modified to evolve 2D geometric structures.
	## phy = phylogeny of class phylo
	## par = var covar matrix: dim equal to the number of coordinates. Recycle var and covar values. Each x,y par need to evolve at the same rate and with the same covariation in this model.
	## nsim = number of simulations
	## model = BM or speciational
	## tps = the procrustres coordinates for the shape of the root.

	## Things to include: return the shapes at the nodes.
	## Things to expand: different rates for x and y coordinates (Ben's idea). Maybe elongated structures, for example, might evolves in different rates in one specific direction.

    model <- match.arg(model, c("BM", "speciational"))
    model.matrix <- geiger:::.make.modelmatrix(par, model)
    nbranches <- nrow(phy$edge)
    nspecies <- Ntip(phy)

#    if (length(root) > 1) 
#        stop("'root' should be a single value")
    if (ncol(tps) != 2)
        stop("'tps' should be a 2-column matrix")

		## m here is the C matrix from the phylogeny.
    m <- geiger:::.get.simulation.matrix(phy)
    if (model == "speciational") {
        m[m > 0] <- 1 ## See? To get the speciational model you need to set branch lengths all equal to 1.
    }
    nchar <- nrow(tps) * 2 ## Number of traits to be simulated. 2 = ncol(tps)
	## This will use the 'mvrnorm' function to simulate a vector of multivariate normal distributions all traits with mean 0 and varcovar matrix equal to the model.matrix (the matrix given as the par for the function).
    rnd <- t(mvrnorm(nsim * nbranches, mu = rep(0, nchar), Sigma = model.matrix))
	## This vector will be used for the changes of the BM model on each branch.
    rnd <- array(rnd, dim = c(nchar, nbranches, nsim))
	## This is the simulation function. The product of the C matrix and the changes from the multivariate normal PLUS the value of the root. Then, we are changing (by addition) the value of the root (by a plus or a minus signal), proportional to the branch lenghts at each branch a time.
	## This need to change because root is not a single value. However, I can sum the value to a root equal zero after all calculations are done.
    simulate <- function(v) (m %*% as.matrix(v)) + 0 ## Root value fixed to 0
    result <- apply(rnd, 1, simulate) ## Apply the 'simulate' over rows.


    result <- aperm(array(result, dim = c(nspecies, nsim, nchar)), c(1, 3, 2))
	## aperm here is used to change the positions of indexes 3 and 2. (just like transpose).
    rownames(result) <- phy$tip.label ## That is all!
 	
	return(result)
}
