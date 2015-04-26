sim.geo.char <- function (phy, par, tps, model = c("BM", "speciational"), nsim = 1) {
	## Version of 'sim.char' from geiger. Modified to evolve 2D geometric structures.
	## phy = phylogeny of class phylo
	## par = var covar matrix: dim equal to the number of coordinates. Recycle var and covar values. Each x,y par need to evolve at the same rate and with the same covariation in this model.
	## nsim = number of simulations
	## model = BM or speciational
	## tps = the procrustres coordinates for the shape of the root.

	## Known bugs: will break if nsim != 1.
	## Things to include: return the shapes at the nodes.
	## Things to expand: different rates for x and y coordinates (Ben's idea). Maybe elongated structures, for example, might evolves in different rates in one specific direction.
	## Function need to work with vectors of selection. Now the change in the landmarks have random directions. Vectors of selection make it possible to have changes bias towards specific directions.

    model <- match.arg(model, c("BM", "speciational"))
    model.matrix <- geiger:::.make.modelmatrix(par, model)
    nbranches <- nrow(phy$edge)
    nspecies <- Ntip(phy)

	## Need to improve the data format test here.
	## If user gives tps data in 'geomorph' format the error is uninformative.
    if (ncol(tps) != 2)
        stop("'tps' should be a 2-column matrix")

    m <- geiger:::.get.simulation.matrix(phy)
    if (model == "speciational") {
        m[m > 0] <- 1
    }

	nchar <- nrow(tps)
	## This will use the 'mvrnorm' function to simulate a vector of multivariate normal distributions all traits with mean 0 and varcovar matrix equal to the model.matrix (the matrix given as the par for the function).
    rnd1 <- mvrnorm(nsim * nbranches, mu = rep(0, nchar), Sigma = model.matrix)
	rnd2 <- mvrnorm(nsim * nbranches, mu = rep(0, nchar), Sigma = model.matrix)
	## This vector will be used for the changes of the BM model on each branch.
	rnd1 <- array(rnd1, dim = c(nchar, nbranches, nsim))
	rnd2 <- array(rnd2, dim = c(nchar, nbranches, nsim))
	## This is the simulation function. The product of the C matrix and the changes from the multivariate normal PLUS the value of the root. Then, we are changing (by addition) the value of the root (by a plus or a minus signal), proportional to the branch lenghts at each branch a time.
	## This need to change because root is not a single value. However, I can sum the value to a root equal zero after all calculations are done.
    simulate <- function(v) (m %*% as.matrix(v)) ## Root value fixed to 0. Sum made at the end.

	sim1 <- apply(rnd1, 1, simulate)
	sim2 <- apply(rnd2, 1, simulate)
	
	## Create a list of simulated tps values for each species.
	sp.sim <- lapply(1:nspecies, FUN = function(x) cbind(sim1[x,], sim2[x,]))
	## Reduce is used to make matrix sum. sum BM changes to root states.
	sp.sim <- lapply(sp.sim, FUN = function(x) Reduce('+', list(x, tps) ) )

	## Create the output in the tps format for geomorph package.
	## Output should be already in correct format here. Need improvement.
	result <- array( dim=c(nchar,2,nspecies) )
	for(i in 1:nspecies)  result[,,i] <- sp.sim[[i]]

	dimnames(result)[[3]] <- phy$tip.label	
	return(result)
}
