library(sets)
library(LaplacesDemon)
library(rjson)
setwd("GraphicalModel/edu_GraphicalModel")
data = read.csv("competitor_pack/data/fire-data-2.csv", nrows = 100000, check.names = FALSE)
#not DP, just using this for a test. user would specify domain normally.
mapping = fromJSON(file="competitor_pack/data/fire-data-specs-mapping.json")
data = data[,c(names(data) %in% names(mapping))]

epsilon = 1/3
delta = 1e-6
sensitivity = 1
max_domain_size = 100000


isDP <- function(params, d_g, e_g, err){
  # Function called by update parameters when using 
  # optimal composition approximation
  #
  # Args:
  #	params: a kx2 matrix where the first column corresponds to epsilon_i values and the second 
  # 			corresponds to delta_i values. 
  #	d_g: global delta value
  #	e_g: global epsilon value
  #   err: error parameter
  # Returns:
  #	Boolean indicating whether the given parameters satisfy (e_g,d_g)-DP or not. Since this is 
  #   an approximation (with error governed by err) it is not always correct. However, it will never 
  #   say something satisfies privacy when it doesn't. 
  
  elist <- params[,1] 
  dlist <- params[,2]
  k <- length(elist)
  
  # First check if the KOV composition theorems can decide privacy:
  if(max(elist) == min(elist) & max(dlist) == min(dlist)){
    eps <- KOVhom(params, d_g)
    if(eps <= e_g){
      return(T)
    }
    else {
      return(F)
    }
  }
  else {
    eps <- KOVhet(params, d_g)
    if(eps <= e_g){
      return(T)
    }
  }
  
  # sort parameters by epsilon value (needed for dynamic programming step)
  params <- params[order(params[,1]), ] 
  
  # discretize privacy parameters and e_g
  beta <- err/(k+sum(elist)+1)
  e0 <- log(1+beta)
  as <- ceiling(elist*((1/beta)+1))
  astar <- floor(e_g/e0)
  return(isDPinternal(as, beta, dlist, astar, d_g))
  
}


knapsack <- function(k, B, as, weights){
  # Helper function for isDP and isDPinternal. Runs dynamic programming procedure 
  # described in [MV16] Lemma 5.1.  
  
  prevrow <- rep(1,times=(B+1))
  for(r in 1:length(as)){
    ar <- as[r]
    wr <- weights[r]
    prevrow[(ar+1):(B+1)] <- prevrow[(ar+1):(B+1)] + wr*prevrow[1:(B+1-ar)]	
  }
  return(prevrow[B+1])
}


isDPinternal <- function(as, beta, dlist, astar, d_g){
  # Helper function for approxComp. Decides given the input privacy
  # parameters and a guess for global epsilon, whether or not privacy
  # is satisfied.
  
  k <- length(as)
  sum <- sum(as)
  
  # beta lets us convert between integer as and epsilons, 
  base <- 1+beta
  eterm <- prod(1+base^as)
  dterm <- prod(1-dlist)
  
  #Right hand side of optimal composition condition
  RHS <- eterm*(1-((1-d_g)/dterm))
  coef1 <- base^sum
  coef2 <- base^astar
  B <- floor((sum-astar)/2)
  
  # For dynamic programming, only need as that are less than B
  as <- subset(as, as<=B)
  # If all as are larger than B, privacy is trivially achieved:
  if(length(as)==0){
    #only nonzero set contains all epsilons
    LHS <- coef1-coef2
  }
  else{
    # Set weights for dynamic program as described in [MV16] Lemma 5.2
    weights1 <- base^(-as)
    weights2 <- base^as
    
    # run dynamic programming procedure
    t1 <- knapsack(k, B, as, weights1)
    t2 <- knapsack(k, B, as, weights2)
    
    # outputs of the dynamic programming give us left hand side of
    # optimal composition condition.
    LHS <- coef1*t1 - coef2*t2
  }
  if(LHS <= RHS){
    return(T)
  }
  else{
    return(F)
  }
}


approxComp <- function(params, d_g, err=.01){
  # Yields an approximation of optimal global epsilon for the given privacy
  # parameters with additive error err and multiplicative error of 2^(-err/2) on
  # global delta (Theorem 1.7 in Murtagh, Vadhan '16) 
  #
  # Args:
  #	params: a kx2 matrix where the first column corresponds to epsilon_i values and the second 
  # 			corresponds to delta_i values. 
  #	d_g: global delta value
  #   err: error parameter
  # Returns:
  #	Approximation of optimal global epsilon for the composition of the input parameters
  
  # sort parameters by epsilon value (needed for dynamic programming step)
  params <- params[order(params[,1]), ] 
  elist <- params[,1] 
  dlist <- params[,2]
  k <- length(elist)
  
  # if the privacy parameters are homogeneous (they do not differ across statistics)
  # just use the homogeneous optimal composition theorem KOV15.
  if(max(elist) == min(elist) & max(dlist) == min(dlist)){
    return(KOVhom(params, d_g))
  }
  
  # Constrain binary search by the best available rapidly computable upper bound:
  epsUpBound <- KOVhet(params, d_g)
  
  # Discretize epsilon values as described in [MV16]
  beta <- err/(k+sum(elist)+1)
  e0 <- log(1+beta)
  as <- ceiling(elist*((1/beta)+1))
  
  # begin binary search for optimal epsilon
  l <- 0
  u <- ceiling(epsUpBound/e0)
  done <- F
  
  # After termination astar is minimum integer such that e_g=astar*e0 satisfies privacy
  count <- 0
  while(!done){
    astar <- l+floor(((u-l)/2))
    
    # Check if current astar*e0 satisfies privacy:
    dp <- isDPinternal(as, beta, dlist, astar, d_g)
    
    # If it does, we can lower astar
    if(dp){
      if(u-l==1){
        # astar is the answer
        done <- T
      }
      else{
        u <- astar
      }
    }
    
    # If we are not satisfying privacy, must try a larger astar
    else {
      if(u-l<=2){
        # astar+1 is the answer
        astar <- astar+1
        done <- T
      }
      else{
        l <- astar	
      }
    }
  }
  return(astar*e0)
}



computeLHS <- function(eguess, k, eps){
  # helper function for KOVhom to compute 
  # the left hand side in Theorem 3.3 [KOV15]
  l <- ceiling((eguess+k*eps)/(2*eps))
  sum <- 0
  for(i in l:k){
    sum <- sum+choose(k,i)*(exp(i*eps)-exp(eguess+eps*(k-i)))
  }
  return(sum)
}

KOVhom <- function(params, d_g){
  # Computes the optimal composition theorem for the case when privacy parameters
  # do not differ across statistics. Theorem 3.3 [KOV15]. 
  # 
  # Args:
  #	params: a kx2 matrix where the first column corresponds to epsilon_i values and the second 
  # 			corresponds to delta_i values. 
  #	d_g: global delta value
  #   
  # Returns:
  #	global epsilon value guaranteed from the composition
  
  k <- length(params[,1])
  eps <- params[1,1]
  del <- params[1,2]
  eterm <- (1+exp(eps))^k
  dterm <- (1-del)^k
  RHS <- eterm*(1-((1-d_g)/dterm))
  u <- k*eps
  l <- 0
  LHS <- Inf
  e_g <- u
  
  # Must do binary search to approach optimal value
  count <- 0
  
  # Completely arbitrary cutoff of 50 rounds. Should have a more principled approach.
  while(count<50){
    eguess <- l+((u-l)/2)
    LHS <- computeLHS(eguess,k,eps)
    
    # If eguess satisfies privacy:
    if(LHS <= RHS){
      e_g <- eguess
      u <- eguess
    }
    else{
      l <- eguess
    }
    count <- count + 1
  }
  return(e_g)
}


KOVhet <- function(params, d_g, print=FALSE){
  # Computes an upper bound of optimal composition (Theorem 3.5, [KOV15]). Works for 
  # different settings of epsilons and deltas (heterogenous case). 
  # Args:
  #	params: a kx2 matrix where the first column corresponds to epsilon_i values and the second 
  # 			corresponds to delta_i values. 
  #	d_g: global delta value
  #   print: Boolean that if TRUE will print each individual term of the theorem rather than just
  #          the minimimum.
  #   
  # Returns:
  #	global epsilon value guaranteed from the composition
  
  elist <- params[,1]
  dlist <- params[,2]
  k <- length(elist)
  del_bar <- 1-((1-d_g)/prod(1-dlist))
  sum_of_squares <- sum(elist^2)
  
  # Function that will be applied to the vector of epsilons
  fun <- function(x){
    return(((exp(x) - 1)*x)/(exp(x) + 1))		
  }
  
  first_term <- sum(sapply(elist, FUN=fun))
  
  # a, b, and c will correspond to the first, second, and third expressions in 
  # theorem 3.5 in KOV15. The minimum is returned.
  
  a <- sum(elist)
  b <- first_term + sqrt((2*log(exp(1) + (sqrt(sum_of_squares)/del_bar)))*sum_of_squares)
  c <- first_term + sqrt(2*log(1/del_bar)*sum_of_squares)
  
  # For testing purposes if one wants to print all three terms
  if(print){
    cat("\nFirst term: ", a)
    cat("\nSecond term: ", b)
    cat("\nThird term: ", c)
    cat("\nFinal result: ", min(a,b,c), "\n")
  }
  
  vec <- c(a,b,c)
  
  #If any of the outputs are NaN, return the minimum of the actual numerical results
  if(sum(!is.na(vec)) == 0){
    return(NaN)
  }
  else{
    return(min(vec[which(!is.na(vec))]))
  }	
}

update_parameters <- function(params, hold, eps, del){
  #
  # Args:
  #	params: kx2 matrix of privacy parameters where column one corresponds
  #			to epsilons and column two is deltas.
  #	hold: vector of indices corresponding to rows of params that will not 
  #		   be updated, either because they were just added or because the 
  #		   user has requested that these values stay fixed (Hold feature). 
  #	       If we are to update every parameter, set hold to 0. 
  #	eps: global epsilon
  #	del: global delta
  #
  # Returns:
  #	kx2 matrix of updated parameters
  
  
  # k is the total number of statistics currently selected
  k <- length(params[ , 1])
  
  # If we are only calculating one statistic, allocate the whole budget to it.
  if(k == 1){
    params[1,1] <- eps
    params[1,2] <- del
    return(params)
  }
  
  elist <- as.numeric(params[ ,1])
  dlist <- as.numeric(params[ ,2])
  #hard coded error tolerance for optimal composition approximation. Might do something more clever one day. 
  #err <- eps/10
  err <- .01
  # Check if there are unset epsilon values
  unsetEpsilons <- c(which(is.na(elist)), which(elist==" "), which(elist==""), which(elist==0))
  unsetEpsilons <- unique(unsetEpsilons)
  
  # Get list of unheld (or free) epsilon indices
  indices <- seq(1,k,1)
  free <- indices[! indices %in% hold]
  newEps <- -10
  #if no epsilons are set: 
  if(length(unsetEpsilons)==length(elist)){
    tempElist <- rep(1,times=length(elist))
    newElist <- scale_eps(tempElist, dlist, eps, del, free, err)
    toReturn <- cbind(newElist, dlist)
    return(toReturn)
  }
  else if(length(unsetEpsilons) > 0){
    # collect all of the free epsilons that have been set
    toAvg <- free[! free %in% unsetEpsilons]
    tempElist <- elist
    # replace all free epsilons with their average
    tempElist[free] <- mean(tempElist[toAvg])
    newElist <- scale_eps(tempElist, dlist, eps, del, free, err)
    neweps <- newElist[free][1]
    if(min(elist[toAvg])==max(elist[toAvg])){
      # if all free unset epsilons are the same then newElist the right update
      toReturn <- cbind(newElist, dlist)
      return(toReturn)
    }
    else{
      #set all unset epsilons to neweps and remove them from the free list
      elist[unsetEpsilons] <- neweps
      free <- free[! free %in% unsetEpsilons]
    }
  }
  # if some epsilons are being held, check that the held ones alone don't exceed the budget:
  else if(length(elist)>length(free)){
    tempElist <- elist
    tempElist[free] <- 0
    dp <- isDP(cbind(tempElist,dlist), del, eps, err)
    if(!dp){
      return("error")
    }
  }
  
  newElist <- scale_eps(elist, dlist, eps, del, free, err)
  toReturn <- cbind(newElist, dlist)
  return(toReturn)
}


scale_eps <- function(elist, dlist, eps, del, free, err){
  # This function returns a list of epsilon values
  # where each epsilon in elist that is not being held
  # is scaled by the same multiplicative factor until 
  # composition is satisfied
  
  # Initialize parameters for binary search. 
  l <- 0
  u <- eps/max(elist)  # no epsilon value in the composition can exceed global eps
  dp <- F
  
  goodlist <- c()
  
  # is there a better stopping condition?
  while(u-l>.0001){	
    # scaling factor		
    r <- l + ((u - l)/2)
    testElist <- elist
    testElist[free] <- r*testElist[free]
    test_params <- cbind(testElist, dlist)
    dp <- isDP(test_params, del, eps, err)
    
    # Reset upper and lower bounds of the search depending on outcome of isDP
    if(dp){
      l <- r
      goodlist <- testElist
    }
    else {
      u <- r
    }
  }
  # If result does not beat simple summing composition
  total <- sum(goodlist)
  if(total < eps){
    toadd <- (eps-total)/length(free)
    toaddlist <- rep(0,times=length(goodlist))
    toaddlist[free] <- toadd
    goodlist <- goodlist + toaddlist
  }
  return(goodlist)	
}


rlap = function(mu=0, b=1, size=1) {
  p <- runif(size) - 0.5
  draws <- mu - b * sgn(p) * log(1 - 2 * abs(p))
  return(draws)
}

mechanism <- setRefClass(
  Class = 'mechanism',
  fields = list(
    mechanism = 'character',
    name = 'character',
    variable = 'character',
    var.type = 'character',
    var.type.orig = 'character',
    n = 'numeric',
    epsilon = 'numeric',
    delta = 'numeric',
    rng = 'ANY',
    result = 'ANY',
    alpha = 'numeric',
    accuracy = 'numeric',
    bins = 'ANY',
    n.bins = 'ANY',
    k = 'numeric',
    error = 'numeric',
    n.boot = 'ANY',
    boot.fun = 'function',
    impute.rng = 'ANY',
    impute = 'logical',
    formula = 'ANY',
    columns = 'ANY',
    intercept = 'logical', 
    stability = 'logical',
    objective = 'function',
    gran = 'numeric',
    percentiles = 'ANY',
    tree.data = 'ANY'
  ))

mechanism$methods(
  getFields = function() {
    f <- names(getRefClass()$fields())
    out <- setNames(vector('list', length(f)), f)
    for (fd in f) {
      out[[fd]] <- .self[[fd]]
    }
    return(out)
  })

mechanism$methods(
  getFunArgs = function(fun) {
    f <- .self$getFields()
    spec <- list()
    for (arg in names(f)) {
      if (arg %in% names(formals(fun))) {
        spec[[arg]] <- f[[arg]]
      }
    }
    return(spec)
  })

mechanismExponential <- setRefClass(
  Class = 'mechanismExponential',
  contains = 'mechanism'
)

#mechanismExponential$methods(
#  getFunArgs = function(fun) {
#    callSuper(fun)
#  })

mechanismExponential$methods(
  evaluate = function(fun, x, sens, postFun, ...) {
    x <- censordata(x, .self$var.type, levels=.self$bins)
    x <- fillMissing(x, .self$var.type, categories=.self$bins)   # Problem that this needs to work over other types than categorical
    field.vals <- .self$getFunArgs(fun)
    ellipsis.vals <- getFuncArgs(list(...), fun)
    #print(field.vals)
    #print(ellipsis.vals)
    true.val <- do.call(fun, c(list(x=x), field.vals, ellipsis.vals))
    quality <- true.val - max(true.val)
    likelihoods <- exp((.self$epsilon * quality) / (2 * sens))   # Problem that this zeroes out for quality differences > 800
    probs <- likelihoods/sum(likelihoods)
    #print(names(true.val))
    release <- sample(names(true.val), size=.self$k, prob=probs) # Problem that this needs to use openssl randomness
    out <- list('release' = release)
    out <- postFun(out, ...)
    return(out)
  })


dpUnif <- function(n, seed=NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
    return(runif(n))
  }
  return(openssl::rand_num(n))
}


censordata <- function(x, var_type, rng=NULL, levels=NULL) {
  if (var_type %in% c('character', 'factor')) {
    if (is.null(levels)) {
      x <- factor(x, exclude=NULL)
    } else {
      x <- factor(x, levels=levels, exclude=NULL)
    }
  } else {
    if (is.null(rng)) {
      stop('range `rng` is required for numeric types')
    }
    if (NCOL(x) > 1) {
      for (j in 1:ncol(x)) {
        rng[j, ] <- checkrange(rng[j, ])
        x[, j][x[, j] < rng[j, 1]] <- rng[j, 1]
        x[, j][x[, j] > rng[j, 2]] <- rng[j, 2]
      }
    } else {
      rng <- checkrange(rng)
      x[x < rng[1]] <- rng[1]
      x[x > rng[2]] <- rng[2]
    }
  }
  return(x)
}

fillMissing1D <- function(x, var.type, lower=NULL, upper=NULL, categories=NULL) {
  naIndices <- is.na(x)         # indices of NA values in x
  nMissing <- sum(naIndices)    # number of missing values
  
  if (nMissing == 0) { 
    return(x) 
  }
  
  u <- dpUnif(nMissing) # array of uniform random numbers of length nMissing
  scaledVals <- scaleValues(u, var.type, lower, upper, categories) # scale uniform vals
  x[naIndices] <- scaledVals #assign to NAs in input array
  return(x)
}

fillMissing2D <- function(x, var.type, impute.rng=NULL) {
  for (j in 1:ncol(x)) {
    x[, j] <- fillMissing1D(x[, j], var.type, lower=impute.rng[j, 1], upper=impute.rng[j, 2])
  }
  return(x)
}

fillMissing = function(x, var.type, impute.rng=NULL, categories=NULL) {
  if (var.type %in% c('numeric', 'integer', 'logical')) {
    if (NCOL(x) > 1) {
      x <- fillMissing2D(x, var.type, impute.rng)
    } else {
      x <- fillMissing1D(x, var.type, impute.rng[1], impute.rng[2])
    }
  } else {
    x <- fillMissing1D(x, var.type, categories=categories)
  }
  return(x)
}

scaleValues = function(vals, var.type, lower=NULL, upper=NULL, categories=NULL) {
  if (var.type %in% c('character', 'factor')) { 
    lower <- 1
    upper <- length(categories)
  }
  
  if (var.type == 'logical') {       # logical values can only be 0 or 1 so set bounds accordingly
    lower <- 0
    upper <- 2                       # upper bound of 2 not 1 because as.integer always rounds down.
  }
  
  out <- vals * (upper - lower) + lower  # scale uniform random numbers based on upper and lower bounds
  
  if (var.type %in% c('logical', 'integer')) { # if logical or integer, trim output to integer values
    out <- as.integer(out)
  } else if(var.type == 'logical') {
    
  } else if (var.type %in% c('character', 'factor')) { # if character or factor, assign output to categories.
    out <- categories[as.integer(out)]
  }
  
  return(out)
}

getFuncArgs <- function(output, target.func) {
  spec <- list()
  for (element in names(output)) {
    if (element %in% names(formals(target.func))) {
      spec[[element]] <- output[[element]]
    }
  }
  return(spec)
}

mutual_information_vector = function(data, x, type = "L1") {
  vec_type = vector(mode = "numeric", length = 1)
  mutual_informations = vapply(x, mutual_information, vec_type, data, type)
  names(mutual_informations) = x
  #print(mutual_informations)
  return(mutual_informations)
}

mutual_information = function(pair, data, type = "L1") {
  pair = as.character(pair)
  #print(pair)
  if (pair == "empty") {
    return(0.2)
  } else{
    these_elements = strsplit(pair, ';')[[1]]
    setA = as.set(strsplit(these_elements[1],',')[[1]])
    setB = as.set(strsplit(these_elements[2],',')[[1]])
    #intersection = set_intersection(setA,setB)
    #setA = set_complement(intersection, setA)
    #setB = set_complement(intersection, setB)
    a = as.vector(setA, mode = "integer")
    b = as.vector(setB, mode = "integer")
    #int = as.vector(intersection, mode = "integer")
    #print(these_elements)
    #print(a)
    #print(b)
    n = dim(data)[1]
    if_independent = outer(as.vector(table(data[,a])), as.vector(table(data[,b])))/n^2
    dim = c(length(table(data[,a])),length(table(data[,b])))
    true_val = array(table(data[,c(a,b)])/n, dim = dim, dimnames = NULL)
    if (type == "L1") {
      return(sum(abs(if_independent-true_val)))
    }
    
    else if (type == "KL") {
      return(KLD(if_independent,true_val)$sum.KLD.py.px)
    }
  }
  
}


#trim whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

identity = function(a,...) {
  return(a)
}

get_domains = function(index, domain, att) {
  return(length(domain[[att[index]]]))
}


select_queries = function(data, max_domain_size, domain, epsilon, delta, queries = NULL) {
  n = dim(data)[1]
  total_mutual_information = 0
  d = length(names(data))
  att = names(data)
  adj = list()
  for (i in 1:d) {
    adj[[as.character(i)]] = as.set(as.character(i))
  }
  type = "KL"
  sensitivity = 2/n*log(n+1/2)+(n-1)/n*log((n+1)/(n-1))
  Q = as.set(queries)
  c = length(queries)
  list_of_pairs = vector(mode = "character")
  #utility = vector()
  #start with just the points and cliques already in the graph
  for (i in 1:d) {
    for (j in 1:i){
      if (j < i && as.double(length(domain[[att[i]]]))*as.double(length(domain[[att[j]]]))< max_domain_size) {
        #utility = c(utility, compute_mutual_information(data, i, j))
        list_of_pairs = c(list_of_pairs, paste(as.character(c(i,j)), collapse = ";"))
      }
    }
  }
  for (q1 in queries) {
    for (q2 in Q) {
      if (q1 != q2) {
        c1 = as.numeric(strsplit(q1, ",")[[1]])
        domain_size_1 = prod(sapply(c1, get_domains, domain, att))
        elements_1 = as.set(c1)
        c2 = as.numeric(strsplit(q2, ",")[[1]])
        domain_size_2 = prod(sapply(c2, get_domains, domain, att))
        elements_2 = as.set(c2)
        if (domain_size_1*domain_size_2 < max_domain_size && set_is_empty(set_intersection(elements_1,elements_2))) {
          new_pair = paste(as.character(c(q1, q2)), collapse = ";")
          list_of_pairs = c(list_of_pairs, paste(new_pair))
        }
      }
      
    }
  }
  list_of_pairs = c(list_of_pairs, "empty")
  set_of_pairs = as.set(list_of_pairs)
  init = rep(c(epsilon/(2*(d-c)),0),2*(d-c))
  params = matrix(init, nrow = 2*(d-c), ncol = 2, byrow = TRUE)
  optimal_values = update_parameters(params = params, hold=0, eps = epsilon, del=delta)
  best_epsilon = optimal_values[1,1]
  best_delta = optimal_values[1,2]
  #print(list_of_pairs)
  #Now choose queries with exponential mechanism, we take pairs of cliques
  #cliques build up over time so we keep adding to the set that we're looking at
  for (i in 1:(2*(d-c))) {
    #print(i)
    myMech = mechanismExponential()
    myMech$k = 1
    myMech$var.type = "character"
    myMech$n = n
    myMech$epsilon = best_epsilon
    myMech$delta = best_delta
    myMech$bins = list_of_pairs
    pair = myMech$evaluate(fun = mutual_information_vector, x = list_of_pairs, sens = sensitivity, postFun = identity, data = data, type = type)$release
    #print(pair)
    total_mutual_information = total_mutual_information + mutual_information(pair, data, type = "KL")
    #print(mutual_information(pair, data, type = "KL"))
    if (pair != "empty") {
      list_of_pairs = list_of_pairs[list_of_pairs != pair]
      S = as.set(strsplit(pair, ';')[[1]])
      #Make sure each element of S is actually a single attribute
      S = paste(as.character(S), collapse = ",")
      S = as.set(strsplit(S, ',')[[1]])
    
    #print(4*d/(best_epsilon*n))
    
    #u = mutual_information(data, pair, "KL") + rnorm(n=1, mean=0, sd = 2*d/(best_epsilon*n))
      #If they all share a common neighbor then adding this pairing completes the clique so we want to only consider that
      #for the mutual information gain (and thus probability) to be accurate.
      for (j in S) {
        adj[[j]] = set_union(adj[[j]], S)
      }
      for (j in S) {
        for (pairing in list_of_pairs) {
          these_elements = strsplit(pairing, ';')[[1]]
          setA = as.set(strsplit(these_elements[1],',')[[1]])
          setB = as.set(strsplit(these_elements[2],',')[[1]])
          these_elements = set_union(setA, setB)
          count = 0
          #The above happens exactly when we're one edge off from a complete graph, so we can count how far off the degree of this
          #subgraph is from being complete and we get rid of the pair iff we're off by 2
          for (k in these_elements) {
            this_length = length(set_intersection(adj[[j]], adj[[k]]))
            count = count + length(these_elements) - this_length + 1
          }
          if (count == 2) {
            list_of_pairs = list_of_pairs[list_of_pairs != pairing]
          }
        }
      }
      set_of_pairs = as.set(list_of_pairs)
      #remove all pairings that only use elements in pair
      #for (pairing in list_of_pairs) {
        #these_elements = strsplit(pairing, ';')[[1]]
        #setA = as.set(strsplit(these_elements[1],',')[[1]])
        #setB = as.set(strsplit(these_elements[2],',')[[1]])
        #these_elements = set_union(setA, setB)
        #if (set_is_subset(these_elements, S)) {
        #  list_of_pairs = list_of_pairs[list_of_pairs != pairing]
        #}

      #}
      new_clique = paste(as.character(S), collapse = ",")
      #Add all pairs with the new set under the domain size bound into consideration
      #We don't include pairs that intersect since those are redundant.
      for (clique in Q) {
        cliqueSet = as.set(strsplit(clique,',')[[1]])
        indices = as.vector(set_union(cliqueSet,S), mode = "integer")
        #print(index)
        domain_size = 1
        for (i in indices) {
          domain_size = domain_size * length(domain[[att[i]]])
        }
        #print(domain_size)
        if (set_is_empty(set_intersection(cliqueSet, S)) && domain_size < max_domain_size) {
          #print("chosen")
          new_pair = paste(as.character(c(new_clique, clique)),collapse = ";")
          if (!set_contains_element(set_of_pairs, new_pair)) {
            list_of_pairs = c(list_of_pairs, new_pair)
            set_of_pairs = set_union(set_of_pairs, as.set(new_pair))
          }
          
        }
      }
      #print(list_of_pairs)
      Q = set_union(Q, as.set(new_clique))
      #print(Q)
      
    }
  }
  #Add 1-way marginals on isolated points to the final set of queries
  for (adjacencies in adj) {
    if (length(element == 1)) {
      Q = set_union(Q, element)
    } 
  }
  output = list()
  output$Q = Q
  output$MI = total_mutual_information
  return(output)
}