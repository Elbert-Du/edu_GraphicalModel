library(sets)
library(LaplacesDemon)
#setwd("GraphicalModel/edu_GraphicalModel")
data = read.csv("competitor_pack/data/fire-data-2.csv", nrows = 100000)
#not DP, just using this for a test. user would specify domain normally.
domain = list()
for (i in 1:dim(data)[2]) {
  domain[[names(data)[i]]] = c(unique(data[,i]))
}

epsilon = 1
delta = 1e-6
sensitivity = 10
max_domain_size = 100000

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
  mutual_informations = vector(mode = "numeric")
  for (pair in x) {
    mutual_informations = c(mutual_informations, mutual_information(data, pair, type))
  }
  names(mutual_informations) = x
  return(mutual_informations)
}

mutual_information = function(data, pair, type = "L1") {
  these_elements = strsplit(pair, ';')[[1]]
  setA = as.set(strsplit(these_elements[1],',')[[1]])
  setB = as.set(strsplit(these_elements[2],',')[[1]])
  intersection = set_intersection(setA,setB)
  setA = set_complement(intersection, setA)
  setB = set_complement(intersection, setB)
  a = as.vector(setA, mode = "integer")
  b = as.vector(setB, mode = "integer")
  print(these_elements)
  print(a)
  print(b)
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


#trim whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

identity = function(a,...) {
  return(a)
}


select_queries = function(data, max_domain_size, sensitivity, domain, epsilon, delta) {
  n = dim(data)[1]
  d = length(names(data))
  Q = as.set(as.vector(1:d, mode = "character"))
  list_of_pairs = vector(mode = "character")
  #utility = vector()
  #start with just the mutual information of pairs of points
  for (i in 1:d) {
    for (j in 1:i){
      if (j < i && as.double(length(table(data[i])))*as.double(length(table(data[j])))< max_domain_size) {
        #utility = c(utility, compute_mutual_information(data, i, j))
        list_of_pairs = c(list_of_pairs, paste(as.character(c(i,j)), collapse = ";"))
      }
    }
  }
  #Now choose queries with exponential mechanism, we take pairs of cliques
  #cliques build up over time so we keep adding to the set that we're looking at
  for (i in 1:d) {
    #TO DO:Make this part work
    myMech = mechanismExponential()
    myMech$k = 1
    myMech$var.type = "character"
    myMech$n = n
    myMech$epsilon = epsilon/(2*d)
    myMech$delta = delta
    myMech$bins = list_of_pairs
    pair = myMech$evaluate(fun = mutual_information_vector, x = list_of_pairs, sens = sensitivity, postFun = identity, data = data, type = "L1")$release
    print(pair)
    S = as.set(strsplit(pair, ';')[[1]])
    u = mutual_information(data, pair, "KL") + rnorm(n=1, mean=0, sd = 2*d*epsilon/n)
    if (u > 4*d*epsilon/n) {
      #remove all pairings that only use elements in pair
      for (pairing in list_of_pairs) {
        these_elements = as.set(strsplit(pairing, ';')[[1]])
        if (set_is_subset(these_elements, S)) {
          list_of_pairs = list_of_pairs[list_of_pairs != pairing]
        }
      }
      new_clique = paste(as.character(S), collapse = ",")
      for (clique in Q) {
        index = as.vector(set_union(clique,S), mode = "integer")
        domain_size = 1
        for (attribute in index) {
          domain_size = domain_size * length(domain[attribute])
        }
        cliqueSet = as.set(strsplit(clique,',')[[1]])
        if (!set_is_subset(cliqueSet, S) && domain_size < max_domain_size) {
          
          new_pair = paste(as.character(c(new_clique, clique)),collapse = ";")
          list_of_pairs = c(list_of_pairs, new_pair)
        }
      }
      Q = set_union(Q, as.set(new_clique))
      print(Q)
      
    }
  }
  return(Q)
}