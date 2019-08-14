library(reticulate)
library(rjson)
#setwd("GraphicalModel/edu_GraphicalModel")
n.bins <- fromJSON(file='data/adult-domain.json')
bins <- fromJSON(file = 'data/domain.json')
var.type <- fromJSON(file = "data/specs.json")
data = data = read.csv("data/adult.csv", check.names = FALSE)
source("query_selection.R")
pd<-import("pandas")

for (att in names(var.type)) {
  var.type[[att]] = var.type[[att]][["type"]]
}
var.type = as.vector(unname(var.type), mode = "character")

mechanismSynthetic <- setRefClass(
  Class = 'mechanismSynthetic',
  contains = 'mechanism',
  fields = list(
    "mech" = "ANY",
    "epsilon1" = "numeric",
    "epsilon2" = "numeric",
    "epsilon3" = "numeric",
    "usefulness" = "numeric"
  )
)
make_gm <- import_from_path("match3")

make_specs = function(domain, var.type, num_bins = NULL) {
  specs = list()
  i=1
  for (att in names(domain)) {
    specs[[att]] = list()
    specs[[att]][["type"]] = var.type[[i]]
    if (var.type[[i]] == "character") {
      specs[[att]][["num_bins"]] = length(domain[[att]])
    } else {
      specs[[att]][["lower_bound"]] = min(domain[[att]])
      specs[[att]][["upper_bound"]] = max(domain[[att]])
      specs[[att]][["num_bins"]] = num_bins[[att]]
    }
    i = i+1
  }
  return(specs)
}

make_mapping = function(domain, specs, num_bins = NULL) {
  mapping = list()
  for (att in names(specs)) {
    mapping[[att]] = list()
    if (specs[[att]][["type"]] == "character") {
      for (i in 1:specs[[att]][["num_bins"]]) {
        mapping[[att]][[as.character(i)]] = domain[[att]][i]
      }
    } else {
      a = specs[[att]][["lower bound"]]
      b = specs[[att]][["upper bound"]]
      if (!is.null(num_bins)) {
        for (i in 1:num_bins[[att]]) {
          mapping[[att]][[as.character(i)]] = i/num_bins[[att]]*(b-a)+a
        }
      }
      
    }
  }
  return(mapping)
}

mechanismSynthetic$methods(
  initMech = function(data) {
    domain = .self$bins
    specs = make_specs(domain, .self$var.type)
    mapping = make_mapping(domain, specs, .self$n.bins)
    my.mech = make_gm$Match3(data, specs, domain, mapping, .self$delta)
    .self$mech = my.mech
  }
)

mechanismSynthetic$methods(
  shrinkDomain = function(data) {
    num_attributes = length(names(data))
    bound = 2
    new_data <- mech$shrink_domain(.self$epsilon1, .self$delta, bound)
    compressed_domain = new_data[[2]]
    compressed_data = py_to_r(new_data[[1]])
    rownames(compressed_data) = c()
    #make sure we don't have any NAs by mapping them to the attribute which is all of the low counts
    for (att in names(compressed_data)) {
      na_indices = which(is.na(compressed_data[[att]]))
      compressed_data[[att]][na_indices] = compressed_domain[[att]]-1
    }
    .self$n.bins = compressed_domain
    return(compressed_data)
  })

mechanismSynthetic$methods(
  selectQueries = function(data) {
    queries = select_queries(data, .self$n, domain = .self$n.bins, epsilon = .self$epsilon2, delta = .self$delta)#change NAN in compressed_data
    return(queries)
  })

mechanismSynthetic$methods(
  release = function(queries) {
    .self$mech$measure(queries, from_r = TRUE)
    .self$mech$postprocess()
    .self$mech$write_output()
    return(.self$mech$measurements)
  })

myMech = mechanismSynthetic()
myMech$n = dim(data)[1]
myMech$var.type = var.type
myMech$n.bins = n.bins
myMech$bins = bins
myMech$delta = 1e-6
myMech$initMech(data)
myMech$epsilon = 1
myMech$epsilon1 = 2/3
myMech$epsilon2 = 1/3
myMech$usefulness = 2
if (myMech$n > 100000/myMech$epsilon1) {
  compressed_data = myMech$shrinkDomain(data)
  
} else {
  sigma = make_gm$moments_calibration(1e-30, 1, myMech$epsilon1, myMech$delta)
  myMech$mech$sigma = sigma
  queries = myMech$selectQueries(data)
}

measurements = myMech$release(queries)
