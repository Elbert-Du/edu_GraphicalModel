library(reticulate)
library(rjson)
library("methods")
pd<-import("pandas")
#library("PSIlence")
#setwd("GraphicalModel/edu_GraphicalModel")
#use_virtualenv("/home/home5/yz488/virtualenvs/PyEktelo/venv/bin/py3env",required = TRUE)

source("query_selection.R")
make_gm <- import_from_path("match3")
#When we actually use this, read these values from PSI interface



path_to_data <- "data/adult.csv"
#path_to_data <- "competitor_pack/data/fire-data.csv"
path_to_specs <- "data/specs.json"
path_to_domain <- "data/domain.json"
path_to_mapping <- "data/mapping.json"

save <- vector(mode = "character")
for (i in 1:100) {
  save = c(save, paste(as.character(i), ".csv"))
}

bound = c(2:11)
howSplit = list()
for (i in 1:10) {
  howSplit[[i]] = c(1,i)
}

epsilon = 1
delta = 1e-6


#data = read.csv(path_to_data, nrows = 100000)
data = read.csv(path_to_data, nrows = 10000, check.names = FALSE)
specs = fromJSON(file=path_to_specs)
domain = fromJSON(file=path_to_domain)
mapping = fromJSON(file=path_to_mapping)


num_iters <- 1000
wrapper <- function(data, specs, domain, mapping, save, epsilon, delta, num_iters,query_selection = "new", bound = NULL, howSplit = c(1,1), encoded = FALSE) {
  n = dim(data)[1]
  total_mutual_information = 0
  if (is.null(bound)) {
    bound = log2(n)
  }
  mech <- make_gm$Match3(data, specs, domain, mapping, save, iters = num_iters, warmup = FALSE, is_encoded = encoded)
  epsilon1 = 2/3*howSplit[1]/sum(howSplit)*epsilon
  #In the future we will have to change epsilon2 to account for less privacy loss when user specifies queries
  epsilon2 = 2/3*howSplit[2]/sum(howSplit)*epsilon
  new_data <- mech$shrink_domain(epsilon1, delta, bound)  # ? epsilon/3 it should be /2?
  if (query_selection == "new"){
    compressed_domain = new_data[[2]]
    compressed_data = py_to_r(new_data[[1]])
    rownames(compressed_data) = c()
    for (att in names(compressed_data)) {
      na_indices = which(is.na(compressed_data[[att]]))
      compressed_data[[att]][na_indices] = compressed_domain[[att]]-1
    }
    output = select_queries(data, n, domain = compressed_domain, epsilon = epsilon/3, delta = delta)#change NAN in compressed_data
    queries = output$Q
    total_mutual_information = output$MI
    from_r = TRUE
  }
  else if(query_selection == "privbayes"){
    queries = mech$privbayes_query_selection(epsilon/2, seed=0)
    from_r = FALSE
  }
  else if(query_selection == "user"){
    queries = list()
  }
  #todo weights from user, network from user
  mech$measure(queries, epsilon2, from_r)
  mech$postprocess()
  mech$write_output()
  return(total_mutual_information)
}

#for (i in 1:100) {
#  print(i)
#  this_save = save[i]
#  this_bound = bound[ceiling(i/10)]
#  this_split = howSplit[[i%%10]]
#  wrapper(data, specs, domain, mapping, this_save, epsilon, delta, num_iters, query_selection = "new", bound = this_bound, howSplit = this_split, encoded = TRUE)
#}
#wrapper(data, specs, domain, mapping, save, epsilon, delta, num_iters, query_selection = "new")
bound = c(2:100)
these_domains = list()
these_data = list()
for (n in (1:10)*10000) {
  these_data[[n/10000]] = list()
  these_domains[[n/10000]] = list()
  for (m in 1:99) {
    my_bound = bound[m]
    bigger_data = data[sample(1:dim(data)[1], size = n, replace = TRUE),]
    encoded = TRUE
    mech <- make_gm$Match3(bigger_data, specs, domain, mapping, save, iters = num_iters, warmup = FALSE, is_encoded = encoded)
    epsilon1 = 2/3*howSplit[[1]][1]/sum(howSplit[[1]])*epsilon
    new_data <- mech$shrink_domain(epsilon1, delta, my_bound)
    compressed_data = new_data[[1]]
    these_data[[n/10000]][[m]] = py_to_r(compressed_data)
    rownames(these_data[[n/10000]][[m]]) = c()
    compressed_domain = new_data[[2]]
    these_domains[[n/10000]][[m]] = compressed_domain
  }
  
  #queries = select_queries(bigger_data, n, domain = compressed_domain, epsilon = epsilon/3, delta = delta)
}

bad_indices = matrix(nrow = 10, ncol = 99)
prop_dropped = matrix(nrow = 10, ncol = 99)
for (n in 1:10) {
  for (j in 1:99) {
    dropped = 0
    for (att in names(data)) {
      m = these_domains[[n]][[j]][[att]]-1
      num_bad = length(which(these_data[[n]][[j]][[att]] == m))
      print(c(num_bad, dropped))
      if (m+1 < length(domain[[att]])) {
        dropped = dropped + num_bad/(n*10000*length(names(data)))
        if (num_bad > n*10000/length(names(data))) {
          bad_indices[n,j] = 1
        }
      }
    }
    prop_dropped[n,j] = dropped
  }
}

domain_reduction = matrix(nrow = 10, ncol = 99)
for (n in 1:10) {
  for (m in 1:99) {
    num_removed = 1
    for (att in names(data)) {
      num_removed = num_removed * length(domain[[att]]) / these_domains[[n]][[m]][[att]]
    }
    domain_reduction[n,m] = log(num_removed)
  }
}
x = vector()
y = vector()
z1 = vector()
z2 = vector()
for (i in 1:10) {
  for (j in 1:99) {
    x = c(x,i)
    y = c(y,j)
    z1 = c(z1, domain_reduction[i,j])
    z2 = c(z2, prop_dropped[i,j])
  }
}
library(scatterplot3d)
scatterplot3d(x,y,z1, xlab = "n/10000", ylab = "bound-2", zlab = "log domain reduction", main = "domain reduction")
dev.copy(png, 'bound_plot_domain_reduction.png')
dev.off()

scatterplot3d(x,y,z2, xlab = "n/10000", ylab = "bound-2", zlab = "proportion data lost", main = "data loss")
dev.copy(png, 'bound_plot_data_loss.png')
dev.off()