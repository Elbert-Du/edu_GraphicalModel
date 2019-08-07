library(reticulate)
library(rjson)
library("methods")
pd<-import("pandas")
#library("PSIlence")
#setwd("GraphicalModel/edu_GraphicalModel")
use_virtualenv("/home/home5/yz488/virtualenvs/PyEktelo/venv/bin/py3env",required = TRUE)

source("query_selection.R")
make_gm <- import_from_path("match3")
#When we actually use this, read these values from PSI interface



path_to_data <- "data/adult.csv"
#path_to_data <- "competitor_pack/data/fire-data.csv"
path_to_specs <- "data/specs.json"
path_to_domain <- "data/domain.json"
path_to_mapping <- "data/mapping.json"
domaininfo <- "data/adult-domain.json"

save <- vector(mode = "character")
for (i in 1:10) {
  save = c(save, paste("pg+new",as.character(i),".csv"))
}

bound = c(3:12)
howSplit = list()
for (i in 1:10) {
  howSplit[[i]] = c(1,i)
}

epsilon = 1
delta = 0.001


#data = read.csv(path_to_data, nrows = 100000)
data = pd$read_csv(path_to_data)
specs = fromJSON(file=path_to_specs)
domain = fromJSON(file=path_to_domain)
mapping = fromJSON(file=path_to_mapping)
domainininfo=fromJSON(file=domaininfo)

num_iters <- 10000
wrapper <- function(data, specs, domain, mapping, save, epsilon, delta, num_iters,domainininfo,query_selection = "new", bound = NULL, howSplit = c(1,1), encoded = FALSE,seed) {
  n = dim(data)[1]
  if (is.null(bound)) {
    bound = log2(n)
  }
  mech <- make_gm$Match3(data, specs, domain, mapping,delta, save, iters = num_iters, warmup = FALSE, is_encoded = encoded)
  epsilon1 = 2/3*howSplit[1]/sum(howSplit)*epsilon
  #In the future we will have to change epsilon2 to account for less privacy loss when user specifies queries
  epsilon2 = 2/3*howSplit[2]/sum(howSplit)*epsilon
 # new_data <- mech$shrink_domain(epsilon1, delta, bound)  # ? epsilon/3 it should be /2?
  if (query_selection == "new"){
   # compressed_data = new_data[1]
   # compressed_domain = new_data[2]
   # compressed_domain=compressed_domain[[1]]
      output = select_queries(data, n, domainininfo, epsilon = epsilon/2, delta = delta)
      queries = output$Q
      total_mutual_information = output$MI
      from_r = TRUE
  }
  else if(query_selection == "privbayes"){
      queries = mech$privbayes_query_selection(epsilon/2, seed)
      from_r = FALSE
  }
  else if(query_selection == "user"){
      queries = list()
  }
  #todo weights from user, network from user
  mech$measure(queries, epsilon/2, from_r)
  mech$postprocess()
  mech$write_output()
}

for (i in 1:10) {
    print(i)
    this_save = save[i]
    this_bound = bound[ceiling(i/10)]
    this_split = i
    seed=i*10
    wrapper(data, specs, domain, mapping, this_save, epsilon, delta, num_iters,domainininfo, query_selection = "new", bound = this_bound, howSplit = this_split, encoded = TRUE,seed)
}
                                        #wrapper(data, specs, domain, mapping, save, epsilon, delta, num_iters, query_selection = "new")

