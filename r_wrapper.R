library(reticulate)
library(rjson)
library("methods")
pd<-import("pandas")
#library("PSIlence")
#setwd("GraphicalModel/edu_GraphicalModel")
use_virtualenv("/home/home5/yz488/virtualenvs/PyEktelo/venv/bin/py3env",required = TRUE)

source("query_selection.R")
make_gm <- import("match3")
#When we actually use this, read these values from PSI interface



#path_to_data <- "competitor_pack/data/fire-data-2.csv"
path_to_data <- "competitor_pack/data/fire-data.csv"
path_to_specs <- "competitor_pack/data/fire-data-specs.json"
path_to_domain <- "domain.json"
path_to_mapping <- "competitor_pack/data/fire-data-specs-mapping.json"

save <- "output.csv" #todo where to save?

epslion = 1
delta = 1e-6


#data = read.csv(path_to_data, nrows = 100000)
data = pd$read_csv(path_to_data, nrows = 1000)
specs = fromJSON(file=path_to_specs)
domain = fromJSON(file=path_to_domain)
mapping = fromJSON(file=path_to_mapping)


num_iters <- 1000
wrapper <- function(data, specs, domain, mapping, save, epsilon, delta, num_iters,query_selection = "new") {
  n = dim(data)[1]
  mech <- make_gm$Match3(data, specs, domain, mapping, save, iters = num_iters, warmup = FALSE)
  new_data <- mech$shrink_domain(epsilon/3, delta)  # ? epsilon/3 it should be /2?
  if (query_selection == "new"){
      compressed_data = new_data[1]
      compressed_domain = new_data[2]
      compressed_domain=compressed_domain[[1]]
      for(name in names(compressed_domain)){
          compressed_domain[[name]]=c(0,seq(as.numeric(compressed_domain[[name]])-1))
      }
      queries = select_queries(data, n, sensitivity = 1, domain = compressed_domain, epsilon = epsilon/3, delta = delta)#change NAN in compressed_data
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
  mech$measure(queries, from_r)
  mech$postprocess()
  mech$write_output()
}
wrapper(data, specs, domain, mapping, save, epsilon, delta, num_iters, query_selection = "privbayes" )

