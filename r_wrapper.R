library(reticulate)
library(rjson)
#library("PSIlence")
#setwd("GraphicalModel/edu_GraphicalModel")
source("query_selection.R")
make_gm <- import("match3")

#When we actually use this, read these values from PSI interface
path_to_data <- "competitor_pack/data/fire-data-2.csv"
path_to_specs <- "competitor_pack/data/fire-data-specs.json"
path_to_domain <- "domain.json"
<<<<<<< HEAD
path_to_mapping <- "competitor_pack/data/fire-data-specs-mapping.json"
out_file_name <- "out.csv"
=======
path_to_mapping <- "competitor_pack/data/fire-data-specs-mapping.csv"
save <- "output.csv" #todo where to save?
>>>>>>> 9aa3deb62e9b60ea4da805cbdf0e9b3b509ddf14

epslion = 1
delta = 1e-6


data = read.csv(path_to_data, nrows = 100000)
specs = fromJSON(paste(readLines(path_to_specs)))
domain = fromJSON(paste(readLines(path_to_domain)))
mapping = fromJSON(paste(readLines(path_to_mapping)))


num_iters <- 1000
wrapper <- function(data, specs, domain, mapping, epsilon, delta, num_iters, out) {
  n = dim(data)[1]
  mech <- make_gm$Match3(data, specs, domain, mapping, iters = num_iters, warmup = TRUE)
  new_data <- mech$shrink_domain(data, epsilon/3, delta)
  compressed_data = first_measurements[1]
  compressed_domain = first_measurements[2]
  queries = select_queries(compressed_data, n, sensitivity = 1, domain = compressed_domain, epsilon = epsilon/3, delta = delta)
  mech$measure(queries)
  mech$postprocess()
<<<<<<< HEAD
  mech$write_output(out)
=======
  mech$write_output(save)
>>>>>>> 9aa3deb62e9b60ea4da805cbdf0e9b3b509ddf14
}
