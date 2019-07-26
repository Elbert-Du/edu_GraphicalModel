library(reticulate)
library(rjson)
#library("PSIlence")
setwd("GraphicalModel/edu_GraphicalModel")
make_gm <- import("match3")

#When we actually use this, read these values from PSI interface
path_to_data <- "competitor_pack/data/fire-data-2.csv"
path_to_specs <- "competitor_pack/data/fire-data-specs.csv"
path_to_domain <- "domain.json"
path_to_mapping <- "competitor_pack/data/fire-data-specs-mapping.csv"

epslion = 1
delta = 1e-6


data = read.csv(path_to_data)
specs = fromJSON(path_to_specs)
domain = fromJSON(path_to_domain)
mapping = fromJSON(path_to_mapping)


wrapper <- function(data, specs, domain, mapping, epsilon, delta) {
  n = dim(data)[1]
  first_measurements <- make_gm$shrink_domain(data, epsilon/3, delta)
  measurements = first_measurements[1]
  supports = first_measurements[2]
  new_data = make_gm$transform_data(data, domain, supports)
  compressed_data = new_data[1]
  compressed_domain = new_data[2]
  queries = select_queries(compressed_data, n, sensitivity = 1, domain = compressed_domain, epsilon = epsilon/3, delta = delta)
  mech <- make_gm$Match3(data, specs, domain, mapping, measurements, supports, queries)
  mech$measure()
  mech$postprocess()
  mech$write_output()
}
