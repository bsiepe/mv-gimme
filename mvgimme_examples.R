library(tidyverse)
data <- mvgimme::simData

# multicore
mv_res <- multiverse.gimme(data = data,
                 groupcutoffs = c(.5, .6, .7),
                 subcutoffs = c(.5),
                 n.cores = 3,
                 subgroup = TRUE)

ref_fit <- gimme::gimme(data = data, 
                 groupcutoff = .75,
                 subcutoff   = .75,
                 subgroup = TRUE)

# saveRDS(mv_res, "mv_res.RDS")
# saveRDS(ref_fit, "ref_fit.RDS")

# single core
before_gimme <- Sys.time()
mv_res <- multiverse.gimme(data = data,
                           n.cores = 1)
after_gimme <- Sys.time() - before_gimme



test_combs <- expand.grid(
  groupcutoffs = c(.5, .6, .7),
  subcutoffs = .51,
  rmsea.cuts = .05,
  srmr.cuts = .05,
  nnfi.cuts = .95,
  cfi.cuts = .95,
  n.excellent = 2
)


test_all <- multiverse.compare(l_res = mv_res, 
                               ref_model = ref_fit)

test_group <- multiverse.compare.group(l_res = mv_res, 
                                 ref_model = ref_fit)


test_subgroup <- multiverse.compare.subgroup(l_res = mv_res, 
                                             ref_model = ref_fit)

test_individual <- multiverse.compare.individual(l_res = mv_res, 
                                                 ref_model = ref_fit)

