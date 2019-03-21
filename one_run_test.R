source('source_data.R')




batch_folder = "over_100"
solving_method = "lsoda"


all_good_pars <- get(load(file = paste0(getwd(), "/114_fits.Rdata")))
all_good_pars <- all_good_pars[[1]]

run_model_for_longer <- function(x) {
  
  # cleaning up the parameters list such that each element only has a length of 1 (those parameters which exist as vectors, e.g. PrEP efficacy)
  x = lapply(x[rownames(ranges)], function(y) {
    if(length(y) == 9)
      y = y[1] else y
  })
  
  combined_ranges = cbind(unlist(x[rownames(ranges)]), unlist(x[rownames(ranges)]))
  
  res_after_prep = cotonou::run_model_with_fit(number_simulations = 1,
                                               par_seq = par_seq, condom_seq = condom_seq,
                                               groups_seq = groups_seq, years_seq = years_seq,
                                               best_set = best_set,
                                               time = c(seq(1986, 2035)),
                                               ranges = combined_ranges, outputs = CEA_outputs,
                                               prev_points = prev_points,
                                               frac_N_discard_points = frac_N_discard_points,
                                               Ntot_data_points = Ntot_data_points,
                                               ART_data_points = ART_data_points,
                                               PrEP_fitting = PrEP_fitting)
  # return(list(res_after_prep[[1]][[1]], res_after_prep[[3]][[1]]))
  return(res_after_prep[[3]][[1]])
}

res = lapply(all_good_pars, run_model_for_longer)
