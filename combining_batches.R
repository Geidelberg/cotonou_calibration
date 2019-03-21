#    devtools::install_github("geidelberg/cotonou")

#     setwd("Q:/cotonou_cluster")

# setwd("C:\\Users\\eg1012\\Desktop")

source('Q:/cotonou_cluster/source_data.R')



batch_folder = "over_100"



# batch_folder = "final_111"

solving_method = "lsoda"


# fits
# res_best_runs <- get(load(file = "C:\\Users\\eg1012\\Desktop/final_111/res_best_runs_final_111.Rdata"))
time <- c(seq(1986, 2035))


# everything below here commented out together can be used to sort through the varioous fits from different batches

if(!exists("res_best_runs")) {
  file_list = list.files(path = paste0("Q:/cotonou_cluster/runs/", batch_folder))
  file_list_res = file_list[grep("res",file_list)]
  file_list_pars = file_list[grep("pars",file_list)]

  all_good_runs = list()
  all_good_pars = list()
  for(i in 1:length(file_list_res))
  {

    # print(length(batch))
    batch = get(load(paste0("Q:/cotonou_cluster/runs/", batch_folder, "/",file_list_res[i])))
    all_good_runs[(length(all_good_runs) + 1):(length(all_good_runs) + length(batch))] <- batch

    pars = get(load(paste0("Q:/cotonou_cluster/runs/", batch_folder, "/",file_list_pars[i])))
    all_good_pars[(length(all_good_pars) + 1):(length(all_good_pars) + length(pars))] <- pars

    print(file_list_res[i])
    print(length(batch))

  }



  # removing the one run which doesn't fit GPM
  all_good_pars <- all_good_pars[unlist(lapply(res_best_runs, function(x) {
    x$prev_men[which(time == 2011)] < 2.854672
  }))]

  # removing the one run which doesn't fit GPM
  all_good_runs <- all_good_runs[unlist(lapply(all_good_runs, function(x) {
    x$prev_men[which(time == 2011)] < 2.854672
  }))]



  # removing the runs that don't fit % women FSW based on lower bound of mapping on numerator
  all_good_pars <- all_good_pars[which(unlist(lapply(all_good_runs, function(x) {
    return(all(x$frac_N[,1] > 0.001944823*0.515666224))})))]

  all_good_runs <- all_good_runs[which(unlist(lapply(all_good_runs, function(x) {
    return(all(x$frac_N[,1] > 0.001944823*0.515666224))})))]
  
  
  
  # # removing parameters not used in the model
  # all_good_pars <- lapply(all_good_pars, function(x) {
  #   x <- x[!names(x) %in% c(unlist(strsplit("c_comm_1993_LowFSW, c_comm_1995_ProFSW, c_comm_2002_Client, c_noncomm_1985_Client, c_noncomm_1985_LowFSW, c_noncomm_1985_ProFSW, fc_y_comm_1985_LowFSW_Client, fc_y_comm_1998_ProFSW_Client, fc_y_comm_2008_ProFSW_Client, fc_y_comm_2015_LowFSW_Client, fc_y_noncomm_1985_LowFSW_Client, fc_y_noncomm_2015_LowFSW_Client, n_y_noncomm_1998_GPF_GPM, n_y_noncomm_2011_GPF_GPM, n_y_noncomm_1998_GPF_Client, n_y_noncomm_2011_GPF_Client", ", ")))]
  #   return(x)
  # })


  length(all_good_runs)
  res_best_runs = all_good_runs

} else {time <- c(seq(1986, 2035))}




# removing the ones that don't fit total number on ART in 2017 for all
all_good_pars <- all_good_pars[unlist(lapply(res_best_runs, function(x) {
  GP_ART_2017 = sum(x$HIV_positive_On_ART[which(time == 2017), c(1:8)])
  return(GP_ART_2017 > 2369+6155 && GP_ART_2017 < 5223+13050)
}))]

res_best_runs <- res_best_runs[unlist(lapply(res_best_runs, function(x) {
  GP_ART_2017 = sum(x$HIV_positive_On_ART[which(time == 2017), c(1:8)])
  return(GP_ART_2017 > 2369+6155 && GP_ART_2017 < 5223+13050)
}))]



length(res_best_runs)
length(all_good_pars)





#  ____  _   _ _   _ _   _ ___ _   _  ____   _____ ___ _       ____   ___ _________
# |  _ \| | | | \ | | \ | |_ _| \ | |/ ___| |_   _|_ _| |     |___ \ / _ \___ / ___|
# | |_) | | | |  \| |  \| || ||  \| | |  _    | |  | || |       __) | | | ||_ \___ \
# |  _ <| |_| | |\  | |\  || || |\  | |_| |   | |  | || |___   / __/| |_| |__) |__) |
# |_| \_\\___/|_| \_|_| \_|___|_| \_|\____|   |_| |___|_____| |_____|\___/____/____/
#

## !!! redefining the ART eligibility and putting it to 0 for all CD4 count...
# this is because we're assuming eligibility won't increase unless TasP happens...

all_good_pars <- lapply(all_good_pars, function(x) {
  # x$above_500_by_group <- c(0,0,0,0,0,0,0,0,0)

  FSW_eligible = 0
  GP_eligible = 1

  x$TasP_testing <- 0

  x$PrEPOnOff =  0
  x$rho_intervention_y = x$rho_intervention_y * 0
  x$tau_intervention_y = x$tau_intervention_y * 0
  x$ART_eligible_CD4_above_500_y <- c(0,0,0,0,1)

  return(x)
})

all_good_pars <- all_good_pars[!duplicated(unlist(lapply(all_good_pars, function(x) x$rate_leave_client)))]



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


#!!!!!!!!!
# 24012019


# length(unique(unlist(lapply(all_good_pars, function(x) x$rate_leave_client))))


# res_both <- res[!duplicated(unlist(lapply(all_good_pars, function(x) x$rate_leave_client)))]

# check for duplicates!


#
res_best_runs <- res
time <- c(seq(1986, 2035))

# above




# dir.create(path = paste0("C:\\Users\\eg1012\\Google Drive\\Imperial\\Cotonou\\Models\\Transmission Model\\Outputs",
#                          "/", batch_folder))




### ! below now sure if wanna do still!




remove_run = function(x) {

  # x$meet_criteria = ifelse(x$lambda_sum_0[which(time == 2012), 1] < 0.05,
  #                          ifelse(x$prev[which(time == 2012), 5] < 5.160244,
  #                                 ifelse(x$prev[which(time == 2011), 6] < 2.854672,
  #                                        ifelse(x$prev[which(time == 1993), 1] < 58.48,
  #                                               ifelse(x$prev[which(time == 1993), 1] > 48.0200000,
  #                                               ifelse(x$prev[which(time == 2015), 1] > 15.7100000,
  #                                        ifelse(x$N[which(time == 2012), 1] < 1300, 1, 0), 0), 0), 0), 0), 0), 0)

  x$meet_criteria = ifelse(x$lambda_sum_0[which(time == 2015), 1] < 0.03,

                           # prevalence
                           ifelse(all(x$prev[, 1]  < 65),
                                  ifelse(x$prev[which(time == 2015), 1]  < 22,
                                         ifelse(x$prev[which(time == 2015), 1]  > 13.79,
                                                ifelse(x$prev[which(time == 1993), 1]  < 58.48,
                                                       ifelse(x$prev[which(time == 1993), 1]  > 48.02,
                                                              ifelse(x$prev[which(time == 2005), 5]  < 10.5, # client
                                                                     ifelse(x$prev[which(time == 2005), 5]  > 4.29,  # client
                                                                            ifelse(x$prev[which(time == 2011), 3]  > 0, # GPF
                                                                                   ifelse(x$prev_women[which(time == 2011)]  < 3.529629, # GPF
                                                                                          ifelse(x$prev_men[which(time == 2011)]  < 2.854672, # GPM
                                                                                                 # ifelse(x$prev[which(time == 2011), 6]  > 2.854672, # GPM

                                                                                                 ifelse(x$N[which(time == 2012), 1]  < 1391,
                                                                                                        ifelse(x$N[which(time == 2012), 1]  > 889,
                                                                                                               ifelse(all(x$HIV_positive_On_ART[which(time == 2015), 1]  < 56),
                                                                                                                      ifelse(all(x$HIV_positive_On_ART[which(time == 2015), 1]  > 42),
                                                                                                                             ifelse(sum(x$HIV_positive_On_ART[which(time == 2017), c(1,2,3,4)])/
                                                                                                                                      sum(x$HIV_positive_On_ART[which(time == 2017), c(5,6)]) < 3,
                                                                                                                                    ifelse(sum(x$HIV_positive_On_ART[which(time == 2017), c(1,2,3,4)])/
                                                                                                                                             sum(x$HIV_positive_On_ART[which(time == 2017), c(5,6)]) > 1.4,
                                                                                                                                           ifelse(rowSums(data.frame(x$I32 + x$I33 + x$I34 + x$I35)[,c(1:8)])[which(time == 2017)]< (5223+13050),
                                                                                                                                                  ifelse(rowSums(data.frame(x$I32 + x$I33 + x$I34 + x$I35)[,c(1:8)])[which(time == 2017)]> (2369+6155),





                                                                                                                                                         1, 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0), 0)

  return(x$meet_criteria)
}


all_good_pars = all_good_pars[unlist(lapply(res_best_runs, remove_run)) == 1]

res_best_runs_screened = res_best_runs[unlist(lapply(res_best_runs, remove_run)) == 1]


length(res_best_runs_screened)

res_best_runs <- res_best_runs_screened




length(all_good_pars)
length(res_best_runs)


save(res_best_runs, file = "res_best_runs_final_111.Rdata")
save(all_good_pars, file = "all_good_pars_final_111.Rdata")


# save(res_best_runs, file = paste0("runs/","final_111","/","best_runs_day_", "final_114","_", "combined", ".Rdata"))
# 
# save(all_good_pars, file = paste0("runs/","final_111","/","best_prams_day_", "final_114","_", "combined", ".Rdata"))
# save(all_good_pars, file = paste0("Q:/CEAcotonou/","final_111","/","best_prams_day_", "final_114","_", "combined", ".Rdata"))
# 
# 





# ignore these ######################################
# frac_ProFSW_F = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N), function(x) (x[,1]/(x[,1] + x[,2] + x[,3] + x[,4] + x[,7])))), 2, cotonou::quantile_95)))
# frac_ProFSW_F = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N/x$frac_F), function(x) x[,1])), 2, cotonou::quantile_95)))








annual_client_volume_pro_FSW = data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x)  {x$c_comm[,1]}))))
colnames(annual_client_volume_pro_FSW) = c("time", as.character(seq(1, (length(annual_client_volume_pro_FSW[1,])-1))))
annual_client_volume_pro_FSW_melted = reshape2::melt(annual_client_volume_pro_FSW, id.vars = c("time"))



N_Pro_FSW = data.frame(time, t(do.call(rbind, lapply(lapply(res_best_runs, function(x)  {x$N}), function(x) {return(x[,c(1)])}))))
colnames(N_Pro_FSW) = c("time", as.character(seq(1, (length(N_Pro_FSW[1,])-1))))
N_Pro_FSW_melted = reshape2::melt(N_Pro_FSW, id.vars = c("time"))
colnames(N_Pro_FSW_melted) = c("time", "run", "value")
N_Low_FSW = data.frame(time, t(do.call(rbind, lapply(lapply(res_best_runs, function(x)  {x$N}), function(x) {return(x[,c(2)])}))))
colnames(N_Low_FSW) = c("time", as.character(seq(1, (length(N_Low_FSW[1,])-1))))
N_Low_FSW_melted = reshape2::melt(N_Low_FSW, id.vars = c("time"))
colnames(N_Low_FSW_melted) = c("time", "run", "value")
N_Client = data.frame(time, t(do.call(rbind, lapply(lapply(res_best_runs, function(x)  {x$N}), function(x) {return(x[,c(5)])}))))
colnames(N_Client) = c("time", as.character(seq(1, (length(N_Client[1,])-1))))
N_Client_melted = reshape2::melt(N_Client, id.vars = c("time"))
colnames(N_Client_melted) = c("time", "run", "value")



Fraction_F =   do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$N}), function(x) {return(x[,c(1, 2, 3, 4, 7)])}), rowSums)) / 
  (do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$N}), function(x) {return(x[,c(5, 6)])}), rowSums)) + do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$N}), function(x) {return(x[,c(1, 2, 3, 4, 7)])}), rowSums)))


lambda_sum_0_ProFSW = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,1]))))
lambda_sum_0_LowFSW = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,2]))))
lambda_sum_0_GPF = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,3]))))
lambda_sum_0_FormerFSW = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,4]))))
lambda_sum_0_Client = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,5]))))
lambda_sum_0_GPM = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,6]))))
lambda_sum_0_Virgin_Female = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,7]))))
lambda_sum_0_Virgin_Male = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,8]))))
lambda_sum_0_Former_FSW_Outside = data.frame(t(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$lambda_sum_0), function(x) x[,9]))))




lambda_sum_0_indiv = rbind(lambda_sum_0_ProFSW, lambda_sum_0_LowFSW, lambda_sum_0_GPF, 
                           lambda_sum_0_FormerFSW, lambda_sum_0_Client, lambda_sum_0_GPM)


lambda_sum_0_indiv = data.frame(time, rep(c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients", "GPM"), each = length(time)), lambda_sum_0_indiv)


colnames(lambda_sum_0_indiv) = c("time", "variable", as.character(seq(1, length(lambda_sum_0_ProFSW[1,]))))
lambda_sum_0_indiv_melted = reshape2::melt(lambda_sum_0_indiv, id.vars = c("time", "variable"))
colnames(lambda_sum_0_indiv_melted) = c("time", "variable", "run", "value")
lambda_sum_0_indiv_melted$variable = factor(lambda_sum_0_indiv_melted$variable, levels = c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients", "GPM"))






frac_ProFSW = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,1])), 2, cotonou::quantile_95)))
frac_LowFSW = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,2])), 2, cotonou::quantile_95)))
frac_GPF = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,3])), 2, cotonou::quantile_95)))
frac_FormerFSW = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,4])), 2, cotonou::quantile_95)))
frac_Client = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,5])), 2, cotonou::quantile_95)))
frac_GPM = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,6])), 2, cotonou::quantile_95)))
frac_Virgin_Female = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,7])), 2, cotonou::quantile_95)))
frac_Virgin_Male = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,8])), 2, cotonou::quantile_95)))
frac_Former_FSW_Outside = data.frame(time, t(apply(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$frac_N*100), function(x) x[,9])), 2, cotonou::quantile_95)))
frac_Active_FSW = data.frame(time, t(apply(do.call(rbind, lapply(res_best_runs, function(x) {100*(x$frac_N[,1] + x$frac_N[,2])})), 2, cotonou::quantile_95)))
Ratio_Low_Pro = data.frame(time, t(apply(do.call(rbind, lapply(res_best_runs, function(x) {x$frac_N[,2]/ x$frac_N[,1]})), 2, cotonou::quantile_95)))

frac = rbind(frac_ProFSW, frac_LowFSW, frac_GPF, frac_FormerFSW, frac_Client, frac_GPM, frac_Virgin_Female, frac_Virgin_Male, frac_Former_FSW_Outside, frac_Active_FSW, Ratio_Low_Pro)
frac = data.frame(frac, group = rep(c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients", "GPM", "Virgin female", "Virgin male", "Former FSW outside Cotonou", "Active FSW", "Low Pro Ratio"), each = length(time)))
colnames(frac) = c("time", "Lower", "Median", "Upper", "variable")
frac$variable = factor(frac$variable, levels = c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients", "GPM", "Virgin female", "Virgin male", "Former FSW outside Cotonou", "Active FSW", "Low Pro Ratio"))

prev_FSW = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$prev_FSW)), 2, cotonou::quantile_95))
prev_LowFSW = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$prev_LowFSW)), 2, cotonou::quantile_95))
prev_client = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$prev_client)), 2, cotonou::quantile_95))
prev_women = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$prev_women)), 2, cotonou::quantile_95))
prev_men = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$prev_men)), 2, cotonou::quantile_95))
prev = rbind(prev_FSW, prev_LowFSW, prev_client, prev_women, prev_men)
prev = data.frame(time, prev, rep(c("Pro FSW", "Low-level FSW", "Clients", "Women", "Men"), each = length(time)))
colnames(prev) = c("time", "Lower", "Median", "Upper", "variable")
prev$variable = factor(prev$variable, levels = c("Pro FSW", "Low-level FSW", "Clients", "Women", "Men"))


prev_FSW_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$prev_FSW)))
prev_LowFSW_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$prev_LowFSW)))
prev_client_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$prev_client)))
prev_women_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$prev_women)))
prev_men_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$prev_men)))
prev_indiv = rbind(prev_FSW_indiv, prev_LowFSW_indiv, prev_client_indiv, prev_women_indiv, prev_men_indiv)

prev_indiv = data.frame(time, rep(c("Pro FSW", "Low-level FSW", "Clients", "Women", "Men"), each = length(time)), prev_indiv)






colnames(prev_indiv) = c("time", "variable", as.character(seq(1, length(prev_FSW_indiv[1,]))))

prev_indiv_melted = reshape2::melt(prev_indiv, id.vars = c("time", "variable"))

colnames(prev_indiv_melted) = c("time", "variable", "run", "value")

prev_indiv_melted$variable = factor(prev_indiv_melted$variable, levels = c("Pro FSW", "Low-level FSW", "Clients", "Women", "Men"))




Ntot = data.frame(time, t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$Ntot)), 2, cotonou::quantile_95)))
colnames(Ntot) = c("time", "Lower", "Median", "Upper")


ART_coverage_women = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_women)), 2, cotonou::quantile_95))
ART_coverage_men = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_men)), 2, cotonou::quantile_95))
ART_coverage_FSW = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_FSW)), 2, cotonou::quantile_95))
ART_coverage_all = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_all)), 2, cotonou::quantile_95))
ART_coverage = rbind(ART_coverage_women, ART_coverage_men, ART_coverage_FSW, ART_coverage_all)
ART_coverage = data.frame(time, ART_coverage, rep(c("Women", "Men", "Pro FSW", "All"), each = length(time)))
colnames(ART_coverage) = c("time", "Lower", "Median", "Upper", "variable")
ART_coverage$variable = factor(ART_coverage$variable, levels = c("Pro FSW", "Women", "Men", "All"))
ART_coverage = ART_coverage[ART_coverage$variable == "All" | ART_coverage$variable == "Pro FSW",]



ART_FSW_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_FSW)))
ART_women_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_women)))
ART_men_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_men)))

ART_all_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$ART_coverage_all)))


ART_indiv = rbind(ART_FSW_indiv, ART_all_indiv)

ART_indiv = data.frame(time, rep(c("Pro FSW", "All"), each = length(time)), ART_indiv)


colnames(ART_indiv) = c("time", "variable", as.character(seq(1, length(ART_FSW_indiv[1,]))))

ART_indiv_melted = reshape2::melt(ART_indiv, id.vars = c("time", "variable"))

colnames(ART_indiv_melted) = c("time", "variable", "run", "value")

ART_indiv_melted$variable = factor(ART_indiv_melted$variable, levels = c("Pro FSW", "All"))


# N of FSW on ART
N_ART_FSW = t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$HIV_positive_On_ART[,1])), 2, cotonou::quantile_95))
N_ART_FSW = data.frame(time, N_ART_FSW)
colnames(N_ART_FSW) = c("time", "Lower", "Median", "Upper")



N_ART_FSW_indiv = t(do.call(rbind, lapply(res_best_runs, function(x) x$HIV_positive_On_ART[,1])))


N_ART_FSW_indiv = data.frame(time, rep(c("Pro FSW"), each = length(time)), N_ART_FSW_indiv)


colnames(N_ART_FSW_indiv) = c("time", "variable", as.character(seq(1, length(N_ART_FSW_indiv[1,])-2)))

N_ART_FSW_indiv_melted = reshape2::melt(N_ART_FSW_indiv, id.vars = c("time", "variable"))

colnames(N_ART_FSW_indiv_melted) = c("time", "variable", "run", "value")



# N on ART and N off ART diagnosed
Diagnosed_Off_ART_FSW = data.frame(time, t(do.call(rbind, lapply(lapply(res_best_runs, function(x) {x$I22 + x$I23 + x$I24 + x$I25}), function(x) x[,1]))))
Diagnosed_On_ART_FSW = data.frame(time, t(do.call(rbind, lapply(lapply(res_best_runs, function(x) {x$I32 + x$I33 + x$I34 + x$I35}), function(x) x[,1]))))
Diagnosed_Dropout_ART_FSW = data.frame(time, t(do.call(rbind, lapply(lapply(res_best_runs, function(x) {x$I42 + x$I43 + x$I44 + x$I45}), function(x) x[,1]))))
Diagnosed_ALL_FSW = data.frame(time, t(do.call(rbind, lapply(lapply(res_best_runs, function(x) {x$I22 + x$I23 + x$I24 + x$I25 + x$I32 + x$I33 + x$I34 + x$I35 + x$I42 + x$I43 + x$I44 + x$I45}), function(x) x[,1]))))

Diagnosed_FSW = data.frame(rbind(Diagnosed_Off_ART_FSW, Diagnosed_On_ART_FSW, Diagnosed_Dropout_ART_FSW, Diagnosed_ALL_FSW), 
                           rep(c("Diagnosed Off ART", "Diagnosed On ART", "Dropout", "All diagnosed"), each = length(time)))
colnames(Diagnosed_FSW) = c("time", as.character(seq(1, (length(Diagnosed_Off_ART_FSW[1,])-1))), "variable")

Diagnosed_FSW_melted = reshape2::melt(Diagnosed_FSW, id.vars = c("time", "variable"))

Diagnosed_FSW_melted$group = "FSW"

colnames(Diagnosed_FSW_melted) = c("time", "variable", "run", "value", "group")


Diagnosed_Off_ART_All = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I22 + x$I23 + x$I24 + x$I25}), function(x) {return(x[,c(1, 2, 3, 4, 5, 6, 7, 8)])}), rowSums))))
Diagnosed_On_ART_All = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I32 + x$I33 + x$I34 + x$I35}), function(x) {return(x[,c(1, 2, 3, 4, 5, 6, 7, 8)])}), rowSums))))
Diagnosed_Dropout_ART_All = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I42 + x$I43 + x$I44 + x$I45}), function(x) {return(x[,c(1, 2, 3, 4, 5, 6, 7, 8)])}), rowSums))))

Diagnosed_All = data.frame(rbind(Diagnosed_Off_ART_All, Diagnosed_On_ART_All, Diagnosed_Dropout_ART_All), 
                           rep(c("Diagnosed Off ART", "Diagnosed On ART", "Dropout"), each = length(time)))
colnames(Diagnosed_All) = c("time", as.character(seq(1, (length(Diagnosed_Off_ART_All[1,])-1))), "variable")

Diagnosed_All_melted = reshape2::melt(Diagnosed_All, id.vars = c("time", "variable"))

Diagnosed_All_melted$group = "All"


colnames(Diagnosed_All_melted) = c("time", "variable", "run", "value", "group")


Diagnosed = rbind(Diagnosed_FSW_melted, Diagnosed_All_melted)



Diagnosed_Off_ART_Men = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I22 + x$I23 + x$I24 + x$I25}), function(x) {return(x[,c(5, 6)])}), rowSums))))
Diagnosed_On_ART_Men = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I32 + x$I33 + x$I34 + x$I35}), function(x) {return(x[,c(5, 6)])}), rowSums))))
Diagnosed_Dropout_ART_Men = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I42 + x$I43 + x$I44 + x$I45}), function(x) {return(x[,c(5, 6)])}), rowSums))))

Diagnosed_Men = data.frame(rbind(Diagnosed_Off_ART_Men, Diagnosed_On_ART_Men, Diagnosed_Dropout_ART_Men), 
                           rep(c("Diagnosed Off ART", "Diagnosed On ART", "Dropout"), each = length(time)))
colnames(Diagnosed_Men) = c("time", as.character(seq(1, (length(Diagnosed_Off_ART_Men[1,])-1))), "variable")

Diagnosed_Men_melted = reshape2::melt(Diagnosed_Men, id.vars = c("time", "variable"))

Diagnosed_Men_melted$group = "Men"


colnames(Diagnosed_Men_melted) = c("time", "variable", "run", "value", "group")



Diagnosed_Off_ART_Women = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I22 + x$I23 + x$I24 + x$I25}), function(x) {return(x[,c(1, 2, 3, 4)])}), rowSums))))
Diagnosed_On_ART_Women = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I32 + x$I33 + x$I34 + x$I35}), function(x) {return(x[,c(1, 2, 3, 4)])}), rowSums))))
Diagnosed_Dropout_ART_Women = data.frame(time, t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x)  {x$I42 + x$I43 + x$I44 + x$I45}), function(x) {return(x[,c(1, 2, 3, 4)])}), rowSums))))

Diagnosed_Women = data.frame(rbind(Diagnosed_Off_ART_Women, Diagnosed_On_ART_Women, Diagnosed_Dropout_ART_Women), 
                             rep(c("Diagnosed Off ART", "Diagnosed On ART", "Dropout"), each = length(time)))
colnames(Diagnosed_Women) = c("time", as.character(seq(1, (length(Diagnosed_Off_ART_Women[1,])-1))), "variable")

Diagnosed_Women_melted = reshape2::melt(Diagnosed_Women, id.vars = c("time", "variable"))

Diagnosed_Women_melted$group = "Women"


colnames(Diagnosed_Women_melted) = c("time", "variable", "run", "value", "group")

Diagnosed_Women_Men = rbind(Diagnosed_Women_melted, Diagnosed_Men_melted)

Diagnosed_Women_Men_ratio = data.frame(Diagnosed_Women_melted[,c("time", "variable", "run")],
                                       value = Diagnosed_Women_melted$value/Diagnosed_Men_melted$value)


HIV_deaths = data.frame(time[-1], t(do.call(rbind, lapply(lapply(lapply(res_best_runs, function(x) x$cumuHIVDeaths), rowSums), diff) )))

colnames(HIV_deaths) = c("time", as.character(seq(1, (length(HIV_deaths[1,])-1))))

HIV_deaths_melted = reshape2::melt(HIV_deaths, id.vars = c("time"))







pc_OfWomen_ProFSW = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1]))/
                                         (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1])) +
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2])) + 
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,3])) +
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,4]))+
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,7]))
                                          
                                          
                                         )))

pc_OfWomen_LowFSW = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2]))/
                                         (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1])) +
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2])) + 
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,3])) +
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,4]))+
                                            do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,7])))))



pc_OfWomen_Active_FSW = data.frame(time, t(100*(do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1])) + do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2])))/
                                             (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1])) +
                                                do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2])) + 
                                                do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,3])) +
                                                do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,4]))+
                                                do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,7])))))


pc_OfWomen_GPF = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,3]))/
                                      (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1])) +
                                         do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2])) + 
                                         do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,3])) +
                                         do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,4]))+
                                         do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,7])))))

pc_OfWomen_FormerFSW = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,4]))/
                                            (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1])) +
                                               do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2])) + 
                                               do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,3])) +
                                               do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,4]))+
                                               do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,7])))))



pc_OfMen_Client = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,5]))/
                                       (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,5])) +
                                          do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,6]))+
                                          do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,8])))))


pc_OfMen_GPM = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,6]))/
                                    (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,5])) +
                                       do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,6]))+
                                       do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,8])))))


pc_OfWomen_VF = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,7]))/
                                     (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,1])) +
                                        do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,2])) + 
                                        do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,3])) +
                                        do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,4]))+
                                        do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,7])))))


pc_OfMen_VM = data.frame(time, t(100*do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,8]))/
                                   (do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,5])) +
                                      do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,6]))+
                                      do.call(rbind, lapply(lapply(res_best_runs, function(x) x$N), function(x) x[,8])))))

Ratio_Low_Pro = data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x) {x$frac_N[,2]/ x$frac_N[,1]}))))


frac_by_gender = rbind(pc_OfWomen_ProFSW,
                       pc_OfWomen_LowFSW,
                       pc_OfWomen_GPF,
                       pc_OfWomen_FormerFSW,
                       pc_OfMen_Client,
                       pc_OfMen_GPM,
                       pc_OfWomen_VF,
                       pc_OfMen_VM,
                       pc_OfWomen_Active_FSW,
                       Ratio_Low_Pro
)
frac_by_gender = data.frame(frac_by_gender, group = rep(c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients", "GPM", "Virgin female", "Virgin male", "Active FSW", "Low Pro Ratio"), each = length(time)))
colnames(frac_by_gender) = c("time",as.character(seq(1, (length(pc_OfWomen_ProFSW[1,])-1))),  "variable")
frac_by_gender$variable = factor(frac_by_gender$variable, levels = c("Pro FSW", "Low-level FSW", "GPF", "Former FSW in Cotonou", "Clients", "GPM", "Virgin female", "Virgin male", "Active FSW", "Low Pro Ratio"))


frac_by_gender_melted = reshape2::melt(frac_by_gender, id.vars = c("time", "variable"))

colnames(frac_by_gender_melted) = c("time", "variable", "run", "value")

colnames(pc_OfWomen_ProFSW) = c("time", as.character(seq(1, (length(pc_OfWomen_ProFSW[1,])-1))))
colnames(pc_OfWomen_LowFSW) = c("time", as.character(seq(1, (length(pc_OfWomen_LowFSW[1,])-1))))
colnames(pc_OfWomen_Active_FSW) = c("time", as.character(seq(1, (length(pc_OfWomen_Active_FSW[1,])-1))))
colnames(pc_OfMen_Client) = c("time", as.character(seq(1, (length(pc_OfMen_Client[1,])-1))))


colnames(pc_OfWomen_VF) = c("time", as.character(seq(1, (length(pc_OfWomen_VF[1,])-1))))
colnames(pc_OfMen_VM) = c("time", as.character(seq(1, (length(pc_OfMen_VM[1,])-1))))

pc_OfWomen_ProFSW_melted = reshape2::melt(pc_OfWomen_ProFSW, id.vars = c("time"))
pc_OfWomen_LowFSW_melted = reshape2::melt(pc_OfWomen_LowFSW, id.vars = c("time"))
pc_OfWomen_Active_FSW_melted = reshape2::melt(pc_OfWomen_Active_FSW, id.vars = c("time"))
pc_OfMen_Client_melted = reshape2::melt(pc_OfMen_Client, id.vars = c("time"))



pc_OfWomen_VF_melted = reshape2::melt(pc_OfWomen_VF, id.vars = c("time"))
pc_OfMen_VM_melted = reshape2::melt(pc_OfMen_VM, id.vars = c("time"))







condom_Pro_FSW_comm = t((do.call(rbind, lapply(res_best_runs, function(x)  {x$fc_comm[,1,][,5]}))))
condom_Pro_FSW_noncomm = t((do.call(rbind, lapply(res_best_runs, function(x)  {x$fc_noncomm[,1,][,5]}))))

condom_Pro_FSW = data.frame(time, rbind(condom_Pro_FSW_comm, condom_Pro_FSW_noncomm), rep(c("Commercial", "Non commercial"), each = length(time)))
colnames(condom_Pro_FSW) = c("time", as.character(seq(1, (length(condom_Pro_FSW[1,])-2))), "variable")

condom_Pro_FSW_melted = reshape2::melt(condom_Pro_FSW, id.vars = c("time", "variable"))
colnames(condom_Pro_FSW_melted) = c("time", "variable", "run", "value")


condom_Low_FSW_comm = t((do.call(rbind, lapply(res_best_runs, function(x)  {x$fc_comm[,2,][,5]}))))
condom_Low_FSW_noncomm = t((do.call(rbind, lapply(res_best_runs, function(x)  {x$fc_noncomm[,2,][,5]}))))

condom_Low_FSW = data.frame(time, rbind(condom_Low_FSW_comm, condom_Low_FSW_noncomm), rep(c("Commercial", "Non commercial"), each = length(time)))
colnames(condom_Low_FSW) = c("time", as.character(seq(1, (length(condom_Low_FSW[1,])-2))), "variable")

condom_Low_FSW_melted = reshape2::melt(condom_Low_FSW, id.vars = c("time", "variable"))
colnames(condom_Low_FSW_melted) = c("time", "variable", "run", "value")

condom_GPF_noncomm = data.frame(time, t((do.call(rbind, lapply(res_best_runs, function(x)  {x$fc_noncomm[,3,][,6]})))))
colnames(condom_GPF_noncomm) = c("time", as.character(seq(1, (length(condom_GPF_noncomm[1,])-1))))
condom_GPF_noncomm_melted = reshape2::melt(condom_GPF_noncomm, id.vars = c("time"))

condom_GPM_noncomm = data.frame(time, t((do.call(rbind, lapply(res_best_runs, function(x)  {x$fc_noncomm[,6,][,3]})))))
colnames(condom_GPM_noncomm) = c("time", as.character(seq(1, (length(condom_GPM_noncomm[1,])-1))))
condom_GPM_noncomm_melted = reshape2::melt(condom_GPM_noncomm, id.vars = c("time"))






testing_rate_ratio_F_M = data.frame(time, t((do.call(rbind, lapply(lapply(res_best_runs, function(x)  {x$testing_prob}), function(x) {return(x[,c(3)])}))/do.call(rbind, lapply(lapply(res_best_runs, function(x)  {x$testing_prob}), function(x) {return(x[,c(6)])})))))

colnames(testing_rate_ratio_F_M) = c("time", as.character(seq(1, (length(testing_rate_ratio_F_M[1,])-1))))

testing_rate_ratio_F_M_melted = reshape2::melt(testing_rate_ratio_F_M, id.vars = c("time"))

colnames(testing_rate_ratio_F_M_melted) = c("time", "run", "value")

Diagnosed_Women_Men_ratio_On_ART = Diagnosed_Women_Men_ratio[Diagnosed_Women_Men_ratio$variable == "Diagnosed On ART",]


testing_rate_ratio_F_M_Women_Men_ratio_On_ART = data.frame(x=testing_rate_ratio_F_M_melted$value, y=Diagnosed_Women_Men_ratio_On_ART$value)


pfFSW = data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x)  {x$pfFSW[,1]}))))
colnames(pfFSW) = c("time", as.character(seq(1, (length(pfFSW[1,])-1))))
pfFSW_melted = reshape2::melt(pfFSW, id.vars = c("time"))
colnames(pfFSW_melted) = c("time", "run", "value")



ART_inits = data.frame(time = time[-length(time)], t(do.call(rbind,lapply(res_best_runs, function(x)  {diff(rowSums(x$cumuARTinitiations))}))))
colnames(ART_inits) = c("time", as.character(seq(1, (length(ART_inits[1,])-1))))
ART_inits_melted = reshape2::melt(ART_inits, id.vars = c("time"))
colnames(ART_inits_melted) = c("time", "run", "value")



# Ratio female:male of rate of initiation per infected = 8.0/4.68 = 1.7
# to compare to the data, I am taking the number of initiations that years as a fraction of the total nubmer of PLHIV that year



# DENOMINTOR = ALL POSTIIVES
female_ART_init_rate_from_all_HIV_pos = data.frame(time = time[-length(time)], t(do.call(rbind,lapply(res_best_runs, function(x)  {
  a = diff(rowSums(x$cumuARTinitiations[,c(1,2,3,4)])) # number of female ART initiations per year
  b = rowSums(x$HIV_positive[,c(1,2,3,4)])[-length(time)] # number of HIV positive females at the beginning of the year
  return(a/b)
}))))


male_ART_init_rate_from_all_HIV_pos = data.frame(time = time[-length(time)], t(do.call(rbind,lapply(res_best_runs, function(x)  {
  c = diff(rowSums(x$cumuARTinitiations[,c(5,6)])) # number of male ART initiations per year
  d = rowSums(x$HIV_positive[,c(5,6)])[-length(time)] # number of HIV positive males at the beginning of the year
  return(c/d)
}))))







# DENOMINTOR = ALL POSTIIVES NOT ON ART

female_ART_init_rate_from_all_HIV_pos_NOT_ON_ART = data.frame(time = time[-length(time)], t(do.call(rbind,lapply(res_best_runs, function(x)  {
  a = diff(rowSums(x$cumuARTinitiations[,c(1,2,3,4)])) # number of female ART initiations per year
  b = rowSums(x$HIV_positive[,c(1,2,3,4)] - x$HIV_positive_On_ART[,c(1,2,3,4)])[-length(time)] # number of HIV positive females not on ART at the beginning of the year
  return(a/b)
}))))
# female_ART_init_rate_from_all_HIV_pos_NOT_ON_ART



male_ART_init_rate_from_all_HIV_pos_NOT_ON_ART = data.frame(time = time[-length(time)], t(do.call(rbind,lapply(res_best_runs, function(x)  {
  a = diff(rowSums(x$cumuARTinitiations[,c(5,6)])) # number of female ART initiations per year
  b = rowSums(x$HIV_positive[,c(5,6)] - x$HIV_positive_On_ART[,c(5,6)])[-length(time)] # number of HIV positive females not on ART at the beginning of the year
  return(a/b)
}))))
# male_ART_init_rate_from_all_HIV_pos_NOT_ON_ART





# FEMALE TO MALE RATIO OF INITIATION RATES WHERE DENOMINATOR IS ART POSITIVE NOT ON ART

FtM_ratio_ART_initiation_rates_from_HIVpos_not_on_ART = data.frame(time = time[-length(time)], t(do.call(rbind, lapply(res_best_runs, function(x)  {
  a = diff(rowSums(x$cumuARTinitiations[,c(1,2,3,4)])) # number of female ART initiations per year
  b = rowSums(x$HIV_positive[,c(1,2,3,4)] - x$HIV_positive_On_ART[,c(1,2,3,4)])[-length(time)] # number of HIV positive females not on ART at the beginning of the year
  c = diff(rowSums(x$cumuARTinitiations[,c(5,6)])) # number of female ART initiations per year
  d = rowSums(x$HIV_positive[,c(5,6)] - x$HIV_positive_On_ART[,c(5,6)])[-length(time)] # number of HIV positive females not on ART at the beginning of the year
  return((a/b)/(c/d))
}))))

colnames(FtM_ratio_ART_initiation_rates_from_HIVpos_not_on_ART) = c("time", as.character(seq(1, (length(FtM_ratio_ART_initiation_rates_from_HIVpos_not_on_ART[1,])-1))))

FtM_ratio_ART_initiation_rates_from_HIVpos_not_on_ART_melted = reshape2::melt(FtM_ratio_ART_initiation_rates_from_HIVpos_not_on_ART, id.vars = c("time"))



ART_inits_ratio_F_over_M =  data.frame(time = time[-length(time)], t(do.call(rbind,lapply(res_best_runs, function(x)  {diff(rowSums(x$cumuARTinitiations[,c(1,2,3,4)])) / 
    diff(rowSums(x$cumuARTinitiations[,c(5,6)]))}))))
ART_inits_ratio_F_over_M_melted = reshape2::melt(ART_inits_ratio_F_over_M, id.vars = c("time"))




ART_inits_cumu = data.frame(time = time, t(do.call(rbind,lapply(res_best_runs, function(x)  {rowSums(x$cumuARTinitiations)}))))
colnames(ART_inits_cumu) = c("time", as.character(seq(1, (length(ART_inits_cumu[1,])-1))))
ART_inits_cumu_melted = reshape2::melt(ART_inits_cumu, id.vars = c("time"))
colnames(ART_inits_cumu_melted) = c("time", "run", "value")



ART_REinits_cumu = data.frame(time = time, t(do.call(rbind,lapply(res_best_runs, function(x)  {rowSums(x$cumuARTREinitiations)}))))
colnames(ART_REinits_cumu) = c("time", as.character(seq(1, (length(ART_REinits_cumu[1,])-1))))
ART_REinits_cumu_melted = reshape2::melt(ART_REinits_cumu, id.vars = c("time"))
colnames(ART_REinits_cumu_melted) = c("time", "run", "value")



ART_REinits = data.frame(time = time[-length(time)], t(do.call(rbind,lapply(res_best_runs, function(x)  {diff(rowSums(x$cumuARTREinitiations))}))))
colnames(ART_REinits) = c("time", as.character(seq(1, (length(ART_REinits[1,])-1))))
ART_REinits_melted = reshape2::melt(ART_REinits, id.vars = c("time"))
colnames(ART_REinits_melted) = c("time", "run", "value")




frac_ART_inits_REinits = data.frame(time = time[-length(time)])
for(i in 2:(length(res_best_runs)+1))
{
  frac_ART_inits_REinits[,i] = ART_REinits[,i]/(ART_inits[,i] + ART_REinits[,i]) * 100
}
colnames(frac_ART_inits_REinits) = c("time", as.character(seq(1, (length(frac_ART_inits_REinits[1,])-1))))
frac_ART_inits_REinits_melted = reshape2::melt(frac_ART_inits_REinits, id.vars = c("time"))
colnames(frac_ART_inits_REinits_melted) = c("time", "run", "value")




ART_init_rate_from_all_HIV_pos_NOT_ON_ART_by_sex = data.frame(
  Women_2014 = unlist(female_ART_init_rate_from_all_HIV_pos_NOT_ON_ART[which(time == 2014),c(2:(length(res_best_runs)+1))]),
  Women_2015 = unlist(female_ART_init_rate_from_all_HIV_pos_NOT_ON_ART[which(time == 2015),c(2:(length(res_best_runs)+1))]),
  
  Men_2014 = unlist(male_ART_init_rate_from_all_HIV_pos_NOT_ON_ART[which(time == 2014),c(2:(length(res_best_runs)+1))]),
  Men_2015 = unlist(male_ART_init_rate_from_all_HIV_pos_NOT_ON_ART[which(time == 2015),c(2:(length(res_best_runs)+1))])
)


ART_init_rate_from_all_HIV_pos_NOT_ON_ART_by_sex_melted = reshape2::melt(ART_init_rate_from_all_HIV_pos_NOT_ON_ART_by_sex)

levels(ART_init_rate_from_all_HIV_pos_NOT_ON_ART_by_sex_melted) = c("Women_2014","Women_2015", "Men_2014","Men_2015")


###################






require(ggplot2)




# want plots on number of initiations of ART, dropout, reinitiations! compare to data points from the reports!

# need to do difference between years for cumulative art inits etc


prev_axes = data.frame(variable = c(rep("Pro FSW", 2),
                                    rep("Clients", 2),
                                    rep("Women", 2),
                                    rep("Men", 2),
                                    rep("Low-level FSW", 2)),
                       time = c(rep(c(1986, 2025), 5)),
                       value = c(0, 70, 0, 70, 0, 15, 0, 15, 0, 70)
)
prev_points_80s = prev_points_all[c(1,2,3),]



prev_points$fitted = c("Yes", "No", "No", "Yes", "No", "Yes", "No", "Yes", 
                       "No", "No", "Yes", "No", "No", "No",
                       "No", "No", "Yes", 
                       "No", "No", "Yes", 
                       "No", "No")

prev_points[prev_points$time == "2015", "lower"][1] = 13.79


# plot prevalence in each group indiv runs
g1=ggplot() + geom_line(data = prev_indiv_melted, aes(x = time, y = value, factor = run), alpha = 0.123) + theme_bw() + facet_wrap(~variable, scales = "free") +#labs(y = "HIV Prevalence (%)") +
  geom_point(data = prev_points, aes(x = time, y = value), size = 2)+ geom_errorbar(data = prev_points, aes(x = time, ymin = lower, ymax = upper, col = fitted), size = 1)+
  scale_colour_manual(values = c("blue", "red"))+
  scale_x_continuous(breaks = seq(1986, 2017, 1), limits = c(1986,2017))+
  
  geom_point(data = prev_points_80s, aes(x = time, y = value), colour = "blue", size = 2)+
  geom_blank(data = prev_axes, aes(x = time, y = value))+
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5))+
  # scale_x_continuous(breaks = seq(1986, 2020, 2), limits = c(1986,2020))+
  ggtitle("FITTING: HIV PREVALENCE (%)")


g1


# plot prevalence in each group indiv runs
g1_paper = ggplot() + geom_line(data = prev_indiv_melted[prev_indiv_melted$variable != "Low-level FSW",], aes(x = time, y = value, factor = run), alpha = 0.08123) + theme_bw() + facet_wrap(~variable, scales = "free") +#labs(y = "HIV Prevalence (%)") +
  geom_point(data = prev_points[prev_points$variable != "Low-level FSW",], aes(x = time, y = value), shape = 15, size = 2)+
  geom_errorbar(data = prev_points[prev_points$variable != "Low-level FSW",], aes(x = time, ymin = lower, ymax = upper, col = fitted), size = 1)+
  scale_colour_manual(values = c("blue", "red"))+
  scale_x_continuous(breaks = seq(1986, 2018, 4), limits = c(1986, 2018))+
  labs(y="HIV Prevalence (%)", x="Year")+
  geom_line(data = prev[prev$variable != "Low-level FSW",], aes(x = time, y = Median))+
  geom_ribbon(data = prev[prev$variable != "Low-level FSW",], aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.2)+
  geom_point(data = prev_points_80s, aes(x = time, y = value), colour = "blue", size = 2)+
  geom_blank(data = prev_axes[prev_axes$variable != "Low-level FSW",], aes(x = time, y = value))+
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5))+
  # scale_x_continuous(breaks = seq(1986, 2020, 2), limits = c(1986,2020))+
  ggtitle("")


g1_paper
gc()

# plot prevalence in each group indiv runs
g1b_paper = ggplot() + geom_line(data = prev_indiv_melted[prev_indiv_melted$variable != "Low-level FSW",], aes(x = time, y = value, factor = run), alpha = 0.08123) + theme_bw() + facet_wrap(~variable, scales = "free") +#labs(y = "HIV Prevalence (%)") +
  geom_point(data = prev_points[prev_points$variable != "Low-level FSW" & prev_points$fitted == "Yes",], aes(x = time, y = value), shape = 15, size = 2)+
  geom_errorbar(data = prev_points[prev_points$variable != "Low-level FSW" & prev_points$fitted == "Yes",], aes(x = time, ymin = lower, ymax = upper, col = fitted), size = 1, width = 1)+
  scale_colour_manual(values = c("black", "black"))+
  scale_x_continuous(breaks = seq(1986, 2018, 8), limits = c(1986, 2018))+
  labs(y="HIV Prevalence (%)", x="Year")+
  geom_line(data = prev[prev$variable != "Low-level FSW",], aes(x = time, y = Median))+
  geom_line(data = prev[prev$variable != "Low-level FSW",], aes(x = time, y = Lower), linetype = "longdash")+
  geom_line(data = prev[prev$variable != "Low-level FSW",], aes(x = time, y = Upper), linetype = "longdash")+
  geom_ribbon(data = prev[prev$variable != "Low-level FSW",], aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.2)+
  geom_blank(data = prev_axes[prev_axes$variable != "Low-level FSW",], aes(x = time, y = value))+
  theme(text = element_text(size=30),plot.title = element_text(hjust = 0.5), legend.position = "")+
  # scale_x_continuous(breaks = seq(1986, 2020, 2), limits = c(1986,2020))+
  ggtitle("")


g1b_paper
gc()


prev_indiv_melted_FSW <- prev_indiv_melted[prev_indiv_melted$variable == "Pro FSW",]
prev_FSW <- prev[prev$variable == "Pro FSW",]
prev_points_FSW <- prev_points[prev_points$variable == "Pro FSW",]



# plot prevalence in each group indiv runs
g1c_paper = ggplot() + geom_line(data = prev_indiv_melted_FSW, aes(x = time, y = value, factor = run), alpha = 0.08123) +
  theme_bw() + #facet_wrap(~variable, scales = "free") +#labs(y = "HIV Prevalence (%)") +
  geom_point(data = prev_points_FSW, aes(x = time, y = value, col = fitted, size = fitted), shape = 15)+
  # geom_errorbar(data = prev_points_FSW[prev_points_FSW$variable != "Low-level FSW" & prev_points_FSW$fitted == "Yes",], aes(x = time, ymin = lower, ymax = upper, col = fitted), size = 1, width = 1)+
  geom_errorbar(data = prev_points_FSW, aes(x = time, ymin = lower, ymax = upper, col = fitted, size = fitted),  width = 1)+
  scale_colour_manual(values = c("gray20", "black"))+
  scale_size_manual(values = c(1, 2.5))+
  scale_x_continuous(breaks = seq(1986, 2018, 8), limits = c(1986, 2018))+
  labs(y="", x="Year")+
  geom_line(data = prev_FSW[prev_FSW$variable != "Low-level FSW",], aes(x = time, y = Median))+
  geom_line(data = prev_FSW[prev_FSW$variable != "Low-level FSW",], aes(x = time, y = Lower), linetype = "longdash")+
  geom_line(data = prev_FSW[prev_FSW$variable != "Low-level FSW",], aes(x = time, y = Upper), linetype = "longdash")+
  geom_ribbon(data = prev_FSW[prev_FSW$variable != "Low-level FSW",], aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.2)+
  # geom_blank(data = prev_axes[prev_axes$variable != "Low-level FSW",], aes(x = time, y = value))+
  theme(text = element_text(size=50),plot.title = element_text(hjust = 0.5), legend.position = "")+
  # scale_x_continuous(breaks = seq(1986, 2020, 2), limits = c(1986,2020))+
  ggtitle("")


g1c_paper
gc()


prev_indiv_melted_GP <- prev_indiv_melted[prev_indiv_melted$variable == "Women" | prev_indiv_melted$variable == "Men",]
prev_GP <- prev[prev$variable == "Women" | prev$variable == "Men",]
prev_points_GP <- prev_points[prev_points$variable == "Women" | prev_points$variable == "Men",]


# plot prevalence in each group indiv runs
g1d_paper = ggplot() + geom_line(data = prev_indiv_melted_GP, aes(x = time, y = value, factor = run), alpha = 0.08123) +
  theme_bw() + facet_wrap(~variable, scales = "free") +#labs(y = "HIV Prevalence (%)") +
  geom_point(data = prev_points_GP[prev_points_GP$variable != "Low-level FSW" ,], aes(x = time, y = value, size = fitted), shape = 15)+
  geom_errorbar(data = prev_points_GP[prev_points_GP$variable != "Low-level FSW" ,], aes(x = time, ymin = lower, ymax = upper, col = fitted, size = fitted), width = 1)+
  scale_colour_manual(values = c("gray20", "black"))+
  scale_size_manual(values = c(1, 2.5))+
  
  scale_x_continuous(breaks = seq(1986, 2018, 8), limits = c(1986, 2018))+
  labs(y="", x="Year")+
  geom_line(data = prev_GP[prev_GP$variable != "Low-level FSW",], aes(x = time, y = Median))+
  geom_line(data = prev_GP[prev_GP$variable != "Low-level FSW",], aes(x = time, y = Lower), linetype = "longdash")+
  geom_line(data = prev_GP[prev_GP$variable != "Low-level FSW",], aes(x = time, y = Upper), linetype = "longdash")+
  geom_ribbon(data = prev_GP[prev_GP$variable != "Low-level FSW",], aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.2)+
  # geom_blank(data = prev_axes[prev_axes$variable != "Low-level FSW",], aes(x = time, y = value))+
  theme(text = element_text(size=50),plot.title = element_text(hjust = 0.5), legend.position = "")+
  # scale_x_continuous(breaks = seq(1986, 2020, 2), limits = c(1986,2020))+
  ggtitle("")


g1d_paper
gc()



ART_data_points_with_numbers_FSW = data.frame(time = c(2012, 2014, 2015, 2016, 2017),
                                              Lower = c(27,	34,	42, 53*0.85,	73*0.85),
                                              Upper = c(37, 46, 56, 53*1.15, 73*1.15),
                                              variable = rep("Pro FSW", 5),
                                              datatype = c("black", "black", "red", "black", "black"))

ART_data_points_with_numbers_FSW_points = data.frame(time = c(2012, 2014, 2015, 2016, 2017),
                                                     value = c(32,	40,	49, 53,	73),
                                                     variable = rep("Pro FSW", 5),
                                                     datatype = c("black", "black", "red", "black", "black"))


ART_upper_bound = data.frame(time = c(2016, 2017),
                                                     value = c(122,	137),
                                                     variable = rep("Pro FSW", 2)
)

ART_data_points_with_numbers_FSW_full = ART_data_points_with_numbers_FSW

ART_data_points_with_numbers_FSW = ART_data_points_with_numbers_FSW[1:3,]

g7=ggplot() + geom_line(data = N_ART_FSW_indiv_melted, aes(x = time, y = value, factor = run), alpha = 0.2) +
  ggtitle("FITTING: Number of professional FSW on ART; fitting 2015 only") +theme_bw() +
  geom_errorbar(data = ART_data_points_with_numbers_FSW, size = 1.5, aes(x = time, ymin = Lower, ymax = Upper, col = datatype)) + labs(x = "Year") +
  scale_x_continuous(breaks = seq(2000, 2034, 1), limits = c(2000,2034))+
  labs(y = "")+
  geom_point(data = ART_data_points_with_numbers_FSW_points,aes(x = time, y = value, col = datatype))+
  theme(legend.position = "none", text = element_text(size=20),plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=24),plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=18),
        legend.key.size = unit(1.3, 'lines'))+ annotate("text", x = 2008.8, y = 120, size = 7,col = "black",
                                                        label = "Data is # FSW file active +- 15% \n Fitting to 2015 estimate only \n 2016&2017 are counterfactuals \n triangles are # on ART with TasP (upper bound)")+
  geom_point(data = ART_upper_bound, aes(x = time, y = value), size = I(3), shape = 2, fill="black", col = "black")+
  scale_color_manual(values = c('blue', 'red', 'darkred'))

g7


g7_paper =ggplot() + geom_line(data = N_ART_FSW_indiv_melted, aes(x = time, y = value, factor = run), alpha = 0.1) +
  theme_bw() +
  geom_errorbar(data = ART_data_points_with_numbers_FSW_full[,], width = 0.5, aes(x = time, ymin = Lower, ymax = Upper, size =datatype, col = datatype)) + labs(x = "Year") +
  scale_x_continuous(breaks = seq(2000, 2018, 3), limits = c(2000, 2018))+
  labs(y = "")+
  geom_point(aes(x = 2015, y = 49), shape = 15, size  =2)+
  geom_ribbon(data = N_ART_FSW, aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.2)+
  geom_line(data = N_ART_FSW, aes(x = time, y = Median), alpha = 1)+
  geom_line(data = N_ART_FSW, aes(x = time, y = Lower), linetype = "longdash", alpha = 1)+
  geom_line(data = N_ART_FSW, aes(x = time, y = Upper), linetype = "longdash", alpha = 1)+

  theme(legend.position = "none", text = element_text(size=30),plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size=50),plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=30),
        legend.key.size = unit(1.3, 'lines'))+ 
  scale_colour_manual(values = c("gray20", "black"))+
  scale_size_manual(values = c(1, 2.5))+
  scale_linetype_manual(values = c("solid"))
# scale_linetype_manual(values = c("longdash", "solid"))

g7_paper
gc()




ART_data_points = data.frame(time = c(2011, 2012, 2013, 2014, 2015, 2016, 2017,
                                      
                                      2012, 2014, 2015, 2016, 2017
),
Lower = c(0.33, 0.42, 0.44, 0.39, 0.44, 0.51, 0.57,
          
          0.09,	0.14,	0.2,	0.22,	0.3
          
),
Upper = c(0.52, 0.66, 0.69, 0.61, 0.69, 0.8, 0.8,
          0.13, 0.2, 0.28, 0.3, 0.42
          
),
variable = c("All", "All", "All", "All", "All", "All", "All",
             "Pro FSW", "Pro FSW", "Pro FSW", "Pro FSW", "Pro FSW"))

# ART_data_points$Data = c(rep("UNAIDS (fitted 2011)", 7), rep("File active FSW (X-validation)", 5))

art_cov = data.frame(time = 2017, Lower = 0.6, Upper = 0.91, variable = "All" )#, Data = "Audit report 2017 (fitting)")
# 
ART_data_points = rbind(ART_data_points, art_cov)
# 
# ART_data_points=ART_data_points[-as.numeric(rownames(ART_data_points[ART_data_points$time==2017 & ART_data_points$Data=="UNAIDS (fitted 2011)", ])),]

ART_text = data.frame(time = 2010, variable = "Pro FSW", value = 50)

ART_data_points$Data = c("Yes", rep("No", 11), "Yes")


g8=ggplot() + geom_line(data = ART_indiv_melted, aes(x = time, y = value*100,  factor = run), alpha = 0.3) + theme_bw() + facet_wrap(~variable, scales = "free") + labs(y = "ART coverage (%)") +
  ggtitle("FITTING: ART coverage (%) in total population, 2011 & 2017 ONLY. \n NOT fitting to FSW coverage by %")+
  geom_errorbar(data = ART_data_points, aes(x = time, ymin = Lower*100, ymax = Upper*100, col = Data), size = I(1))+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20)) +   scale_colour_manual(values = c("blue", "red"))+
  theme(legend.position = "top") + geom_text(data = ART_text, aes(x = time, y = value, factor = variable),
                                             label = "2016 & 2017 are \n NOT including FSW on TasP/PrEP \n so we expect a higher coverage \n n.b. we fit to # FSW on ART not %", size = 5)
g8



ART_data_points = ART_data_points[!(ART_data_points$time == 2017 & ART_data_points$Data == "No"),]


g8_paper=ggplot() + geom_line(data = ART_indiv_melted[ART_indiv_melted$variable == "All",], aes(x = time, y = value*100,  factor = run), alpha = 0.1) + theme_bw() + #facet_wrap(~variable, scales = "free") +
  labs(y = "", x= "Year") +
  geom_errorbar(data = ART_data_points[ART_data_points$variable == "All",], aes(x = time, ymin = Lower*100, ymax = Upper*100,col = Data, size = Data),  width = 0.7)+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=50)) + 
  scale_colour_manual(values = c("gray20", "black"))+
  scale_size_manual(values = c(1, 2.5))+
  theme(legend.position = "") + 
  scale_x_continuous(breaks = seq(2000, 2018, 3), limits = c(2000, 2018))+
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100))+
  geom_ribbon(data = ART_coverage[ART_coverage$variable == "All",], aes(x = time, ymin = 100*Lower, ymax = 100*Upper), alpha = 0.3)+
  geom_line(data = ART_coverage[ART_coverage$variable == "All",], aes(x = time, y = 100*Median), alpha = 1)+
  geom_line(data = ART_coverage[ART_coverage$variable == "All",], aes(x = time, y = 100*Lower), linetype = "longdash", alpha = 1)+
  geom_line(data = ART_coverage[ART_coverage$variable == "All",], aes(x = time, y = 100*Upper), linetype = "longdash", alpha = 1)
g8_paper


jpeg(filename=paste0(getwd(), "/", batch_folder, "/", "paper_1",".jpg"), width = 1200, height = 700)
print(g1b_paper)
dev.off()
jpeg(filename=paste0(getwd(), "/", batch_folder, "/", "paper_2",".jpg"), width = 1200, height = 700)
print(g7_paper)
dev.off()
jpeg(filename=paste0(getwd(),  "/",batch_folder, "/", "paper_3",".jpg"), width = 1200, height = 700)
print(g8_paper)
dev.off()



Ntot = data.frame(c(seq(1986, 2035)), t(apply(do.call(rbind, lapply(res_best_runs, function(x) x$Ntot)), 2, cotonou::quantile_95)))
colnames(Ntot) = c("time", "Lower", "Median", "Upper")

# plot total population size
g6=ggplot(Ntot) + geom_line(aes(x = time, y = Median)) + geom_ribbon(aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.5) +
  ggtitle("FITTING: Total population size of Grand Cotonou; red points fitted")+
  theme_bw() + labs(y = "") +
  geom_point(data = Ntot_data_points, aes(x = time, y = point, color = colour), size = I(2), shape = 15) + geom_errorbar(data = Ntot_data_points, aes(x = time, ymax = upper, ymin = lower, color = colour), width = 2)+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20))
g6












frac_by_gender_discard_points_no_FSW_LB = data.frame(variable = c("Pro FSW", "Clients", "Virgin female", "Virgin male", "Active FSW", "Low Pro Ratio"),
                                                     min = c(0, 0.074, 0.07896475, 0.07039551, 0.0048, .01),
                                                     max = c(0.0143/2, 0.3 , 0.2, 0.17,  0.0143, .05))



# plot fraction in each group
g2=ggplot(frac_by_gender_melted) + geom_line(aes(x = time, y = value, factor = run))  +
  ggtitle("FITTING: Summary of all demography fitting")+
  theme_bw() + labs(y = "Percent in each group in their respective gender (%)") +  facet_wrap(~variable, scales = "free") +
  # geom_point(data = frac_N_data_points, aes(x = time, y = point), size = I(2), color = "red", shape = 15) +
  geom_hline(data = frac_by_gender_discard_points_no_FSW_LB, aes(yintercept = 100*min), size = I(0.5), color = "red", linetype = 1, alpha = 0.7) +
  geom_hline(data = frac_by_gender_discard_points_no_FSW_LB, aes(yintercept = 100*max), size = I(0.5), color = "red", linetype = 1, alpha = 0.7)+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20))
g2






g3=ggplot() + geom_line(data = pc_OfWomen_ProFSW_melted, aes(x = time, y = value, factor = variable), alpha = 1) + theme_bw() + labs(y = "% of all women 15-59 that are pro FSW")+
  ggtitle("FITTING: Fraction professional FSW")+
  
  theme(text = element_text(size=20)) + geom_hline(yintercept = 0.24,col = "orange", linetype="dashed", size = 2) + geom_hline(yintercept = 0.72, linetype="dashed",col = "red", size = 2)+
  theme(text = element_text(size=24),
        legend.text=element_text(size=18),
        legend.key.size = unit(1.3, 'lines'))+ annotate("text", x = 2002.8, y = 0.7, size = 7,col = "red", label = "Fitted to below this line")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=24),
        legend.text=element_text(size=18),
        legend.key.size = unit(1.3, 'lines'))+ annotate("text", x = 2003.8, y = 0.5, size = 7,col = "orange", label = "Original estimate lower bound - \n I relaxed this in order to get fits \n Using point estimate # pro FSW from mapping") +
  geom_hline(aes( yintercept = 0.001944823*100), size =2,linetype="dashed", col = "red", fill ="red")+
  # geom_point(aes(x = 2012, y = 0.001944823*100), shape = 17, size =I(3), col = "red", fill ="red")
  annotate("text", x = 2003.8, y = 0.15, size = 7,col = "red", label = "Using lower bound # pro FSW from mapping") 
  
g3







g4=ggplot() + geom_line(data = pc_OfWomen_LowFSW_melted, aes(x = time, y = value, factor = variable), alpha = 0.3) + theme_bw() + labs(y = "% of all women 15-59 that are low level FSW")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20))


g4b = ggplot() + geom_line(data = pc_OfWomen_VF_melted, aes(x = time, y = value, factor = variable), alpha = 0.3) + theme_bw() + labs(y = "% of all women 15-59 that are virgins")+
  ggtitle("FITTING: Fraction women not yet sexually active")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20))+ geom_hline(yintercept = c(7.9, 20), col = "red", size = 2) + 
  theme(text = element_text(size=24), legend.text=element_text(size=18),
        legend.key.size = unit(1.3, 'lines'))+ annotate("text", x = 2010, y = 19, size = 7,col = "red", label = "Screening between lines")

g4b


g5=ggplot() + geom_line(data = pc_OfMen_Client_melted, aes(x = time, y = value, factor = variable), alpha = 1) + theme_bw() + labs(y = "% of all men 15-59 that are clients")+
  ggtitle("FITTING: Fraction men clients")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20)) + geom_hline(yintercept = c(7.4, 30), col = "red", size = 2) + 
  theme(text = element_text(size=24), legend.text=element_text(size=18),
        legend.key.size = unit(1.3, 'lines'))+ annotate("text", x = 2003.8, y = 28, size = 7,col = "red", label = "Screening between lines")

g5

g5b = ggplot() + geom_line(data = pc_OfMen_VM_melted, aes(x = time, y = value, factor = variable), alpha = 0.3) + theme_bw() + labs(y = "% of all men 15-59 that are virgins")+
  ggtitle("FITTING: Fraction men not yet sexually active")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20))+ geom_hline(yintercept = c(7, 17), col = "red", size = 2)+ 
  theme(text = element_text(size=24), legend.text=element_text(size=18),
        legend.key.size = unit(1.3, 'lines'))+ annotate("text", x = 2010, y = 15, size = 7,col = "red", label = "Screening between lines")

g5b
# # plot prevalence in each group
# ggplot() + geom_line(data = prev, aes(x = time, y = Median))+ geom_ribbon(data = prev, aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.5) + theme_bw() + facet_wrap(~variable, scales = "free") + labs(y = "prevalence (%)") +
#   geom_point(data = prev_points, aes(x = time, y = value))+ geom_errorbar(data = prev_points, aes(x = time, ymin = lower, ymax = upper)) +
#   geom_point(data = prev_points_80s, aes(x = time, y = value), colour = "red")+
#   geom_blank(data = prev_axes, aes(x = time, y = value))
# 



# ggplot(ART_coverage) +
#   geom_line(aes(x = time, y = Median))+ geom_ribbon(aes(x = time, ymin = Lower, ymax = Upper), alpha = 0.5) + theme_bw() +
#   facet_wrap(~variable) + labs(y = "ART coverage ") +
#   geom_errorbar(data = ART_data_points, aes(x = time, ymin = Lower, ymax = Upper), colour = "darkred")



# ART_indiv_melted=ART_indiv_melted[ART_indiv_melted$variable == "All",]
# ART_data_points=ART_data_points[ART_data_points$variable == "All",]
# 
# g8=ggplot() + geom_line(data = ART_indiv_melted, aes(x = time, y = value*100,  factor = run), alpha = 1) +
#   theme_bw() +  labs(y = "ART coverage (%)", col = "Fitted") +
#   ggtitle("FITTING: ART coverage (%) in total population, 2011 & 2017 ONLY. \n NOT fitting to FSW coverage by %")+
#   geom_errorbar(data = ART_data_points, aes(x = time, ymin = Lower*100, ymax = Upper*100, col = Data), size = 2, width = 2)+
#   theme(plot.title = element_text(hjust = 0.5),text = element_text(size=40)) +   scale_colour_manual(values = c("blue", "red"))+
#   theme(legend.position = "top")
# g8

FOI_point = data.frame(time = 2015, variable = "Pro FSW", value = 3)
g9=ggplot() + geom_line(data = lambda_sum_0_indiv_melted[lambda_sum_0_indiv_melted$time > 2009,], aes(x = time, y = value*100,  factor = run), alpha = 0.3) + 
  ggtitle("FITTING: susceptible professional FSW FOI below 3% in 2015")+
  geom_point(data= FOI_point, aes(x = time, y = value), col = "red", fill = "red", size = 3, shape = 25)+ theme_bw() + facet_wrap(~variable, scales = "free") + 
  labs(y = "Force of infection on susceptibles off PrEP (%)")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20)) 
g9

Diagnosed_FSW_melted$col = Diagnosed_FSW_melted$variable
Diagnosed_FSW_melted$col = ifelse(Diagnosed_FSW_melted$col == "Diagnosed Off ART", "#7fc97f", 
                                  ifelse(Diagnosed_FSW_melted$col == "Diagnosed On ART", "#beaed4",
                                         ifelse(Diagnosed_FSW_melted$col == "Dropout", "#fdc086", 0)))


# 2012 FSW IBBA ever tested 0.75 - 0.83
# EQS mapping 2012 N = 889 - 1391
# 2012 prevalence = 27.4% [0.23 - 0.322]
# lower N diagnosed = 0.75 * 889 * 0.23 =
# upper N diagnosed = 0.83 * 1391 * 0.322 =


# bounds of infected ever tested = 304 * 0.75, 304 * 0.83 = 228 - 252
pc_diagnosed = data.frame(time = 2012, Lower = 153, Upper = 372, variable = "Pro FSW")

g10=ggplot() + geom_line(data = Diagnosed_FSW_melted, aes(x = time, y = value,  factor = run, colour = variable), alpha = 0.25) + theme_bw() + labs(y = "Numbers of pro FSW in Grand Cotonou")+
  ggtitle("CROSS-VALIDATION")+
  scale_x_continuous(breaks = seq(min(time), max(time), 4) ) +
  scale_y_continuous(breaks = seq(min(Diagnosed$value), #max(Diagnosed_FSW_melted$value), 20) )+
                                  380, 40) )+
  theme(text = element_text(size=20), legend.position = "top")+  
  geom_errorbar(data = ART_data_points_with_numbers_FSW, size=1.5,aes(x = time, ymin = Lower, ymax = Upper), col ="#4daf4a") +
  labs(x = "Year", colour = "Care state")+  
  geom_errorbar(data = pc_diagnosed, size=1.5, aes(x = time, ymin = Lower, ymax = Upper), col ="#e41a1c") +
  scale_colour_brewer(type = "qual", palette = 6) + 
  guides(color = guide_legend(override.aes = list(size=5, alpha = 1)))+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=24),
        legend.text=element_text(size=18))+ annotate("text", x = 2003.8, y = 222, size = 7,col = "#e41a1c", label = "DATA: Number of FSW * \n prevalence * \n % ever tested \n (cross-validation)")

g10
# # 
# ggplot() + geom_line(data = Diagnosed_All_melted, aes(x = time, y = value, factor = variable, factor = run, colour = variable), alpha = 1) + theme_bw() + labs(y = "Numbers of all people in Grand Cotonou")+
#   scale_x_continuous(breaks = seq(min(time), max(time)+1, 4) ) +
#   scale_y_continuous(breaks = seq(min(Diagnosed$value), max(Diagnosed$value), 1000) ) +
#   geom_vline(xintercept=2005, col = "black", alpha = 0.5, linetype = "dashed") +
#   geom_text(aes(x=2005, label="CD4 < 200", y=-400))+
#   geom_vline(xintercept=2012, col = "black", alpha = 0.5, linetype = "dashed") +
#   geom_text(aes(x=2012, label="CD4 < 350", y=-800))+
#   geom_vline(xintercept=2015, col = "black", alpha = 0.5, linetype = "dashed") +
#   geom_text(aes(x=2015, label="CD4 < 500", y=-400))+
#   geom_vline(xintercept=2016, col = "black", alpha = 0.5, linetype = "dashed") +
#   geom_text(aes(x=2016, label="Everyone", y=-800)) 
# # 


# I2x= Diagnosed_All_melted[Diagnosed_All_melted$variable == "Diagnosed Off ART",]
# I3x= Diagnosed_All_melted[Diagnosed_All_melted$variable == "Diagnosed On ART",]
# I4x= Diagnosed_All_melted[Diagnosed_All_melted$variable == "Dropout",]


I2xF= Diagnosed_Women_melted[Diagnosed_Women_melted$variable == "Diagnosed Off ART",]
I3xF= Diagnosed_Women_melted[Diagnosed_Women_melted$variable == "Diagnosed On ART",]
I4xF= Diagnosed_Women_melted[Diagnosed_Women_melted$variable == "Dropout",]


I2xM= Diagnosed_Men_melted[Diagnosed_Men_melted$variable == "Diagnosed Off ART",]
I3xM= Diagnosed_Men_melted[Diagnosed_Men_melted$variable == "Diagnosed On ART",]
I4xM= Diagnosed_Men_melted[Diagnosed_Men_melted$variable == "Dropout",]







diag_off_art_auditF = data.frame(time = 2017, Lower = 195, Upper = 209, variable = "Women")
diag_off_art_auditM = data.frame(time = 2017, Lower = 100, Upper = 109, variable = "Men")



g10a = ggplot() + geom_line(data = I2xF, aes(x = time, y = value,  factor = run), colour = "black", alpha = 1) + theme_bw() +
  labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), max(time)+1, 4) ) +
  scale_y_continuous(breaks = seq(min(Diagnosed$value), max(Diagnosed$value), 1000) ) +
  geom_vline(xintercept=2005, col = "black", alpha = 0.5, linetype = "dashed") +
  geom_text(aes(x=2005, label="CD4 < 200", y=-400))+
  geom_vline(xintercept=2012, col = "black", alpha = 0.5, linetype = "dashed") +
  geom_text(aes(x=2012, label="CD4 < 350", y=-800))+
  geom_vline(xintercept=2015, col = "black", alpha = 0.5, linetype = "dashed") +
  geom_text(aes(x=2015, label="CD4 < 500", y=-400))+
  geom_vline(xintercept=2016, col = "black", alpha = 0.5, linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20), legend.position = "none")+
  geom_errorbar(data = diag_off_art_auditF, size=1, aes(x = time, ymin = Lower, ymax = Upper), col ="#e41a1c") +
  ggtitle("CROSS:VALIDATION: Numbers of all women in Grand Cotonou diagnosed NOT on ART compared to data in red (Pre-ART audit report 2017)")
g10a 




g10b = ggplot() + geom_line(data = I2xM, aes(x = time, y = value, factor = run), colour = "black", alpha = 1) + theme_bw() +
  labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), max(time)+1, 4) ) +
  scale_y_continuous(breaks = seq(min(Diagnosed$value), max(Diagnosed$value), 1000) ) +
  geom_vline(xintercept=2005, col = "black", alpha = 0.5, linetype = "dashed") +
  geom_text(aes(x=2005, label="CD4 < 200", y=-400))+
  geom_vline(xintercept=2012, col = "black", alpha = 0.5, linetype = "dashed") +
  geom_text(aes(x=2012, label="CD4 < 350", y=-800))+
  geom_vline(xintercept=2015, col = "black", alpha = 0.5, linetype = "dashed") +
  geom_text(aes(x=2015, label="CD4 < 500", y=-400))+
  geom_vline(xintercept=2016, col = "black", alpha = 0.5, linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20), legend.position = "none")+
  geom_errorbar(data = diag_off_art_auditM, size=1, aes(x = time, ymin = Lower, ymax = Upper), col ="#e41a1c") +
  ggtitle("CROSS:VALIDATION: Numbers of all men in Grand Cotonou diagnosed NOT on ART compared to data in red (Pre-ART audit report 2017)")
g10b




I2xF_2017 = data.frame(
  # Women_2016=I2xF[I2xF$time == 2016, "value"],
  Women_2017=I2xF[I2xF$time == 2017, "value"],
  # Men_2016=I2xM[I2xM$time == 2016, "value"], 
  Men_2017=I2xM[I2xM$time == 2017, "value"]
  )
I2xF_2017_melted = reshape2::melt(I2xF_2017)
  
g10c = ggplot(I2xF_2017_melted, aes(x = variable, y=value)) + geom_boxplot() + theme_bw() +
  geom_errorbar(aes(x = "Women_2017", ymax = 201, ymin = 215), col = "red", width = I(0.2))+
  geom_errorbar(aes(x = "Men_2017", ymax = 108, ymin = 117), col = "red", width = I(0.2))+
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) + 
  labs(x = "Sex", y = "")+ ggtitle("CROSS VALIDATION: Across the difference runs, the number of women and men in Grand Cotonou \n diagnosed but NOT on ART in 2017 vs the data in red for 2017")
g10c

I3xF_2017 = data.frame(Women=I3xF[I3xF$time == 2017, "value"], Men=I3xM[I3xM$time == 2017, "value"])
I3xF_2017_melted = reshape2::melt(I3xF_2017)

I3xF_2017_melted = rbind(I3xF_2017_melted, data.frame(variable = c("Both"), value = rowSums(I3xF_2017)))



g10d = ggplot(I3xF_2017_melted, aes(x = variable, y=value)) + geom_boxplot() + theme_bw() +
  geom_errorbar(aes(x = "Women", ymax = 6155, ymin = 13050), col = "red", width = I(0.2))+
  geom_errorbar(aes(x = "Men", ymax = 2369, ymin = 5223), col = "red", width = I(0.2))+
  geom_errorbar(aes(x = "Both", ymin = 2369+6155, ymax = 5223+13050), col = "red", width = I(0.2))+
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) + 
  labs(x = "Sex", y = "")+
  ggtitle("FITTING: Number of women and men in Grand Cotonou diagnosed \n ON ART in 2017 vs the data in red \n FITTING TO TOTAL ON ART (MEN & WOMEN TOGETHER)")
g10d









I4xF_2015 = data.frame(Women=I4xF[I4xF$time == 2015, "value"], Men=I4xM[I4xM$time == 2015, "value"])
I4xF_2015_melted = data.frame(variable = "All", value = rowSums(I4xF_2015))



ART_sex_ratio = data.frame(time = 2017, Lower = 1.5, Upper = 3)



g12=ggplot() + geom_line(data = Diagnosed_Women_Men_ratio[Diagnosed_Women_Men_ratio$variable == "Diagnosed On ART",], aes(x = time, y = value, factor = run), alpha = 1) + 
  theme_bw() + labs(y = "")+
  ggtitle("FITTING: Ratio of Women:Men on ART compared to data (2017) in red from audit 2017")+
  scale_x_continuous(breaks = seq(min(time), (max(time)+1), 4) ) +
  scale_y_continuous(breaks = seq(1, 5, 0.2) )+
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) + 
  geom_point(aes(x = 2017, y = 2.5), col = "red", size = I(03))+
  # geom_hline(yintercept = c(1,2), col = "orange") + 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=24),
        legend.text=element_text(size=18))+
  geom_errorbar(aes(x = 2017, ymax = 3, ymin = 1.5), col = "red", width = I(1))
  #+ annotate("text", x = 1994, y = 2.5, size = 7,col = "orange")


g12


g10e = ggplot(I4xF_2015_melted, aes(x = variable, y=value)) + geom_boxplot() + theme_bw() +
  # geom_errorbar(aes(x = "Women", ymax = 6155, ymin = 13050), col = "red", width = I(0.2))+
  geom_errorbar(aes(x = "All", ymax = 56, ymin = 56), col = "red", width = I(0.2))+
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) + 
  labs(x = "Sex", y = "")+
  ggtitle("CROSS VALIDATION: Number of women and men in Grand Cotonou in dropout \n from ART in 2015 vs the data in red (Lower bound)")
g10e


g10f = ggplot() + geom_line(data = ART_inits_cumu_melted, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() +
  labs(y = "")+
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5))+
  ggtitle("Cumulative ART initiations of all groups combined")

g10f

g10g = ggplot() + geom_line(data = ART_inits_melted, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() +
  labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), max(time), 2) ) +
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) +
  geom_errorbar(aes(x = 2015, ymax = 2051, ymin = 1318), col = "red", width = I(0.5), size = I(02))+
  ggtitle("CROSS VALIDATION: Numbers of ART initiations each year of all groups combined vs. the data in red")
  
  
g10g




g10h=ggplot() + geom_line(data = ART_inits_ratio_F_over_M_melted, aes(x = time, y = value, factor = variable), alpha = 1) +
  theme_bw() + labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), (max(time)+3), 2) ) +
  scale_y_continuous(breaks = seq(1, 5, 0.2) )+
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5)) +
  geom_point(aes(x = 2015, y = 2.4), col = "red", size = I(03))+
  ggtitle("CROSS VALIDATION: Ratio of ART initiations women:Men, \n compared with data from audit 2017 in red")
g10h




g10i=ggplot() + geom_line(data = FtM_ratio_ART_initiation_rates_from_HIVpos_not_on_ART_melted, aes(x = time, y = value, factor = variable), alpha = 1) +
  theme_bw() + labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), (max(time)+3), 2) ) +
  scale_y_continuous(breaks = seq(1, 5, 0.2) )+
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5)) +
  geom_point(aes(x = 2015, y = 1.7), col = "red", size = I(03))+
  ggtitle("CROSS VALIDATION: Ratio of ART initiation rate from HIV+ not on ART Women:Men, \n compared with data from audit 2017 in red")
g10i

g10j = ggplot(ART_init_rate_from_all_HIV_pos_NOT_ON_ART_by_sex_melted, aes(x = variable, y = value)) + geom_boxplot() + theme_bw() +
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5)) +
  labs(y = "")+
  geom_point(aes(x = "Women_2015", y = 0.169), col = "red", size = 3)+
  geom_point(aes(x = "Men_2015", y = 0.1), col = "red", size = 3)+
  ggtitle("CROSS VALIDATION: ART initiation rate from HIV+ not on ART Women & Men in 2014 and 2015, \n compared with data for 2015 in red")
g10j


# Ratio female:male of rate of initiation per infected = 8.0/4.68 = 1.7




g10k = ggplot() + geom_line(data = ART_REinits_cumu_melted, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() +
  labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), max(time), 2) ) +
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) +
  # geom_errorbar(aes(x = 2015, ymax = 2051, ymin = 1318), col = "red", width = I(0.5), size = I(02))+
  ggtitle("Cumulative numbers of ART RE initiations each year of all groups combined")


g10k


g10l = ggplot() + geom_line(data = ART_REinits_melted, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() +
  labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), max(time), 2) ) +
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) +
  geom_point(aes(x = 2015, y = 213), size = 5, col = "red") +
  ggtitle("CROSS VALIDATION: Numbers of ART RE initiations each year of all groups combined vs. the data in red") + 
  geom_text(col = "red",aes(x=2002, label="Data from monitoring reports for 2015 \n for Littoral + Atlantique + Oueme", y=300), size = 8)
  


g10l




g10m = ggplot() + geom_line(data = frac_ART_inits_REinits_melted[frac_ART_inits_REinits_melted$time > 1999,], aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() +
  labs(y = "")+
  scale_x_continuous(breaks = seq(min(time), max(time), 2) ) +
  theme(text = element_text(size=20), legend.position = "none",plot.title = element_text(hjust = 0.5)) +
  geom_point(aes(x = 2013.5, y = 8), size = 5, col = "red") +
  geom_point(aes(x = 2014.5, y = 7), size = 5, col = "red") +
  geom_point(aes(x = 2015, y = 7), size = 5, col = "red") +
  geom_point(aes(x = 2015.5, y = 9), size = 5, col = "red") +
  geom_point(aes(x = 2016, y = 12), size = 5, col = "red") +
  geom_point(aes(x = 2017, y = 11), size = 5, col = "red") +
  ggtitle("Percent of ART initiations which are RE initiations (%) vs. the data in red")  +
  geom_text(col = "red",aes(x=2005, label="Data from monitoring reports by semester for BENIN", y=30), size = 8)

g10m
# ggplot() + geom_line(data = I2x, aes(x = time, y = value, factor = variable, factor = run, colour = variable), alpha = 1) + theme_bw() + labs(y = "Numbers of all people in Grand Cotonou")+
#   scale_x_continuous(breaks = seq(min(time), max(time)+1, 4) ) +
#   scale_y_continuous(breaks = seq(min(Diagnosed$value), max(Diagnosed$value), 1000) ) +
#   geom_vline(xintercept=2005, col = "black", alpha = 0.5, linetype = "dashed") +
#   geom_text(aes(x=2005, label="CD4 < 200", y=-400))+
#   geom_vline(xintercept=2012, col = "black", alpha = 0.5, linetype = "dashed") +
#   geom_text(aes(x=2012, label="CD4 < 350", y=-800))+
#   geom_vline(xintercept=2015, col = "black", alpha = 0.5, linetype = "dashed") +
#   geom_text(aes(x=2015, label="CD4 < 500", y=-400))+
#   geom_vline(xintercept=2016, col = "black", alpha = 0.5, linetype = "dashed")




# 
# 
# ggplot() + geom_line(data = Diagnosed_Women_Men, aes(x = time, y = value, factor = variable, factor = run, colour = variable, linetype = group), alpha = 1) + theme_bw() + labs(y = "Numbers of all people in Grand Cotonou")+
#   scale_x_continuous(breaks = seq(min(time), max(time), 4) ) +
#   scale_y_continuous(breaks = seq(min(Diagnosed$value), max(Diagnosed$value), 1000) )
# 
# ggplot() + geom_line(data = Diagnosed_Women_Men_ratio, aes(x = time, y = value, factor = variable, factor = run, colour = variable), alpha = 1) + 
#   theme_bw() + labs(y = "Ratio of Women:Men")+
#   scale_x_continuous(breaks = seq(min(time), max(time), 4) ) +
#   scale_y_continuous(breaks = seq(1, 5, 0.2) )

# ggplot() + geom_line(data = Diagnosed, aes(x = time, y = value, factor = variable, factor = run, colour = group), alpha = 1) + theme_bw() + facet_wrap(~variable, scales = "fixed") + labs(y = "Numbers of all people in Grand Cotonou") +
#   scale_x_continuous(breaks = seq(min(time), max(time), 4) ) +
#   scale_y_continuous(breaks = seq(min(Diagnosed$value), max(Diagnosed$value), 1000) )
# 
# #

g11=ggplot() + geom_line(data = Diagnosed_Women_Men_ratio[Diagnosed_Women_Men_ratio$variable == "Diagnosed Off ART",], aes(x = time, y = value, factor = run, colour = variable), alpha = 1) +
  theme_bw() + labs(y = "Ratio of Women:Men")+
  scale_x_continuous(breaks = seq(min(time), max(time), 4) ) +
  scale_y_continuous(breaks = seq(1, 5, 0.2) )+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20))



# cross validation to ratio female:male on ART from audit report 2017


HIV_deaths_aidsinfo = data.frame(time = seq(1990, 2016),
                                 lower = c(100,100,100,100,100,100,100,100,200,200,200,500,500,500,500,500,500,500,500,200,200,200,200,200,200,200,100),
                                 upper = c(300,300,300,300,300,400,400,500,600,900,1200,1500,1500,1500,1500,2000,2000,2000,2000,900,900,800,1300,1300,1300,1300,700))

g13=ggplot() + geom_line(data = HIV_deaths_melted, aes(x = time, y = value, factor = variable), alpha = 1) + theme_bw() + 
  labs(y = "HIV deaths") +
  ggtitle("CROSS VALIDATION: HIV deaths") +
  geom_errorbar(data = HIV_deaths_aidsinfo, aes(x = time, ymin = lower, ymax = upper), col = "orange")+
  theme(text = element_text(size=20)) +
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=24),
        legend.text=element_text(size=18))+ annotate("text", x = 1992, y = 1700, size = 7,col = "orange",
                                                                                           label = "Cross validation")

g13



# EQS 2012 mapping 889 – 1391
EQS_mapping = data.frame(time = 2012, lower = 889, upper = 1391)
g14=ggplot() + geom_line(data = N_Pro_FSW_melted, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() + 
  labs(y = "Number of Pro FSW") +
  ggtitle("FITTING: Number of Pro FSW") +
  theme(text = element_text(size=20)) + 
  geom_errorbar(data = EQS_mapping, aes(x = time, ymin = lower, ymax = upper), col = "red")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=24),
        legend.text=element_text(size=18))+ annotate("text", x = 1992, y = 1200, size = 7,col = "red",
                                                     label = "Fitting to EQS mapping,\n  CI given by Michel")


g14

N_low_level = data.frame(time = 2012, lower = 2000, upper = 3000)

g15=ggplot() + geom_line(data = N_Low_FSW_melted, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() + 
  labs(y = "Number of Low level FSW") +
  ggtitle("CROSS VALIDATION: Number of Low level FSW") +
  theme(text = element_text(size=20)) + 
  geom_errorbar(data = N_low_level, aes(x = time, ymin = lower, ymax = upper), col = "orange")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=24),
        legend.text=element_text(size=18))+ annotate("text", x = 1992, y = 2750, size = 7,col = "orange",
                                                     label = "Cross validating")
  

low_pro_ratio = N_Low_FSW_melted
low_pro_ratio$value = N_Low_FSW_melted$value / N_Pro_FSW_melted$value 

g22 = ggplot() + geom_line(data = low_pro_ratio, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() + 
  labs(y = "Ratio of low level FSW to pro FSW") +
  ggtitle("FITTING: Ratio of low level FSW to pro FSW") +
  theme(text = element_text(size=20)) + 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=24),
        legend.text=element_text(size=18))+ annotate("text", x = 1992, y = 4.5, size = 7,col = "red",
                                                     label = "Fitting") +
  geom_hline(yintercept = c(1, 5), col = "red")


# lowndes 2002 19970
N_clients = data.frame(time = 1998,
                       value = 19970)

g16=ggplot() + geom_line(data = N_Client_melted, aes(x = time, y = value, factor = run), alpha = 1) + theme_bw() + 
  ggtitle("CROSS VALIDATION: Number of clients (4 cities study)") +
  labs(y = "Number of Clients")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=24)) + 
  geom_point(data = N_clients, aes(x = time, y = value), col = "orange", size = 3)+
  annotate("text", x = 1992, y = 90000, size = 7,col = "orange",
           label = "Note that their study area (Cotonou) is smaller than ours (Grand Cotonou)")
g16







g17=ggplot() + geom_line(data = condom_Pro_FSW_melted, aes(x = time, y = value,  factor = run, colour = variable), alpha = 1) +
  theme_bw() + labs(y = "Condom use of Pro FSW")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20), legend.position = "top")

g18=ggplot() + geom_line(data = condom_Low_FSW_melted, aes(x = time, y = value,  factor = run, colour = variable), alpha = 1) +
  theme_bw() + labs(y = "Condom use of Low level FSW")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20), legend.position = "top")

g19=ggplot() + geom_line(data = condom_GPF_noncomm_melted, aes(x = time, y = value, factor = variable), alpha = 1) +
  theme_bw() + labs(y = "Condom use of GPF")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20))

# ggplot() + geom_line(data = condom_GPM_noncomm_melted, aes(x = time, y = value, factor = variable), alpha = 1) +
#   theme_bw() + labs(y = "Condom use of GPM")


#



g20=ggplot() + geom_line(data = annual_client_volume_pro_FSW_melted, aes(x = time, y = value, factor = variable), alpha = 1) +
  theme_bw() + labs(y = "Annual client volume for Pro FSW")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20))

g21=ggplot() + geom_line(data = pfFSW_melted, aes(x = time, y = value*100, factor = run), alpha = 1) +
  theme_bw() + labs(y = "HIV prevalence of incoming FSWs from neighbouring countries(%)")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20))

### VIRALLY SUPPRESSED!!

viral_supp_vec <- rep(unlist(lapply(res_best_runs, function(x) x$viral_supp[1,3])), each = length(time))

ART_indiv_melted_ALL <- ART_indiv_melted[ART_indiv_melted$variable == "All", ]


ART_indiv_melted_ALL$value <- ART_indiv_melted_ALL$value * viral_supp_vec

viral_supp_data_aidsinfo <- data.frame(time = c(2015, 2016, 2017),
                                       Lower = c(14, 15, 28),
                                       Upper = c(30, 32, 60),
                                       variable = "All")


g22 <- ggplot() + geom_line(data = ART_indiv_melted_ALL, aes(x = time, y = value*100,  factor = run), alpha = 0.3) + theme_bw()  + labs(y = "Viral suppression (%)") +
  ggtitle("CROSS-VALIDATION: People living with HIV who have suppressed viral loads (%) vs data from UNAIDS AIDSINFO")+
  geom_errorbar(data = viral_supp_data_aidsinfo, aes(x = time, ymin = Lower, ymax = Upper), col = "orange",size = I(1))+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=20)) + theme(legend.position = "top")+
  scale_x_continuous(breaks = seq(2000, 2018, 1), limits = c(2000,2018))
  



plotlist_paper = list(g1b_paper, )


plotlist = list(g1, g7, g8, g10d, g12, g3, g4b, g5,g5b,  g9, g6, 
                g14, g22,g2,
                # g10a, g10b, g10f, 
                g10, g10c, g10e, 
                g10g, g10h, g10i, g10j, 
                  g13, g15, g16, g17, g18, g19, g20, g21 ,g4, g22)

# outputfolder = "0103"

  




# for(i in 1:length(plotlist))
# 
# {
#   # jpeg(filename=paste0("C:\\Users\\eg1012\\Google Drive\\Imperial\\Cotonou\\Models\\Transmission Model\\Outputs/",
#   #                      batch_folder, "/","all", i,".jpg"), width = 1200, height = 700)
#   print(plotlist[i][[1]])
#   # dev.off()
# 
# 
#   # save(plotlist[i][[1]][[1]], file = paste0("C:\\Users\\eg1012\\Google Drive\\Imperial\\Cotonou\\Models\\Transmission Model\\Outputs/",
#   #                                      batch_folder, "/","test_obj.Rda", i))
# 
# }








# par_table = data.frame(do.call(rbind, lapply(all_good_pars, function(x) {
#   return(x$rate_leave_pro_FSW)
# })),
# 
# do.call(rbind, lapply(all_good_pars, function(x) {
#   return(x$prev_non_ben_fsw_1993)
# })),
# 
# do.call(rbind, lapply(all_good_pars, function(x) {
#   return(x$prev_non_ben_fsw_2015)
# })),
# 
# do.call(rbind, lapply(all_good_pars, function(x) {
#   return(x$fraction_FSW_foreign)
# })),
# 
# do.call(rbind, lapply(all_good_runs, function(x) {
#   return(x$lambda_sum_0[which(time == 1993),1])
# })),
# do.call(rbind, lapply(all_good_runs, function(x) {
#   return(x$lambda_sum_0[which(time == 2012),1])
# })),
# do.call(rbind, lapply(all_good_runs, function(x) {
#   return(x$lambda_sum_0[which(time == 2015),1])
# })),
# do.call(rbind, lapply(all_good_runs, function(x) {
#   return(x$prev[which(time == 1993),1])
# })),
# do.call(rbind, lapply(all_good_runs, function(x) {
#   return(x$prev[which(time == 2012),1])
# })),
# do.call(rbind, lapply(all_good_runs, function(x) {
#   return(x$prev[which(time == 2015),1])
# }))
# 
# 
# )
# colnames(par_table) <- c("rate of leaving sex work",
#                          "prevalence of incoming professional FSW 1993",
#                          "prevalence of incoming professional FSW 2015",
#                          "fraction of FSW that are foreign",
#                          
#                          "incidence professional FSW 1993",
#                          "incidence professional FSW 2012",
#                          "incidence professional FSW 2015",
#                          
#                          
#                          "prevalence professional FSW 1993",
#                          "prevalence professional FSW 2012",
#                          
#                          "prevalence professional FSW 2015"
# )
# 
# par_table <- cbind(run = c(1: length(all_good_runs)), par_table)
# 
# write.table(par_table,sep = ",", row.names = F,file = paste0("C:\\Users\\eg1012\\Google Drive\\Imperial\\Cotonou\\Models\\Transmission Model\\Outputs/",
#                                                              batch_folder, "/", "pars_and_outcomes.csv"))
# 
# 
# 
# 
# # Analysing the fits to ART:
# 
# # How does the number of FSW in 2012 affect how well each simulation fits to ART coverage (%)?
# 
# cor.test(N_Pro_FSW_melted[N_Pro_FSW_melted$time == 2012, "value"],
#          ART_indiv_melted[ART_indiv_melted$time == 2017 & ART_indiv_melted$variable == "Pro FSW", "value"])
# 
# plot(N_Pro_FSW_melted[N_Pro_FSW_melted$time == 2012, "value"],
#          ART_indiv_melted[ART_indiv_melted$time == 2017 & ART_indiv_melted$variable == "Pro FSW", "value"],
#      ylab = "FSW ART Coverage 2017", xlab = "N FSW in 2012")
# 
# 
# # 2.	What are the posteriors generally for the testing rate, 
# # ART initiation rate, dropout rate, reinitiation rate?
# 
# # FSW
# testing_rate_FSW = data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x) {
#   x$tau[,1]#[which(time == 2015)]
# }))))
# 
# # GP
# testing_rate_men = data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x) {
#   x$tau[,5]#[which(time == 2015)]
# }))))
# 
# testing_rate_women = data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x) {
#   x$tau[,3]#[which(time == 2015)]
# }))))
# 
# # data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x) {
# #   tot = x$rho * x$ART_eligible_CD4_above_500  # *x$above_500_by_group[i] # NOTE THIS IS ALWAYS 1 BEFORE INTERVENTION
# #   return(tot[,1])
# #   }))))
# # 
# # data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x) {
# #   tot = x$rho * x$ART_eligible_CD4_350_500 
# #   return(tot[,1])
# # }))))
# # 
# # data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x) {
# #   tot = x$rho * x$ART_eligible_CD4_200_349 
# #   return(tot[,1])
# # }))))
# # 
# # data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x) {
# #   tot = x$rho * x$ART_eligible_CD4_below_200
# #   return(tot[,1])
# # }))))
# 
# 
# # weighted average
# # FSW
# average_ART_init_rate_FSW = data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x) {
#   tot = x$rho * x$ART_eligible_CD4_above_500 * x$prop_FSW_I0_2/(x$prop_FSW_I0_2 + x$prop_FSW_I0_3 + x$prop_FSW_I0_4 + x$prop_FSW_I0_5)+
#     x$rho * x$ART_eligible_CD4_350_500 * x$prop_FSW_I0_3/(x$prop_FSW_I0_2 + x$prop_FSW_I0_3 + x$prop_FSW_I0_4 + x$prop_FSW_I0_5)+
#     x$rho * x$ART_eligible_CD4_200_349 * x$prop_FSW_I0_4/(x$prop_FSW_I0_2 + x$prop_FSW_I0_3 + x$prop_FSW_I0_4 + x$prop_FSW_I0_5)+
#     x$rho * x$ART_eligible_CD4_below_200 * x$prop_FSW_I0_5/(x$prop_FSW_I0_2 + x$prop_FSW_I0_3 + x$prop_FSW_I0_4 + x$prop_FSW_I0_5)
#   return(tot[,1])
# }))))
# 
# # NOTE THAT FROM NEW RUNS THERE WILL BE DIFFERENCE BETWEEN MEN AND WOMEN
# average_ART_init_rate_men_and_women = data.frame(time, t(do.call(rbind, lapply(res_best_runs, function(x) {
#   tot = x$rho * x$ART_eligible_CD4_above_500 * x$prop_FSW_I0_2/(x$prop_FSW_I0_2 + x$prop_FSW_I0_3 + x$prop_FSW_I0_4 + x$prop_FSW_I0_5)+
#     x$rho * x$ART_eligible_CD4_350_500 * x$prop_FSW_I0_3/(x$prop_FSW_I0_2 + x$prop_FSW_I0_3 + x$prop_FSW_I0_4 + x$prop_FSW_I0_5)+
#     x$rho * x$ART_eligible_CD4_200_349 * x$prop_FSW_I0_4/(x$prop_FSW_I0_2 + x$prop_FSW_I0_3 + x$prop_FSW_I0_4 + x$prop_FSW_I0_5)+
#     x$rho * x$ART_eligible_CD4_below_200 * x$prop_FSW_I0_5/(x$prop_FSW_I0_2 + x$prop_FSW_I0_3 + x$prop_FSW_I0_4 + x$prop_FSW_I0_5)
#   return(tot[,2])
# }))))
# 
# 
# 
# # FSW
# dropout_rate_FSW = data.frame(time, t(do.call(rbind, lapply(all_good_pars, function(x) {
#   rep(x$phi2[1], length(time))
# }))))
# 
# # GP
# dropout_rate_men_and_women = data.frame(time, t(do.call(rbind, lapply(all_good_pars, function(x) {
#   rep(x$phi2[2], length(time))
# }))))
# 
# # FSW
# reinit_rate_FSW = data.frame(time, t(do.call(rbind, lapply(all_good_pars, function(x) {
#   rep(x$iota[1], length(time))
# }))))
# 
# # GP
# reinit_rate_men_and_women = data.frame(time, t(do.call(rbind, lapply(all_good_pars, function(x) {
#   rep(x$iota[2], length(time))
# }))))
# 
# 
# 
# 
# 
# FSW_ART_rates = data.frame(parameter = c(rep("testing rate", length(time)),
#                          rep("average ART initiation rate", length(time)),
#                          rep("dropout rate", length(time)),
#                          rep("ART reinitiation rate", length(time))), 
#            group = "Pro FSW",
#            rbind(testing_rate_FSW,
#                  average_ART_init_rate_FSW,
#                  dropout_rate_FSW,
#                  reinit_rate_FSW))
# 
# FSW_ART_rates_melted = reshape2::melt(FSW_ART_rates, id.vars = c("parameter", "group", "time"))
# 
# FSW_ART_rates_melted$parameter = factor(FSW_ART_rates_melted$parameter, levels = c("testing rate",
#                                                                                    "average ART initiation rate",
#                                                                                    "dropout rate",
#                                                                                    "ART reinitiation rate"))
# 
# 
# ggplot(FSW_ART_rates_melted[FSW_ART_rates_melted$time > 1999,]) + geom_line(aes(x = time, y = value, factor = variable, colour = parameter), size = 1.5) + 
#   theme_bw()+ scale_colour_brewer(type = "qual", palette = 2) +
#   ggtitle("FSW parameters relevant to ART") +   theme(text = element_text(size=20), legend.position = "top",plot.title = element_text(hjust = 0.5))
# 
# 
# 
# 
# 
# 
# GP_ART_rates = data.frame(parameter = c(rep("testing rate", length(time)),
#                                         rep("testing rate", length(time)),
#                                         rep("average ART initiation rate", length(time)),
#                                         rep("average ART initiation rate", length(time)),
#                                         rep("dropout rate", length(time)),
#                                         rep("dropout rate", length(time)),
#                                         rep("ART reinitiation rate", length(time)), 
#                                         rep("ART reinitiation rate", length(time))), 
# 
#                           group = c(rep("Men", length(time)), 
#                                     rep("Women", length(time)),
#                                     
#                                     rep("Men", length(time)),
#                                     rep("Women", length(time)),
#                                     
#                                     rep("Men", length(time)),
#                                     rep("Women", length(time)),
#                                     
#                                     rep("Men", length(time)),
#                                     rep("Women", length(time))),
#                           
#                           rbind(testing_rate_men,
#                                 testing_rate_women,
#                                 
#                                 average_ART_init_rate_men_and_women,
#                                 average_ART_init_rate_men_and_women,
#                                 
#                                 dropout_rate_men_and_women,
#                                 dropout_rate_men_and_women,
#                                 
#                                 reinit_rate_men_and_women,
#                                 reinit_rate_men_and_women))
# 
# GP_ART_rates_melted = reshape2::melt(GP_ART_rates, id.vars = c("parameter", "group", "time"))
# 
# GP_ART_rates_melted$parameter = factor(GP_ART_rates_melted$parameter, levels = c("testing rate",
#                                                                                    "average ART initiation rate",
#                                                                                    "dropout rate",
#                                                                                    "ART reinitiation rate"))
# 
# 
# ggplot(GP_ART_rates_melted[GP_ART_rates_melted$time > 1999 & GP_ART_rates_melted$parameter != "average ART initiation rate" ,]) + geom_line(aes(x = time, y = value, factor = variable, colour = parameter), size = 1.5) + 
#   facet_wrap(~group)+
#   theme_bw()+ scale_colour_brewer(type = "qual", palette = 2) +
#   ggtitle("Whole population parameters relevant to ART") +   theme(text = element_text(size=20), legend.position = "top", plot.title = element_text(hjust = 0.5))
# 
# 
# 
# lapply(res_best_runs, function(x) {
#   # tot = x$rho * x$ART_eligible_CD4_above_500 + # *x$above_500_by_group[i] # NOTE THIS IS ALWAYS 1 BEFORE INTERVENTION
#     # x$rho * x$ART_eligible_CD4_350_500 +
#     # x$rho * x$ART_eligible_CD4_200_349 +
#   x$ART_eligible_CD4_above_500# * x$ART_eligible_CD4_below_200
# })
# 
# 
# lapply(all_good_pars, function(x) {
#   # tot = x$rho * x$ART_eligible_CD4_above_500 + # *x$above_500_by_group[i] # NOTE THIS IS ALWAYS 1 BEFORE INTERVENTION
#   # x$rho * x$ART_eligible_CD4_350_500 +
#   # x$rho * x$ART_eligible_CD4_200_349 +
#  x$phi5[2]# * x$ART_eligible_CD4_below_200
# })

# prevalence 2016
# range(prev_indiv_melted[prev_indiv_melted$time == 2016 & prev_indiv_melted$variable == "Pro FSW", "value"])
# 
# 
# range(Diagnosed_FSW_melted[Diagnosed_FSW_melted$group == "FSW" & Diagnosed_FSW_melted$variable == "Diagnosed On ART" & Diagnosed_FSW_melted$time == 2016, "value"])
# range(Diagnosed_FSW_melted[Diagnosed_FSW_melted$group == "FSW" & Diagnosed_FSW_melted$variable == "Diagnosed On ART" & Diagnosed_FSW_melted$time == 2017, "value"])
# 
# 
# 
# range(lapply(res_best_runs, function(x) sum(x$N[which(time == 2016),c(1:8)])))
# range(lapply(res_best_runs, function(x) sum(x$N[which(time == 2017),c(1:8)])))
# range(lapply(all_good_pars, function(x) x$fraction_F))
# range(lapply(all_good_pars, function(x) x$frac_women_ProFSW))
# range(lapply(res_best_runs, function(x) x$prev[which(time == 2016), 1]))
# range(lapply(res_best_runs, function(x) x$prev[which(time == 2017), 1]))


