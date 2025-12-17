# list of functions for plotting data

plot_effect_estimates = function(effect_ests, #list of stan outputs
                                 plot_models, # indices of models to plot in list
                                 my_pch=1,
                                 mod_cols = NULL,
                                 study_threshold){
  
  if (length(my_pch)==1) my_pch = (1:length(plot_models))+15
  if (length(my_pch)!=length(plot_models)) stop('length of my_pch needs to be the same as the number of input models')
  if (is.null(mod_cols)){
    mod_cols = brewer.pal(n = length(plot_models), name = 'Set1')
    if(length(plot_models)>9) writeLines('too many models - supply user defined colors')
  }
  if(length(mod_cols)!=length(plot_models)) stop('number of colors needs to be equal to number of models to be plotted')
  
  K_treatments = nrow(effect_ests[[plot_models[1]]])
  
  xlims = (exp(range(c(0, range(sapply(effect_ests[plot_models], rbind)))) )-1)*100
  x_points = pretty((xlims),6)
  plot(NA, NA, xlim = range(x_points),
       ylim = c(0.75,K_treatments+.25),
       panel.first=grid(), ylab='', yaxt='n', type='n',
       xlab = 'Change in rate of clearance (%)',
       xaxt = 'n')
  
  index_p = rev(seq(-.2,.2, length.out = length(plot_models)))
  abline(v=0,lwd=2)
  polygon(c(-1000, 100*(study_threshold-1), 100*(study_threshold-1), -1000),
          c(-100, -100, 100, 100), border = NA,
          col = adjustcolor('grey',.4))
  
  sort_ind = order(effect_ests[[plot_models[1]]][,'50%'])
  for(i in 1:length(plot_models)){
    points((exp(effect_ests[[plot_models[i]]][sort_ind,'50%'])-1)*100,
           1:K_treatments+index_p[i],pch=my_pch[i],
           col=mod_cols[i],cex=1.5)
    for(j in 1:length(sort_ind)){
      kk = sort_ind[j]
      lines((exp(c(effect_ests[[plot_models[i]]][kk,'2.5%'],
                   effect_ests[[plot_models[i]]][kk,'97.5%']))-1)*100,
            rep(j+index_p[i],2),col=mod_cols[i],lwd=1)
      lines((exp(c(effect_ests[[plot_models[i]]][kk,'10%'],
                   effect_ests[[plot_models[i]]][kk,'90%']))-1)*100,
            rep(j+index_p[i],2),col=mod_cols[i],lwd=3)
    }
  }
  axis(1, at = x_points)
  axis(2, at = 1:K_treatments, labels = rownames(effect_ests[[plot_models[1]]])[sort_ind])
  
}



plot_baseline_data = function(input_data){
  
  baseline_ind = input_data$Timepoint_ID==0
  bvl = aggregate(log10_viral_load ~ ID, input_data[baseline_ind, ], median)
  
  hist(bvl$log10_viral_load,
       breaks = seq(1,9,by=.5),
       xlab='Baseline viral load (SARS CoV2 genomes/mL)',
       ylab ='Number of patients',xlim=c(1,9),
       main='', xaxt ='n')
  axis(1, at = c(2,4,6,8), labels = c(expression(10^2),
                                      expression(10^4),
                                      expression(10^6),
                                      expression(10^8)))
  grid(); 
}


plot_serial_data = function(xx, xlims=c(0,7)){
  
  daily_VL_data = xx %>% group_by(ID, Time) %>%
    mutate(daily_VL = mean(log10_viral_load,na.rm = T))
  
  summary_VL_data = daily_VL_data %>% group_by(Timepoint_ID, Trt) %>%
    summarise(daily_VL = mean(daily_VL,na.rm = T),
              trt_color=unique(trt_color))
  
  summary_dat = daily_VL_data %>% ungroup() %>% distinct(ID, .keep_all = T) %>%
    group_by(Trt) %>%
    summarise(n = n(),
              trt_color = unique(trt_color)) %>% ungroup() %>%
    mutate(legend = paste0(Trt, ' (n=',n,')'))
  
  par(las=1)
  plot(summary_VL_data$Timepoint_ID, summary_VL_data$daily_VL,
       ylab = 'SARS CoV2 genomes/mL', panel.first=grid(),
       xlab = 'Time since randomization (days)',
       xlim = xlims, yaxt='n',type='n',
       ylim = range(xx$log10_viral_load))
  axis(2, at = c(2,4,6), labels = c(expression(10^2),
                                    expression(10^4),
                                    expression(10^6)))
  
  points(daily_VL_data$Time, daily_VL_data$daily_VL,
         col= adjustcolor(daily_VL_data$trt_color,.3),
         pch = 1)
  for(tt in unique(summary_VL_data$Trt)){
    ind = summary_VL_data$Trt==tt
    lines(summary_VL_data$Timepoint_ID[ind],
          summary_VL_data$daily_VL[ind], pch=15,
          type='b',col = summary_VL_data$trt_color[ind],lwd=3)
  }
  
  legend('topright', col = summary_dat$trt_color, 
         legend = summary_dat$legend,
         lwd=2, pch=15, cex=1, inset = 0.03)
}


bayes_R2 = function(mod_preds, mod_residuals) {
  var_pred = apply(mod_preds, 1, var)
  var_res = apply(mod_residuals, 1, var)
  var_pred / (var_pred + var_res)
}



make_stan_inputs = function(input_data_fit, 
                            int_covs_base,
                            int_covs_full,
                            slope_covs_base,
                            slope_covs_full,
                            trt_frmla,
                            Dmax
){
  
  ## check censored values come last
  if(!all(diff(input_data_fit$log10_viral_load == input_data_fit$log10_cens_vl)>=0)) stop()
  ind_dup = !duplicated(input_data_fit$ID)
  
  for(ll in unique(input_data_fit$Lab)){
    ind = input_data_fit$Lab==ll
    med_val = median(input_data_fit$CT_RNaseP[ind], na.rm = T)
    input_data_fit$CT_RNaseP[ind & is.na(input_data_fit$CT_RNaseP)]=med_val
    input_data_fit$CT_RNaseP[ind] = input_data_fit$CT_RNaseP[ind]-med_val
  }
  input_data_fit$RnaseP_scaled = input_data_fit$CT_RNaseP
  
  # make the covariate matrix
  # check no missing data
  if(!all(!apply(input_data_fit[, union(int_covs_full,slope_covs_full)], 2, function(x) any(is.na(x))))){
    stop('Missing data in covariate matrix!')
  }
  
  ind_contr = which(apply(input_data_fit[, int_covs_base,drop=F], 2, 
                          function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_intcpt_1 = model.matrix( ~ ., 
                               data = input_data_fit[, int_covs_base[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_intcpt_1 = array(dim = c(nrow(input_data_fit),0))
  }
  
  ind_contr = which(apply(input_data_fit[, int_covs_full,drop=F], 2, function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_intcpt_2 = model.matrix( ~ ., 
                               data = input_data_fit[, int_covs_full[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_intcpt_2 = array(dim = c(nrow(input_data_fit),0))
  }
  
  ind_contr = which(apply(input_data_fit[, slope_covs_base,drop=F], 2, function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_slope_1 = model.matrix( ~ ., 
                              data = input_data_fit[, slope_covs_base[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_slope_1 = array(dim = c(nrow(input_data_fit),0))
  }
  
  ind_contr = which(apply(input_data_fit[, slope_covs_full,drop=F], 2, function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_slope_2 = model.matrix( ~ ., 
                              data = input_data_fit[, slope_covs_full[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_slope_2 = array(dim = c(nrow(input_data_fit),0))
  }
  
  if(!nrow(X_intcpt_1) == nrow(input_data_fit)) stop()
  if(!nrow(X_slope_1) == nrow(input_data_fit)) stop()
  
  cov_matrices = list(X_int=list(X_intcpt_1, X_intcpt_2),
                      X_slope=list(X_slope_1, X_slope_2))
  
  ID_map = data.frame(ID_key = input_data_fit$ID,
                      ID_stan = as.numeric(as.factor(input_data_fit$ID)))
  writeLines(sprintf('There are a total of %s patients in the database with a total of %s PCRs analysable',
                     max(ID_map$ID_stan),
                     nrow(input_data_fit)))
  
  ind_cens = !input_data_fit$log10_viral_load>
    input_data_fit$log10_cens_vl
  
  writeLines(sprintf('%s%% of samples are below LOD',
                     round(100*mean(ind_cens),digits = 2)))
  
  analysis_data_stan = list(Ntot = nrow(input_data_fit),
                            N_obs = sum(!ind_cens),
                            n_id = max(ID_map$ID_stan),
                            id = ID_map$ID_stan,
                            ind_start = which(!duplicated(ID_map$ID_stan)),
                            obs_day = input_data_fit$Time,
                            log_10_vl = input_data_fit$log10_viral_load,
                            log10_cens_vl = input_data_fit$log10_cens_vl,
                            RNaseP = input_data_fit$RnaseP_scaled,
                            Time_max = Dmax)
  
  ID_map = ID_map[!duplicated(ID_map$ID_key), ]
  
  writeLines('check stan data formatting:')
  all(analysis_data_stan$log_10_vl[1:analysis_data_stan$N_obs]>
        analysis_data_stan$log10_cens_vl[1:analysis_data_stan$N_obs]) &
    all(analysis_data_stan$log_10_vl[(1+analysis_data_stan$N_obs):analysis_data_stan$Ntot] ==
          analysis_data_stan$log10_cens_vl[(1+analysis_data_stan$N_obs):analysis_data_stan$Ntot])
  
  Treatment_matrix = model.matrix(trt_frmla, data = input_data_fit)
  Treatment_matrix[,1]=0 # first column is dummy
  
  analysis_inputs = list(cov_matrices=cov_matrices,
                         analysis_data_stan=analysis_data_stan,
                         Treatment_matrix=Treatment_matrix,
                         ID_map=ID_map)
  return(analysis_inputs)
}






get_rates_linear_mod = function(mod_out, # single model fit - not a list
                                analysis_data_stan){
  
  # get the indices of first datapoint for each individual
  ind_id = which(!duplicated(analysis_data_stan$id))
  thetas = 
    rstan::extract(mod_out, 
                   pars = c('beta_0','beta_cov','theta_rand_id','trt_effect'))
  beta_cov = x_slope*slope_coefs;
  
  beta_0*exp(trt_slope[i]+theta_rand_id[id[i]][2]+beta_cov[i])
}

# plot_data_model_fits = function(mod_out, # model fits
#                              models_plot, # which models to plot
#                              K_plots,
#                              mod_cols,
#                              ID_map,
#                              analysis_data_stan
# ){
#   
#   # extract posterior parameters and outputs
#   thetas = list()
#   for(mm in 1:length(mod_out)){
#     thetas[[mm]] = rstan::extract(mod_out[[mm]])
#   }
#   counter = 1
#   
#   ID_map$Treatment = gsub(pattern = '\n',
#                     replacement = '',
#                     x = ID_map$Treatment,fixed = T)
#   while(counter <= nrow(ID_map)){
#     
#     # draw individual model fit with data
#     ind = analysis_data_stan$id==ID_map$ID_stan[counter]
#     plot(analysis_data_stan$obs_day[ind],
#          analysis_data_stan$log_10_vl[ind],
#          xlab='', ylab='', 
#          xaxt='n', yaxt='n',
#          panel.first=grid(), xlim=c(0,7),
#          ylim = range(analysis_data_stan$log_10_vl))
#     # if(counter %% sqrt(K_plots) == 1){
#     #   mtext(text = 'RNA copies per mL',side = 2,
#     #         line = 3,las = 3)
#     # }
#     axis(1, at = c(0,3,7))
#     axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
#                                         expression(10^4),
#                                         expression(10^6),
#                                         expression(10^8)))
#     # if((counter%%K_plots) >= K_plots - sqrt(K_plots)){
#     #   mtext(text = 'Days',side = 1,line = 2)
#     # }
#     for(mm in models_plot){
#       ix = order(analysis_data_stan$obs_day[ind])
#       my_xs = analysis_data_stan$obs_day[ind][ix]
#       polygon(x = c(my_xs, rev(my_xs)),
#               y = c(apply(thetas[[mm]]$preds[,ind],2,
#                           quantile,probs=0.025)[ix],
#                     rev(apply(thetas[[mm]]$preds[,ind],2,
#                               quantile,probs=0.975)[ix])),
#               border = NA, 
#               col = adjustcolor(mod_cols[mm],alpha.f = .3))
#       lines(my_xs,
#             colMeans(thetas[[mm]]$preds[,ind])[ix],
#             col = mod_cols[mm],lwd=2)
#     }
#     points(analysis_data_stan$obs_day[ind],
#            analysis_data_stan$log_10_vl[ind],pch=16)
#     
#     mtext(text = paste0(ID_map$ID[counter],
#                         '\n',
#                         ID_map$Treatment[counter]),
#           side = 3, line = -0.5, cex=0.8)
#     counter=counter+1
#   }
#   
# }

plot_data_model_fits = 
  function(model_list, # list of model fits
           models_to_plot, # which models to plot
           mod_cols,
           ID_map,
           analysis_data_stan
  ){
    
    # extract posterior parameters and outputs
    thetas = list()
    for(mm in 1:length(model_list)){
      thetas[[mm]] = rstan::extract(model_list[[mm]], pars='preds')
    }
    counter = 1
    
    # ID_map$Treatment = gsub(pattern = '\n',
    #                   replacement = '',
    #                   x = ID_map$Treatment,fixed = T)
    xlims = range(analysis_data_stan$obs_day)
    while(counter <= nrow(ID_map)){
      
      # draw individual model fit with data
      ind = analysis_data_stan$id==ID_map$ID_stan[counter]
      cens = as.numeric(analysis_data_stan$log10_cens_vl[ind]==analysis_data_stan$log_10_vl[ind])
      plot(analysis_data_stan$obs_day[ind],
           analysis_data_stan$log_10_vl[ind],
           xlab='', ylab='', 
           xaxt='n', yaxt='n',
           panel.first=grid(), xlim=xlims,
           ylim = range(analysis_data_stan$log_10_vl),pch=cens+1)
      
      axis(1, at = c(0,3,5))
      axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
                                          expression(10^4),
                                          expression(10^6),
                                          expression(10^8)))
      
      for(mm in models_to_plot){
        ix = order(analysis_data_stan$obs_day[ind])
        my_xs = analysis_data_stan$obs_day[ind][ix]
        polygon(x = c(my_xs, rev(my_xs)),
                y = c(apply(thetas[[mm]]$preds[,ind],2,
                            quantile,probs=0.025)[ix],
                      rev(apply(thetas[[mm]]$preds[,ind],2,
                                quantile,probs=0.975)[ix])),
                border = NA, 
                col = adjustcolor(mod_cols[mm],alpha.f = .3))
        lines(my_xs,
              colMeans(thetas[[mm]]$preds[,ind])[ix],
              col = mod_cols[mm],lwd=2)
      }
      points(analysis_data_stan$obs_day[ind],
             analysis_data_stan$log_10_vl[ind],pch=cens+16)
      
      mtext(text = ID_map$ID_key[counter],
            side = 3, line = 0.5, cex=0.7)
      counter=counter+1
    }
    
  }



plot_individ_curves = function(xx, IDs, xlims, mITT_threshold){
  
  if(!all(IDs %in% xx$ID)) stop('missing IDs!!')
  ylims = range(xx$log10_viral_load)
  xx$Day = round(xx$Time)
  
  xx = xx%>% group_by(ID, Day) %>%
    mutate(Daily_VL = mean(log10_viral_load))
  
  for(id in unique(IDs)){
    xx_id = xx%>% filter(ID==id)
    plot(xx_id$Time,
         xx_id$log10_viral_load,
         pch = as.numeric(xx_id$log10_cens_vl==xx_id$log10_viral_load)+1,
         xlab='', ylab='', 
         xaxt='n', yaxt='n',
         panel.first=grid(), xlim=xlims,
         ylim = ylims)
    title(id)
    points(xx_id$Time,
           xx_id$log10_viral_load,
           pch = as.numeric(xx_id$log10_cens_vl==xx_id$log10_viral_load)+16)
    xx_id = xx_id %>% distinct(Day, .keep_all = T)
    lines(xx_id$Day, xx_id$Daily_VL, type='b',pch='+',col='darkred')
    
    axis(1, at = c(0,7,14,21))
    axis(2, at = c(0,2,4,6,8), labels = c(1,
                                          expression(10^2),
                                          expression(10^4),
                                          expression(10^6),
                                          expression(10^8)))
    abline(h=log10(mITT_threshold), lty=2)
    
  }
  #legend('topright', legend = c('>LLOQ','<LLOQ'), pch = 16:17, inset = 0.03)
}


plot_coef_effects = function(stan_out, cov_mat, stan_inputs){
  
  thetas = rstan::extract(stan_out)
  alpha_coefs = apply(thetas$intercept_coefs,2,
                      quantile,probs=c(0.025,.1,.5,.9,0.975))
  xlims=range(alpha_coefs)
  
  cov_names_intercept =
    plyr::mapvalues(x = colnames(stan_inputs$cov_matrices$X_int[[cov_mat]]),
                    from = c('Age_scaled','Antibody_test',
                             'Symptom_onset','N_dose'),
                    to = c('Age','Serology RDT',
                           'Days since\nsymptom onset',
                           'Number of\nvaccine doses'))
  
  cov_names_slope =
    plyr::mapvalues(x = colnames(stan_inputs$cov_matrices$X_slope[[cov_mat]]),
                    from = c('Age_scaled','Antibody_test',
                             'Symptom_onset','N_dose'),
                    to = c('Age','Serology RDT',
                           'Days since\nsymptom onset',
                           'Number of\nvaccine doses'))
  
  plot(alpha_coefs['50%', ], 1:ncol(alpha_coefs),
       xlim=xlims,yaxt='n',ylab='',bty='n',xaxt='n',
       panel.first=grid(), xlab='Intercept (fold change)')
  abline(v=0,lty=2,lwd=2)
  for(i in 1:ncol(alpha_coefs)){
    lines(c(alpha_coefs['10%',i], alpha_coefs['90%',i]),
          c(i,i), lwd=3)
    lines(c(alpha_coefs['2.5%',i], alpha_coefs['97.5%',i]),
          c(i,i), lwd=1)
  }
  axis(2, at =1:ncol(alpha_coefs), labels = cov_names_intercept,tick = F)
  
  x_points = signif(10^seq(xlims[1], xlims[2],length.out = 5),2)
  axis(1, at = log10(x_points), labels = x_points)
  
  beta_coefs = apply(thetas$slope_coefs,2,quantile,
                     probs=c(0.025,.1,.5,.9,0.975))
  xlims=range(beta_coefs)
  plot(beta_coefs['50%', ], 1:ncol(beta_coefs),
       xlim=xlims,yaxt='n',ylab='',bty='n',xaxt='n',
       panel.first=grid(), xlab='Slope (multiplicative effect)')
  abline(v=0,lty=2,lwd=2)
  for(i in 1:ncol(beta_coefs)){
    lines(c(beta_coefs['10%',i], beta_coefs['90%',i]),
          c(i,i), lwd=3)
    lines(c(beta_coefs['2.5%',i], beta_coefs['97.5%',i]),
          c(i,i), lwd=1)
  }
  axis(2, at =1:ncol(beta_coefs), labels = cov_names_slope, tick = F)
  x_points = signif(10^seq(xlims[1], xlims[2],length.out = 5),2)
  axis(1, at = log10(x_points), labels = x_points)
}

# Checks a function for use of global variables
# Returns TRUE if ok, FALSE if globals were found.
checkStrict <- function(f, silent=FALSE) {
  vars <- codetools::findGlobals(f)
  found <- !vapply(vars, exists, logical(1), envir=as.environment(2))
  if (!silent && any(found)) {
    warning("global variables used: ", paste(names(found)[found], collapse=', '))
    return(invisible(FALSE))
  }
  
  !any(found)
}

calculate_fever_clearance = function(temp_dat,
                                     window_clear = 24/24, # look ahead window to define "fever clearance"
                                     threshold=37){
  
  if(!'ax_temperature' %in% colnames(temp_dat)) stop('needs to contain a ax_temperature column')
  
  temp_dat$clearance_time = NA
  # For interval censored data, the status indicator is 0=right censored, 1=event at time, 2=left censored, 3=interval censored. 
  temp_dat$clearance_time_cens = 1
  
  temp_dat$fever_binary = temp_dat$ax_temperature>threshold
  temp_dat = dplyr::arrange(temp_dat, ID, Time) 
  temp_dat = temp_dat[!is.na(temp_dat$ax_temperature), ]
  
  for(id in unique(temp_dat$ID)){
    ind = temp_dat$ID==id
    if(all(!temp_dat$fever_binary[ind])){ # never fever
      temp_dat$clearance_time[ind]=0
    } else if(all(temp_dat$fever_binary[ind])){ # always fever
      writeLines(sprintf('all fever for %s with %s FUP points',id,sum(ind)))
      temp_dat$clearance_time[ind] = max(temp_dat$Time[ind])
      temp_dat$clearance_time_cens[ind] = 0 #censored obs
    } else { # fever cleared
      j_cleared = which(ind & !temp_dat$fever_binary)
      check_ahead=F
      for(j in j_cleared){
        if(!check_ahead){
          ind_check = 
            which(ind & 
                    temp_dat$Time>temp_dat$Time[j] &
                    temp_dat$Time<temp_dat$Time[j] + window_clear)
          if(length(ind_check)>0 & all(!temp_dat$fever_binary[ind_check])){
            temp_dat$clearance_time[ind]=temp_dat$Time[j]
            check_ahead=T
          }
        }
      }
      if(!check_ahead){
        temp_dat$clearance_time[ind]=tail(temp_dat$Time[ind],1)
        temp_dat$clearance_time_cens[ind]=0
      }
    }
  }
  
  return(temp_dat[!duplicated(temp_dat$ID), ])
}



make_slopes_plot = function(stan_out, 
                            analysis_data_stan,
                            ID_map, 
                            data_summary,
                            my_lims = c(5, 72), # hours
                            my_vals = c(7,24,48,72)){
  
  slopes = rstan::extract(stan_out, pars='slope')$slope
  
  t12_output = data.frame(t_12_med = 24*log10(.5)/apply(slopes,2,mean),
                          t_12_up = 24*log10(.5)/apply(slopes,2,quantile,.9),
                          t_12_low = 24*log10(.5)/apply(slopes,2,quantile,.1),
                          slope_median = apply(slopes,2,median),
                          ID_stan = analysis_data_stan$id[analysis_data_stan$ind_start])
  t12_output = merge(t12_output, ID_map, by = 'ID_stan')
  data_summary = merge(data_summary, t12_output, by.x = 'ID', by.y = 'ID_key')
  
  data_summary = dplyr::arrange(data_summary, Treatment, t_12_med)
  
  par(mar=c(5,2,2,2))
  plot(data_summary$t_12_med, 1:nrow(data_summary),
       yaxt='n', xaxt='n',
       xlim=my_lims, panel.first=grid(), xlab='', 
       pch=15, ylab='', col=data_summary$trt_color)
  mtext(text = expression('t'[1/2] ~ ' (hours)'),side = 1,line=3)
  axis(1, at = my_vals,labels = my_vals)
  
  for(i in 1:nrow(data_summary)){
    lines(c(data_summary$t_12_low[i],
            data_summary$t_12_up[i]),
          rep(i,2),
          col=adjustcolor(data_summary$trt_color[i],alpha.f = .5))
  }
  
  for(kk in which(!duplicated(data_summary$Treatment))){
    ind = data_summary$Treatment==data_summary$Treatment[kk]
    writeLines(sprintf('In %s the median clearance half life was %s (range %s to %s)',
                       data_summary$Treatment[kk],
                       round(median(data_summary$t_12_med[ind]),1),
                       round(min(data_summary$t_12_med[ind]),1),
                       round(max(data_summary$t_12_med[ind]),1)))
    abline(v = median(data_summary$t_12_med[ind]), col=data_summary$trt_color[kk],lty=2,lwd=2)
  }
  
  legend('bottomright', legend = unique(data_summary$Treatment),
         col = unique(data_summary$trt_color),
         pch=15,lwd=2,inset=0.03)
  return(data_summary)
}

get_itt_population = function(){
  
  require(tidyverse)
  # TH58 and TH57 used envelopes so allocation is not recorded in the data-XXX.csv files
  rand.TH58 <- read.csv("~/Dropbox/PLATCOV/rand-TH58.csv")[0:9, ]
  rand.TH57 <- read.csv("~/Dropbox/PLATCOV/rand-TH57.csv")[0:10, ]
  rand.TH58$ID = paste('PLT-TH58-',rand.TH58$RandomisationID,sep='')
  rand.TH57$ID = paste('PLT-TH57-',rand.TH57$RandomisationID,sep='')
  
  ff_names = list.files(path = "~/Dropbox/PLATCOV", pattern = 'data',full.names = T)
  data_list = list()
  for(i in 1:length(ff_names)){
    data_list[[i]] = read.csv(ff_names[i])
    data_list[[i]]$Date = as.POSIXct(data_list[[i]]$Date,format='%a %b %d %H:%M:%S %Y')
    my_prefix=gsub(x = strsplit(ff_names[i], split = 'data-')[[1]][2], pattern = '.csv',replacement = '')
    print(my_prefix)
    
    data_list[[i]]$ID = paste('PLT-', my_prefix, '-', data_list[[i]]$randomizationID, sep='')
    data_list[[i]] = data_list[[i]][, c('ID', 'Treatment')]
  }
  data_list[[length(ff_names)+1]]=rand.TH57[, c('ID', 'Treatment')]
  data_list[[length(ff_names)+2]]=rand.TH58[, c('ID', 'Treatment')]
  
  xx = bind_rows(data_list)
  
  library(stringr)
  for(i in 1:nrow(xx)){
    id = unlist(strsplit(xx$ID[i],split = '-'))
    id[3] = str_pad(id[3], 3, pad = "0")
    id = paste(id, collapse = '-')
    xx$ID[i]=id
  }
  
  return(xx)
}




get_trt_colors = function(){
  trt_cols = array(dim = 5)
  names(trt_cols) = 
    c("Baloxavir",
      'No Study Drug',
      "Favipiravir",
      "Molnupiravir",
      "Oseltamivir")
  trt_cols['No Study Drug'] = viridis::viridis(n = 10)[8]
  trt_cols['Oseltamivir'] = viridis::magma(n = 10)[8]
  trt_cols['Favipiravir'] = viridis::plasma(n = 100)[92]
  trt_cols['Oseltamivir'] = viridis::plasma(n = 10)[1]
  trt_cols['Baloxavir'] = viridis::inferno(n = 10)[5]
  trt_cols['Molnupiravir'] = viridis::inferno(n = 10)[7]
  
  
  return(trt_cols)
}



plot_baseline_vl <- function(Baseline_data){
  
  Baseline_data$symptomDay <- as.factor(Baseline_data$symptomDay)
  
  G1 <- ggplot(Baseline_data, aes(x = fluType, y = Baseline.viral.load)) +
    geom_jitter(width = 0.1, size = 3, alpha = 0.75, aes(col = fluType)) +
    geom_boxplot(width = 0.5, alpha = 0.2, aes(fill = fluType), size = 0.8) +
    theme_bw() +
    scale_color_manual(values = c("#EE4266", "#387ADF"), name = "Influenza type")  +
    scale_fill_manual(values = c("#EE4266", "#387ADF"), guide = "none") +
    guides(shape = guide_legend(override.aes = list(size = 4))) +
    xlab("") +
    ylab("Baseline viral densities (log10 genomes/mL)") +
    theme(axis.text = element_text(size = 11),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          panel.spacing = unit(1, "lines")
    )
  
  G2 <- ggplot(Baseline_data, aes(x = symptomDay, y = Baseline.viral.load)) +
    geom_point(position=position_jitterdodge(dodge.width=0.8, jitter.width = 0.2), size = 3, alpha = 0.75, aes(col = fluType)) +
    geom_boxplot(width = 0.5, alpha = 0.2, aes(fill = fluType), size = 0.8, position=position_dodge(width=0.8)) +
    theme_bw() +
    scale_color_manual(values = c("#EE4266", "#387ADF"), name = "Influenza type")  +
    scale_fill_manual(values = c("#EE4266", "#387ADF"), guide = "none") +
    guides(shape = guide_legend(override.aes = list(size = 4))) +
    xlab("Time since symptom onset (days)") +
    ylab("Baseline viral densities (log10 genomes/mL)") +
    theme(axis.text = element_text(size = 11),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold"),
          panel.spacing = unit(1, "lines")
    )
  
  ggarrange(G1, G2, ncol = 2, labels = "AUTO", common.legend = T, legend = "right")
  
}


plot_vl_box <- function(dataplot, trt_colors, fluType = F, trt_order = NULL){
  dataplot <- dataplot %>%
    distinct(ID, Timepoint_ID, daily_VL, .keep_all = T) %>%
    filter(Timepoint_ID <= 5) 
  
  if(!is.null(trt_order)) {
    dataplot$Trt <- factor(dataplot$Trt, levels = trt_order)
  } else {
    dataplot$Trt <- factor(dataplot$Trt)
  }
  
  dataplot$fluType <- as.factor(dataplot$fluType)
  levels(dataplot$fluType) <- paste0("Influenza ", levels(dataplot$fluType))
  
  dataplot$Timepoint_ID <- as.factor(dataplot$Timepoint_ID)
  
  if(fluType){
    dataplot_median <- dataplot %>%
      group_by(Trt, Timepoint_ID, fluType) %>%
      summarise(median_VL = median(daily_VL), .groups = 'drop');
    f_tab <- dataplot %>%
      distinct(ID, Trt, fluType) %>%
      group_by(Trt, fluType) %>%
      summarise(n = n()) %>%
      as.data.frame()
  } else {
    dataplot_median <- dataplot %>%
      group_by(Trt, Timepoint_ID) %>%
      summarise(median_VL = median(daily_VL), .groups = 'drop');
    f_tab <- dataplot %>%
      distinct(ID, Trt) %>%
      group_by(Trt) %>%
      summarise(n = n()) %>%
      as.data.frame()
  }
  
  f_tab$lab <- paste0("n = ", f_tab$n)
  
  colors <- trt_colors[names(trt_colors) %in% unique(dataplot$Trt)]
  colors <- colors[levels(dataplot$Trt)]
  labels <- names(colors)
  
  # numeric x positions for boxplots/medians (1..5); continuous Time for points
  G <- ggplot(dataplot, aes(y = daily_VL, fill = Trt)) +
    # points placed at continuous Time
    geom_point(data = dataplot, aes(x = Time, y = daily_VL, shape = censor, fill = Trt),
               alpha = 0.25, size = 1.5, stroke = 0.5) +
    # boxplots centered at integer positions (as.numeric(Timepoint_ID)); group by Timepoint_ID
    geom_boxplot(data = dataplot, aes(x = as.numeric(Timepoint_ID)-1, y = daily_VL, group = Timepoint_ID, fill = Trt),
                 width = 0.65, size = 0.5, outlier.shape = NA, coef = 0, alpha = 0.3) +
    # median line/points at integer positions
    geom_line(data = dataplot_median, aes(x = as.numeric(Timepoint_ID)-1, y = median_VL, group = Trt, col = Trt),
              linewidth = 1.2, linetype = 1) +
    geom_point(data = dataplot_median, aes(x = as.numeric(Timepoint_ID)-1, y = median_VL, fill = Trt),
               size = 3.5, shape = 24, col = "black") +
    scale_shape_manual(values = c(25, 21), guide = NULL) +
    scale_color_manual(label = labels, values = colors, name = "") +
    scale_fill_manual(label = labels, values = colors, name = "") +
    theme_bw() +
    # set x axis to continuous with integer breaks labeled by Timepoint_ID levels
    scale_x_continuous(breaks = seq_along(levels(dataplot$Timepoint_ID))-1,
                       labels = levels(dataplot$Timepoint_ID)) +
    scale_y_continuous(labels = label_math(), breaks = seq(0,10,2), limits = c(0,9)) +
    xlab("Time since randomisation (days)") +
    ylab("Viral densities (genomes/mL)") + 
    theme(axis.title  = element_text(face = "bold"),
          plot.title = element_text(face = "bold"),
          legend.position = "none",
          axis.text = element_text(size = 10),
          strip.text = element_text(size = 10, face = 'bold')) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed", linewidth = 0.75) +
    geom_text(data = f_tab, x = 5, y = 9, aes(label = lab),
              hjust = 1, vjust = 1) +
    ggtitle("Viral density dynamics") 
  
  if(fluType){G <- G +  facet_grid(fluType~Trt)} else {G <- G +  facet_grid(.~Trt)}
  G
}


slope_to_hl  <- function(slope){
  24*log10(2)/(-(slope)) 
}

formatter <- function(x){  
  (x-1)*100 
}

plot_hl <- function(Half_life, trt_colors, trt_order){
  Half_life$fluType <- as.factor(Half_life$fluType)
  levels(Half_life$fluType) <- paste0("Influenza ", levels(Half_life$fluType))
  
  Half_life <- Half_life %>%
    mutate(Trt = factor(Trt, levels = rev(trt_order))) %>%
    arrange(Trt, t_12_med) %>%
    mutate(ID = factor(ID, levels = unique(ID)))
  
  
  
  Half_life_med <- Half_life %>%
    group_by(Trt) %>%
    summarise(med_hl = median(t_12_med))  %>%
    as.data.frame()
  
  colors <- trt_colors[names(trt_colors) %in% unique(Half_life$Trt)]
  colors <- colors[levels(Half_life$Trt)]
  labels <- names(colors)
  
  f_tab <- Half_life %>%
    distinct(ID, Trt) %>%
    group_by(Trt) %>%
    summarise(n = n()) %>%
    as.data.frame()
  f_tab$med_hl <- Half_life_med$med_hl
  
  f_tab$lab <- paste0(f_tab$Trt, " (n = ", f_tab$n, "): ", sprintf("%.1f", f_tab$med_hl), " h")
  f_tab <- f_tab %>% arrange(match(Trt, (trt_order)))
  freq_lab <- paste(f_tab$lab, collapse = '\n')
  
  
  G <- ggplot(Half_life, aes(x = t_12_med, y = as.factor(as.numeric(ID)), col = Trt)) +
    geom_errorbar(aes(xmin = t_12_low, xmax = t_12_up),width = 0, alpha = 0.4) +
    geom_point(size = 2.5, aes(shape = fluType), alpha = 0.8) +
    geom_vline(data = Half_life_med, aes(xintercept = med_hl, col = Trt),linewidth = 1) +
    theme_bw(base_size = 18) +  
    scale_y_discrete(expand = c(0.01,0.01), breaks = NULL) +
    scale_color_manual(label = labels, values = colors, name = "") +
    scale_shape_manual(values = c(16, 17), name = "") +
    scale_x_continuous(breaks = seq(0,40,5), expand = c(0,0), minor_breaks = NULL) +
    guides(color = guide_legend(override.aes=list(linetype = rep(0, length(unique(Half_life$Trt)))))) +
    theme(axis.ticks.y = element_blank(),
          axis.title = element_text(face = "bold"),
          plot.title = element_text( face = "bold"),
          legend.position = "bottom",
          legend.box="vertical",
          legend.margin=margin(),
          legend.title = element_text(face = "bold")) +
    
    # Dynamically update xlim on the right and the legend size
    coord_cartesian(xlim = c(0, max(Half_life$t_12_med) + 5)) +
    xlab("Estimated viral clearance half-life (h)") +
    ylab("Individual participants") +
    ggtitle("A) Individual viral clearance half-life\n") +
    annotate("text", x = 27, y = nrow(Half_life)/8, label = freq_lab, hjust = 0, vjust = 1, size = 4) 
  G
  
  
}


plot_hl_flutype <- function(Half_life, trt_colors, trt_order){
  Half_life <- Half_life |>
    mutate(
      fluType = as.factor(fluType),
      fluType = factor(paste0("Influenza ", levels(fluType))[as.integer(fluType)], levels = paste0("Influenza ", unique(fluType)))
    )
  
  Half_life <- Half_life |>
    mutate(Trt = factor(Trt, levels = rev(trt_order))) |>
    arrange(fluType, Trt, t_12_med) |>
    mutate(ID = factor(ID, levels = unique(ID)))
  
  Half_life_med <- Half_life |>
    group_by(Trt, fluType) |>
    summarise(med_hl = median(t_12_med, na.rm = TRUE), .groups = "drop") |>
    as.data.frame()
  
  colors <- trt_colors[names(trt_colors) %in% unique(Half_life$Trt)]
  colors <- colors[levels(Half_life$Trt)]
  labels <- names(colors)
  
  f_tab <- Half_life |>
    distinct(ID, Trt, fluType) |>
    group_by(Trt, fluType) |>
    summarise(n = n(), .groups = "drop") |>
    as.data.frame()
  
  f_tab <- f_tab |>
    left_join(Half_life_med, by = c("Trt", "fluType"))
  
  f_tab$lab <- paste0(f_tab$Trt, " (n = ", f_tab$n, "): ", sprintf("%.1f", f_tab$med_hl), " h")
  f_tab <- f_tab %>% arrange(match(Trt, (trt_order)))
  
  freq_lab <- f_tab |>
    group_by(fluType) |>
    summarise(freq_lab = paste(lab, collapse = "\n"), .groups = "drop")
  
  # X position for per-facet labels (relative)
  x_max <- max(Half_life$t_12_med, na.rm = TRUE)
  x_pos <- 33.5
  freq_lab <- freq_lab |> mutate(x_pos = x_pos)
  
  G <- ggplot(Half_life, aes(x = t_12_med, y = ID)) +
    geom_errorbar(aes(xmin = t_12_low, xmax = t_12_up, colour = Trt), width = 0, alpha = 0.4) +
    geom_point(size = 2.5, aes(shape = fluType, colour = Trt)) +
    geom_vline(data = Half_life_med, aes(xintercept = med_hl, colour = Trt), linewidth = 1) +
    theme_bw(base_size = 14) +
    scale_y_discrete(expand = c(0.01, 0.01), breaks = NULL) +
    scale_color_manual(label = labels, values = colors, name = "") +
    scale_shape_manual(values = c(16, 17), name = "") +
    scale_x_continuous(breaks = seq(0, 40, 5), expand = c(0, 0)) +
    guides(color = guide_legend(override.aes = list(linetype = rep(0, length(unique(Half_life$Trt)))))) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text( face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(),
      legend.title = element_text(face = "bold")
    ) +
    coord_cartesian(xlim = c(0, x_max + 5)) +
    xlab("Estimated viral clearance half-life (h)") +
    ylab("") +
    ggtitle("Individual viral clearance half-life\n") +
    geom_text(
      data = freq_lab,
      mapping = aes(x = x_pos, label = freq_lab),
      inherit.aes = FALSE,
      y = Inf,
      hjust = 0,
      vjust = 1.2,
      size = 3.5
    ) +   facet_grid(fluType ~ ., scales = "free_y")
  
  G
}

plot_trt_effs <- function(effect_ests, model_cols, trt_order){
  effect_ests_plots <- NULL
  trt_order2 <- trt_order[trt_order != ref_arm] %>% rev()
  
  for(i in 1:length(effect_ests)){
    effect_ests_plot <- as.data.frame(effect_ests[[i]])
    effect_ests_plot <- exp(effect_ests_plot)
    colnames(effect_ests_plot)[1:5] <- c("L95", "L80", "med", "U80", "U95")
    effect_ests_plot$arm <- row.names(effect_ests_plot)
    effect_ests_plot$arm <- as.factor(effect_ests_plot$arm)
    effect_ests_plot$arm <- factor(effect_ests_plot$arm, levels = trt_order2)
    
    effect_ests_plot$model <- names(effect_ests)[i]
    effect_ests_plots <- rbind(effect_ests_plots, effect_ests_plot)
  }
  
  effect_ests_plots$model <- as.factor(effect_ests_plots$model)
  
  #Labeling reference arm
  lab_ref <- ref_arm
  #Labeling intervention arm
  my.labs <- levels(effect_ests_plot$arm)
  
  title <- paste0("B) Estimated treatment effects relative to \n", tolower(lab_ref), " arm")
  
  names(model_cols) <- levels(effect_ests_plots$model)
  
  G <- ggplot(effect_ests_plots, 
              aes(x = arm, y = med, col = model, group = model)) +
    geom_rect(aes(ymin = min(0.75, min(L95)-0.05), ymax = study_threshold, xmin = 0, xmax = length(my.labs)+1), fill = "#7D7C7C", alpha = 0.2, col = NA) +
    geom_point(position = position_dodge(width = 0.5), size = 4.5) +
    geom_errorbar(aes(x = arm, ymin = L95, ymax = U95),position = position_dodge(width = 0.5), width = 0, linewidth = 1) +
    geom_errorbar(aes(x = arm, ymin = L80, ymax = U80),position = position_dodge(width = 0.5), width = 0, linewidth = 2.25) +
    scale_color_manual(values = model_cols, guide = NULL) +
    coord_flip() +
    theme_bw(base_size = 18) +
    geom_hline(yintercept = 1, col = "red", linetype = "32") +
    scale_y_continuous(
      labels = formatter,
      limits = c(min(0.75, min(effect_ests_plots$L95) - 0.05), max(effect_ests_plots$U95) + .25),
      expand = c(0, 0),
      breaks = (function(rng){
        rng_pct <- rng * 100
        seq(from = floor(rng_pct[1] / 20) * 20, to = ceiling(rng_pct[2] / 20) * 20, by = 20) / 100
      })(c(min(0.75, min(effect_ests_plots$L95) - 0.05), max(effect_ests_plots$U95) + .25))
    ) +
    scale_x_discrete(labels= my.labs) +
    ylab("Change in viral clearance rate (%)") +
    xlab("") +
    ggtitle(title)  + 
    theme( axis.title  = element_text(face = "bold"),
           plot.title  = element_text(face = "bold"),
           legend.position = "bottom",
           legend.title = element_text(face = "bold"))
  G
}



plot_trt_effs_2 <- function(effect_ests, model_cols, lower, upper){
  
  effect_ests_plots <- NULL
  for(i in 1:length(effect_ests)){
    effect_ests_plot <- as.data.frame(effect_ests[[i]])
    effect_ests_plot <- (exp(effect_ests_plot) - 1) * 100
    
    ## Here: already in percent change units
    colnames(effect_ests_plot)[1:5] <- c("L95", "L80", "med", "U80", "U95")
    effect_ests_plot$arm <- row.names(effect_ests_plot)
    effect_ests_plot$arm <- as.factor(effect_ests_plot$arm)
    effect_ests_plot$model <- names(effect_ests)[i]
    effect_ests_plots <- rbind(effect_ests_plots, effect_ests_plot)
  }
  
  effect_ests_plots$model <- as.factor(effect_ests_plots$model)
  
  lab_ref <- ref_arm
  my.labs <- levels(effect_ests_plot$arm)
  title <- paste0("Estimated treatment effects relative to ", tolower(lab_ref), " arm")
  names(model_cols) <- levels(effect_ests_plots$model)
  
  ## Calculate y-limits
  ymin_all <- min(lower, min(effect_ests_plots$L95) - 5)
  ymax_all <- max(upper, max(effect_ests_plots$U95) + 5)
  
  G <- ggplot(effect_ests_plots, aes(x = arm, y = med, col = model, group = model)) +
    
    ## Shaded regions:
    geom_rect(aes(ymin = ymin_all, ymax = -10, xmin = 0, xmax = length(my.labs)+1),
              fill = "#7D7C7C", alpha = 0.25, col = NA, inherit.aes = FALSE) +
    geom_rect(aes(ymin = -10, ymax = 20, xmin = 0, xmax = length(my.labs)+1),
              fill = "#E6E6E6", alpha = 0.25, col = NA, inherit.aes = FALSE) +
    # geom_rect(aes(ymin = ymin_all, ymax = -10, xmin = 0, xmax = length(my.labs)+1),
    #           fill = "#EE4266", alpha = 0.5, col = NA, inherit.aes = FALSE) +
    # geom_rect(aes(ymin = -10, ymax = 20, xmin = 0, xmax = length(my.labs)+1),
    #           fill = "#D6FFD6", alpha = 0.2, col = NA, inherit.aes = FALSE) +
    
    ## Points and error bars:
    geom_point(position = position_dodge(width = 0.5), size = 4) +
    geom_errorbar(aes(ymin = L95, ymax = U95), position = position_dodge(width = 0.5),
                  width = 0, linewidth = 0.65) +
    geom_errorbar(aes(ymin = L80, ymax = U80), position = position_dodge(width = 0.5),
                  width = 0, linewidth = 1.5) +
    
    scale_color_manual(values = model_cols) +
    coord_flip() +
    theme_bw() +
    
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    scale_y_continuous(limits = c(ymin_all, ymax_all),
                       breaks = seq(lower, upper, 10),
                       expand = c(0,0)) +
    scale_x_discrete(labels = my.labs) +
    ylab("Change in viral clearance rate (%)") +
    xlab("") +
    ggtitle(title)  + 
    theme(axis.title  = element_text(face = "bold"),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom",
          axis.text = element_text(size = 10))
  
  G
}


plot_randomisation <- function(Baseline_data, trt_order){
  Baseline_data <- Baseline_data %>%
    mutate(Trt = factor(Trt, levels = rev(trt_order)))
  
  Baseline_data$Rand_date <- as.Date(Baseline_data$Rand_date)
  
  ggplot(Baseline_data, aes(x = Rand_date, y = Trt, col = fluType)) +
    geom_jitter(size = 3, alpha = 0.5, height = 0.1) +
    theme_bw() +
    xlab("Randomisation date") +
    ylab("") +
    scale_x_date(date_breaks = "6 months", date_minor_breaks = "3 months",
                 date_labels = "%b %y") +
    scale_color_manual(values = c("#EE4266", "#387ADF"), name = "Influenza type")  +
    theme(axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 11))
  
  
}


plot_inds <- function(model_selects, Baseline_data){
  model_list = list()
  
  stan_inputs_i <- stan_inputs[[model_settings$trt_formulas[model_selects[1]]]]
  adastra_dat_analysis <- data.frame("ID" = stan_inputs_i$analysis_data_stan$id,
                                     "log10_viral_load" = stan_inputs_i$analysis_data_stan$log_10_vl,
                                     "Time" = stan_inputs_i$analysis_data_stan$obs_day)
  ID_map <- stan_inputs_i$ID_map
  for(i in model_selects){
    load(ff[i])
    model_list[[i]] = out
  }
  preds <- lapply(model_list, rstan::extract, "preds")
  preds <- lapply(preds, `[[`, 1)
  preds_list <- lapply(preds, function(x) as.data.frame(t(apply(x, 2, quantile, c(0.025, 0.5, 0.975)))))
  preds_list <- lapply(preds_list, function(x) cbind(adastra_dat_analysis, x))
  models <- c("Linear")#, "Non-linear")
  preds_list <- Map(cbind, preds_list, models)
  col_names <- c(colnames(adastra_dat_analysis), c("low", "med", "up", "model"))
  preds_list <- lapply(preds_list, setNames, col_names)
  preds_list <- lapply(preds_list, function(x) data.frame(x, "censor" = c(rep("none", stan_inputs_i$analysis_data_stan$N_obs), 
                                                                          rep("left", (stan_inputs_i$analysis_data_stan$Ntot - stan_inputs_i$analysis_data_stan$N_obs)))))
  preds_list <- lapply(preds_list, function(x) merge(x, ID_map, by.x = "ID", by.y = "ID_stan"))
  preds_list <- lapply(preds_list, function(x) merge(x, as.data.frame(Baseline_data[,c("ID", "Trt", "fluType")]), by.x = "ID_key", by.y = "ID"))
  preds_dat <- do.call("rbind", preds_list)
  preds_dat$Trt <- as.character(preds_dat$Trt)
  preds_dat$Trt <- as.factor(preds_dat$Trt)
  preds_dat$model <- as.factor(preds_dat$model)
  preds_dat$censor <- as.factor(preds_dat$censor)
  ID_map = merge(stan_inputs_i$ID_map, Baseline_data, by.x = 'ID_key',by.y = 'ID')
  ID_map <- ID_map[order(ID_map$ID_key),]
  
  ind_plot_list <- list()
  #resid_dat <- NULL
  for(i in 1:nrow(ID_map)){
    plot_data <- preds_dat %>% filter(ID_key == ID_map$ID_key[i])
    plot_data$resid <- plot_data$log10_viral_load - plot_data$med
    plot_data$Timepoint_ID <- round(plot_data$Time)
    #resid_dat <- rbind(resid_dat, plot_data)
    
    lab <- paste0(plot_data$ID_key[1], "\n", plot_data$Trt[1], ": Influenza ", plot_data$fluType)
    
    ind_plot_list[[i]] <- ggplot() +
      geom_point(data = plot_data[plot_data$model == "Linear",], aes(x = Time, y = log10_viral_load, shape = censor),
                 size = 1.75, alpha = 0.8) +
      # geom_ribbon(data = plot_data, aes(x = Time, ymin = low, ymax = up, fill = model), alpha = 0.2) +
      geom_line(data = plot_data, aes(x = Time, y = med, col = model), linewidth = 1, alpha = 0.8) +
      theme_bw() +
      scale_y_continuous(labels=label_math(), breaks = seq(0,10,2)) +
      coord_cartesian(ylim = c(0,9), xlim = c(0,6))+
      scale_x_continuous(breaks = 0:14) +
      ylab("") +
      xlab("") +
      theme(
        axis.title  = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
        legend.position = "none",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), 'lines')) +
      scale_color_manual(values = c("#1640D6", "#BE3144"), name = "Model") +
      scale_fill_manual(values = c("#1640D6", "#BE3144"), name = "Model") +
      scale_shape_manual(values = c(6, 1), guide = "none", drop=FALSE) +
      geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
      ggtitle(lab)
    
  }
  return(ind_plot_list)
}

plot_vl_base <- function(dataplot, fluType = F){
  dataplot <- dataplot %>%
    distinct(ID, Timepoint_ID, daily_VL, .keep_all = T) %>%
    filter(Timepoint_ID <= 5) 
  
  dataplot$fluType <- as.factor(dataplot$fluType)
  levels(dataplot$fluType) <- paste0("Influenza ", levels(dataplot$fluType))
  
  dataplot$Timepoint_ID <- as.factor(dataplot$Timepoint_ID)
  
  if(fluType){
    dataplot_median <- dataplot %>%
      group_by(Timepoint_ID, fluType) %>%
      summarise(median_VL = median(daily_VL), .groups = 'drop');
    f_tab <- dataplot %>%
      distinct(ID, fluType) %>%
      group_by(fluType) %>%
      summarise(n = n()) %>%
      as.data.frame()
  } else {
    dataplot_median <- dataplot %>%
      group_by( Timepoint_ID) %>%
      summarise(median_VL = median(daily_VL), .groups = 'drop');
    f_tab <- dataplot %>%
      distinct(ID) %>%
      summarise(n = n()) %>%
      as.data.frame()
  }
  
  f_tab$lab <- paste0("n = ", f_tab$n)
  
  
  G <- ggplot(dataplot, aes(x = Timepoint_ID, y = daily_VL)) +
    geom_point(alpha = 0.25, size = 2, aes(shape = censor),
               stroke = 0.5) +
    geom_line(aes(group = ID), alpha = 0.25, linewidth = 0.4) +
    geom_line(data = dataplot_median, aes(x =  Timepoint_ID, y = median_VL, group = 1, col = fluType), 
              linewidth = 1.2, linetype = 1) +
    geom_point(data = dataplot_median, aes(x = Timepoint_ID, y = median_VL, fill = fluType), 
               size = 3.5, shape = 24, col = "black", stroke = 0.75) +
    scale_shape_manual(values = c(25, 21), guide = NULL) +
    scale_color_manual(label = labels, values = c("#EE4266", "#387ADF"), name = "") +
    scale_fill_manual(label = labels, values = c("#EE4266", "#387ADF"), name = "") +
    theme_bw() +
    scale_y_continuous(labels=label_math(), breaks = seq(0,10,2), limits = c(0,9)) +
    xlab("Time since randomisation (days)") +
    ylab("Viral densitites (log10 genomes/mL)") + 
    theme(axis.title  = element_text(face = "bold"),
          plot.title = element_text(face = "bold"),
          legend.position = "none",
          axis.text = element_text(size = 10),
          strip.text = element_text(size = 10, face = 'bold')) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed", linewidth = 0.75) +
    geom_text(data = f_tab, x = 6, y = 9, aes(label = lab),
              hjust = 1, vjust = 1) +
    ggtitle("Viral density dynamics")
  
  if(fluType){G <- G +  facet_grid(.~fluType)}
  G
  
}

plot_mITT <- function(adastra_dat, start_dt = "2023-04-01"){
  G <- adastra_dat %>%
    mutate(symptomDay = ceiling(symptomDay)) %>% # round up to nearest integer
    distinct(ID,.keep_all = T) %>%
    group_by(symptomDay, fluType) %>%
    mutate(n = n(),
           n_mITT_F = sum(!mITT),
           lab = paste0(n_mITT_F, "/",n, " (", round(100*n_mITT_F/n), "%)"),
           fluType_lab = paste0("Influenza ", fluType),
           fluType_lab = as.factor(fluType_lab),
           symptomDay_lab = ifelse(symptomDay == 1, paste0(symptomDay, " day"), paste0(symptomDay, " days")),
           symptomDay_lab = as.factor(symptomDay_lab)) %>%
    ungroup()
  
  # Dynamically update to number of treatments
  max_y <- length(unique(G$Trt))
  
  G <- G %>%
    ggplot(aes(x = Rand_date, y = Trt)) +
    geom_point(aes(col = mITT, size = mITT, shape = all_negative), alpha = 0.5) +
    theme_bw(base_size = 13) +
    scale_size_manual(values = c(5,2)) +
    scale_shape_manual(values = c(16,17), name = "All negative qPCR") +
    facet_grid(symptomDay_lab~fluType_lab) +
    geom_text(aes(x = as.POSIXct(start_dt), y = max_y+0.3, label = lab, hjust = 0)) +
    ylab("") +
    xlab("Randomisation date")
  
  G
}

cal_fever_clearance <- function(temp_dat_for_plot, threshold, window_clear){
  
  temp_dat_for_plot = temp_dat_for_plot %>%
    filter(Fever_Baseline==1) %>%
    mutate(fever_binary = temp > threshold,
           clearance_time = NA,
           clearance_time_cens = 1,# For interval censored data, the status indicator is 0=right censored, 1=event at time, 2=left censored, 3=interval censored. 
           Time = ifelse(Time<0, 0, Time),
    ) %>% group_by(Label) %>%
    mutate(n=n()) %>% filter(n>1)
  
  ############################################################################################
  for(id in unique(temp_dat_for_plot$Label)){
    ind = temp_dat_for_plot$Label==id
    
    if(all(temp_dat_for_plot$fever_binary[ind])){ # always fever
      
      writeLines(sprintf('all fever for %s with %s FUP points',id,sum(ind)))
      temp_dat_for_plot$clearance_time[ind] = max(temp_dat_for_plot$Time[ind])
      temp_dat_for_plot$clearance_time_cens[ind] = 0 #censored obs
      
    } else { # fever cleared
      j_cleared = which(ind & !temp_dat_for_plot$fever_binary & temp_dat_for_plot$Time>0)
      check_ahead=F
      for(j in j_cleared){
        if(!check_ahead){
          ind_check = 
            which(ind & 
                    temp_dat_for_plot$Time>temp_dat_for_plot$Time[j] &
                    temp_dat_for_plot$Time<temp_dat_for_plot$Time[j] + window_clear)
          if(length(ind_check)>0 & all(!temp_dat_for_plot$fever_binary[ind_check])){
            temp_dat_for_plot$clearance_time[ind]=temp_dat_for_plot$Time[j]
            check_ahead=T
          }
        }
      }
      if(!check_ahead){
        temp_dat_for_plot$clearance_time[ind]=tail(temp_dat_for_plot$Time[ind],1)
        temp_dat_for_plot$clearance_time_cens[ind]=0
      }
    }
  }
  return(temp_dat_for_plot)
}


plot_inds_obs <- function(adastra_dat, adastra_dat_analysis){
  ID_map <- adastra_dat %>% distinct(ID, .keep_all = T) %>%
    arrange(Trt, fluType)
  
  ind_plot_list <- list()
  for(i in 1:nrow(ID_map)){
    plot_data <- adastra_dat %>% filter(ID == ID_map$ID[i], Timepoint_ID<=5)
    lab <- paste0(plot_data$ID[1], "\n", plot_data$Trt[1], ": Influenza ", plot_data$fluType[1])
    
    if(!plot_data$ID[1] %in% adastra_dat_analysis$ID){col = "red"}else{col = "black"}
    
    ind_plot_list[[i]] <- ggplot() +
      geom_point(data = plot_data, aes(x = Time, y = log10_viral_load, shape = censor),
                 size = 1.75, alpha = 0.8, fill = col) +
      theme_bw() +
      scale_y_continuous(labels=label_math(), breaks = seq(0,10,2)) +
      coord_cartesian(ylim = c(0,9), xlim = c(0,6))+
      scale_x_continuous(breaks = 0:14) +
      ylab("") +
      xlab("") +
      theme(
        axis.title  = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
        legend.position = "none",
        plot.margin = unit(c(0.1,0.1,0.1,0.1), 'lines')) +
      scale_shape_manual(values = (c(25, 21)), guide = "none", drop=FALSE) +
      geom_hline(yintercept = 0, col = "black", linetype = "dashed") +
      geom_hline(yintercept = log10(250), col = "red", linetype = "dashed") +
      ggtitle(lab) 
    
  }
  return(ind_plot_list)
}

checkStrict(make_stan_inputs)
checkStrict(plot_serial_data)
checkStrict(plot_effect_estimates)
checkStrict(plot_data_model_fits)
checkStrict(plot_coef_effects)
checkStrict(calculate_fever_clearance)
