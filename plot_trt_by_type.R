effect_ests <- list()

for(i in 1:length(ff)){
  load(ff[i])
  
  #All variants
  if(i == 1){
    effect_ests[[i]] = 
      as.data.frame(summary(out, pars='trt_effect',use_cache=F,probs=my_probs)$summary[,c('2.5%','10%','50%','90%','97.5%'),drop=F])
    arms <- levels(adastra_dat_analysis$Trt)
    rownames(effect_ests[[i]]) <-  arms[-1]
    rownames(effect_ests[[i]]) <- paste0("All types:", row.names(effect_ests[[i]]))
  } 
  
  #Interactions with variants
  if(i == 2){
    posterior <- rstan::extract(out, pars = 'trt_effect')$trt_effect %>% as.data.frame()
    colnames(posterior) <- colnames(stan_inputs[[model_settings$trt_formulas[i]]]$Treatment_matrix)[-1]

    print(colnames(posterior))
    
    vars <-   levels(as.factor(adastra_dat_analysis$fluType))
    ref_var <- vars[1]
    
    arms <- levels(adastra_dat_analysis$Trt)
    ref_arm <- arms[1]
    
    effect <- array(NA, c(nrow(posterior), length(vars), length(arms)-1))
    for(j in 1:length(vars)){
      for(k in 1:(length(arms)-1)){
        colnames(posterior)
        inds <- (colnames(posterior) == paste0("Trt", arms[k+1])) |
          (grepl( paste0("Trt", arms[k+1]), colnames(posterior)) & grepl(paste0("fluType", vars[j]), colnames(posterior)))
        
        if(sum(inds) == 1){effect[,j,k] <- posterior[,which(inds)]} else {effect[,j,k] <- rowSums(posterior[,which(inds)])}
        
      }
     
    }
    
    effect <- as.data.frame(effect)
    colnames(effect) <- apply(expand.grid(paste0("fluType",vars), arms[-1]), 1, paste, collapse=":")
    effect_ests[[i]] = t(apply(effect, 2, quantile, my_probs))
  }
}

effect_ests <- do.call("rbind", effect_ests)

effect_ests$gr <- row.names(effect_ests)
effect_ests <- effect_ests %>% 
  separate(gr, c("fluType", "Trt"), sep = ":") 
effect_ests$fluType[effect_ests$fluType == "fluTypeA"] <- "Influenza A"
effect_ests$fluType[effect_ests$fluType == "fluTypeB"] <- "Influenza B"

  effect_ests_plots <- as.data.frame(effect_ests)
  effect_ests_plots <- cbind(exp(effect_ests_plots[1:5]),effect_ests_plots[6:7])
  colnames(effect_ests_plots)[1:5] <- c("L95", "L80", "med", "U80", "U95")
  
effect_ests_plots$Trt <- factor(effect_ests_plots$Trt, levels = rev(trt_order))
effect_ests_plots$fluType <- factor(effect_ests_plots$fluType, levels = rev(c("All types", "Influenza B", "Influenza A")))

  lab_ref <- ref_arm

  title <- paste0("Estimated treatment effects relative to \n", tolower(lab_ref), " arm")
  
A <-  ggplot(effect_ests_plots,  aes(x = Trt, y = med, col = fluType, group = fluType)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax =  1.2, 
             fill = "#7D7C7C", alpha = 0.5, col = NA) +
    geom_point(position = position_dodge(width = 0.5), size = 4) +
    geom_errorbar(aes(x = Trt, ymin = L95, ymax = U95),position = position_dodge(width = 0.5), width = 0, linewidth = 0.65) +
    geom_errorbar(aes(x = Trt, ymin = L80, ymax = U80),position = position_dodge(width = 0.5), width = 0, linewidth = 1.5) +
    scale_color_manual(values = rev(c("#002379", "#387ADF", "#EE4266")), name = "Influenza type") +
    theme_bw(base_size = 18) +
    geom_hline(yintercept = 1, col = "red", linetype = "dashed") +
    scale_y_continuous(labels = formatter, expand = c(0,0), 
                       breaks = seq(0.2,6, 0.4)) +
    coord_flip(ylim = c(min(0.75, min(effect_ests_plots$L95)-0.05), min(c(4, max(effect_ests_plots$U95) + .25)))) +
    ylab("Change in viral clearance rate (%)") +
    xlab("") +
    ggtitle(title)  + 
    theme(axis.title  = element_text(face = "bold"),
          plot.title = element_blank(),
          legend.position = "bottom"
    )


A 

plot_name <- here("Plots", intervention, paste0(intervention,  "_mITT_", mITT_threshold, "_esimated_effect_vs_by_type.png"))

png(plot_name, width = 8, height = 6, unit = 'in', res = 350)
A
dev.off()


effect_ests_plots %>%
  mutate(lab = paste0(
    round((med-1)*100),
    "% [95%CrI: ",
    round((L95-1)*100),
    " to ",
    round((U95-1)*100),
    "%]"))
