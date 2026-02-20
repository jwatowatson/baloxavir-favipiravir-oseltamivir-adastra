intervention = 'Final'
ref_arm = 'No Study Drug'

ref_code <- ifelse(ref_arm == "No Study Drug", "NSD", ref_arm)

ff = list.files(
  'Rout/',
  pattern = paste0(
    intervention, "(_mITT)?(_[0-9]+)?_vs_", ref_code,
    "\\.RData$"
  )
)
ff = ff[grep(pattern = 'model_fits_1_', x = ff)]
# if(!length(ff)==nrow(model_settings)) stop('not all outputs are ready for all model settings')
ff = paste0('Rout/', ff)



#main_mod = 1
model_cols = c("#1640D6", "#16B2D6")
names(model_cols) = paste('model', 1:2)#nrow(model_settings))

effect_ests=list()
for(i in 1:length(ff)){
  load(ff[i])
  effect_ests[[i]] = 
    summary(out, pars='trt_effect',use_cache=F,probs=my_probs)$summary[,c('2.5%','10%','50%','90%','97.5%'),drop=F]
  rownames(effect_ests[[i]]) = trts
  names(effect_ests)[i] <- names(model_cols)[i]
}




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

effect_ests_plots$model <- factor(effect_ests_plots$model, levels = rev(names(model_cols)),
                                  labels = rev(c("mITT population", "ITT population")))

## Build label: "XX% [XX to XX]" using the same formatter as the y-axis
effect_ests_plots$lbl <- paste0(
  round(formatter(effect_ests_plots$med)), "% [",
  round(formatter(effect_ests_plots$L95)), " to ",
  round(formatter(effect_ests_plots$U95)), "]"
)

## Put mITT label above its CI, ITT label below its CI
effect_ests_plots$lbl_vjust <- ifelse(effect_ests_plots$model == "mITT population", -0.9, 1.8)

#Labeling reference arm
lab_ref <- ref_arm
#Labeling intervention arm
my.labs <- levels(effect_ests_plot$arm)

title <- paste0("B) Estimated treatment effects relative to \n", tolower(lab_ref), " arm")

names(model_cols) <- rev(levels(effect_ests_plots$model))

G <- ggplot(effect_ests_plots, 
            aes(x = arm, y = med, col = model, group = model)) +
  geom_rect(aes(ymin = min(0.75, min(L95)-0.05), ymax = study_threshold, xmin = 0, xmax = length(my.labs)+1), fill = "#7D7C7C", alpha = 0.2, col = NA) +
  geom_point(position = position_dodge(width = 0.5), size = 6) +
  geom_errorbar(aes(x = arm, ymin = L95, ymax = U95),position = position_dodge(width = 0.5), width = 0, linewidth = 1.25) +
  geom_errorbar(aes(x = arm, ymin = L80, ymax = U80),position = position_dodge(width = 0.5), width = 0, linewidth = 3) +
  geom_text(aes(label = lbl, vjust = lbl_vjust), colour = "black",
            position = position_dodge(width = 0.5),
            size = 4.2, show.legend = FALSE) +
  scale_color_manual(values = model_cols, name = NULL) +
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
  #ggtitle(title)  + 
  theme( axis.title  = element_text(face = "bold"),
         plot.title  = element_text(face = "bold"),
         legend.position = "bottom",
         legend.title = element_text(face = "bold"))
G

plot_name <- here("Plots", intervention, paste0( "sensitivity_trt_effect.png"))
png(plot_name, width = 8, height = 6, units = "in", res = 350)
G
dev.off()
