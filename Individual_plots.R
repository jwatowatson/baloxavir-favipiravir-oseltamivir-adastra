library(ggplot2)
library(ggforce)

df = readr::read_csv('Analysis_Data/interim_firstResults_analysis.csv')

mITT_threshold <- 250
use_threshold = T

df = df %>% group_by(ID) %>%
  mutate(
    baseline_vl = mean(log10_viral_load[which(Timepoint_ID==0)]),
    mITT = ifelse(use_threshold, baseline_vl > log10(mITT_threshold), T),
    all_negative = ifelse(all(CT == 40), T, F),
    ID_Trt = paste(ID, Trt, sep = " - "),
    ID = factor(ID)) %>%
  group_by(ID, Timepoint_ID) %>%
    mutate(CTmean = mean(CT)) %>% filter(Timepoint_ID<6) %>% 
  arrange(Site, Trt, ID)

# ensure ID is a factor (important for ordering)

# how many panels per page
nrow <- 4
ncol <- 4
per_page <- nrow * ncol
ct_min = min(df$CT)

ids <- levels(df$ID_Trt)
id_pages <- split(ids, ceiling(seq_along(ids) / per_page))

# open pdf device
pdf("individual_curves.pdf", width = 11, height = 8.5)

for (pg in seq_along(id_pages)) {
  ids_this <- id_pages[[pg]]
  
  dpg <- df[df$ID_Trt %in% ids_this, , drop = FALSE]

    p <- dpg %>% ggplot( aes(x = Time, y = CTmean, colour = mITT)) +
    geom_line() + ylim(40, 15)+
    geom_point(aes(x=Time, y=CT))+
      facet_wrap(~ID_Trt, nrow = nrow, ncol = ncol, drop = TRUE)+
    scale_color_manual(
      values = c(
        "TRUE" = "red",
        "FALSE" = "blue"
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 8),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 8)
    )
  
  print(p)
}

dev.off()