rm(list = ls())
pacman::p_load(dplyr, magrittr, ggplot2, drake, knitr, ggpubr, RSvgDevice)
loadd(rna_kr_data)
loadd(dna_kr_data)

save_path <- "/Users/adnaniazi/Documents/phd/code/krauseNiazi2019Analyses/figures"

# CAP the data
rna_kr_data %<>% mutate(tail_length_tf = ifelse(tail_length_tf > 300, 300, tail_length_tf))
dna_kr_data %<>% mutate(tail_length_st = ifelse(tail_length_st > 300, 300, tail_length_st)) %>%
  mutate(tail_length_ff = ifelse(tail_length_ff > 300, 300, tail_length_ff))


# RNA ALONE
devtools::source_gist("2a1bb0133ff568cbe28d", filename = "geom_flat_violin.R")
p <- ggplot(data = rna_kr_data, aes(x = barcode, y = tail_length_tf, color=barcode, fill=barcode)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .5) +
  facet_grid(~barcode, scales = 'free') +
  geom_hline(aes(yintercept=as.numeric(as.character(barcode)))) +
  theme(panel.background = element_blank(),
        panel.grid.minor = ggplot2::element_line(colour="grey", size=0.5),
        panel.grid.major = ggplot2::element_line(colour="grey", size=0.5),
        panel.grid.major.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggplot2::scale_y_continuous(minor_breaks = seq(0, 300, 25), breaks = seq(0, 300, 100)) +
  ylab('tailfindr tail length estimate [nt]') 
RSvgDevice::devSVG(file.path(save_path, 'rna-alone.svg'), 
                   width = 10, height = 6.5,
                   bg = "white", fg = "black", onefile=TRUE, xmlHeader=TRUE)
print(p)
dev.off()

# RNA001 vs RNA002
devtools::source_gist("44fa0fdc48ddc495e1ed76e7524e4164", filename = "geom_flat_violin.R")
p <- ggplot(data = rna_kr_data, aes(x = barcode, y = tail_length_tf, color = kit)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .5) +
  facet_grid(~barcode, scales = 'free') +
  geom_hline(aes(yintercept=as.numeric(as.character(barcode)))) +
  theme(panel.background = element_blank(),
        panel.grid.minor = ggplot2::element_line(colour="grey", size=0.5),
        panel.grid.major = ggplot2::element_line(colour="grey", size=0.5),
        panel.grid.major.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggplot2::scale_y_continuous(minor_breaks = seq(0, 300, 25), breaks = seq(0, 300, 100)) +
  ylab('tailfindr tail length estimate [nt]') 
RSvgDevice::devSVG(file.path(save_path, 'rna001-vs-rna002.svg'), 
                   width = 10, height = 6.5,
                   bg = "white", fg = "black", onefile=TRUE, xmlHeader=TRUE)
print(p)
dev.off()


# 3) polyA vs PolyT 
devtools::source_gist("44fa0fdc48ddc495e1ed76e7524e4164", filename = "geom_flat_violin.R")
p <- ggplot(data = dna_kr_data, aes(x = barcode, y = tail_length_st, color = read_type)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .5) +
  facet_grid(~barcode, scales = 'free') +
  geom_hline(aes(yintercept=as.numeric(as.character(barcode)))) +
  theme(panel.background = element_blank(),
        panel.grid.minor = ggplot2::element_line(colour="grey", size=0.5),
        panel.grid.major = ggplot2::element_line(colour="grey", size=0.5),
        panel.grid.major.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggplot2::scale_y_continuous(minor_breaks = seq(0, 300, 25), breaks = seq(0, 300, 100)) +
  ylab('tailfindr tail length estimate [nt]') 
RSvgDevice::devSVG(file.path(save_path, 'polya-vs-polyt.svg'), 
                   width = 10, height = 6.5,
                   bg = "white", fg = "black", onefile=TRUE, xmlHeader=TRUE)
print(p)
dev.off()

# and 4) RNA vs DNA
polyat_dna_ab <- dplyr::select(dna_kr_data,
                               tail_length_st, read_type, barcode) %>% 
  mutate(barcode = as.numeric(as.character(barcode)))
polya_rna_ab <- dplyr::select(rna_kr_data,
                              tail_length_tf, barcode) %>% 
  mutate(barcode = as.numeric(as.character(barcode)))

# tabulate read in each DNA barcode
dna_barcode_nums <- as.data.frame(table(polyat_dna_ab$barcode))

dna_barcode_nums <- dplyr::rename(dna_barcode_nums,
                                  barcode=Var1,
                                  n_desired=Freq)

dna_barcode_nums <- dplyr::mutate(dna_barcode_nums,
                                  barcode=as.numeric(as.character(barcode)))

# randomize rna data
polya_rna_ab <- polya_rna_ab[sample(nrow(polya_rna_ab)),]

# inlclude only complete cases
polya_rna_ab <- polya_rna_ab[complete.cases(polya_rna_ab), ]

library(dplyr)
polya_rna_ab_subset <- dplyr::left_join(polya_rna_ab, dna_barcode_nums, by = "barcode") %>%
  dplyr::group_by(barcode) %>% 
  dplyr::slice(seq(dplyr::first(n_desired))) %>%
  dplyr::select(-n_desired)
table(polya_rna_ab_subset$barcode)

# combine rna and dna datasets
polya_rna_ab_subset <- dplyr::rename(polya_rna_ab_subset, tail_length = tail_length_tf)
polya_rna_ab_subset <- dplyr::mutate(polya_rna_ab_subset, data_type='RNA')
polyat_dna_ab <- dplyr::mutate(polyat_dna_ab, data_type='DNA') %>% 
  select(-read_type) %>% 
  rename(tail_length = tail_length_st)

polyat_dna_ab$barcode <- as.numeric(polyat_dna_ab$barcode)

df <- bind_rows(polyat_dna_ab, polya_rna_ab_subset)
df <- dplyr::mutate(df, 
                    data_type=as.factor(data_type),
                    barcode=as.factor(barcode)
)

devtools::source_gist("44fa0fdc48ddc495e1ed76e7524e4164", filename = "geom_flat_violin.R")
p <- ggplot(data = df, aes(x = barcode, y = tail_length, color = data_type)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .5) +
  facet_grid(~barcode, scales = 'free') +
  geom_hline(aes(yintercept=as.numeric(as.character(barcode)))) +
  theme(panel.background = element_blank(),
        panel.grid.minor = ggplot2::element_line(colour="grey", size=0.5),
        panel.grid.major = ggplot2::element_line(colour="grey", size=0.5),
        panel.grid.major.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggplot2::scale_y_continuous(minor_breaks = seq(0, 300, 25), breaks = seq(0, 300, 100)) +
  ylab('tailfindr tail length estimate [nt]') 
RSvgDevice::devSVG(file.path(save_path, 'rna-vs-dna.svg'), 
                   width = 10, height = 6.5,
                   bg = "white", fg = "black", onefile=TRUE, xmlHeader=TRUE)
print(p)
dev.off()


## 6. Nanopolish's tail start estimate vs. tailfindr tail start estimate
data <- mutate(rna_kr_data, diff = tail_start_np - tail_start_tf)
p <- ggplot(data, aes(x = diff)) + 
  geom_histogram(binwidth = 1) + 
  theme_bw() + 
  coord_cartesian(xlim = c(-500, 500)) +
  scale_x_continuous(minor_breaks = seq(-500, 500, 50)) +
  xlab('Nanopolish tail start - tailfindr tail start') +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_blank(),
                 panel.background = element_rect(fill = "transparent") # bg of the panel
                 , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
                 , panel.grid.minor = element_line(size = 0.7)
                 , panel.grid.major = element_line(size = 0.7)
                 , legend.background = element_rect(fill = "transparent") # get rid of legend bg
                 , legend.box.background = element_rect(fill = "transparent")) # get rid of legend panel bg) 
ggplot2::ggsave(plot=p, filename=file.path(save_path, 'np_start_vs_tf_start.png'), height  = 7, width = 7, units = 'in', dpi=600)


## 6. Nanopolish's tail end estimate vs. tailfindr tail end estimate
data <- mutate(rna_kr_data, diff = tail_end_np - tail_end_tf)
p <- ggplot(data, aes(x = diff)) + 
  geom_histogram(binwidth = 1) + 
  theme_bw() + 
  coord_cartesian(xlim = c(-500, 500)) +
  scale_x_continuous(minor_breaks = seq(-500, 500, 50)) +
  xlab('Nanopolish tail start - tailfindr tail start') +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                 panel.border = ggplot2::element_blank(),
                 panel.background = element_rect(fill = "transparent") # bg of the panel
                 , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
                 , panel.grid.minor = element_line(size = 0.7)
                 , panel.grid.major = element_line(size = 0.7)
                 , legend.background = element_rect(fill = "transparent") # get rid of legend bg
                 , legend.box.background = element_rect(fill = "transparent")) # get rid of legend panel bg) 
ggplot2::ggsave(plot=p, filename=file.path(save_path, 'np_end_vs_tf_end.png'), height  = 7, width = 7, units = 'in', dpi=600)




