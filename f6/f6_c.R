library(tidyverse)
library(rtracklayer)
library(ggtext)

bins <- import.bed('../data/1kb.bed')
bl <- import.bed('../data/blacklist.bed')
k <- !overlapsAny(bins, bl) & 
  as.character(seqnames(bins)) %in% paste0('chr', c(1:19, 'X'))

cts <- cbind(readRDS('../data/cts.rds')$`1kb`[k,],
             readRDS('../data/ferrari.cts.rds')$`1kb`[k,])
lsz <- rbind(read_csv('../data/libsz.csv'),
             read_csv('../data/ferrari.libsz.csv')) 

clrs2 <- c('WT' = '#0C8140', 'EED-KO' = '#F26522', 'sgNSD1' = '#EE2E2E')

grep('27ac', colnames(cts), value = T) %>%
  tibble(samp = .) %>%
  merge(lsz, by = 'samp', all.x = T) %>%
  rowwise() %>%
  mutate(v = mean(cts[[samp]][cts[[samp]] > quantile(cts[[samp]], .99)])) %>%
  ungroup() %>%
  mutate(study = ifelse(grepl('E36', samp), 'Ferrari 2014', 'This study') %>%
           factor(c('This study', 'Ferrari 2014')),
         cond = case_when(grepl('NSD1', samp) ~ 'sgNSD1',
                          grepl('EedKO', samp) ~ 'EED-KO',
                          T ~ 'WT'),
         v = 1e6 * v / endo.chip,
         clr = clrs2[cond],
         nm = sprintf("<span style='color:%s'>%s</span>",
                      clr, cond)) %>%
  ggplot(aes(y = nm, x = v, fill = clr)) +
  geom_col() +
  facet_grid(study ~ ., scales = 'free_y') +
  scale_fill_identity() +
  scale_x_continuous('H3K27ac peakiness', expand = expansion(c(0, .05))) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text.y = element_markdown(),
        legend.direction = 'vertical',
        axis.title.y = element_blank(),
        strip.text = element_text(color = 'black', size = 13),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = 'grey80', linetype = 'dashed')) -> p

ggsave('f6_c.pdf', p, height = 4, width = 3, bg = 'transparent')

