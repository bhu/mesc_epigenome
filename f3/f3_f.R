library(tidyverse)
library(rtracklayer)
library(ggtext)

bins <- import.bed('../data/1kb.bed')
bl <- import.bed('../data/blacklist.bed')
k <- !overlapsAny(bins, bl) & 
  as.character(seqnames(bins)) %in% paste0('chr', c(1:19, 'X'))

cts <- readRDS('../data/cts.rds')$`1kb`[k,]
lsz <- read_csv('../data/libsz.csv') 

clrs2 <- c('WT' = '#0C8140','DNMT-TKO' = '#3953A4', 'sgNSD1' = '#EE2E2E')

grep('27me3', colnames(cts), value = T) %>%
  grep('WT|NSD1_KO|rep2', ., value = T) %>%
  tibble(samp = .) %>%
  merge(lsz, by = 'samp', all.x = T) %>%
  rowwise() %>%
  mutate(v = mean(cts[[samp]][cts[[samp]] > quantile(cts[[samp]], .99)])) %>%
  ungroup() %>%
  mutate(line = ifelse(grepl('NSD1|A1', samp), 'A1', 'J1') %>%
           paste('mESC'),
         cond = case_when(grepl('NSD1', samp) ~ 'sgNSD1',
                          grepl('TKO', samp) ~ 'DNMT-TKO',
                          T ~ 'WT'),
         v = 1e6 * v / endo.chip,
         clr = clrs2[cond],
         nm = sprintf("<span style='color:%s'>%s H3K27me3</span>",
                      clr, cond)) %>%
  ggplot(aes(y = nm, x = v, fill = clr)) +
  geom_col() +
  facet_grid(line ~ ., scales = 'free_y') +
  scale_fill_identity() +
  scale_x_continuous('H3K27me3 peakiness', expand = expansion(c(0, .05))) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = 'grey80', linetype = 'dashed'),
        axis.text = element_text(size = 11, color = 'black'),
        axis.title.y = element_blank(),
        axis.text.y = element_markdown(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 13, color = 'black'),
        legend.position = 'none') -> p
ggsave('f3_f.pdf', p, height = 3.6, width = 5.5, bg = 'transparent')
