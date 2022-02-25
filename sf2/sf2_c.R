library(tidyverse)
library(rtracklayer)
library(ggtext)
library(ggh4x)

bins <- import.bed('../data/1kb.bed')
bl <- import.bed('../data/blacklist.bed')
k <- !overlapsAny(bins, bl) & 
  as.character(seqnames(bins)) %in% paste0('chr', c(1:19, 'X'))

cts <- readRDS('../data/cts.rds')$`1kb`[k,]
lsz <- read_csv('../data/libsz.csv') 

clrs2 <- c('WT' = '#0C8140','DNMT-TKO' = '#3953A4', 'sgNSD1' = '#EE2E2E')

grep('27me3', colnames(cts), value = T) %>%
  grep('Par|sgNSD1|rep1', ., value = T) %>%
  tibble(samp = .) %>%
  merge(lsz, by = 'samp', all.x = T) %>%
  rowwise() %>%
  mutate(v = mean(cts[[samp]][cts[[samp]] > quantile(cts[[samp]], .99)])) %>%
  ungroup() %>%
  dplyr::select(samp, v, tot = endo.chip) %>%
  mutate(cond = case_when(grepl('NSD1', samp) ~ 'sgNSD1',
                          grepl('TKO', samp) ~ 'DNMT-TKO',
                          T ~ 'WT'),
         v = 1e6 * v / tot,
         clr = clrs2[cond],
         nm = sprintf("<span style='color:%s'>%s H3K27me3</span>",
                      clr, cond),
         proj = 'This study') -> d1
  

load('../data/TKO.K27me3.public.rda')

bins <- bins$`1kb` %>%
  mutate(start = start + 1) %>%
  makeGRangesFromDataFrame()
k <- !overlapsAny(bins, bl) & 
  as.character(seqnames(bins)) %in% paste0('chr', c(1:19, 'X'))

cts <- cts$`1kb`[k,]


colnames(cts) %>%
  tibble(samp = .) %>%
  mutate(tot = lsz[samp]) %>%
  rowwise() %>%
  mutate(v = mean(cts[[samp]][cts[[samp]] > quantile(cts[[samp]], .99)])) %>%
  ungroup() %>%
  dplyr::select(samp, v, tot) %>%
  mutate(cond = case_when(grepl('TKO', samp) ~ 'DNMT-TKO',
                          T ~ 'WT'),
         v = 1e6 * v / tot,
         clr = clrs2[cond],
         nm = sprintf("<span style='color:%s'>%s H3K27me3</span>",
                      clr, cond),
         proj = ifelse(grepl('rep', samp), 'Hagarman', 'Brinkman')) -> d2


rbind(d1, d2) %>%
  mutate(line = case_when(proj == 'This study' & !grepl('DNMT', samp) ~'A1',
                          proj == 'This study' ~ 'J1',
                          T ~ 'V6.5')) %>%
  filter(line != 'J1') %>%
  ggplot(aes(x = nm, y = v, fill = clr)) +
  stat_summary(geom = 'errorbar', fun.data = mean_se, width = .25) +
  stat_summary(geom = 'col', fun = mean) +
  facet_nested(.~proj, space = 'free', scales = 'free_x',
               nest_line = element_line(color = 'black')) +
  scale_fill_identity() +
  scale_y_continuous('H3K27me3 peakiness', expand = expansion(c(0, .05))) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = 'grey80', linetype = 'dashed'),
        axis.text = element_text(size = 11, color = 'black'),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(angle = 45, hjust = 1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 13, color = 'black'),
        legend.position = 'none') -> p
ggsave('sf2_c.pdf', p, height = 3.5, width = 5.8, bg = 'transparent')
