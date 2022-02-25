library(tidyverse)
library(rtracklayer)
library(ggtext)

load('../data/lavarone.cts.rda')

bins <- bins$`1kb` %>%
  mutate(start = start + 1) %>%
  makeGRangesFromDataFrame()
bl <- import.bed('../data/blacklist.bed')
k <- !overlapsAny(bins, bl) & 
  as.character(seqnames(bins)) %in% paste0('chr', c(1:19, 'X'))

cts <- cts$`1kb`[k,]
clrs2 <- c('WT' = '#0C8140', 'EZH2-KO' = '#EC7723', 'EZH1/2-DKO' = '#CC6828')

grep('27ac', colnames(cts), value = T) %>%
  tibble(samp = .) %>%
  mutate(tot = lsz[samp]) %>%
  rowwise() %>%
  mutate(v = mean(cts[[samp]][cts[[samp]] > quantile(cts[[samp]], .99)])) %>%
  ungroup() %>%
  mutate(cond = case_when(grepl('Ezh1', samp) ~ 'EZH1/2-DKO',
                          grepl('Ezh2', samp) ~ 'EZH2-KO', 
                          T ~ 'WT') %>%
           factor(names(clrs2)),
         v = 1e6 * v / tot,
         clr = clrs2[cond],
         nm = sprintf("<span style='color:%s'>%s H3K27ac</span>",
                      clr, cond)) %>%
  ggplot(aes(y = nm, x = v, fill = clr)) +
  geom_col() +
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
p
ggsave('sf4_c.pdf', p, height = 2, width = 6.2, bg = 'transparent')
