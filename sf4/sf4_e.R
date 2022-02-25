library(tidyverse)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

ann <- import('../data/mESC_par_H3K27ac_peaks.narrowPeak') %>%
  annotatePeak(TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
               annoDb = 'org.Mm.eg.db')

enh <- ann@anno[ann@anno$annotation == 'Distal Intergenic']

bins <- import('../data/1kb.bed')
bl <- import('../data/blacklist.bed')
k <- !overlapsAny(bins, bl) & 
  as.character(seqnames(bins)) %in% paste0('chr', c(1:19, 'X'))
cts <- readRDS('../data/cts.rds')$`1kb`[k,]
bins <- bins[k,]
lsz <- read_csv('../data/libsz.csv')
ms <- read_csv('../data/massspec.facs.csv') %>%
  filter(grepl('A1_par|sgNSD1', type) & grepl('36me2', type)) %>%
  mutate(cond = ifelse(grepl('NSD1', type), 'sgNSD1', 'Par'))

cts[overlapsAny(bins, enh), grepl('36me2|27ac', names(cts)) & grepl('Par|sg', names(cts))] %>%
  mutate(idx = 1:n()) %>%
  pivot_longer(-idx, names_to = 'samp', values_to = 'v') %>%
  merge(lsz) %>%
  mutate(rx = (endo.chip / exo.chip) / (endo.input / exo.input),
         v = v / endo.chip) %>%
  separate(samp, c(NA, 'cond', 'mark')) %>%
  merge(ms, all.x = T) %>%
  mutate(v = case_when(is.na(a) ~ v * rx,
                       T ~ v * a)) %>%
  dplyr::select(idx, cond, mark, v) %>%
  pivot_wider(names_from = 'cond', values_from = 'v') %>%
  mutate(LFC = log2((sgNSD1) / (Par))) %>%
  dplyr::select(idx, LFC, mark) %>%
  filter(is.finite(LFC)) %>%
  pivot_wider(names_from = 'mark', values_from = 'LFC') %>%
  ggplot(aes(x = H3K27ac, y = H3K36me2)) +
  geom_point(shape = 16, size = .9, alpha = .9) +
  geom_smooth(method = 'lm') +
  facet_wrap(~'Epigenetic changes at distal enhancers') +
  labs(x = expression(log[2]*frac('sgNSD1', 'WT')*' H3K27ac'),
       y = expression(log[2]*frac('sgNSD1', 'WT')*' H3K6me2'),) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        strip.text = element_text(size = 13, color = 'black'),
        axis.text = element_text(size = 11, color = 'black'),
        strip.background = element_rect(fill = NA)) -> p

ggsave('sf4_e.pdf', p, height = 4, width = 4.2, bg = 'transparent')

