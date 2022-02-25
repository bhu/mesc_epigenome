library(tidyverse)
library(rtracklayer)
library(ggdist)
library(ggsignif)
library(ggh4x)
library(ggtext)

bl <- import.bed('../data/blacklist.bed')
bins <- import.bed('../data/10kb.bed')
k <- !overlapsAny(bins, bl) & 
  as.character(seqnames(bins)) %in% paste0('chr', c(1:19, 'X'))

cts <- readRDS('../data/cts.rds')$`10kb`[k,]
bins <- bins[k,]

bins$cg <- read_tsv('../data/10kb.cg.bed', 
               col_names = c('seqnames', 'start', 'end', 'cg')) %>%
  mutate(start = start + 1) %>%
  merge(as_tibble(bins), ., all.x = T, by = c('seqnames', 'start', 'end')) %>%
  mutate(seqnames = factor(seqnames, seqlevels(bins))) %>%
  arrange(seqnames, start) %>%
  pull(cg)

lsz <- read_csv('../data/libsz.csv') 

clrs2 <- c('WT' = '#0C8140','DNMT-TKO' = '#3953A4', 'sgNSD1' = '#EE2E2E')
nms <- clrs2 %>%
  mapply(function(clr, nm) {
    sprintf("<span style='color:%s'>%s</span>", clr, nm)
  }, ., names(.))


grep('Par|sgNSD1|J1|rep2', colnames(cts), value = T) %>%
  grep('27me3', ., value = T) %>%
  tibble(samp = .) %>%
  merge(lsz, by = 'samp', all.x = T) %>%
  group_by(samp) %>%
  do(., tibble(v = cts[[.$samp]] * 100 / (.$endo.chip / length(bins))) %>%
       mutate(idx = 1:n())) %>%
  ungroup() %>%
  mutate(cond = case_when(grepl('TKO', samp) ~ 'DNMT-TKO',
                          grepl('NSD1', samp) ~ 'sgNSD1',
                          T ~ 'WT'),
         line = ifelse(grepl('J1|TKO', samp), 'J1', 'A1')) %>%
  split(., .$line) %>%
  lapply(function(x) {
    x %>%
      mutate(cond = ifelse(cond == 'WT', 'WT', 'KO')) %>%
      dplyr::select(cond, idx, v) %>%
      pivot_wider(names_from = 'cond', values_from = 'v') %>%
      mutate(d = log((KO + 1)/(WT + 1)),
             perc = ntile(d, 100)) %>%
      filter(perc %in% c(1, 100)) 
  }) %>%
  bind_rows(.id = 'line') %>%
  mutate(comp = ifelse(line == 'J1', 'DNMT-TKO', 'sgNSD1'),
         nm = nms[comp],
         cg = bins$cg[idx],
         grp = ifelse(perc == 1, '< WT (loss)', '> WT (gain)') %>%
           factor(c('> WT (gain)', '< WT (loss)'))) %>%
  ggplot(aes(x = nm, y = cg, fill = comp, color = comp)) +
  geom_violin(color = NA, alpha = .5) +
  stat_pointinterval() +
  geom_signif(comparisons = list(nms[c('DNMT-TKO', 'sgNSD1')]),
              tip_length = 0, y_position = 500, color = 'black') +
  facet_nested(.~ 'H3K27me3' + grp,
               nest_line = element_line(color = 'black')) +
  coord_cartesian(ylim = c(0, 550)) +
  scale_fill_manual(values = clrs2) +
  scale_color_manual(values = clrs2) +
  ylab('CpG density / 10kb') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 13, color = 'black'),
        panel.grid = element_blank(),
        legend.position = 'none',
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        axis.text = element_text(size = 11, color = 'black'),
        axis.text.x = element_markdown(angle = 45, hjust = 1),
        axis.title.x = element_blank()) -> p

ggsave('f4_b.pdf', p, height = 4.3, width = 3, bg = 'transparent')
