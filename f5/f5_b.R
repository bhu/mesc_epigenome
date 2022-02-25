library(tidyverse)
library(rtracklayer)
library(gggenes)
library(patchwork)
library(pals)
library(ggh4x)
library(ggrepel)

# region
reg <- GRanges('chr4', IRanges(32.5e6,33.6e6))
ttl <- 'chr4: 32.5 - 33.6mb'
xbrks <- seq(32.5e6, 33.6e6, 5e5)
xlabs <- sprintf('%.1fmb', xbrks / 1e6)
xlims <- c(start(reg), end(reg))

# genes
g2u <- read_tsv('../data/GCF_000001635.26_GRCm38.p6_assembly_report.txt', 
                comment = '#', col_names = F) %>%
  select(X7, X10) %>%
  filter(X10 %in% paste0('chr', c(1:22, 'X', 'Y'))) %>%
  deframe()
g <- import.gff3('../data/GCF_000001635.26_GRCm38.p6_genomic.gff.gz') %>%
  {.[!is.na(.$tag) & .$tag == 'RefSeq Select']} %>%
  {.[seqnames(.) %in% names(g2u)]} %>%
  renameSeqlevels(g2u) %>%
  keepSeqlevels(g2u)

# import bigwigs
odr <- c('H3K36me2', 'H3K27me3', 'mCG')
clrs2 <- c('WT' = '#0C8140','DNMT-TKO' = '#3953A4',
           'EZH2-KO' = '#F26522', 'EED-KO' = '#F26522')

bws <- list.files('../tracks', full.names = T, pattern = '27me3|36me2|mCG') %>%
  grep('NSD|Par_H3|rep1|GSM|E36', ., value = T, invert = T) %>%
  setNames(., sub('\\..*', '', basename(.))) %>%
  lapply(function(x) {
    import.bw(x, selection = BigWigSelection(reg))
  })
bins <- tile(reg, 300)[[1]]
xs <- mid(bins)
scs <- lapply(bws, function(x) {
  if (all(width(x) > 1)) {
    keepSeqlevels(x, seqnames(reg)) %>%
      mcolAsRleList('score') %>%
      binnedAverage(bins, ., 'score') %>%
      score() %>%
      tibble(score = .) %>%
      mutate(idx = xs)
  } else {
    findOverlaps(bins, x) %>%
      as("List") %>%
      extractList(x$score, .) %>%
      mean(na.rm = T) %>%
      tibble(score = .) %>%
      mutate(idx = xs)
  }
}) %>% 
  bind_rows(.id = 'track') %>%
  mutate(score = replace_na(score, 0),
         cnd = case_when(grepl('Par|A1', track) ~ 'WT_A1',
                         grepl('TKO', track) ~ 'DNMT-TKO',
                         grepl('NSD', track) ~ 'sgNSD1',
                         grepl('Eed', track) ~ 'EED-KO',
                         grepl('EZH2', track) ~ 'EZH2-KO',
                         T ~ 'WT_J1'),
         mark = case_when(grepl('36me2', track) ~ 'H3K36me2',
                          grepl('27me3', track) ~ 'H3K27me3',
                          T ~ 'mCG'))

ms <- read_csv('../../massspec.facs.csv') %>%
  mutate(cnd = c('mESC_A1_par' = 'WT_A1',
                 'mESC_A1_sgNSD1' = 'sgNSD1',
                 'mESC_J1_par' = 'WT_J1',
                 'mESC_J1-DNMT-TKO' = 'DNMT-TKO',
                 'mESC_A1_EZH2-KO' = 'EZH2-KO')[type])

pd <- read_csv('../data/libsz.csv') %>%
  mutate(rx = (endo.chip / endo.input) / (exo.chip / exo.input)) %>%
  merge(scs, all.y = T, by.x = 'samp', by.y = 'track') %>%
  mutate(ttl = ttl,
         cpm = replace_na(1e6 / endo.chip, 1) * score,
         dm = replace_na(1e6 / exo.chip, 1) * score,
         rx = replace_na(1e6 * rx / endo.chip, 1) * score,
         cond = case_when(grepl('A1', samp) ~ 'A1 WT',
                          grepl('J1', samp) ~ 'J1 WT',
                          T ~ cnd),
         cndd = sub('_.*', '', cnd),
         track = paste(cnd, mark),
         samp = sub('_rep[12]$', '', samp),
         chip = factor(samp, c('mESC_WT_J1_H3K36me2',
                               'mESC_DNMT_TKO_H3K36me2',
                               'mESC_WT_A1_H3K36me2',
                               'mESC_EZH2_KO_H3K36me2',
                               'mESC_WT_J1_H3K27me3',
                               'mESC_DNMT_TKO_H3K27me3',
                               'mESC_WT_A1_H3K27me3',
                               'mESC_EZH2_KO_H3K27me3',
                               'mESC_Par_mCG',
                               'mESC_EedKO_mCG'))) %>%
  merge(ms, by = c('cnd', 'mark'), all.x = T) %>%
  arrange(chip) %>%
  mutate(track = fct_inorder(track),
         ms = replace_na(1e6 * f / endo.chip, 1) * score)

yls <- c('cpm', 'dm', 'rx', 'ms') %>%
  setNames(.,.) %>%
  lapply(function(s) {
    list(1:4, 5:8, 9:10) %>%
      lapply(function(i) {
        s <- pd[[s]][pd$track %in% levels(pd$track)[i]]
        s <- s[is.finite(s)]
        m <- max(s)
        l <- if(m >= 10) {
          floor(m)
        } else if (m >= 1) {
          floor(m * 10) / 10
        } else {
          floor(m * 100) / 100
        }
        scale_y_continuous(limits = c(0, m),
                           breaks = c(0, l),
                           labels = c('', l),
                           expand = expansion(c(0, .1)),
                           position = 'right') %>%
          list(.) %>%
          rep(length(i))
      }) %>%
      unlist(recursive = F, use.names = F)
  })

# tracks
nrm <- 'ms'
pd %>%
  dplyr::rename(sc = !!nrm) %>%
  ggplot(aes(x = idx, y = sc)) +
  scale_fill_manual(values = clrs2) +
  geom_area(aes(fill = cndd), alpha = 1, color = NA) +
  scale_x_continuous(expand = expansion(0),
                     breaks = xbrks,
                     limits = xlims) +
  facet_grid2(track ~ ., scales = 'free_y', switch = 'y',
              strip = strip_themed(
                text_y = elem_list_text(color = unname(clrs2[c(1,2,1,3,1,2,1,3, 1,3)]))
              )) +
  facetted_pos_scales(y = yls[[nrm]]) +
  ylab('MS norm coverage or mCG%') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        plot.margin = margin(),
        axis.title.y.right = element_text(margin = margin(0,0,0,10)),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 11),
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        axis.text = element_text(color = 'black', size  = 11)) -> p1

# cgis
read_tsv('../data/cpgIslandExt.txt.gz', col_names = F) %>% 
  dplyr::select(chr = X2, start = X3, end = X4) %>% 
  filter(chr == as.character(seqnames(reg)) &
           start > start(reg) &
           end < end(reg)) %>%
  mutate(mid = (start + end) / 2) %>%
  ggplot() +
  geom_linerange(aes(x = mid, ymin = 0, ymax = 1), size = 2) +
  #geom_rect(size = 2, fill = 'black', color = 'black') +
  scale_x_continuous(expand = expansion(0),
                     breaks = xbrks,
                     limits = xlims) +
  scale_y_continuous(position = 'right', expand = expansion(.2)) +
  facet_grid('CpG island'~., switch = 'y') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        plot.margin = margin(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color = 'black', size = 11),
        strip.text.y.left = element_text(angle = 0, hjust = 1)) -> p2

# genes
subsetByOverlaps(g, reg) %>% 
  as_tibble() %>%
  filter(type %in% c('mRNA', 'exon')) %>%
  mutate(forward = strand == '+',
         position = case_when(strand == '-' ~ end, T ~ start),
         mid = (start + end) / 2,
         alt = as.numeric(fct_inorder(gene)) %% 2 == 0) %>%
  ggplot(aes(xmin = start, xmax = end, y = 1, color = strand,
             fill = strand, forward = forward)) +
  geom_text_repel(aes(x = position, y = 1, label = gene), 
                  force_pull = 0,
                  nudge_y = -5,
                  direction = 'x',
                  hjust = .5,
                  size = 3,
                  segment.color = 'grey70',
                  segment.linetype = 'dashed',
                  min.segment.length = 0,
                  data = ~subset(., type == 'mRNA' & alt)) +
  geom_text_repel(aes(x = position, y = 1, label = gene), 
                  force_pull = 0,
                  nudge_y = 5,
                  direction = 'x',
                  hjust = .5,
                  size = 3,
                  segment.color = 'grey70',
                  segment.linetype = 'dashed',
                  min.segment.length = 0,
                  data = ~subset(., type == 'mRNA' & !alt)) +
  geom_linerange(data = ~subset(., type == 'mRNA')) +
  geom_gene_arrow(color = NA, data = ~subset(., type != 'mRNA')) +
  geom_feature(aes(x = position, y = 1, forward = forward, color = strand),
               data = ~subset(., type == 'mRNA')) +
  scale_color_manual(values = c('black', 'black')) +
  scale_fill_manual(values = c('black', 'black')) +
  scale_x_continuous(expand = expansion(0),
                     breaks = xbrks,
                     labels = xlabs) +
  xlab(seqnames(reg)) +
  scale_y_continuous(limits = c(.8,1.2)) +
  coord_cartesian(xlim = xlims, clip = 'on') +
  facet_grid('Gene'~., switch = 'y') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color= 'black', size = 11),
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        axis.line.x = element_line(color = 'black'),
        panel.grid = element_blank(),
        axis.text =element_text(color = 'black', size = 11),
        plot.margin = margin(t = 0)) -> p3

{ wrap_plots(p1, p2, p3, ncol = 1, heights = c(8, 1, 2))  &
    theme(plot.background = element_blank()) } %>%
  ggsave('f5_b.pdf', ., height = 6.5, width = 7, bg = 'transparent')
