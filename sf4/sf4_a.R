library(tidyverse)
library(rtracklayer)
library(gggenes)
library(patchwork)
library(pals)
library(ggh4x)
library(ggrepel)

# region
reg <- GRanges('chr6', IRanges(107.81e6,109.26e6))
ttl <- 'chr6: 107.81 - 109.26mb'
xbrks <- seq(108e6, 109e6, 5e5)
xlabs <- sprintf('%.1fmb', xbrks / 1e6)
xlims <- c(start(reg), end(reg))
rcts <- tibble(start = c(107.81e6,108553727,108864728),
               end = c(108067494,108783800,109.26e6),
               kind = c('inactive', 'inactive', 'inactive'))

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
odr <- c('H3K27ac','H3K27me1','H3K27me2')
clrs2 <- c('WT' = '#0C8140', 'EZH2-KO' = '#EC7723', 'EZH1/2-DKO' = '#CC6828')
clrs <- c('intergenic' = '#FBB040',
          'active' = '#1C75BC',
          'inactive' = '#7F3F98')

bws <- list.files('../tracks', full.names = T, pattern = 'K27') %>%
  grep('E14_H3K27ac|Ezh1KOEzh2KO_H3K27ac|_Ezh2KO_', ., value = T) %>%
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
         cond = case_when(grepl('Ezh1', track) ~ 'EZH1/2-DKO',
                          grepl('Ezh2', track) ~ 'EZH2-KO', 
                          T ~ 'WT') %>%
           factor(names(clrs2)), 
         mark = sub('.*_', '', track) %>%
           factor(odr)) %>%
  arrange(mark, cond) %>%
  mutate(track = paste(cond, mark) %>% fct_inorder())

yl <- split(scs$score, scs$track) %>%
  lapply(function(s) {
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
                       position = 'right') 
  })

# tracks
scs %>%
  ggplot(aes(x = idx, y = score)) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = start, xmax = end, fill = kind),
            data = rcts, inherit.aes = F, alpha = .2) +
  scale_fill_manual(values = c(clrs, clrs2)) +
  geom_area(aes(fill = cond), alpha = 1, color = NA) +
  scale_x_continuous(expand = expansion(0),
                     breaks = xbrks,
                     limits = xlims) +
  facet_grid2(track ~ ., scales = 'free_y', switch = 'y',
              strip = strip_themed(
                text_y = elem_list_text(color = unname(clrs2[c(1,2,3,2,2)]))
              )) +
  facetted_pos_scales(y = yl) +
  ylab('Coverage') +
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
         mid = (start + end) / 2) %>%
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
                  data = ~subset(., type == 'mRNA')) +
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
  scale_y_continuous(limits = c(.8,1.1)) +
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

{ wrap_plots(p1, p2, p3, ncol = 1, heights = c(11, 1, 2))  &
    theme(plot.background = element_blank()) } %>%
  ggsave('sf4_a.pdf', ., height = 5, width = 5.5, bg = 'transparent')

