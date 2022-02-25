library(data.table)
library(tidyverse)
library(rtracklayer)
library(matrixStats)
library(isoband)
library(sf)
library(MASS)
library(lwgeom)
library(hexbin)
library(pals)
library(patchwork)

cts <- readRDS('../data/cts.rds')$`10kb`
bins <- fread('../data/10kb.bed', col.names = c('chr', 'start', 'end')) %>%
  mutate(start = start + 1) %>%
  makeGRangesFromDataFrame()

greg <- import.bed('../data/gene.bed')
ireg <- import.bed('../data/intergenic.bed')
i1 <- overlapsAny(bins, greg) & !overlapsAny(bins, ireg)
i2 <- !overlapsAny(bins, greg) & overlapsAny(bins, ireg)
olap <- tibble(gene = overlapsAny(bins, greg),
               igr = overlapsAny(bins, ireg)) %>%
  mutate(out = case_when(
    gene & !igr ~ 1,
    !gene & igr ~ -1,
    TRUE ~ 0
  )) %>%
  pull(out)

cfunc <- colorRampPalette(c("#93003a","#ad1042","#c32a49",
                            "#d4444d","#e35e50","#ee7851",
                            "#f59350","#faad4c","#fbc846",
                            "#f9e43c","#f3ff2c","#cbed92",
                            "#b2daa9","#9cc6b3","#88b2b8",
                            "#759fb8","#638bb6","#5178b2",
                            "#3e66ac","#2854a5","#00429d"))
lineclr <- "black"
horz <- F
gradientn1 <- cfunc(21)
cramp <- colorRampPalette(c("#000000ff","#ffffff00"), alpha = T)(5)
leg.brks <- seq(-1, 1, length.out = 19)[seq(2, 18, by = 2)]
leg.labs <- c(sprintf('Genic\u25bc'), rep('', 3), '50%',
              rep('', 3), sprintf('Genic\u25b2'))
len <- 9
pal <- cfunc(len)
cmat <- seq(0, 255, length.out = len + 1) %>%
  {.[-1]} %>%
  round() %>%
  as.hexmode() %>%
  format(width = 2, upper.case = T) %>%
  lapply(function(x) {
    paste0(pal, x)
  }) %>%
  do.call(cbind, .) %>%
  data.frame() %>%
  `colnames<-`(1:len) %>%
  mutate(clr = 1:dplyr::n()) %>%
  reshape2::melt(id.vars = "clr", variable.name = "opa") %>%
  mutate(opa = as.integer(opa))
leg <- ggplot() +
  geom_tile(aes(x = opa, y = clr, fill = value),
            data = cmat) +
  scale_fill_identity() +
  labs(x = "# of bins \u25ba",
       y = "% genic \u25ba") +
  coord_fixed(expand = F) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Arial',
                                  color = "black"))
leg

ggsave('f5_c.leg.pdf', leg, device = cairo_pdf, height = 1.2, width = 1.5)


lsz <- read_csv('../data/libsz.csv') 

ms <- read_csv('../../massspec.facs.csv') %>%
  mutate(cond = c('mESC_A1_par' = 'Par',
                  'mESC_A1_EZH2-KO' = 'EZH2-KO')[type]) %>%
  na.omit()

grep('A1|EZH2', colnames(cts), value = T) %>%
  grep('27me3|36me2', ., value = T) %>%
  tibble(samp = .) %>%
  mutate(cond = ifelse(grepl('EZH', samp), 'EZH2-KO', 'Par'),
         mark = sub('.*_', '', samp)) %>%
  merge(lsz, by = 'samp', all.x = T) %>%
  mutate(rx = (endo.chip / endo.input) / (exo.chip / exo.input)) %>%
  merge(ms, by = c('cond', 'mark'), all.X = T) %>%
  group_by(cond, mark, samp) %>%
  do(., {
    x <- .
    #tibble(v = 1e6 * cts[[.$samp]] * .$f / .$endo.chip) %>%
    tibble(v = cts[[.$samp]] * .$a * 100 / mean(cts[[.$samp]])) %>%
      mutate(idx = 1:n())
  }) %>%
  ungroup() %>%
  dplyr::select(cond, v, mark, idx) %>%
  pivot_wider(names_from = 'cond', values_from = 'v') %>%
  dplyr::select(x = Par, y = `EZH2-KO`, idx, mark) %>%
  mutate(mark = factor(mark, c('H3K27me3', 'H3K36me2'))) %>%
  split(., .$mark) %>%
  lapply(function(d) {
    pdat <- d %>%
      filter(x != 0 | y != 0) %>%
      mutate(across(c(x, y), ~log10(.x))) %>%
      filter(is.finite(x) & is.finite(y)) %>%
      filter(x < 2 & y < 2) %>%
      mutate(r = olap[idx]) %>%
      dplyr::select(x, y, r) %>%
      rbind(tibble(x = c(0, 0, 2, 2),
                   y = c(0, 0, 2, 2),
                   r = rep(0, 4)))
    hex <- hexbin(pdat$x, pdat$y, xbins = 50, IDs = T)
    pdat$cell <- hex@cID
    hex <- data.frame(hcell2xy(hex),
                      cell = hex@cell,
                      count = hex@count)
    pdat2 <- pdat %>%
      group_by(cell) %>%
      summarise(prop = mean(r, na.rm = T)) %>%
      ungroup() %>%
      right_join(hex, by = "cell") %>%
      mutate(logcount = log10(count),
             mark = d$mark[1])
    
    pow <- .5
    ggplot() +
      geom_hex(aes(x = x, y = y, fill = prop, alpha = count^pow, color = prop),
               stat = "identity", color = NA, data = pdat2, size = 5) +
      scale_fill_gradientn(colors = gradientn1, name = "% Genic",
                           breaks = leg.brks, labels = leg.labs,
                           limits = c(-1, 1)) +
      scale_color_gradientn(colors = gradientn1, name = "% Genic",
                            breaks = leg.brks, labels = leg.labs,
                            limits = c(-1, 1)) +
      scale_alpha(range = c(0.01, 1), name = "Number of bins", guide = F) +
      coord_cartesian(expand = F) +
      scale_y_continuous('EZH2-KO', breaks = -5:2,
                         labels = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10,  100)) +
      scale_x_continuous('WT', breaks = 0:2,
                         labels = c(1,  10, 100)) +
      annotate('segment', x = -Inf, xend =Inf, y = Inf,
               yend = Inf, color = 'white') +
      geom_abline(slope = 1, intercept = 0, color = 'grey50', alpha = .5) +
      facet_wrap(~mark) +
      theme(panel.background = element_rect(fill =NA, color = 'black', size = 1),
            legend.position = "none",
            panel.grid = element_blank(),
            plot.background = element_blank(),
            legend.background = element_blank(),
            axis.text = element_text(size = 11, color = 'black'),
            axis.line = element_blank(),
            plot.title = element_blank(),
            strip.text = element_text(color = "black", size = 13),
            strip.background = element_rect(fill = NA))
  }) -> ps1


ms2 <- read_csv('../../massspec.facs.csv') %>%
  mutate(cond = c('mESC_J1_par' = 'Par',
                  'mESC_J1-DNMT-TKO' = 'TKO')[type])

grep('WT|J1|TKO', colnames(cts), value = T) %>%
  grep('me[123]', ., value = T) %>%
  grep('J1_H3K27me3|J1_H3K36me2|27me[12]|TKO', ., value = T) %>%
  grep('36me3|rep1|_Par_', ., value = T, invert = T) %>%
  tibble(samp = .) %>%
  mutate(cond = ifelse(grepl('TKO', samp), 'TKO', 'Par'),
         mark = sub('_rep[12]$', '', samp) %>%
           sub('.*_', '', .)) %>%
  merge(lsz, by = 'samp', all.x = T) %>%
  mutate(rx = (endo.chip / endo.input) / (exo.chip / exo.input)) %>%
  merge(ms2, by = c('cond', 'mark'), all.X = T) %>%
  group_by(cond, mark, samp) %>%
  do(., tibble(v = cts[[.$samp]] * .$a * 100 / (.$endo.chip / length(bins))) %>%
       mutate(idx = 1:n())) %>%
  ungroup() %>%
  dplyr::select(cond, v, mark, idx) %>%
  pivot_wider(names_from = 'cond', values_from = 'v') %>%
  dplyr::select(x = Par, y = TKO, idx, mark) %>%
  mutate(mark = factor(mark, c(paste0('H3K27me', 3:1), 'H3K36me2'))) %>%
  split(., .$mark) %>%
  lapply(function(d) {
    pdat <- d %>%
      filter(x != 0 | y != 0) %>%
      mutate(across(c(x, y), ~log10(.x))) %>%
      filter(is.finite(x) & is.finite(y)) %>%
      filter(between(x, 0, 2) & between(y, 0, 2)) %>%
      mutate(r = olap[idx]) %>%
      dplyr::select(x, y, r) %>%
      rbind(tibble(x = c(0, 0, 2, 2),
                   y = c(0, 0, 2, 2),
                   r = rep(0, 4)))
    hex <- hexbin(pdat$x, pdat$y, xbins = 50, IDs = T)
    pdat$cell <- hex@cID
    hex <- data.frame(hcell2xy(hex),
                      cell = hex@cell,
                      count = hex@count)
    pdat2 <- pdat %>%
      group_by(cell) %>%
      summarise(prop = mean(r, na.rm = T)) %>%
      ungroup() %>%
      right_join(hex, by = "cell") %>%
      mutate(logcount = log10(count),
             mark = d$mark[1])
    
    pow <- .4
    ggplot() +
      geom_hex(aes(x = x, y = y, fill = prop, alpha = count^pow, color = prop),
               stat = "identity", color = NA, data = pdat2, size = 5) +
      scale_fill_gradientn(colors = gradientn1, name = "% Genic",
                           breaks = leg.brks, labels = leg.labs,
                           limits = c(-1, 1)) +
      scale_color_gradientn(colors = gradientn1, name = "% Genic",
                            breaks = leg.brks, labels = leg.labs,
                            limits = c(-1, 1)) +
      scale_alpha(range = c(0.01, 1), name = "Number of bins", guide = F) +
      coord_cartesian(expand = F) +
      scale_x_continuous('WT', breaks = 0:2,
                         labels = c(0, 10,  100)) +
      scale_y_continuous('DNMT-TKO', breaks = 0:2,
                         labels = c(0,  10, 100)) +
      annotate('segment', x = -Inf, xend =Inf, y = Inf,
               yend = Inf, color = 'white') +
      geom_abline(slope = 1, intercept = 0, color = 'grey50', alpha = .5) +
      facet_wrap(~mark) +
      theme(panel.background = element_rect(fill =NA, color = 'black', size = 1),
            legend.position = "none",
            panel.grid = element_blank(),
            plot.background = element_blank(),
            legend.background = element_blank(),
            axis.text = element_text(size = 11, color = 'black'),
            axis.line = element_blank(),
            plot.title = element_blank(),
            strip.text = element_text(color = "black", size = 13),
            strip.background = element_rect(fill = NA))
  })  -> ps2



{wrap_plots(ps1$H3K36me2, ps2$H3K36me2, nrow = 1) &
    theme(plot.background = element_blank())} %>%
  ggsave('f5_c.pdf', ., height = 2.7, width = 5.1, bg = 'transparent')

