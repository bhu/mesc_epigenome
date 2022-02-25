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

bws <- list(Parental = 'mESC_Par_mCG.bw',
            sgNSD1 = 'mESC_sgNSD1_mCG.bw') %>%
  lapply(function(x) {
    import.bw(sprintf('../tracks/%s', x))
  })

scs <- lapply(bws, function(x) {
  findOverlaps(bins, x) %>%
    as("List") %>%
    extractList(x$score, .) %>%
    mean(na.rm = T)
}) %>% bind_cols()

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
ggsave('f1_d.leg.pdf', leg, device = cairo_pdf, height = 1.2, width = 1.5)


lsz <- read_csv('../data/libsz.csv') 

ms <- read_csv('../../massspec.facs.csv') %>%
  filter(grepl('A1_par|sgNSD1', type)) %>%
  mutate(cond = c('mESC_A1_par' = 'Par',
                  'mESC_A1_sgNSD1' = 'sgNSD1')[type])


grep('Par|sgNSD1', colnames(cts), value = T) %>%
  grep('me[123]|mCG', ., value = T) %>%
  tibble(samp = .) %>%
  separate(samp, c(NA, 'cond', 'mark'), '_', F) %>%
  merge(lsz, by = 'samp', all.x = T) %>%
  mutate(rx = (endo.chip / endo.input) / (exo.chip / exo.input)) %>%
  merge(ms, by = c('cond', 'mark'), all.X = T) %>%
  group_by(cond, mark, samp) %>%
  do(., tibble(v = cts[[.$samp]] * .$a * 100 / (.$endo.chip / length(bins))) %>%
       mutate(idx = 1:n())) %>%
  ungroup() %>%
  dplyr::select(cond, v, mark, idx) %>%
  pivot_wider(names_from = 'cond', values_from = 'v') %>%
  dplyr::select(x = Par, y = sgNSD1, idx, mark) %>%
  rbind(mutate(dplyr::rename(scs, x = Parental, y = sgNSD1) , idx = 1:n(), mark = 'DNAme')) %>%
  mutate(mark = factor(mark, c(paste0('H3K27me', 3:1),
                               'H3K36me3', 'H3K36me2', 'DNAme'))) %>%
  split(., .$mark) %>%
  lapply(function(d) {
    pdat <- if (d$mark[1] != 'DNAme') {
      d %>%
        filter(x != 0 | y != 0) %>%
        mutate(across(c(x, y), ~log10(.x))) %>%
        filter(is.finite(x) & is.finite(y)) %>%
        filter(between(x, 0, 2) & between(y, 0, 2)) %>%
        mutate(r = olap[idx]) %>%
        dplyr::select(x, y, r) %>%
        rbind(tibble(x = c(0, 0, 2, 2),
                     y = c(0, 0, 2, 2),
                     r = rep(0, 4)))
    } else {
      d %>%
        filter(is.finite(x) & is.finite(y)) %>%
        mutate(r = olap[idx]) %>%
        dplyr::select(x, y, r) %>%
        rbind(tibble(x = c(0, 0, 100, 100),
                     y = c(0, 0, 100, 100),
                     r = rep(0, 4)))
    }
    
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
    
    p <- ggplot() +
      geom_hex(aes(x = x, y = y, fill = prop, alpha = count^pow, color = prop),
               stat = "identity", color = NA, data = pdat2, size = 5) +
      scale_fill_gradientn(colors = gradientn1, name = "% Genic",
                           breaks = leg.brks, labels = leg.labs,
                           limits = c(-1, 1)) +
      scale_color_gradientn(colors = gradientn1, name = "% Genic",
                            breaks = leg.brks, labels = leg.labs,
                            limits = c(-1, 1)) +
      scale_alpha(range = c(0.01, 1), name = "Number of bins", guide = F) +
      geom_abline(slope = 1, intercept = 0, color = 'grey50', alpha = .5) +
      facet_wrap(~mark) +
      theme(panel.background = element_rect(fill = NA, color = 'black', size = 1),
            legend.position = "none",
            panel.grid = element_blank(),
            plot.background = element_blank(),
            legend.background = element_blank(),
            axis.text = element_text(size = 11, color = 'black'),
            axis.line = element_blank(),
            plot.title = element_blank(),
            strip.text = element_text(color = "black", size = 13),
            strip.background = element_rect(fill = NA))
    
    p <- if (d$mark[1] != 'DNAme') {
      p + 
        coord_cartesian(expand = F, xlim = c(0,2), ylim = c(0, 2)) +
        scale_x_continuous('WT', breaks = 0:2,
                           labels = c(0, 10,  100)) +
        scale_y_continuous('sgNSD1', breaks = 0:2,
                           labels = c(0,  10, 100))
    } else {
      p +
        coord_cartesian(expand = F, xlim = c(0, 100), ylim = c(0, 100)) +
        scale_x_continuous('WT', breaks = c(0, 50, 100),
                           labels = c(0, 50, 100)) +
        scale_y_continuous('sgNSD1', breaks = c(0, 50, 100),
                           labels = c(0, 50, 100)) 
    }
    
  }) %>%
  {wrap_plots(., nrow = 2) & theme(plot.background = element_blank())} %>%
  ggsave('f1_d.pdf', ., height = 5, width = 6.7, bg = 'transparent')

