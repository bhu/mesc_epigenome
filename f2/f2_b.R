library(data.table)
library(tidyverse)
library(matrixStats)
library(ggh4x)
library(ggnewscale)
library(patchwork)

rx <- read_csv('../data/libsz.csv') %>%
  mutate(fac = 1e6 / endo.chip) %>%
  dplyr::select(samp, fac) %>%
  deframe()

list.files('../aggr', full.names = T, pattern = 'intergenic|activeParNsd1') %>%
  grep('DNMT|mCHH', ., value = T) %>%
  grep('_Par|_sgNSD1', ., value = T) %>%
  tibble(p = .) %>%
  mutate(f = sub('.scale.mat.gz', '', basename(p))) %>%
  separate(f, c('track','reg'), '\\.', F) %>%
  group_by(track, reg) %>%
  do(., {
    x <- .
    m <-  fread(x$p, skip = 1, drop = 1:6) %>% 
      as.matrix()
    m[!is.na(m) & m > quantile(m, .9999, na.rm = T)] <- NA
    if (x$track %in% names(rx)) {
      m <- m * rx[x$track]
    }
    m %>%
      {tibble(mu = colMeans(., na.rm = T),
              sd = colSds(., na.rm = T),
              num = colSums(is.finite(m)))} %>%
      mutate(idx = 1:n())
  }) %>%
  ungroup() -> aggr

clrs2 <- c('WT' = '#0C8140', 'sgNSD1' = '#EE2E2E')
clrs <- c('Intergenic regions' = '#FBB040',
          'Active genes' = '#1C75BC',
          'Inactive genes' = '#7F3F98')

bg <- '#1C75BC'
aggr %>%
  mutate(cond = ifelse(grepl('NSD1', track), 'sgNSD1', 'WT') %>%
           factor(names(clrs2)),
         track = sub('mESC_[A-Za-z0-9]*_', '', track) %>%
           sub('mCHH', 'DNAme (mCHH)', .) %>%
           sub('DNMT3A2', 'DNMT3A', .) %>%
           factor(c('DNMT3A','DNAme (mCHH)')),
         reg = c('activeParNsd1' = 'Active genes',
                 'intergenic' = 'Intergenic regions')[reg],
         se = sd / sqrt(num),
         xmin = ifelse(reg == 'Active genes', 40, 20),
         xmax = ifelse(reg == 'Active genes', 80, 40)) %>%
  split(., .$track) %>%
  lapply(function(d) {
    if (d$reg[1] == 'Intergenic regions') {
      k <- 60
      l1 <- 'Start'
      l2 <- 'End'
      bg <- '#FBB040'
    } else {
      k <- 120
      l1 <- 'TSS'
      l2 <- 'TES'
      bg <- '#1C75BC'
    }
    if (d$track[1] == 'DNMT3A') {
      yl <- 'Depth norm coverage'
    } else {
      yl <- 'mCHH/CHH (%)'
    }
    p <- ggplot(d, aes(x = idx, y = mu, ymin = mu - se, ymax = mu + se)) + 
      geom_rect(aes(fill = reg, xmin = xmin, xmax = xmax), alpha = .2,
                ymin = -Inf, ymax = Inf, inherit.aes = F, show.legend = F,
                data = ~distinct(., track, reg, .keep_all = T)) +
      scale_fill_manual(values = clrs) +
      new_scale_fill() +
      geom_ribbon(aes(fill = cond), alpha = .5) +
      geom_line(aes(color = cond), size = 1) + 
      scale_color_manual(values = clrs2) +
      scale_fill_manual(values = clrs2) +
      facet_nested(. ~ track + reg, scales = 'free', independent = 'y',
                   nest_line = element_line(color = 'black')) +
      facetted_pos_scales(x = list(
        scale_x_continuous(breaks = c(1, 40, 80, 120),
                           labels = c('-20kb','TSS','TES','+20kb')),
        scale_x_continuous(breaks = c(1, 20, 40, 60),
                           labels = c('-20kb','Start','End','+20kb'))
        
      )) +
      scale_y_continuous(yl, breaks = scales::breaks_pretty(3)) +
      theme(plot.background = element_blank(),
            panel.background = element_rect(fill = NA, color = 'black', size = 1),
            axis.text = element_text(color = 'black', size = 11),
            strip.background = element_blank(),
            legend.background = element_blank(),
            legend.key = element_blank(),
            legend.title = element_blank(),
            legend.position = c(0,0),
            legend.justification = c(0,1.4),
            legend.direction = 'horizontal',
            strip.text = element_text(color = 'black', size = 13),
            axis.title.x = element_blank(),
            panel.grid = element_blank())
    if (d$track[1] == 'DNMT3A') {
      p <- p + theme(legend.position = 'none',
                     strip.background.y = element_blank(),
                     strip.text.y = element_blank())
    } else {
      p <- p + theme(legend.position = c(1,0),
                     legend.justification = c(1,0),
                     legend.direction = 'vertical')
    }
    p
  }) %>%
  wrap_plots(ncol = 2) -> p

ggsave('f2_b.pdf', p & theme(plot.background = element_blank()),
       height = 2.2, width = 10, bg = 'transparent')

