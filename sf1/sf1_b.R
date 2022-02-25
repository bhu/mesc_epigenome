library(data.table)
library(tidyverse)
library(matrixStats)
library(ggh4x)
library(patchwork)
library(ggnewscale)

libsz <- read_csv('../data/libsz.csv') %>%
  dplyr::select(samp, endo.chip) %>%
  deframe()

ms <- read_csv('../../massspec.facs.csv') %>%
  filter(grepl('A1_par|sgNSD1', type)) %>%
  mutate(cond = c('mESC_A1_par' = 'WT',
                  'mESC_A1_sgNSD1' = 'sgNSD1')[type])

list.files('../aggr', full.names = T, pattern = '36me2|27me3')  %>%
  grep('_Par_|sgNSD1|J1|EZH2|DNMT|Eed', ., invert = T, value = T) %>%
  grep('ParNsd|intergenic', ., value = T) %>%
  grep('scale', ., value = T) %>%
  tibble(p = .) %>%
  mutate(f = sub('.mat.gz', '', basename(p))) %>%
  separate(f, c('track','reg', NA), '\\.', F) %>%
  separate(track, c(NA, 'cond', NA, 'mark'), '_', F) %>%
  mutate(cond = sub('NSD1', 'sgNSD1', cond)) %>%
  group_by(track, reg, cond, mark) %>%
  do(., {
    x <- .
    m <-  fread(x$p, skip = 1, drop = 1:6) %>% 
      as.matrix()
    m[!is.na(m) & m > quantile(m, .9999, na.rm = T)] <- NA
    if (x$track %in% names(libsz)) {
      m <- 1e6 * m * ms$f[ms$cond == x$cond & ms$mark == x$mark] / libsz[x$track]
    }
    m %>%
      {tibble(mu = colMeans(., na.rm = T),
              sd = colSds(., na.rm = T),
              num = colSums(is.finite(m)))} %>%
      mutate(idx = 1:n())
  }) %>%
  ungroup() -> aggr

clrs <- c('Intergenic\nregions' = '#FBB040',
          'Active\ngenes' = '#1C75BC',
          'Inactive\ngenes' = '#7F3F98')

clrs2 <- c('sgNSD1' = '#EE2E2E', 'WT' = '#0C8140')


aggr %>%
  mutate(cond = ifelse(grepl('NSD', track), 'sgNSD1', 'WT'),
         mark = sub('.*_', '', track) %>%
           factor(c('H3K36me3','H3K36me2','DNAme',
                    'H3K27me3', 'H3K27me2','H3K27me1')),
         rr = c('activeParNsd1' = 'Active\ngenes',
                'intergenic' = 'Intergenic\nregions',
                'silentParNsd1' = 'Inactive\ngenes')[reg],
         se = sd / sqrt(num),
         idx = case_when(reg == 'intergenic' ~ as.integer(idx * 2 - 1), T ~ idx)) %>%
  ggplot(aes(x = idx, y = mu, ymin = mu - se, ymax = mu + se)) + 
  geom_rect(aes(fill = rr), xmin = 40, xmax = 80, ymin = -Inf, ymax = Inf,
            data = ~distinct(., mark, rr), inherit.aes = F,
            alpha = .2, show.legend = F) +
  scale_fill_manual(values = clrs) +
  new_scale_fill() +
  geom_ribbon(aes(fill = cond), color = NA, alpha = .5) +
  geom_line(aes(color = cond), size = 1) + 
  scale_fill_manual(values = clrs2) +
  scale_color_manual(values = clrs2) +
  facet_nested(mark~ 'Replicate 2' + rr, scales=  'free',
               nest_line = element_line(color = 'black')) +
  facetted_pos_scales(x = list(
    rr == 'Intergenic\nregions' ~
      scale_x_continuous(breaks = c(1, 40, 80, 120),
                         labels = c('-20kb','Start','End','+20kb')),
    T ~
      scale_x_continuous(breaks = c(1, 40, 80, 120),
                         labels = c('-20kb','TSS','TES', '+20kb'))
  )) +
  scale_y_continuous('MS norm coverage') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_blank(),
        legend.title = element_blank(),
        strip.text = element_text(color = 'black', size = 13),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = 'bottom',
        legend.justification = 'left',
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 11),
        axis.title.x = element_blank()) -> p2
p2
ggsave('sf1_b.pdf', p2, height = 4.0, width = 6, bg = 'transparent')


