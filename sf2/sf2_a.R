library(data.table)
library(tidyverse)
library(matrixStats)
library(ggh4x)
library(ggnewscale)
library(patchwork)

libsz <- read_csv('../data/libsz.csv') %>%
  dplyr::select(samp, endo.chip) %>%
  deframe()

ms <- read_csv('../../massspec.facs.csv') %>%
  filter(grepl('TKO|par|sgNSD1', type)) %>%
  mutate(cond = c('mESC_A1_par' = 'WT_A1',
                  'mESC_A1_sgNSD1' = 'sgNSD1',
                  'mESC_J1_par' = 'WT_J1',
                  'mESC_J1-DNMT-TKO' = 'DNMT-TKO')[type])

list.files('../aggr', full.names = T, pattern = 'activeParNsd1') %>%
  grep('27me', ., value = T) %>%
  grep('me[12]|rep2|J1', ., value = T) %>%
  grep('NSD1|_Par', ., value = T, invert = T) %>%
  tibble(p = .) %>%
  mutate(f = sub('.scale.mat.gz', '', basename(p))) %>%
  separate(f, c('track','reg'), '\\.', F) %>%
  mutate(cnd = case_when(grepl('Par|A1', track) ~ 'WT_A1',
                          grepl('TKO', track) ~ 'DNMT-TKO',
                          grepl('NSD', track) ~ 'sgNSD1',
                          T ~ 'WT_J1'),
         mark = sub('_rep[12]$', '', track) %>% sub('.*_', '', .)) %>% 
  group_by(track, reg, cnd, mark) %>%
  do(., {
    x <- .
    m <-  fread(x$p, skip = 1, drop = 1:6) %>% 
      as.matrix()
    m[!is.na(m) & m > quantile(m, .9999, na.rm = T)] <- NA
    if (x$track %in% names(libsz)) {
      m <- 1e6 * m * ms$f[ms$mark == x$mark & ms$cond == x$cnd] / libsz[x$track]
    }
    m %>%
      {tibble(mu = colMeans(., na.rm = T),
              sd = colSds(., na.rm = T),
              num = colSums(is.finite(m)))} %>%
      mutate(idx = 1:n())
  }) %>%
  ungroup() %>%
  mutate(cond = sub('_.*', '', cnd)) -> aggr

clrs2 <- c('WT' = '#0C8140','DNMT-TKO' = '#3953A4')
clrs <- c('Intergenic regions' = '#FBB040',
          'Active genes' = '#1C75BC',
          'Inactive genes' = '#7F3F98')

bg <- '#1C75BC'
aggr %>%
  mutate(cond = factor(cond, names(clrs2)),
         mark = factor(mark, paste0('H3K27me', 3:1)),
         se = sd / sqrt(num)) %>%
  ggplot(aes(x = idx, y = mu, ymin = mu - se, ymax = mu + se)) + 
  geom_rect(xmin = 40, xmax = 80, ymin = -Inf, ymax = Inf,
            data = ~distinct(., mark), inherit.aes = F,
            fill = bg, alpha = .2) +
  geom_ribbon(aes(fill = cond), alpha = .5) +
  geom_line(aes(color = cond), size = 1) + 
  scale_color_manual(values = clrs2) +
  scale_fill_manual(values = clrs2) +
  facet_nested(.~'Active genes' + mark, scales=  'free_x',
               nest_line = element_line(color = 'black')) +
  scale_x_continuous('Active genes', breaks = c(1, 40, 80, 120),
                     labels = c('-20kb','TSS','TES','+20kb')) +
  scale_y_continuous('MS norm coverage',
                     breaks = scales::breaks_pretty(3),
                     limits = c(.025, .105)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        legend.background = element_blank(),
        axis.title.x = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0,0),
        legend.justification = c(0,1.4),
        legend.direction = 'horizontal',
        strip.text = element_text(color = 'black', size = 13),
        panel.grid = element_blank()) -> p

ggsave('sf2_a.pdf', p, height = 2.2, width = 6, bg = 'transparent')

