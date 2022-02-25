library(data.table)
library(tidyverse)
library(matrixStats)
library(ggh4x)
library(ggnewscale)

list.files('../aggr', full.names = T, pattern = '27me1') %>%
  grep('Ezh|E14', ., value = T) %>%
  grep('Lost|Y726D', ., value = T, invert = T) %>%
  tibble(p = .) %>%
  mutate(f = sub('.scale.mat.gz', '', basename(p))) %>%
  separate(f, c('track','reg'), '\\.', F) %>%
  group_by(track, reg) %>%
  do(., {
    x <- .
    m <-  fread(x$p, skip = 1, drop = 1:6) %>% 
      as.matrix()
    m[!is.na(m) & m > quantile(m, .9999, na.rm = T)] <- NA
    m %>%
      {tibble(mu = colMeans(., na.rm = T),
              sd = colSds(., na.rm = T),
              num = colSums(is.finite(m)))} %>%
      mutate(idx = 1:n())
  }) %>%
  ungroup() -> aggr

clrs2 <- c('WT' = '#0C8140', 'EZH2-KO' = '#EC7723', 'EZH1/2-DKO' = '#CC6828')
clrs <- c('Intergenic\nregions' = '#FBB040',
          'Active\ngenes' = '#1C75BC',
          'Inactive\ngenes' = '#7F3F98')

bg <- '#7F3F98'
aggr %>%
  mutate(cond = case_when(grepl('Ezh1', track) ~ 'EZH1/2-DKO',
                          grepl('Ezh2', track) ~ 'EZH2-KO', 
                          T ~ 'WT') %>%
           factor(names(clrs2)),
         idx = case_when(reg == 'intergenic' ~ as.integer(idx * 2 - 1), T ~ idx),
         reg = ifelse(reg == '2iactive', 'Active\ngenes', 'Intergenic\nregions'),
         se = sd / sqrt(num)) %>%
  ggplot(aes(x = idx, y = mu, ymin = mu - se, ymax = mu + se)) + 
  geom_rect(aes(fill = reg), xmin = 40, xmax = 80, ymin = -Inf, ymax = Inf,
            data = ~distinct(., reg, cond), inherit.aes = F,
            alpha = .2, show.legend = F) +
  scale_fill_manual(values = clrs) +
  new_scale_fill() +
  geom_ribbon(aes(fill = cond, group = track), alpha = .5) +
  geom_line(aes(color = cond, group = track), size = 1) + 
  scale_color_manual(values = clrs2) +
  scale_fill_manual(values = clrs2) +
  facet_grid2(reg ~ cond, scales =  'free', independent = 'all') +
  facetted_pos_scales(x = list(
    reg == 'Intergenic\nregions' ~
      scale_x_continuous(breaks = c(1, 40, 80, 120),
                         labels = c('-20kb','Start','End','+20kb')),
    T ~
      scale_x_continuous(breaks = c(1, 40, 80, 120),
                         labels = c('-20kb','TSS','TES', '+20kb'))
  )) +
  scale_y_continuous('H3K27me1 (CPM)',
                     breaks = scales::breaks_pretty(3)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        #legend.direction = 'hoirizon',
        axis.title.x = element_blank(),
        strip.text = element_text(color = 'black', size = 13),
        panel.grid = element_blank()) -> p
p
ggsave('sf4_d.pdf', p, height = 4.2, width = 6.6, bg = 'transparent')

