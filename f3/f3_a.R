library(data.table)
library(tidyverse)
library(readxl)
library(ggh4x)
library(ggtext)
library(ggpubr)

signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
         symbols = c("****", "***", "**", "*", "ns"))
}

clrs2 <- c('WT' = '#0C8140', 'DNMT-TKO' = '#3953A4')

pdat <- read_xlsx('../data/histone_ratios.xlsx') %>%
  {.[,colSums(is.na(.)) != nrow(.)]} %>% 
  na.omit() %>% 
  {.[,c(1, which(.[1,] == 'Area'))]} %>%
  tail(-1) %>%
  dplyr::rename(p = 1) %>%
  `names<-`(., sub('\\..*', '', names(.))) %>%
  separate(p, c('pep', 'mod'), '\\ ') %>%
  mutate(pep = sub('H33', 'H3', pep)) %>%
  group_by(pep) %>%
  mutate(across(-mod, function(x) as.numeric(x) / sum(as.numeric(x)))) %>%
  ungroup() %>% 
  filter(grepl('^H[34]', pep) & mod != 'unmod') %>%
  mutate(mod = strsplit(mod, '(?<=.)(?=K)', perl = T)) %>% 
  unnest(mod) %>%
  mutate(mod = paste0(substr(pep, 1, 2), mod)) %>%
  dplyr::select(-pep) %>%
  group_by(mod) %>%
  summarise(across(everything(), sum), .groups = 'drop') %>%
  pivot_longer(-mod, names_to = 'samp', values_to = 'a') %>%
  filter(grepl('J1_par|TKO', samp) &
           grepl('H3K(27|36)me', mod)) %>%
  mutate(cond = ifelse(grepl('TKO', samp), 'DNMT-TKO', 'WT') %>%
           factor(names(clrs2)),
         a = a * 100) %>%
  separate(mod, c('resid', 'mod'), 'me') %>%
  arrange(cond) %>%
  mutate(mod = paste0('me', mod),
         x = sprintf("<span style='color:%s'>%s</span>", clrs2[cond], cond) %>%
           fct_inorder())
ann <- pdat %>%
  group_by(resid, mod) %>%
  do(., {
    ym <- max(.$a)
    compare_means(a ~ x, ., method = 't.test') %>%
      mutate(y = ym)
  }) %>%
  ungroup() %>%
  mutate(sig = as.character(signif.num(p)),
         y = y + 1) 

ggplot(pdat, aes(x = x, y = a, fill = cond, color = cond)) +
  stat_summary(geom = 'errorbar', fun.data = mean_se, width = .25) +
  stat_summary(geom = 'col', fun = mean) +
  stat_pvalue_manual(ann, label = 'sig', y.position = 'y', inherit.aes = F,
                     tip.length = 0) +
  facet_nested(.~resid + mod,
               nest_line = element_line(color = 'black')) +
  scale_fill_manual(values = clrs2) +
  scale_color_manual(values = clrs2) +
  scale_y_continuous('Abundance (%)', expand = expansion(c(0, .1))) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        legend.position = 'none',
        strip.background = element_rect(fill = NA),
        strip.text = element_text(color= 'black', size = 13),
        axis.line.x = element_line(color = 'black'),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(color = 'grey80', linetype = 'dashed'),
        axis.text.x = element_markdown(angle = 45, hjust = 1),
        axis.text =element_text(color = 'black', size = 11)) -> p
p
ggsave('f3_a.pdf', p, height = 2.7, width = 4.5, bg = 'transparent')
