library(data.table)
library(tidyverse)
library(matrixStats)
library(pals)

fread('../tracks/mESC_Parental_BS_1-bismethyl-mm10_CGI.bedGraph') %>% 
  filter(V4 != '.') %>% 
  mutate(beta = as.numeric(V4) / (as.numeric(V4) + as.numeric(V5)),
         q = cut(beta, c(0,.01,.2,.4,.8,1), include.lowest = T),
         crd = sprintf('%s:%d-%d', V1, V2, V3)) -> dd

div <- dd %>%
  dplyr::select(crd, q) %>%
  distinct() %>%
  deframe()

rr <- dd %>%
  group_by(q) %>%
  summarise(r = q) %>%
  #summarise(r = sprintf('%.3f - %.3f', min(beta), max(beta))) %>%
  ungroup()

fread('../data/massspec.facs.csv') %>%
  filter(mark == 'H3K27me3') %>%
  mutate(cond = case_when(grepl('TKO', type) ~ 'DNMT-TKO',
                          grepl('NSD', type) ~ 'sgNSD1',
                          grepl('par', type) ~ 'WT',
                          grepl('EZH2', type) ~ 'EZH2-KO'),
         line = ifelse(grepl('J1', type), 'J1', 'A1')) -> ms

fread('../data/libsz.csv') %>%
  filter(grepl('27me3', samp) & !grepl('EZH2', samp)) %>%
  mutate(cond = case_when(grepl('TKO', samp) ~ 'DNMT-TKO',
                          grepl('NSD', samp) ~ 'sgNSD1',
                          T ~ 'WT'),
         line = case_when(grepl('Par|NSD|A1', samp) ~ 'A1',
                          T ~ 'J1')) %>%
  merge(ms, by = c('cond', 'line')) %>%
  split(., .$samp) %>%
  lapply(function(x) {
    
    crd <- fread(sprintf('../aggr/%s.cgi.point.mat.gz', x$samp), 
                 skip = 1, select = 1:6)
    
    mc <- fread(sprintf('../aggr/%s.cgi.point.mat.gz', x$samp), 
                skip = 1, drop = 1:6)
    mi <- fread(sprintf('../aggr/%s.cgi.point.mat.gz', x$inp), 
                skip = 1, drop = 1:6)
    
    rx <-  (x$endo.chip / x$exo.chip) / (x$endo.input / x$exo.input)
    fc <- min(x$endo.chip, x$endo.input) / x$endo.chip
    fi <- min(x$endo.chip, x$endo.input) / x$endo.input
    
    crd %>% 
      mutate(idx = 1:n(), 
             q = div[V4]) %>% 
      na.omit() %>%
      split(., .$q) %>%
      lapply(function(y) {
         log2((mc[y$idx,] * fc * x$f + 1) / (mi[y$idx,] * fi + 1)) %>%
          as.matrix() %>%
          {tibble(mu = colMeans(., na.rm = T),
                  std = colSds(., na.rm = T),
                  num = colSums(is.finite(.)))} %>%
          mutate(idx = 1:n())
      }) %>%
      bind_rows(.id = 'q') %>%
      mutate(cond = x$cond[1], 
             line = x$line[1])
  }) %>%
  bind_rows(.id = 'samp') %>%
  mutate(se = std / sqrt(num)) -> aggr

clrs2 <- c('WT' = '#0C8140','DNMT-TKO' = '#3953A4', 'sgNSD1' = '#EE2E2E')
clrs <- tableau20()[c(1,9,5,3,7)]
aggr %>%
  #mutate(q = as.integer(q)) %>%
  #merge(rr, by = 'q') %>%
  #mutate(qnt = sprintf('%d (%s)', q, r) %>% fct_inorder()) %>%
  #mutate(qnt = q) %>%
  mutate(qnt = factor(q, levels(dd$q)),
         rep = ifelse(grepl('Par|sg|rep1', samp), 'Replicate 1', 'Replicate 2'),
         cond = case_when(grepl('TKO', samp) ~ 'DNMT-TKO', 
                          grepl('NSD', samp) ~ 'sgNSD1',
                          grepl('Par|A1', samp) ~ 'WT A1',
                          T ~ 'WT J1') %>%
           factor(c('WT A1', 'WT J1', 'DNMT-TKO', 'sgNSD1'))) %>%
  split(., .$rep) %>%
  lapply(function(x) {
    ggplot(x, aes(x = idx, y = mu, ymin = mu - se, ymax = mu + se)) +
      geom_ribbon(aes(fill = qnt), alpha = .5) +
      geom_line(aes(color = qnt), size = 1) + 
      facet_nested(. ~ rep + cond,
                   nest_line = element_line(color = 'black'),
                   strip = strip_nested(
                     text_x = elem_list_text(
                       color = c('black', clrs2[sub(' .*', '', (sort(unique(x$cond))))])
                     )
                   )) +
      scale_color_manual('Quantile of CGI mCG/CG',
                         values = clrs) +
      scale_fill_manual('Quantile of CGI mCG/CG',
                         values = clrs) +
      ylab('H3K27me3\nMS norm enrichment') +
      scale_x_continuous(breaks = c(1, 40, 80),
                         labels = c('-20kb', 'CGI', '+20kb')) +
      theme(plot.background = element_blank(),
            panel.background = element_rect(fill = NA, color = 'black', size = 1),
            panel.grid = element_blank(),
            strip.background = element_rect(fill = NA),
            strip.text = element_text(size = 13, color = 'black'),
            axis.text = element_text(color = 'black', size = 11),
            axis.title.x = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(size = 11),
            legend.background = element_blank(),
            panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
            axis.text.x = element_text(angle = 45, hjust = 1))
  }) -> ps


ggsave('sf4_d.top.pdf', ps$`Replicate 1`, height = 4.3, width = 9, bg = 'transparent')

ggsave('sf4_d.bot.pdf', ps$`Replicate 2`, height = 4.3, width = 12.5, bg = 'transparent')


