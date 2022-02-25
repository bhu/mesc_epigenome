library(tidyverse)
library(ComplexHeatmap)
library(pals)

cts <- readRDS('../data/cts.rds')$`10kb`
cts <- log10(cts + 1)

md <- tibble(samp = colnames(cts)) %>%
  mutate(cond = case_when(grepl('TKO', samp) ~ 'DNMT-TKO',
                          grepl('EZH2', samp) ~ 'EZH2-KO',
                          grepl('NSD1', samp) ~ 'sgNSD1',
                          T ~ 'WT'),
         mark = sub('_rep[0-9]$', '', samp) %>%
           sub('.*_', '', .) %>%
           sub('HA', 'input', .) %>%
           sub('^K36me$', 'input', .),
         rep = case_when(grepl('sgNSD1|rep1|Par', samp) ~ '1',
                         cond == 'WT' & grepl('A1', samp) ~ '3',
                         T ~ '2'),
         nm = sprintf('%s %s %s', cond, mark, rep))

ss <- md %>%
  filter(mark != 'input') %>%
  add_count(cond, mark) %>%
  filter(n > 1)

cm <- cts[,ss$samp] %>%
  `colnames<-`(ss$nm) %>%
  cor(method = 'spearman') -> cm

hc <- hclust(as.dist(1 - cm), method = 'average')
plot(hc)

ann <- md %>%
  filter(samp %in% ss$samp) %>%
  dplyr::select(Condition = cond, Mark = mark, nm) %>%
  column_to_rownames('nm') %>%
  .[ss$nm,] 

aclrs <- list(Condition = c('WT' = '#0C8140','DNMT-TKO' = '#3953A4', 'sgNSD1' = '#EE2E2E'),
              Mark = setNames(tableau20()[c(3,9,11,15)], sort(unique(ann$Mark))))

pdf('sf1_d.pdf', height = 4, width = 6.8)
pheatmap(cm, cluster_rows = hc, cluster_col = hc, show_row_dend  = F,
         show_colnames = F, border_color = NA, name = 'Correlation',
         col = rev(brewer.rdylbu(11)), breaks = seq(0, 1, .1),
         annotation_col = ann, annotation_colors = aclrs)
dev.off()

