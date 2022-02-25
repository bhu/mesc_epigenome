library(data.table)
library(tidyverse)
library(matrixStats)
library(ggh4x)

rx <- read_csv('../data/libsz.csv') %>%
  mutate(rx = (endo.chip / endo.input) / (exo.chip / exo.input),
         fac = 1e6 / endo.chip) %>%
  dplyr::select(samp, fac) %>%
  deframe()

rx <- read_csv('../data/ferrari.libsz.csv') %>%
  mutate(fac = 1e6 / endo.chip) %>%
  dplyr::select(samp, fac) %>%
  deframe() %>%
  c(rx)

resampler <- function(data) {
  n <- nrow(data)
  resample.rows <- sample(1:n,size=n,replace=TRUE)
  return(data[resample.rows,])
}

spline.estimator <- function(data,m=300, spar=.5) {
  fit <- smooth.spline(x=data[,1],y=data[,2],cv=TRUE, spar=spar)
  eval.grid <- seq(from=min(data[,1]),to=max(data[,1]),length.out=m)
  return(predict(fit,x=eval.grid)$y) # We only want the predicted values
}

spline.cis <- function(data,B,alpha=0.05,m=300) {
  spline.main <- spline.estimator(data,m=m)
  spline.boots <- replicate(B,spline.estimator(resampler(data),m=m))
  cis.lower <- 2*spline.main - apply(spline.boots,1,quantile,probs=1-alpha/2)
  cis.upper <- 2*spline.main - apply(spline.boots,1,quantile,probs=alpha/2)
  return(list(main.curve=spline.main,lower.ci=cis.lower,upper.ci=cis.upper,
              x=seq(from=min(data[,1]),to=max(data[,1]),length.out=m)))
}

list.files('../aggr', full.names = T, pattern = '27ac.K36me2Lost') %>%
  grep('E14|Ezh', ., value = T, invert = T) %>%
  tibble(p = .) %>%
  mutate(f = sub('.scale.mat.gz', '', basename(p))) %>%
  separate(f, c('track','reg'), '\\.', F) %>%
  group_by(track, reg) %>%
  do(., {
    x <- .
    m <-  fread(x$p, skip = 1, drop = 1:6) %>% 
      as.matrix()
    #m[!is.na(m) & m > quantile(m, .9999, na.rm = T)] <- NA
    if (x$track %in% names(rx)) {
      m <- m * rx[x$track]
    }
    #m %>%
    #  {tibble(mu = colMeans(., na.rm = T),
    #          sd = colSds(., na.rm = T),
    #          num = colSums(is.finite(m)))} %>%
    #  mutate(idx = 1:n())
    as.data.frame(m) %>%
      mutate(reg = 1:n()) %>%
      pivot_longer(-reg, names_to = 'pos') %>%
      mutate(pos = sub('^V', '', pos),
             pos = as.integer(pos) - 6) %>%
      filter(is.finite(value)) %>%
      dplyr::select(x = pos, y = value) %>%
      as.data.frame() -> data
    sp.cis <- spline.cis(data, B = 100, alpha = 0.05, m = 120)
    data.frame(sp.cis)
  }) %>%
  ungroup() -> aggr

clrs2 <- c('WT' = '#0C8140', 'EED-KO' = '#F26522', 'sgNSD1' = '#EE2E2E')
clrs <- c('Intergenic\nregions' = '#FBB040',
          'Active\ngenes' = '#1C75BC',
          'Inactive\ngenes' = '#7F3F98')

bg <- '#7F3F98'
aggr %>%
  mutate(cond = case_when(grepl('EedKO', track) ~ 'EED-KO', 
                          grepl('NSD1', track) ~ 'sgNSD1',
                          T ~ 'WT'),
         mark = 'H3K27ac',
         study = ifelse(grepl('E36', track), 'Ferrari 2014', 'This study') %>%
           factor(c('This study', 'Ferrari 2014'))) %>%
  ggplot(aes(x = x, y = main.curve, ymin = lower.ci, ymax = upper.ci)) + 
  geom_rect(xmin = 40, xmax = 80, ymin = -Inf, ymax = Inf, inherit.aes = F,
            fill = bg, alpha = .2, data = ~distinct(., mark, study)) +
  geom_ribbon(aes(fill = cond, group = track), alpha = .5) +
  geom_line(aes(color = cond, group = track), size = 1) + 
  scale_color_manual(values = clrs2) +
  scale_fill_manual(values = clrs2) +
  facet_grid2(study ~ mark, scales =  'free') +
  scale_x_continuous('Regions depleted of\nH3K36me2 in sgNSD1',
                     breaks = c(1, 40, 80, 120),
                     labels = c('-20kb','Start','End','+20kb'),
                     expand = expansion(0)) +
  scale_y_continuous('Depth norm coverage',
                     breaks = scales::breaks_pretty(3)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        strip.text = element_text(color = 'black', size = 13),
        panel.grid = element_blank()) -> p

p
ggsave('f6_b.pdf', p, height = 5.5, width = 4.4, bg = 'transparent')

