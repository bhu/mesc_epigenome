library(data.table)
library(tidyverse)
library(matrixStats)
library(ggh4x)

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
  grep('Ezh|E14', ., value = T) %>%
  tibble(p = .) %>%
  mutate(f = sub('.scale.mat.gz', '', basename(p))) %>%
  separate(f, c('track','reg'), '\\.', F) %>%
  group_by(track, reg) %>%
  do(., {
    x <- .
    m <-  fread(x$p, skip = 1, drop = 1:6) %>% 
      as.matrix()
    # m[!is.na(m) & m > quantile(m, .9999, na.rm = T)] <- NA
    # m %>%
    #   {tibble(mu = colMeans(., na.rm = T),
    #           sd = colSds(., na.rm = T),
    #           num = colSums(is.finite(m)))} %>%
    #   mutate(idx = 1:n()) %>%
    #   arrange(-mu)
    #   {smooth.spline(.$idx, .$mu)} %>% plot()
    # 
    # as.data.frame(m) %>%
    #   mutate(reg = 1:n()) %>%
    #   pivot_longer(-reg, names_to = 'pos') %>%
    #   mutate(pos = sub('^V', '', pos),
    #          pos = as.integer(pos) - 6) %>%
    #   filter(is.finite(value)) %>%
    #   gam(value ~ s(pos, k = 30), data = ., method = 'REML') -> modd
    
    as.data.frame(m) %>%
      mutate(reg = 1:n()) %>%
      pivot_longer(-reg, names_to = 'pos') %>%
      mutate(pos = sub('^V', '', pos),
             pos = as.integer(pos) - 6) %>%
      filter(is.finite(value)) %>%
      dplyr::select(x = pos, y = value) %>%
      as.data.frame() -> data
    
    #smooth.spline(data[,1], data[,2], spar = .5) %>%
    #  predict() %>%
    #  {tibble(x = .$x, y = .$y)}
    
    sp.cis <- spline.cis(data, B = 100, alpha = 0.05, m = 120)
    data.frame(sp.cis)
  }) %>%
  ungroup() -> aggr

clrs2 <- c('WT' = '#0C8140', 'EZH2-KO' = '#EC7723', 'EZH1/2-DKO' = '#CC6828')
bg <- '#7F3F98'
aggr %>%
  mutate(cond = case_when(grepl('Ezh1', track) ~ 'EZH1/2-DKO',
                          grepl('Ezh2', track) ~ 'EZH2-KO', 
                          T ~ 'WT') %>%
           factor(names(clrs2)),
         mark = 'H3K27ac') %>%
  ggplot(aes(x = x, y = main.curve, ymin = lower.ci, ymax = upper.ci)) + 
  geom_rect(xmin = 40, xmax = 80, ymin = -Inf, ymax = Inf, inherit.aes = F,
            fill = bg, alpha = .2, data = ~distinct(., mark)) +
  geom_ribbon(aes(fill = cond, group = track), alpha = .5) +
  geom_line(aes(color = cond, group = track), size = 1) + 
  scale_color_manual(values = clrs2) +
  scale_fill_manual(values = clrs2) +
  facet_grid2(. ~ mark, scales =  'free') +
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
        #legend.direction = 'hoirizon',
        strip.text = element_text(color = 'black', size = 13),
        panel.grid = element_blank()) -> p
p
ggsave('sf4_b.pdf', p, height = 2.8, width = 6.2, bg = 'transparent')

