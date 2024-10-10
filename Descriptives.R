#####################################################################################
### Analyse and plot data on first positive test results in Ho Chi Minh City (incidence.csv)
#####################################################################################

library(tidyverse)
library(ggplot2)
library(scales)
library(incidence)
library(RColorBrewer)
library(ggpattern)

## Read data
linelist <- read.csv("Data/Incidence.csv") #
linelist$first_date_test_positive <- as.Date(linelist$first_date_test_positive)

## compute incidence of new positive test results
inc_obj <- incidence(linelist$first_date_test_positive, 1, groups = linelist$set)
inc_obj$inclusion <- ifelse(inc_obj$dates < inc_peak, "exponential growth", "delay in reporting")
inc_peak <- find_peak(inc_obj)


inc_obj_new <- with(linelist, incidence(first_date_test_positive[first_date_test_positive < inc_peak],1))
inc_fit <- fit(inc_obj_new)
r <- inc_fit$info$r
r025 <- inc_fit$info$r.conf[1]
r975 <- inc_fit$info$r.conf[2]
est_se <- (r975 - r025)/(2*1.96)

cols <- brewer_pal(palette = "Greys")(9)[c(4,7)]
p_epicurve <-
    ggplot() +
    geom_col(aes( x = inc_obj$dates, y = rowSums(inc_obj$counts), fill = inc_obj$inclusion), col = "black") +
    geom_col_pattern(aes( x = inc_obj$dates, y = inc_obj$counts[,1]), fill = NA, col = "black",
                   ## For the pattern:
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.005,
                   pattern_spacing = 0.005) +
    geom_line(aes(x = inc_fit$info$pred$dates, y = inc_fit$info$pred$fit,
                  linetype = "doubling time (6.55 days)")) +
    geom_line(aes(x = inc_fit$info$pred$dates, y = inc_fit$info$pred$lwr,
                  linetype = "confidence limits (95%, 6.05-7.14 days)"), ) +
    geom_line(aes(x = inc_fit$info$pred$dates, y = inc_fit$info$pred$upr,
                  linetype = "confidence limits (95%, 6.05-7.14 days)")) +
    theme_bw() + coord_cartesian(ylim = c(0,2500)) +
    labs(x = "First positive test day (2021)", y = "Individuals (n)", title = "",
         fill = "Epidemic phase", linetype = "Estimated growth") +
    scale_fill_manual(values = cols[c(2,1)]) +
    theme(legend.position = "inside", legend.position.inside = c(0.3, 0.55), panel.border = element_rect(linewidth = 0.5, linetype = "solid", colour = "black"))




