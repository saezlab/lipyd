#!/usr/bin/env Rscript

# Dénes Türei 2017 EMBL
# turei.denes@gmail.com

library(ggplot2)
library(gtable)
library(grid)
library(RColorBrewer)
library(colorspace)
library(dplyr)

source('theme_black.R')

onlyclass1 <- FALSE
screen1file <- 'antonella_final.csv'
screen2file <- 'enric_processed.csv'
barpdf      <- 'intensities_bar.pdf'
ccunspdf    <- 'cc_uns_heatmap.pdf'

a <- read.table(screen1file, sep = '\t', header = TRUE)
aI <- a %>% filter(cls == 'I' & uhgroup != 'P40')

e <- read.table(screen2file, sep = '\t', header = TRUE)
eI <- e %>% filter(cls == 'I' & uhgroup != 'P40')

aeI <- rbind(aI, eI)

aeIhg <- aeI %>%
    group_by(protein, ionm, screen, uhgroup) %>%
    mutate(ihg = sum(as.numeric(intensity))) %>%
    group_by(protein, ionm, screen) %>%
    mutate(itotal = sum(as.numeric(ihg))) %>%
    group_by(protein, ionm, screen, uhgroup) %>%
    mutate(irel = ihg / itotal, scr_ionm = paste0(screen, ionm)) %>%
    summarize_all(first) %>%
    mutate(uhgroup = as.character(uhgroup))

p <- ggplot(aeIhg, aes(x = uhgroup, y = irel * 100.0 + 1.0)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    scale_y_log10() +
    facet_grid(
        protein ~ scr_ionm,
        labeller = labeller(
            scr_ionm = c(
                Apos = 'HEK cells / positive',
                Aneg = 'HEK cells / negative',
                Epos = 'E. coli / positive',
                Eneg = 'E. coli / negative'
            )
        )
    ) +
    xlab('Screening and ion mode') +
    ylab('Protein') +
    theme_bw() +
    theme(
        text = element_text(family = 'DINPro'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7, hjust = 1)
    )

cairo_pdf(barpdf, width = 18, height = 99)
    g <- ggplotGrob(p)
dev.off()

panels <- grep('panel', g$layout$name)
top <- unique(g$layout$t[panels])
ytitle <- grep('ylab-l', g$layout$name)

all <- gtable:::rbind_gtable(g[seq.int(min(top)),],
                            g[max(top)+1,], 'first')

for(i in seq(2, length(top) - 1)){
    all <- gtable:::rbind_gtable(
        all,
        gtable:::rbind_gtable(
            g[c(top[i]-1, top[i]),],
            g[max(top)+1,], 'first'
        ),
        'first'
    )
}

all <- gtable:::rbind_gtable(all, g[(max(top)-1):nrow(g),], 'first')
ally <- gtable:::gtable_add_grob(all, g$grobs[ytitle][[1]],
                        t = g$layout[ytitle,]$t,
                        b = g$layout[ytitle,]$b + length(top) - 1,
                        l = g$layout[ytitle,]$l,
                        r = g$layout[ytitle,]$r)

cairo_pdf(barpdf, width = 18, height = 99)

    grid.newpage()
    grid.draw(ally)

dev.off()

### 


aeI <- rbind(aI, eI)

aeIcc <- rbind(
        aeI %>% select(protein, ionm, screen, uhgroup, cc = fa1c, uns = fa1u, intensity),
        aeI %>% select(protein, ionm, screen, uhgroup, cc = fa2c, uns = fa2u, intensity),
        aeI %>% select(protein, ionm, screen, uhgroup, cc = fa3c, uns = fa3u, intensity)
    ) %>%
    filter(!is.nan(cc) & !is.nan(uns)) %>%
    group_by(protein, ionm, screen) %>%
    mutate(itotal = sum(as.numeric(intensity))) %>%
    group_by(protein, ionm, screen, cc, uns) %>%
    mutate(iccuns = sum(as.numeric(intensity))) %>%
    mutate(irel   = iccuns / itotal) %>%
    group_by(ionm, screen, cc, uns) %>%
    mutate(irels  = sum(irel)) %>%
    summarize_all(first)

p <- ggplot(aeIcc, aes(x = uns, y = cc, fill = log10(irels))) +
    geom_tile() +
    scale_fill_continuous(name = 'Relative\nintensity') +
    facet_grid(screen ~ ionm,
               labeller = labeller(screen = c(A = 'HEK cells', E = 'E. coli'),
                                   ionm = c(neg = 'negative', pos = 'positive'))) +
    theme_bw() +
    xlab('Unsaturation') +
    ylab('Carbon count') +
    theme_bw() +
    theme(
        text = element_text(family = 'DINPro'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7, hjust = 1)
        )

cairo_pdf(ccunspdf, width = 6, height = 9)
    p
dev.off()
