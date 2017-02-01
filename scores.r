#!/usr/bin/Rscript

# (c) Denes Turei 2017 EMBL

require(ggplot2)

infile <- 'scores.csv'

data <- read.csv(infile, sep = '\t', header = TRUE)

p <- ggplot(data, aes(x = boolenr, y = log(prs))) +
    geom_boxplot() +
    facet_wrap(protein ~ mode, scales = 'free_y') +
    ggtitle('Peak ratio score vs. Enric`s filter') +
    xlab('Enric`s filter (final)') +
    ylab('Peak ratio score (mean; log)')

ggsave('gg_prs_enric_box.pdf', device = cairo_pdf, width = 10, height = 20)


p <- ggplot(data, aes(x = boolenr, y = log(prsf))) +
geom_boxplot() +
facet_wrap(protein ~ mode, scales = 'free_y') +
ggtitle('Peak ratio score vs. Enric`s filter') +
xlab('Enric`s filter (final)') +
ylab('Peak ratio score (first fractions; log)')

ggsave('gg_prsf_enric_box.pdf', device = cairo_pdf, width = 10, height = 20)

p <- ggplot(data, aes(x = boolenr, y = log(prsh))) +
geom_boxplot() +
facet_wrap(protein ~ mode, scales = 'free_y') +
ggtitle('Peak ratio score vs. Enric`s filter') +
xlab('Enric`s filter (final)') +
ylab('Peak ratio score (highest diff; log)')

ggsave('gg_prsh_enric_box.pdf', device = cairo_pdf, width = 10, height = 20)

p <- ggplot(data, aes(x = boolppr, y = log(prs))) +
geom_boxplot() +
facet_wrap(protein ~ mode, scales = 'free_y') +
ggtitle('Peak ratio score vs. Protein ratio') +
xlab('Protein ratio (first fractions)') +
ylab('Protein ratio score (mean; log)')

ggsave('gg_prs_ppr_box.pdf', device = cairo_pdf, width = 10, height = 20)
