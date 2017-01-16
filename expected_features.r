#!/usr/bin/Rscript

# (c) Denes Turei 2016 EMBL

require(ggplot2)

infile <- 'expected_features.tab'
clinfile <- 'expected_features_closest.tab'

data <- read.csv(infile, sep = '\t', header = TRUE)
cldata <- read.csv(clinfile, sep = '\t', header = TRUE)

data$fr <- factor(data$fr, levels = levels(data$fr)[order(substr(levels(data$fr), 1, 1),
                                         as.numeric(substr(levels(data$fr), 2, 3)))])

cldata$fr <- factor(cldata$fr, levels = levels(cldata$fr)[order(substr(levels(cldata$fr), 1, 1),
                                             as.numeric(substr(levels(cldata$fr), 2, 3)))])

p <- ggplot(data, aes(y = i, x = fr)) +
    facet_wrap(protein ~ mode, scales = 'free') +
    geom_line(aes(group = fe), alpha = 0.4, color = '#049AA7') +
    xlab('Fraction') +
    ylab('Intensity')

ggsave('gg_fr_int_protein_mode_line.pdf', device = cairo_pdf, width = 12, height = 12)

p <- ggplot(cldata, aes(y = i, x = fr)) +
    facet_wrap(protein ~ mode, scales = 'free') +
    geom_line(aes(group = fe), alpha = 0.4, color = '#049AA7') +
    xlab('Fraction') +
    ylab('Intensity')

ggsave('gg_fr_int_protein_mode_line_closest.pdf', device = cairo_pdf, width = 10, height = 10)

###

infile <- 'nonzero_counts.tab'
data <- read.csv(infile, sep = '\t', header = TRUE)

data$fr <- factor(data$fr, levels = levels(data$fr)[order(substr(levels(data$fr), 1, 1),
                                         as.numeric(substr(levels(data$fr), 2, 3)))])

p <- ggplot(data, aes(y = cnt, x = fr, fill = peak)) +
    facet_wrap(~ protein + mode + set, scales = 'free') +
    geom_bar(stat = 'identity') +
    xlab('Fraction') +
    ylab('Non-zero values') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave('gg_nonzero_count_protein_mode_bar.pdf', device = cairo_pdf, width = 16, height = 24)
