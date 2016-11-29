#!/usr/bin/Rscript

library(ggplot2)
library(RColorBrewer)

cat("\t:: Reading data.\n")

getPalette = colorRampPalette(brewer.pal(12, "Paired"))

data <- read.csv('final_results.csv', sep = '\t', header = TRUE)

data['best'] <- data$cls %in% c('I', 'II')

adata <- aggregate(x = data, by = data[c('protein', 'id', 'ionm')], FUN = head, 1)

best_lab <- function(labs){
    lapply(labs, function(lab) ifelse(lab, 'Class I & II', 'Other classes'))
}

ionm_lab <- c(
    'pos' = 'Positive',
    'neg' = 'Negative'
)

# class I and II results
data12 <- data[data$cls %in% c('I', 'II'),]
adata12 <- adata[adata$cls %in% c('I', 'II'),]

# filter empty:
nonempty <- NULL
for(pr in levels(adata12$protein)){
    if(dim(adata12[adata12$protein == pr,])[1] != 0){
        nonempty <- c(nonempty, pr)
    }
}

adata12v <- adata12[adata12$protein %in% nonempty,]

for(col in colnames(adata12v)){
    if(is.factor(adata12v[[col]])){
        adata12v[[col]] <- factor(adata12v[[col]])
    }
}

# those having carbon counts:
adata12cc <- adata12v[!is.na(adata12v$cc),]
for(col in colnames(adata12cc)){
    if(is.factor(adata12cc[[col]])){
        adata12cc[[col]] <- factor(adata12cc[[col]])
    }
}

# those having fatty acids:
adata12ccfa <- adata12v[!is.na(adata12v$ccfa),]
for(col in colnames(adata12ccfa)){
    if(is.factor(adata12ccfa[[col]])){
        adata12ccfa[[col]] <- factor(adata12ccfa[[col]])
    }
}

datacc <- data[((!is.na(data$cc)) & data$cls != ""),]
datacc$cls <- factor(datacc$cls)

#
# plots:
#

cat("\t:: Plotting `Mean retention time by ion mode and class`.\n")
ggplot(adata, aes(x = 1, y = rtmean)) +
    geom_violin() +
    facet_grid(ionm ~ best, labeller = labeller(ionm = ionm_lab, best = best_lab)) +
    xlab('Category') +
    ylab('Mean RT [min]') +
    ggtitle('Mean retention time\nby ion mode and class') +
    theme(axis.title = element_text(size = 24),
          axis.text = element_text(size = 18),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 21),
          strip.text = element_text(size = 18),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())

ggsave('gg_rt_cat_viol.pdf', device = cairo_pdf, width = 6, height = 12)

cat("\t:: Plotting `Mean retention time by ion mode, class and headgroup`.\n")
ggplot(adata, aes(x = 1, y = rtmean, fill = uhgroup)) +
geom_dotplot(binaxis = 'y', stackdir = 'centerwhole', dotsize = 1, binwidth = 0.25, method = 'histodot') +
facet_grid(ionm ~ best, labeller = labeller(ionm = ionm_lab, best = best_lab)) +
scale_fill_discrete(name = 'Headgroup') +
xlab('Category') +
ylab('Mean RT [min]') +
ggtitle('Mean retention time\nby ion mode, class and headgroup') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_rt_cat_hg_dot.pdf', device = cairo_pdf, width = 12, height = 12)

cat("\t:: Plotting `Mean retention time by ion mode and class`.\n")
ggplot(adata, aes(x = 1, y = rtmean)) +
geom_boxplot() +
facet_grid(ionm ~ best, labeller = labeller(ionm = ionm_lab, best = best_lab)) +
xlab('Category') +
ylab('Mean RT [min]') +
ggtitle('Mean retention time\nby ion mode and class') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_rt_cat_box.pdf', device = cairo_pdf, width = 6, height = 12)

cat("\t:: Plotting `Mean retention time by protein`.\n")
ggplot(adata12, aes(x = protein, y = rtmean)) +
geom_boxplot() +
xlab('Protein') +
ylab('Mean RT [min]') +
ggtitle('Mean retention time by protein') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(size = 36))

ggsave('gg_rt_protein.pdf', device = cairo_pdf, width = 24, height = 9)

cat("\t:: Plotting `Carbon count and unsaturation by protein and ion mode`.\n")
ggplot(adata12, aes(x = unsat, y = carb)) +
geom_point(alpha = 0.5) +
scale_x_continuous(breaks = seq(0, max(adata12$unsat[!is.nan(adata12$unsat)]), 2)) +
facet_grid(protein ~ ionm, labeller = labeller(ionm = ionm_lab)) +
xlab('Unsaturation') +
ylab('Carbon count') +
ggtitle('Carbon count and unsaturation\nby protein and ion mode') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      plot.title = element_text(size = 18),
      strip.text = element_text(size = 18))

ggsave('gg_cc_unsat_protein_mode.pdf', device = cairo_pdf, width = 4, height = 48)

cat("\t:: Plotting `Carbon count and unsaturation`.\n")
ggplot(adata12[!is.nan(adata12$unsat),], aes(x = as.factor(unsat), y = carb)) +
geom_violin() +
xlab('Unsaturation') +
ylab('Carbon count') +
ggtitle('Carbon count and unsaturation') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      plot.title = element_text(size = 24),
      strip.text = element_text(size = 18))

ggsave('gg_cc_unsat_viol.pdf', device = cairo_pdf, width = 12, height = 6)

cat("\t:: Plotting `Carbon counts`.\n")
ggplot(adata12[!is.nan(adata12$unsat),], aes(carb, fill = ionm)) +
geom_histogram(binwidth = 1) +
geom_density(kernel = 'gaussian') +
scale_x_discrete(breaks = seq(0, max(adata12$carb[!is.nan(adata12$carb)]), 5)) +
scale_fill_discrete(name = 'Ion mode', breaks = c('neg', 'pos'), labels = c('Negative', 'Positive')) +
xlab('Carbon counts') +
ylab('Frequency') +
ggtitle('Carbon counts') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      plot.title = element_text(size = 24),
      strip.text = element_text(size = 18),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 21),
      legend.key = element_blank())

ggsave('gg_cc_mode_hist.pdf', device = cairo_pdf, width = 12, height = 6)

cat("\t:: Plotting `Carbon counts`.\n")
ggplot(adata12[!is.nan(adata12$unsat) & adata12$unsat < 8,], aes(carb, fill = as.factor(unsat))) +
#geom_histogram(binwidth = 1) +
geom_density(kernel = 'gaussian', alpha = 0.5) +
scale_x_discrete(breaks = seq(0, max(adata12$carb[!is.nan(adata12$carb)]), 5)) +
scale_fill_discrete(name = 'Unsaturation', breaks = seq(0,7), labels = seq(0, 7)) +
xlab('Carbon counts') +
ylab('Frequency') +
ggtitle('Carbon counts') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      plot.title = element_text(size = 24),
      strip.text = element_text(size = 18),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 21),
      legend.key = element_blank())

ggsave('gg_cc_unsat_dens.pdf', device = cairo_pdf, width = 12, height = 4)

cat("\t:: Plotting `Carbon count and unsaturation by protein and ion mode`.\n")
p <- ggplot(adata12, aes(x = unsat, y = carb, size = log10(intensity), color = uhgroup)) +
geom_point(alpha = 0.33, stroke = 0) +
scale_color_discrete(guide = guide_legend(title = 'Headgroup')) +
scale_size_continuous(guide = guide_legend(title = 'Intensity (log)')) +
scale_x_continuous(breaks = seq(0, max(adata12$unsat[!is.nan(adata12$unsat)]), 2)) +
scale_radius(range=c(1,6)) +
facet_grid(protein ~ ionm, labeller = labeller(ionm = ionm_lab)) +
xlab('Unsaturation') +
ylab('Carbon count') +
ggtitle('Carbon count and unsaturation\nby protein and ion mode') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      plot.title = element_text(size = 18),
      strip.text = element_text(size = 18))

ggsave('gg_cc_unsat_int_protein_mode.pdf', device = cairo_pdf, width = 5.33, height = 48)

cat("\t:: Plotting `Carbon count, unsaturation, intensity and headgroup\nby protein and ion mode`.\n")
p <- ggplot(adata12, aes(x = unsat, y = carb, size = log10(intensity), color = uhgroup)) +
geom_text(aes(label=uhgroup), alpha = 0.3, nudge_x = 0, nudge_y = 0) +
geom_point(size = 1, alpha = 1.0, shape = 43) +
scale_x_continuous(breaks = seq(0, max(adata12$unsat[!is.nan(adata12$unsat)]), 2)) +
scale_color_discrete(guide = guide_legend(title = 'Headgroup', title.position = 'top', ncol = 4, order = 2)) +
scale_size_continuous(guide = guide_legend(title = 'Intensity (log)', title.position = 'top', order = 1)) +
#scale_radius(range=c(1,6)) +
facet_grid(protein ~ ionm, labeller = labeller(ionm = ionm_lab)) +
xlab('Unsaturation') +
ylab('Carbon count') +
ggtitle('Carbon count, unsaturation,\nintensity and headgroup\nby protein and ion mode') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      plot.title = element_text(size = 18),
      strip.text = element_text(size = 18),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 12),
      legend.position = 'bottom')

ggsave('gg_cc_unsat_int_protein_hg_mode.pdf', device = cairo_pdf, width = 4.33, height = 48)

cat("\t:: Plotting `Intensity and headgroup by protein and ion mode`.\n")
p <- ggplot(adata12, aes(x = factor(1), y = intensity, fill = uhgroup, label = uhgroup)) +
geom_bar(width = 1, stat = 'sum', position = 'fill') +
scale_size_continuous(guide = FALSE) +
#annotate(geom = "text") +
coord_polar(theta = 'y') +
scale_fill_manual(values = getPalette(length(levels(as.factor(adata12$uhgroup)))),
                  guide = guide_legend(title = 'Headgroup', title.position = 'top', ncol = 3, order = 2)) +
#geom_text(aes(label = uhgroup, y = intensity, x = factor(12)), stat = 'sum', position = position_dodge(1)) +
facet_grid(protein ~ ionm, labeller = labeller(ionm = ionm_lab), scales = 'free') +
xlab(NULL) +
ylab('Intensity') +
ggtitle('Intensity and headgroup\nby protein and ion mode') +
theme(axis.text.x = element_blank(),
      axis.title = element_text(size = 18),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = 16),
      strip.text = element_text(size = 18),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 12),
      legend.position = 'bottom')

ggsave('gg_int_protein_hg_mode.pdf', device = cairo_pdf, width = 3, height = 49)

cat("\t:: Plotting `Mean retention time by headgroup`.\n")
ggplot(adata12, aes(x = uhgroup, y = rtmean)) +
geom_boxplot() +
#facet_grid(ionm, labeller = labeller(ionm = ionm_lab, best = best_lab)) +
xlab('Headgroup') +
ylab('Mean RT [min]') +
ggtitle('Mean retention time\nby headgroup') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      #axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      #panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_rt_hg_boxplot.pdf', device = cairo_pdf, width = 12, height = 6)

cat("\t:: Plotting `Mean retention time by carbon count`.\n")
ggplot(adata12, aes(x = carb, y = rtmean)) +
geom_point(shape = 1) +
geom_smooth(method = lm, se = FALSE) +
#facet_grid(ionm, labeller = labeller(ionm = ionm_lab, best = best_lab)) +
xlab('Carbon count') +
ylab('Mean RT [min]') +
ggtitle('Mean retention time by carbon count') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.minor.x = element_blank())

ggsave('gg_rt_cc_dot.pdf', device = cairo_pdf, width = 12, height = 6)

cat("\t:: Plotting `Mean retention time by carbon count and headgroup`.\n")
ggplot(adata12cc, aes(x = hgcc, y = rtmean)) +
geom_boxplot() +
xlab('Headgroup') +
ylab('Mean RT [min]') +
ggtitle('Mean retention time by carbon count and headgroup') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.minor.x = element_blank())

ggsave('gg_rt_hg_cc_boxplot.pdf', device = cairo_pdf, width = 36, height = 6)

cat("\t:: Plotting `Mean retention time by carbon count, headgroup and class`.\n")
ggplot(adata12cc, aes(x = hgcc, y = rtmean, color = cls)) +
geom_point(shape = 'o', size = 3) +
scale_color_discrete(guide = guide_legend(title = 'Class')) +
xlab('Headgroup') +
ylab('Mean RT [min]') +
ggtitle('Mean retention time by carbon count and headgroup') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.minor.x = element_blank())

ggsave('gg_rt_hg_cc_class_point.pdf', device = cairo_pdf, width = 36, height = 6)

cat("\t:: Plotting `Mean retention time by carbon count, headgroup, class and protein`.\n")
ggplot(adata12cc, aes(x = hgcc, y = rtmean, color = cls, label = protein)) +
geom_text(size = 1) +
scale_color_discrete(guide = guide_legend(title = 'Class')) +
xlab('Headgroup') +
ylab('Mean RT [min]') +
ggtitle('Mean retention time by carbon count, headgroup and protein') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.minor.x = element_blank())

ggsave('gg_rt_hg_cc_class_protein_point.pdf', device = cairo_pdf, width = 36, height = 6)

cat("\t:: Plotting `Mean retention time by carbon count, headgroup and class`.\n")
ggplot(datacc, aes(x = hgcc, y = rtmean, color = cls)) +
geom_point(shape = 'o', size = 2) +
scale_color_discrete(guide = guide_legend(title = 'Class')) +
xlab('Headgroup') +
ylab('Mean RT [min]') +
ggtitle('Mean retention time by carbon count, headgroup and class') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.minor.x = element_blank())

ggsave('gg_all_rt_hg_cc_class_point.pdf', device = cairo_pdf, width = 49, height = 6)


cat("\t:: Plotting `Mean retention time by carbon count and headgroup`.\n")
ggplot(adata12, aes(x = carb, y = rtmean, color = uhgroup, label = uhgroup)) +
geom_point(shape = '+', alpha = 1.0) +
geom_text(alpha = 0.5) +
scale_color_discrete(guide = guide_legend(title = 'Headgroup')) +
xlab('Carbon count') +
ylab('Mean RT [min]') +
ggtitle('Mean retention time by carbon count and headgroup') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      #axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      #panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_rt_cc_hg_dot.pdf', device = cairo_pdf, width = 12, height = 6)

##############################
#
cat("\t:: Plotting `Headgroup vs. carbon count & unsaturation`.\n")
ggplot(adata12cc, aes(x = cc, y = uhgroup)) +
geom_count() +
#scale_color_discrete(guide = guide_legend(title = 'Headgroup')) +
#facet_grid(ionm, labeller = labeller(ionm = ionm_lab, best = best_lab)) +
scale_size_continuous(guide = guide_legend(title = 'Count')) +
xlab('Carbon count & unsaturation') +
ylab('Headgroup') +
ggtitle('Headgroup vs. carbon count & unsaturation') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8, hjust = 1),
      #axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      #panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_cc_hg_size.pdf', device = cairo_pdf, width = 18, height = 9)

#
cat("\t:: Plotting `Headgroup vs. fatty acids`.\n")
ggplot(adata12ccfa, aes(x = ccfa, y = uhgroup)) +
geom_count() +
#scale_color_discrete(guide = guide_legend(title = 'Headgroup')) +
#facet_grid(ionm, labeller = labeller(ionm = ionm_lab, best = best_lab)) +
scale_size_continuous(guide = guide_legend(title = 'Count')) +
xlab('Fatty acids') +
ylab('Headgroup') +
ggtitle('Headgroup vs. fatty acids') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14, hjust = 1),
      #axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      #panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_facc_hg_size.pdf', device = cairo_pdf, width = 18, height = 9)

#
cat("\t:: Plotting `Protein vs. carbon count & unsaturation'`.\n")
ggplot(adata12cc, aes(x = cc, y = protein)) +
geom_count() +
#scale_color_discrete(guide = guide_legend(title = 'Headgroup')) +
#facet_grid(ionm, labeller = labeller(ionm = ionm_lab, best = best_lab)) +
scale_size_continuous(guide = guide_legend(title = 'Count')) +
xlab('Carbon count & unsaturation') +
ylab('Protein') +
ggtitle('Protein vs. carbon count & unsaturation') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
      #axis.text.y = element_text(size = 9),
      #axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      #panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_cc_protein_size.pdf', device = cairo_pdf, width = 18, height = 12)

cat("\t:: Plotting `Protein vs. fatty acids`.\n")
ggplot(adata12ccfa, aes(x = ccfa, y = protein)) +
geom_count() +
#scale_color_discrete(guide = guide_legend(title = 'Headgroup')) +
#facet_grid(ionm, labeller = labeller(ionm = ionm_lab, best = best_lab)) +
scale_size_continuous(guide = guide_legend(title = 'Count')) +
xlab('Fatty acids') +
ylab('Protein') +
ggtitle('Protein vs. fatty acids') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14, hjust = 1),
      axis.text.y = element_text(size = 14),
      #axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      #panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_facc_protein_size.pdf', device = cairo_pdf, width = 18, height = 6)

#
cat("\t:: Plotting `Protein vs. headgroup, carbon count & unsaturation`.\n")
ggplot(adata12cc, aes(x = hgcc, y = protein, color = uhgroup)) +
geom_count() +
#scale_color_discrete(guide = guide_legend(title = 'Headgroup')) +
#facet_grid(ionm, labeller = labeller(ionm = ionm_lab, best = best_lab)) +
scale_size_continuous(guide = guide_legend(title = 'Count')) +
scale_color_manual(values = getPalette(length(levels(as.factor(adata12$uhgroup)))),
                  guide = guide_legend(title = 'Headgroup')) +
xlab('Headgroup, carbon count & unsaturation') +
ylab('Protein') +
ggtitle('Protein vs. headgroup, carbon count & unsaturation') +
theme(axis.title = element_text(size = 24),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9, hjust = 1),
      axis.text.y = element_text(size = 14),
      #axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      #panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_hg_cc_protein_size.pdf', device = cairo_pdf, width = 30, height = 9)

#
cat("\t:: Plotting `Protein vs. headgroup & fatty acids`.\n")
ggplot(adata12ccfa, aes(x = hgfa, y = protein, color = uhgroup)) +
geom_count() +
scale_size_continuous(guide = guide_legend(title = 'Count')) +
scale_color_discrete(guide = guide_legend(title = 'Headgroup')) +
xlab('Headgroup & fatty acids') +
ylab('Protein') +
ggtitle('Protein vs. headgroup & fatty acids') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14, hjust = 1),
      axis.text.y = element_text(size = 14),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.minor.x = element_blank(),
      legend.position = 'bottom')

ggsave('gg_hg_facc_protein_size.pdf', device = cairo_pdf, width = 18, height = 6)

##############################
#
cat("\t:: Plotting `Headgroup vs. carbon count & unsaturation vs. intensity`.\n")
ggplot(adata12cc, aes(x = cc, y = uhgroup, size = log10(intensity))) +
geom_point(stat = 'sum') +
scale_radius(guide = guide_legend(title = 'Intensity (log)')) +
xlab('Carbon count & unsaturation') +
ylab('Headgroup') +
ggtitle('Headgroup vs. carbon count & unsaturation vs. intensity') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8, hjust = 1),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.minor.x = element_blank())

ggsave('gg_cc_hg_int_size.pdf', device = cairo_pdf, width = 18, height = 9)

#
cat("\t:: Plotting `Headgroup vs. fatty acids vs. intensity`.\n")
ggplot(adata12ccfa, aes(x = ccfa, y = uhgroup, size = log10(intensity))) +
geom_point(stat = 'sum') +
scale_radius(guide = guide_legend(title = 'Intensity (log)')) +
xlab('Fatty acids') +
ylab('Headgroup') +
ggtitle('Headgroup vs. fatty acids vs. intensity') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14, hjust = 1),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.minor.x = element_blank())

ggsave('gg_facc_hg_int_size.pdf', device = cairo_pdf, width = 18, height = 9)

#
cat("\t:: Plotting `Protein vs. carbon count & unsaturation vs. intensity`.\n")
ggplot(adata12cc, aes(x = cc, y = protein, size = log10(intensity))) +
geom_point(stat = 'sum') +
scale_radius(guide = guide_legend(title = 'Intensity (log)')) +
xlab('Carbon count & unsaturation') +
ylab('Protein') +
ggtitle('Protein vs. carbon count & unsaturation vs. intensity') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.minor.x = element_blank())

ggsave('gg_cc_protein_int_size.pdf', device = cairo_pdf, width = 18, height = 12)

###
#
cat("\t:: Plotting `Protein vs. carbon count & unsaturation vs. intensity`.\n")
p <- ggplot(datacc, aes(x = cc, y = protein, size = log10(intensity), color = cls)) +
geom_point(stat = 'sum', alpha = 0.4) +
scale_radius(guide = guide_legend(title = 'Intensity (log)')) +
xlab('Carbon count & unsaturation') +
ylab('Protein') +
ggtitle('Protein vs. carbon count & unsaturation vs. intensity') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.minor.x = element_blank())

ggsave('gg_all_cc_protein_int_class_size.pdf', device = cairo_pdf, width = 18, height = 12)


#
cat("\t:: Plotting `Protein vs. carbon count & unsaturation vs. intensity`.\n")
ggplot(adata12cc, aes(x = cc, y = protein, size = log10(intensity), color = cls)) +
geom_point(stat = 'sum', alpha = 0.4) +
scale_radius(guide = guide_legend(title = 'Intensity (log)')) +
xlab('Carbon count & unsaturation') +
ylab('Protein') +
ggtitle('Protein vs. carbon count & unsaturation vs. intensity') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.minor.x = element_blank())

ggsave('gg_cc_protein_int_class_size.pdf', device = cairo_pdf, width = 18, height = 12)

#
cat("\t:: Plotting `Protein vs. carbon count & unsaturation vs. intensity`.\n")
ggplot(adata12cc, aes(x = cc, y = protein, size = log10(intensity), color = cls)) +
geom_point(stat = 'sum', alpha = 0.4) +
scale_radius(guide = guide_legend(title = 'Intensity (log)')) +
xlab('Carbon count & unsaturation') +
ylab('Protein') +
ggtitle('Protein vs. carbon count & unsaturation vs. intensity') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.minor.x = element_blank())

ggsave('gg_cc_protein_int_class_hg_size.pdf', device = cairo_pdf, width = 18, height = 12)

#
cat("\t:: Plotting `Protein vs. carbon count & unsaturation vs. intensity`.\n")
p <- ggplot(adata12cc, aes(x = cc, y = protein, size = log10(intensity), color = cls, label = uhgroup)) +
geom_text(stat = 'sum', alpha = 0.4) +
scale_radius(guide = guide_legend(title = 'Intensity (log)')) +
xlab('Carbon count & unsaturation') +
ylab('Protein') +
ggtitle('Protein vs. carbon count & unsaturation vs. intensity') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.minor.x = element_blank())

ggsave('gg_cc_protein_int_class_hg_text.pdf', device = cairo_pdf, width = 18, height = 12)

###
#
cat("\t:: Plotting `Protein vs. fatty acids vs. intensity`.\n")
ggplot(adata12ccfa, aes(x = ccfa, y = protein, size = log10(intensity))) +
geom_point(stat = 'sum') +
scale_radius(guide = guide_legend(title = 'Intensity (log)')) +
xlab('Fatty acids') +
ylab('Protein') +
ggtitle('Protein vs. fatty acids vs. intensity') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14, hjust = 1),
      axis.text.y = element_text(size = 14),
      #axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      #panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_facc_protein_int_size.pdf', device = cairo_pdf, width = 18, height = 6)

#
cat("\t:: Plotting `Protein vs. headgroup, carbon count & unsaturation vs. intensity`.\n")
ggplot(adata12cc, aes(x = hgcc, y = protein, color = uhgroup, size = log10(intensity))) +
geom_point(stat = 'sum') +
scale_radius(guide = guide_legend(title = 'Intensity (log)')) +
scale_color_manual(values = getPalette(length(levels(as.factor(adata12$uhgroup)))),
                   guide = guide_legend(title = 'Headgroup')) +
xlab('Headgroup, carbon count & unsaturation') +
ylab('Protein') +
ggtitle('Protein vs. headgroup, carbon count & unsaturation vs. intensity') +
theme(axis.title = element_text(size = 24),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9, hjust = 1),
      axis.text.y = element_text(size = 14),
      #axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      #panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_hg_cc_protein_int_size.pdf', device = cairo_pdf, width = 30, height = 9)

###
#
cat("\t:: Plotting `Protein vs. headgroup, carbon count & unsaturation vs. intensity vs. class`.\n")
p <- ggplot(adata12cc, aes(x = hgcc, y = protein, color = cls, size = log10(intensity))) +
geom_point(stat = 'sum', alpha = 0.33) +
scale_radius(guide = guide_legend(title = 'Intensity (log)')) +
scale_color_discrete(guide = guide_legend(title = 'Class')) +
xlab('Headgroup, carbon count & unsaturation') +
ylab('Protein') +
ggtitle('Protein vs. headgroup, carbon count & unsaturation vs. intensity') +
theme(axis.title = element_text(size = 24),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9, hjust = 1),
      axis.text.y = element_text(size = 14),
      #axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      #panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_hg_cc_protein_int_class_size.pdf', device = cairo_pdf, width = 30, height = 9)

#
cat("\t:: Plotting `Protein vs. headgroup, carbon count & unsaturation vs. intensity vs. class`.\n")
p <- ggplot(datacc, aes(x = hgcc, y = protein, color = cls, size = log10(intensity))) +
geom_point(stat = 'sum', alpha = 0.33) +
scale_radius(guide = guide_legend(title = 'Intensity (log)')) +
scale_color_discrete(guide = guide_legend(title = 'Class')) +
xlab('Headgroup, carbon count & unsaturation') +
ylab('Protein') +
ggtitle('Protein vs. headgroup, carbon count & unsaturation vs. intensity') +
theme(axis.title = element_text(size = 24),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5, hjust = 1),
      axis.text.y = element_text(size = 11),
      #axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      #panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())

ggsave('gg_all_hg_cc_protein_int_class_size.pdf', device = cairo_pdf, width = 49, height = 9)
###

#
cat("\t:: Plotting `Protein vs. headgroup & fatty acids vs. intensity`.\n")
ggplot(adata12ccfa, aes(x = hgfa, y = protein, color = uhgroup, size = log10(intensity))) +
geom_point(stat = 'sum') +
scale_radius(guide = guide_legend(title = 'Intensity (log)')) +
scale_color_discrete(guide = guide_legend(title = 'Headgroup', ncol = 13)) +
xlab('Headgroup & fatty acids') +
ylab('Protein') +
ggtitle('Protein vs. headgroup & fatty acids vs. intensity') +
theme(axis.title = element_text(size = 24),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14, hjust = 1),
      axis.text.y = element_text(size = 14),
      #axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      #panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = 'bottom')

ggsave('gg_hg_facc_protein_int_size.pdf', device = cairo_pdf, width = 18, height = 8)

#
cat("\t:: Plotting `Protein vs. headgroup vs. intensity`.\n")
ggplot(adata12cc, aes(x = uhgroup, y = protein, fill = log10(intensity))) +
geom_raster(stat = 'sum') +
scale_fill_continuous(guide = guide_legend(title = 'Intensity (log)')) +
xlab('Headgroup') +
ylab('Protein') +
ggtitle('Protein vs. headgroup vs. intensity') +
theme_bw() +
theme(axis.title = element_text(size = 24),
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 11, hjust = 1),
      axis.text.y = element_text(size = 11),
      #axis.ticks.x = element_blank(),
      plot.title = element_text(size = 21),
      strip.text = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank())

ggsave('gg_hg_protein_int_raster.pdf', device = cairo_pdf, width = 9, height = 9)

###
#
my.formula <- y ~ order(as.factor(x))

cat("\t:: Plotting `RT vs. carbon count & unsaturation by headgroup`.\n")
p <- ggplot(adata12cc, aes(x = as.integer(carb), y = rtmean, color = cls, size = log10(intensity))) +
geom_point(stat = 'sum', alpha = 0.33) +
stat_smooth(fullrange = TRUE, show.legend = FALSE, se = TRUE, method = lm, linetype = 2, geom = 'smooth') +
facet_wrap( ~ uhgroup, scales = 'free') +
scale_radius(guide = guide_legend(title = 'Intensity (log)')) +
scale_color_discrete(guide = guide_legend(title = 'Class')) +
xlab('Carbon count unsaturation') +
ylab('Mean RT [min]') +
ggtitle('RT vs. carbon count & unsaturation by headgroup') +
theme(axis.title = element_text(size = 24),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 9, hjust = 1),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 21),
    strip.text = element_text(size = 18),
    panel.grid.minor.x = element_blank())

ggsave('gg_rt_hg_cc_unsat_int_class_lm_size.pdf', device = cairo_pdf, width = 30, height = 30)
