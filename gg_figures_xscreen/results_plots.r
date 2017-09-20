#!/usr/bin/Rscript

###
cat("\t:: Populating workspace.\n")

library(ggplot2)
library(RColorBrewer)
library(colorspace)
library(dplyr)

source('theme_black.R')

onlyclass1 <- TRUE
screen1file <- 'antonella_final.csv'
screen2file <- 'enric_processed.csv'

getPalette = colorRampPalette(brewer.pal(12, "Paired"))

# code the colorkey palette
pal <- colorRampPalette(c("white", "dark blue"),space="rgb")

# here we just get a big pool of distinct colors
# later we use as many as we need
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
distinct_cols <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors,
                               rownames(qual_col_pals)))

distinct_cols2 <- c(
    "#000000",
    #"#FFFF00",
    "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    #"#FFDBE5",
    "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693",
    #"#FEFFE6",
    "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
    # "#DDEFFF",
    # "#000035",
    "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C"
)

best_lab <- function(labs){
    lapply(labs, function(lab) ifelse(lab, ifelse(onlyclass1, 'Class I', 'Class I & II'), 'Other classes'))
}

ionm_lab <- c(
    'pos' = 'Positive',
    'neg' = 'Negative'
)

screen_lab <- c(
    'A' = 'HEK cells',
    'AE' = 'Both',
    'E' = 'E. coli & liposomes'
)

interleave <- function(v1,v2){
    ord1 <- 2*(1:length(v1))-1
    ord2 <- 2*(1:length(v2))
    c(v1,v2)[order(c(ord1,ord2))]
}

protein_colors <- function(proteins){
    return(
        as.vector(
            sapply(
                as.vector(
                    sort(
                        unique(proteins)
                    )
                ),
                function(p) cprotein[[p]]
            )
        )
    )
}

sort_factors <- function(df){
    for(col in colnames(df)){
        if(is.factor(df[[col]])){
            df[[col]] <- factor(df[[col]],
                                levels = sort(levels(df[[col]])))
        }
    }
    return(df)
}

a <- function(...){
    print(hasArg(i))
    cat(...)
}

merge_screens <- function(df, ...){
    df_g <-
        df %>%
        group_by(...) %>%
        mutate(screens = paste(sort(unique(screen)), collapse = "")) %>%
        summarise_all(first)
    df_g$screens <- factor(df_g$screens)
    return(sort_factors(df_g))
}

merge_intensities <- function(df, ...){
    df_g <-
        df %>%
        group_by(..., protein, ionm, id, screen) %>%
        summarise_all(first) %>%
        group_by(..., screen) %>%
        mutate(isum = sum(nintensity)) %>%
        mutate(cnt  = n()) %>%
        summarise_all(first)
    return(sort_factors(df_g))
}

merge_screens2 <- function(df, ...){
    # aggregate the data frame by grouping vars plus `screen`
    df_g             <- df %>%
        group_by(..., screen)
    # sorting all factors
    df_g             <- sort_factors(df_g)
    # grouping by grouping variables only and concatenating values of `screen`
    df_screen        <- df_g %>%
        group_by(...) %>%
        summarise(screen = paste(screen, collapse = ""))
    # casting `screen` to factor
    df_screen$screen <- factor(df_screen$screen)
    # do the same grouping on the remaining columns here taking the first value within groups
    df_othercols     <- df_g %>%
        group_by(...) %>%
        summarise_all(first)
    # changing the `screen` col of this dataframe with the combined one
    df_merged        <- bind_cols(df_othercols %>% select(-screen),
                                  as.data.frame(df_screen['screen']))
    return(df_merged)
}

norm_intensity <- function(df){
    isum <- df %>%
        group_by(protein, ionm) %>%
        summarise(nintensity = sum(intensity))
    df <- df %>%
        inner_join(isum, by = c('protein', 'ionm'))
    df$nintensity <- df$intensity / df$nintensity
    return(df)
}

###
cat("\t:: Reading data.\n")

data1 <- read.csv(screen1file, sep = '\t', header = TRUE)
data1 <- norm_intensity(data1)
data2 <- read.csv(screen2file, sep = '\t', header = TRUE)
data2 <- norm_intensity(data2)

cat("\t:: Preprocessing data.\n")

data <- rbind(data1, data2)

data$protein_ionm <- paste(data$protein, data$ionm, sep = "_")

cprotein <- list()
for(protein in unique(data$protein)){
    screens <- unique(data[data$protein == protein,]$screen)
    cprotein[[protein]] <- ifelse('A' %in% screens, ifelse('E' %in% screens, 'black', 'cyan'), 'magenta')
}

if(onlyclass1){
    data['best'] <- data$cls == 'I'
    pname <- 'onlyI_'
}else{
    data['best'] <- data$cls %in% c('I', 'II')
    pname <- ''
}

data['inscreen'] <- as.vector(sapply(data$protein, function(p) cprotein[[p]]))

# aggregate by protein, ion mode, feature and screening
adata <- merge_screens(data, protein, ionm, id)

# class I and II results
data12 <- data %>% filter(best)

# class I & II with screens merged
adata12_s <- merge_screens(data12, protein, ionm, id)
adata12 <- data12 %>% group_by(protein, ionm, id, screen) %>% summarise_all(first)
adata12 <- sort_factors(adata12)

datacc       <- data   %>% filter(!is.na(cc))
data12cc     <- data12 %>% filter(!is.na(cc))
dataccfa     <- data   %>% filter(!is.na(ccfa))
data12ccfa   <- data12 %>% filter(!is.na(ccfa))

adata12hg      <- merge_screens(data12, protein, uhgroup)
adata12ihg     <- merge_screens(data12, protein, ionm, uhgroup)
adata12hgi     <- merge_intensities(data12, protein, uhgroup)
adata12hgcc    <- merge_screens(data12cc, protein, uhgroup, cc)
adata12ihgcc   <- merge_screens(data12cc, protein, ionm, uhgroup, cc)
adata12hgcci   <- merge_intensities(data12, protein, uhgroup, cc)
adata12cc      <- merge_screens(data12cc, protein, cc)
adata12icc     <- merge_screens(data12cc, protein, ionm, cc)
adata12hgccfa  <- merge_screens(data12ccfa, protein, uhgroup, cc, ccfa)
adata12ihgccfa <- merge_screens(data12ccfa, protein, ionm, uhgroup, cc, ccfa)
adata12ccfa    <- merge_screens(data12ccfa, protein, cc, ccfa)
adata12iccfa   <- merge_screens(data12ccfa, protein, ionm, cc, ccfa)

adatahg     <- merge_screens(data, protein, ionm, uhgroup)
adatahgcc   <- merge_screens(datacc, protein, ionm, uhgroup, cc)
adatahgccfa <- merge_screens(dataccfa, protein, ionm, uhgroup, cc, ccfa)

#
# plots:
#

cat("\t:: Plotting `Headgroup vs. carbon count & unsaturation vs. intensity`.\n")
p <- ggplot(adata12hgi, aes(x = uhgroup, y = protein,
                            color = screen, size = log10(isum))) +
geom_point(alpha = 0.33) +
scale_radius(guide = guide_legend(title = 'Intensity (log)')) +
scale_color_manual(values = c('cyan', 'magenta'),
                   guide = guide_legend(title = 'Screening'),
                   labels = c('HEK cells', 'E. coli & liposomes')) +
#scale_shape_manual(values = c(3, 4), guide = FALSE) +
xlab('Headgroup') +
ylab('Protein') +
ggtitle('Intensities in screenings by proteins vs. headgroups') +
theme_bw() +
theme(
    text = element_text(family = "DINPro"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, hjust = 1),
    axis.text.y = element_text(size = 14, color = protein_colors(adata12hgi$protein)),
    plot.title = element_text(size = 21),
    strip.text = element_text(size = 18),
    panel.grid.minor.x = element_blank())

ggsave(paste('gg_', pname, 'protein_hg_int_screen_size.pdf', sep = ''), device = cairo_pdf, width = 12, height = 12)

cat("\t:: Plotting `Headgroup vs. carbon count & unsaturation vs. count`.\n")
p <- ggplot(adata12hgi, aes(x = uhgroup, y = protein,
                            color = screen, size = cnt)) +
geom_point(alpha = 0.33) +
scale_radius(guide = guide_legend(title = 'Count of distinct species')) +
scale_color_manual(values = c('cyan', 'magenta'),
                   guide = guide_legend(title = 'Screening'),
                   labels = c('HEK cells', 'E. coli & liposomes')) +
#scale_shape_manual(values = c(3, 4), guide = FALSE) +
xlab('Headgroup') +
ylab('Protein') +
ggtitle('Counts in screenings by proteins vs. headgroups') +
theme_bw() +
theme(
    text = element_text(family = "DINPro"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, hjust = 1),
    axis.text.y = element_text(size = 14, color = protein_colors(adata12hgi$protein)),
    plot.title = element_text(size = 21),
    strip.text = element_text(size = 18),
    panel.grid.minor.x = element_blank())

ggsave(paste('gg_', pname, 'protein_hg_cnt_screen_size.pdf', sep = ''), device = cairo_pdf, width = 12, height = 12)


cat("\t:: Plotting `Protein vs. fatty acids by screenings`.\n")
p <- ggplot(adata12ccfa, aes(x = ccfa, y = protein, shape = screens)) +
    geom_point(
        #alpha = .5
        ) +
    #scale_color_discrete(guide = guide_legend(title = 'Headgroup')) +
    #facet_grid(ionm, labeller = labeller(ionm = ionm_lab, best = best_lab)) +
    #scale_size_continuous(guide = guide_legend(title = 'Count')) +
    scale_shape_manual(
                    values = c(65, 13, 69),
                    labels = as.vector(screen_lab),
                    guide = guide_legend(title = 'Screening')) +
    xlab('Fatty acids') +
    ylab('Protein') +
    ggtitle('Protein vs. fatty acids by screenings') +
    theme_bw() +
    theme(
        text = element_text(family = "DINPro"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
        axis.text.y = element_text(size = 9,
                                   color = protein_colors(adata12ccfa$protein)),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 21),
        strip.text = element_text(size = 18),
        #panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggsave(paste('gg_', pname, 'facc_protein_screen_point.pdf', sep = ''), device = cairo_pdf, width = 20, height = 8)

#

cat("\t:: Plotting `Protein vs. headgroup, carbon count & unsaturation by screening` with facets.\n")
p <- ggplot(adata12hgcc, aes(x = cc, y = protein, shape = screens)) +
    geom_point(size = 3) +
    #scale_color_discrete(guide = guide_legend(title = 'Headgroup')) +
    facet_grid(~ uhgroup, scales = 'free_x', space = 'free_x') +
    # scale_size_continuous(guide = guide_legend(title = 'Count')) +
    scale_shape_manual(values = c(65, 13, 69),
                    labels = as.vector(screen_lab),
                    guide = guide_legend(title = 'Screening')) +
    xlab('Headgroup, carbon count & unsaturation') +
    ylab('Protein') +
    ggtitle('Protein vs. headgroup, carbon count & unsaturation') +
    theme_bw() +
    theme(
        text = element_text(family = "DINPro"),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
        axis.text.y = element_text(size = 14, color = protein_colors(adata12cc$protein)),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 21),
        strip.text = element_text(size = 18, angle = 90),
        #panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggsave(paste('gg_', pname, 'hg_cc_protein_screen_point_facet.pdf', sep = ''),
       device = cairo_pdf, width = 64, height = 12, limitsize = FALSE)

cat("\t:: Plotting `Protein vs. headgroup, carbon count & unsaturation by screening` with facets.\n")
p <- ggplot(adata12hgcci, aes(x = cc, y = protein, shape = screen, color = cls, size = log10(isum))) +
    geom_point(alpha = 0.7) +
    #scale_color_discrete(guide = guide_legend(title = 'Headgroup')) +
    facet_grid(~ uhgroup, scales = 'free_x', space = 'free_x') +
    # scale_size_continuous(guide = guide_legend(title = 'Count')) +
    scale_shape_manual(values = c(65, 69),
                    labels = c(A = 'HEK cells', E = 'E. coli & liposomes'),
                    guide = guide_legend(title = 'Screening')) +
    scale_color_manual(values = c('black', 'red'), guide = guide_legend(title = 'Class')) +
    scale_size_continuous(guide = guide_legend(title = 'Intensity (log)')) +
    xlab('Headgroup, carbon count & unsaturation') +
    ylab('Protein') +
    ggtitle('Protein vs. headgroup, carbon count & unsaturation') +
    theme_bw() +
    theme(
        text = element_text(family = "DINPro"),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
        axis.text.y = element_text(size = 14, color = protein_colors(adata12hgcci$protein)),
        #axis.ticks.x = element_blank(),
        plot.title = element_text(size = 21),
        strip.text = element_text(size = 18, angle = 90),
        #panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggsave(paste('gg_', pname, 'hg_cc_protein_screen_int_class_point_facet.pdf', sep = ''),
       device = cairo_pdf, width = 64, height = 12, limitsize = FALSE)
##

cat("\t:: Plotting `Carbon count, unsaturation, intensity and headgroup\nby protein and ion mode`.\n")
p <- ggplot(adata12hgi, aes(x = unsat, y = carb, size = log10(isum), color = uhgroup)) +
    geom_text(aes(label=uhgroup), alpha = 0.35, nudge_x = 0, nudge_y = 0) +
    geom_point(size = 1, alpha = 1.0, shape = 43) +
    scale_x_continuous(breaks = seq(0, max(adata12$unsat[!is.nan(adata12$unsat)]), 2)) +
    scale_color_manual(guide = guide_legend(title = 'Headgroup', title.position = 'top', ncol = 2, order = 2),
                    values = distinct_cols2[1:length(levels(as.factor(adata12hgi$uhgroup)))]) +
    scale_size_continuous(guide = guide_legend(title = 'Intensity (log)', title.position = 'top', order = 1)) +
    facet_grid(protein ~ screen, labeller = labeller(screen = c(A = 'HEK cells', E = 'E.coli & liposomes'))) +
    xlab('Unsaturation') +
    ylab('Carbon count') +
    ggtitle('Carbon count, unsaturation,\nintensity and headgroup\nby protein and ion mode') +
    theme_bw() +
    theme(
        text = element_text(family = "DINPro"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 18),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 9),
        legend.position = 'right',
        strip.background = element_rect(fill="white")
    )

ggsave(paste('gg_', pname, 'cc_unsat_int_protein_hg_mode.pdf', sep = ''),
       device = cairo_pdf, width = 5, height = 64, limitsize = FALSE)
