#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(ggplot2)
require(dplyr)
require(lazyeval)
#require(RColorBrewer)

domain_order <- function(df){
    
    df$protein <- factor(df$protein,
                    levels = unique(df[order(as.character(df$domain),
                                             as.character(df$protein)), 'protein']),
                    ordered = TRUE)
    return(df)
}

#
# constants
#

# fname.e <- 'combined_network_df.csv'
fname.e <- 'combined_network_classI_df.csv'

data <- read.csv(fname.e, header = TRUE, sep = '\t')

#

data$protein <- factor(data$protein,
                       levels = unique(data[order(as.character(data$domain), as.character(data$protein)), 'protein']),
                       ordered = TRUE)

data$reltpsa <- data$tpsa / data$heavyatoms
data$charge <- factor(data$charge, levels = c('----', '---', '--', '-', '+--', '+-', '0'), ordered = TRUE)
data$charge[is.na(data$charge)] <- '0'

screensym <- list(
    MSA = '{',
    MSE = '}',
    LITB = '@',
    LIMAL = '#',
    HPTLC = '%'
)

screens <- c('MSA', 'MSE', 'LITB', 'LIMAL', 'HPTLC')

for(sc in screens){
    inscreen <- unique(as.character((data %>% filter(datasource == sc))$protein))
    for(pr in inscreen){
        levels(data$protein)[levels(data$protein) == pr] <- paste(screensym[[sc]], pr)
    }
}

dcolors <- c('#3A7AB3', '#608784', '#03928C', '#CF5836', '#7575BE',
             '#D6708B', '#65B9B9', '#69B3D6', '#C441B3', '#9B998D')

names(dcolors) <- sort(unique(levels(data$domain)))

dprotein <- list()
for(i in 1:dim(data)[1]){
    dprotein[[as.character(data$protein[i])]] <- as.character(data$domain[i])
}

pcolors <- as.character(dcolors[as.character(dprotein)])
names(pcolors) <- names(dprotein)


# compare cargo screenings and literature
this_view <- data %>%
    filter(etype == 'cargo_binding') %>%
    group_by(protein, other, datasource) %>%
    summarise_all(first) %>%
    as.data.frame()

this_view$protein <- factor(this_view$protein,
                            levels = unique(this_view[order(as.character(this_view$domain),
                                                            as.character(this_view$protein)), 'protein']),
                            ordered = TRUE)

p <- ggplot(this_view,
       aes(y = protein, x = other, color = category, group = datasource,
           shape = datasource)) +
    geom_point(position = position_dodge(width = 0.8)) +
    annotate(xmin = 46.6, xmax = 51.4, ymin = -Inf, ymax = Inf, fill = '#CCCCCC', geom = 'rect', alpha = 0.2) +
    scale_color_manual(values = c('#1B7D8D', '#F9382D'),
                       labels = c('Literature', 'Our experiments'),
                       guide = guide_legend(title = 'Data source')) +
    scale_shape_manual(
        values = c(17, 16, 18, 8),
        labels = c('HPTCL' = 'HPTLC', 'MSA' = 'MS, HEK cells', 'MSE' = 'MS, E. coli & liposomes', 'LITB' = 'Literature'),
        guide = guide_legend(title = '')) +
    xlab('Lipids') +
    ylab('Lipid transfer proteins') +
    ggtitle('Cargo binding specificities: only class I') +
    theme_bw() +
    theme(
        text = element_text(family = "DINPro"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
        axis.text.y = element_text(color = pcolors[levels(this_view$protein)])
    )

ggsave('gg_combined_cargo_classI_our-lit.pdf', device = cairo_pdf, width = 12, height = 12)

this_view <- this_view %>%
    filter(datasource %in% c('MSA', 'MSE', 'HPTLC'))

p <- ggplot(this_view,
       aes(y = protein, x = other, color = category, group = datasource,
           shape = datasource)) +
    geom_point(position = position_dodge(width = 0.8)) +
    #annotate(xmin = 50.6, xmax = 56.4, ymin = -Inf, ymax = Inf, fill = '#CCCCCC', geom = 'rect', alpha = 0.2) +
    scale_color_manual(values = c('#1B7D8D', '#F9382D'),
                       labels = c('Our experiments'),
                       guide = guide_legend(title = 'Data source')) +
    scale_shape_manual(
        values = c(17, 16, 18, 8),
        labels = c('HPTCL' = 'HPTLC', 'MSA' = 'MS, HEK cells', 'MSE' = 'MS, E. coli & liposomes', 'LITB' = 'Literature'),
        guide = guide_legend(title = '')) +
    xlab('Lipids') +
    ylab('Lipid transfer proteins') +
    ggtitle('Cargo binding specificities: our screenings, only class I') +
    theme_bw() +
    theme(
        text = element_text(family = "DINPro"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
        axis.text.y = element_text(color = pcolors[levels(this_view$protein)])
    )

ggsave('gg_combined_cargo_classI_our.pdf', device = cairo_pdf, width = 12, height = 12)


# membrane lipid affinities vs. cargo
this_view <- data %>%
    filter(etype == 'cargo_binding' | (etype == 'membrane_affinity' & otype == 'lipid')) %>%
    group_by(protein, other, etype, category) %>%
    summarise_all(first) %>%
    arrange(domain) %>%
    as.data.frame()

this_view$protein <- factor(this_view$protein,
                       levels = unique(this_view[order(as.character(this_view$domain),
                                                       as.character(this_view$protein)), 'protein']),
                       ordered = TRUE)

p <- ggplot() +
    geom_point(
        data = this_view,
        aes(
            y = protein,
            x = other,
            color = etype,
            group = datasource,
            shape = category
        ),
        stat = 'identity',
        position = position_dodge(width = 0.75)
    ) +
    annotate(xmin = 51.6, xmax = 61.4, ymin = -Inf, ymax = Inf, fill = '#CCCCCC', geom = 'rect', alpha = 0.2) +
        #data = data.frame(xmin = c(47), xmax = c(57), ymin = c(1), ymax = c(length(unique(this_view$protein)))),
        #        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        #        fill = 'gray80') +
    scale_color_manual(values = c('cargo_binding' = '#1B7D8D', 'membrane_affinity' = '#F9382D'),
                       labels = c('cargo_binding' = 'Cargo binding', 'membrane_affinity' = 'Membrane affinity'),
                        guide = guide_legend(title = 'Interaction type')) +
    scale_shape_manual(
        values = c('literature' = 18, 'own_experiment' = 8),
        labels = c('literature' = 'Literature', 'own_experiment' = 'Our experiments'),
        guide = guide_legend(title = 'Data source')) +
    #facet_grid(domain ~ ., scales = 'free', space = 'free') +
    xlab('Lipids') +
    ylab('Lipid transfer proteins') +
    ggtitle('Cargo binding vs. membrane affinity specificities in our screenings and the literature') +
    theme_bw() +
    theme(
        text = element_text(family = "DINPro"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
        axis.text.y = element_text(color = pcolors[levels(this_view$protein)])
    )

ggsave('gg_combined_cargo-affy_our-lit.pdf', device = cairo_pdf, width = 12, height = 12)

### charge

this_view <- data %>%
    filter(etype == 'cargo_binding') %>%
    group_by(protein, other) %>%
    summarise_all(first) %>%
    group_by(protein, charge) %>%
    mutate(cnt = n()) %>%
    as.data.frame()

this_view$protein <- factor(this_view$protein,
                            levels = unique(this_view[order(as.character(this_view$domain),
                                                            as.character(this_view$protein)), 'protein']),
                            ordered = TRUE)

p <- ggplot(this_view,
       aes(x = protein, y = charge, size = cnt)) +
    geom_point() +
    xlab('Lipid transfer proteins') +
    ylab('Charges') +
    ggtitle('Charge states of cargos') +
    scale_size_area(guide = guide_legend(title = 'Count of\nlipid classes\n(e.g. PC = 1 class)')) +
    theme_bw() +
    theme(
        text = element_text(family = "DINPro"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1,
                                   color = pcolors[levels(this_view$protein)]),
    )

ggsave('gg_combined_cargo_charge_our-lit.pdf', device = cairo_pdf, width = 16, height = 6)

#  charge by family

this_view <- data %>%
    filter(etype == 'cargo_binding') %>%
    group_by(domain, other) %>%
    summarise_all(first) %>%
    group_by(domain, charge) %>%
    mutate(cnt = n()) %>%
    as.data.frame()

this_view$domain <- factor(this_view$domain,
                            levels = unique(this_view[order(as.character(this_view$domain)), 'domain']),
                            ordered = TRUE)

p <- ggplot(this_view,
       aes(x = domain, y = charge, size = cnt)) +
    geom_point() +
    xlab('Lipid transfer domains') +
    ylab('Charges') +
    ggtitle('Charge states of cargos') +
    scale_size_area(guide = guide_legend(title = 'Count of\nlipid classes\n(e.g. PC = 1 class)')) +
    theme_bw() +
    theme(
        text = element_text(family = "DINPro"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
    )

ggsave('gg_combined_cargo_charge_by-domain_our-lit.pdf', device = cairo_pdf, width = 9, height = 6)


### BY LIPID PHYS & CHEM PROPERTIES ###

physchemplots <- list(
    list(
        ttl = 'Polar surface area of cargos',
        ylab = parse(text = 'Polar~Surface~Area~~bgroup("[", \u212b^{2}, "]")'),
        name = 'tpsa',
        yvar = as.symbol('tpsa')
    ),
    list(
        ttl = 'Relative polar surface area of cargos',
        ylab = parse(text = 'Relative~Polar~Surface~Area~~bgroup("[", \u212b^{2}/heavy~atoms, "]")'),
        name = 'reltpsa',
        yvar = as.symbol('reltpsa')
        ),
    list(
        ttl = 'Partitioning coefficient of cargos',
        ylab = 'Partitioning coefficient [XlogP]',
        name = 'xlogp',
        yvar = as.symbol('xlogp')
    ),
    list(
        ttl = 'Structural complexity of cargos',
        ylab = 'Complexity (Bertz-Hendrickson-Ihlenfeldt formula)',
        name = 'complexity',
        yvar = as.symbol('complexity')
    ),
    list(
        ttl = 'Heavy atom count of cargos',
        ylab = 'Number of heavy atoms',
        name = 'heavyatoms',
        yvar = as.symbol('heavyatoms')
    ),
    list(
        ttl = 'Hydrogen donors on cargos',
        ylab = 'Number of hydrogen donors',
        name = 'hdonors',
        yvar = as.symbol('hdonor')
        ),
    list(
        ttl = 'Hydrogen acceptors on cargos',
        ylab = 'Number of hydrogen acceptors',
        name = 'hacceptors',
        yvar = as.symbol('haccept')
        )
    #,
    #list(
    #    ttl = 'Net charge of cargos',
    #    ylab = 'Net charge',
    #    name = 'netcharge',
    #    yvar = quote(netcharge)
    #)
)

# compare cargo screenings and literature

for(param in physchemplots){
    
    this_view <- data %>%
        filter(etype == 'cargo_binding') %>%
        group_by(protein, other, datasource) %>%
        summarise_all(first) %>%
        as.data.frame()
    
    ordr <- as.character((this_view %>%
        group_by(protein) %>%
        #mutate(meanval = max(tpsa))
        mutate_(.dots = list(meanval = interp(~ f(as.name(as.character(param$yvar))), f = quote(`max`)))) %>%
        summarise_all(first) %>%
        arrange(meanval) %>%
        mutate(protein = factor(protein, protein)) %>%
        as.data.frame())$protein)
    
    print(ordr)
    
    this_view$protein <- factor(this_view$protein, levels = ordr, ordered = TRUE)
    
    p <- ggplot(this_view,
        aes(y = param$yvar, x = protein, color = category, group = datasource,
            shape = datasource)) +
        geom_point(position = position_dodge(width = 0.8)) +
        #annotate(xmin = 51.6, xmax = 57.4, ymin = -Inf, ymax = Inf, fill = '#CCCCCC', geom = 'rect', alpha = 0.2) +
        scale_color_manual(values = c('#1B7D8D', '#F9382D'),
                        labels = c('Literature', 'Our experiments'),
                        guide = guide_legend(title = 'Data source')) +
        scale_shape_manual(
            values = c(17, 16, 18, 8),
            labels = c('HPTCL' = 'HPTLC', 'MSA' = 'MS, HEK cells', 'MSE' = 'MS, E. coli & liposomes', 'LITB' = 'Literature'),
            guide = guide_legend(title = '')) +
        xlab('Proteins') +
        ylab(param$ylab) +
        ggtitle(param$ttl) +
        theme_bw() +
        theme(
            text = element_text(family = "DINPro"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1,
                                    color = pcolors[levels(this_view$protein)])
        )
    
    ggsave(paste('gg_combined_cargo_', param$name, '_y-asc_our-lit.pdf', sep = ''),
        device = cairo_pdf, width = 18, height = 12)

    this_view$protein <- factor(this_view$protein,
                                levels = unique(this_view[order(as.character(this_view$domain),
                                                                as.character(this_view$protein)), 'protein']),
                                ordered = TRUE)
    
    p <- ggplot(this_view,
        aes(y = this_view[[as.character(param$yvar)]], x = protein, color = category, group = datasource,
            shape = datasource)) +
        geom_point(position = position_dodge(width = 0.8)) +
        #annotate(xmin = 51.6, xmax = 57.4, ymin = -Inf, ymax = Inf, fill = '#CCCCCC', geom = 'rect', alpha = 0.2) +
        scale_color_manual(values = c('#1B7D8D', '#F9382D'),
                        labels = c('Literature', 'Our experiments'),
                        guide = guide_legend(title = 'Data source')) +
        scale_shape_manual(
            values = c(17, 16, 18, 8),
            labels = c('HPTCL' = 'HPTLC', 'MSA' = 'MS, HEK cells', 'MSE' = 'MS, E. coli & liposomes', 'LITB' = 'Literature'),
            guide = guide_legend(title = '')) +
        xlab('Proteins') +
        ylab(param$ylab) +
        ggtitle(param$ttl) +
        theme_bw() +
        theme(
            text = element_text(family = "DINPro"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1,
                                    color = pcolors[levels(this_view$protein)])
        )
    
    ggsave(paste('gg_combined_cargo_', param$name, '_our-lit.pdf', sep = ''),
        device = cairo_pdf, width = 18, height = 12)

    # TPSA by family:

    this_view <- data %>%
        filter(etype == 'cargo_binding') %>%
        group_by(domain, other, datasource) %>%
        summarise_all(first) %>%
        as.data.frame()

    this_view$protein <- factor(this_view$protein,
                                levels = unique(this_view[order(as.character(this_view$domain)), 'domain']),
                                ordered = TRUE)

    p <- ggplot(this_view,
            ) +
    geom_boxplot(aes(y = this_view[[as.character(param$yvar)]], x = domain)) +
    geom_point(aes(y = this_view[[as.character(param$yvar)]], x = domain, color = category, group = datasource,
                   shape = datasource), position = position_dodge(width = 0.5)) +
        #annotate(xmin = 51.6, xmax = 57.4, ymin = -Inf, ymax = Inf, fill = '#CCCCCC', geom = 'rect', alpha = 0.2) +
        scale_color_manual(values = c('#1B7D8D', '#F9382D'),
                        labels = c('Literature', 'Our experiments'),
                        guide = guide_legend(title = 'Data source')) +
        scale_shape_manual(
            values = c(17, 16, 18, 8),
            labels = c('HPTCL' = 'HPTLC', 'MSA' = 'MS, HEK cells', 'MSE' = 'MS, E. coli & liposomes', 'LITB' = 'Literature'),
            guide = guide_legend(title = '')) +
        xlab('Lipid transfer domain') +
        ylab(param$ylab) +
        ggtitle(param$ttl) +
        theme_bw() +
        theme(
            text = element_text(family = "DINPro"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1)
        )

    ggsave(paste('gg_combined_cargo_', param$name, '_by-domain', '_our-lit.pdf', sep = ''),
        device = cairo_pdf, width = 9, height = 12)

}

# XLogP by protein

this_view <- data %>%
    filter(etype == 'cargo_binding') %>%
    group_by(protein, other, datasource) %>%
    summarise_all(first) %>%
    as.data.frame()

this_view$protein <- factor(this_view$protein,
                            levels = unique(this_view[order(as.character(this_view$domain),
                                                            as.character(this_view$protein)), 'protein']),
                            ordered = TRUE)

p <- ggplot(this_view,
       aes(y = xlogp, x = protein, color = category, group = datasource,
           shape = datasource)) +
    geom_point(position = position_dodge(width = 0.8)) +
    #annotate(xmin = 51.6, xmax = 57.4, ymin = -Inf, ymax = Inf, fill = '#CCCCCC', geom = 'rect', alpha = 0.2) +
    scale_color_manual(values = c('#1B7D8D', '#F9382D'),
                       labels = c('Literature', 'Our experiments'),
                       guide = guide_legend(title = 'Data source')) +
    scale_shape_manual(
        values = c(17, 16, 18, 8),
        labels = c('HPTCL' = 'HPTLC', 'MSA' = 'MS, HEK cells', 'MSE' = 'MS, E. coli & liposomes', 'LITB' = 'Literature'),
        guide = guide_legend(title = '')) +
    xlab('Proteins') +
    ylab('Partitioning coefficient [XlogP]') +
    ggtitle('Partitioning coefficient of cargos') +
    theme_bw() +
    theme(
        text = element_text(family = "DINPro"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1,
                                   color = pcolors[levels(this_view$protein)])
    )

ggsave(paste('gg_combined_cargo_', 'xlogp', '_our-lit.pdf', sep = ''),
       device = cairo_pdf, width = 18, height = 12)

# XLogP by family:

this_view <- data %>%
    filter(etype == 'cargo_binding') %>%
    group_by(domain, other, datasource) %>%
    summarise_all(first) %>%
    as.data.frame()

this_view$protein <- factor(this_view$protein,
                            levels = unique(this_view[order(as.character(this_view$domain)), 'domain']),
                            ordered = TRUE)

p <- ggplot(this_view,
       aes(y = xlogp, x = domain, color = category, group = datasource,
           shape = datasource)) +
    geom_point(position = position_dodge(width = 0.5)) +
    #annotate(xmin = 51.6, xmax = 57.4, ymin = -Inf, ymax = Inf, fill = '#CCCCCC', geom = 'rect', alpha = 0.2) +
    scale_color_manual(values = c('#1B7D8D', '#F9382D'),
                       labels = c('Literature', 'Our experiments'),
                       guide = guide_legend(title = 'Data source')) +
    scale_shape_manual(
        values = c(17, 16, 18, 8),
        labels = c('HPTCL' = 'HPTLC', 'MSA' = 'MS, HEK cells', 'MSE' = 'MS, E. coli & liposomes', 'LITB' = 'Literature'),
        guide = guide_legend(title = '')) +
    xlab('Lipid transfer domain') +
    ylab('Partitioning coefficient [XlogP]') +
    ggtitle('Partitioning coefficient of cargos') +
    theme_bw() +
    theme(
        text = element_text(family = "DINPro"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1)
    )

ggsave(paste('gg_combined_cargo_', 'xlogp_by-domain', '_our-lit.pdf', sep = ''),
       device = cairo_pdf, width = 9, height = 12)

# complexity by protein



this_view <- data %>%
    filter(etype == 'cargo_binding') %>%
    group_by(protein, other, datasource) %>%
    summarise_all(first) %>%
    as.data.frame()

this_view$protein <- factor(this_view$protein,
                            levels = unique(this_view[order(as.character(this_view$domain),
                                                            as.character(this_view$protein)), 'protein']),
                            ordered = TRUE)

p <- ggplot(this_view,
       aes(y = xlogp, x = protein, color = category, group = datasource,
           shape = datasource)) +
    geom_point(position = position_dodge(width = 0.8)) +
    #annotate(xmin = 51.6, xmax = 57.4, ymin = -Inf, ymax = Inf, fill = '#CCCCCC', geom = 'rect', alpha = 0.2) +
    scale_color_manual(values = c('#1B7D8D', '#F9382D'),
                       labels = c('Literature', 'Our experiments'),
                       guide = guide_legend(title = 'Data source')) +
    scale_shape_manual(
        values = c(17, 16, 18, 8),
        labels = c('HPTCL' = 'HPTLC', 'MSA' = 'MS, HEK cells', 'MSE' = 'MS, E. coli & liposomes', 'LITB' = 'Literature'),
        guide = guide_legend(title = '')) +
    xlab('Proteins') +
    ylab('Complexity [XlogP]') +
    ggtitle('Partitioning coefficient of cargos') +
    theme_bw() +
    theme(
        text = element_text(family = "DINPro"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1,
                                   color = pcolors[levels(this_view$protein)])
    )

ggsave(paste('gg_combined_cargo_', 'xlogp', '_our-lit.pdf', sep = ''),
       device = cairo_pdf, width = 18, height = 12)

# complexity by family:

this_view <- data %>%
    filter(etype == 'cargo_binding') %>%
    group_by(domain, other, datasource) %>%
    summarise_all(first) %>%
    as.data.frame()

this_view$protein <- factor(this_view$protein,
                            levels = unique(this_view[order(as.character(this_view$domain)), 'domain']),
                            ordered = TRUE)

p <- ggplot(this_view,
       aes(y = xlogp, x = domain, color = category, group = datasource,
           shape = datasource)) +
    geom_point(position = position_dodge(width = 0.5)) +
    #annotate(xmin = 51.6, xmax = 57.4, ymin = -Inf, ymax = Inf, fill = '#CCCCCC', geom = 'rect', alpha = 0.2) +
    scale_color_manual(values = c('#1B7D8D', '#F9382D'),
                       labels = c('Literature', 'Our experiments'),
                       guide = guide_legend(title = 'Data source')) +
    scale_shape_manual(
        values = c(17, 16, 18, 8),
        labels = c('HPTCL' = 'HPTLC', 'MSA' = 'MS, HEK cells', 'MSE' = 'MS, E. coli & liposomes', 'LITB' = 'Literature'),
        guide = guide_legend(title = '')) +
    xlab('Lipid transfer domain') +
    ylab('Partitioning coefficient [XlogP]') +
    ggtitle('Partitioning coefficient of cargos') +
    theme_bw() +
    theme(
        text = element_text(family = "DINPro"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1)
    )

ggsave(paste('gg_combined_cargo_', 'xlogp_by-domain', '_our-lit.pdf', sep = ''),
       device = cairo_pdf, width = 9, height = 12)


###

this_view <- data %>%
    filter(otype == 'location') %>%
    group_by(protein, other, etype, category) %>%
    summarise_all(first) %>%
    arrange(domain) %>%
    as.data.frame()

this_view$protein <- factor(this_view$protein,
                            levels = unique(this_view[order(as.character(this_view$domain),
                                                            as.character(this_view$protein)), 'protein']),
                            ordered = TRUE)

p <- ggplot() +
    geom_point(
        data = this_view,
        aes(
            y = protein,
            x = other,
            color = category,
            group = datasource,
            shape = datasource
            ),
        stat = 'identity',
        position = position_dodge(width = 0.75)
    ) +
    #annotate(xmin = 46.6, xmax = 56.4, ymin = -Inf, ymax = Inf, fill = '#CCCCCC', geom = 'rect', alpha = 0.2) +
    #data = data.frame(xmin = c(47), xmax = c(57), ymin = c(1), ymax = c(length(unique(this_view$protein)))),
    #        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    #        fill = 'gray80') +
    scale_color_manual(
        values = c('own_experiment' = '#F9382D', 'database' = '#1B7D8D'),
        labels = c('own_experiment' = 'Membrane affinity (LiMA screening)',
                   'database' = 'Localization (literature)'),
        guide = guide_legend(title = 'Interaction type')
    ) +
    scale_shape_manual(
        values = c('LIMAM' = 20, 'LIMCM' = 18, 'CCLOC' = 6, 'HPALOC' = 2),
        labels = c('LIMAM' = 'LiMA simple membranes',
                   'LIMCM' = 'LiMA complex membranes',
                   'CCLOC' = 'UniProt localization',
                   'HPALOC' = 'HPA localization'),
        guide = guide_legend(title = 'Data source')
    ) +
    #facet_grid(domain ~ ., scales = 'free', space = 'free') +
    xlab('Localizations') +
    ylab('Lipid transfer proteins') +
    ggtitle('Localization vs. membrane affinity in the literature and the LiMA screening') +
    theme_bw() +
    theme(
        text = element_text(family = "DINPro"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
        axis.text.y = element_text(color = pcolors[levels(this_view$protein)])
    )

ggsave('gg_combined_loc-affy_lit-our.pdf', device = cairo_pdf, width = 12, height = 12)

############################ BY DOMAIN FAMILY ###



############################ THIS IS A DROPPED ATTEMPT ###
if(FALSE){
    
    call_with_args <- function(t, method, args){
        
        #lapply(args, eval)
        if(!is.character(method)){
            method <- deparse(substitute(method))
        }
        if(tail(strsplit(deparse(substitute(group_by)), '')[[1]], 1) != '_'){
            method <- paste(method, '_', sep = '')
        }
        
        return(do.call(method, list(t, .dots = args)))
    }



    plot_generic <- function(data, name1, name2, title, xlab, ylab,
                                group_args,
                                filter_args,
                                aes_args,
                                shape_param,
                                color_param,
                                annot_param = NULL,
                                sum_method = first,
                                typeface = 'DINPro'){
        
        group_args <- append(group_args, ~other)
        group_args <- append(group_args, ~datasource)
        
        # compare cargo screenings and literature
        this_view <- data %>%
            filter_(filter_args) %>%
            call_with_args(group_by, group_args) %>%
            summarise_all(sum_method) %>%
            as.data.frame()

        this_view <- domain_order(this_view)
        
        p <- ggplot(this_view,
            do.call(aes, aes_args)) +
            geom_point(position = position_dodge(width = 0.8)) +
            do.call('scale_color_manual', color_param) +
            do.call('scale_shape_manual', shape_param) +
            xlab(xlab) +
            ylab(ylab) +
            ggtitle(title) +
            theme_bw() +
            theme(
                text = element_text(family = typeface),
                axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
                axis.text.y = element_text(color = pcolors[levels(this_view$protein)])
            )
        
        if(length(annot_param)){
            p <- p + annotate(xmin = annot_param[1], xmax = annot_param[2],
                            ymin = -Inf, ymax = Inf,
                            fill = '#CCCCCC', geom = 'rect', alpha = 0.2)
        }
        
        ggsave(paste('gg_combined_cargo_', name1, '_', name2, '.pdf', sep = ''),
            device = cairo_pdf, width = 12, height = 12)
    }

    plot_generic(data, 'cargo_our', 'by-protein',
                    'Cargo binding specificities',
                    'Lipids',
                    'Lipid transfer proteins',
                    group_args = list(~protein),
                    filter_args = quote(etype == 'cargo_binding'),
                    aes_args = list(y = quote(protein), x = quote(other),
                                    color = quote(category),
                                    group = quote(datasource),
                                    shape = quote(datasource)),
                    shape_param = list(
                            values = c(HPTLC = 17, LITB = 16, MSA = 18, MSE = 8),
                            labels = c(HPTCL = 'HPTLC', MSA = 'MS, HEK cells', MSE = 'MS, E. coli & liposomes', LITB = 'Literature'),
                            guide = guide_legend(title = '')),
                    color_param = list(
                            values = c(literature = '#1B7D8D', own_experiment = '#F9382D'),
                            labels = c(literature = 'Literature', own_experiment = 'Our experiments'),
                            guide = guide_legend(title = 'Data source')),
                    annot_param = c(51.6, 57.4)
                        )

}
