#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(ggplot2)
require(dplyr)
require(lazyeval)
require(stringr)
require(ggrepel)

#
# constants
#

infile_c12 <- 'combined_network_df.csv'
infile_c1  <- 'combined_network_classI_df.csv'
infile_a   <- 'antonella_final.csv'
infile_e   <- 'enric_processed.csv'
gpl <- c('PE', 'PC', 'PG', 'PA', 'PI', 'PS', 'DAG', 'FA')
mss <- c('MSA', 'MSE')

#
# functions
#

get_domain_colors <- function(df){
    
    dcolors <- c('#3A7AB3', '#608784', '#03928C', '#CF5836',  '#7575BE',
                 '#D6708B', '#65B9B9', '#69B3D6', '#C441B3', '#9B998D')
    
    names(dcolors) <- sort(unique(levels(df$domain)))
    
    dprotein <- list()
    
    for(i in 1:dim(df)[1]){
        
        dprotein[[as.character(df$protein[i])]] <- as.character(df$domain[i])
        
    }
    
    pcolors <- as.character(dcolors[as.character(dprotein)])
    names(pcolors) <- names(dprotein)
    
    return(pcolors)
    
}


preprocess <- function(infile, by_hg = TRUE){
    
    data <- read.csv(infile, header = TRUE, sep = '\t')

    # order by domain and protein name
    data <- data %>%
        filter(datasource %in% mss) %>%
        mutate(gsub('MAG', 'LysoDAG', uhgroup)) %>%
        arrange(domain, protein) %>%
        mutate(
            protein = factor(protein, unique(protein)),
            lyso    = startsWith(as.character(other), 'Lyso'),
            ether   = endsWith(as.character(other), '-O'),
            hg0     = gsub('-O$', '', gsub('^Lyso', '', as.character(other)))
        ) %>%
        mutate(
            lysoether = factor(ifelse(lyso & ether, 'Lyso & ether',
                                      ifelse(lyso, 'Lyso',
                                             ifelse(ether, 'Ether', '2x Ester')
                                    )
                            ),
                            levels = c('2x Ester', 'Ether', 'Lyso', 'Lyso & ether')
                        ),
            gpl       = hg0 %in% gpl
        )

    #

    data$reltpsa <- data$tpsa / data$heavyatoms
    data$charge <- factor(data$charge,
                          levels = c('----', '---', '--', '-', '+--', '+-', '0'),
                          ordered = TRUE)
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

    this_view <- data %>%
        filter(gpl) %>%
        droplevels()
    
    if(by_hg){
        this_view <- this_view %>%
            group_by(protein, datasource, lysoether, hg0) %>%
            summarise_all(first)
    }else{
        this_view <- this_view %>%
            group_by(protein, datasource, lysoether) %>%
            summarise_all(first)
    }
    
    return(this_view)

}

preprocess2 <- function(by_hg = TRUE, only_class1 = TRUE){
    
    a  <- read.table(infile_a, sep = '\t', header = TRUE)
    e  <- read.table(infile_e, sep = '\t', header = TRUE)
    ae <- rbind(a, e)
    
    ae <- ae %>%
        mutate(gsub('MAG', 'LysoDAG', uhgroup)) %>%
        group_by(protein) %>%
        mutate(
            screens = paste0(unique(sort(screen)), collapse = ''),
            protein0 = protein
        ) %>%
        ungroup() %>%
        mutate(
            protein = ifelse(
                screens == 'AE',
                paste('}{', protein),
                ifelse(screens == 'A',
                    paste('{', protein),
                    paste('}', protein)
                )
            )
        )
    
    if(only_class1){
    
        aeg <- ae %>%
            filter(cls == 'I' & uhgroup != 'P40')
    
    }else{
        
        aeg <- ae %>%
            filter(cls %in% c('I', 'II') & uhgroup != 'P40')
    }
    
    aeg <- aeg %>%
        mutate(
            datasource = paste0(screen, ionm),
            lyso2    = startsWith(as.character(uhgroup), 'Lyso') | uhgroup == 'FA',
            ether   = endsWith(as.character(uhgroup), '-O'),
            hg0     = gsub('-O$', '', gsub('^Lyso', '', as.character(uhgroup))),
            scr_ionm = paste0(screen, ionm)
        ) %>%
        mutate(
            lysoether =
                factor(
                    ifelse(
                        lyso2 & ether,
                        'Lyso & Ether',
                        ifelse(
                            lyso2,
                            'Lyso',
                            ifelse(
                                ether,
                                'Ether', '2x Ester'
                            )
                        )
                    ),
                    levels = c('2x Ester', 'Ether', 'Lyso', 'Lyso & Ether')
                ),
            gpl = hg0 %in% gpl
        ) %>%
        group_by(protein, ionm, screen, uhgroup, cls) %>%
        mutate(ihg = sum(as.numeric(intensity))) %>%
        ungroup() %>%
        group_by(protein, ionm, screen) %>%
        mutate(itotal = sum(as.numeric(ihg))) %>%
        ungroup() %>%
        group_by(protein, ionm, screen, uhgroup, cls) %>%
        mutate(irel = ihg / itotal) %>%
        summarize_all(first) %>%
        ungroup() %>%
        filter(gpl) %>%
        ungroup() %>%
        #group_by(protein, scr_ionm, hg0, lysoether) %>%
        #summarise_all(first) %>%
        arrange(domain, sapply(strsplit(protein, ' '), function(x){last(x)})) %>%
        mutate(
            protein = factor(protein, unique(protein)),
            datasource = factor(datasource)
        )
    
    if(by_hg){
        
        aeg <- aeg %>%
            group_by(protein, datasource, lysoether, hg0) %>%
            summarise_all(first) %>%
            ungroup()
        
    }else{
        
        aeg <- aeg %>%
            group_by(protein, datasource, lysoether) %>%
            summarise_all(first) %>%
            ungroup()
        }
    
    return(aeg)
    
}

intensity_bubble <- function(df, pdfname, title, by_hg = TRUE, intensities = FALSE){
    
    width <- ifelse(by_hg, 12, 6)
    
    pcolors <- get_domain_colors(df)
    
    if(intensities){
        
        p <- ggplot(
                df,
                aes(
                    y = protein,
                    x = lysoether,
                    group = datasource,
                    size  = irel,
                    color = datasource,
                    shape = cls
                )
            ) +
            geom_point(position = position_dodge(width = 0.6), alpha = 0.7) +
            scale_color_manual(
                values = c(
                    'Apos' = '#1B7D8D',
                    'Aneg' = '#1B618D',
                    'Epos' = '#F9382D',
                    'Eneg' = '#F96B2D'
                ),
                labels = c(
                    'Apos' = 'HEK cells, + ion mode',
                    'Aneg' = 'HEK cells, - ion mode',
                    'Epos' = 'E. coli & liposomes,\n+ ion mode',
                    'Eneg' = 'E. coli & liposomes,\n- ion mode'
                ),
                guide = guide_legend(title = 'Data source\nand ion mode')
            ) +
            scale_size(
                limits = c(0.0, 1.0),
                guide = guide_legend(title = 'Relative intensity')
            ) +
            scale_shape_manual(
                values = c(
                    'I'  = 16,
                    'II' = 10
                ),
                labels = c(
                    'I'  = 'Class I',
                    'II' = 'Class II'
                ),
                guide = guide_legend(title = 'Result class')
            )
        
    }else{
        
        p <- ggplot(
                df,
                aes(
                    y = protein,
                    x = lysoether,
                    group = datasource,
                    shape = datasource
                )
            ) +
            geom_point(position = position_dodge(width = 0.6)) +
        #     scale_color_manual(values = c('#1B7D8D', '#F9382D'),
        #                        labels = c('Literature', 'Our experiments'),
        #                        guide = guide_legend(title = 'Data source')) +
            scale_shape_manual(
                values = c(17, 16, 18, 8),
                labels = c(
                    'MSA' = 'MS, HEK cells\n(in vivo)',
                    'MSE' = 'MS, E. coli & liposomes\n(in vitro)'
                    ),
                guide = guide_legend(title = ''),
                drop = FALSE
            )
        
    }
    
    if(by_hg){
        
        p <- p + facet_grid(.~hg0)
        
    }
    
    p <- p +
        scale_x_discrete(drop = FALSE) +
        xlab('Lipids') +
        ylab('Lipid transfer proteins') +
        ggtitle(title) +
        theme_bw() +
        theme(
            text = element_text(family = "DINPro"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
            axis.text.y = element_text(color = pcolors[levels(df$protein)])
        )
        # +
        #annotate(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, fill = '#CCCCCC', geom = 'rect', alpha = .2)+
        #annotate(xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf, fill = '#CCCCCC', geom = 'rect', alpha = .2)

    ggsave(pdfname, device = cairo_pdf, width = width, height = 8)
}

align_by_hg0 <- function(d){
    
    d0 <- d %>%
        select(
            protein0, datasource, lysoether, lyso, lyso2, ether,
            screen, screens, cls, ionm, uhgroup, hg0, irel
        ) %>%
        droplevels() %>%
        group_by(protein0, datasource, cls, uhgroup) %>%
        mutate(clsh = first(sort(cls))) %>%
        summarise_all(first) %>%
        ungroup()
    
    lyso  <- d0 %>% filter(lyso2 & !ether)
    ether <- d0 %>% filter(ether & !lyso2)
    lyset <- d0 %>% filter(ether & lyso2)
    ester <- d0 %>% filter(!lyso2 & !ether)
    
    joincols <- c('protein0', 'datasource', 'hg0')
    
    full  <- ester %>%
        full_join(
            ether,
            by = joincols,
            suffix = (c('.ester', '.ether'))
        ) %>%
        full_join(
            lyso,
            by = joincols,
            suffix = (c('', '.lyso'))
        ) %>%
        full_join(
            lyset,
            by = joincols,
            suffix = (c('.lyso', '.lyset'))
        ) %>%
        mutate(
            irel.ester = ifelse(is.na(irel.ester), 0.0, irel.ester),
            irel.ether = ifelse(is.na(irel.ether), 0.0, irel.ether),
            irel.lyso  = ifelse(is.na(irel.lyso ), 0.0, irel.lyso ),
            irel.lyset = ifelse(is.na(irel.lyset), 0.0, irel.lyset)
        ) %>%
        mutate(
            ifc.ether  = ifelse(irel.ester == 0 & irel.ether == 0, 0.0, (irel.ether / (irel.ester + irel.ether) - 0.5) * 2.0),
            ifc.lyso   = ifelse(irel.ester == 0 & irel.lyso  == 0, 0.0, (irel.lyso  / (irel.ester + irel.lyso ) - 0.5) * 2.0),
            ifc.lyset  = ifelse(irel.ester == 0 & irel.lyset == 0, 0.0, (irel.lyset / (irel.ester + irel.lyset) - 0.5) * 2.0)
        )
    
    return(full)
    
}

pref_bar <- function(data, y){
    
    label <- strsplit(deparse(substitute(y)), '\\.')[[1]][2]
    pdfname <- sprintf('pref-%s.pdf', label)
    
    p <- ggplot(data, aes_(x = quote(datasource), substitute(y), fill = substitute(y) < 0)) +
        geom_col() +
        facet_grid(protein0 ~ hg0) +
        scale_fill_manual(
            values = c(
                'FALSE' = '#1B7D8D',
                'TRUE'  = '#F9382D'
            ),
            labels = c(
                'FALSE' = sprintf('Preferring %s', label),
                'TRUE'  = 'Preferring 2x ester'
            ),
            guide = guide_legend(title = 'Preference')
        ) +
        ggtitle(sprintf('Preference towards %s species', label)) +
        xlab('Lipid') +
        ylab('Protein') +
        theme_bw() +
        theme(
            text = element_text(family = "DINPro"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
        )
    
    ggsave(pdfname, device = cairo_pdf, width = 6, height = 24)
    
}

pref_point <- function(data, y, labsize = 2, labangle = 45){
    
    label <- strsplit(deparse(substitute(y)), '\\.')[[1]][2]
    pdfname <- sprintf('pref-dot-%s.pdf', label)
    
    p <- ggplot(data,
            aes_(
                x = quote(irel.ester + 1e-9),
                y = substitute(y + 1e-9),
                color = quote(hg0),
                label = quote(protein0)
            )
        ) +
        scale_x_log10() +
        scale_y_log10() +
        scale_color_discrete(guide = guide_legend(title = 'Headgroup')) +
        geom_text(family = 'DINPro', angle = labangle, size = labsize) +
        ggtitle(sprintf('Preference towards %s species', label)) +
        xlab('2x ester relative intensity (log)') +
        ylab(sprintf('%s relative intensity (log)', str_to_title(label))) +
        theme_bw() +
        theme(
            text = element_text(family = "DINPro")
        )
    
    ggsave(pdfname, device = cairo_pdf, width = 8, height = 8)
    
}

# data_1     <- preprocess(infile_c1,  by_hg = FALSE)
# data_12    <- preprocess(infile_c12, by_hg = FALSE)
# data_1_hg  <- preprocess(infile_c1,  by_hg = TRUE)
# data_12_hg <- preprocess(infile_c12, by_hg = TRUE)
# 
# intensity_bubble(data_1,  'gg_combined_cargo_classI_lyso_ether_1.pdf',
#         'Lyso and ether species: only class I', by_hg = FALSE)
# intensity_bubble(data_12, 'gg_combined_cargo_lyso_ether_1.pdf',
#         'Lyso and ether species: with class II from in vivo', by_hg = FALSE)
# intensity_bubble(data_1_hg,  'gg_combined_cargo_classI_lyso_ether.pdf',
#         'Lyso and ether species: only class I', by_hg = TRUE)
# intensity_bubble(data_12_hg, 'gg_combined_cargo_lyso_ether.pdf',
#         'Lyso and ether species: with class II from in vivo', by_hg = TRUE)

data_1     <- preprocess2(by_hg = FALSE)
data_12    <- preprocess2(by_hg = FALSE, only_class1 = FALSE)
data_1_hg  <- preprocess2(by_hg = TRUE)
data_12_hg <- preprocess2(by_hg = TRUE, only_class1 = FALSE)

intensity_bubble(data_1,
        'gg_combined_cargo_classI_lyso_ether_int_1.pdf',
        'Lyso and ether species: only class I',
        by_hg = FALSE, intensities = TRUE)

intensity_bubble(data_12,
        'gg_combined_cargo_lyso_ether_int_1.pdf',
        'Lyso and ether species: with class II',
        by_hg = FALSE, intensities = TRUE)

intensity_bubble(data_1_hg,
        'gg_combined_cargo_classI_lyso_ether_int.pdf',
        'Lyso and ether species: only class I',
        by_hg = TRUE, intensities = TRUE)

intensity_bubble(data_12_hg,
        'gg_combined_cargo_lyso_ether_int.pdf',
        'Lyso and ether species: with class II',
        by_hg = TRUE, intensities = TRUE)

full <- align_by_hg0(data_12_hg)
pref_bar(full, irel.ether)
pref_bar(full, irel.lyso)

pref_point(full, irel.ether)
pref_point(full, irel.lyso)
