#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(ggplot2)
require(dplyr)

#
# constants
#

infile_c12 <- 'combined_network_df.csv'
infile_c1  <- 'combined_network_classI_df.csv'
infile_a   <- 'antonella_final.csv'
infile_e   <- 'enric_processed.csv'
gpl <- c('PE', 'PC', 'PG', 'PA', 'PI', 'PIP', 'PS')
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
    
    pcolors <- get_domain_colors(ae)
    
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
            lyso2    = startsWith(as.character(uhgroup), 'Lyso'),
            ether   = endsWith(as.character(uhgroup), '-O'),
            hg0     = gsub('-O$', '', gsub('^Lyso', '', as.character(uhgroup)))
        ) %>%
        group_by(protein, ionm, screen, uhgroup) %>%
        mutate(ihg = sum(as.numeric(intensity))) %>%
        group_by(protein, ionm, screen) %>%
        mutate(itotal = sum(as.numeric(ihg))) %>%
        group_by(protein, ionm, screen, uhgroup) %>%
        mutate(irel = ihg / itotal, scr_ionm = paste0(screen, ionm)) %>%
        summarize_all(first) %>%
        mutate(uhgroup = as.character(uhgroup))
    
    return(aeg)
    
}

do_plot <- function(df, pdfname, title, by_hg = TRUE){
    
    width <- ifelse(by_hg, 12, 6)
    
    dcolors <- c('#3A7AB3', '#608784', '#03928C', '#CF5836', '#7575BE',
                 '#D6708B', '#65B9B9', '#69B3D6', '#C441B3', '#9B998D')
    
    names(dcolors) <- sort(unique(levels(df$domain)))
    
    dprotein <- list()
    for(i in 1:dim(df)[1]){
        dprotein[[as.character(df$protein[i])]] <- as.character(df$domain[i])
        }
    
    pcolors <- as.character(dcolors[as.character(dprotein)])
    names(pcolors) <- names(dprotein)
    
    p <- ggplot(df,
        aes(y = protein, x = lysoether, group = datasource,
            shape = datasource)) +
        geom_point(position = position_dodge(width = 0.6))
    #     scale_color_manual(values = c('#1B7D8D', '#F9382D'),
    #                        labels = c('Literature', 'Our experiments'),
    #                        guide = guide_legend(title = 'Data source'))
    
    if(by_hg){
        
        p <- p + facet_grid(.~hg0)
        
    }
    
    p <- p +
        scale_x_discrete(drop = FALSE) +
        scale_shape_manual(
            values = c(17, 16, 18, 8),
            labels = c(
                'MSA' = 'MS, HEK cells\n(in vivo)',
                'MSE' = 'MS, E. coli & liposomes\n(in vitro)'
            ),
            guide = guide_legend(title = ''),
            drop = FALSE
        ) +
        xlab('Lipids') +
        ylab('Lipid transfer proteins') +
        ggtitle(title) +
        theme_bw() +
        theme(
            text = element_text(family = "DINPro"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
            axis.text.y = element_text(color = pcolors[levels(df$protein)])
        )

    ggsave(pdfname, device = cairo_pdf, width = width, height = 6)
}


data_1     <- preprocess(infile_c1,  by_hg = FALSE)
data_12    <- preprocess(infile_c12, by_hg = FALSE)
data_1_hg  <- preprocess(infile_c1,  by_hg = TRUE)
data_12_hg <- preprocess(infile_c12, by_hg = TRUE)

do_plot(data_1,  'gg_combined_cargo_classI_lyso_ether_1.pdf',
        'Lyso and ether species: only class I', by_hg = FALSE)
do_plot(data_12, 'gg_combined_cargo_lyso_ether_1.pdf',
        'Lyso and ether species: with class II from in vivo', by_hg = FALSE)
do_plot(data_1_hg,  'gg_combined_cargo_classI_lyso_ether.pdf',
        'Lyso and ether species: only class I', by_hg = TRUE)
do_plot(data_12_hg, 'gg_combined_cargo_lyso_ether.pdf',
        'Lyso and ether species: with class II from in vivo', by_hg = TRUE)
