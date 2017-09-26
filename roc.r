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

infile_a   <- 'antonella_final.csv'
infile_e   <- 'enric_processed.csv'
infile_lit_ligands <- 'binding_properties_plain.csv'

#
# functions
#

preprocess_minimal <- function(){
    
    a  <- read.table(infile_a, sep = '\t', header = TRUE)
    e  <- read.table(infile_e, sep = '\t', header = TRUE)
    ae <- rbind(a, e) %>% mutate(lit = lit == 'True')
    
    return(ae)

}

read_literature_ligands <- function(){
    
    d <- read.csv(infile_lit_ligands, header = FALSE, sep = '\t')
    names(d) <- c('carrier', 'ligand')
    
    litlig <- sapply(levels(d$carrier),
                     function(cr){as.character(d$ligand[d$carrier == cr])})
    names(litlig) <- levels(d$carrier)
    
    return(litlig)
    
}

sens_spec <- function(pos, neg){
    
    tp <- dim(pos %>% filter( lit))[1]
    fp <- dim(pos %>% filter(!lit))[1]
    
    fn <- dim(neg %>% filter( lit))[1]
    tn <- dim(neg %>% filter(!lit))[1]
    
    return(list(
        sens = tp / (tp + fn),
        spec = tn / (tn + fp),
        prec = tp / (tp + fp),
        neg  = tn + fn,
        pos  = tp + fp,
        fpr  = fp / (fp + tn),
        fnr  = fn / (tp + fn),
        fdr  = fp / (tp + fp),
        tp   = tp,
        fp   = fp,
        tn   = tn,
        fn   = fn,
        n    = tp + fp + tn + fn
    ))
    
}

roc_conditions <- function(){
    
    in_screen <- function(d, sc, cl){
        
        return(
            d %>%
            filter(cls == cl & screen == sc) %>%
            select(protein, screen, uhgroup) %>%
            group_by(protein, screen, uhgroup) %>%
            summarise_all(first) %>%
            ungroup() %>%
            mutate(inn = TRUE)
        )
        
    }
    
    by_cols <- c('protein', 'screen', 'uhgroup')
    
    d <- preprocess_minimal() %>%
        filter(cls %in% c('I', 'II')) %>%
        group_by(protein, ionm, screen, id) %>%
        summarise_all(first) %>%
        ungroup()
    
    inAI  <- in_screen(d, 'A', 'I')  %>% rename(in_ai  = inn)
    inAII <- in_screen(d, 'A', 'II') %>% rename(in_aii = inn)
    inEI  <- in_screen(d, 'E', 'I')  %>% rename(in_ei  = inn)
    inEII <- in_screen(d, 'E', 'II') %>% rename(in_eii = inn)
    
    d <- d %>%
        left_join(inAI,  by = by_cols) %>%
        left_join(inAII, by = by_cols) %>%
        left_join(inEI,  by = by_cols) %>%
        left_join(inEII, by = by_cols) %>%
        mutate(
            in_ai  = !is.na(in_ai ),
            in_aii = !is.na(in_aii),
            in_ei  = !is.na(in_ei ),
            in_eii = !is.na(in_eii)
        )
    
    return(list(
        list(
            name = 'Only class I',
            vals = sens_spec(
                d %>% filter(cls == 'I'),
                d %>% filter(cls != 'I')
            )
        ),
        list(
            name = 'In vivo, class I',
            vals = sens_spec(
                d %>% filter(cls == 'I' & screen == 'A'),
                d %>% filter(cls != 'I' | screen != 'A')
            )
        ),
        list(
            name = 'In vitro, class I',
            vals = sens_spec(
                d %>% filter(cls == 'I' & screen == 'E'),
                d %>% filter(cls != 'I' | screen != 'E')
            )
        ),
        list(
            name = 'In vitro, class I & II',
            vals = sens_spec(
                d %>% filter(screen == 'E'),
                d %>% filter(screen != 'E')
            )
        ),
        list(
            name = 'In vivo, class I & II',
            vals = sens_spec(
                d %>% filter(screen == 'A'),
                d %>% filter(screen != 'A')
            )
        ),
        list(
            name = 'Class I and in vivo class II',
            vals = sens_spec(
                d %>% filter(cls == 'I' | (screen == 'A' & cls =='II')),
                d %>% filter(cls == 'II' & screen != 'A')
            )
        ),
        list(
            name = 'Class I & in vivo class II confirmed by in vitro',
            vals = sens_spec(
                d %>% filter(cls == 'I'  | (screen == 'A' & (in_ei | in_eii))),
                d %>% filter(cls == 'II' & (screen != 'A' | (!in_ei & !in_eii)))
            )
        ),
        list(
            name = 'Class I & class II confirmed by other screen',
            vals = sens_spec(
                d %>% filter(
                    cls == 'I'  |
                    (screen == 'A' & (in_ei | in_eii)) |
                    (screen == 'E' & (in_ai | in_aii))
                ),
                d %>% filter(
                    cls != 'I'  & (
                        (screen == 'A' & (!in_ei & !in_eii)) |
                        (screen == 'E' & (!in_ai & !in_aii))
                    )
                )
            )
        )
    ))
    
}

roc_df <- function(){
    
    roc <- roc_conditions()
    
    rocdf <- do.call(rbind, lapply(roc, function(r){data.frame(c(list(name = r$name), r$vals))}))
    
    return(as.data.frame(rocdf))
    
}

roc_plot <- function(){
    
    rocdf <- roc_df()
    
    rocdf <- rocdf[1:6,]
    rocdf <- rocdf %>% filter(sens > 0.0)
    
    p <- ggplot(rocdf, aes(x = sens, y = 1 - spec, label = name)) +
        geom_abline(intercept = 0, slope = 1, color = 'red') +
        geom_point() +
        geom_text_repel(family = 'DINPro') +
        ggtitle('Various conditions in ROC space') +
        xlab('1 - specificity') +
        ylab('Sensitivity') +
        xlim(0.0, 1.0) +
        ylim(0.0, 1.0) +
        theme_bw() +
        theme(
            text = element_text(family = 'DINPro')
        )
    
    ggsave('roc_classes.pdf', device = cairo_pdf, width = 6, height = 6)
    
}
