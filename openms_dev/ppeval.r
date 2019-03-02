require(ggplot2)
require(viridis)
require(dplyr)
require(readr)

infiles <- c(
    'STARD10_ppeval__iterative__no-smoothing.tsv',
    'STARD10_ppeval__iterative__gaussian-20.tsv',
    'STARD10_ppeval__iterative__gaussian-10.tsv'
)


for(infile in infiles){
    
    all_plots(infile)
    
}

ppm <- function(theoretical, measured){
    
    (measured - theoretical) / theoretical * 1e6
    
}


all_plots <- function(infile){
    
    label <- gsub('.*__(.*)\\.tsv', '\\1', infile)
    
    d <- suppressMessages(read_tsv(infile)) %>%
        mutate(
            ppm_feature_peaks = ppm(mz_theoretical, mz_feature_peaks),
            ppm_feature_progenesis = ppm(mz_theoretical, mz_feature_progenesis),
            ppm_feature_oms = ppm(mz_theoretical, mz_feature_oms),
            ppm_scan_oms = ppm(mz_theoretical, mz_scan_oms)
        )

    features <- unique(d$feature_id_oms)

    for(ife in features){
        
        this_data <- d %>% filter(feature_id_oms == ife)
        
        mz_theoretical <- this_data$mz_theoretical[1]
        
        peaks_ppm <- this_data$ppm_feature_peaks[1]
        progenesis_ppm <- this_data$ppm_feature_progenesis[1]
        oms_ppm <- this_data$ppm_feature_oms[1]
        
        p <- ggplot(
                this_data %>% filter(isotope == 0),
                aes(x = rt_scan_oms, y = ppm_scan_oms)
            ) +
            geom_line(color = '#EF3A43', lwd = .3) +
            geom_point(size = .3, color = '#EF3A43') +
            # PEAKS A10
            geom_hline(yintercept = peaks_ppm, color = '#EC6224', lwd = .2) +
            # Progenesis A10
            geom_hline(yintercept = progenesis_ppm, color = '#4C4B6B', lwd = .2) +
            # OpenMS A10
            geom_hline(yintercept = oms_ppm, color = '#B7CA54', lwd = .2) +
            # instrument error
            geom_hline(
                yintercept = -1.778,
                color = '#529986',
                lwd = .2,
                linetype = 'dashed',
            ) +
            xlab('RT [s]') +
            ylim(c(NA, 0)) +
            ylab('Error [ppm]') +
            ggtitle(sprintf(
                '[%s+H]+ (%.04f) %s',
                this_data$lipid_species[1],
                mz_theoretical,
                label
            )) +
            theme_linedraw()
        
        ggsave(
            sprintf(
                'img/STARD10_A10_pos_%.1f_ppms_%s.pdf',
                mz_theoretical,
                label
            ),
            device = cairo_pdf,
            height = 2.5,
            width = 5
        )
        
        p <- ggplot(
                this_data,
                aes(x = rt_scan_oms, y = ppm_scan_oms)
            ) +
            facet_grid(rows = vars(isotope)) +
            geom_line(color = '#EF3A43', lwd = .3) +
            geom_point(size = .3, color = '#EF3A43') +
            # PEAKS A10
            geom_hline(yintercept = peaks_ppm, color = '#EC6224', lwd = .2) +
            # Progenesis A10
            geom_hline(yintercept = progenesis_ppm, color = '#4C4B6B', lwd = .2) +
            # OpenMS A10
            geom_hline(yintercept = oms_ppm, color = '#B7CA54', lwd = .2) +
            # instrument error
            geom_hline(
                yintercept = -1.778,
                color = '#529986',
                lwd = .2,
                linetype = 'dashed',
            ) +
            xlab('RT [s]') +
            ylim(c(NA, 0)) +
            ylab('Error [ppm]') +
            ggtitle(sprintf(
                '[%s+H]+ (%.04f) %s',
                this_data$lipid_species[1],
                mz_theoretical,
                label
            )) +
            theme_linedraw()
        
        ggsave(
            sprintf(
                'img/STARD10_A10_pos_%.1f_ppms_%s_isotopes.pdf',
                mz_theoretical,
                label
            ),
            device = cairo_pdf,
            height = 5,
            width = 5
        )
        
        
    }
    
}
