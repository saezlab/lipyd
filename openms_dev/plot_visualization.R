require(ggplot2)
require(viridis)
require(dplyr)
require(readr)


reqplot1 <- ggplot(df %>% filter(V2 > 1 & V5 > 1), aes(x = V2, y = V5))+
  geom_hex() +
  scale_fill_viridis_c() +
  #geom_point(alpha = .1)+
  xlab("RT1")+
  ylab("RT2")+
  scale_size_continuous(name="1/4 mile time")+
  ggtitle("Difference between RT")+
  theme(plot.title = element_text(hjust = 0.5))

ggsave('rt-rt_peaks-vs-openms_STARD10_A10_pos', device = "jpeg", height = 4, width = 5)
print(reqplot1)


V3log10 = log10(df$V3)
V6log10 = log10(df$V6)

df <- df %>%
  mutate(V3log = log10(V3), V6log = log10(V6))

cor.test(df$V3log, df$V6log)

Areaplot <- ggplot(df %>% filter(V3log10 > 1 & V6log10 > 1), aes(x = V3, y = V6))+
  geom_hex()+
  scale_fill_viridis_c() + 
  scale_x_log10() +
  scale_y_log10() +
  #geom_point(alpha = .1)+
  xlab("Log10 Intensity")+
  ylab("Log10 Intensity")+
  scale_size_continuous(name="1/4 mile time")+
  ggtitle("Difference between intensities")+
  theme(plot.title = element_text(hjust = 0.5))

print(Areaplot)
ggsave('Int-Int_peaks-vs-openms_STARD10_A10_pos', device = "jpeg", height = 4, width = 5)

Density <- ggplot(df, aes(x = V7)) + 
  geom_density() +
  #geom_vline(aes(xintercept= 0),
  #color="blue", linetype="dashed", size=1) +
  xlab("ppm")
ggtitle("Ppm density")+
  theme(plot.title = element_text(hjust = 0.5))
print(Density)
ggsave('Density-ppm-openms_STARD10_A10_pos', device = "jpeg", height = 4, width = 5)

print(log10(82900000))


teo <- suppressWarnings(suppressMessages(
    read_table2('teo.txt', col_names = c('rt', 'mz_measured', 'nothing'))
)) %>%
    filter(!is.na(rt)) %>%
    mutate(
        isotope = ifelse(
        mz_measured < 786.7,
        0,
        ifelse(
            mz_measured < 787.7,
            1,
            ifelse(
            mz_measured < 788.7,
            2,
            3
            )
        )
        ),
        mz_theoretical = ifelse(
        mz_measured < 786.7,
        786.60073213211,
        ifelse(
            mz_measured < 787.7,
            787.60939704799,
            ifelse(
            mz_measured < 788.7,
            788.61806196387,
            789.62672687975
            )
        )
        ),
        ppm = (mz_measured - mz_theoretical) / mz_theoretical * 1e06
    )

facet_plot <- ggplot(PPIter3, aes(V5, V6)) +
    geom_line(color = '#EF3A43', lwd = .3) +
    geom_point(size = .3, color = '#EF3A43') +
    # PEAKS A10
    geom_hline(yintercept = -1.6935, color = '#B7CA54', lwd = .1) +
    # instrument error
    geom_hline(yintercept = -1.778, color = '#529986', lwd = .1) +
    facet_grid(rows = vars(V2)) +
    xlab('RT [s]') +
    ylim(c(NA, 0)) +
    ylab('Error (measured vs. theoretical) [ppm]') +
    theme_linedraw()

ggsave(
    'STARD10_A10_pos_786.5_isotopes_ppms.pdf',
    device = cairo_pdf,
    height = 4,
    width = 5
)

facet_plot <- ggplot(teo %>% filter(isotope == 0), aes(rt, ppm)) +
    geom_line(color = '#EF3A43', lwd = .3) +
    geom_point(size = .3, color = '#EF3A43') +
    # PEAKS A10
    geom_hline(yintercept = -1.6935, color = '#B7CA54', lwd = .2) +
    # instrument error
    geom_hline(yintercept = -1.778, color = '#529986', lwd = .2) +
    xlab('RT [s]') +
    ylim(c(NA, 0)) +
    ylab('Error [ppm]') +
    theme_linedraw()

ggsave(
    'STARD10_A10_pos_786.5_ppms.pdf',
    device = cairo_pdf,
    height = 2,
    width = 5
)

#28.02.2019
PPIterative <- read.table(
    'A10_pos_best_samples.txt',
    col.names = c('RT', 'mz_measured', 'mz_theoretical', 'Lipid name', 'ppm')
)


PPIterative$V3 <- ifelse(mz_measured < 786.7, 0, ifelse(mz_measured < 787.7, 1, ifelse(mz_measured < 788.7, 2, 3)))
delta <- 0.01

PPIterative$V3 = ifelse(
  abs(PPIterative$mz_measured - mz_database$mz_theoretical) < 0.01,
  mz_database$mz_theoretical)

PPIterative$V3 <- mz_database$mz_theoretical


PPIterative2 <- suppressWarnings(suppressMessages(
  read_table2('/home/igor/Documents/Lipyd/openms_dev/pc_best_samples_A10_pos.txt')
)) %>%
  mutate(
    PPIterative$V4 <- ifelse(
      abs(PPIterative$mz_measured - mz_database$mz_theoretical < 0.01),
      PPIterative$V4 <-  mz_database$isotope))
    
 
apply(PPIterative, MARGIN = 2, FUN = ifelse(abs(PPIterative - mz_database$mz_theoretical < 0.01, PPIterative$V3 == mz_database$mz_theoretical )))

  
PPIterative$V3 <- PPIterative %>% mutate(ifelse(mz_measured < 704.7, 0))

for (i in PPIterative$mz_measured)
  for (k in mz_database$mz_theoretical)
    ifelse (abs(i - k) < 0.01, PPIterative$V3 == mz_database$mz_theoretical, no = NA)
  
PPIter3 <- read.table('/home/igor/Documents/Lipyd/openms_dev/result_iterative.txt', sep = " ")
PPIter3$V6 = (PPIter3$V4 - PPIter3$V3) / PPIter3$V3 * 1e06


facet_plot <- ggplot(PPIter3, aes(V5, V6)) +
  geom_line(color = '#EF3A43', lwd = .3) +
  geom_point(size = .3, color = '#EF3A43') +
  # PEAKS A10
  geom_hline(yintercept = -1.6935, color = '#B7CA54', lwd = .1) +
  # instrument error
  geom_hline(yintercept = -1.778, color = '#529986', lwd = .1) +
  facet_grid(rows = vars(V2)) +
  xlab('RT [s]') +
  ylim(c(NA, 0)) +
  ylab('Error (measured vs. theoretical) [ppm]') +
  theme_linedraw()


