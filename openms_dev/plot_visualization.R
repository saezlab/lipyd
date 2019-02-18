require(ggplot2)
require(viridis)
require(dplyr)

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

require(readr)
require(dplyr)

help(read_table)

teo <- suppressWarnings(suppressMessages(read_table2(
  "/home/igor/Documents/Black_scripts/Comparison/teo.txt",
  col_names = c('rt', 'mz_measured', 'nothing')
))) %>%
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
help(facet_grid)
facet_plot <- ggplot(teo, aes(rt, ppm)) +
  geom_line(color = 'black') +
  geom_point(size = .3) +
  facet_grid(rows = vars(isotope)) +
  theme_minimal()

facet_plot + facet_grid(rows = vars('RT'))
facet_plot + facet_grid(cols = vars('ppm'))
facet_plot + facet_grid(vars('RT'), vars('ppm'))

head(teo)
ggsave('Facet_plot_STARD10_A10_pos', device = "jpeg", height = 4, width = 5)
