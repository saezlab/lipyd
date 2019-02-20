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