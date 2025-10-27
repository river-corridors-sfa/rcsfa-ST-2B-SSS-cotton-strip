# ==============================================================================
#
# Make figures for ERsed manuscript
#
# Status: Complete
#
# ==============================================================================
#
# Author: Brieanne Forbes 
# 23 January 2025
#
# ==============================================================================
library(tidyverse) 
# library(segmented) #only needed if using segmented regression which we are not anymore
library(ggpubr)
library(ggpmisc)

rm(list=ls(all=T))

# Setting wd to parent folder
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
setwd("./..")

# =================================== find files ===============================

ER <- './v2_SSS_Water_Sediment_Total_Respiration_GPP.csv' %>%
  read_csv(comment = '#', na = '-9999') %>%
  mutate(Total_Ecosystem_Respiration= case_when(Total_Ecosystem_Respiration > 0 ~ NA, # set positive ERtot values to NA
                                                TRUE ~ Total_Ecosystem_Respiration))


# pulls gap-filled (QAQCd) data from Bernhardt et al. (2022). Repo: https://github.com/streampulse/metabolism_synthesis/tree/master
ER_lit <- readRDS(url("https://raw.githubusercontent.com/streampulse/metabolism_synthesis/master/output_data/lotic_gap_filled.rds")) %>%
  bind_rows() %>%
  group_by(Site_ID) %>%
  summarise(mean_ER = mean(ER_filled, na.rm = T))

ER_hz_model <- './v2_SSS_ER_d50_TotalOxygenConsumed.csv' %>%
  read_csv(comment = '#', na = '-9999')

# =============================== set theme ====================================
theme_set(
  theme(
    text = element_text(family = 'serif', face = 'plain'),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    line = element_line(size = 0.05),
    axis.line = element_line(size = 0.5),
    panel.background = element_rect(color = 'white'),
    panel.border = element_rect(
      colour = 'black',
      fill = NA,
      size = 0.5
    ),
    plot.title = element_text(size = 16, face = 'bold', hjust = 0, vjust = 1.5),
    axis.ticks.length = unit(.25, 'cm'),
    plot.subtitle = element_text(size = 14),
    legend.title = element_blank()
  )
)

# ====================== Density:ERtot, ERsed, ERwc =============================

p2 <- ggplot(data=ER) + 
  geom_segment(aes( x = min(Water_Column_Respiration), xend = max(Water_Column_Respiration), 
                    y = 0.08, yend = 0.08, color = 'wc'), 
               size = 1, linetype = 'dashed') +
  geom_segment(aes(x = min(Water_Column_Respiration), xend = min(Water_Column_Respiration), 
                   y = 0.075, yend = 0.085, color = 'wc'), 
               size = 1) +
  geom_segment(aes(x = max(Water_Column_Respiration), xend = max(Water_Column_Respiration), 
                   y = 0.075, yend = 0.085, color = 'wc'),
               size = 1) + #ER WC
  geom_density(data = ER_lit, aes(x=mean_ER,colour="lit",fill='lit'),adjust = 6)+
  geom_density(aes(x=Sediment_Respiration,colour="sed",fill='sed'),adjust = 6,alpha=0.8)+
  geom_density(aes(x=Total_Ecosystem_Respiration,colour="tot",fill='tot'),adjust = 6,alpha=0.6)+
  geom_vline(data = ER_lit, aes(xintercept=median(mean_ER, na.rm = T)),color="black",  size=1)+
  geom_vline(aes(xintercept=median(Total_Ecosystem_Respiration, na.rm = T)),color="grey",  size=1)+
  geom_vline(aes(xintercept=median(Sediment_Respiration, na.rm = T)),color="coral4",  size=1, linetype = "dashed")+
  labs(x = expression(paste("Ecosystem Respiration"*" (g O"[2]*" m"^-2*" day"^-1*")")), y = 'Density')+
  scale_fill_manual("",breaks = c("tot",'sed', 'lit'),labels = c(expression("ER"[tot]*" (YRB)"),expression("ER"[sed]*" (YRB)"),expression("ER"[tot]*" (Lit)")),
                    values = c('grey','coral4', "black"))+
  scale_colour_manual("",breaks = c("tot","sed",'wc', 'lit'),labels = c(expression("ER"[tot]*" (YRB)"),expression("ER"[sed]*" (YRB)"),expression("ER"[wc]*" (YRB) range"),expression("ER"[tot]*" (Lit)")),
                      values = c('darkgrey','coral4','blue', "black")    
  )+
  theme(
    legend.position = c(.2, .95),
    legend.justification = c( "top"),
    legend.text = element_text(size=12,hjust = 0), #, margin = margin(l = 0, r = 5, unit = "pt")
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
    legend.key = element_rect(fill = "white", color = "black", linewidth = 0.2),
    legend.box.just = "right",
    legend.title = element_blank()
  )


# ====================== Density: ERwc =============================

p3 <- ggplot(data=ER) +
  geom_density(aes(x=Water_Column_Respiration,colour="wc",fill='wc'),adjust = 6)+
  geom_vline(aes(xintercept=median(Water_Column_Respiration, na.rm = T)),color="blue",  size=1)+
  labs(x = expression(paste("Ecosystem Respiration"*" (g O"[2]*" m"^-2*" day"^-1*")")), y = 'Density')+
  scale_fill_manual("",breaks = c('wc'),labels = c(expression("ER"[wc]*" (YRB)")),
                    values = c('skyblue'))+
  scale_colour_manual("",breaks = c('wc'),labels = c(expression("ER"[wc]*" (YRB)")),
                      values = c('blue')
  )+ theme(
    legend.position = c(.2, .95),
    legend.justification = c( "top"),
    legend.text = element_text(size=12,hjust = 0), #, margin = margin(l = 0, r = 5, unit = "pt")
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
    legend.key = element_rect(fill = "white", color = "black", linewidth = 0.2),
    legend.box.just = "right"
  )

density <- ggarrange(
  p2,
  p3,
  ncol = 2,
  nrow = 1,
  widths = c(5, 5),
  heights = 3
)

ggsave('./Figures/Figure3_ERtot_ERsed_ERlit_ERwc_Density.pdf',
       density,
       device = 'pdf',
       width = 10.5,
       height = 3,
       units = 'in',
       dpi = 300
)

# ====================== Scatter: GPP vs ERtot+ERsed =============================
# compare linear and segmented regression
# there is barely a difference between the AIC and BIC between the two regressions. 
# Both have a (very) slightly lower value for linear so going with this for simplicity

# fit <- lm(Sediment_Respiration~Gross_Primary_Production, data=ER)
# summary(fit)
# AIC(fit)
# BIC(fit)
# 
# segmented.fit <- segmented(fit, seg.Z = ~Gross_Primary_Production, psi=9,
#                            control =seg.control(display = TRUE, maxit.glm=3))
# summary(segmented.fit)
# AIC(segmented.fit)
# BIC(segmented.fit)


p4 <- ggplot(ER, aes(x = Gross_Primary_Production)) +
  geom_abline(slope = -1, intercept = 0, color = 'darkgrey', linetype = 'dashed') +
  geom_point(aes(y = Total_Ecosystem_Respiration), color = "grey32", size = 2.5) +
  geom_smooth(aes(y = Total_Ecosystem_Respiration), method = "lm", 
              color = "grey32", fill = "grey32", alpha = 0.2) +
  geom_point(aes(y = Sediment_Respiration), color = "coral4", size = 2.5) +
  geom_smooth(aes(y = Sediment_Respiration), method = "lm", 
              color = "coral4", fill = "coral4", alpha = 0.2) +
  annotate("text", x = 0, y = -16, 
           label = paste(
             "ERtot: y =",
             sprintf("%.2f", coef(lm(Total_Ecosystem_Respiration ~ Gross_Primary_Production, data = ER))[2]),
             "*x +",
             sprintf("%.2f", coef(lm(Total_Ecosystem_Respiration ~ Gross_Primary_Production, data = ER))[1]),
             "\nR² =",
             sprintf("%.2f", summary(lm(Total_Ecosystem_Respiration ~ Gross_Primary_Production, data = ER))$r.squared),
             "\n",
             ifelse(summary(lm(Total_Ecosystem_Respiration ~ Gross_Primary_Production, data = ER))$coefficients[8] < 0.001,
                    "p < 0.001",
                    paste("p =", sprintf("%.3f", summary(lm(Total_Ecosystem_Respiration ~ Gross_Primary_Production, data = ER))$coefficients[8]))
             )),
           color = "grey32", vjust = -0.5, hjust = 0, family = 'serif') +
  annotate("text", x = 0, y = -19, 
           label = paste(
             "ERsed: y =",
             sprintf("%.2f", coef(lm(Sediment_Respiration ~ Gross_Primary_Production, data = ER))[2]),
             "*x +",
             sprintf("%.2f", coef(lm(Sediment_Respiration ~ Gross_Primary_Production, data = ER))[1]),
             "\nR² =",
             sprintf("%.2f", summary(lm(Sediment_Respiration ~ Gross_Primary_Production, data = ER))$r.squared),
             "\n",
             ifelse(summary(lm(Sediment_Respiration ~ Gross_Primary_Production, data = ER))$coefficients[8] < 0.001,
                    "p < 0.001",
                    paste("p =", sprintf("%.3f", summary(lm(Sediment_Respiration ~ Gross_Primary_Production, data = ER))$coefficients[8]))
             )),
           color = "coral4", vjust = 0.5, hjust = 0, family = 'serif') +
  labs(x = expression(paste("Gross Primary Productivity"*" (g O"[2]*" m"^-2*" day"^-1*")")), 
       y = expression(paste("Ecosystem Respiration"*" (g O"[2]*" m"^-2*" day"^-1*")")))



ggsave('./Figures/Figure5_ERtot_vs_GPP_Regression.pdf',
       p4,
       device = 'pdf',
       width = 5,
       height = 4.5,
       units = 'in',
       dpi = 300
)

# =========== Scatter: Normalized + Rank order: ERtot vs ERhz =================

combine_ER <- ER %>%
  full_join(ER_hz_model)%>%
  mutate(Total_Ecosystem_Respiration_Z = c(scale(Total_Ecosystem_Respiration, center = TRUE, scale = TRUE)),
         HZ_Respiration_Z =  c(scale(Total_Oxygen_Consumed_g_per_m2_per_day, center = TRUE, scale = TRUE)),
         Total_Ecosystem_Respiration_rank = rank(Total_Ecosystem_Respiration),
         HZ_Respiration_rank = rank(Total_Oxygen_Consumed_g_per_m2_per_day))

p5 <- ggplot(data = combine_ER, aes(x = Total_Ecosystem_Respiration_Z, y = HZ_Respiration_Z)) +
  geom_point(alpha = 0.5, size = 3)+
  xlab(expression(paste("Normalized Field-Estimated Total Ecosystem Respiration")))+
  ylab(expression(paste("Normalized Predicted Hyporheic Zone Respiration")))+
  ggtitle('(a)')+
  xlim(-3.5, 1)+
  ylim(-3.5, 1)

p6 <- ggplot(data = combine_ER, aes(x = Total_Ecosystem_Respiration_rank, y = HZ_Respiration_rank)) +
  geom_point(alpha = 0.5, size = 3)+
  xlab(expression(paste("Rank Order - Field-Estimated Total Ecosystem Respiration")))+
  ylab(expression(paste("Rank Order - Predicted Hyporheic Zone Respiration")))+
  stat_cor(method = "spearman",cor.coef.name = c( "rho"),
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 30,label.y = 47,color='black',size= 5, family = 'serif')+
  ggtitle('(b)')

norm_rank <- ggarrange(
  p5,
  p6,
  ncol = 2,
  nrow = 1,
  widths = c(5, 5),
  heights = 5
)

ggsave('./Figures/Figure2_ERtot_ERhz_ZScore_RankOrder.pdf',
       norm_rank,
       device = 'pdf',
       width = 10.5,
       height = 5,
       units = 'in',
       dpi = 300
)

# =========== Scatter: Normalized + Rank order: ERsed vs ERhz =================


combine_ER <- combine_ER %>%
  mutate(Sediment_Respiration_Z = c(scale(Sediment_Respiration, center = TRUE, scale = TRUE)),
         Sediment_Respiration_rank = rank(Sediment_Respiration))

p7 <- ggplot(data = combine_ER, aes(x = Sediment_Respiration_Z, y = HZ_Respiration_Z)) +
  geom_point(alpha = 0.5, size = 3)+
  xlab(expression(paste("Normalized Field-Estimated Sediment Respiration")))+
  ylab(expression(paste("Normalized Predicted Hyporheic Zone Respiration")))+
  ggtitle('(a)')+
  xlim(-4.5, 1)+
  ylim(-4.5, 1)

p8 <- ggplot(data = combine_ER, aes(x = Sediment_Respiration_rank, y = HZ_Respiration_rank)) +
  geom_point(alpha = 0.5, size = 3)+
  xlab(expression(paste("Rank Order - Field-Estimated Sediment Respiration")))+
  ylab(expression(paste("Rank Order - Predicted Hyporheic Zone Respiration")))+
  stat_cor(method = "spearman",cor.coef.name = c( "rho"),
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 30,label.y = 47,color='black',size= 5, family = 'serif')+
  ggtitle('(b)')

ersed_norm_rank <- ggarrange(
  p7,
  p8,
  ncol = 2,
  nrow = 1,
  widths = c(5, 5),
  heights = 5
)

ggsave('./Figures/FigureS2_ERsed_ERhz_ZScore_RankOrder.pdf',
       ersed_norm_rank,
       device = 'pdf',
       width = 10.5,
       height = 5,
       units = 'in',
       dpi = 300
)

# =========== Scatter:  ERsed vs ERtot =================

p9 <- ggplot(data = ER, aes(x = Sediment_Respiration, y = Total_Ecosystem_Respiration)) +
  geom_abline(slope = 1, intercept = 0, color = 'grey32', linetype = 'dashed') +
  geom_point(size = 2.5, alpha = 0.6)+
  labs(x = expression(paste("Sediment Respiration"*" (g O"[2]*" m"^-2*" day"^-1*")")), 
       y = expression(paste("Total Ecosystem Respiration"*" (g O"[2]*" m"^-2*" day"^-1*")")))

summary(lm(Total_Ecosystem_Respiration ~ Sediment_Respiration, data = ER))

ggsave('./Figures/Intermediate_Files/ERsed_ERtot_Scatter.pdf',
       p9,
       device = 'pdf',
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300
)

# =========== Scatter:  ERwc vs ERtot =================

p10 <- ggplot(data = ER, aes(x = Water_Column_Respiration, y = Total_Ecosystem_Respiration)) +
  geom_abline(slope = 1, intercept = 0, color = 'grey32', linetype = 'dashed') +
  geom_point(size = 2.5, alpha = 0.6)+
  labs(x = expression(paste("Water Column Respiration"*" (g O"[2]*" m"^-2*" day"^-1*")")), 
       y = expression(paste("Total Ecosystem Respiration"*" (g O"[2]*" m"^-2*" day"^-1*")")))

ggsave('./Figures/Intermediate_Files/ERwc_ERtot_Scatter.pdf',
       p10,
       device = 'pdf',
       width = 5,
       height = 5,
       units = 'in',
       dpi = 300
)

