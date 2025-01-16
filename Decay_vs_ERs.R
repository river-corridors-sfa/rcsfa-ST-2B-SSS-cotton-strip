require(readr)

# read in ER data
ER.data = read_csv("v2_SSS_Water_Sediment_Total_Respiration_GPP.csv",comment = '#')
head(ER.data)

# read in decay rate data
Decay.data = read_csv("Outputs/Decay_Data.csv")
head(Decay.data)
head(Decay.data[,c('Parent_ID','mean_degree_decay_rate')])
length((unique(Decay.data$Parent_ID))) # 46

# get mean per day decay rate
mean.decay.per.day = as.data.frame(tapply(X = Decay.data$Decay_Rate_per_day,INDEX = Decay.data$Parent_ID,FUN = 'mean'))
colnames(mean.decay.per.day) = c('mean_day_decay_rate')
mean.decay.per.day$Parent_ID = rownames(mean.decay.per.day)
head(mean.decay.per.day)

# merge mean decay with decay data
Decay.data = merge(x = Decay.data,y = mean.decay.per.day,by = 'Parent_ID')
head(Decay.data)
dim(Decay.data)

# drop unneeded columns and reduce to a single row per site (parent id)
Decay.data = unique(Decay.data[,c('Parent_ID','mean_day_decay_rate','mean_degree_decay_rate')])
head(Decay.data)
length((unique(Decay.data$Parent_ID))) # should be 46 (same as length above on line 11)

# merge decay rates with ERs
merged.data = merge(ER.data,Decay.data,by = 'Parent_ID',all = T)
head(merged.data)

# remove sites with missing data
merged.data = merged.data[-grep(pattern = -9999,x = merged.data$Sediment_Respiration),]
merged.data = merged.data[-which(is.na(merged.data$mean_day_decay_rate)==T),]
head(merged.data)
dim(merged.data)

# cube root transform
cube_root <- function(x) sign(x) * (abs(x))^(1/3)
merged.data$cube.root_mean_day_decay_rate = cube_root(merged.data$mean_day_decay_rate)
merged.data$cube.root_mean_degree_decay_rate = cube_root(merged.data$mean_degree_decay_rate)
merged.data$cube.root_Total_Ecosystem_Respiration = cube_root(merged.data$Total_Ecosystem_Respiration)
merged.data$cube.root_Sediment_Respiration = cube_root(merged.data$Sediment_Respiration)
merged.data$cube.root_Water_Column_Respiration = cube_root(merged.data$Water_Column_Respiration)

# generate the 6 panel plot

pdf("Outputs/Decay_vs_ERs.pdf", width = 14, height = 9)

# Set plotting parameters
par(pty = "s", mfrow = c(2, 3), oma = c(4, 4, 4, 4), mar = c(5, 5, 2, 1), mgp = c(3.5, 1, 0))

# per day decay plots
mod.to.plot = merged.data$cube.root_mean_day_decay_rate ~ merged.data$cube.root_Total_Ecosystem_Respiration
mod.sum = summary(lm(mod.to.plot))
plot(mod.to.plot, 
     xlab = "", # expression(ER[tot]~(g~O[2]~m^-2~day^-1)^{1/3}), 
     ylab = expression(Decay~Rate~(day^-1)^{1/3}), 
     cex.lab = 2, cex.axis = 1.5)
abline(mod.sum,lwd=3)
mtext(bquote(R^2 == .(round(mod.sum$r.squared, digits = 2))), side = 3, adj = 0.95, line = -2)
mtext(paste0("p = ",round(mod.sum$coefficients[2,4],digits = 3)),side = 3,adj = 0.95,line = -3.5)
mtext("A",cex=2,side = 1,adj = 0.05,line=-1.5)

mod.to.plot = merged.data$cube.root_mean_day_decay_rate ~ merged.data$cube.root_Sediment_Respiration
mod.sum = summary(lm(mod.to.plot))
plot(mod.to.plot, 
     xlab = "", # expression(ER[sed]~(g~O[2]~m^-2~day^-1)^{1/3}), 
     ylab = "", # expression(Decay~Rate~(day^-1)^{1/3}), 
     cex.lab = 2, cex.axis = 1.5)
abline(mod.sum,lwd=3)
mtext(bquote(R^2 == .(round(mod.sum$r.squared, digits = 2))), side = 3, adj = 0.95, line = -2)
mtext(paste0("p = ",round(mod.sum$coefficients[2,4],digits = 3)),side = 3,adj = 0.95,line = -3.5)
mtext("B",cex=2,side = 1,adj = 0.05,line=-1.5)

mod.to.plot = merged.data$cube.root_mean_day_decay_rate ~ merged.data$cube.root_Water_Column_Respiration
mod.sum = summary(lm(mod.to.plot))
plot(mod.to.plot, 
     xlab = "", # expression(ER[wc]~(g~O[2]~m^-2~day^-1)^{1/3}), 
     ylab = "", #expression(Decay~Rate~(day^-1)^{1/3}), 
     cex.lab = 2, cex.axis = 1.5)
abline(mod.sum,lwd=3)
mtext(bquote(R^2 == .(round(mod.sum$r.squared, digits = 2))), side = 3, adj = 0.90, line = -2)
mtext(paste0("p = ",round(mod.sum$coefficients[2,4],digits = 3)),side = 3,adj = 0.90,line = -3.5)
mtext("C",cex=2,side = 1,adj = 0.05,line=-1.5)

# per degree day decay plots
mod.to.plot = merged.data$cube.root_mean_degree_decay_rate ~ merged.data$cube.root_Total_Ecosystem_Respiration
mod.sum = summary(lm(mod.to.plot))
plot(mod.to.plot, 
     xlab = expression(ER[tot]~(g~O[2]~m^-2~day^-1)^{1/3}), 
     ylab = bquote(Decay ~ Rate ~ ('degree ' ~ day^-1)^{1/3}), 
     cex.lab = 2, cex.axis = 1.5)
abline(mod.sum,lwd=3)
mtext(bquote(R^2 == .(round(mod.sum$r.squared, digits = 2))), side = 3, adj = 0.95, line = -2)
mtext(paste0("p = ",round(mod.sum$coefficients[2,4],digits = 3)),side = 3,adj = 0.95,line = -3.5)
mtext("D",cex=2,side = 1,adj = 0.05,line=-1.5)

mod.to.plot = merged.data$cube.root_mean_degree_decay_rate ~ merged.data$cube.root_Sediment_Respiration
mod.sum = summary(lm(mod.to.plot))
plot(mod.to.plot, 
     xlab = expression(ER[sed]~(g~O[2]~m^-2~day^-1)^{1/3}), 
     ylab = "", # bquote(Decay ~ Rate ~ ('degree ' ~ day^-1)^{1/3}), 
     cex.lab = 2, cex.axis = 1.5)
abline(mod.sum,lwd=3)
mtext(bquote(R^2 == .(round(mod.sum$r.squared, digits = 2))), side = 3, adj = 0.95, line = -2)
mtext(paste0("p = ",round(mod.sum$coefficients[2,4],digits = 3)),side = 3,adj = 0.95,line = -3.5)
mtext("E",cex=2,side = 1,adj = 0.05,line=-1.5)

mod.to.plot = merged.data$cube.root_mean_degree_decay_rate ~ merged.data$cube.root_Water_Column_Respiration
mod.sum = summary(lm(mod.to.plot))
plot(mod.to.plot, 
     xlab = expression(ER[wc]~(g~O[2]~m^-2~day^-1)^{1/3}), 
     ylab = "", #bquote(Decay ~ Rate ~ ('degree ' ~ day^-1)^{1/3}), 
     cex.lab = 2, cex.axis = 1.5)
abline(mod.sum,lwd=3)
mtext(bquote(R^2 == .(round(mod.sum$r.squared, digits = 2))), side = 3, adj = 0.9, line = -2)
mtext(paste0("p = ",round(mod.sum$coefficients[2,4],digits = 3)),side = 3,adj = 0.9,line = -3.5)
mtext("F",cex=2,side = 1,adj = 0.05,line=-1.5)

dev.off()


# do multiple regression compared to univariate

day_mult_mod = (lm(merged.data$cube.root_mean_day_decay_rate ~ merged.data$cube.root_Water_Column_Respiration + merged.data$cube.root_Sediment_Respiration))
day_mult_inter_mod = (lm(merged.data$cube.root_mean_day_decay_rate ~ merged.data$cube.root_Water_Column_Respiration + merged.data$cube.root_Sediment_Respiration + merged.data$cube.root_Water_Column_Respiration*merged.data$cube.root_Sediment_Respiration))
degree_mult_mod = (lm(merged.data$cube.root_mean_degree_decay_rate ~ merged.data$cube.root_Water_Column_Respiration + merged.data$cube.root_Sediment_Respiration))
degree_mult_inter_mod = (lm(merged.data$cube.root_mean_degree_decay_rate ~ merged.data$cube.root_Water_Column_Respiration + merged.data$cube.root_Sediment_Respiration + merged.data$cube.root_Water_Column_Respiration*merged.data$cube.root_Sediment_Respiration))

day_ERtot_mod = (lm(merged.data$cube.root_mean_day_decay_rate ~ merged.data$cube.root_Total_Ecosystem_Respiration))
degree_ERtot_mod = (lm(merged.data$cube.root_mean_degree_decay_rate ~ merged.data$cube.root_Total_Ecosystem_Respiration))

# outcome is that univariate is a better model than multivariate for both rates
AIC(day_ERtot_mod) # [1] -59.01036
AIC(day_mult_mod) # [1] -55.24735
AIC(day_mult_inter_mod) # -53.27437

AIC(degree_ERtot_mod) # [1] -113.1095
AIC(degree_mult_mod) # [1] -110.0843
AIC(degree_mult_inter_mod) # -108.3438


