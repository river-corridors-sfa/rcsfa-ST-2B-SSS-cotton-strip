# compare statistical distributions of cotton strip decay rates between those observed in the YRB and those from a global dataset
# the global dataset is from: https://zenodo.org/records/14166081

global.cotton = read.csv(file = "~/GitHub/rcsfa-ST-2B-SSS-cotton-strip/data/LeRoy.ExpandedDataset.Kvalues.csv" )
str(global.cotton)

# calculate decay per degree day
global.cotton$Kdd = global.cotton$k..d./as.numeric(global.cotton$mean..ºC.)

# cleaning data due to missing data
global.trim = global.cotton[,c('Kdd','k..d.')]
global.trim = global.trim[-which(is.na(global.trim$Kdd)==T |  is.na(global.trim$k..d.)==T | global.trim$Kdd==Inf),]
head(global.trim)
range(global.trim$Kdd)
range(global.trim$k..d.)

# just a check
plot(global.trim$Kdd ~ global.trim$k..d.)

#
yrb.cotton = read.csv(file = "~/GitHub/rcsfa-ST-2B-SSS-cotton-strip/Outputs/Decay_Data.csv")
yrb.trim = yrb.cotton[,c("degree_decay_rate","Decay_Rate_per_day")]
range(yrb.trim$degree_decay_rate)
range(yrb.trim$Decay_Rate_per_day)

# generate kernels
global.Kdd.kern = density(global.trim$Kdd,from = 0)
global.Kcd.kern = density(global.trim$k..d.,from = 0)
yrb.Kdd.kern = density(yrb.trim$degree_decay_rate,from = 0)
yrb.Kcd.kern = density(yrb.trim$Decay_Rate_per_day,from = 0)

# make plots

pdf("~/GitHub/rcsfa-ST-2B-SSS-cotton-strip/Outputs/Global_YRB_Kern.pdf",width = 10,height = 5)

par(pty="s",mfrow=c(1,2))

quantile(x = global.trim$k..d.,probs = 0.99)
plot(global.Kcd.kern,typ="l",xlim=c(0,0.25),main="",xlab = "Kcd",lwd=2) # setting x range to be over 99% of global data, based on the previous line
points(yrb.Kcd.kern,typ="l",col=2,lwd=2)
legend(x = 0.1,y=40,col = c(1,2),lwd=2,legend = c('Global','YRB'),cex = 1.25)

quantile(x = global.trim$Kdd,probs = 0.95)
plot(global.Kdd.kern,typ="l",xlim=c(0,0.02),main="",xlab = "Kdd",lwd=2) # setting x range to be over 95% of global data, based on the previous line
points(yrb.Kdd.kern,typ="l",col=2,lwd=2)

dev.off()
