# compare statistical distributions of cotton strip decay rates between those observed in the YRB and those from a global dataset

global.cotton.Kcd = read.csv(file = "~/GitHub/rcsfa-ST-2B-SSS-cotton-strip/data/str_k.csv" )
str(global.cotton.Kcd)
range(global.cotton.Kcd$k)

global.cotton.Kdd = read.csv(file = "~/GitHub/rcsfa-ST-2B-SSS-cotton-strip/data/str_k_dd.csv" )
str(global.cotton.Kdd)
range(global.cotton.Kdd$k)

#
yrb.cotton = read.csv(file = "~/GitHub/rcsfa-ST-2B-SSS-cotton-strip/Outputs/Decay_Data.csv")
yrb.trim = yrb.cotton[,c("degree_decay_rate","Decay_Rate_per_day")]
range(yrb.trim$degree_decay_rate)
range(yrb.trim$Decay_Rate_per_day)

# generate kernels
global.Kcd.kern = density(global.cotton.Kcd$k,from = 0)
global.Kdd.kern = density(global.cotton.Kdd$k,from = 0)
yrb.Kdd.kern = density(yrb.trim$degree_decay_rate,from = 0)
yrb.Kcd.kern = density(yrb.trim$Decay_Rate_per_day,from = 0)

# make plots

pdf("~/GitHub/rcsfa-ST-2B-SSS-cotton-strip/Outputs/Global_YRB_Kern.pdf",width = 10,height = 5)

par(pty="s",mfrow=c(1,2))

plot(global.Kcd.kern,typ="l",main="",xlab = "Kcd",lwd=2,xlim=c(0,0.25),cex.lab=2,cex.axis=1.5)
points(yrb.Kcd.kern,typ="l",col=2,lwd=2)
legend(x = 0.1,y=20,col = c(1,2),lwd=2,legend = c('Global','YRB'),cex = 1.25)

plot(global.Kdd.kern,typ="l",main="",xlab = "Kdd",lwd=2,xlim=c(0,0.025),ylim=c(0,250),cex.lab=2,cex.axis=1.5)
points(yrb.Kdd.kern,typ="l",col=2,lwd=2)

dev.off()
