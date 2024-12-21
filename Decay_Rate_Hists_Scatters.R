data = read.csv("Outputs/Decay_Data.csv")
head(data)
str(data)

# function to generate scatter plot with adjacent histograms
# original version modified and pulled from https://www.r-bloggers.com/2011/06/example-8-41-scatterplot-with-marginal-histograms/
# some hard coding in this function to tweak the plot aesthetics

scatterhist = function(x, y, xlab="", ylab=""){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE,breaks=20)
  yhist = hist(y, plot=FALSE,breaks=20)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(5,5,1,1))
  plot(x,y,cex.axis=1.5,cex=1.5,xlab="",ylab="")
  points(lowess(y~x),typ="l",col=4,lwd=2)
  par(mar=c(0,5,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(5,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, cex = 2,
        at=.35 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0,cex = 2, 
        at=(.6 * (mean(y) - min(y))/(max(y) - min(y))))
}

# generate the plot

pdf("Outputs/Decay_Scatter_Hists.pdf")
  scatterhist(x = data$degree_decay_rate,y = data$Decay_Rate_per_day,xlab = "Decay Rate (per degree day)",ylab = "Decay Rate (per day)")
dev.off()
