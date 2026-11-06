####
data <- read.csv("metadata.csv")

data <- read.csv("PlottingUpdateJun25.csv")


status_colors <- c(
  "_"   = "grey60",
  "_S"  = "steelblue",
  "L_"  = "orange",
  "L_S" = "darkgreen"
)

pdf("figures/statusJUN2025.pdf",width = 8,height = 3)
plot(data$Est.Cal.YrBP,as.factor(data$CoreID),ylim=c(0.5,3.5),xlim=c(12000,0),pch="|",col="grey",cex=1.5,yaxt = "n",ylab="",xlab="Estimated Cal. Years BP")
points(data$Est.Cal.YrBP[grep("L",data$Status)],as.factor(data$CoreID)[grep("L",data$Status)],ylim=c(0.5,3.5),xlim=c(12000,0),pch="-",col="grey25")
axis(2, at = 1:3, labels = levels(as.factor(data$CoreID)), las = 1)
dev.off()

status_colors <- c(
  "_"   = "grey60",
  "_S"  = "steelblue",
  "L_"  = "orange",
  "L_S" = "darkgreen"
)

status_labels <- c(
  "_"   = "Being sequenced",
  "_S"  = "Sequenced",
  "L_"  = "Lipids only",
  "L_S" = "Sequenced + Lipids"
)

# Create a standalone legend plot
pdf("figures/statusJUN2025_legend.pdf", width = 4, height = 3)

# Empty plot
plot.new()
par(mar = c(0, 0, 0, 0))  # remove margins
legend("center",
       legend = status_labels,
       col = status_colors,
       pch = 15,
       pt.cex = 2,
       cex = 1.2,
       bty = "n",
)

dev.off()


# Example usage in a plot:
plot(x, y, col = status_colors[as.factor(data$Status)])



datin <- data[data$CoreID=="PC012",]

pdf("figures/PC012samples.pdf",height = 9,width = 2)
plot(rep(1,length(datin$age_mean_.cal_BP.)),datin$age_mean_.cal_BP.,pch="-",cex=3,ylab="Cal Yr BP",xlab="",xaxt='n',col="purple4")
dev.off()


datin <- data[data$CoreID=="PC014",]

pdf("figures/PC014samples.pdf",height = 14,width = 2)
plot(rep(1,length(datin$age_interpolated_May2024)),datin$age_interpolated_May2024,pch="-",cex=3,ylab="Cal Yr BP",xlab="",xaxt='n')
dev.off()


datin <- data[data$CoreID=="GC017",]

pdf("figures/GC017samples.pdf",height = 14,width = 2)
plot(rep(1,length(datin$age_interpolated_May2024)),datin$age_interpolated_May2024,pch="-",cex=3,ylab="Cal Yr BP",xlab="",xaxt='n')
dev.off()


data2 <- read.csv("LIPIDv2.csv")

datin <- data2[data2$Core=="GC017",]
palette(brewer.pal(6, "Set2"))


plot(data2$Age,as.numeric(as.factor(data2$Core)),pch="|",col=as.numeric(as.factor(data2$eDNAorLipid))+2,yaxt="n",cex=2,ylim=c(0,4))
axis(2,1:3,labels = levels(as.factor(data2$Core)))


pdf("figures/GC017lipid.pdf",height =11,width = 2)
plot(rep(1,length(datin$Age)),1950-datin$Age,pch="-",cex=3,ylab="Year",xlab="",xaxt='n',yaxt='n',ylim=c(-9000,2000),bty="n",col=factor(datin$eDNAorLipid,levels=c('B','E','L')))
axis(2,at=seq(-8000,2000,2000),labels=paste0(sqrt(seq(-8000,2000,2000)^2),c("BCE","BCE","BCE","BCE","","CE")))
dev.off()

datin <- data2[data2$Core=="PC014",]

pdf("figures/PC014lipid.pdf",height = 11,width = 2)
plot(rep(2,length(datin$Age)),1950-datin$Age,pch="-",cex=3,ylab="year",xlab="",xaxt='n',yaxt='n',ylim=c(-9000,2000),bty="n",col=factor(datin$eDNAorLipid,levels=c('B','E','L')))
axis(2,at=seq(-8000,2000,2000),labels=paste0(sqrt(seq(-8000,2000,2000)^2),c("BCE","BCE","BCE","BCE","","CE")))
dev.off()


datin <- data2[data2$Core=="PC012",]

pdf("figures/PC012lipid.pdf",height =11,width = 2)
plot(rep(2,length(datin$Age)),1950-datin$Age,pch="-",cex=3,ylab="year",xlab="",xaxt='n',yaxt='n',ylim=c(-9000,2000),bty="n",col=factor(datin$eDNAorLipid,levels=c('B','E','L')))
axis(2,at=seq(-8000,2000,2000),labels=paste0(sqrt(seq(-8000,2000,2000)^2),c("BCE","BCE","BCE","BCE","","CE")))
dev.off()


datin <- data[data$CoreID=="GC017",]

pdf("figures/GC017allsamples.pdf",height = 14,width = 3)
plot(rep(1,length(datin$age_interpolated_May2024)),datin$age_interpolated_May2024,pch="-",cex=3,ylab="Cal Yr BP",xlab="",xaxt='n',xlim=c(0,3))
datin <- data2[data2$Core=="GC017",]
points(rep(2,length(datin$Age)),datin$Age,pch="-",col="darkred",cex=3)
dev.off()


datin <- data[data$CoreID=="PC014",]

pdf("figures/PC014allsamples.pdf",height = 14,width = 3)
plot(rep(1,length(datin$age_interpolated_May2024)),datin$age_interpolated_May2024,pch="-",cex=3,ylab="Cal Yr BP",xlab="",xaxt='n',xlim=c(0,3))
datin <- data2[data2$Core=="PC014",]
points(rep(2,length(datin$Age)),datin$Age,pch="-",col="darkred",cex=3)
dev.off()





