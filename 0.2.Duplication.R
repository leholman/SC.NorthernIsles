####Northern isles duplication

metadata <- read.csv("metadata.csv")

plot(metadata$Dct,metadata$duplication.rate,pch=16,col=as.factor(metadata$CoreID),xlab="Ct difference",ylab="Duplication Rate (%)")
legend("topleft",legend = c("Controls","GC017","PC012","PC014"),pch=16,col=c(1,2,3,4))
