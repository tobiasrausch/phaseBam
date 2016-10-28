library(ggplot2)
library(scales)
library(RColorBrewer)

args=commandArgs(trailingOnly=TRUE)
swr =read.table(args[1], header=T)
swr = swr[!apply(is.na(swr), 1, any),]
swr$distance = swr$distance + 1


png("swr.png", width=800, height=600)
p1=ggplot(data=swr, aes(x=acPre, y=acSuc)) + geom_point(aes(color=distance))
p1=p1 + scale_colour_gradientn(colours=brewer.pal(3, "Set1"), trans="log", breaks=c(1,10,100,1000,10000,100000,1000000,10000000), labels=c("1","10","100","1Kbp","10Kbp","100Kbp","1Mbp","10Mbp"))
p1=p1 + labs(color="Distance to previous\nswitch error")
p1=p1 + xlab("Allele Count Preceeding SNP")
p1=p1 + ylab("Allele Count Succeeding SNP")
#p1=p1 + scale_x_log10() + scale_y_log10()
p1
dev.off()

ac = c(swr$acPre, swr$acSuc)
png("histAC.png")
hist(log(ac), main="Log Allele Count Distribution of\nPreceeding and Succeeding SNP\nadjacent to Switch Error", xlab="log(AC)")
dev.off()

png("histDistance.png")
hist(log(swr$distance), main="Log Distance to previous Switch Error (in bp)", xlab="log(distance)")
dev.off()

print(warnings())

