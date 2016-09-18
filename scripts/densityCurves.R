args=(commandArgs(TRUE))
label=args[1]
table1=paste("known",label,sep=".")
table2=paste("novel",label,sep=".")


library(ggplot2)

dataA=read.table(table1)
dataB=read.table(table2)

dataA$V1="known"
dataB$V1="novel"

allData <- rbind(dataA, dataB)

outputPDF <- paste(label,"pdf",sep=".")
pdf(outputPDF)
ggplot(allData, aes(allData$V2, fill = V1)) + geom_density(alpha = 0.2)
dev.off()


## SNP.FS
#outputPDF <- paste(label,"2","pdf",sep=".")
#pdf(outputPDF)
#ggplot(allData, aes(allData$V2, fill = V1)) + geom_density(alpha = 0.2) + coord_cartesian(xlim = c(-20, 120))
#dev.off()

## SNP.MG
#outputPDF <- paste(label,"2","pdf",sep=".")
#pdf(outputPDF)
#ggplot(allData, aes(allData$V2, fill = V1)) + geom_density(alpha = 0.2) + coord_cartesian(xlim = c(50, 63))
#dev.off()

## SNP.MQRankSum
#outputPDF <- paste(label,"2","pdf",sep=".")
#pdf(outputPDF)
#ggplot(allData, aes(allData$V2, fill = V1)) + geom_density(alpha = 0.2) + coord_cartesian(xlim = c(-4, 4))
#dev.off()

## SNP.MQRankSum
#outputPDF <- paste(label,"2","pdf",sep=".")
#pdf(outputPDF)
#ggplot(allData, aes(allData$V2, fill = V1)) + geom_density(alpha = 0.2) + coord_cartesian(xlim = c(-4, 4))
#dev.off()

## SNP.DP
#outputPDF <- paste(label,"2","pdf",sep=".")
#pdf(outputPDF)
#ggplot(allData, aes(allData$V2, fill = V1)) + geom_density(alpha = 0.2) + coord_cartesian(xlim = c(-100, 5000))
#dev.off()


