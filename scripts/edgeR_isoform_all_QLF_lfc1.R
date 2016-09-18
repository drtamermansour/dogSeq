library(edgeR, lib.loc = "~/R/v3.0.1/library")
library(Biobase)
source("/mnt/ls15/scratch/users/mansourt/Tamer/Hus_plant/scripts/trinity_util/R/rnaseq_plot_funcs.R")

data_count = read.table("eXpress_isoform_matrix.counts.matrix", header=T, row.names=1)            ## 89632
data_fpkm = read.table("eXpress_isoform_matrix.TMM.fpkm.matrix", header=T, row.names=1)

dir.create("All_DE_dir/isoform_edgeR", recursive = TRUE)
setwd("All_DE_dir/isoform_edgeR")

## re-write the input file to correct the table header (to fix the issue of columns starting with no)
write.table(data_count, file='eXpress_isoform_matrix.counts.matrix', sep='\t', quote=F, row.names=T, col.names=NA)
write.table(data_fpkm, file='eXpress_isoform_matrix.TMM.fpkm.matrix', sep='\t', quote=F, row.names=T, col.names=NA)

## edgeR expects integer values so convert the float no into intgers
round_count=round(data_count)
## exclude isoforms with low expression level
round_count_sel = round_count[rowMax(as.matrix(data_fpkm))>=10,]

## initiate the edgeR DGEList
y <- DGEList(counts=round_count_sel)

## insert the group names according to columns names
NoIrr_1h_id=which(names(data_count)%in%c("X1.express", "X2.express", "X3.express"))
#NoIrr_6h_id=which(names(data_count)%in%c("X5.express", "X6.express"))
#NoIrr_12h_id=which(names(data_count)%in%c("X7.express", "X9.express"))
NoIrr_6h_id=which(names(data_count)%in%c("X4.express", "X5.express", "X6.express"))
NoIrr_12h_id=which(names(data_count)%in%c("X7.express"))
Irr_1h_id=which(names(data_count)%in%c("X10.express", "X11.express", "X12.express"))
Irr_6h_id=which(names(data_count)%in%c("X13.express", "X14.express", "X15.express"))
Irr_12h_id=which(names(data_count)%in%c("X16.express", "X17.express"))

levels(y$samples$group)=c("NoIrr_1h","NoIrr_6h","NoIrr_12h","Irr_1h","Irr_6h","Irr_12h")
y$samples$group[NoIrr_1h_id]=c(rep("NoIrr_1h",3))
#y$samples$group[NoIrr_6h_id]=c(rep("NoIrr_6h",2))
#y$samples$group[NoIrr_12h_id]=c(rep("NoIrr_12h",2))
y$samples$group[NoIrr_6h_id]=c(rep("NoIrr_6h",3))
y$samples$group[NoIrr_12h_id]=c(rep("NoIrr_12h",1))
y$samples$group[Irr_1h_id]=c(rep("Irr_1h",3))
y$samples$group[Irr_6h_id]=c(rep("Irr_6h",3))
y$samples$group[Irr_12h_id]=c(rep("Irr_12h",2))
levels(y$samples$group)

## create design matrix for your expermint
design <- model.matrix(~group, data=y$samples)
design

## normalization
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)

## mds plots
logcpm <- cpm(y, prior.count=2, log=TRUE)
pdf("mdsplot_pairwise.pdf")
plotMDS(logcpm, top = 1000, labels=y$samples$group)
dev.off()

pdf("mdsplot_common.pdf")
plotMDS(logcpm, top = 1000, labels=y$samples$group, gene.selection = "common")
dev.off()

## pairwise comparisons
for(i in 1:5)
{
    y2 <- y
    y2$samples$group <- relevel(y2$samples$group, ref=levels(y$samples$group)[i])
    print(levels(y2$samples$group))
    design2 <- model.matrix(~group, data=y2$samples)
    fit <- glmQLFit(y2,design2)
    datasetA=levels(y2$samples$group)[1]
    for(z in (i+1):6)
    {
        lrt <- glmTreat(fit,coef=z,lfc=1)
        FD <- topTags(lrt, n=NULL)
        datasetB=levels(y$samples$group)[z]
        table.name = paste("eXpress_isoform_matrix.",datasetA,"_vs_",datasetB,".edgeR.DE_results",sep="")
        print(table.name)
        write.table(FD, file=table.name, sep='\t', quote=F, row.names=T)
        FDf=FD$table
        
        logFC_cutoff=2
        FDR_cutoff=0.05
        plotfile = paste("eXpress_isoform_matrix.",datasetA,"_vs_",datasetB,".MA_n_Volcano_P",FDR_cutoff,"_C",logFC_cutoff,".pdf",sep="")
        pdf(plotfile)
        ## I edit the script to select for target FDR and absolute log fold change
        plot_MA_and_Volcano(FDf$logCPM, FDf$logFC, FDf$FDR, logFC_cutoff, FDR_cutoff)
        dev.off()
        
        logFC_cutoff=2
        FDR_cutoff=0.001
        plotfile = paste("eXpress_isoform_matrix.",datasetA,"_vs_",datasetB,".MA_n_Volcano_P",FDR_cutoff,"_C",logFC_cutoff,".pdf",sep="")
        pdf(plotfile)
        plot_MA_and_Volcano(FDf$logCPM, FDf$logFC, FDf$FDR, logFC_cutoff, FDR_cutoff)
        dev.off()
        
        logFC_cutoff=3
        FDR_cutoff=0.001
        plotfile = paste("eXpress_isoform_matrix.",datasetA,"_vs_",datasetB,".MA_n_Volcano_P",FDR_cutoff,"_C",logFC_cutoff,".pdf",sep="")
        pdf(plotfile)
        plot_MA_and_Volcano(FDf$logCPM, FDf$logFC, FDf$FDR, logFC_cutoff, FDR_cutoff)
        dev.off()
    }
}

dir.create("ANOVA")
setwd("ANOVA")

## ANOVA for effect of sunrise
fit <- glmQLFit(y,design)
lrt <- glmQLFTest(fit,coef=2:3)
FD <- topTags(lrt, n=NULL)
datasetA=levels(y$samples$group)[1]
datasetB=levels(y$samples$group)[2]
datasetC=levels(y$samples$group)[3]
table.name = paste("eXpress_isoform_matrix.",datasetA,"_vs_",datasetB,"_vs_",datasetC,".edgeR.DE_results",sep="")
print(table.name)
write.table(FD, file=table.name, sep='\t', quote=F, row.names=T)

## anova for effect of Irrgation and sunrise
y2 <- y
y2$samples$group <- relevel(y2$samples$group, ref=levels(y$samples$group)[4])
levels(y2$samples$group)
design2 <- model.matrix(~group, data=y2$samples)
fit <- glmQLFit(y2,design2)
lrt <- glmQLFTest(fit,coef=5:6)
FD <- topTags(lrt, n=NULL)
datasetA=levels(y2$samples$group)[1]
datasetB=levels(y2$samples$group)[5]
datasetC=levels(y2$samples$group)[6]
table.name = paste("eXpress_isoform_matrix.",datasetA,"_vs_",datasetB,"_vs_",datasetC,".edgeR.DE_results",sep="")
print(table.name)
write.table(FD, file=table.name, sep='\t', quote=F, row.names=T)


## model for effect of irregation alone
dir.create("../interaction")
setwd("../interaction")

design <- model.matrix(~0+group, data=y$samples)
design
colnames(design) <- levels(y$samples$group)
design

fit <- glmQLFit(y,design)
my.contrasts <- makeContrasts(
 NoIrr_6vs1 = NoIrr_6h-NoIrr_1h,
 NoIrr_12vs1 = NoIrr_12h-NoIrr_1h,
 Irr_6vs1 = Irr_6h-Irr_1h,
 Irr_12vs1 = Irr_12h-Irr_1h,
 IrrEffect_1h = Irr_1h-NoIrr_1h,
 IrrEffect_6h = (Irr_6h-Irr_1h)-(NoIrr_6h-NoIrr_1h),
 IrrEffect_12h = (Irr_12h-Irr_1h)-(NoIrr_12h-NoIrr_1h),
 levels=design)
lrt <- glmTreat(fit,contrast=my.contrasts[,"IrrEffect_1h"],lfc=1)
FD <- topTags(lrt, n=NULL)
table.name = paste("eXpress_isoform_matrix.","IrrEffect_1h",".edgeR.DE_results",sep="")
print(table.name)
write.table(FD, file=table.name, sep='\t', quote=F, row.names=T)

lrt <- glmTreat(fit, contrast=my.contrasts[,"IrrEffect_6h"],lfc=1)
FD <- topTags(lrt, n=NULL)
table.name = paste("eXpress_isoform_matrix.","IrrEffect_6h",".edgeR.DE_results",sep="")
print(table.name)
write.table(FD, file=table.name, sep='\t', quote=F, row.names=T)

lrt <- glmTreat(fit, contrast=my.contrasts[,"IrrEffect_12h"],lfc=1)
FD <- topTags(lrt, n=NULL)
table.name = paste("eXpress_isoform_matrix.","IrrEffect_12h",".edgeR.DE_results",sep="")
print(table.name)
write.table(FD, file=table.name, sep='\t', quote=F, row.names=T)




























