
source("F:/wangj/Desktop/temp/brain asymmetry/functions/set_theme.R")

library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)


symmetry_gene = read.table("F:/wangj/Desktop/temp/brain asymmetry/data/HGNC_method/all_GENE_in_brain_clump_instance2.txt",header=FALSE)[,1]

ac_GO<-enrichGO(gene=symmetry_gene, OrgDb = org.Hs.eg.db,ont='ALL',keyType = 'ENSEMBL')
ac_GO_BP<-enrichGO(gene=symmetry_gene, OrgDb = org.Hs.eg.db,ont = 'BP',keyType = 'ENSEMBL')
ac_GO_MF<-enrichGO(gene=symmetry_gene, OrgDb = org.Hs.eg.db,ont = 'MF',keyType = 'ENSEMBL')
ac_GO_CC<-enrichGO(gene=symmetry_gene, OrgDb = org.Hs.eg.db,ont = 'CC',keyType = 'ENSEMBL')

save(ac_GO,file="F:/wangj/Desktop/temp/brain asymmetry/data/ac_GO.RData")
save(ac_GO_BP,file="F:/wangj/Desktop/temp/brain asymmetry/data/ac_GO_BP.RData")
save(ac_GO_MF,file="F:/wangj/Desktop/temp/brain asymmetry/data/ac_GO_MF.RData")
save(ac_GO_CC,file="F:/wangj/Desktop/temp/brain asymmetry/data/ac_GO_CC.RData")

#

dataTemp<-data.frame(Description=ac_GO_BP$Description,
	p.adjust=ac_GO_BP$p.adjust,
	Count=ac_GO_BP$Count,
	GeneRatio=sapply(ac_GO_BP$GeneRatio,function(x){
		temp<-as.numeric(unlist(strsplit(x,"/")))
		return(temp[1]/temp[2])
	})
)

dataTemp<-dataTemp[order(dataTemp$GeneRatio),]
dataTemp$Description<-factor(dataTemp$Description,level=dataTemp$Description)

ggplot(dataTemp[(dim(dataTemp)[1]-9):dim(dataTemp)[1],],aes(x=GeneRatio,y=Description,color=p.adjust,size=Count))+
geom_point(show.legend=T)+
scale_color_gradient(low = "red", high = "blue")+
scale_size_area(max_size=2)+
scale_x_continuous(breaks=c(0.05,0.1))+
theme(axis.text.y=element_text(angle=15,hjust=1),axis.text=element_text(size=rel(6/11)))+
labs(x="Gene ratio",y="GO",title="Biological process")


ggsave("F:/wangj/Desktop/temp/brain asymmetry/figs/GO_enrich_clusterProfiler_BP.pdf",width=unit(3.89,"in"),height=unit(1.5,"in"))

#


dataTemp<-data.frame(Description=ac_GO_CC$Description,
	p.adjust=ac_GO_CC$p.adjust,
	Count=ac_GO_CC$Count,
	GeneRatio=sapply(ac_GO_CC$GeneRatio,function(x){
		temp<-as.numeric(unlist(strsplit(x,"/")))
		return(temp[1]/temp[2])
	})
)

dataTemp<-dataTemp[order(dataTemp$GeneRatio),]
dataTemp$Description<-factor(dataTemp$Description,level=dataTemp$Description)


ggplot(dataTemp[(dim(dataTemp)[1]-9):dim(dataTemp)[1],],aes(x=GeneRatio,y=Description,color=p.adjust,size=Count))+
geom_point(show.legend=T)+
scale_color_gradient(low = "red", high = "blue")+
scale_size_area(max_size=2)+
scale_x_continuous(breaks=c(0.05,0.08))+
theme(axis.text.y=element_text(angle=15,hjust=1),axis.text=element_text(size=rel(6/11)))+
labs(x="Gene ratio",y="GO",title="Cellular component")


ggsave("F:/wangj/Desktop/temp/brain asymmetry/figs/GO_enrich_clusterProfiler_CC.pdf",width=unit(3.04,"in"),height=unit(1.5,"in"))



###

DEG_enrich<-read.table("F:/wangj/Desktop/temp/brain asymmetry/data/HGNC_method/FUMA_gene2func105843/gtex_v8_ts_DEG.txt",header=T,fill=T)

dataTemp<-DEG_enrich

brain_tissue<-setdiff(names(table(DEG_enrich$GeneSet))[grepl("Brain",names(table(DEG_enrich$GeneSet)))],"Brain_Spinal_cord_cervical_c-1")
non_brain_tissue<-setdiff(names(table(DEG_enrich$GeneSet)),brain_tissue)

dataTemp$GeneSet<-factor(dataTemp$GeneSet,level=c(brain_tissue,non_brain_tissue))

dataTemp$color1<-factor(ifelse(dataTemp$adjP<0.05,"red","blue"),level=c("red","blue"))

#

ggplot(dataTemp[dataTemp$Category=="DEG.up",],aes(x=GeneSet,y=-log10(p)))+
geom_bar(stat="identity",width=0.8,fill=dataTemp[dataTemp$Category=="DEG.up",]$color1)+
scale_y_continuous(breaks=c(0,4,8))+
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),axis.text=element_text(size=rel(6/11)))+
labs(x="Tissue",y=expression(paste("-",log[10],"(",italic(p),")")),title="")

ggsave("F:/wangj/Desktop/temp/brain asymmetry/figs/DEG_enrich_up.pdf",width=unit(5.8,"in"),height=unit(2.8,"in"))

#

ggplot(dataTemp[dataTemp$Category=="DEG.down",],aes(x=GeneSet,y=-log10(p)))+
geom_bar(stat="identity",width=0.8,fill=dataTemp[dataTemp$Category=="DEG.down",]$color1)+
scale_y_continuous(breaks=c(0,5,10))+
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),axis.text=element_text(size=rel(6/11)))+
labs(x="Tissue",y=expression(paste("-",log[10],"(",italic(p),")")),title="")

ggsave("F:/wangj/Desktop/temp/brain asymmetry/figs/DEG_enrich_down.pdf",width=unit(5.8,"in"),height=unit(2.8,"in"))



###


gene_tiss_exp<-read.table("F:/wangj/Desktop/temp/brain asymmetry/data/HGNC_method/FUMA_gene2func105843/gtex_v8_ts_avg_log2TPM_exp.txt",header=T)

#GTEx_tissue_name<-sapply(names(gene_tiss_exp),function(x){paste0(unlist(strsplit(x,"_")),collapse=" ")})
#write.csv(GTEx_tissue_name,"F:/wangj/Desktop/temp/brain asymmetry/MS/brain_paper/whole gene fuma/GTEx_tissue_name.csv")
#brain_tissue<-setdiff(names(gene_tiss_exp)[grepl("Brain",names(gene_tiss_exp))],"Brain_Spinal_cord_cervical_c.1")

library(pheatmap)
library(extrafont)
library(ggplot2)

#font_import()
loadfonts()



dataTemp<-gene_tiss_exp[,match(gsub("-",".",c(brain_tissue,non_brain_tissue)),names(gene_tiss_exp))]
rownames(dataTemp)<-gene_tiss_exp[,1]

tissue_group<-data.frame(Group=c(rep("Brain",length(brain_tissue)),rep("Non_brain",length(non_brain_tissue))))
rownames(tissue_group)<-gsub("-",".",c(brain_tissue,non_brain_tissue))

annotation_gene<-rep("Insignificant",dim(gene_tiss_exp)[1])
annotation_gene[sapply(1:dim(gene_tiss_exp)[1],function(i){t.test(dataTemp[i,1:12],gene_tiss_exp[i,13:54],alternative="less")$p.value})<0.05/dim(gene_tiss_exp)[1]]<-"Down_regulated"
annotation_gene[sapply(1:dim(gene_tiss_exp)[1],function(i){t.test(dataTemp[i,1:12],gene_tiss_exp[i,13:54],alternative="greater")$p.value})<0.05/dim(gene_tiss_exp)[1]]<-"Up_regulated"
annotation_gene<-data.frame(Type=annotation_gene)
rownames(annotation_gene)<-rownames(dataTemp)

annotation_colors_1<-list(Group=c(Brain="#800080",Non_brain="#87CCEB"),Type=c(Down_regulated="#00BFC4",Insignificant="yellow",Up_regulated="#F8766D"))

pdf("F:/wangj/Desktop/temp/brain asymmetry/figs/gene_tiss_exp.pdf",width=7.36,height=4)

pheatmap(dataTemp,cluster_rows=T,cluster_cols=F,angle_col=315,fontsize_col=unit(5,"pt"),treeheight_row=10,
color=colorRampPalette(c("navy","white","firebrick3"))(100),
show_colnames=T,border_color=NA,scale="none",
show_rownames=F,annotation_row=annotation_gene,annotation_col=tissue_group,annotation_colors=annotation_colors_1)

dev.off()

write.csv(data.frame(gene_tiss_exp[,1:2],annotation_gene),"F:/wangj/Desktop/temp/brain asymmetry/data/gene_up_down.csv",row.names=F)





















