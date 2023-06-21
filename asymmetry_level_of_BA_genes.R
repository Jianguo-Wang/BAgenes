

source("F:/wangj/Desktop/temp/brain asymmetry/functions/set_theme.R")
source("F:/wangj/Desktop/temp/brain asymmetry/functions/rowColRecovery.R")   

brain_asymmetry_gene<-read.table("F:/wangj/Desktop/temp/brain asymmetry/data/HGNC_method/all_GENE_in_brain_clump_instance2.txt",header=FALSE)[,1]
human_gene_list<-read.csv("F:/wangj/Desktop/temp/brain asymmetry/data/hgnc_complete_set_2023-03-02.csv")
#lines <- readLines("F:/wangj/Desktop/temp/brain asymmetry/data/hgnc_complete_set_2023-03-02.txt")


### gene expression data and annotation files
# donor 1 

d1_Genes<-read.csv("F:/wangj/Desktop/temp/brain asymmetry/data/rnaseq_donor9861/Genes.csv")
d1_Ontology<-read.csv("F:/wangj/Desktop/temp/brain asymmetry/data/rnaseq_donor9861/Ontology.csv")
d1_TPM<-read.csv("F:/wangj/Desktop/temp/brain asymmetry/data/rnaseq_donor9861/RNAseqTPM.csv",header=FALSE)
d1_annot<-read.csv("F:/wangj/Desktop/temp/brain asymmetry/data/rnaseq_donor9861/SampleAnnot.csv")
d1_TPM_rownames<-d1_Genes$gene_symbol
d1_TPM_colnames<-d1_annot$RNAseq_sample_name
d1_TPM<-d1_TPM[,-1]
rownames(d1_TPM)<-d1_TPM_rownames
names(d1_TPM)<-d1_TPM_colnames

# donor 2

d2_Genes<-read.csv("F:/wangj/Desktop/temp/brain asymmetry/data/rnaseq_donor10021/Genes.csv")
d2_Ontology<-read.csv("F:/wangj/Desktop/temp/brain asymmetry/data/rnaseq_donor10021/Ontology.csv")
d2_TPM<-read.csv("F:/wangj/Desktop/temp/brain asymmetry/data/rnaseq_donor10021/RNAseqTPM.csv",header=FALSE)
d2_annot<-read.csv("F:/wangj/Desktop/temp/brain asymmetry/data/rnaseq_donor10021/SampleAnnot.csv")
d2_TPM_rownames<-d2_Genes$gene_symbol
d2_TPM_colnames<-d2_annot$RNAseq_sample_name
d2_TPM<-d2_TPM[,-1]
rownames(d2_TPM)<-d2_TPM_rownames
names(d2_TPM)<-d2_TPM_colnames


### overlapped genes
#174 genes overlapped
brain_overlap_gene<-human_gene_list$entrez_id[na.omit(match(brain_asymmetry_gene,human_gene_list$ensembl_gene_id))]
brain_overlap_ind<-sort(na.omit(match(brain_overlap_gene,d1_Genes$entrez_id)))


### samples including both L and R regions

# donor 1 

d1_ontology_structure_acronym_table<-table(d1_annot$ontology_structure_acronym)
d1_ontology_structure_acronym_table_names<-names(d1_ontology_structure_acronym_table)[which(d1_ontology_structure_acronym_table>1)]
d1_ontology_structure_acronym_table_names_LR<-d1_ontology_structure_acronym_table_names[sapply(d1_ontology_structure_acronym_table_names,function(x){
	ifelse(length(unique(d1_annot$hemisphere[d1_annot$ontology_structure_acronym==x]))==2,TRUE,FALSE)
})]

d1_LR_sample_ind<-lapply(d1_ontology_structure_acronym_table_names_LR,function(x){
	list(L=which(d1_annot$ontology_structure_acronym==x&d1_annot$hemisphere=="L"),R=which(d1_annot$ontology_structure_acronym==x&d1_annot$hemisphere=="R"))
})

# donor 2

d2_ontology_structure_acronym_table<-table(d2_annot$ontology_structure_acronym)
d2_ontology_structure_acronym_table_names<-names(d2_ontology_structure_acronym_table)[which(d2_ontology_structure_acronym_table>1)]
d2_ontology_structure_acronym_table_names_LR<-d2_ontology_structure_acronym_table_names[sapply(d2_ontology_structure_acronym_table_names,function(x){
	ifelse(length(unique(d2_annot$hemisphere[d2_annot$ontology_structure_acronym==x]))==2,TRUE,FALSE)
})]

d2_LR_sample_ind<-lapply(d2_ontology_structure_acronym_table_names_LR,function(x){
	list(L=which(d2_annot$ontology_structure_acronym==x&d2_annot$hemisphere=="L"),R=which(d2_annot$ontology_structure_acronym==x&d2_annot$hemisphere=="R"))
})

### cal AL 

## function 

cal_FA<-function(dat,genInd,LR_indList,cutoff=0){

	L_values<-dat[genInd,LR_indList[[1]]]
	R_values<-dat[genInd,LR_indList[[2]]]
	L_values[which(L_values<cutoff)]<-NA
	R_values[which(R_values<cutoff)]<-NA
	L=mean(L_values,na.rm=T)
	R=mean(R_values,na.rm=T)
	
	return((L-R)/sqrt(2*(L^2+R^2)))
}


## donor 1
# withour cutoff 

d1_AL<-sapply(1:dim(d1_TPM)[1],function(i){
	print(i)
	return(sapply(1:length(d1_LR_sample_ind),function(j){cal_FA(d1_TPM,i,d1_LR_sample_ind[[j]])}))
})

colnames(d1_AL)<-rownames(d1_TPM)
rownames(d1_AL)<-d1_ontology_structure_acronym_table_names_LR

# with cutoff 

d1_AL_filter<-sapply(1:dim(d1_TPM)[1],function(i){
	print(i)
	return(sapply(1:length(d1_LR_sample_ind),function(j){cal_FA(d1_TPM,i,d1_LR_sample_ind[[j]],cutoff=1e-7)}))
})


colnames(d1_AL_filter)<-rownames(d1_TPM)
rownames(d1_AL_filter)<-d1_ontology_structure_acronym_table_names_LR

length(which(as.matrix(d1_TPM)<1e-7))/length(as.matrix(d1_TPM))


## donor 2
# withour cutoff 

d2_AL<-sapply(1:dim(d2_TPM)[1],function(i){
	print(i)
	return(sapply(1:length(d2_LR_sample_ind),function(j){cal_FA(d2_TPM,i,d2_LR_sample_ind[[j]])}))
})

colnames(d2_AL)<-rownames(d2_TPM)
rownames(d2_AL)<-d2_ontology_structure_acronym_table_names_LR

# with cutoff 

d2_AL_filter<-sapply(1:dim(d2_TPM)[1],function(i){
	print(i)
	return(sapply(1:length(d2_LR_sample_ind),function(j){cal_FA(d2_TPM,i,d2_LR_sample_ind[[j]],cutoff=1e-7)}))
})

colnames(d2_AL_filter)<-rownames(d2_TPM)
rownames(d2_AL_filter)<-d2_ontology_structure_acronym_table_names_LR

length(which(as.matrix(d2_TPM)<1e-7))/length(as.matrix(d2_TPM))


#### compare between BA genes and other genes 

genetype<-rep("other genes",dim(d1_TPM)[1])
genetype[brain_overlap_ind]<-"brain asymmetry"


d1_AL_df<-data.frame(AL=c(d1_AL),AL_filter=c(d1_AL_filter),type=c(sapply(genetype,rep,dim(d1_AL)[1])),gene=c(sapply(colnames(d1_AL),rep,dim(d1_AL)[1])),region=rep(rownames(d1_AL),dim(d1_AL)[2]))
d2_AL_df<-data.frame(AL=c(d2_AL),AL_filter=c(d2_AL_filter),type=c(sapply(genetype,rep,dim(d2_AL)[1])),gene=c(sapply(colnames(d2_AL),rep,dim(d2_AL)[1])),region=rep(rownames(d2_AL),dim(d2_AL)[2]))

d1_AL_df_nonNA<-d1_AL_df[!is.na(d1_AL_df$AL)&d1_AL_df$AL!=Inf&d1_AL_df$AL!=-Inf,]
d2_AL_df_nonNA<-d2_AL_df[!is.na(d2_AL_df$AL)&d2_AL_df$AL!=Inf&d2_AL_df$AL!=-Inf,]


common_d1_d2_nonNA_identifier<-intersect(apply(d1_AL_df_nonNA[,4:5],1,paste0,collapse="_"),apply(d2_AL_df_nonNA[,4:5],1,paste0,collapse="_"))
d1_AL_df_nonNA_com<-d1_AL_df_nonNA[match(common_d1_d2_nonNA_identifier,apply(d1_AL_df_nonNA[,4:5],1,paste0,collapse="_")),]
d2_AL_df_nonNA_com<-d2_AL_df_nonNA[match(common_d1_d2_nonNA_identifier,apply(d2_AL_df_nonNA[,4:5],1,paste0,collapse="_")),]

cor.test(d1_AL_df_nonNA_com[,1],d2_AL_df_nonNA_com[,1])
length(which(sign(d1_AL_df_nonNA_com[,1]*d2_AL_df_nonNA_com[,1])==1))/length(d1_AL_df_nonNA_com[,1])

cor.test(d1_AL_df_nonNA_com[,2],d2_AL_df_nonNA_com[,2],use="pairwise.complete.obs")
length(which(sign(d1_AL_df_nonNA_com[,2]*d2_AL_df_nonNA_com[,2])==1))/length(which(!is.na(d1_AL_df_nonNA_com[,2])&!is.na(d2_AL_df_nonNA_com[,2])))


### donor 1

## without cutoff
#total distribution

ggplot(d1_AL_df_nonNA,aes(x=abs(AL),col=type,fill=type))+geom_density(alpha=0.5,show.legend=FALSE)+
labs(x="AL of genes for donor 1",y="Density",title="")

wilcox.test(abs(d1_AL_df_nonNA[which(d1_AL_df_nonNA$type=="other genes"),1]),abs(d1_AL_df_nonNA[which(d1_AL_df_nonNA$type=="brain asymmetry"),1]),alternative="less")

# respective regions 

sapply(names(table(d1_AL_df_nonNA$region)),function(x){
	temp<-d1_AL_df_nonNA[d1_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),1]),abs(temp[which(temp$type=="brain asymmetry"),1]),alternative="less")$p.value)
	
})

## with cutoff 
#total distribution

ggplot(d1_AL_df_nonNA,aes(x=abs(AL_filter),col=type,fill=type))+geom_density(alpha=0.5,show.legend=FALSE)+
labs(x="AL of genes for donor 1",y="Density",title="")

ggsave("F:/wangj/Desktop/temp/brain asymmetry/figs/d1_AL_nonNA.pdf",width=unit(3,"in"),height=unit(3.1,"in"))


wilcox.test(abs(d1_AL_df_nonNA[which(d1_AL_df_nonNA$type=="other genes"),2]),abs(d1_AL_df_nonNA[which(d1_AL_df_nonNA$type=="brain asymmetry"),2]),alternative="less")

# respective regions 

sapply(names(table(d1_AL_df_nonNA$region)),function(x){
	temp<-d1_AL_df_nonNA[d1_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),2]),abs(temp[which(temp$type=="brain asymmetry"),2]),alternative="less")$p.value)
	
})

#multiple testing 

p.adjust(sapply(names(table(d1_AL_df_nonNA$region)),function(x){
	temp<-d1_AL_df_nonNA[d1_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),2]),abs(temp[which(temp$type=="brain asymmetry"),2]),alternative="less")$p.value)
	
}),method="BH")


### donor 2
## without cutoff
#total distribution

ggplot(d2_AL_df_nonNA,aes(x=abs(AL),col=type,fill=type))+geom_density(alpha=0.5,show.legend=FALSE)+
labs(x="AL of genes for donor 2",y="Density",title="")

wilcox.test(abs(d2_AL_df_nonNA[which(d2_AL_df_nonNA$type=="other genes"),1]),abs(d2_AL_df_nonNA[which(d2_AL_df_nonNA$type=="brain asymmetry"),1]),alternative="less")

# respective regions 

sapply(names(table(d2_AL_df_nonNA$region)),function(x){
	temp<-d2_AL_df_nonNA[d2_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),1]),abs(temp[which(temp$type=="brain asymmetry"),1]),alternative="less")$p.value)
	
})

## with cutoff 
#total distribution

ggplot(d2_AL_df_nonNA,aes(x=abs(AL_filter),col=type,fill=type))+geom_density(alpha=0.5,show.legend=FALSE)+
labs(x="AL of genes for donor 2",y="Density",title="")

ggsave("F:/wangj/Desktop/temp/brain asymmetry/figs/d2_AL_nonNA.pdf",width=unit(3,"in"),height=unit(3.1,"in"))

wilcox.test(abs(d2_AL_df_nonNA[which(d2_AL_df_nonNA$type=="other genes"),2]),abs(d2_AL_df_nonNA[which(d2_AL_df_nonNA$type=="brain asymmetry"),2]),alternative="less")

# respective regions 

sapply(names(table(d2_AL_df_nonNA$region)),function(x){
	temp<-d2_AL_df_nonNA[d2_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),2]),abs(temp[which(temp$type=="brain asymmetry"),2]),alternative="less")$p.value)
	
})

#multiple testing 

p.adjust(sapply(names(table(d2_AL_df_nonNA$region)),function(x){
	temp<-d2_AL_df_nonNA[d2_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),2]),abs(temp[which(temp$type=="brain asymmetry"),2]),alternative="less")$p.value)
	
}),method="BH")


### average of donors

com_AL_df_nonNA<-data.frame(AL=(d1_AL_df_nonNA_com[,1]+d2_AL_df_nonNA_com[,1])/2,
	AL_filter=(d1_AL_df_nonNA_com[,2]+d2_AL_df_nonNA_com[,2])/2)
com_AL_df_nonNA$type=d1_AL_df_nonNA_com$type 
com_AL_df_nonNA$gene=d1_AL_df_nonNA_com$gene
com_AL_df_nonNA$region=d1_AL_df_nonNA_com$region

## without cutoff 
#total distribution

ggplot(com_AL_df_nonNA,aes(x=abs(AL),col=type,fill=type))+geom_density(alpha=0.5,show.legend=FALSE)+
labs(x="Mean AL of genes",y="Density",title="")

ggsave("F:/wangj/Desktop/temp/brain asymmetry/figs/d1_d2_mean_AL_nonNA.pdf",width=unit(3,"in"),height=unit(3.1,"in"))

wilcox.test(abs(com_AL_df_nonNA[which(com_AL_df_nonNA$type=="other genes"),1]),abs(com_AL_df_nonNA[which(com_AL_df_nonNA$type=="brain asymmetry"),1]),alternative="less")

# respective regions 

sapply(names(table(com_AL_df_nonNA$region)),function(x){
	temp<-com_AL_df_nonNA[com_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),1]),abs(temp[which(temp$type=="brain asymmetry"),1]),alternative="less")$p.value)
	
})

## with cutoff 
#total distribution

ggplot(com_AL_df_nonNA,aes(x=abs(AL_filter),col=type,fill=type))+geom_density(alpha=0.5,show.legend=FALSE)+
labs(x="Mean AL of genes",y="Density",title="")

wilcox.test(abs(com_AL_df_nonNA[which(com_AL_df_nonNA$type=="other genes"),2]),abs(com_AL_df_nonNA[which(com_AL_df_nonNA$type=="brain asymmetry"),2]),alternative="less")

# respective regions 

sapply(names(table(com_AL_df_nonNA$region)),function(x){
	temp<-com_AL_df_nonNA[com_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),2]),abs(temp[which(temp$type=="brain asymmetry"),2]),alternative="less")$p.value)
	
})

#multiple testing 

p.adjust(sapply(names(table(com_AL_df_nonNA$region)),function(x){
	temp<-com_AL_df_nonNA[com_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),2]),abs(temp[which(temp$type=="brain asymmetry"),2]),alternative="less")$p.value)
	
}),method="BH")




##########################################

save.image("F:/wangj/Desktop/temp/brain asymmetry/figs/asymmetrical_gene_expression.RData")

d1_Ontology[match(c("LiG-pest","LOrG","ITG-mts","MTG-s"),d1_Ontology$acronym),]


### output supplementary tables 
#d1_AL_filter
#d2_AL_filter
#d1_AL_df_nonNA
#d2_AL_df_nonNA
#com_AL_df_nonNA

## pvalues for d1 

d1_p_total<-wilcox.test(abs(d1_AL_df_nonNA[which(d1_AL_df_nonNA$type=="other genes"),2]),abs(d1_AL_df_nonNA[which(d1_AL_df_nonNA$type=="brain asymmetry"),2]),alternative="less")$p.value

d1_p_region<-sapply(names(table(d1_AL_df_nonNA$region)),function(x){
	temp<-d1_AL_df_nonNA[d1_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),2]),abs(temp[which(temp$type=="brain asymmetry"),2]),alternative="less")$p.value)
	
})

d1_p_region_mltiple_testing<-p.adjust(sapply(names(table(d1_AL_df_nonNA$region)),function(x){
	temp<-d1_AL_df_nonNA[d1_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),2]),abs(temp[which(temp$type=="brain asymmetry"),2]),alternative="less")$p.value)
	
}),method="BH")

d1_region<-names(table(d1_AL_df_nonNA$region))

d1_pvalue<-data.frame(region=c("Total",d1_region),pvalue=c(d1_p_total,d1_p_region),p.adjust=c(NA,d1_p_region_mltiple_testing))


## pvalues for d2

d2_p_total<-wilcox.test(abs(d2_AL_df_nonNA[which(d2_AL_df_nonNA$type=="other genes"),2]),abs(d2_AL_df_nonNA[which(d2_AL_df_nonNA$type=="brain asymmetry"),2]),alternative="less")$p.value

d2_p_region<-sapply(names(table(d2_AL_df_nonNA$region)),function(x){
	temp<-d2_AL_df_nonNA[d2_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),2]),abs(temp[which(temp$type=="brain asymmetry"),2]),alternative="less")$p.value)
	
})

d2_p_region_mltiple_testing<-p.adjust(sapply(names(table(d2_AL_df_nonNA$region)),function(x){
	temp<-d2_AL_df_nonNA[d2_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),2]),abs(temp[which(temp$type=="brain asymmetry"),2]),alternative="less")$p.value)
	
}),method="BH")

d2_region<-names(table(d2_AL_df_nonNA$region))

d2_pvalue<-data.frame(region=c("Total",d2_region),pvalue=c(d2_p_total,d2_p_region),p.adjust=c(NA,d2_p_region_mltiple_testing))


## pvalues for average of d1 and d2

d1d2_p_total<-wilcox.test(abs(com_AL_df_nonNA[which(com_AL_df_nonNA$type=="other genes"),2]),abs(com_AL_df_nonNA[which(com_AL_df_nonNA$type=="brain asymmetry"),2]),alternative="less")$p.value

d1d2_p_region<-sapply(names(table(com_AL_df_nonNA$region)),function(x){
	temp<-com_AL_df_nonNA[com_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),2]),abs(temp[which(temp$type=="brain asymmetry"),2]),alternative="less")$p.value)
	
})

d1d2_p_region_mltiple_testing<-p.adjust(sapply(names(table(com_AL_df_nonNA$region)),function(x){
	temp<-com_AL_df_nonNA[com_AL_df_nonNA$region==x,]
	return(wilcox.test(abs(temp[which(temp$type=="other genes"),2]),abs(temp[which(temp$type=="brain asymmetry"),2]),alternative="less")$p.value)
	
}),method="BH")


d1d2_region<-names(table(com_AL_df_nonNA$region))

d1d2_pvalue<-data.frame(region=c("Total",d1d2_region),pvalue=c(d1d2_p_total,d1d2_p_region),p.adjust=c(NA,d1d2_p_region_mltiple_testing))


# Load the 'openxlsx' package
library(openxlsx)

# Create a new Excel workbook
wb <- createWorkbook()

# Add data frames to separate worksheets
addWorksheet(wb, "donor1_AL")
writeData(wb, sheet = 1, x = t(d1_AL_filter),rowNames=TRUE)

addWorksheet(wb, "donor2_AL")
writeData(wb, sheet = 2, x = t(d2_AL_filter),rowNames=TRUE)

addWorksheet(wb, "donor1_AL_non_NA")
writeData(wb, sheet = 3, x = d1_AL_df_nonNA[,-1])

addWorksheet(wb, "donor2_AL_non_NA")
writeData(wb, sheet = 4, x = d2_AL_df_nonNA[,-1])

addWorksheet(wb, "donor1_donor2_AL_common")
writeData(wb, sheet = 5, x = com_AL_df_nonNA[,-1])

addWorksheet(wb, "donor1_pvalue")
writeData(wb, sheet = 6, x = d1_pvalue,keepNA=TRUE)

addWorksheet(wb, "donor2_pvalue")
writeData(wb, sheet = 7, x = d2_pvalue,keepNA=TRUE)

addWorksheet(wb, "donors_average_pvalue")
writeData(wb, sheet = 8, x = d1d2_pvalue,keepNA=TRUE)

# Save the workbook to a file
saveWorkbook(wb, file = "F:/wangj/Desktop/temp/brain asymmetry/MS/Tab.S7-Asymmetry_Level_of_genes_for_two_donors.xlsx", overwrite = TRUE)


data<-list(t(d1_AL_filter),t(d2_AL_filter),d1_AL_df_nonNA[,-1],d2_AL_df_nonNA[,-1],com_AL_df_nonNA[,-1],d1_pvalue,d2_pvalue,d1d2_pvalue)
names(data)<-c("d1_AL_filter","d2_AL_filter","d1_AL_df_nonNA","d2_AL_df_nonNA","com_AL_df_nonNA","d1_pvalue","d2_pvalue","d1d2_pvalue")
saveRDS(data, file = "F:/wangj/Desktop/temp/brain asymmetry/MS/Tab.S7-Asymmetry_Level_of_genes_for_two_donors.rds")



