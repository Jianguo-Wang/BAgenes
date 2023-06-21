
source("F:/wangj/Desktop/temp/brain asymmetry/functions/set_theme.R")
source("F:/wangj/Desktop/temp/brain asymmetry/functions/rowColRecovery.R")   

setwd("F:/wangj/Desktop/temp/brain asymmetry/data/species-organ-expression")

human<-read.table("Human_rpkm.txt",header=T)
chicken<-read.table("Chicken_rpkm.txt",header=T)
macaque<-read.table("Macaque_rpkm.txt",header=T)
mouse<-read.table("Mouse_rpkm.txt",header=T)
opossum<-read.table("Opossum_rpkm.txt",header=T)
rabbit<-read.table("Rabbit_rpkm.txt",header=T)
rat<-read.table("Rat_rpkm.txt",header=T)

brain_asymmetry_gene<-read.table("F:/wangj/Desktop/temp/brain asymmetry/data/HGNC_method/all_GENE_in_brain_clump_instance2.txt",header=FALSE)[,1]

cerebellum_annot<-read.csv("F:/wangj/Desktop/temp/brain asymmetry/data/cerebellum_annot.csv")
cerebellum_gene<-setdiff(unique(union(cerebellum_annot$emsembl_phenoscanner,cerebellum_annot$Annovar_ensembl)),c("-",""))

cerebrum_annot<-read.csv("F:/wangj/Desktop/temp/brain asymmetry/data/cerebrum_annot.csv")
cerebrum_gene<-setdiff(unique(union(cerebrum_annot$emsembl_phenoscanner,cerebrum_annot$Annovar_ensembl)),c("-",""))

cerebellum_specific_gene<-setdiff(cerebellum_gene,cerebrum_gene)

brain_asymmetry_gene<-setdiff(brain_asymmetry_gene,cerebellum_specific_gene)

gene_tiss_exp<-read.table("F:/wangj/Desktop/temp/brain asymmetry/data/HGNC_method/FUMA_gene2func105843/gtex_v8_ts_avg_log2TPM_exp.txt",header=T)

stage_inf<-read.csv("table_development_stage_2_wjg.csv")

#

human_chicken<-read.table("humanVSchicken.txt",header=T,sep ="\t")
human_macaque<-read.table("humanVSmacaque.txt",header=T,sep ="\t")
human_mouse<-read.table("humanVSmouse.txt",header=T,sep ="\t")
human_opossum<-read.table("humanVSopossum.txt",header=T,sep ="\t")
human_rabbit<-read.table("humanVSrabbit.txt",header=T,sep ="\t")
human_rat<-read.table("humanVSrat.txt",header=T,sep ="\t")

get_1v1<-function(human_species){
	human_species_conf<-human_species[which(human_species[,9]==1),]
	human_species_conf_uniq<-human_species_conf[intersect(intersect(match(unique(human_species_conf[,3]),human_species_conf[,3]),match(unique(human_species_conf[,1]),human_species_conf[,1])),match(unique(human_species_conf[,2]),human_species_conf[,2])),]
	return(human_species_conf_uniq)
}

human_macaque_1v1<-get_1v1(human_macaque)
human_mouse_1v1<-get_1v1(human_mouse)
human_opossum_1v1<-get_1v1(human_opossum)
human_rabbit_1v1<-get_1v1(human_rabbit)
human_rat_1v1<-get_1v1(human_rat)


#

human_ensembl_union<-sort(unique(c(human_macaque_1v1[,1],human_mouse_1v1[,1],human_opossum_1v1[,1],human_rabbit_1v1[,1],human_rat_1v1[,1])))

human_ensembl_correspond<-data.frame(human_ensembl=human_ensembl_union)
human_ensembl_correspond$macaque_ensembl<-rep(NA,length(human_ensembl_union))
human_ensembl_correspond$macaque_ensembl[which(!is.na(match(human_ensembl_union,human_macaque_1v1[,1])))]<-human_macaque_1v1[na.omit(match(human_ensembl_union,human_macaque_1v1[,1])),3]

human_ensembl_correspond$mouse_ensembl<-rep(NA,length(human_ensembl_union))
human_ensembl_correspond$mouse_ensembl[which(!is.na(match(human_ensembl_union,human_mouse_1v1[,1])))]<-human_mouse_1v1[na.omit(match(human_ensembl_union,human_mouse_1v1[,1])),3]

human_ensembl_correspond$rat_ensembl<-rep(NA,length(human_ensembl_union))
human_ensembl_correspond$rat_ensembl[which(!is.na(match(human_ensembl_union,human_rat_1v1[,1])))]<-human_rat_1v1[na.omit(match(human_ensembl_union,human_rat_1v1[,1])),3]

human_ensembl_correspond$opossum_ensembl<-rep(NA,length(human_ensembl_union))
human_ensembl_correspond$opossum_ensembl[which(!is.na(match(human_ensembl_union,human_opossum_1v1[,1])))]<-human_opossum_1v1[na.omit(match(human_ensembl_union,human_opossum_1v1[,1])),3]

human_ensembl_correspond$rabbit_ensembl<-rep(NA,length(human_ensembl_union))
human_ensembl_correspond$rabbit_ensembl[which(!is.na(match(human_ensembl_union,human_rabbit_1v1[,1])))]<-human_rabbit_1v1[na.omit(match(human_ensembl_union,human_rabbit_1v1[,1])),3]


#

extract_exp_correspond<-function(organ,stage){
	
	temp<-data.frame(human=rep(NA,dim(human_ensembl_correspond)[1]),
		macaque=rep(NA,dim(human_ensembl_correspond)[1]),
		mouse=rep(NA,dim(human_ensembl_correspond)[1]),
		rat=rep(NA,dim(human_ensembl_correspond)[1]),
		opossum=rep(NA,dim(human_ensembl_correspond)[1]),
		rabbit=rep(NA,dim(human_ensembl_correspond)[1]))
	
	human_ind<-which(grepl(paste0(organ,"\\.",stage_inf[match(stage,stage_inf$stage),"human"],"\\."),names(human)))
	macaque_ind<-which(grepl(paste0(organ,"\\.",stage_inf[match(stage,stage_inf$stage),"macaque"],"\\."),names(macaque)))
	mouse_ind<-which(grepl(paste0(organ,"\\.",stage_inf[match(stage,stage_inf$stage),"mouse"],"\\."),names(mouse)))
	rat_ind<-which(grepl(paste0(organ,"\\.",stage_inf[match(stage,stage_inf$stage),"rat"],"\\."),names(rat)))
	opossum_ind<-which(grepl(paste0(organ,"\\.",stage_inf[match(stage,stage_inf$stage),"opossum"],"\\."),names(opossum)))
	rabbit_ind<-which(grepl(paste0(organ,"\\.",stage_inf[match(stage,stage_inf$stage),"rabbit"],"\\."),names(rabbit)))
	
	if(length(human_ind)>1){
		temp$human[which(!is.na(match(human_ensembl_correspond$human_ensembl,human[,1])))]<-apply(human[na.omit(match(human_ensembl_correspond$human_ensembl,human[,1])),human_ind],1,mean,na.rm=T)
	}
	if(length(human_ind)==1){
		temp$human[which(!is.na(match(human_ensembl_correspond$human_ensembl,human[,1])))]<-human[na.omit(match(human_ensembl_correspond$human_ensembl,human[,1])),human_ind]
	}

	if(length(macaque_ind)>1){
		temp$macaque[which(!is.na(match(human_ensembl_correspond$macaque_ensembl,macaque[,1])))]<-apply(macaque[na.omit(match(human_ensembl_correspond$macaque_ensembl,macaque[,1])),macaque_ind],1,mean,na.rm=T)
	}
	if(length(macaque_ind)==1){
		temp$macaque[which(!is.na(match(human_ensembl_correspond$macaque_ensembl,macaque[,1])))]<-macaque[na.omit(match(human_ensembl_correspond$macaque_ensembl,macaque[,1])),macaque_ind]
	}
	
	if(length(mouse_ind)>1){
		temp$mouse[which(!is.na(match(human_ensembl_correspond$mouse_ensembl,mouse[,1])))]<-apply(mouse[na.omit(match(human_ensembl_correspond$mouse_ensembl,mouse[,1])),mouse_ind],1,mean,na.rm=T)
	}
	if(length(mouse_ind)==1){
		temp$mouse[which(!is.na(match(human_ensembl_correspond$mouse_ensembl,mouse[,1])))]<-mouse[na.omit(match(human_ensembl_correspond$mouse_ensembl,mouse[,1])),mouse_ind]
	}
	
	if(length(rat_ind)>1){
		temp$rat[which(!is.na(match(human_ensembl_correspond$rat_ensembl,rat[,1])))]<-apply(rat[na.omit(match(human_ensembl_correspond$rat_ensembl,rat[,1])),rat_ind],1,mean,na.rm=T)	
	}
	if(length(rat_ind)==1){
		temp$rat[which(!is.na(match(human_ensembl_correspond$rat_ensembl,rat[,1])))]<-rat[na.omit(match(human_ensembl_correspond$rat_ensembl,rat[,1])),rat_ind]
	}
	
	if(length(opossum_ind)>1){
		temp$opossum[which(!is.na(match(human_ensembl_correspond$opossum_ensembl,opossum[,1])))]<-apply(opossum[na.omit(match(human_ensembl_correspond$opossum_ensembl,opossum[,1])),opossum_ind],1,mean,na.rm=T)	
	}
	if(length(opossum_ind)==1){
		temp$opossum[which(!is.na(match(human_ensembl_correspond$opossum_ensembl,opossum[,1])))]<-opossum[na.omit(match(human_ensembl_correspond$opossum_ensembl,opossum[,1])),opossum_ind]	
	}

	if(length(rabbit_ind)>1){
		temp$rabbit[which(!is.na(match(human_ensembl_correspond$rabbit_ensembl,rabbit[,1])))]<-apply(rabbit[na.omit(match(human_ensembl_correspond$rabbit_ensembl,rabbit[,1])),rabbit_ind],1,mean,na.rm=T)
	}
	if(length(rabbit_ind)==1){
		temp$rabbit[which(!is.na(match(human_ensembl_correspond$rabbit_ensembl,rabbit[,1])))]<-rabbit[na.omit(match(human_ensembl_correspond$rabbit_ensembl,rabbit[,1])),rabbit_ind]
	}
	
	row.names(temp)<-human_ensembl_correspond$human_ensembl
	
	return(temp)
}


#

extract_exp_correspond2<-function(sta_org_list){
	
	return(lapply(sta_org_list,function(x){temp<-lapply(1:6,function(i){temp<-data.frame(sapply(1:7,function(j){x[[j]][,i]}));names(temp)<-organ_names;return(temp)});names(temp)<-species_names;return(temp)}))
	
}


###

organ_names<-setdiff(unique(sapply(names(human),function(x){unlist(strsplit(x,"\\."))[1]})),"Names")
stage_names<-paste0("s",c(1:14))
species_names<-c("human","macaque","mouse","rat","rabbit","opossum")


###

stage_organ_list<-list()

for(i in 1:14){
	print(i)
	stage_organ<-list()
	for(j in 1:7){
		stage_organ[[j]]<-extract_exp_correspond(organ_names[j],stage_names[i])
	}
	names(stage_organ)<-organ_names
	stage_organ_list[[i]]<-stage_organ
}
names(stage_organ_list)<-stage_names

stage_species_list<-extract_exp_correspond2(stage_organ_list)
stage_species_list<-lapply(stage_species_list,function(x){lapply(x,function(y){rownames(y)<-human_ensembl_union;return(y)})})

#

brain_asymmetry_gene_overlap_ind<-match(intersect(human_ensembl_union,brain_asymmetry_gene),human_ensembl_union)
stage_species_list_overlapInd<-lapply(stage_species_list,function(x){lapply(x,function(y){y[brain_asymmetry_gene_overlap_ind,c(1,3:7)]})})

#

cal_BSI_pairSpecies<-function(x,y,lowExp=0){
	tempInd<-which(!is.na(x[-1])&!is.na(y[-1])&x[-1]>lowExp&y[-1]>lowExp)
	if(length(tempInd)==0){
		return(c(NA,NA))
	}
	
	if(length(tempInd)==1){
		return(c(log2(x[1]/x[tempInd+1]),log2(y[1]/y[tempInd+1])))
	}
	
	if(length(tempInd)>1){
		return(c(log2(x[1]/mean(x[tempInd+1])),log2(y[1]/mean(y[tempInd+1]))))
	}	
	
}


remove_na_inf<-function(dat){
	indTemp<-which(is.na(dat)|dat==Inf|dat==-Inf)
	if(length(indTemp)==0){
		return(dat)
	}
	if(length(indTemp)>0){
		naInd<-unique(sapply(indTemp,function(x){rowColRec(x,dim(dat)[1],dim(dat)[2])[1]}))
		if(length(naInd)==dim(dat)[1]){
			return(cbind(NA,NA))
		}else{
			return(dat[-naInd,])
		}		
	}
	

}



# lowly expressed genes are filtered with RPKM<1 

BSI_pairSpecies<-list()

for(i in 1:14){
		
	BSI_pairSpecies_focal_tempStage<-list()
	
	for(j in 1:5){
	
		print(paste0("i:",i,";j:",j))

		BSI_focal_pairSpecies<-t(sapply(1:158,function(k){cal_BSI_pairSpecies(as.numeric(stage_species_list_overlapInd[[i]][[1]][k,]),as.numeric(stage_species_list_overlapInd[[i]][[j+1]][k,]),lowExp=1)}))
		rownames(BSI_focal_pairSpecies)<-intersect(human_ensembl_union,brain_asymmetry_gene)
		colnames(BSI_focal_pairSpecies)<-c("BSI_human",paste0("BSI_",species_names[j+1]))
		
		BSI_pairSpecies_focal_tempStage[[j]]<-BSI_focal_pairSpecies
	}
	names(BSI_pairSpecies_focal_tempStage)<-paste0("human_",species_names[2:6])
	
	BSI_pairSpecies[[i]]<-BSI_pairSpecies_focal_tempStage
	
}
names(BSI_pairSpecies)<-stage_names


#

pvalues2<-t(sapply(1:14,function(i){

	sapply(1:5,function(j){
		print(paste0("i:",i,";j:",j))
		focalDF<-BSI_pairSpecies[[i]][[j]]
		focal_temp<-remove_na_inf(focalDF)
		if(!is.na(focal_temp[1,1])){
			return(wilcox.test(focal_temp[,1],focal_temp[,2],paired=T,alternative="greater")$p.value)
		}else{
			return(NA)
		}		

	})
}))

rownames(pvalues2)<-stage_names
colnames(pvalues2)<-paste0(species_names[1],"_",species_names[2:6])

pvalues2_df<-data.frame(pvalue=c(pvalues2),
	species_pair=factor(c(sapply(colnames(pvalues2),rep,dim(pvalues2)[1])),level=paste0(species_names[1],"_",species_names[2:6])),
	N=c(t(sapply(1:14,function(i){sapply(1:5,function(j){
		focalDF<-BSI_pairSpecies[[i]][[j]]
		focal_temp<-remove_na_inf(focalDF)		
		return(ifelse(is.na(focal_temp[1,1]),0,dim(focal_temp)[1]))
		})}))),
	stage=rep(1:14,dim(pvalues2)[2]))
	

	
library(ggplot2)
	
## comparison of BSI between species at stage 13

lapply(stage_species_list_overlapInd[[13]],function(x){apply(x,2,function(y){length(which(!is.na(y)))})})

dataTemp<-data.frame(BSI=c(apply(sapply(1:5,function(i){BSI_pairSpecies[[13]][[i]][,1]}),1,mean,na.rm=T),do.call(rbind,BSI_pairSpecies[[13]][c(1,2,3,5,4)])[,2]),
	species=factor(c(sapply(species_names,rep,158)),level=rev(species_names)),
	gene=rep(rownames(BSI_pairSpecies[[13]][[1]]),6))


ggplot(dataTemp,aes(y=BSI,x=species))+
geom_boxplot(size=0.3,width=0.65,outlier.size=0.25)+
geom_hline(yintercept=0,lty=2,col="#F8766D",lwd=0.3)+
scale_y_continuous(breaks=c(-5,0,5),limits=c(min(dataTemp$BSI,na.rm=T),-min(dataTemp$BSI,na.rm=T)))+
labs(y=expression(paste("BSI at stage 13 (",log[2],"FC)")),x="Species",title="Comparison of BSI between human and other species")+
theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.5))

ggsave("F:/wangj/Desktop/temp/brain asymmetry/figs/BSI_comparison_between_species_stage_13.pdf",width=unit(4,"in"),height=unit(3,"in"))

pvalues2_df[pvalues2_df$stage==13,]
	
## across stages

ggplot(pvalues2_df,aes(x=stage,y=-log10(pvalue),col=species_pair,fill=species_pair,shape=species_pair))+geom_point()+geom_line(lty=1)+
geom_hline(yintercept=-log10(0.05),lty=2,col="#F8766D")+
scale_x_continuous(breaks=c(1:14))+
labs(x="Stage",y=expression(paste("-",log[10],italic(P)," for BSI difference")),title="Comparison of BSI difference across stages")+
theme(legend.position=c(0.4,0.85))

ggsave("F:/wangj/Desktop/temp/brain asymmetry/figs/BSI_comparison_across_stages.pdf",width=unit(4,"in"),height=unit(3,"in"))



# examples of BA genes 

cor(sapply(1:5,function(i){BSI_pairSpecies[[13]][[i]][,1]}),use="pairwise.complete.obs")

BSI_stage_13<-data.frame(apply(sapply(1:5,function(i){BSI_pairSpecies[[13]][[i]][,1]}),1,mean,na.rm=T),sapply(1:5,function(i){BSI_pairSpecies[[13]][[i]][,2]}))
names(BSI_stage_13)<-species_names

cor(BSI_stage_13,use="pairwise.complete.obs")

example_genes_index<-which(apply(BSI_stage_13,1,function(x){ifelse(length(which(!is.na(x)))==6,TRUE,FALSE)&x[1]==max(x,na.rm=T)}))
length(example_genes_index)

pdf("F:/wangj/Desktop/temp/brain asymmetry/figs/BSI_stage_13_example_genes.pdf")

par(mfrow=c(2,2))
for(i in 1:length(example_genes_index)){
	plot(as.numeric(BSI_stage_13[example_genes_index[i],]),type="l",main=paste0(example_genes_index[i],":",rownames(BSI_stage_13)[example_genes_index[i]]))
	abline(h=0,lty=2)
}

dev.off()


geneLabel<-paste(rownames(BSI_stage_13),as.character(gene_tiss_exp[match(rownames(BSI_stage_13),gene_tiss_exp[,1]),2]),sep=" | ")


dataTemp<-data.frame(BSI=c(as.matrix(BSI_stage_13)),
	species=factor(c(sapply(names(BSI_stage_13),rep,dim(BSI_stage_13)[1])),level=rev(species_names)),
	gene=rep(geneLabel,6),
	geneIndex=rep(1:dim(BSI_stage_13)[1],6))
	
dataTemp2<-data.frame(BSI=c(as.matrix(BSI_stage_13)),
	species=c(sapply(6:1,rep,dim(BSI_stage_13)[1])),
	gene=rep(geneLabel,6),
	geneIndex=rep(1:dim(BSI_stage_13)[1],6))
	
select_gene_index<-c(12,14,16,17,50,73,103)

multi_match<-function(x,y){
	return(sort(unlist(lapply(x,function(k){which(y==k)}))))
}

ggplot(dataTemp[multi_match(select_gene_index,dataTemp$geneIndex),],aes(x=species,y=BSI))+
geom_point()+
geom_line(data=dataTemp2[multi_match(select_gene_index,dataTemp$geneIndex),],aes(x=species,y=BSI))+
labs(x="Species",y=expression(paste("BSI of individual genes (",log[2],"FC)")),title="Example genes at stage 13")+
facet_wrap(vars(gene),nrow=7,scales="free_y",strip.position = "top")+
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
	strip.text.x=element_text(size=rel(8/11)),
	strip.background=element_blank())

ggsave("F:/wangj/Desktop/temp/brain asymmetry/figs/BSI_stage_13_example_genes.pdf",width=unit(2.4,"in"),height=unit(7.04,"in"))


### robustness of RPKM>1 : compare with RPKM>0.5 and RPKM>1.5 



# lowly expressed genes are filtered with RPKM<1 

get_pvalues2_df<-function(thresh){

BSI_pairSpecies<-list()

for(i in 1:14){
		
	BSI_pairSpecies_focal_tempStage<-list()
	
	for(j in 1:5){
	
		print(paste0("i:",i,";j:",j))

		BSI_focal_pairSpecies<-t(sapply(1:158,function(k){cal_BSI_pairSpecies(as.numeric(stage_species_list_overlapInd[[i]][[1]][k,]),as.numeric(stage_species_list_overlapInd[[i]][[j+1]][k,]),lowExp=thresh)}))
		rownames(BSI_focal_pairSpecies)<-intersect(human_ensembl_union,brain_asymmetry_gene)
		colnames(BSI_focal_pairSpecies)<-c("BSI_human",paste0("BSI_",species_names[j+1]))
		
		BSI_pairSpecies_focal_tempStage[[j]]<-BSI_focal_pairSpecies
	}
	names(BSI_pairSpecies_focal_tempStage)<-paste0("human_",species_names[2:6])
	
	BSI_pairSpecies[[i]]<-BSI_pairSpecies_focal_tempStage
	
}
names(BSI_pairSpecies)<-stage_names


#

pvalues2<-t(sapply(1:14,function(i){

	sapply(1:5,function(j){
		print(paste0("i:",i,";j:",j))
		focalDF<-BSI_pairSpecies[[i]][[j]]
		focal_temp<-remove_na_inf(focalDF)
		if(!is.na(focal_temp[1,1])){
			return(wilcox.test(focal_temp[,1],focal_temp[,2],paired=T,alternative="greater")$p.value)
		}else{
			return(NA)
		}		

	})
}))

rownames(pvalues2)<-stage_names
colnames(pvalues2)<-paste0(species_names[1],"_",species_names[2:6])

pvalues2_df<-data.frame(pvalue=c(pvalues2),
	species_pair=factor(c(sapply(colnames(pvalues2),rep,dim(pvalues2)[1])),level=paste0(species_names[1],"_",species_names[2:6])),
	N=c(t(sapply(1:14,function(i){sapply(1:5,function(j){
		focalDF<-BSI_pairSpecies[[i]][[j]]
		focal_temp<-remove_na_inf(focalDF)		
		return(ifelse(is.na(focal_temp[1,1]),0,dim(focal_temp)[1]))
		})}))),
	stage=rep(1:14,dim(pvalues2)[2]))
	
return(pvalues2_df)

}


pvalues2_df_RPKM0p5<-get_pvalues2_df(0.5)
pvalues2_df_RPKM1p5<-get_pvalues2_df(1.5)

dataTemp<-data.frame(p_rpkm1=pvalues2_df$pvalue,p_rpkm0p5=pvalues2_df_RPKM0p5$pvalue,p_rpkm1p5=pvalues2_df_RPKM1p5$pvalue)

ggplot(dataTemp,aes(x=p_rpkm1,y=p_rpkm0p5))+
geom_point()+
geom_smooth(method="lm")+
scale_x_log10()+
scale_y_log10()+
labs(x=expression(paste("-",log[10],"(",italic(P),") under RPKM>1")),y=expression(paste("-",log[10],"(",italic(P),") under RPKM>0.5")),title="")

ggsave("F:/wangj/Desktop/temp/brain asymmetry/figs/BSI_robustness_RPKM0p5.pdf",width=unit(3,"in"),height=unit(3.1,"in"))


ggplot(dataTemp,aes(x=p_rpkm1,y=p_rpkm1p5))+
geom_point()+
geom_smooth(method="lm")+
scale_x_log10()+
scale_y_log10()+
labs(x=expression(paste("-",log[10],"(",italic(P),") under RPKM>1")),y=expression(paste("-",log[10],"(",italic(P),") under RPKM>1.5")),title="")


ggsave("F:/wangj/Desktop/temp/brain asymmetry/figs/BSI_robustness_RPKM1p5.pdf",width=unit(3,"in"),height=unit(3.1,"in"))


save.image("F:/wangj/Desktop/temp/brain asymmetry/figs/gene_exp_analysis.RData")


################################################################################

### output 
#BSI_pairSpecies
#pvalues2_df

BSI_pairSpecies_refine<-lapply(1:length(BSI_pairSpecies),function(i){lapply(BSI_pairSpecies[[i]],function(y){
	if(class(y)[1]=="matrix"){
		temp<-rep(paste0(do.call(rbind,strsplit(colnames(y),"BSI_"))[,2],collapse="_"),dim(y)[1])
		temp2<-rep(stage_names[i],dim(y)[1])
		return(cbind(y,species_pair=temp,stage=temp2))
	}
})})


refine_BSI<-do.call(rbind,lapply(BSI_pairSpecies,function(x){do.call(rbind,x)}))

refine_stage<-c(sapply(stage_names,rep,158*5))
refine_species_pair<-rep(c(sapply(names(BSI_pairSpecies[[1]]),rep,158)),14)
BSI_pairSpecies_df<-data.frame(refine_BSI,stage=refine_stage,species_pair=refine_species_pair)
names(BSI_pairSpecies_df)[2]<-"BSI_anoter_species"
BSI_pairSpecies_df$Ensembl=rep(rownames(BSI_pairSpecies[[1]][[1]]),5*14)


# Load the 'openxlsx' package
library(openxlsx)

# Create a new Excel workbook
wb <- createWorkbook()

# Add data frames to separate worksheets
addWorksheet(wb, "BSI_species_stage")
writeData(wb, sheet = 1, x = BSI_pairSpecies_df,keepNA=TRUE)

addWorksheet(wb, "pvalue")
writeData(wb, sheet = 2, x = pvalues2_df,keepNA=TRUE)

# Save the workbook to a file
saveWorkbook(wb, file = "F:/wangj/Desktop/temp/brain asymmetry/MS/Tab.S5-Brain_Specificity_Index_of_Species_Pairs_Across_Stages.xlsx", overwrite = TRUE)











