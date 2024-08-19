cd PATH/TO/ega_"$snam"/vep
module load R/4.0.4
R
library(vcfR)
library(tidyverse)

rfnam<-dir(pattern="vcf.gz_soma_dp15_qual20_vep_pic_tab.txt$")
vfnam<-dir(pattern="vcf.gz_soma_dp15_qual20_vep_pic.vcf$")

rsnam=gsub("_1pc|.vcf.gz_soma_dp15_qual20_vep_pic_tab.txt","",rfnam)
vsnam=gsub("_1pc|.vcf.gz_soma_dp15_qual20_vep_pic.vcf","",vfnam)
table(rsnam==vsnam)

trimdat<-lapply(1:length(rfnam), function(i){
	tryCatch({
	message(i," ",rsnam[i])
	tab<-read.delim(rfnam[i],skip=63)
	vcf<-read.vcfR(vfnam[i])
	return(
	tibble(tab) %>% filter(grepl("PICK=",Extra)) %>% mutate(
		PolyPhen= ifelse(grepl("PolyPhen",Extra),gsub("PolyPhen=","",grep("PolyPhen",unlist(strsplit(Extra,";")),value=T)),NA),
		SIFT=ifelse(grepl("SIFT",Extra),gsub("SIFT=","",grep("SIFT",unlist(strsplit(Extra,";")),value=T)),NA),
		IMPACT=ifelse(grepl("IMPACT",Extra),gsub("IMPACT=","",grep("IMPACT",unlist(strsplit(Extra,";")),value=T)),NA),
		MAX_AF=ifelse(grepl("MAX_AF",Extra),gsub("MAX_AF=","",grep("MAX_AF=",unlist(strsplit(Extra,";")),value=T)),NA)
		)  %>% mutate(
		polyn= tryCatch(unlist(lapply(strsplit(PolyPhen,"\\(|\\)"), function(X) X[2])), error=function(e) NA),
		polyc= tryCatch(unlist(lapply(strsplit(PolyPhen,"\\(|\\)"), function(X) X[1])), error=function(e) NA),
		siftn= tryCatch(unlist(lapply(strsplit(SIFT,"\\(|\\)"), function(X) X[2])), error=function(e) NA),
		siftc= tryCatch(unlist(lapply(strsplit(SIFT,"\\(|\\)"), function(X) X[1])), error=function(e) NA)
		) %>% left_join(data.frame(vcf@fix[,1:6]) %>% mutate(X.Uploaded_variation=paste0(CHROM,"_",POS,"_",REF,"/",ALT)) %>% cbind(.,extract_info_tidy(vcf)) %>% tibble)
	)
	}
	,error=function(e) NA)
}
)
names(trimdat)<-rsnam


library(tidyverse)
library(data.table)
library(ggpubr)
library(patchwork)

sift_1<-lapply(trimdat, function(X){
group_by(X,siftc,IMPACT,Consequence) %>% summarise(count=n()) %>% tibble
})
for (i in 1:length(sift_1)) sift_1[[i]]$list=names(sift_1)[i]

poly_1<-lapply(trimdat, function(X){
group_by(X,polyc,,IMPACT,Consequence) %>% summarise(count=n()) %>% tibble
})
for (i in 1:length(poly_1)) poly_1[[i]]$list=names(poly_1)[i]

# bin, sum and calculate; DEV (id,sample,comparator) COZ NAME HAS EXTRA UNDERSCORE
sift<-group_by(rbindlist(sift_1),list)  %>% summarize(
	total=sum(count),
	nonSynon=total-sum(count[is.na(siftc) & !grepl("HIGH|MODERATE",IMPACT)]),
	benign=sum(count[grepl("tolerated",siftc)]),
	deleterious=sum(count[grepl("deleterious",siftc) | (is.na(siftc) & grepl("HIGH",IMPACT)) ]),
	undefined=sum(count[(is.na(siftc) & grepl("MODERATE",IMPACT)) ])
	) %>% mutate(
		perc_total_benign=100*benign/total,
		perc_total_DPM=100*deleterious/total,
		perc_NonSyn_benign=100*benign/nonSynon,
		perc_NonSyn_DPM=100*deleterious/nonSynon,
		benignVSDPM=benign/deleterious,
		id=unlist(lapply(strsplit(list,"_"),function(X) paste0(X[1]))),
		sample=unlist(lapply(strsplit(list,"_"),function(X) X[2])),
		comparator=gsub("^v","",unlist(lapply(strsplit(list,"_"),function(X) X[4])))
	) #changeID

poly<-group_by(rbindlist(poly_1),list) %>% summarize(
		total=sum(count),
		nonSynon=total-sum(count[is.na(polyc) & !grepl("HIGH|MODERATE",IMPACT)]),
		benign=sum(count[grepl("benign",polyc) | grepl("unknown",polyc)]),
		deleterious=sum(count[grepl("damaging",polyc) | (is.na(polyc) & grepl("HIGH",IMPACT))]),
		undefined=sum(count[(is.na(polyc) & grepl("MODERATE",IMPACT))])
		) %>% mutate(
			perc_total_benign=100*benign/total,
			perc_total_DPM=100*deleterious/total,
			perc_NonSyn_benign=100*benign/nonSynon,
			perc_NonSyn_DPM=100*deleterious/nonSynon,
			benignVSDPM=benign/deleterious,
			id=unlist(lapply(strsplit(list,"_"),function(X) paste0(X[1]))),
			sample=unlist(lapply(strsplit(list,"_"),function(X) X[2])),
			comparator=gsub("^v","",unlist(lapply(strsplit(list,"_"),function(X) X[4])))
		)#changeID



# metadata different coz metadata is trash
 samd<-  read.csv("PATH/TO/Intersection_dataframe.csv")[,-1]  %>% tibble() 
 md<- samd %>% group_by(id) %>% summarise(gender=unique(donor_gender),condition=paste0(unique(cond), collapse=",")) %>% mutate(id=as.character(id))

subset(md,!md$id %in% poly$id) 
subset(poly,!poly$id %in% md$id)


forSumm=list(
	sift_bigtab=rbindlist(sift_1) %>% pivot_wider(values_from=count,names_from=list) %>% rename_with("make.unique"),
	sift_summary=left_join(sift,md)%>% rename_with("make.unique"),
	poly_bigtab=rbindlist(poly_1) %>% pivot_wider(values_from=count,names_from=list)%>% rename_with("make.unique"),
	poly_summary=left_join(poly,md)%>% rename_with("make.unique"),
	criteria=data.frame(procedure=eval('total=all var that arent intergenic
nonSynon=total - (no sift/poly score but not High or Moderate IMPACT)
benign= "tolerated" in sift, "benign" or "unknown" in polyphen
deleterious= "deleterious" in sift, "damaging" in polyphen ; PLUS (no sift/poly score but High IMPACT)
undefined= "no sift/poly score but moderate IMPACT'))
)

openxlsx::write.xlsx(forSumm,file="PATH/TO/ega_"$snam"_1pc_soma_DP15_QUAL20_VEP_PICK_ANNO_COMB_COUNT_CALCULATE.xlsx")
saveRDS(forThom,"PATH/TO/ega_"$snam"_1pc_soma_DP15_QUAL20_VEP_PICK_ANNO_COMB_COUNT_CALCULATE.RDS")


# adding GTEX expression
# equal sized bins
 gt<-data.table::fread("data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz") %>% dplyr::select(Name,Description,Liver)
summary(log2(gt$Liver[gt$Liver>0]))

# cut
gg<-gt %>% mutate(gtex_bin=cut(log2(Liver+1e-5),breaks=c(-Inf,-10.69,-2.69,0.1,2.75,Inf),labels=c("0_undetected","1st_quartile","2nd_quartile","3rd_quartile","4th_quartile"))) %>% mutate(ensid=gsub("\\..+$","",Name)) 

#annotate
library(biomaRt)
mouse<- useMart(biomart='ensembl', dataset = "hsapiens_gene_ensembl")
hasEntrez2<-getBM(attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol"),filters="ensembl_gene_id",values=gg$ensid ,mart=mouse) 
ggg<-left_join(gg,hasEntrez2 %>% group_by(ensembl_gene_id) %>% summarize(symbols=paste0(hgnc_symbol,collapse=";"),ENTREZID=entrezgene_id[1]),by=c("ensid"="ensembl_gene_id")) %>% group_by(ensid) %>% arrange(gtex_bin) %>% summarize(ids=paste0(Description,collapse=","),gtex_bin=gtex_bin[length(gtex_bin)],ENTREZID=ENTREZID[1])


# combining for later
listofComb<-lapply(fnam, function(exp){
    print(grep(exp,fnam))
nam<-gsub("_VEP_PICK_ANNO.RDS|_soma_DP15_QUAL20|r_objects/[0-9]{6}_","",exp)
 dbgap<-readRDS(exp)

 if (! "IMPACT" %in% colnames(dbgap[[1]])) dbgap<-lapply(dbgap[!is.na(dbgap)], function(X) mutate(X, IMPACT=ifelse(grepl("IMPACT",Extra),gsub("IMPACT=","",grep("IMPACT",unlist(strsplit(Extra,";")),value=T)),NA)))

comb1<-lapply(dbgap[!is.na(dbgap)], function(X){
tmp<-X %>% left_join(ggg,by=c('Gene'='ensid')) %>% mutate(
	sift_nonSynon=!(is.na(siftc) & !grepl("HIGH|MODERATE",IMPACT)),
    poly_nonSynon=!(is.na(polyc) & !grepl("HIGH|MODERATE",IMPACT)),
	sift_mut_class=ifelse(grepl("tolerated",siftc),"tolerated",ifelse(grepl("deleterious",siftc) | (is.na(siftc) & grepl("HIGH",IMPACT)),"deleterious",ifelse((is.na(siftc) & grepl("MODERATE",IMPACT)) ,"undetermined","not_assigned"))),
    poly_mut_class=ifelse(grepl("benign",polyc) | grepl("unknown",polyc),"benign",ifelse(grepl("damaging",polyc) | (is.na(polyc) & grepl("HIGH",IMPACT)),"deleterious",ifelse((is.na(polyc) & grepl("MODERATE",IMPACT)) ,"undefined","not_assigned"))))
})
for (i in 1:length(comb1)) comb1[[i]]$list=names(comb1)[i]
comb<-rbindlist(comb1)


sift_go=tryCatch(compareCluster(ENTREZID~sift_mut_class,data=subset(comb,sift_nonSynon),fun="enrichGO",OrgDb="org.Hs.eg.db",universe=as.character((subset(comb,!sift_nonSynon))$ENTREZID),ont="BP") %>% setReadable(OrgDb="org.Hs.eg.db"),error=function(e) NA)

sift_go_bin=tryCatch(compareCluster(ENTREZID~sift_mut_class+gtex_bin,data=subset(comb,sift_nonSynon),fun="enrichGO",OrgDb="org.Hs.eg.db",universe=as.character((subset(comb,!sift_nonSynon))$ENTREZID),ont="BP") %>% setReadable(OrgDb="org.Hs.eg.db"),error=function(e) NA)

poly_go=tryCatch(compareCluster(ENTREZID~poly_mut_class,data=subset(comb,poly_nonSynon),fun="enrichGO",OrgDb="org.Hs.eg.db",universe=as.character((subset(comb,!poly_nonSynon))$ENTREZID),ont="BP") %>% setReadable(OrgDb="org.Hs.eg.db"),error=function(e) NA)

poly_go_bin=tryCatch(compareCluster(ENTREZID~poly_mut_class+gtex_bin,data=subset(comb,poly_nonSynon),fun="enrichGO",OrgDb="org.Hs.eg.db",universe=as.character((subset(comb,!poly_nonSynon))$ENTREZID),ont="BP") %>% setReadable(OrgDb="org.Hs.eg.db"),error=function(e) NA)

list(
	exp_name=nam,
    data=comb,
	sift_go=sift_go,
	sift_go_bin=sift_go_bin,
	poly_go=poly_go,
    poly_go_bin=poly_go_bin
)
})

## adding alphafold
af<-data.table::fread("data/AlphaMissense_hg38.tsv.gz")
colnames(af)=paste0("alphaFold_",colnames(af))

# to sort to single value per variant (like pick i guess)

alphaFold_data=lapply(all_dat, function(X) {
unique_anno=inner_join(dplyr::select(X$data %>% mutate(POS=as.integer(POS)),CHROM,POS,REF,ALT),af,by=c("CHROM"="alphaFold_#CHROM","POS"="alphaFold_POS","REF"="alphaFold_REF","ALT"="alphaFold_ALT"),relationship = "many-to-many") %>% group_by(CHROM,POS,REF,ALT) %>% arrange(-alphaFold_am_pathogenicity) %>% slice_head(n=1) 

tes<-list(
    data=left_join(X$data %>% mutate(POS=as.integer(POS)),unique_anno,relationship = "many-to-one") %>% mutate(AF_nonSynon=!(is.na(alphaFold_am_class) & !grepl("HIGH|MODERATE",IMPACT)), AF_mut_class=ifelse(grepl("likely_benign|ambiguous",alphaFold_am_class),"benign",ifelse(grepl("likely_pathogenic",alphaFold_am_class) | (is.na(alphaFold_am_class) & grepl("HIGH",IMPACT)),"deleterious",ifelse((is.na(alphaFold_am_class) & grepl("MODERATE",IMPACT)) ,"undetermined","not_assigned")))),
    exp_name=X$exp_name,
    manuscript_label=X$manuscript_label
)
})


sdat<-data.table::rbindlist(lapply(names(alphaFold_data), function(X) alphaFold_data[[X]]$data %>% group_by(list) %>% summarize(nvar=n(),med_cov=median(DP,na.rm=T),study=X)))

# summarise with trimdata

proc_vars=function(annodf){
    list(sift=group_by(annodf,list,sift_mut_class)  %>% summarize(
	total=n()) %>% pivot_wider(names_from="sift_mut_class",values_from="total"),
    poly=group_by(annodf,list,poly_mut_class)  %>% summarize(
	total=n()) %>% pivot_wider(names_from="poly_mut_class",values_from="total"),
    alphaFold=group_by(annodf,list,AF_mut_class)  %>% summarize(total=n()) %>% pivot_wider(names_from="AF_mut_class",values_from="total")
    )
    }

lapply(names(alphaFold_data), function(X) proc_vars(alphaFold_data[[X]]$data)


### enrichment

# get postitional enriechment
library(msigdbr)
m_df = msigdbr(species = "Homo sapiens",category="C1")
term2gene_pos=dplyr::select(m_df,gs_name,entrez_gene)  %>% rename(c("ENTREZID"="entrez_gene"))

term2gene_bin=dplyr::select(ggg,gtex_bin,ENTREZID) 

library(clusterProfiler)

annos<-lapply(alphaFold_data,function(X){
    print(X$exp_name)
    list(
        exp_name=X$exp_name,
        manuscript_label=X$manuscript_label,
        sift_go=tryCatch(compareCluster(ENTREZID~sift_mut_class,data=subset(X$data,sift_nonSynon),fun="enrichGO",OrgDb="org.Hs.eg.db",universe=as.character((subset(X$data,!sift_nonSynon))$ENTREZID),ont="BP") %>% setReadable(OrgDb="org.Hs.eg.db"),error=function(e) NA),
        sift_bin=tryCatch(compareCluster(ENTREZID~sift_mut_class,data=subset(X$data,sift_nonSynon),fun="enricher",TERM2GENE=term2gene_bin,universe=as.character((subset(X$data,!sift_nonSynon))$ENTREZID)) %>% setReadable(OrgDb="org.Hs.eg.db",keyType="ENTREZID"),error=function(e) NA),
        sift_pos=tryCatch(compareCluster(ENTREZID~sift_mut_class,data=subset(X$data,sift_nonSynon),fun="enricher",TERM2GENE=term2gene_pos,universe=as.character((subset(X$data,!sift_nonSynon))$ENTREZID)) %>% setReadable(OrgDb="org.Hs.eg.db",keyType="ENTREZID"),error=function(e) NA),
        poly_go=tryCatch(compareCluster(ENTREZID~poly_mut_class,data=subset(X$data,poly_nonSynon),fun="enrichGO",OrgDb="org.Hs.eg.db",universe=as.character((subset(X$data,!poly_nonSynon))$ENTREZID),ont="BP") %>% setReadable(OrgDb="org.Hs.eg.db"),error=function(e) NA),
        poly_bin=tryCatch(compareCluster(ENTREZID~poly_mut_class,data=subset(X$data,poly_nonSynon),fun="enricher",TERM2GENE=term2gene_bin,universe=as.character((subset(X$data,!poly_nonSynon))$ENTREZID)) %>% setReadable(OrgDb="org.Hs.eg.db",keyType="ENTREZID"),error=function(e) NA),
        poly_pos=tryCatch(compareCluster(ENTREZID~poly_mut_class,data=subset(X$data,poly_nonSynon),fun="enricher",TERM2GENE=term2gene_pos,universe=as.character((subset(X$data,!poly_nonSynon))$ENTREZID)) %>% setReadable(OrgDb="org.Hs.eg.db",keyType="ENTREZID"),error=function(e) NA),
         alphaFold_go=tryCatch(compareCluster(ENTREZID~AF_mut_class,data=subset(X$data,AF_nonSynon),fun="enrichGO",OrgDb="org.Hs.eg.db",universe=as.character((subset(X$data,!AF_nonSynon))$ENTREZID),ont="BP") %>% setReadable(OrgDb="org.Hs.eg.db"),error=function(e) NA),
        alphaFold_bin=tryCatch(compareCluster(ENTREZID~AF_mut_class,data=subset(X$data,AF_nonSynon),fun="enricher",TERM2GENE=term2gene_bin,universe=as.character((subset(X$data,!AF_nonSynon))$ENTREZID)) %>% setReadable(OrgDb="org.Hs.eg.db",keyType="ENTREZID"),error=function(e) NA),
        alphaFold_pos=tryCatch(compareCluster(ENTREZID~AF_mut_class,data=subset(X$data,AF_nonSynon),fun="enricher",TERM2GENE=term2gene_pos,universe=as.character((subset(X$data,!AF_nonSynon))$ENTREZID)) %>% setReadable(OrgDb="org.Hs.eg.db",keyType="ENTREZID"),error=function(e) NA)
    )
})


# diver mutations from https://www.nature.com/articles/s41568-020-0290-x
# from https://www.intogen.org/search?cancer=HCC
library(data.table)
drivers<-fread("data/240527_IntOGen-DriverGenes_HCC.tsv")

library(survminer)
library(survival)
library(tidyverse)
library(patchwork)
library(pROC)
library(gridExtra)
library(data.table)


all_dat<-readRDS("PATH/TO_allVars_8_studies_plus_AlphaFold.RDS")

for (i in 1:length(all_dat)) {
    ind<-all_dat[[i]]$data$id %in% drivers$Symbol
    all_dat[[i]]$data$isLiverDriver=ifelse(ind & all_dat[[i]]$data$IMPACT=="HIGH","driver_damaging",ifelse(ind,"driver_notDamaging","not_driver"))
}

allvarsum<-rbindlist(lapply(all_dat, function(X) select(X$data,IMPACT,polyc,siftc,alphaFold_am_class,poly_mut_class,sift_mut_class,AF_mut_class,isLiverDriver)))


table(paste(allvarsum$IMPACT,allvarsum$polyc),allvarsum$poly_mut_class)
table(paste(allvarsum$IMPACT,allvarsum$siftc),allvarsum$sift_mut_class)
table(paste(allvarsum$IMPACT,allvarsum$alphaFold_am_class),allvarsum$AF_mut_class)
table(allvarsum$poly_mut_class,allvarsum$isLiverDriver)

#grab followup data

md<-lapply(sapply(names(all_dat), function(x) grep(x,dir("2022_liver_variants_project_info/",pattern="ANNO",full.names=T),value=T)[1]),function(y) {
    print(y) 
    openxlsx::read.xlsx(y,sheet="poly_summary") %>% dplyr::select(-c(total,nonSynon,benign,deleterious,undefined,perc_total_benign,perc_total_DPM,perc_NonSyn_benign,perc_NonSyn_DPM,benignVSDPM)) %>% distinct(list,id,sample,comparator) 
 })
names(md)<-names(all_dat)

fudat<-list(
 #followup data curated by coauthor
)
# had some issues (HC23 e.g. with multiple values, picked worst)

lapply(fudat,function(X){
    print (length(X$id)) 
    print(length(unique(X$id)))
})

# samples with folowup data
proc_vars=function(nam){
    list(sift=group_by(all_dat[[nam]]$data,list,sift_mut_class)  %>% summarize(total=n()) %>% pivot_wider(names_from="sift_mut_class",values_from="total") %>% mutate(nonSyn=sum(del
eterious,tolerated,undetermined),perc_DPM_nonSyn=100*deleterious/nonSyn,perc_benign_nonSyn=100*tolerated/nonSyn) %>% left_join(md[[nam]],by="list") %>% left_join(fudat[[nam]],by="i
d"),
    poly=group_by(all_dat[[nam]]$data,list,poly_mut_class)  %>% summarize(
total=n()) %>% pivot_wider(names_from="poly_mut_class",values_from="total")  %>% mutate(nonSyn=sum(deleterious,benign,undefined),perc_DPM_nonSyn=100*deleterious/nonSyn,perc_benign_
nonSyn=100*benign/nonSyn) %>% left_join(md[[nam]],by="list") %>% left_join(fudat[[nam]],by="id"),
    alphaFold=group_by(all_dat[[nam]]$data,list,AF_mut_class)  %>% summarize(total=n()) %>% pivot_wider(names_from="AF_mut_class",values_from="total")  %>% mutate(nonSyn=sum(delete
rious,benign,undetermined),perc_DPM_nonSyn=100*deleterious/nonSyn,perc_benign_nonSyn=100*benign/nonSyn) %>% left_join(md[[nam]],by="list") %>% left_join(fudat[[nam]],by="id"),
    driver=group_by(all_dat[[nam]]$data ,list,isLiverDriver)  %>% summarize(total=n()) %>% pivot_wider(names_from="isLiverDriver",values_from="total") %>% mutate(totalvars=sum(driv
er_notDamaging,not_driver,driver_damaging,na.rm=T))  %>% left_join(md[[nam]],by="list") %>% left_join(fudat[[nam]],by="id")
    )
    }
dat<-lapply(names(fudat), function(X) proc_vars(X))
names(dat)<-names(fudat)


# The rest
fproc_vars=function(nam){
    list(sift=group_by(all_dat[[nam]]$data,list,sift_mut_class)  %>% summarize(total=n()) %>% pivot_wider(names_from="sift_mut_class",values_from="total") %>% mutate(nonSyn=sum(deleterious,tolerated,undetermined),perc_DPM_nonSyn=100*deleterious/nonSyn,perc_benign_nonSyn=100*tolerated/nonSyn) %>% left_join(md[[nam]],by="list")%>% left_join(md[[nam]],by="list"),
    poly=group_by(all_dat[[nam]]$data,list,poly_mut_class)  %>% summarize(
	total=n()) %>% pivot_wider(names_from="poly_mut_class",values_from="total")  %>% mutate(nonSyn=sum(deleterious,benign,undefined),perc_DPM_nonSyn=100*deleterious/nonSyn,perc_benign_nonSyn=100*benign/nonSyn)  ,
    alphaFold=group_by(all_dat[[nam]]$data,list,AF_mut_class)  %>% summarize(total=n()) %>% pivot_wider(names_from="AF_mut_class",values_from="total")  %>% mutate(nonSyn=sum(deleterious,benign,undetermined),perc_DPM_nonSyn=100*deleterious/nonSyn,perc_benign_nonSyn=100*benign/nonSyn) %>% left_join(md[[nam]],by="list") ,
    driver=group_by(all_dat[[nam]]$data ,list,isLiverDriver)  %>% summarize(total=n()) %>% pivot_wider(names_from="isLiverDriver",values_from="total")  %>% mutate(totalvars=sum(driver_notDamaging,not_driver,driver_damaging,na.rm=T)) %>%left_join(md[[nam]],by="list")
    )
    }

fdat<-lapply(names(all_dat)[!names(all_dat) %in% names(fudat)], function(X) fproc_vars(X))
names(fdat)<-names(all_dat)[!names(all_dat) %in% names(fudat)]
#

ccdat<-c(dat,fdat)
for (i in 1:names(ccdat)){
    openxlsx::write.xlsx(ccdat[[i]]],file=paste0("PATH/",.date(),"_",all_dat[[i]]$manuscript_label,"_scores_PLUS_driversAndNtotal_survWhereAvail.xlsx"),
asTable=T)
}
for (i in 1:names(ccdat)){
    openxlsx::write.xlsx(ccdat[[i]],file=paste0("PATH/",.date(),"_",all_dat[[i]]$manuscript_label,"_scores_PLUS_driversAndNtotal_survWhereAvail.xlsx"),a
sTable=T)
}
i in 1:names(ccdat)
names(ccdat)
for (i in names(ccdat)){
    openxlsx::write.xlsx(ccdat[[i]],file=paste0("PATH/",.date(),"_",all_dat[[i]]$manuscript_label,"_scores_PLUS_driversAndNtotal_survWhereAvail.xlsx"),a
sTable=T)
}

ccdat<-c(dat,fdat)

for (i in names(ccdat)){
    openxlsx::write.xlsx(ccdat[[i]],file=paste0("PATH/",.date(),"_",all_dat[[i]]$manuscript_label,"_scores_PLUS_driversAndNtotal_survWhereAvail.xlsx"),asTable=T)
}
