# download sample metadata from ega into excel file
sdat<-openxlsx::read.xlsx("excel_file.xlsx",sheet="PROJECT_sample_meta") %>% mutate(sample_id=SAMPLE_ALIAS,phenotype=gsub(" ","_",phenotype)) 
ss<-unique(unlist(lapply(strsplit(dir("~/mtpt/fastq",pattern="gz$"),"_S|\\-R[1-2]"),function(X) X[1])))
#check
 sdat$sample_id[!sdat$sample_id %in% ss]
 ss[!ss %in% sdat$sample_id]
 
#arrange
sdat<-arrange(sdat,sample_id)
table(sdat$sample_id == ss)
write.csv(sdat,"PATH/TO/project/source/sample_metadata.csv",row.names=F)