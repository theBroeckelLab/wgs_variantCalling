## Annotation of Variants/Loci using VEP and BiomaRt
## 10.04.2023


## load biomaRt and openxlsx
library(openxlsx)
library(biomaRt)

## read-in variant/loci data
snps=read.xlsx("Z:/Projects/Project Management/Analysis/12.16.2019/Data Analysis/2022-03-28 topMed KO network grant/snps_from_bing_oct2023.xlsx", sheet=1)

## sort variants by genomic chr and coordinate
snps=snps[with(snps, order(CHR_LVM_2DMM, POS_LVM_2DMM)), ]




######################################################################
### Enselbls VEP (https://useast.ensembl.org/Tools/VEP) ##############
######################################################################
## snp data into required txt format
vep2write=paste0(snps$CHR_LVM_2DMM, " ", snps$POS_LVM_2DMM, " ", snps$POS_LVM_2DMM, " ", paste0(snps$REF_LVM_2DMM, "/", snps$ALT_LVM_2DMM), " +")
vep2write

## write table to txt file - then import to VEP to identify any known/deleterious effects of variants
write.table(vep2write, "Z:/Projects/Project Management/Analysis/12.16.2019/Data Analysis/2022-03-28 topMed KO network grant/data_for_vep_test.txt", row.names = F, col.names = F, quote=F)







############################################################
## Ensembl biomart Query using BiomaRt Package #############
#############################################################
## for each variant we want to search up/downstream window for nearby genes
window.size=500000                                                            ##set 500kbp upstream and downstream threshold
paste0("chr", snps$CHR_LVM_2DMM[1], ":", snps$POS_LVM_2DMM[1])                ##location of first variant
paste0("chr", snps$CHR_LVM_2DMM[1], ":", snps$POS_LVM_2DMM[1]-window.size)    ##500kbp upstream of first variant 
paste0("chr", snps$CHR_LVM_2DMM[1], ":", snps$POS_LVM_2DMM[1]+window.size)    ##500kbp downstream of first variant 
## window region for first variant to input into biomaRt
paste0("chr", snps$CHR_LVM_2DMM[1], ":", snps$POS_LVM_2DMM[1]-window.size, "-", snps$POS_LVM_2DMM[1]+window.size)


#pull all window coordinates into one vector
filterlist=paste0(snps$CHR_LVM_2DMM, ":", snps$POS_LVM_2DMM-window.size, ":", snps$POS_LVM_2DMM+window.size)
filterlist

## load human (ie hsapiens_gene_ensembl) biomart database 
hg=useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")


## results for each variant will be saved as a separate dataframe in list 'list2write'
list2write=list()

## loop through each variant location, pull all genes 500kb up/downstream with gene annotations
for (i in 1:length(filterlist)) {
  
  #progress message
  print(paste0("processing variant ", i, " of ", length(filterlist)))
  
  #query biomaRt- search database based on chromosomal_region and return hgnc_symbol chromosome_name start_position etc
  full.results=getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position","gene_biotype","entrezgene_description"),
                filters=c("chromosomal_region", "biotype"),
                values=list(chromosomal_region=filterlist[i], biotype="protein_coding"), mart=hg)
  
  #if there are no genes in the window, save an empty data frame and skip to next iteration
  if(nrow(full.results)==0) {
    list2write[[i]]=data.frame(loci=snps$SNP[i],cis.gene=NA,cis.gene.TSS=NA,cis.gene.TES=NA, intragenic=NA, 
                               dist.to.TSS=NA, cis.gene.type=NA, cis.gene.description=NA, cis.gene.annotation=NA)
    next
  }
  
  #for each gene in window query biomaRt again- search database based on chromosomal_region and return hgnc_symbol and name_1006 (ie GO terms)
  full.results$go.terms=NA
  go.pull=getBM(attributes=c("hgnc_symbol","name_1006"),
                filters=c("chromosomal_region", "biotype"),
                values=list(chromosomal_region=filterlist[i], biotype="protein_coding"), mart=hg)
  
  #collapse Go terms into string, add to initial biomaRt query results
  go.pull=split(go.pull, go.pull$hgnc_symbol)
  for (k in 1:length(go.pull)) {
    full.results$go.terms[match(names(go.pull)[k], full.results$hgnc_symbol)]=sort(paste0(go.pull[[k]]$name_1006, collapse="; "))
  }
  
  #create a dataframe to compile biomaRt results
  df2write=data.frame(loci=rep(NA,nrow(full.results)), cis.gene=NA, cis.gene.TSS=NA, cis.gene.TES=NA, intragenic=NA, 
                      dist.to.TSS=NA, cis.gene.type=NA, cis.gene.description=NA, cis.gene.annotation=NA)
  df2write$loci=snps$SNP[i]                             ##save snp name
  df2write$cis.gene=full.results$hgnc_symbol            ##save gene symbol
  df2write$cis.gene.TSS=full.results$start_position     ##save gene start
  df2write$cis.gene.TES=full.results$end_position       ##save gene end
  #determine if variant is within gene
  for (j in 1:nrow(full.results)) {
    df2write$intragenic[j]="No"
    if(snps$POS_LVM_2DMM[i]>full.results$start_position[j] & snps$POS_LVM_2DMM[i]<full.results$end_position[j]) {df2write$intragenic[j]="Yes"; next}
    #if variant is not within gene, calculate distance from variant to gene start
    df2write$dist.to.TSS[j]=snps$POS_LVM_2DMM[i]-full.results$start_position[j]
  }
  df2write$cis.gene.type=full.results$gene_biotype                         ##save gene type
  df2write$cis.gene.description=full.results$entrezgene_description        ##save gene description
  df2write$cis.gene.annotation=full.results$go.terms                       ##save gene GO terms
  
  #write df to list2write
  list2write[[i]]=df2write
}

## use variant coords to name each dataframe in list2write 
names(list2write)=paste0("chr", paste0(snps$CHR_LVM_2DMM,"-",snps$POS_LVM_2DMM))
View(list2write$`chr2-16147224`)

## write list to excel (each dataframe in list2write is saved as a seperate sheet in the excel workbook)
write.xlsx(list2write, "Z:/Projects/Project Management/Analysis/12.16.2019/Data Analysis/2022-03-28 topMed KO network grant/snps_from_bing_oct2023_annotation.xlsx")


## for any variants that are intra-genic, use VEP (https://useast.ensembl.org/Tools/VEP) to determine if they are intronic, exonic, regulatory, deleterious, etc
list2write$`chr2-167466782`[,1:8]   #variant is within B3GALT1 gene

## entering the below coordinates into VEP indicates an intronic variant in B3GALT1
paste0(c(snps$CHR_LVM_2DMM[9],  snps$POS_LVM_2DMM[9], snps$POS_LVM_2DMM[9], paste0(snps$REF_LVM_2DMM[9], "/", snps$ALT_LVM_2DMM[9]), "+"), collapse =" ")







