library(tidyverse)
library(DESeq2)
library(ggrepel)
library(hrbrthemes)
library(patchwork)

#gig genes 
wgd<-data.table::fread("files/gig_wd_LORe_AORe_td.tsv",header=T) %>% mutate(gig = gsub("_(?=[a-g])", "", gig,perl=TRUE),gfam=gsub('_.*','',gig))

gig1<-wgd %>% transmute(gid=gene,gig) %>% filter(grepl('gig1',gig))
gig2<-wgd %>% transmute(gid=gene,gig) %>% filter(grepl('gig2',gig))

gigs<-rbind(gig1,gig2)

#load RNAseq expression from BMC study 
counts.tab<-data.table::fread("files/sav3.salmon.merged.gene_counts.tsv",header = T)  %>% select(-gene_name) %>% as.data.frame()

#metadata
sample=c("SRR10482313","SRR10482314","SRR10482315","SRR10482316","SRR10482317","SRR10482318","SRR10482319","SRR10482320","SRR10482321","SRR10482322","SRR10482323","SRR10482324")
pheno=c("survived","survived","survived","survived","died_11_dpi","died_11_dpi","died_11_dpi","died_11_dpi","survived","survived","died_11_dpi","died_11_dpi") 

meta<-cbind(sample,pheno) %>% as.data.frame() %>% mutate(rows=sample) %>% tibble::column_to_rownames(var="rows")
meta$pheno<-as.factor(pheno)

counts.fixed<-as.data.frame(counts.tab [,-1])
rownames(counts.fixed)<-counts.tab %>% pull(gene_id)
counts.fixed<-round(counts.fixed)

xpr<- counts.fixed[colnames(counts.fixed) %in% meta$sample]
xpr<-counts.fixed[,meta$sample]

#sanity check
all(colnames(xpr) %in% rownames(meta))
all(colnames(xpr) == rownames(meta))

#create deseq2 object
dds<-DESeqDataSetFromMatrix(countData = xpr,
                            colData = meta,
                            design = ~pheno)

#remove genes with 0 expression throughout
keep<-rowSums(counts(dds)) > 0
table(keep)
dds<-dds[keep,]

#estimate the Deseq model
dds<-DESeq(dds)

#estimate the stats for all genes
res<-results(dds,alpha=0.05,contrast = c("pheno","died_11_dpi","survived")) %>% as.data.frame() %>% tibble::rownames_to_column(var="gid")
gres<-res %>% inner_join(gigs,by="gid") %>% mutate(sample="SAV3: Heart tissue")

data.table::fwrite(hktgres,"files/deg.sav3.gene_tpm.tsv",col.names = T,row.names = F,sep="\t",quote = F)

###################################
###################################
###  aquafaang RNAseq data
gcounts<-data.table::fread("files/polyIC.salmon.merged.gene_counts.tsv",header = T)  %>% 
  select(-gene_name) %>% as.data.frame()

gmeta<-data.table::fread("files/polyIC.pbs.metadata") %>% 
  separate(sample_alias,c("ss","uni","vivo","pheno","sample"),"_") %>% transmute(sample=run_accession,pheno,vivo) %>% filter(pheno != "vibrio") %>% 
  mutate(rows=sample,pheno=as_factor(pheno)) %>% tibble::column_to_rownames(var="rows") %>% 
  mutate(experiment=ifelse(vivo=="invitro","HK_primary_cells","HK_tissue"))

gcounts.fixed<-as.data.frame(gcounts [,-1]) 
gcounts.fixed<-round(gcounts.fixed)
rownames(gcounts.fixed)<-gcounts %>% pull(gene_id)


#######################
#Head kidney tissue DEG
hktmeta<-gmeta %>% filter(experiment=="HK_tissue")

hktgxpr<- gcounts.fixed[colnames(gcounts.fixed) %in% hktmeta$sample]
hktgxpr<-gcounts.fixed[,hktmeta$sample]

#sanity check
all(colnames(hktgxpr) %in% rownames(hktmeta))
all(colnames(hktgxpr) == rownames(hktmeta))

#create deseq2 object
hktdds<-DESeqDataSetFromMatrix(countData = hktgxpr,
                            colData = hktmeta,
                            design = ~ pheno)

#remove genes with 0 expression throughout
gkeep<-rowSums(counts(hktdds)) > 0
table(gkeep)
hktdds<-hktdds[gkeep,]

#estimate the Deseq model
hktdds<-DESeq(hktdds)

#what are the possible contrasts?
resultsNames(hktdds)

#estimate the stats for all genes between HK cell line PBS and poly:IC
hktres<-results(hktdds,alpha=0.05,contrast = c("pheno","polyIC","pbs")) %>% as.data.frame() %>% tibble::rownames_to_column(var="gid")
hktgres<-hktres %>% inner_join(gigs,by="gid") %>% mutate(sample="polyIC: HK tissue")

#write the DEG
data.table::fwrite(hktgres,"files/deg.hktissue.gene_tpm.tsv",col.names = T,row.names = F,sep="\t",quote = F)


#######################
#Head kidney cell line
hkcmeta<-gmeta %>% filter(experiment=="HK_primary_cells")

hkcgxpr<- gcounts.fixed[colnames(gcounts.fixed) %in% hkcmeta$sample]
hkcgxpr<-gcounts.fixed[,hkcmeta$sample]

#sanity check
all(colnames(hkcgxpr) %in% rownames(hkcmeta))
all(colnames(hkcgxpr) == rownames(hkcmeta))

#create deseq2 object
hkcdds<-DESeqDataSetFromMatrix(countData = hkcgxpr,
                               colData = hkcmeta,
                               design = ~ pheno)

#remove genes with 0 expression throughout
gkeep<-rowSums(counts(hkcdds)) > 0
table(gkeep)
hkcdds<-hkcdds[gkeep,]

#estimate the Deseq model
hkcdds<-DESeq(hkcdds)

#what are the possible contrasts?
resultsNames(hkcdds)

#estimate the stats for all genes between HK cell line PBS and poly:IC
hkcres<-results(hkcdds,alpha=0.05,contrast = c("pheno","polyIC","pbs")) %>% as.data.frame() %>% tibble::rownames_to_column(var="gid")
hkcgres<-hkcres %>% inner_join(gigs,by="gid") %>% mutate(sample="polyIC:: HK primary cell line")

data.table::fwrite(hkcgres,"files/deg.hkcells.gene_tpm.tsv",col.names = T,row.names = F,sep="\t",quote = F)
