library(dotplot)
library(ggplot2)
library(seqinr)
library(tidyverse)

gigs<-data.table::fread("files/Ssal_Gig1_Gig2_tss.txt",header=T) %>% 
  mutate(gclade=gsub("_.*","",ggroup),gfam=ifelse(grepl('gig1',gclade),"gig1","gig2")) %>% arrange(gfam,gclade,chrom) %>% 
  filter(tentative_nickname != "gig1_ssa04_e2",ggroup != "gig2_ssa04f")


gig1.inf<-read.fasta("files/gig1_cds_sequence.fa",seqtype = "DNA", set.attributes = FALSE)
gig1<-read.fasta("files/gig1_cds_sequence.fa",seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
gig2.inf<-read.fasta("files/gig2_cds_sequence.fa",seqtype = "DNA", set.attributes = FALSE)
gig2<-read.fasta("files/gig2_cds_sequence.fa",seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)

gigs<-data.table::fread("files/Ssal_Gig1_Gig2_tss.txt",header=T)
gigs$tentative_nickname<-factor(gigs$tentative_nickname,levels = gigs$tentative_nickname)

ssal_gig1<-gig1.inf[c(which(grepl('ENSSSAG',names(gig1.inf))))]
ssal_gig2<-gig2.inf[c(which(grepl('ENSSSAG',names(gig2.inf))))]

#Salmon gig1
for (i in c(1:length(names(ssal_gig1)))){
  genseq<-ssal_gig1[[i]]
  genid<-names(ssal_gig1[i])
  gigid<-gigs[gigs$ensembl_id %in% genid,9] %>% pull()
  dotPlot(genseq,genseq, wsize = 10,wstep=1,nmatch=9, xlab=gigid,ylab=gigid)
}

#Salmon gig2
for (z in c(1:length(names(ssal_gig2)))){
  genseq<-ssal_gig2[[z]]
  genid<-names(ssal_gig2[z])
  gigid<-gigs[gigs$ensembl_id %in% genid,9] %>% pull()
  dotPlot(genseq,genseq, wsize = 10,wstep=1,nmatch=9, xlab=gigid,ylab=gigid)
}

