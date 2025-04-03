library(tidyverse)
library(ggrepel)
library(gggenomes)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(UpSetR)
library(aplot)

#load all gigs 
gigs<-data.table::fread("files/Ssal_Gig1_Gig2_tss.txt",header=T) %>% 
  mutate(gclade=gsub("_.*","",ggroup),gfam=ifelse(grepl('gig1',gclade),"gig1","gig2")) %>% arrange(gfam,gclade,chrom) %>% 
  filter(tentative_nickname != "gig1_ssa04_e2",ggroup != "gig2_ssa04f")

gigs$gclade<-factor(gigs$gclade,levels=unique(gigs$gclade))
gigs$gfam<-factor(gigs$gfam,levels=unique(gigs$gfam))

#load the pfam results 
pfc<-data.table::fread("files/gig_interpro_domains.tsv") %>% 
  mutate(gene=gsub("_1","",V1),
         V9=ifelse(V9=="-",1,V9)) %>% select(-V1) %>%
  mutate(V9=as.numeric(V9)) 


#grab Atlantic salmon genes only and incorporate the phylogenetic gene clusters
clusters<-pfc %>% filter(grepl('ENSSSAG',gene)) %>% inner_join(.,gigs,by=c('gene'='ensembl_id')) 

#gig1 genes
g1cpf<-clusters %>% filter(gfam=="gig1") %>% 
  transmute(gene=ggroup,cluster=gclade,chrom,
            db=V4,db.id=V5,db.desc=V6,pval=as.numeric(V9),aa.start=V7,aa.stop=V8,aa=V3) %>% 
  arrange(cluster,chrom,gene)

g2cpf<-clusters %>% filter(gfam=="gig2") %>% 
  transmute(gene=ggroup,cluster=gclade,chrom,
            db=V4,db.id=V5,db.desc=V6,pval=as.numeric(V9),aa.start=V7,aa.stop=V8,aa=V3) %>% 
  arrange(cluster,chrom,gene)

#how many in gig1, how many in gig2?
g1cpf %>% distinct(db.desc,db.id,db)
g2cpf %>% distinct(db.desc,db.id,db)


#assign functional clusters 
filtgig1<-g1cpf %>%
  mutate(db.desc=gsub("Region of a membrane-bound protein predicted to be outside .*","Membrane-bound protein region (outside)",db.desc),
         db.desc=gsub("Region of a membrane-bound protein predicted to be embedded .*","Membrane-bound protein region (embedded)",db.desc),
         db.desc=gsub(".*region of a signal peptide.","Signal peptide region",db.desc),db.desc=gsub("Signa.*","Signal peptide region",db.desc),
         db.desc=gsub("Hydrophobic.*","Signal peptide region",db.desc))

filtgig2<-g2cpf %>%
  mutate(db.desc=gsub("Region of a membrane-bound protein predicted to be outside .*","Membrane-bound protein region (outside)",db.desc),
         db.desc=gsub("Region of a membrane-bound protein predicted to be embedded .*","Membrane-bound protein region (embedded)",db.desc),
         db.desc=gsub(".*region of a signal peptide.","Signal peptide region",db.desc),
         db.desc=gsub("SignalP-noTM","Signal peptide region",db.desc),
         db.desc=gsub("Grass carp reovirus.*","Gig2-(like) protein",db.desc),
         db.desc=gsub("GIG2-LIKE PROTEIN DRED-RELATED","Gig2-(like) protein",db.desc),
         db.desc=gsub("PARP catalytic domain profile.","ADP-ribosylation",db.desc),
         db.desc=gsub("Diphtheria Toxin, domain 1","ADP-ribosylation",db.desc),
         db.desc=gsub(".*polymerase catalytic domain","ADP-ribosylation",db.desc))

## gggenome plot preparation
# Function to merge intervals with improved handling for edge cases
merge_intervals <- function(data) {
  if (nrow(data) == 1) {
    return(data)
  }
  
  data <- data %>% arrange(aa.start)
  merged <- data.frame(aa.start = integer(0), aa.stop = integer(0))
  
  current_start <- data$aa.start[1]
  current_stop <- data$aa.stop[1]
  
  for (i in 2:nrow(data)) {
    if (!is.na(data$aa.start[i]) && !is.na(current_stop) && data$aa.start[i] <= current_stop) {
      # If overlapping or adjacent, merge intervals
      current_stop <- max(current_stop, data$aa.stop[i], na.rm = TRUE)
    } else {
      # If no overlap, save the current interval and start a new one
      merged <- rbind(merged, data.frame(aa.start = current_start, aa.stop = current_stop))
      current_start <- data$aa.start[i]
      current_stop <- data$aa.stop[i]
    }
  }
  # Save the last interval
  merged <- rbind(merged, data.frame(aa.start = current_start, aa.stop = current_stop))
  return(merged)
}

# Apply the function by grouping by gene and db.desc
mfg1<- filtgig1 %>% select(gene,db.desc,aa.start,aa.stop) %>% 
  group_by(gene, db.desc) %>%
  nest() %>%
  mutate(merged_intervals = map(data, merge_intervals)) %>%
  unnest(merged_intervals) %>%
  select(-data) %>% arrange(gene)
fmfg1<-mfg1 %>% inner_join(.,filtgig1[,c(1:3,10)],relationship="many-to-many") %>% distinct() %>% arrange(gene,aa.start)

mfg2<- filtgig2 %>% select(gene,db.desc,aa.start,aa.stop) %>% 
  group_by(gene, db.desc) %>%
  nest() %>%
  mutate(merged_intervals = map(data, merge_intervals)) %>%
  unnest(merged_intervals) %>%
  select(-data)
fmfg2<-mfg2 %>% inner_join(.,filtgig2[,c(1:3,10)],relationship="many-to-many") %>% distinct()


#make custom palette
df_fg2 <- mfg2 %>% group_by(db.desc) %>% dplyr::count() 
df_fg1 <- mfg1 %>% group_by(db.desc) %>% dplyr::count()
palette <- c("#fca311","#108720","#a7c957","#005f99","#fb6f92","#48cae4","#c77dff","#d00000","grey","yellow3")


levels <- union(levels(factor(df_fg1$db.desc)), levels(factor(df_fg2$db.desc)))
names(palette) <- levels


ssag1gn <- fmfg1 %>% ungroup() %>%  
  transmute(file_id="Salmo salar",seq_id=gene,start=aa.start,end=aa.stop,
            strand="+",feat_id=db.desc,source=db.desc,width=aa,name=db.desc,geom_id=db.desc) %>% as.data.frame()
ssag1s<-ssag1gn %>% transmute(file_id="Salmo salar",seq_id,seq_desc=name,length=width) %>% filter(!duplicated(seq_id)) %>% 
  arrange(seq_id) %>% as.data.frame()

ssag2gn<-fmfg2 %>% ungroup() %>% 
  transmute(file_id="Salmo salar",seq_id=gene,start=aa.start,end=aa.stop,
            strand="+",feat_id=db.desc,source=db.desc,width=aa,name=db.desc,geom_id=db.desc) %>% as.data.frame()
ssag2s<-ssag2gn %>% transmute(file_id="Salmo salar",seq_id,seq_desc=name,length=width) %>% filter(!duplicated(seq_id)) %>% 
  arrange(seq_id) %>% as.data.frame()

#make sure the genes are in proper arrangement
g1gn<-clusters %>% filter(gfam=="gig1") %>% distinct(ggroup,gclade,chrom) %>% arrange(gclade,chrom,ggroup)

ssag1gn$seq_id<-factor(ssag1gn$seq_id,levels = unique(rev(g1gn$ggroup)))
ssag1s$seq_id<-factor(ssag1s$seq_id,levels = unique(rev(g1gn$ggroup)))

sg1g<-ssag1gn %>% arrange(seq_id) %>% 
  mutate(galpha=ifelse(name=="SI:CH211-198C19.1-RELATED",1.7,1))
sg1s<-ssag1s %>% arrange(seq_id)

ssag1<-gggenomes(genes = sg1g,seqs=sg1s) + 
  geom_seq() + 
  geom_gene(aes(fill=name),position="pile") +
  geom_bin_label(size=4) +
  scale_fill_manual(values = palette) +
  labs(fill="Predicted functional clusters") 


#gig2
g2gn<-clusters %>% filter(gfam=="gig2") %>% distinct(ggroup,gclade,chrom) %>% arrange(gclade,chrom,ggroup)
ssag2gn$seq_id<-factor(ssag2gn$seq_id,levels = unique(rev(g2gn$ggroup)))
ssag2s$seq_id<-factor(ssag2s$seq_id,levels = unique(rev(g2gn$ggroup)))

sg2g<-ssag2gn %>% arrange(seq_id)
sg2s<-ssag2s %>% arrange(seq_id)

ssag2<-gggenomes(genes = sg2g,seqs=sg2s) + 
  geom_seq() + 
  geom_gene(aes(fill=name),position="pile") + 
  geom_bin_label(size=4) +
  scale_fill_manual(values = palette) +
  labs(fill="Predicted functional cluster") +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank())  


#plots
ssag1
ssag2



#make the gene plots and dotplots in figure 4  

dtpg<-ssag2gn %>% filter(seq_id=="gig2B_ssa02d" | seq_id=="gig2B_ssa12a" | seq_id=="gig2B_ssa12b") %>% arrange(seq_id)
dtps<-ssag2s %>% filter(seq_id=="gig2B_ssa02d" | seq_id=="gig2B_ssa12a" | seq_id=="gig2B_ssa12b") %>% arrange(seq_id)


gggenomes(genes = dtpg,seqs=dtps) + geom_seq() + geom_gene(aes(fill=name),position="pile") +
  geom_bin_label(size=4) +  scale_fill_manual(values = palette) + labs(fill="Predicted functional clusters")

library(dotplot)
library(seqinr)

gig2.inf<-read.fasta("files/gig2_cds_sequence.fa",seqtype = "DNA", set.attributes = FALSE)
gig2<-read.fasta("files/gig2_cds_sequence.fa",seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
this<-gig2.inf[c(which(grepl('ENSSSAG',names(gig2.inf))))]

dotPlot(this$ENSSSAG00000009405,this$ENSSSAG00000009405, wsize = 10,wstep=1,nmatch=9, xlab="gig2B_ssa02d",ylab="gig2B_ssa02d")
dotPlot(this$ENSSSAG00000104005,this$ENSSSAG00000104005, wsize = 10,wstep=1,nmatch=9, xlab="gig2B_ssa12a",ylab="gig2B_ssa12a")
dotPlot(this$ENSSSAG00000083114,this$ENSSSAG00000083114, wsize = 10,wstep=1,nmatch=9, xlab="gig2B_ssa12b",ylab="gig2B_ssa12b")

