#expression heatmaps
library(tidyverse)
library(tools)
library(viridis)
library(ggrepel)
library(patchwork)
library(pheatmap)
library(aplot)
## heatmap annotation
library(ggh4x)
library(ggnewscale)

###############################
#load the functional clusters for the gig genes
cg1s1<-data.table::fread("files/gig1_prottree_split1.txt") %>% mutate(cluster="gig1A") %>% transmute(ensembl_id=V1,cluster)
cg1s2<-data.table::fread("files/gig1_prottree_split2.txt") %>% mutate(cluster="gig1B") %>% transmute(ensembl_id=V1,cluster)
cg2s1<-data.table::fread("files/gig2_prottree_split1.txt") %>% mutate(cluster="gig2A") %>% transmute(ensembl_id=V1,cluster)
cg2s2<-data.table::fread("files/gig2_prottree_split2.txt") %>% mutate(cluster="gig2B") %>% transmute(ensembl_id=V1,cluster)
cg2s3<-data.table::fread("files/gig2_prottree_split3.txt") %>% mutate(cluster="gig2C") %>% transmute(ensembl_id=V1,cluster)

#retain only salmon genes 
clusters<-rbind(cg1s1,cg1s2,cg2s1,cg2s2,cg2s3) %>% filter(grepl('ENSSSAG',ensembl_id)) 

#load the gig genes
gigs<-data.table::fread("files/Ssal_Gig1_Gig2_tss.txt",header=T) %>% 
  mutate(gclade=gsub("_.*","",ggroup),gfam=ifelse(grepl('gig1',gclade),"gig1","gig2")) %>% arrange(gfam,gclade,chrom) %>% 
  filter(tentative_nickname != "gig1_ssa04_e2",ggroup != "gig2_ssa04f")

gigs$gclade<-factor(gigs$gclade,levels=unique(gigs$gclade))
gigs$gfam<-factor(gigs$gfam,levels=unique(gigs$gfam))


gig1<-gigs %>% filter(grepl("gig1",gfam)) %>% transmute(gid=ensembl_id,gig=ggroup)
gig2<-gigs %>% filter(grepl("gig2",gfam)) %>% transmute(gid=ensembl_id,gig=ggroup)


##################################
#### use the SAV3 expression data

counts.tab<-data.table::fread("files/sav3.salmon.merged.gene_counts.tsv") %>% select(-gene_name) %>% as.data.frame()

#add the metadata
meta<-tibble(sample=c("SRR10482313","SRR10482314","SRR10482315","SRR10482316","SRR10482317","SRR10482318","SRR10482319","SRR10482320","SRR10482321","SRR10482322","SRR10482323","SRR10482324"),
             pheno=c("survived","survived","survived","survived","died_11_dpi","died_11_dpi","died_11_dpi","died_11_dpi","survived","survived","died_11_dpi","died_11_dpi") )
##gig1
ctxp<-counts.tab[counts.tab$gene_id %in% gig1$gid,] %>% pivot_longer(contains('SRR'),names_to = "sample",values_to="exp") %>% 
  inner_join(meta,by="sample",relationship="many-to-many") %>% transmute(gene=gene_id,sample,exp,pheno) %>% 
  inner_join(gig1,by=c('gene'='gid'),relationship="many-to-many")
ctdf<-counts.tab[counts.tab$gene_id %in% gig1$gid,] %>% inner_join(gig1,by=c('gene_id'='gid'),relationship="many-to-many") %>% 
  select(-'gene_id') %>% tibble::column_to_rownames(var="gig") %>% as.matrix()
lctdfg1<-log2(ctdf+1)
mtg1<-meta[meta$sample %in% colnames(lctdfg1),] %>% tibble::column_to_rownames(var="sample") %>% arrange(pheno)

##gig2
ctxp<-counts.tab[counts.tab$gene_id %in% gig2$gid,] %>% pivot_longer(contains('SRR'),names_to = "sample",values_to="exp") %>% 
  inner_join(meta,by="sample",relationship="many-to-many") %>% transmute(gene=gene_id,sample,exp,pheno) %>% 
  inner_join(gig2,by=c('gene'='gid'),relationship="many-to-many")
ctdf<-counts.tab[counts.tab$gene_id %in% gig2$gid,] %>% inner_join(gig2,by=c('gene_id'='gid'),relationship="many-to-many") %>% 
  select(-'gene_id') %>% tibble::column_to_rownames(var="gig") %>% as.matrix()
lctdfg2<-log2(ctdf+1)
mtg2<-meta[meta$sample %in% colnames(lctdfg2),] %>% tibble::column_to_rownames(var="sample") %>% arrange(pheno)

#combine everything
bmctab<-rbind(lctdfg1[,rownames(mtg1)],lctdfg2[,rownames(mtg2)])
bmcmt<-rbind(mtg1,mtg2)

#average the expression
dmeta<-bmcmt %>% filter(pheno=="died_11_dpi")
smeta<-bmcmt %>% filter(pheno=="survived")
dxp<-rowMeans(bmctab[,colnames(bmctab) %in% rownames(dmeta)]) %>% as.data.frame() %>% tibble::rownames_to_column(var="gig") %>% rename(.,dead_11_dpi='.')
sxp<-rowMeans(bmctab[,colnames(bmctab) %in% rownames(smeta)]) %>% as.data.frame() %>% tibble::rownames_to_column(var="gig") %>% rename(.,survived='.')

#final table of mean expressions (log2(TPM)+1)
bmc.mexp<-full_join(dxp,sxp,by="gig") #%>% tibble::column_to_rownames(var="gig")


##########################

## AquaFAANG data

gcounts<-data.table::fread("files/polyIC.salmon.merged.gene_counts.tsv") %>% 
  select(-gene_name) %>% as.data.frame()

meta<-data.table::fread("files/polyIC.pbs.metadata") %>% 
  separate(sample_alias,c("ss","uni","vivo","pheno","sample"),"_") %>%
  #filter(grepl("invitro",sample_title)) %>% 
  transmute(sample=run_accession,pheno,vivo) %>% filter(pheno != "vibrio")
pic.counts<-gcounts[,colnames(gcounts) %in% meta$sample] %>% cbind(.,gene=gcounts$gene_id)

#gig1
ctdf<-pic.counts[pic.counts$gene %in% gig1$gid,] %>% inner_join(gig1,by=c('gene'='gid'),relationship="many-to-many") %>% 
  select(-'gene') %>% tibble::column_to_rownames(var="gig") %>% as.matrix()

hkg1df<-log2(ctdf+1)
hkg1mt<-meta[meta$sample %in% colnames(hkg1df),] %>% tibble::column_to_rownames(var="sample") %>% arrange(pheno)

#gig2
ctdf<-pic.counts[pic.counts$gene %in% gig2$gid,] %>% inner_join(gig2,by=c('gene'='gid'),relationship="many-to-many") %>% 
  select(-'gene') %>% tibble::column_to_rownames(var="gig") %>% as.matrix()
hkg2df<-log2(ctdf+1)
hkg2mt<-meta[meta$sample %in% colnames(hkg2df),] %>% tibble::column_to_rownames(var="sample") %>% arrange(pheno)

hktab<-rbind(hkg1df[rownames(hkg1df) %in% gig1$gig,rownames(hkg1mt)],hkg2df[ rownames(hkg2df) %in% gig2$gig,rownames(hkg2mt)]) %>% as.data.frame() 
hkmeta<-rbind(hkg1mt,hkg2mt)


##average the expression::here you need to make four tables
pvv<-hkmeta %>% filter(pheno=="pbs" & vivo=="invivo")
icvv<-hkmeta %>% filter(pheno=="polyIC" & vivo=="invivo")
pvt<-hkmeta %>% filter(pheno=="pbs" & vivo=="invitro")
icvt<-hkmeta %>% filter(pheno=="polyIC" & vivo=="invitro")

pvvxp<-rowMeans(hktab[,colnames(hktab) %in% rownames(pvv)]) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var="gig") %>% rename(.,pbs_invivo='.')
icvvxp<-rowMeans(hktab[,colnames(hktab) %in% rownames(icvv)]) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var="gig") %>% rename(.,polyIC_invivo='.')
pvtxp<-rowMeans(hktab[,colnames(hktab) %in% rownames(pvt)]) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var="gig") %>% rename(.,pbs_invitro='.')
icvtxp<-rowMeans(hktab[,colnames(hktab) %in% rownames(icvt)]) %>% as.data.frame() %>% 
  tibble::rownames_to_column(var="gig") %>% rename(.,polyIC_invitro='.')

#final table of mean expressions (log2(TPM)+1)
hk.mexp<-full_join(pvvxp,icvvxp,by="gig") %>% full_join(pvtxp,by="gig") %>% full_join(icvtxp,by="gig") 



########################################
#assemble everything together
atab<-full_join(bmc.mexp,hk.mexp,by="gig") %>% pivot_longer(-gig,names_to = "pheno",values_to = "mxp") %>% 
  mutate(infection=ifelse(pheno=="dead_11_dpi" | pheno=="survived","SAV3 infection",
                          ifelse(pheno=="pbs_invivo" | pheno== "polyIC_invivo" ,"In-vivo infection","In-vitro infection")),
         pheno=gsub("_.*","",pheno)) %>% 
  mutate(pheno=ifelse(pheno=="pbs","control",pheno))
atab$pheno=factor(atab$pheno,levels=c("survived","dead","control","polyIC"))#,"vibrio"))

########################################
#add regulatory information and adjust the gig names
reg<-data.table::fread("files/gig_2kb_stat_allmtf_deg_dap_master.tsv",header=T)
reg.acc<-reg %>% inner_join(gigs[,c(5,9)],by="ensembl_id",relationship="many-to-many") %>% mutate(gig=ggroup) %>% select(-ggroup)


#isre regulation
isre.regs<-reg.acc %>% filter(isre.mtf!="") %>% select(gig,isre.mtf,isre.dap) %>% mutate(isre.dap=ifelse(isre.dap=="FALSE","no","yes")) %>% 
  group_by(gig) %>% 
  dplyr::count(isre.dap) %>% 
  mutate(isre.motif="yes") %>% select(-n) %>%  group_by(gig) %>%
  arrange(desc(isre.dap)) %>% # Sort rows so "yes" appears first in isre.dap
  dplyr::slice(1) %>% # Take the first instance of each group
  ungroup() %>% transmute(gig,motif=isre.motif,dap=isre.dap) %>%
  right_join(gigs,by=c('gig'='ggroup')) %>% 
  mutate(motif=ifelse(is.na(motif),"no",motif),dap=ifelse(is.na(dap),"no",dap),db="ISRE") %>% 
  select(gig,motif,dap,db) %>% 
  mutate(mdp=ifelse(motif=="yes" & dap=="yes","dap",
                    ifelse(motif=="yes" & dap=="no","yes","no")))

#JASPAR motif regulation
jaspar.regs<-reg.acc %>% filter(motif.seq !="") %>% select(gig,motif.seq,dap) %>% 
  group_by(gig) %>% 
  dplyr::count(dap) %>% 
  mutate(jaspar.motif="yes") %>%
  mutate(dap=ifelse(dap=="FALSE","no","yes")) %>% 
  select(-n) %>% 
  transmute(gig,motif=jaspar.motif,dap) %>% 
  right_join(gigs,by=c('gig'='ggroup')) %>% 
  mutate(motif=ifelse(is.na(motif),"no",motif),dap=ifelse(is.na(dap),"no",dap),db="Jaspar") %>% 
  select(gig,motif,dap,db) %>% 
  mutate(mdp=ifelse(motif=="yes" & dap=="yes","dap",ifelse(motif=="yes" & dap=="no","yes","no")))

regs.all<-bind_rows(isre.regs,jaspar.regs)
regs.all$db<-factor(regs.all$db,levels=c("ISRE","Jaspar")) 
regs.all$mdp<-factor(regs.all$mdp,levels = c("no","yes","dap"))

#ggplot(regs.all) + geom_tile(aes(db,gig,fill=mdp)) + scale_fill_manual(values=c("#ffba08","#e85d04","#d00000"))

regg1<-regs.all %>% filter(grepl('gig1',gig))
regg2<-regs.all %>% filter(grepl('gig2',gig))


###################################
## make the plots

#get the regulatory accessory
reg1p<-ggplot(regg1 %>% filter(!is.na(gig))) + geom_tile(aes(db,gig,fill=mdp),colour="white") + theme_minimal() + 
  scale_fill_manual(values=c("gray","#ff9100","#d00000")) +
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(),axis.text.y = element_text(face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) + labs(fill="motif")

#get the regulatory accessory
reg2p<-ggplot(regg2 %>% filter(!is.na(gig))) + geom_tile(aes(db,gig,fill=mdp),colour="white") + theme_minimal() + 
  scale_fill_manual(values=c("gray","#ff9100","#d00000")) +
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(),axis.text.y = element_text(face="bold"),
      axis.text.x = element_text(angle = 90, hjust = 1)) + labs(fill="motif")


##############################
# plot the average expression for gig1 and gig2 

g1.avg<-atab[atab$gig %in% regg1$gig,] %>% group_by(gig) %>% mutate(mnv =mean(mxp)) %>% distinct(gig,mnv)
g2.avg<-atab[atab$gig %in% regg2$gig,] %>% group_by(gig) %>% mutate(mnv =mean(mxp)) %>% distinct(gig,mnv)

g1.avg.pg<-ggplot(g1.avg) + 
  geom_tile(aes(x="",gig,fill=mnv),color="#F5F5F1") + scale_fill_viridis(option="A",limits=c(0,8)) + 
  labs(x="",fill="avg.expr") + theme_bw() +
  theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())

g2.avg.pg<-ggplot(g2.avg) + 
  geom_tile(aes(x="",gig,fill=mnv),color="#F5F5F1") + scale_fill_viridis(option="A",limits=c(0,8)) + 
  labs(x="",fill="avg.expr") + theme_bw() +  
  theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())


## perform row-scaling: centering (-mean) and dividing by the standard deviation. 
rr.atabg1<-atab[atab$gig %in% regg1$gig,] %>% group_by(gig) %>%
  mutate(scaled_value =(mxp-mean(mxp))/sd(mxp))
rr.atabg1$gig<-factor(rr.atabg1$gig,levels=unique(gigs$ggroup))

rr.atabg2<-atab[atab$gig %in% regg2$gig,] %>% group_by(gig) %>%
  mutate(scaled_value =(mxp-mean(mxp))/sd(mxp)) 
rr.atabg2$gig<-factor(rr.atabg2$gig,levels=unique(gigs$ggroup))



## produce the plot of row scaled expression 
rr.mexp.htmp.gig1<-ggplot(rr.atabg1) + 
  geom_tile(aes(pheno,gig,fill=scaled_value),color="#F5F5F1") + scale_fill_viridis(option="G",limits=c(-3,3)) + xlab(NULL) + ylab(NULL) + 
  scale_x_discrete(position="bottom") + 
  theme( axis.title.x = element_blank(),axis.title.y = element_blank(), axis.ticks.x = element_blank(),
         axis.text.x = element_blank(), legend.title = element_text(angle = 90)) +
  labs(fill="log2(TPM+1)")+
  theme_bw() +
  facet_nested(~ infection + pheno, scales = "free") + 
  theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        panel.spacing = unit(0.05, "lines")) 


rr.mexp.htmp.gig2<-ggplot(rr.atabg2) + 
  geom_tile(aes(pheno,gig,fill=scaled_value),color="#F5F5F1") + scale_fill_viridis(option="G",limits=c(-3,3)) + xlab(NULL) + ylab(NULL) + 
  scale_x_discrete(position="bottom") + 
  theme( axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(), legend.title = element_text(angle = 90)) +
  theme_bw() +
  facet_nested(~ infection + pheno, scales = "free") + 
  theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        panel.spacing = unit(0.05, "lines"))  #strip.text.x = element_text(size = 11, colour = "black",face = "bold" ),strip.background = element_blank(),



#produce all plots 
aplot::plot_list(reg1p,rr.mexp.htmp.gig1, g1.avg.pg,widths = c(0.08,1, 0.04))
aplot::plot_list(reg2p,rr.mexp.htmp.gig2, g2.avg.pg,widths = c(0.08,1, 0.04))



########################################################
##DEG expression of gig1 and gig2 genes

#retrieve DEG information for all genes
degs<-reg.acc %>% select(c(1,20,23,26)) %>% distinct() 

#inspect the conditions under wich each gene is differentialy expressed
degs %>% filter(grepl('gig2',gig)) %>%
  filter(hkt.padj < 0.05)
degs %>% filter(grepl('gig2',gig)) %>%
  filter(hkc.padj < 0.05)
degs %>% filter(grepl('gig2',gig)) %>%
  filter(sav3.padj < 0.05)
degs %>% filter(grepl('gig1',gig)) %>%
  filter(hkt.padj < 0.05)


###############################################
##are all DAP expressing significantly higher?

avg.tab<-atab %>% group_by(gig) %>% mutate(mnv =mean(mxp)) %>% distinct(gig,mnv)
rg<-rbind(regg1,regg2) %>% select(gig,dap) %>% group_by(gig) %>%
  arrange(desc(dap == "yes")) %>%  # Put 'short reads' first if exists
  dplyr::slice(1) %>%                             # Keep the first entry for each group
  ungroup() %>% full_join(avg.tab,by="gig") %>% 
  mutate(dap=ifelse(dap=="yes","overlaps DR region","no DR region overlap"))

deg.rg.tab<-degs %>% mutate(hkt.padj=ifelse(!is.na(hkt.padj),"yes","no"),
                            hkc.padj=ifelse(!is.na(hkc.padj),"yes","no"),
                            sav3.padj=ifelse(!is.na(sav3.padj),"yes","no")) %>%
  mutate(deg=ifelse(hkt.padj=="no" & hkc.padj=="no" & sav3.padj=="no","no DEG","DEG gene")) %>% 
  select(gig,deg) %>% distinct() %>% inner_join(rg,by="gig")

deg.pval<-degs %>% pivot_longer(-gig, names_to = "exp",values_to="pval") %>% 
  mutate(pval=ifelse(!is.na(pval),pval,1)) %>% group_by(gig) %>% 
  filter(pval==min(pval)) %>% 
  filter(!duplicated(gig)) %>% ungroup() %>% 
  mutate(pval=ifelse(pval==1,NA,pval)) %>% 
  select(gig,pval) %>% 
  inner_join(deg.rg.tab,by="gig")
  

#Figure 6: Boxplots of DEG genes, regulation and expression
data <- deg.pval

# Replace missing p-values with 1 (so that -log10(1) = 0)
data <- data %>%
  mutate(pval = ifelse(is.na(pval), 1, as.numeric(pval)),  # Convert missing p-values to 1
         log_pval = -log10(pval),  # Transform p-values to -log10 scale
         mnv = as.numeric(mnv),  # Convert MNV to numeric type
         gig_family = ifelse(str_detect(gig, "gig1"), "Gig1", "Gig2"))  # Classify genes into gig1 and gig2

# Plot creation (Boxplot + Bubble Plot)
ggplot(data, aes(x = gig_family, y = mnv, fill = log_pval, shape = dap)) +
  # Add a background boxplot (without fill and outliers)
  geom_boxplot(aes(group = interaction(gig_family, dap)), 
               outlier.shape = NA, fill = NA, color = "gray30", position = position_dodge(0.8)) +
  # Bubble plot (Jitter to separate by DR Region Overlap)
  geom_jitter(position = position_dodge(0.8), alpha = 0.9,size=3) +
  scale_fill_viridis() +
  scale_shape_manual(values = c(21, 23)) +  # Assign shapes manually
  see::theme_modern()  +  
  labs(title = "Bubble Plot of DEG Genes with Boxplots",
       x = "Gene Family",
       y = "Average gene expression log2(TPM+1)",
       fill = "-log10(p-value)",
       shape = "DR Region Overlap") +
  theme(axis.text.x = element_text(size = 12, face = "bold"))
  


