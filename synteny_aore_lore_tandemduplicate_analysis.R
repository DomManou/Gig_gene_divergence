#gig and Lore/Aore regions taken from the salmobase
library(tidyverse)

gig<-data.table::fread("files/Ssal_Gig1_Gig2_tss.txt",header=T) %>% 
  mutate(gclade=gsub("_.*","",ggroup), 
         gfam=ifelse(grepl('gig1',gclade),"gig1","gig2")) %>% 
  arrange(gfam,gclade,chrom)

# Extract the start and stop of each gene 
gloc<-gig %>% 
  select(chrom, ensembl_start, ensembl_end, defined_tss, ggroup, strand,ensembl_id) %>% 
  transmute(chr=chrom,start=ensembl_start,end=ensembl_end,ggroup,ensembl_id) %>% 
  arrange(chr,start)

#write out the gig file
#data.table::fwrite(gloc,"files/tss_gigs_bedformat.bed",row.names = F,col.names = F,sep="\t")

#fetch synteny data from Salmobasedata 
alore<-data.table::fread("files/synteny_AtlanticSalmon_Ssal_v3.1.tsv",header = T) %>% select(c(2:10)) %>% 
  arrange(chromosome_x,begin_x) %>% mutate(link=row_number())

#make a long list
alorex<-alore %>% select(c(1:3,7:10)) %>% rename_all(~c("chr","start","stop","anchors","naore","nlore","link")) 
alorey<-alore %>% select(c(4:10)) %>% rename_all(~c("chr","start","stop","anchors","naore","nlore","link"))

#write out the data and run BEDtools intersect to locate gigs
#data.table::fwrite(rbind(alorex,alorey),"synteny_AtlanticSalmon_Ssal3.1.bed",col.names = F,row.names = F,sep="\t",quote = F)

###################

#BEDtools intersection

###################

#look at intersection results:: if there are more LORe links, classify as LORe, similarly for AORe
alre<-data.table::fread("files/gig_AORE_LORe_BEDTools_intersection.bed",header=F) %>% 
  mutate(ore=ifelse(V5 > V6,"Aore","Lore"))


#make a table that includes the synteny block ID 
rbind(alorex,alorey) %>% transmute(V1=chr,V2=start,V3=stop,V4=anchors,synteny_blockid=link) %>% inner_join(alre,relationship="many-to-many") %>% 
  transmute(chr=V7,start=V8,stop=V9,gig=V10,gid=V11,ore,synteny_blockid) %>% 
  mutate(gig=gsub("_(?=[a-g])", "", gig,perl=TRUE)) %>% 
  data.table::fwrite(.,"files/gig_lore_aore_syntenyblock.tsv",
                     sep="\t",row.names = F,col.names = T,quote=F)


#from these results you are able to say if gigs are in AORe or LORe regions. 
#filter tandem duplication based on physical proximity (1Mbp from the previous gene)

aldt<-alre %>% select(c(7:12)) %>% 
  rename_all(~c("chr","start","stop","gig","gene","dup"))

check_proximity <- function(df, distance = 1000000) {
  df <- df %>%
    arrange(start) %>%
    mutate(indication = "singleton")
  
  for (i in seq_along(df$start)) {
    for (j in seq_along(df$start)) {
      if (i != j && abs(df$start[i] - df$start[j]) <= distance) {
        df$indication[i] <- "tandem duplicate"
        df$indication[j] <- "tandem duplicate"
      } else if (i != j && abs(df$start[i] - df$stop[j]) <= distance) {
        df$indication[i] <- "tandem duplicate"
        df$indication[j] <- "tandem duplicate"
      } else if (i != j && abs(df$stop[i] - df$start[j]) <= distance) {
        df$indication[i] <- "tandem duplicate"
        df$indication[j] <- "tandem duplicate"
      } else if (i != j && abs(df$stop[i] - df$stop[j]) <= distance) {
        df$indication[i] <- "tandem duplicate"
        df$indication[j] <- "tandem duplicate"
      }
    }
  }
  return(df)
}

# Group by chromosome and apply the proximity check function
result <- aldt %>%
  group_by(chr) %>%
  nest() %>%
  mutate(data = map(data, check_proximity)) %>%
  unnest(data) %>% 
  mutate(gfam=ifelse(grepl('gig1',gig),"gig1","gig2"))

# Print result
print(result)

#data.table::fwrite(result,"files/gig_wd_LORe_AORe_td.tsv",row.names = F,col.names = T,quote=F,sep="\t")



