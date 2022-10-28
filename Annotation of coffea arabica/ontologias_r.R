##load libraries
#

library(readr)
library(tidyverse)
library(janitor)

blast2go_table <- read_delim("blast2go_table.txt",
delim = "\t", escape_double = FALSE,
trim_ws = TRUE)


blast2go_annot <- read_table("blast2go_annot.annot")


ncbi_dataset_coffea <- read_delim("D:/tesis cafe/blast2go/ncbi_dataset.tsv", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

res_37_vs_control.lfc<-read.csv("D:/tesis cafe/DESeq2 coffea/res_37_vs_control.1.lfc.csv",header = T)
#----making an txdb object
library(GenomicFeatures)

# txdb<-makeTxDbFromGFF("C:/Users/monic/OneDrive/Documentos/data tesis piero/GCF_003713225.1_Cara_1.0_genomic.gff",format = "gff3",organism = "Coffea arabica",taxonomyId = 13443)
# 
# saveDb(txdb,file="TxDb.Carabica.cara01.sqlite")

txdb <- loadDb("TxDb.Carabica.cara01.sqlite")

trans<-transcripts(txdb)

##work on gff3 file from ncbi
# library(ape)
# 
# gff_coffea<-read.gff("D:/tesis cafe/blast2go/GCF_003713225.1_Cara_1.0_genomic.gff",GFF3 = T)
# 
# gff_coffea<-as_tibble(gff_coffea)%>%
#   filter(type=="gene")
# 
# get_geneids<-function(x){
#   gene_id<-str_extract(x,"GeneID:(\\w*?);")
#   return(gene_id)
# }
# 
# get_genenames<-function(x){
#   gene_name<-str_extract(x,"ID=(.*?);")
#   return(gene_name)
# }
# 
# gene_ids<-map_dfr(gff_coffea,get_geneids)
# gene_ids<-str_extract(gene_ids$attributes,"[0-9]+")
# 
# gene_names<-map_dfr(gff_coffea,get_genenames)
# 
# gene_names<-gsub(gsub(gene_names$attributes,pattern = "ID=gene-",replacement = ""),pattern = ";",replacement = "")
# 
# 
# gene_ids_names<-bind_cols(ENTREZID=gene_ids,
#                           gene_names=gene_names,
#                           seqid=gff_coffea$seqid,
#                           start=gff_coffea$start,
#                           end=gff_coffea$end,
#                           strand=gff_coffea$strand)
# 
# #Get orgdb from annotationhub
# library(AnnotationHub)
# 
# ah<-AnnotationHub()
# 
# ah_Carabica<-subset(ah,species == "Coffea arabica")
# 
# query(ah_Carabica,"org.C")
# 
# Carabica<-ah[["AH100848"]]
# 
# keys_1<-keys(Carabica)
# 
# #columns(Carabica)
# 
# info_cara<-select(Carabica,keys = keys_1,keytype = "ENTREZID",columns = c("GID","REFSEQ","SYMBOL","GENENAME","ACCNUM"))
# 
# 
# #----extract genenames from results rsem
# tx_and_genes <- read_table("D:/tesis cafe/blast2go/SRR11196524_.genes.results")
# 
# tx_and_genes<-tx_and_genes%>%separate_rows(`transcript_id(s)`,sep = ",")
# 
# tx_and_genes<-tx_and_genes%>%
#   dplyr::select(gene_id,`transcript_id(s)`)
# 
# colnames(tx_and_genes)<-c("gene_names","tx_id")


#---------adding genes to blast2go table
table_with_genes<-blast2go_table%>%
  clean_names()%>%
  distinct()%>%
  dplyr::select(!number_go:inter_pro_go_names)

ncbi_dataset_coffea<-ncbi_dataset_coffea%>%
  clean_names()

table(table_with_genes$seq_name%in%ncbi_dataset_coffea$transcript_accession)



table_with_genes<-table_with_genes%>%
  right_join(ncbi_dataset_coffea,by=c("seq_name"="transcript_accession"))

table_with_genes<-table_with_genes%>%
  dplyr::select(!c(transcript_transcript_name,
                   nomenclature_id,
                   ensembl_gene_i_ds,
                   omim_i_ds,
                   protein_accession,
                   protein_name,
                   transcript_protein_accession,
                   transcript_protein_name,
                   synonyms,
                   swiss_prot_accessions))

#----------blast2go_table
library(data.table)

table_with_genes<-data.table(table_with_genes)



coffea_sym<- table_with_genes%>%
  dplyr::select(ncbi_gene_id,symbol,description.y,gene_type,chromosomes,
                transcript_genomic_accession,transcript_genomic_start,
                transcript_genomic_stop,orientation,seq_name,description.x,
                length)


colnames(coffea_sym)<-c("GID","SYMBOL","GENENAMENCBI","GENETYPE",
                        "CHROMOSOMES","TXGENOMICACCNUM","START",
                        "END","STRAND","TXNAME","GENENAMEBLAST",
                        "GENELENGTH")

coffea_sym<-coffea_sym%>%
  map_dfr(as.character)

apply(coffea_sym,2,function(x){
  sum(is.na(x))
})

coffea_sym<-coffea_sym%>%
  dplyr::select(!c(TXGENOMICACCNUM:GENELENGTH))
  
colnames(coffea_sym)

coffea_sym<-coffea_sym%>%
  distinct()

#--------blast2go_annot

blast2go_annot<- blast2go_annot%>%
  clean_names()

coffea_go<-blast2go_annot%>%
  pivot_longer(name:id,values_to = "GO",names_to = "names")%>%
  dplyr::select(!names)

coffea_go<-coffea_go%>%
  mutate(EVIDENCE="IEA")

dummy_df<-as_tibble(ncbi_dataset_coffea)%>%
  dplyr::select(ncbi_gene_id,transcript_accession)

colnames(coffea_go)<-c("TXNAME","GO","EVIDENCE")

coffea_go<-coffea_go%>%
  right_join(dummy_df,by=c("TXNAME"="transcript_accession"))

coffea_go<-coffea_go%>%
  dplyr::select(GID=ncbi_gene_id,GO,EVIDENCE)
  
apply(coffea_go,2,function(x){
  sum(is.na(x))
})

coffea_go<-coffea_go%>%
  na.omit()

library(AnnotationForge)


coffea_go<-coffea_go%>%
  distinct()%>%
  data.table()

coffea_go<-coffea_go%>%
  map_dfr(as.character)




makeOrgPackage(gene_info=coffea_sym,
               go=coffea_go,
               version = "0.1",
               maintainer = "Piero Palacios Bernuy <p.palacios.bernuy@gmail.com>",
               author = "Piero Palacios Bernuy",
               outputDir = ".",
               tax_id = "13443",
               genus = "Coffea",
               species = "arabica",
               goTable="go")

install.packages("./org.Carabica.eg.db",repos = NULL,type = "source")

library(org.Carabica.eg.db)

columns(org.Carabica.eg.db)
keys<- (keys(org.Carabica.eg.db))

an<-AnnotationDbi::select(org.Carabica.eg.db,keys = keys,columns = c("GID","SYMBOL","CHROMOSOMES"))


#-----------
#Creting a package for Coffea arabica

library(GO.db)
library(OrganismDbi)


gd = list( join1 = c(GO.db="GOID", org.Carabica.eg.db="GO"),
           join2 = c(org.Carabica.eg.db="GID",
                     txdb="GENEID"))
if (!file.exists("Cara1.0")) # don't do twice...
  makeOrganismPackage(pkgname="Cara1.0",  # simplify typing!
                      graphData=gd, organism="Coffea arabica",
                      version="1.0.0", maintainer="Piero Palacios <p.palacios.bernuy@gmail.com>",
                      author="Piero Palacios <p.palacios.bernuy@gmail.com>",
                      destDir=".",
                      license="Artistic-2.0")

install.packages("Cara1.0", repos=NULL, type="source")

library(Cara1.0)

columns(Cara1.0)

n<-(keys(Cara1.0,keytype="GID"))

b<-select(Cara1.0,keys = n,keytype = "GID",columns = c("SYMBOL","TERM"))
b




