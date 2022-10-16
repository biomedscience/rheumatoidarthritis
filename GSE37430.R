library("tidyverse")
library("data.table")
library("metafor")
library("GEOquery")
library("ggplot2")
library("ggbeeswarm")
library("egg")
library("tidyr")
library("tsne")
library("plotly")

GSE37430 <- getGEO("GSE37430", GSEMatrix = TRUE,AnnotGPL=TRUE)
phen <- pData(phenoData(GSE37430[[1]]))

FKPM <- as.data.frame(exprs(GSE37430[[1]]))
head(FKPM)

probe<-fData(GSE37430[[1]]) %>% select(ID,"Gene symbol")
head(probe)

input<-data.frame(FKPM,ID=probe[match(rownames(FKPM),probe$ID),2])
input$MIF=input$GSM919192/input$GSM919191
input$IL8=input$GSM919193/input$GSM919191
input %>% filter(grepl("IL",ID)) %>% arrange(desc(MIF)) %>% head(n=20)
input %>% filter(grepl("VEG",ID)) %>% arrange(desc(MIF)) %>% head(n=20)
input %>% filter(grepl("S1P",ID)) %>% arrange(desc(MIF)) %>% head(n=20)

write.csv(input,file="GSE37430.csv",quote=F)
