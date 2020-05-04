# clusterProfiler of GO analysis for Oryza sativa in R
#安装包
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
#获得OrgDb
#通过AnnotationHub在线检索并抓取OrgDb
require(AnnotationHub)
hub<-AnnotationHub()
query(hub,"Oryza_sativa")
#title                                               
  AH10561 | hom.Oryza_sativa.inp8.sqlite                        
  AH75915 | org.Oryza_sativa_(japonica_cultivar-group).eg.sqlite
  AH75916 | org.Oryza_sativa_Japonica_Group.eg.sqlite           
  AH75917 | org.Oryza_sativa_subsp._japonica.eg.sqlite 
#通过检索， org.Oryza_sativa_Japonica_Group.eg.sqlite就是我们所要的OrgDb，可以通过相应的accession number, AH75916抓取文件，并存入了rice对象中，它包含了51097个基因的注释：
rice<-hub[['AH75916']]
length(keys(rice))[1]
#[1] 35257
columns(rice) 
这个OrgDb，包含有以下一些注释信息:
#[1] "ACCNUM"      "ALIAS"       "CHR"         "ENTREZID"    "EVIDENCE"    "EVIDENCEALL"
#[7] "GENENAME"    "GID"         "GO"          "GOALL"       "ONTOLOGY"    "ONTOLOGYALL"
#[13] "PMID"        "REFSEQ"      "SYMBOL"   
head(keys(rice,keytype = "SYMBOL"))  #查看注释的类型和例子
#导入已经差异分析过的数据
setwd("C:/Users/Jin qiongli/Documents")
getwd()
read.csv("your_csv_file.csv")
d<-read.csv("your_csv_file.csv")
gene=as.character(d[,1])
head(gene)
#[1] "Os07g0684000" "Os04g0598000" "Os09g0433800" "Os02g0519700" "Os02g0735200" "Os07g0615200"

library("riceidconverter")
convert_id<-RiceIDConvert(myID = gene,'RAP',toType = 'SYMBOL')  #通过查看注释的类型，发现不能转换RAPto GO,所以得先得到SYMBOL，再转GO
convert_id  #查看自己转的对不对
#              RAP       SYMBOL
#1    EPlOSAG00000000424         <NA>
#2    EPlOSAG00000003207         <NA>
#3    EPlOSAG00000003214         <NA>
#发现数据中有大量的NA和None SYMBOL，需要去除
  convert_id=na.omit(convert_id)  #去除含有NA的行 
geneList=as.character(convert_id[,2])  
head(geneList)
# 重要错误点：
> require(clusterProfiler)
> bitr(keys(convert_id)[2], 'SYMBOL', c( "GO", "ONTOLOGY"), rice)
# Error in (function (classes, fdef, mtable)  : 
# unable to find an inherited method for function ‘keys’ for signature ‘"data.frame"’
> library(clusterProfiler)
> bitr(keys(geneList), 'SYMBOL', c("REFSEQ", "GO", "ONTOLOGY"), rice) 
#Error in (function (classes, fdef, mtable)  : 
#unable to find an inherited method for function ‘keys’ for signature ‘"character"’
#正确操作：上述的情况的出现是因为keys这个函数不仅再AnnotationHub中有，也用于ClusterProfiler.所以需要安装ClusterProfiler。
install.packages("clusterProfiler")
library("clusterProfiler")
bitr(keys(geneList)[1], 'SYMBOL', c("REFSEQ", "GO", "ONTOLOGY"), rice)
#错误: external pointer is not valid     
#唉，又碰到了新问题。
 
  
library(clusterProfiler)
bitr(keys(rice)[1], 'ENTREZID', c("REFSEQ", "GO", "ONTOLOGY"), rice)
sample_genes <- keys(rice)[1:100]
head(sample_genes)
require(clusterProfiler)
res = enrichGO(sample_genes, OrgDb=rice, pvalueCutoff=1, qvalueCutoff=1)
res
dotplot(res)
barplot(res,showCategory = 15)
getwd()
list.files("C:/Users/Jin qiongli/Documents")
read.csv("your_csv_file.csv")
d<-read.csv("your_csv_file.csv")
geneList<-d[,3]
names(geneList)<-as.character(d[,1])
geneList<-sort(geneList,decreasing = TRUE)
data(geneList,package = "DOSE")
head(geneList)
rice<-hub[['AH75916']]
library(clusterProfiler)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(keys(gene),  "ENTREZID",c( "GO"),rice)
head(gene.df)
