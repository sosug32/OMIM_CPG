library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)
library(tidyverse)  # data manipulation
require(tclust)    # clustering algorithms
library(factoextra) 
setwd("~/Documents/clustering/")

c.input <- fread(inputdata.clustering, data.table = F)
c.input$OR.gene <- as.numeric(c.input$OR.gene)
c.input$OR.LOH <- as.numeric(c.input$OR.LOH)
c.input %>% dplyr::group_by(Group) %>% dplyr::summarise(n = n_distinct(VEP_HGNC))
unique(c.input$VEP_HGNC) %>% length()
c.input$pancan.logor.gene <- log(c.input$OR.gene,2)
c.input$pancan.logor.loh <- log(c.input$OR.LOH,2)

tmpp <- c.input %>% unique() 
dt1 <- tmpp[,c(1,10,7,8)] ###7. OR.LOH + o/e value + TAU expression values
summary(dt1)
dt11 <- dt1[,2:4]
str(dt11)
pca_res  <- prcomp(dt11, scale. = T)
summary(pca_res)
plot(pca_res, type = "l")
vars <- apply(pca_res$x, 2, var) 
props <- vars / sum(vars)
props
dev.off()
pca_dt1r <- pca_res$x %>% data.frame()
pca_dt1r$"VEP_HGNC" <- dt1$VEP_HGNC
head(pca_dt1r)
pca_dt1r <- merge(pca_dt1r, tmpp[,c(1,4)], by = "VEP_HGNC", all.x = T)

set.seed(20)
k.max <- 10
head(pca_dt1r)
data <- pca_dt1r[, 2:4] ##PCA clustering input columns
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50, iter.max = 15 )$tot.withinss})
wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
dev.off()

dff <- scale(pca_dt1r[, 2:4])
e_clust <- eclust(dff, "kmeans", hc_metric="eucliden", k=4)

ggdata2 <- data.frame(pca_dt1r, eCluster=e_clust$cluster) 
ggdata2$eCluster <- as.factor(ggdata2$eCluster)
ncol(ggdata2)
head(ggdata2)

final.rdgv5 <- merge(tmpp, ggdata2[,c(1,2,3,4,6)], by = "VEP_HGNC",all.x = T)
final.rdgv5 %>% group_by(eCluster) %>% dplyr::summarise(gene.count = n_distinct(VEP_HGNC))
write.table(final.rdgv5, "PCA_Clustering_loh_oe_tau.tsv", sep = "\t")
