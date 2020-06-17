install.packages("Rcpp")
library(Rcpp)
library(readxl)
library(edgeR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(UpSetR)
library(ggpubr)


#import the dataset

biopanning_r1 <- read_excel("D:/Dropbox/2-laboratorija_ongoing/1_phd-LU/2019-1-gads/3_ERASMUS_tambet_teesalu/6_protocols-data-results/9-biopanning/Kaarel Kurm - round1_allpeptides.xlsx", 
                            na = "0", skip = 5)

biopanning_r2 <- read_excel("D:/Dropbox/2-laboratorija_ongoing/1_phd-LU/2019-1-gads/3_ERASMUS_tambet_teesalu/6_protocols-data-results/9-biopanning/Kaarel Kurm - round2_allpeptides.xlsx", 
                            na = "0", skip = 5)
biopanning_r2_2 <- read_excel("D:/Dropbox/2-laboratorija_ongoing/1_phd-LU/2019-1-gads/3_ERASMUS_tambet_teesalu/6_protocols-data-results/9-biopanning/Kaarel Kurm - 20191114_Maarja_allpeptides.xlsx", 
                              na = "0", skip = 5)
#replace NAs with 0

biopanning_r1[is.na(biopanning_r1)] <- 0
biopanning_r2[is.na(biopanning_r2)] <- 0
biopanning_r2_2[is.na(biopanning_r2_2)] <- 0

biopanning_R2 <- full_join(biopanning_r2, biopanning_r2_2, by = "Peptides")
biopanning_all <- full_join(biopanning_r1, biopanning_R2, by = "Peptides")

biopanning_all[is.na(biopanning_all)] <- 0
biopanning_R2[is.na(biopanning_R2)] <- 0

View(biopanning_all)

#create DGEList

biopanning_all_dge <- DGEList(counts = biopanning_all[,2:55], genes = biopanning_all[,1])
biopanning_R2_dge_brain <- DGEList(counts = biopanning_R2[2:19], genes = biopanning_R2[,1])
biopanning_R2_dge_lung <- DGEList(counts = biopanning_R2[20:37], genes = biopanning_R2[,1])


keep_R2_brain_dge <- rowSums((biopanning_R2_dge_brain$counts)>5)>= 2
biopanning_R2_dge_brain_filtered <- biopanning_R2_dge_brain[keep_R2_brain_dge, , keep.lib.sizes = FALSE]

biopanning_R2_brain_dge_norm <- calcNormFactors(biopanning_R2_dge_brain_filtered)

keep_R2_lung_dge <- rowSums((biopanning_R2_dge_lung$counts)>5)>= 2
biopanning_R2_dge_lung_filtered <- biopanning_R2_dge_lung[keep_R2_lung_dge, , keep.lib.sizes = FALSE]

biopanning_R2_lung_dge_norm <- calcNormFactors(biopanning_R2_dge_lung_filtered)

plotMDS(biopanning_all_dge)
plotMDS(biopanning_R2_brain_dge_norm, top=20)
plotMDS(biopanning_R2_lung_dge_norm, top=20)
##round 1 analysis
#----------
#create datasets

lib_vs_li <- biopanning_r1[,c(1,2,8,14,5,11,17)]
lib_vs_br <- biopanning_r1[,c(1,2,8,14,3,9,15)]
lib_vs_ki <- biopanning_r1[,c(1,2,8,14,6,12,18)]
lib_vs_lu <- biopanning_r1[,c(1,2,8,14,4,10,16)]
lib_vs_mus <- biopanning_r1[,c(1,2,8,14,7,13,19)]

View(lib_vs_mus)

#create groups
group_biopanning <- factor(c(1,1,1,2,2,2))

#create model matrix
biopanning_design <- model.matrix(~group_biopanning)

#make dge lists

lib_vs_li_dge <- DGEList(counts = lib_vs_li[,2:7], genes = lib_vs_li[,1], group = group_biopanning)
lib_vs_br_dge <- DGEList(counts = lib_vs_br[,2:7], genes = lib_vs_br[,1], group = group_biopanning)
lib_vs_ki_dge <- DGEList(counts = lib_vs_ki[,2:7], genes = lib_vs_ki[,1], group = group_biopanning)
lib_vs_lu_dge <- DGEList(counts = lib_vs_lu[,2:7], genes = lib_vs_lu[,1], group = group_biopanning)
lib_vs_mus_dge <- DGEList(counts = lib_vs_mus[,2:7], genes = lib_vs_mus[,1], group = group_biopanning)

nrow(lib_vs_li_dge$genes)

#filter out low abundance peptides

keep_lib_vs_li_dge <- rowSums((lib_vs_li_dge$counts)>5)>= 2
lib_vs_li_dge_filtered <- lib_vs_li_dge[keep_lib_vs_li_dge, , keep.lib.sizes = FALSE]

keep_lib_vs_br_dge <- rowSums((lib_vs_br_dge$counts)>5)>= 2
lib_vs_br_dge_filtered <- lib_vs_br_dge[keep_lib_vs_br_dge, , keep.lib.sizes = FALSE]

keep_lib_vs_ki_dge <- rowSums((lib_vs_ki_dge$counts)>5)>= 2
lib_vs_ki_dge_filtered <- lib_vs_ki_dge[keep_lib_vs_ki_dge, , keep.lib.sizes = FALSE]

keep_lib_vs_lu_dge <- rowSums((lib_vs_lu_dge$counts)>5)>= 2
lib_vs_lu_dge_filtered <- lib_vs_lu_dge[keep_lib_vs_lu_dge, , keep.lib.sizes = FALSE]

keep_lib_vs_mus_dge <- rowSums((lib_vs_mus_dge$counts)>5)>= 2
lib_vs_mus_dge_filtered <- lib_vs_mus_dge[keep_lib_vs_mus_dge, , keep.lib.sizes = FALSE]

nrow(lib_vs_br_dge_filtered)

plotMDS(lib_vs_br_dge_filtered,top = 20)

#normalize samples
lib_vs_li_dge_filtered_norm <- calcNormFactors(lib_vs_li_dge_filtered)

lib_vs_br_dge_filtered_norm <- calcNormFactors(lib_vs_br_dge_filtered)

lib_vs_ki_dge_filtered_norm <- calcNormFactors(lib_vs_ki_dge_filtered)

lib_vs_lu_dge_filtered_norm <- calcNormFactors(lib_vs_lu_dge_filtered)

lib_vs_mus_dge_filtered_norm <- calcNormFactors(lib_vs_mus_dge_filtered)

#estimate dispersions, fit, perform stat tests and save results

lib_vs_li_de <- estimateGLMCommonDisp(lib_vs_li_dge_filtered_norm, biopanning_design)
lib_vs_li_de <- estimateGLMTrendedDisp(lib_vs_li_de, biopanning_design)
lib_vs_li_de <- estimateGLMTagwiseDisp(lib_vs_li_de, biopanning_design)
lib_vs_li_fit <- glmQLFit(lib_vs_li_de, biopanning_design)
lib_vs_li_test <-glmQLFTest(lib_vs_li_fit)
lib_vs_li_res <- topTags(lib_vs_li_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
lib_vs_li_res_lfc2 <- lib_vs_li_res %>% 
  filter(lib_vs_li_res$logFC > 2)
View(lib_vs_li_res_lfc2)

lib_vs_br_de <- estimateGLMCommonDisp(lib_vs_br_dge_filtered_norm, biopanning_design)
lib_vs_br_de <- estimateGLMTrendedDisp(lib_vs_br_de, biopanning_design)
lib_vs_br_de <- estimateGLMTagwiseDisp(lib_vs_br_de, biopanning_design)
lib_vs_br_fit <- glmQLFit(lib_vs_br_de, biopanning_design)
lib_vs_br_test <-glmQLFTest(lib_vs_br_fit)
lib_vs_br_res <- topTags(lib_vs_br_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
lib_vs_br_res_lfc2 <- lib_vs_br_res %>% 
  filter(lib_vs_br_res$logFC > 2)
View(lib_vs_br_res_lfc2)
View(lib_vs_br_res)

lib_vs_ki_de <- estimateGLMCommonDisp(lib_vs_ki_dge_filtered_norm, biopanning_design)
lib_vs_ki_de <- estimateGLMTrendedDisp(lib_vs_ki_de, biopanning_design)
lib_vs_ki_de <- estimateGLMTagwiseDisp(lib_vs_ki_de, biopanning_design)
lib_vs_ki_fit <- glmQLFit(lib_vs_ki_de, biopanning_design)
lib_vs_ki_test <-glmQLFTest(lib_vs_ki_fit)
lib_vs_ki_res <- topTags(lib_vs_ki_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
lib_vs_ki_res_lfc2 <- lib_vs_ki_res %>% 
  filter(lib_vs_ki_res$logFC > 2)
View(lib_vs_ki_res_lfc2)

lib_vs_lu_de <- estimateGLMCommonDisp(lib_vs_lu_dge_filtered_norm, biopanning_design)
lib_vs_lu_de <- estimateGLMTrendedDisp(lib_vs_lu_de, biopanning_design)
lib_vs_lu_de <- estimateGLMTagwiseDisp(lib_vs_lu_de, biopanning_design)
lib_vs_lu_fit <- glmQLFit(lib_vs_lu_de, biopanning_design)
lib_vs_lu_test <-glmQLFTest(lib_vs_lu_fit)
lib_vs_lu_res <- topTags(lib_vs_lu_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
lib_vs_lu_res_lfc2 <- lib_vs_lu_res %>% 
  filter(lib_vs_lu_res$logFC > 2)
View(lib_vs_lu_res_lfc2)
write.csv(lib_vs_lu_res_lfc2, file = "Round1_lib_vs_lung_enrichment_lfc2.csv")
View(lib_vs_lu_res)

lib_vs_mus_de <- estimateGLMCommonDisp(lib_vs_mus_dge_filtered_norm, biopanning_design)
lib_vs_mus_de <- estimateGLMTrendedDisp(lib_vs_mus_de, biopanning_design)
lib_vs_mus_de <- estimateGLMTagwiseDisp(lib_vs_mus_de, biopanning_design)
lib_vs_mus_fit <- glmQLFit(lib_vs_mus_de, biopanning_design)
lib_vs_mus_test <-glmQLFTest(lib_vs_mus_fit)
lib_vs_mus_res <- topTags(lib_vs_mus_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
lib_vs_mus_res_lfc2 <- lib_vs_mus_res %>% 
  filter(lib_vs_mus_res$logFC > 2)
View(lib_vs_mus_res_lfc2)
write.csv(lib_vs_mus_res)


###DO not collapse----END of 1st round----------


###round 2 analysis

### round 2 brain
###----

###round 2 brain - enrichment

#create datasets

br2_lib_vs_li <- biopanning_all[,c(1,2,8,14,23,29,35)]
br2_lib_vs_br <- biopanning_all[,c(1,2,8,14,21,27,33)]
br2_lib_vs_ki <- biopanning_all[,c(1,2,8,14,24,30,36)]
br2_lib_vs_lu <- biopanning_all[,c(1,2,8,14,22,28,34)]
br2_lib_vs_mus <- biopanning_all[,c(1,2,8,14,25,31,37)]

View(br2_lib_vs_ki)

#create groups
group_biopanning <- factor(c(1,1,1,2,2,2))

#create model matrix
biopanning_design <- model.matrix(~group_biopanning)

#make dge lists

br2_lib_vs_li_dge <- DGEList(counts = br2_lib_vs_li[,2:7], genes = br2_lib_vs_li[,1], group = group_biopanning)
br2_lib_vs_br_dge <- DGEList(counts = br2_lib_vs_br[,2:7], genes = br2_lib_vs_br[,1], group = group_biopanning)
br2_lib_vs_ki_dge <- DGEList(counts = br2_lib_vs_ki[,2:7], genes = br2_lib_vs_ki[,1], group = group_biopanning)
br2_lib_vs_lu_dge <- DGEList(counts = br2_lib_vs_lu[,2:7], genes = br2_lib_vs_lu[,1], group = group_biopanning)
br2_lib_vs_mus_dge <- DGEList(counts = br2_lib_vs_mus[,2:7], genes = br2_lib_vs_mus[,1], group = group_biopanning)

nrow(br2_lib_vs_li_dge$genes)

#filter out low abundance peptides

keep_br2_lib_vs_li_dge <- rowSums((br2_lib_vs_li_dge$counts)>5)>= 2
br2_lib_vs_li_dge_filtered <- br2_lib_vs_li_dge[keep_br2_lib_vs_li_dge, , keep.lib.sizes = FALSE]

keep_br2_lib_vs_br_dge <- rowSums((br2_lib_vs_br_dge$counts)>5)>= 2
br2_lib_vs_br_dge_filtered <- br2_lib_vs_br_dge[keep_br2_lib_vs_br_dge, , keep.lib.sizes = FALSE]

keep_br2_lib_vs_ki_dge <- rowSums((br2_lib_vs_ki_dge$counts)>5)>= 2
br2_lib_vs_ki_dge_filtered <- br2_lib_vs_ki_dge[keep_br2_lib_vs_ki_dge, , keep.lib.sizes = FALSE]

keep_br2_lib_vs_lu_dge <- rowSums((br2_lib_vs_lu_dge$counts)>5)>= 2
br2_lib_vs_lu_dge_filtered <- br2_lib_vs_lu_dge[keep_br2_lib_vs_lu_dge, , keep.lib.sizes = FALSE]

keep_br2_lib_vs_mus_dge <- rowSums((br2_lib_vs_mus_dge$counts)>5)>= 2
br2_lib_vs_mus_dge_filtered <- br2_lib_vs_mus_dge[keep_br2_lib_vs_mus_dge, , keep.lib.sizes = FALSE]

nrow(br2_lib_vs_br_dge_filtered)

plotMDS(br2_lib_vs_br_dge_filtered,top = 20)

#normalize samples
br2_lib_vs_li_dge_filtered_norm <- calcNormFactors(br2_lib_vs_li_dge_filtered)

br2_lib_vs_br_dge_filtered_norm <- calcNormFactors(br2_lib_vs_br_dge_filtered)

br2_lib_vs_ki_dge_filtered_norm <- calcNormFactors(br2_lib_vs_ki_dge_filtered)

br2_lib_vs_lu_dge_filtered_norm <- calcNormFactors(br2_lib_vs_lu_dge_filtered)

br2_lib_vs_mus_dge_filtered_norm <- calcNormFactors(br2_lib_vs_mus_dge_filtered)

#estimate dispersions, fit, perform stat tests and save results

###IT IS NECESSARY TO ADD logFC calculations


br2_lib_vs_li_de <- estimateGLMCommonDisp(br2_lib_vs_li_dge_filtered_norm, biopanning_design)
br2_lib_vs_li_de <- estimateGLMTrendedDisp(br2_lib_vs_li_de, biopanning_design)
br2_lib_vs_li_de <- estimateGLMTagwiseDisp(br2_lib_vs_li_de, biopanning_design)
br2_lib_vs_li_fit <- glmQLFit(br2_lib_vs_li_de, biopanning_design)
br2_lib_vs_li_test <-glmQLFTest(br2_lib_vs_li_fit)
br2_lib_vs_li_res <- topTags(br2_lib_vs_li_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
br2_lib_vs_li_res_lfc2 <- br2_lib_vs_li_res %>% 
  filter(br2_lib_vs_li_res$logFC > 2)
View(br2_lib_vs_li_res_lfc2)

br2_lib_vs_br_de <- estimateGLMCommonDisp(br2_lib_vs_br_dge_filtered_norm, biopanning_design)
br2_lib_vs_br_de <- estimateGLMTrendedDisp(br2_lib_vs_br_de, biopanning_design)
br2_lib_vs_br_de <- estimateGLMTagwiseDisp(br2_lib_vs_br_de, biopanning_design)
br2_lib_vs_br_fit <- glmQLFit(br2_lib_vs_br_de, biopanning_design)
br2_lib_vs_br_test <-glmQLFTest(br2_lib_vs_br_fit)
br2_lib_vs_br_res <- topTags(br2_lib_vs_br_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
br2_lib_vs_br_res_lfc2 <- br2_lib_vs_br_res %>% 
  filter(br2_lib_vs_br_res$logFC > 2)
View(br2_lib_vs_br_res_lfc2)
View(br2_lib_vs_br_res)
write.csv(br2_lib_vs_br_res_lfc2, file = "Round2_lib_vs_brain_enrichment_lfc2.csv")

br2_lib_vs_ki_de <- estimateGLMCommonDisp(br2_lib_vs_ki_dge_filtered_norm, biopanning_design)
br2_lib_vs_ki_de <- estimateGLMTrendedDisp(br2_lib_vs_ki_de, biopanning_design)
br2_lib_vs_ki_de <- estimateGLMTagwiseDisp(br2_lib_vs_ki_de, biopanning_design)
br2_lib_vs_ki_fit <- glmQLFit(br2_lib_vs_ki_de, biopanning_design)
br2_lib_vs_ki_test <-glmQLFTest(br2_lib_vs_ki_fit)
br2_lib_vs_ki_res <- topTags(br2_lib_vs_ki_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
br2_lib_vs_ki_res_lfc2 <- br2_lib_vs_ki_res %>% 
  filter(br2_lib_vs_ki_res$logFC > 2)
View(br2_lib_vs_ki_res_lfc2)

br2_lib_vs_lu_de <- estimateGLMCommonDisp(br2_lib_vs_lu_dge_filtered_norm, biopanning_design)
br2_lib_vs_lu_de <- estimateGLMTrendedDisp(br2_lib_vs_lu_de, biopanning_design)
br2_lib_vs_lu_de <- estimateGLMTagwiseDisp(br2_lib_vs_lu_de, biopanning_design)
br2_lib_vs_lu_fit <- glmQLFit(br2_lib_vs_lu_de, biopanning_design)
br2_lib_vs_lu_test <-glmQLFTest(br2_lib_vs_lu_fit)
br2_lib_vs_lu_res <- topTags(br2_lib_vs_lu_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
br2_lib_vs_lu_res_lfc2 <- br2_lib_vs_lu_res %>% 
  filter(br2_lib_vs_lu_res$logFC > 2)
View(br2_lib_vs_lu_res_lfc2)
View(br2_lib_vs_lu_res)
write.csv(br2_lib_vs_lu_res_lfc2, file = "br2_lib_vs_lu_res_lfc2.csv")

br2_lib_vs_mus_de <- estimateGLMCommonDisp(br2_lib_vs_mus_dge_filtered_norm, biopanning_design)
br2_lib_vs_mus_de <- estimateGLMTrendedDisp(br2_lib_vs_mus_de, biopanning_design)
br2_lib_vs_mus_de <- estimateGLMTagwiseDisp(br2_lib_vs_mus_de, biopanning_design)
br2_lib_vs_mus_fit <- glmQLFit(br2_lib_vs_mus_de, biopanning_design)
br2_lib_vs_mus_test <-glmQLFTest(br2_lib_vs_mus_fit)
br2_lib_vs_mus_res <- topTags(br2_lib_vs_mus_test, sort.by = "logFC",p.value = 0.001, n = 200000)$table
br2_lib_vs_mus_res_lfc2 <- br2_lib_vs_mus_res %>% 
  filter(br2_lib_vs_mus_res$logFC > 2)
View(br2_lib_vs_mus_res_lfc2)
write.csv(br2_lib_vs_mus_res)



####DO not collapse ------ END of 2nd round BRAIN----

####round 2 LUNG
###----
#create datasets

lu2_lib_vs_li <- biopanning_all[,c(1,2,8,14,41,47,53)]
lu2_lib_vs_br <- biopanning_all[,c(1,2,8,14,39,45,51)]
lu2_lib_vs_ki <- biopanning_all[,c(1,2,8,14,42,48,54)]
lu2_lib_vs_lu <- biopanning_all[,c(1,2,8,14,40,46,52)]
lu2_lib_vs_mus <- biopanning_all[,c(1,2,8,14,43,49,55)]

View(lu2_lib_vs_lu)

#create groups
group_biopanning <- factor(c(1,1,1,2,2,2))

#create model matrix
biopanning_design <- model.matrix(~group_biopanning)

#make dge lists

lu2_lib_vs_li_dge <- DGEList(counts = lu2_lib_vs_li[,2:7], genes = lu2_lib_vs_li[,1], group = group_biopanning)
lu2_lib_vs_br_dge <- DGEList(counts = lu2_lib_vs_br[,2:7], genes = lu2_lib_vs_br[,1], group = group_biopanning)
lu2_lib_vs_ki_dge <- DGEList(counts = lu2_lib_vs_ki[,2:7], genes = lu2_lib_vs_ki[,1], group = group_biopanning)
lu2_lib_vs_lu_dge <- DGEList(counts = lu2_lib_vs_lu[,2:7], genes = lu2_lib_vs_lu[,1], group = group_biopanning)
lu2_lib_vs_mus_dge <- DGEList(counts = lu2_lib_vs_mus[,2:7], genes = lu2_lib_vs_mus[,1], group = group_biopanning)

nrow(lu2_lib_vs_li_dge$genes)

#filter out low abundance peptides

keep_lu2_lib_vs_li_dge <- rowSums((lu2_lib_vs_li_dge$counts)>5)>= 2
lu2_lib_vs_li_dge_filtered <- lu2_lib_vs_li_dge[keep_lu2_lib_vs_li_dge, , keep.lib.sizes = FALSE]

keep_lu2_lib_vs_br_dge <- rowSums((lu2_lib_vs_br_dge$counts)>5)>= 2
lu2_lib_vs_br_dge_filtered <- lu2_lib_vs_br_dge[keep_lu2_lib_vs_br_dge, , keep.lib.sizes = FALSE]

keep_lu2_lib_vs_ki_dge <- rowSums((lu2_lib_vs_ki_dge$counts)>5)>= 2
lu2_lib_vs_ki_dge_filtered <- lu2_lib_vs_ki_dge[keep_lu2_lib_vs_ki_dge, , keep.lib.sizes = FALSE]

keep_lu2_lib_vs_lu_dge <- rowSums((lu2_lib_vs_lu_dge$counts)>5)>= 2
lu2_lib_vs_lu_dge_filtered <- lu2_lib_vs_lu_dge[keep_lu2_lib_vs_lu_dge, , keep.lib.sizes = FALSE]

keep_lu2_lib_vs_mus_dge <- rowSums((lu2_lib_vs_mus_dge$counts)>5)>= 2
lu2_lib_vs_mus_dge_filtered <- lu2_lib_vs_mus_dge[keep_lu2_lib_vs_mus_dge, , keep.lib.sizes = FALSE]

nrow(lu2_lib_vs_br_dge_filtered)

plotMDS(lu2_lib_vs_lu_dge_filtered,top = 20)

#normalize samples
lu2_lib_vs_li_dge_filtered_norm <- calcNormFactors(lu2_lib_vs_li_dge_filtered)

lu2_lib_vs_br_dge_filtered_norm <- calcNormFactors(lu2_lib_vs_br_dge_filtered)

lu2_lib_vs_ki_dge_filtered_norm <- calcNormFactors(lu2_lib_vs_ki_dge_filtered)

lu2_lib_vs_lu_dge_filtered_norm <- calcNormFactors(lu2_lib_vs_lu_dge_filtered)

lu2_lib_vs_mus_dge_filtered_norm <- calcNormFactors(lu2_lib_vs_mus_dge_filtered)

#estimate dispersions, fit, perform stat tests and save results

###IT IS NECESSARY TO ADD logFC calculations


lu2_lib_vs_li_de <- estimateGLMCommonDisp(lu2_lib_vs_li_dge_filtered_norm, biopanning_design)
lu2_lib_vs_li_de <- estimateGLMTrendedDisp(lu2_lib_vs_li_de, biopanning_design)
lu2_lib_vs_li_de <- estimateGLMTagwiseDisp(lu2_lib_vs_li_de, biopanning_design)
lu2_lib_vs_li_fit <- glmQLFit(lu2_lib_vs_li_de, biopanning_design)
lu2_lib_vs_li_test <-glmQLFTest(lu2_lib_vs_li_fit)
lu2_lib_vs_li_res <- topTags(lu2_lib_vs_li_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
lu2_lib_vs_li_res_lfc2 <- lu2_lib_vs_li_res %>% 
  filter(lu2_lib_vs_li_res$logFC > 2)
View(lu2_lib_vs_li_res_lfc2)

lu2_lib_vs_br_de <- estimateGLMCommonDisp(lu2_lib_vs_br_dge_filtered_norm, biopanning_design)
lu2_lib_vs_br_de <- estimateGLMTrendedDisp(lu2_lib_vs_br_de, biopanning_design)
lu2_lib_vs_br_de <- estimateGLMTagwiseDisp(lu2_lib_vs_br_de, biopanning_design)
lu2_lib_vs_br_fit <- glmQLFit(lu2_lib_vs_br_de, biopanning_design)
lu2_lib_vs_br_test <-glmQLFTest(lu2_lib_vs_br_fit)
lu2_lib_vs_br_res <- topTags(lu2_lib_vs_br_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
lu2_lib_vs_br_res_lfc2 <- lu2_lib_vs_br_res %>% 
  filter(lu2_lib_vs_br_res$logFC > 2)
View(lu2_lib_vs_br_res_lfc2)
View(lu2_lib_vs_br_res)

lu2_lib_vs_ki_de <- estimateGLMCommonDisp(lu2_lib_vs_ki_dge_filtered_norm, biopanning_design)
lu2_lib_vs_ki_de <- estimateGLMTrendedDisp(lu2_lib_vs_ki_de, biopanning_design)
lu2_lib_vs_ki_de <- estimateGLMTagwiseDisp(lu2_lib_vs_ki_de, biopanning_design)
lu2_lib_vs_ki_fit <- glmQLFit(lu2_lib_vs_ki_de, biopanning_design)
lu2_lib_vs_ki_test <-glmQLFTest(lu2_lib_vs_ki_fit)
lu2_lib_vs_ki_res <- topTags(lu2_lib_vs_ki_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
lu2_lib_vs_ki_res_lfc2 <- lu2_lib_vs_ki_res %>% 
  filter(lu2_lib_vs_ki_res$logFC > 2)
View(lu2_lib_vs_ki_res_lfc2)

lu2_lib_vs_lu_de <- estimateGLMCommonDisp(lu2_lib_vs_lu_dge_filtered_norm, biopanning_design)
lu2_lib_vs_lu_de <- estimateGLMTrendedDisp(lu2_lib_vs_lu_de, biopanning_design)
lu2_lib_vs_lu_de <- estimateGLMTagwiseDisp(lu2_lib_vs_lu_de, biopanning_design)
lu2_lib_vs_lu_fit <- glmQLFit(lu2_lib_vs_lu_de, biopanning_design)
lu2_lib_vs_lu_test <-glmQLFTest(lu2_lib_vs_lu_fit)
lu2_lib_vs_lu_res <- topTags(lu2_lib_vs_lu_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
lu2_lib_vs_lu_res_lfc2 <- lu2_lib_vs_lu_res %>% 
  filter(lu2_lib_vs_lu_res$logFC > 2)
View(lu2_lib_vs_lu_res_lfc2)
View(lu2_lib_vs_lu_res)
write.csv(lu2_lib_vs_lu_res_lfc2, file = "Round2_lib_vs_lung_enrichment_lfc2.csv")

lu2_lib_vs_mus_de <- estimateGLMCommonDisp(lu2_lib_vs_mus_dge_filtered_norm, biopanning_design)
lu2_lib_vs_mus_de <- estimateGLMTrendedDisp(lu2_lib_vs_mus_de, biopanning_design)
lu2_lib_vs_mus_de <- estimateGLMTagwiseDisp(lu2_lib_vs_mus_de, biopanning_design)
lu2_lib_vs_mus_fit <- glmQLFit(lu2_lib_vs_mus_de, biopanning_design)
lu2_lib_vs_mus_test <-glmQLFTest(lu2_lib_vs_mus_fit)
lu2_lib_vs_mus_res <- topTags(lu2_lib_vs_mus_test, sort.by = "logFC",p.value = 0.001, n = 200000)$table
lu2_lib_vs_mus_res_lfc2 <- lu2_lib_vs_mus_res %>% 
  filter(lu2_lib_vs_mus_res$logFC > 2)
View(lu2_lib_vs_mus_res_lfc2)
write.csv(br2_lib_vs_mus_res)

## do not collapes THE END of Round 2 LUNG ----

### diagrams
## for ggplot use area1_2_3 and total_area - area1_2_3
###------

# input lib (1 of 18)
area_input1 <- nrow(subset(biopanning_all, biopanning_all$`442 Input Model Round1` != 0))
area_input2 <- nrow(subset(biopanning_all, biopanning_all$`443 Input Model Round1` != 0))
area_input3 <- nrow(subset(biopanning_all, biopanning_all$`444 Input Model Round1` != 0))
area_input1_2 <- nrow(subset(biopanning_all, biopanning_all$`442 Input Model Round1` != 0 & 
                               biopanning_all$`443 Input Model Round1` != 0))
area_input2_3 <- nrow(subset(biopanning_all, biopanning_all$`443 Input Model Round1` != 0 & 
                               biopanning_all$`444 Input Model Round1` != 0))
area_input1_3 <- nrow(subset(biopanning_all, biopanning_all$`442 Input Model Round1` != 0 & 
                               biopanning_all$`444 Input Model Round1` != 0))
area_input1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`442 Input Model Round1` != 0 & biopanning_all$`443 Input Model Round1` != 0 & biopanning_all$`444 Input Model Round1` != 0))

total_input <- (area_input1 + (area_input2 - area_input1_2) + (area_input3 - area_input1_3 - (area_input2_3 - area_input1_2_3)))

not_reprod_input <- total_input - area_input1_2_3

area_input_perc <- area_input1_2_3/total_input*100

total_input
area_input1_2_3 #this plus next should be total
not_reprod_input

area_input_perc

#brain 1 (2 of 18)

area_br1 <- nrow(subset(biopanning_all, biopanning_all$`442 Br Model Round1` != 0))
area_br2 <- nrow(subset(biopanning_all, biopanning_all$`443 Br Model Round1` != 0))
area_br3 <- nrow(subset(biopanning_all, biopanning_all$`444 Br Model Round1` != 0))
area_br1_2 <- nrow(subset(biopanning_all, biopanning_all$`442 Br Model Round1` != 0 & 
                               biopanning_all$`443 Br Model Round1` != 0))
area_br2_3 <- nrow(subset(biopanning_all, biopanning_all$`443 Br Model Round1` != 0 & 
                               biopanning_all$`444 Br Model Round1` != 0))
area_br1_3 <- nrow(subset(biopanning_all, biopanning_all$`442 Br Model Round1` != 0 & 
                               biopanning_all$`444 Br Model Round1` != 0))
area_br1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`442 Br Model Round1` != 0 & biopanning_all$`443 Br Model Round1` != 0 & biopanning_all$`444 Br Model Round1` != 0))

total_br1 <- (area_br1 + (area_br2 - area_br1_2) + (area_br3 - area_br1_3 - (area_br2_3 - area_br1_2_3)))

area_br1_perc <- area_br1_2_3/total_br1*100

not_reprod_br1 <- total_br1 - area_br1_2_3 

area_br1_2_3 #reproducible

not_reprod_br1 #not reproduc

total_br1 #total peptides

area_br_perc #reproducible perc



#lung 1 (3 of 18)

area_lu1 <- nrow(subset(biopanning_all, biopanning_all$`442 Lu Model Round1` != 0))
area_lu2 <- nrow(subset(biopanning_all, biopanning_all$`443 Lu Model Round1` != 0))
area_lu3 <- nrow(subset(biopanning_all, biopanning_all$`444 Lu Model Round1` != 0))
area_lu1_2 <- nrow(subset(biopanning_all, biopanning_all$`442 Lu Model Round1` != 0 & 
                            biopanning_all$`443 Lu Model Round1` != 0))
area_lu2_3 <- nrow(subset(biopanning_all, biopanning_all$`443 Lu Model Round1` != 0 & 
                            biopanning_all$`444 Lu Model Round1` != 0))
area_lu1_3 <- nrow(subset(biopanning_all, biopanning_all$`442 Lu Model Round1` != 0 & 
                            biopanning_all$`444 Lu Model Round1` != 0))
area_lu1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`442 Lu Model Round1` != 0 & biopanning_all$`443 Lu Model Round1` != 0 & biopanning_all$`444 Lu Model Round1` != 0))

total_lu1 <- (area_lu1 + (area_lu2 - area_lu1_2) + (area_lu3 - area_lu1_3 - (area_lu2_3 - area_lu1_2_3)))

area_lu1_perc <- area_lu1_2_3/total_lu1*100

not_reprod_lu1 <- total_lu1 - area_lu1_2_3

total_lu1 #total peptides

area_lu1_2_3 # reprod pep

not_reprod_lu1 # not reprod pep

area_lu1_perc # % reprod of all

#liver 1 (4 of 18)

area_li1 <- nrow(subset(biopanning_all, biopanning_all$`442 Li Model Round1` != 0))
area_li2 <- nrow(subset(biopanning_all, biopanning_all$`443 Li Model Round1` != 0))
area_li3 <- nrow(subset(biopanning_all, biopanning_all$`444 Li Model Round1` != 0))
area_li1_2 <- nrow(subset(biopanning_all, biopanning_all$`442 Li Model Round1` != 0 & 
                            biopanning_all$`443 Li Model Round1` != 0))
area_li2_3 <- nrow(subset(biopanning_all, biopanning_all$`443 Li Model Round1` != 0 & 
                            biopanning_all$`444 Li Model Round1` != 0))
area_li1_3 <- nrow(subset(biopanning_all, biopanning_all$`442 Li Model Round1` != 0 & 
                            biopanning_all$`444 Li Model Round1` != 0))
area_li1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`442 Li Model Round1` != 0 & biopanning_all$`443 Li Model Round1` != 0 & biopanning_all$`444 Li Model Round1` != 0))

total_li1 <- (area_li1 + (area_li2 - area_li1_2) + (area_li3 - area_li1_3 - (area_li2_3 - area_li1_2_3)))

area_li1_perc <- area_li1_2_3/total_li1*100

not_reprod_li1 <- total_li1 - area_li1_2_3

total_li1 #total peptides

area_li1_2_3 # reprod pep

not_reprod_li1 # not reprod pep

area_li1_perc # % reprod of all

#kidney 1 (5 of 18)

area_ki1 <- nrow(subset(biopanning_all, biopanning_all$`442 Ki Model Round1` != 0))
area_ki2 <- nrow(subset(biopanning_all, biopanning_all$`443 Ki Model Round1` != 0))
area_ki3 <- nrow(subset(biopanning_all, biopanning_all$`444 Ki Model Round1` != 0))
area_ki1_2 <- nrow(subset(biopanning_all, biopanning_all$`442 Ki Model Round1` != 0 & 
                            biopanning_all$`443 Ki Model Round1` != 0))
area_ki2_3 <- nrow(subset(biopanning_all, biopanning_all$`443 Ki Model Round1` != 0 & 
                            biopanning_all$`444 Ki Model Round1` != 0))
area_ki1_3 <- nrow(subset(biopanning_all, biopanning_all$`442 Ki Model Round1` != 0 & 
                            biopanning_all$`444 Ki Model Round1` != 0))
area_ki1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`442 Ki Model Round1` != 0 & biopanning_all$`443 Ki Model Round1` != 0 & biopanning_all$`444 Ki Model Round1` != 0))

total_ki1 <- (area_ki1 + (area_ki2 - area_ki1_2) + (area_ki3 - area_ki1_3 - (area_ki2_3 - area_ki1_2_3)))

area_ki1_perc <- area_ki1_2_3/total_ki1*100

not_reprod_ki1 <- total_ki1 - area_ki1_2_3

total_ki1 #total peptides

area_ki1_2_3 # reprod pep

not_reprod_ki1 # not reprod pep

area_ki1_perc # % reprod of all

#muscle 1 (6 of 18)

area_mus1 <- nrow(subset(biopanning_all, biopanning_all$`442 Mus Model Round1` != 0))
area_mus2 <- nrow(subset(biopanning_all, biopanning_all$`443 Mus Model Round1` != 0))
area_mus3 <- nrow(subset(biopanning_all, biopanning_all$`444 Mus Model Round1` != 0))
area_mus1_2 <- nrow(subset(biopanning_all, biopanning_all$`442 Mus Model Round1` != 0 & 
                            biopanning_all$`443 Mus Model Round1` != 0))
area_mus2_3 <- nrow(subset(biopanning_all, biopanning_all$`443 Mus Model Round1` != 0 & 
                            biopanning_all$`444 Mus Model Round1` != 0))
area_mus1_3 <- nrow(subset(biopanning_all, biopanning_all$`442 Ki Model Round1` != 0 & 
                            biopanning_all$`444 Mus Model Round1` != 0))
area_mus1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`442 Mus Model Round1` != 0 & biopanning_all$`443 Mus Model Round1` != 0 & biopanning_all$`444 Mus Model Round1` != 0))

total_mus1 <- (area_mus1 + (area_mus2 - area_mus1_2) + (area_mus3 - area_mus1_3 - (area_mus2_3 - area_mus1_2_3)))

area_mus1_perc <- area_mus1_2_3/total_mus1*100

not_reprod_mus1 <- total_mus1 - area_mus1_2_3

total_mus1 #total peptides

area_mus1_2_3 # reprod pep

not_reprod_mus1 # not reprod pep

area_mus1_perc # % reprod of all

#brain 2 input (7 of 18)

area_br2in_1 <- nrow(subset(biopanning_all, biopanning_all$`445 Brain Input for R2` != 0))
area_br2in_2 <- nrow(subset(biopanning_all, biopanning_all$`446 Brain input for R2` != 0))
area_br2in_3 <- nrow(subset(biopanning_all, biopanning_all$`447 Brain input for R2` != 0))
area_br2in_1_2 <- nrow(subset(biopanning_all, biopanning_all$`445 Brain Input for R2` != 0 & 
                             biopanning_all$`446 Brain input for R2` != 0))
area_br2in_2_3 <- nrow(subset(biopanning_all, biopanning_all$`446 Brain input for R2` != 0 & 
                             biopanning_all$`447 Brain input for R2` != 0))
area_br2in_1_3 <- nrow(subset(biopanning_all, biopanning_all$`445 Brain Input for R2` != 0 & 
                             biopanning_all$`447 Brain input for R2` != 0))
area_br2in_1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`445 Brain Input for R2` != 0 & biopanning_all$`446 Brain input for R2` != 0 & biopanning_all$`447 Brain input for R2` != 0))

total_br2in <- (area_br2in_1 + (area_br2in_2 - area_br2in_1_2) + (area_br2in_3 - area_br2in_1_3 - (area_br2in_2_3 - area_br2in_1_2_3)))

area_br2in_perc <- area_br2in_1_2_3/total_br2in*100

not_reprod_br2in <- total_br2in - area_br2in_1_2_3

total_br2in #total peptides

area_br2in_1_2_3 # reprod pep

not_reprod_br2in # not reprod pep

area_br2in_perc # % reprod of all

#brain2 brain (8 of 18)

area_br2br_1 <- nrow(subset(biopanning_all, biopanning_all$`445 Br Round2` != 0))
area_br2br_2 <- nrow(subset(biopanning_all, biopanning_all$`446 Br Round2` != 0))
area_br2br_3 <- nrow(subset(biopanning_all, biopanning_all$`447 Br Round2` != 0))
area_br2br_1_2 <- nrow(subset(biopanning_all, biopanning_all$`445 Br Round2` != 0 & 
                                biopanning_all$`446 Br Round2` != 0))
area_br2br_2_3 <- nrow(subset(biopanning_all, biopanning_all$`446 Br Round2` != 0 & 
                                biopanning_all$`447 Br Round2` != 0))
area_br2br_1_3 <- nrow(subset(biopanning_all, biopanning_all$`445 Br Round2` != 0 & 
                                biopanning_all$`447 Br Round2` != 0))
area_br2br_1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`445 Br Round2` != 0 & biopanning_all$`446 Br Round2` != 0 & biopanning_all$`447 Br Round2` != 0))

total_br2br <- (area_br2br_1 + (area_br2br_2 - area_br2br_1_2) + (area_br2br_3 - area_br2br_1_3 - (area_br2br_2_3 - area_br2br_1_2_3)))

area_br2br_perc <- area_br2br_1_2_3/total_br2br*100

not_reprod_br2br <- total_br2br - area_br2br_1_2_3

total_br2br #total peptides

area_br2br_1_2_3 # reprod pep

not_reprod_br2br # not reprod pep

area_br2br_perc # % reprod of all

#brain2 lung (9 of 18)

area_br2lu_1 <- nrow(subset(biopanning_all, biopanning_all$`445 Lu Round2` != 0))
area_br2lu_2 <- nrow(subset(biopanning_all, biopanning_all$`446 Lu Round2` != 0))
area_br2lu_3 <- nrow(subset(biopanning_all, biopanning_all$`447 Round2 Lu` != 0))
area_br2lu_1_2 <- nrow(subset(biopanning_all, biopanning_all$`445 Lu Round2` != 0 & 
                                biopanning_all$`446 Lu Round2` != 0))
area_br2lu_2_3 <- nrow(subset(biopanning_all, biopanning_all$`446 Lu Round2` != 0 & 
                                biopanning_all$`447 Round2 Lu` != 0))
area_br2lu_1_3 <- nrow(subset(biopanning_all, biopanning_all$`445 Lu Round2` != 0 & 
                                biopanning_all$`447 Round2 Lu` != 0))
area_br2lu_1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`445 Lu Round2` != 0 & biopanning_all$`446 Lu Round2` != 0 & biopanning_all$`447 Round2 Lu` != 0))

total_br2lu <- (area_br2lu_1 + (area_br2lu_2 - area_br2lu_1_2) + (area_br2lu_3 - area_br2lu_1_3 - (area_br2lu_2_3 - area_br2lu_1_2_3)))

area_br2lu_perc <- area_br2lu_1_2_3/total_br2lu*100

not_reprod_br2lu <- total_br2lu - area_br2lu_1_2_3

total_br2lu #total peptides

area_br2lu_1_2_3 # reprod pep

not_reprod_br2lu # not reprod pep

area_br2lu_perc # % reprod of all

#brain2 kidney (10 of 18)
area_br2ki_1 <- nrow(subset(biopanning_all, biopanning_all$`445 Ki Round2` != 0))
area_br2ki_2 <- nrow(subset(biopanning_all, biopanning_all$`446 Ki Round2` != 0))
area_br2ki_3 <- nrow(subset(biopanning_all, biopanning_all$`447 Round2 Ki` != 0))
area_br2ki_1_2 <- nrow(subset(biopanning_all, biopanning_all$`445 Ki Round2` != 0 & 
                                biopanning_all$`446 Ki Round2` != 0))
area_br2ki_2_3 <- nrow(subset(biopanning_all, biopanning_all$`446 Ki Round2` != 0 & 
                                biopanning_all$`447 Round2 Ki` != 0))
area_br2ki_1_3 <- nrow(subset(biopanning_all, biopanning_all$`445 Ki Round2` != 0 & 
                                biopanning_all$`447 Round2 Ki` != 0))
area_br2ki_1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`445 Ki Round2` != 0 & biopanning_all$`446 Ki Round2` != 0 & biopanning_all$`447 Round2 Ki` != 0))

total_br2ki <- (area_br2ki_1 + (area_br2ki_2 - area_br2ki_1_2) + (area_br2ki_3 - area_br2ki_1_3 - (area_br2ki_2_3 - area_br2ki_1_2_3)))

area_br2ki_perc <- area_br2ki_1_2_3/total_br2ki*100

not_reprod_br2ki <- total_br2ki - area_br2ki_1_2_3

total_br2ki #total peptides

area_br2ki_1_2_3 # reprod pep

not_reprod_br2ki # not reprod pep

area_br2ki_perc # % reprod of all

#brain2 liver (11 of 18)

area_br2li_1 <- nrow(subset(biopanning_all, biopanning_all$`445 Li Round2` != 0))
area_br2li_2 <- nrow(subset(biopanning_all, biopanning_all$`446 Li Round2` != 0))
area_br2li_3 <- nrow(subset(biopanning_all, biopanning_all$`447 Roun2 Li` != 0))
area_br2li_1_2 <- nrow(subset(biopanning_all, biopanning_all$`445 Li Round2` != 0 & 
                                biopanning_all$`446 Li Round2` != 0))
area_br2li_2_3 <- nrow(subset(biopanning_all, biopanning_all$`446 Li Round2` != 0 & 
                                biopanning_all$`447 Roun2 Li` != 0))
area_br2li_1_3 <- nrow(subset(biopanning_all, biopanning_all$`445 Li Round2` != 0 & 
                                biopanning_all$`447 Roun2 Li` != 0))
area_br2li_1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`445 Li Round2` != 0 & biopanning_all$`446 Li Round2` != 0 & biopanning_all$`447 Roun2 Li` != 0))

total_br2li <- (area_br2li_1 + (area_br2li_2 - area_br2li_1_2) + (area_br2li_3 - area_br2li_1_3 - (area_br2li_2_3 - area_br2li_1_2_3)))

area_br2li_perc <- area_br2li_1_2_3/total_br2li*100

not_reprod_br2li <- total_br2li - area_br2li_1_2_3

total_br2li #total peptides

area_br2li_1_2_3 # reprod pep

not_reprod_br2li # not reprod pep

area_br2li_perc # % reprod of all

#brain2 muscle (12 of 18)

area_br2mus_1 <- nrow(subset(biopanning_all, biopanning_all$`445 Mu Round2` != 0))
area_br2mus_2 <- nrow(subset(biopanning_all, biopanning_all$`446 Mu Round2` != 0))
area_br2mus_3 <- nrow(subset(biopanning_all, biopanning_all$`447 Round2 Mu` != 0))
area_br2mus_1_2 <- nrow(subset(biopanning_all, biopanning_all$`445 Mu Round2` != 0 & 
                                biopanning_all$`446 Mu Round2` != 0))
area_br2mus_2_3 <- nrow(subset(biopanning_all, biopanning_all$`446 Mu Round2` != 0 & 
                                biopanning_all$`447 Round2 Mu` != 0))
area_br2mus_1_3 <- nrow(subset(biopanning_all, biopanning_all$`445 Mu Round2` != 0 & 
                                biopanning_all$`447 Round2 Mu` != 0))
area_br2mus_1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`445 Mu Round2` != 0 & biopanning_all$`446 Mu Round2` != 0 & biopanning_all$`447 Round2 Mu` != 0))

total_br2mus <- (area_br2mus_1 + (area_br2mus_2 - area_br2mus_1_2) + (area_br2mus_3 - area_br2mus_1_3 - (area_br2mus_2_3 - area_br2mus_1_2_3)))

area_br2mus_perc <- area_br2mus_1_2_3/total_br2mus*100

not_reprod_br2mus <- total_br2mus - area_br2mus_1_2_3

total_br2mus #total peptides

area_br2mus_1_2_3 # reprod pep

not_reprod_br2mus # not reprod pep

area_br2mus_perc # % reprod of all

#lung2 input (13 of 18)

area_lu2in_1 <- nrow(subset(biopanning_all, biopanning_all$`448 Lung Input for R2` != 0))
area_lu2in_2 <- nrow(subset(biopanning_all, biopanning_all$`449 Lung Input for R2` != 0))
area_lu2in_3 <- nrow(subset(biopanning_all, biopanning_all$`450 Lung Input for R2` != 0))
area_lu2in_1_2 <- nrow(subset(biopanning_all, biopanning_all$`448 Lung Input for R2` != 0 & 
                                 biopanning_all$`449 Lung Input for R2` != 0))
area_lu2in_2_3 <- nrow(subset(biopanning_all, biopanning_all$`449 Lung Input for R2` != 0 & 
                                 biopanning_all$`450 Lung Input for R2` != 0))
area_lu2in_1_3 <- nrow(subset(biopanning_all, biopanning_all$`448 Lung Input for R2` != 0 & 
                                 biopanning_all$`450 Lung Input for R2` != 0))
area_lu2in_1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`448 Lung Input for R2` != 0 & biopanning_all$`449 Lung Input for R2` != 0 & biopanning_all$`450 Lung Input for R2` != 0))

total_lu2in <- (area_lu2in_1 + (area_lu2in_2 - area_lu2in_1_2) + (area_lu2in_3 - area_lu2in_1_3 - (area_lu2in_2_3 - area_lu2in_1_2_3)))

area_lu2in_perc <- area_lu2in_1_2_3/total_lu2in*100

not_reprod_lu2in <- total_lu2in - area_lu2in_1_2_3

total_lu2in #total peptides

area_lu2in_1_2_3 # reprod pep

not_reprod_lu2in # not reprod pep

area_lu2in_perc # % reprod of all


#lung2 lung (14 of 18)

area_lu2lu_1 <- nrow(subset(biopanning_all, biopanning_all$`448 Lu` != 0))
area_lu2lu_2 <- nrow(subset(biopanning_all, biopanning_all$`449 Lu` != 0))
area_lu2lu_3 <- nrow(subset(biopanning_all, biopanning_all$`450 Lu` != 0))
area_lu2lu_1_2 <- nrow(subset(biopanning_all, biopanning_all$`448 Lu` != 0 & 
                                biopanning_all$`449 Lu` != 0))
area_lu2lu_2_3 <- nrow(subset(biopanning_all, biopanning_all$`449 Lu` != 0 & 
                                biopanning_all$`450 Lu` != 0))
area_lu2lu_1_3 <- nrow(subset(biopanning_all, biopanning_all$`448 Lu` != 0 & 
                                biopanning_all$`450 Lu` != 0))
area_lu2lu_1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`448 Lu` != 0 & biopanning_all$`449 Lu` != 0 & biopanning_all$`450 Lu` != 0))

total_lu2lu <- (area_lu2lu_1 + (area_lu2lu_2 - area_lu2lu_1_2) + (area_lu2lu_3 - area_lu2lu_1_3 - (area_lu2lu_2_3 - area_lu2lu_1_2_3)))

area_lu2lu_perc <- area_lu2lu_1_2_3/total_lu2lu*100

not_reprod_lu2lu <- total_lu2lu - area_lu2lu_1_2_3

total_lu2lu #total peptides

area_lu2lu_1_2_3 # reprod pep

not_reprod_lu2lu # not reprod pep

area_lu2lu_perc # % reprod of all

#lung2 brain (15 of 18)

area_lu2br_1 <- nrow(subset(biopanning_all, biopanning_all$`448 Br` != 0))
area_lu2br_2 <- nrow(subset(biopanning_all, biopanning_all$`449 Br` != 0))
area_lu2br_3 <- nrow(subset(biopanning_all, biopanning_all$`450 Br` != 0))
area_lu2br_1_2 <- nrow(subset(biopanning_all, biopanning_all$`448 Br` != 0 & 
                                biopanning_all$`449 Br` != 0))
area_lu2br_2_3 <- nrow(subset(biopanning_all, biopanning_all$`449 Br` != 0 & 
                                biopanning_all$`450 Br` != 0))
area_lu2br_1_3 <- nrow(subset(biopanning_all, biopanning_all$`448 Br` != 0 & 
                                biopanning_all$`450 Br` != 0))
area_lu2br_1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`448 Br` != 0 & biopanning_all$`449 Br` != 0 & biopanning_all$`450 Br` != 0))

total_lu2br <- (area_lu2br_1 + (area_lu2br_2 - area_lu2br_1_2) + (area_lu2br_3 - area_lu2br_1_3 - (area_lu2br_2_3 - area_lu2br_1_2_3)))

area_lu2br_perc <- area_lu2br_1_2_3/total_lu2br*100

not_reprod_lu2br <- total_lu2br - area_lu2br_1_2_3

total_lu2br #total peptides

area_lu2br_1_2_3 # reprod pep

not_reprod_lu2br # not reprod pep

area_lu2br_perc # % reprod of all

#lung2 kidney(16 of 18)

area_lu2ki_1 <- nrow(subset(biopanning_all, biopanning_all$`448 Ki` != 0))
area_lu2ki_2 <- nrow(subset(biopanning_all, biopanning_all$`449 Ki` != 0))
area_lu2ki_3 <- nrow(subset(biopanning_all, biopanning_all$`450 Ki` != 0))
area_lu2ki_1_2 <- nrow(subset(biopanning_all, biopanning_all$`448 Ki` != 0 & 
                                biopanning_all$`449 Ki` != 0))
area_lu2ki_2_3 <- nrow(subset(biopanning_all, biopanning_all$`449 Ki` != 0 & 
                                biopanning_all$`450 Ki` != 0))
area_lu2ki_1_3 <- nrow(subset(biopanning_all, biopanning_all$`448 Ki` != 0 & 
                                biopanning_all$`450 Ki` != 0))
area_lu2ki_1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`448 Ki` != 0 & biopanning_all$`449 Ki` != 0 & biopanning_all$`450 Ki` != 0))

total_lu2ki <- (area_lu2ki_1 + (area_lu2ki_2 - area_lu2ki_1_2) + (area_lu2ki_3 - area_lu2ki_1_3 - (area_lu2ki_2_3 - area_lu2ki_1_2_3)))

area_lu2ki_perc <- area_lu2ki_1_2_3/total_lu2ki*100

not_reprod_lu2ki <- total_lu2ki - area_lu2ki_1_2_3

total_lu2ki #total peptides

area_lu2ki_1_2_3 # reprod pep

not_reprod_lu2ki # not reprod pep

area_lu2ki_perc # % reprod of all


#lung2 liver (17 of 18)

area_lu2li_1 <- nrow(subset(biopanning_all, biopanning_all$`448 Li` != 0))
area_lu2li_2 <- nrow(subset(biopanning_all, biopanning_all$`449 Li` != 0))
area_lu2li_3 <- nrow(subset(biopanning_all, biopanning_all$`450 Li` != 0))
area_lu2li_1_2 <- nrow(subset(biopanning_all, biopanning_all$`448 Li` != 0 & 
                                biopanning_all$`449 Li` != 0))
area_lu2li_2_3 <- nrow(subset(biopanning_all, biopanning_all$`449 Li` != 0 & 
                                biopanning_all$`450 Li` != 0))
area_lu2li_1_3 <- nrow(subset(biopanning_all, biopanning_all$`448 Li` != 0 & 
                                biopanning_all$`450 Li` != 0))
area_lu2li_1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`448 Li` != 0 & biopanning_all$`449 Li` != 0 & biopanning_all$`450 Li` != 0))

total_lu2li <- (area_lu2li_1 + (area_lu2li_2 - area_lu2li_1_2) + (area_lu2li_3 - area_lu2li_1_3 - (area_lu2li_2_3 - area_lu2li_1_2_3)))

area_lu2li_perc <- area_lu2li_1_2_3/total_lu2li*100

not_reprod_lu2li <- total_lu2li - area_lu2li_1_2_3

total_lu2li #total peptides

area_lu2li_1_2_3 # reprod pep

not_reprod_lu2li # not reprod pep

area_lu2li_perc # % reprod of all



#lung2 muscle (18 of 18)

area_lu2mus_1 <- nrow(subset(biopanning_all, biopanning_all$`448 Mus` != 0))
area_lu2mus_2 <- nrow(subset(biopanning_all, biopanning_all$`449 Mus` != 0))
area_lu2mus_3 <- nrow(subset(biopanning_all, biopanning_all$`450 Mus` != 0))
area_lu2mus_1_2 <- nrow(subset(biopanning_all, biopanning_all$`448 Mus` != 0 & 
                                biopanning_all$`449 Mus` != 0))
area_lu2mus_2_3 <- nrow(subset(biopanning_all, biopanning_all$`449 Mus` != 0 & 
                                biopanning_all$`450 Mus` != 0))
area_lu2mus_1_3 <- nrow(subset(biopanning_all, biopanning_all$`448 Mus` != 0 & 
                                biopanning_all$`450 Mus` != 0))
area_lu2mus_1_2_3 <- nrow(subset(biopanning_all, biopanning_all$`448 Mus` != 0 & biopanning_all$`449 Mus` != 0 & biopanning_all$`450 Mus` != 0))

total_lu2mus <- (area_lu2mus_1 + (area_lu2mus_2 - area_lu2mus_1_2) + (area_lu2mus_3 - area_lu2mus_1_3 - (area_lu2mus_2_3 - area_lu2mus_1_2_3)))

area_lu2mus_perc <- area_lu2mus_1_2_3/total_lu2mus*100

not_reprod_lu2mus <- total_lu2mus - area_lu2mus_1_2_3

total_lu2mus #total peptides

area_lu2mus_1_2_3 # reprod pep

not_reprod_lu2mus # not reprod pep

area_lu2mus_perc # % reprod of all

###create data frame with all info
# selection - reliable/unrealiable - count

# rounds

found_peptides <- data.frame("In","Round_1", area_input1_2_3, not_reprod_input)

found_peptides

colnames(found_peptides) <- c("selection", "round", "rel_no_reads", "unrel_no_reads")

found_peptides <- add_row(found_peptides, selection = "Br", 
                                round = "Round_1", rel_no_reads = area_br1_2_3 , unrel_no_reads = not_reprod_br1)
found_peptides <- add_row(found_peptides, selection = "Lu", 
                          round = "Round_1", rel_no_reads = area_lu1_2_3 , unrel_no_reads = not_reprod_lu1)
found_peptides <- add_row(found_peptides, selection = "Li", 
                          round = "Round_1", rel_no_reads = area_li1_2_3 , unrel_no_reads = not_reprod_li1)
found_peptides <- add_row(found_peptides, selection = "Ki", 
                          round = "Round_1", rel_no_reads = area_ki1_2_3 , unrel_no_reads = not_reprod_ki1)
found_peptides <- add_row(found_peptides, selection = "Mu", 
                          round = "Round_1", rel_no_reads = area_mus1_2_3 , unrel_no_reads = not_reprod_mus1)

found_peptides <- add_row(found_peptides, selection = "In", 
                          round = "Round_2_BR", rel_no_reads = area_br2in_1_2_3 , unrel_no_reads = not_reprod_br2in)
found_peptides <- add_row(found_peptides, selection = "Br", 
                          round = "Round_2_BR", rel_no_reads = area_br2br_1_2_3 , unrel_no_reads = not_reprod_br2br)
found_peptides <- add_row(found_peptides, selection = "Lu", 
                          round = "Round_2_BR", rel_no_reads = area_br2lu_1_2_3 , unrel_no_reads = not_reprod_br2lu)
found_peptides <- add_row(found_peptides, selection = "Li", 
                          round = "Round_2_BR", rel_no_reads = area_br2li_1_2_3, unrel_no_reads = not_reprod_br2li)
found_peptides <- add_row(found_peptides, selection = "Ki", 
                          round = "Round_2_BR", rel_no_reads = area_br2ki_1_2_3 , unrel_no_reads = not_reprod_br2ki)
found_peptides <- add_row(found_peptides, selection = "Mu", 
                          round = "Round_2_BR", rel_no_reads = area_br2mus_1_2_3 , unrel_no_reads = not_reprod_br2mus)

found_peptides <- add_row(found_peptides, selection = "In", 
                          round = "Round_2_LU", rel_no_reads = area_lu2in_1_2_3 , unrel_no_reads = not_reprod_lu2in)
found_peptides <- add_row(found_peptides, selection = "Br", 
                          round = "Round_2_LU", rel_no_reads = area_lu2br_1_2_3 , unrel_no_reads = not_reprod_lu2br)
found_peptides <- add_row(found_peptides, selection = "Lu", 
                          round = "Round_2_LU", rel_no_reads = area_lu2lu_1_2_3 , unrel_no_reads = not_reprod_lu2lu)
found_peptides <- add_row(found_peptides, selection = "Li", 
                          round = "Round_2_LU", rel_no_reads = area_lu2li_1_2_3, unrel_no_reads = not_reprod_lu2li)
found_peptides <- add_row(found_peptides, selection = "Ki", 
                          round = "Round_2_LU", rel_no_reads = area_lu2ki_1_2_3 , unrel_no_reads = not_reprod_lu2ki)
found_peptides <- add_row(found_peptides, selection = "Mu", 
                          round = "Round_2_LU", rel_no_reads = area_lu2mus_1_2_3 , unrel_no_reads = not_reprod_lu2mus)


found_peptides$reliable_perc_of_all <- found_peptides$rel_no_reads/(found_peptides$rel_no_reads + found_peptides$unrel_no_reads) * 100

View(found_peptides)
#all on one graph
#plot_r2 <- ggplot(found_peptides, aes(selection, reliable_perc_of_all, fill = round, colour = round))

found_peptides

plot_r2 <- ggplot(found_peptides, aes(selection, reliable_perc_of_all, fill = round))
plot_r2 + geom_bar(stat = "identity", position = "dodge")+ theme_linedraw()# + facet_grid(rows = vars(round)) + theme_bw()

plot_r2 + geom_bar(stat = "identity", position = "stack")

plot_r1 <- found_peptides %>% 
  filter(found_peptides$round == "Round_1")
View(plot_r2_lu)

plot_r2_br <- found_peptides %>% 
  filter(found_peptides$round == "Round_2_BR")

plot_r2_lu <- found_peptides %>% 
  filter(found_peptides$round == "Round_2_LU")

level_order <- c("In", "Br", "Lu", "Ki", "Li", "Mu")

plot_r1gg <- ggplot(plot_r1, aes(x = factor(selection, level = level_order), reliable_perc_of_all, level = level_order))
round_1_plot <- plot_r1gg + geom_bar(stat = "identity", position = "dodge")+ theme_linedraw() + ylim(0,35) + ylab("% of reliable peptides") + xlab("Round 1")

plot_r2_brgg <- ggplot(plot_r2_br, aes(x = factor(selection, level = level_order), reliable_perc_of_all, level=level_order))
round_2_br_plot <- plot_r2_brgg + geom_bar(stat = "identity", position = "dodge")+ theme_linedraw() + ylim(0,35) + ylab("% of reliable peptides") + xlab("Round 2 Brain")

plot_r2_lugg <- ggplot(plot_r2_lu, aes(x = factor(selection, level = level_order), reliable_perc_of_all, level=level_order))
round_2_lu_plot <- plot_r2_lugg + geom_bar(stat = "identity", position = "dodge")+ theme_linedraw() + ylim(0,35) + ylab("% of reliable peptides") + xlab("Round 2 Lung")

figure <- ggarrange(round_1_plot, round_2_br_plot, round_2_lu_plot, 
                    labels = c("a", "b", "c"), ncol = 3, nrow = 1)

figure

ggsave("reliable-pep.pdf", width = 8, height = 2.5)

# DO NOT COLLAPSE - END OF 

## organ specificity Round 1



##Organ specificity Round 2

#create datasets

lib_vs_li <- biopanning_r1[,c(1,2,8,14,5,11,17)]
lib_vs_br <- biopanning_r1[,c(1,2,8,14,3,9,15)]
lib_vs_ki <- biopanning_r1[,c(1,2,8,14,6,12,18)]
lib_vs_lu <- biopanning_r1[,c(1,2,8,14,4,10,16)]
lib_vs_mus <- biopanning_r1[,c(1,2,8,14,7,13,19)]

View(lib_vs_mus)

#create groups
group_biopanning <- factor(c(1,1,1,2,2,2))

#create model matrix
biopanning_design <- model.matrix(~group_biopanning)

#make dge lists

lib_vs_li_dge <- DGEList(counts = lib_vs_li[,2:7], genes = lib_vs_li[,1], group = group_biopanning)
lib_vs_br_dge <- DGEList(counts = lib_vs_br[,2:7], genes = lib_vs_br[,1], group = group_biopanning)
lib_vs_ki_dge <- DGEList(counts = lib_vs_ki[,2:7], genes = lib_vs_ki[,1], group = group_biopanning)
lib_vs_lu_dge <- DGEList(counts = lib_vs_lu[,2:7], genes = lib_vs_lu[,1], group = group_biopanning)
lib_vs_mus_dge <- DGEList(counts = lib_vs_mus[,2:7], genes = lib_vs_mus[,1], group = group_biopanning)

nrow(lib_vs_li_dge$genes)

#filter out low abundance peptides

keep_lib_vs_li_dge <- rowSums((lib_vs_li_dge$counts)>5)>= 2
lib_vs_li_dge_filtered <- lib_vs_li_dge[keep_lib_vs_li_dge, , keep.lib.sizes = FALSE]

keep_lib_vs_br_dge <- rowSums((lib_vs_br_dge$counts)>5)>= 2
lib_vs_br_dge_filtered <- lib_vs_br_dge[keep_lib_vs_br_dge, , keep.lib.sizes = FALSE]

keep_lib_vs_ki_dge <- rowSums((lib_vs_ki_dge$counts)>5)>= 2
lib_vs_ki_dge_filtered <- lib_vs_ki_dge[keep_lib_vs_ki_dge, , keep.lib.sizes = FALSE]

keep_lib_vs_lu_dge <- rowSums((lib_vs_lu_dge$counts)>5)>= 2
lib_vs_lu_dge_filtered <- lib_vs_lu_dge[keep_lib_vs_lu_dge, , keep.lib.sizes = FALSE]

keep_lib_vs_mus_dge <- rowSums((lib_vs_mus_dge$counts)>5)>= 2
lib_vs_mus_dge_filtered <- lib_vs_mus_dge[keep_lib_vs_mus_dge, , keep.lib.sizes = FALSE]

nrow(lib_vs_br_dge_filtered)

plotMDS(lib_vs_br_dge_filtered,top = 20)

#normalize samples
lib_vs_li_dge_filtered_norm <- calcNormFactors(lib_vs_li_dge_filtered)

lib_vs_br_dge_filtered_norm <- calcNormFactors(lib_vs_br_dge_filtered)

lib_vs_ki_dge_filtered_norm <- calcNormFactors(lib_vs_ki_dge_filtered)

lib_vs_lu_dge_filtered_norm <- calcNormFactors(lib_vs_lu_dge_filtered)

lib_vs_mus_dge_filtered_norm <- calcNormFactors(lib_vs_mus_dge_filtered)

#estimate dispersions, fit, perform stat tests and save results

###IT IS NECESSARY TO ADD logFC calculations


lib_vs_li_de <- estimateGLMCommonDisp(lib_vs_li_dge_filtered_norm, biopanning_design)
lib_vs_li_de <- estimateGLMTrendedDisp(lib_vs_li_de, biopanning_design)
lib_vs_li_de <- estimateGLMTagwiseDisp(lib_vs_li_de, biopanning_design)
lib_vs_li_fit <- glmQLFit(lib_vs_li_de, biopanning_design)
lib_vs_li_test <-glmQLFTest(lib_vs_li_fit)
lib_vs_li_res <- topTags(lib_vs_li_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
lib_vs_li_res_lfc2 <- lib_vs_li_res %>% 
  filter(lib_vs_li_res$logFC > 2)
View(lib_vs_li_res_lfc2)

lib_vs_br_de <- estimateGLMCommonDisp(lib_vs_br_dge_filtered_norm, biopanning_design)
lib_vs_br_de <- estimateGLMTrendedDisp(lib_vs_br_de, biopanning_design)
lib_vs_br_de <- estimateGLMTagwiseDisp(lib_vs_br_de, biopanning_design)
lib_vs_br_fit <- glmQLFit(lib_vs_br_de, biopanning_design)
lib_vs_br_test <-glmQLFTest(lib_vs_br_fit)
lib_vs_br_res <- topTags(lib_vs_br_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
lib_vs_br_res_lfc2 <- lib_vs_br_res %>% 
  filter(lib_vs_br_res$logFC > 2)
View(lib_vs_br_res_lfc2)
View(lib_vs_br_res)

lib_vs_ki_de <- estimateGLMCommonDisp(lib_vs_ki_dge_filtered_norm, biopanning_design)
lib_vs_ki_de <- estimateGLMTrendedDisp(lib_vs_ki_de, biopanning_design)
lib_vs_ki_de <- estimateGLMTagwiseDisp(lib_vs_ki_de, biopanning_design)
lib_vs_ki_fit <- glmQLFit(lib_vs_ki_de, biopanning_design)
lib_vs_ki_test <-glmQLFTest(lib_vs_ki_fit)
lib_vs_ki_res <- topTags(lib_vs_ki_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
lib_vs_ki_res_lfc2 <- lib_vs_ki_res %>% 
  filter(lib_vs_ki_res$logFC > 2)
View(lib_vs_ki_res_lfc2)

lib_vs_lu_de <- estimateGLMCommonDisp(lib_vs_lu_dge_filtered_norm, biopanning_design)
lib_vs_lu_de <- estimateGLMTrendedDisp(lib_vs_lu_de, biopanning_design)
lib_vs_lu_de <- estimateGLMTagwiseDisp(lib_vs_lu_de, biopanning_design)
lib_vs_lu_fit <- glmQLFit(lib_vs_lu_de, biopanning_design)
lib_vs_lu_test <-glmQLFTest(lib_vs_lu_fit)
lib_vs_lu_res <- topTags(lib_vs_lu_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
lib_vs_lu_res_lfc2 <- lib_vs_lu_res %>% 
  filter(lib_vs_lu_res$logFC > 2)
View(lib_vs_lu_res_lfc2)
View(lib_vs_lu_res)

lib_vs_mus_de <- estimateGLMCommonDisp(lib_vs_mus_dge_filtered_norm, biopanning_design)
lib_vs_mus_de <- estimateGLMTrendedDisp(lib_vs_mus_de, biopanning_design)
lib_vs_mus_de <- estimateGLMTagwiseDisp(lib_vs_mus_de, biopanning_design)
lib_vs_mus_fit <- glmQLFit(lib_vs_mus_de, biopanning_design)
lib_vs_mus_test <-glmQLFTest(lib_vs_mus_fit)
lib_vs_mus_res <- topTags(lib_vs_mus_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
lib_vs_mus_res_lfc2 <- lib_vs_mus_res %>% 
  filter(lib_vs_mus_res$logFC > 2)
View(lib_vs_mus_res_lfc2)
write.csv(lib_vs_mus_res)


###DO not collapse----END of 1st round----------


###round 2 analysis

### round 2 brain
###----

###round 2 brain

#create datasets

br2_lib_vs_li <- biopanning_all[,c(1,2,8,14,23,29,35)]
br2_lib_vs_br <- biopanning_all[,c(1,2,8,14,21,27,33)]
br2_lib_vs_ki <- biopanning_all[,c(1,2,8,14,24,30,36)]
br2_lib_vs_lu <- biopanning_all[,c(1,2,8,14,22,28,34)]
br2_lib_vs_mus <- biopanning_all[,c(1,2,8,14,25,31,37)]

View(br2_lib_vs_ki)

#create groups
group_biopanning <- factor(c(1,1,1,2,2,2))

#create model matrix
biopanning_design <- model.matrix(~group_biopanning)

#make dge lists

br2_lib_vs_li_dge <- DGEList(counts = br2_lib_vs_li[,2:7], genes = br2_lib_vs_li[,1], group = group_biopanning)
br2_lib_vs_br_dge <- DGEList(counts = br2_lib_vs_br[,2:7], genes = br2_lib_vs_br[,1], group = group_biopanning)
br2_lib_vs_ki_dge <- DGEList(counts = br2_lib_vs_ki[,2:7], genes = br2_lib_vs_ki[,1], group = group_biopanning)
br2_lib_vs_lu_dge <- DGEList(counts = br2_lib_vs_lu[,2:7], genes = br2_lib_vs_lu[,1], group = group_biopanning)
br2_lib_vs_mus_dge <- DGEList(counts = br2_lib_vs_mus[,2:7], genes = br2_lib_vs_mus[,1], group = group_biopanning)

nrow(br2_lib_vs_li_dge$genes)

#filter out low abundance peptides

keep_br2_lib_vs_li_dge <- rowSums((br2_lib_vs_li_dge$counts)>5)>= 2
br2_lib_vs_li_dge_filtered <- br2_lib_vs_li_dge[keep_br2_lib_vs_li_dge, , keep.lib.sizes = FALSE]

keep_br2_lib_vs_br_dge <- rowSums((br2_lib_vs_br_dge$counts)>5)>= 2
br2_lib_vs_br_dge_filtered <- br2_lib_vs_br_dge[keep_br2_lib_vs_br_dge, , keep.lib.sizes = FALSE]

keep_br2_lib_vs_ki_dge <- rowSums((br2_lib_vs_ki_dge$counts)>5)>= 2
br2_lib_vs_ki_dge_filtered <- br2_lib_vs_ki_dge[keep_br2_lib_vs_ki_dge, , keep.lib.sizes = FALSE]

keep_br2_lib_vs_lu_dge <- rowSums((br2_lib_vs_lu_dge$counts)>5)>= 2
br2_lib_vs_lu_dge_filtered <- br2_lib_vs_lu_dge[keep_br2_lib_vs_lu_dge, , keep.lib.sizes = FALSE]

keep_br2_lib_vs_mus_dge <- rowSums((br2_lib_vs_mus_dge$counts)>5)>= 2
br2_lib_vs_mus_dge_filtered <- br2_lib_vs_mus_dge[keep_br2_lib_vs_mus_dge, , keep.lib.sizes = FALSE]

nrow(br2_lib_vs_br_dge_filtered)

plotMDS(br2_lib_vs_br_dge_filtered,top = 20)

#normalize samples
br2_lib_vs_li_dge_filtered_norm <- calcNormFactors(br2_lib_vs_li_dge_filtered)

br2_lib_vs_br_dge_filtered_norm <- calcNormFactors(br2_lib_vs_br_dge_filtered)

br2_lib_vs_ki_dge_filtered_norm <- calcNormFactors(br2_lib_vs_ki_dge_filtered)

br2_lib_vs_lu_dge_filtered_norm <- calcNormFactors(br2_lib_vs_lu_dge_filtered)

br2_lib_vs_mus_dge_filtered_norm <- calcNormFactors(br2_lib_vs_mus_dge_filtered)

#estimate dispersions, fit, perform stat tests and save results

###IT IS NECESSARY TO ADD logFC calculations


br2_lib_vs_li_de <- estimateGLMCommonDisp(br2_lib_vs_li_dge_filtered_norm, biopanning_design)
br2_lib_vs_li_de <- estimateGLMTrendedDisp(br2_lib_vs_li_de, biopanning_design)
br2_lib_vs_li_de <- estimateGLMTagwiseDisp(br2_lib_vs_li_de, biopanning_design)
br2_lib_vs_li_fit <- glmQLFit(br2_lib_vs_li_de, biopanning_design)
br2_lib_vs_li_test <-glmQLFTest(br2_lib_vs_li_fit)
br2_lib_vs_li_res <- topTags(br2_lib_vs_li_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
br2_lib_vs_li_res_lfc2 <- br2_lib_vs_li_res %>% 
  filter(br2_lib_vs_li_res$logFC > 2)
View(br2_lib_vs_li_res_lfc2)

br2_lib_vs_br_de <- estimateGLMCommonDisp(br2_lib_vs_br_dge_filtered_norm, biopanning_design)
br2_lib_vs_br_de <- estimateGLMTrendedDisp(br2_lib_vs_br_de, biopanning_design)
br2_lib_vs_br_de <- estimateGLMTagwiseDisp(br2_lib_vs_br_de, biopanning_design)
br2_lib_vs_br_fit <- glmQLFit(br2_lib_vs_br_de, biopanning_design)
br2_lib_vs_br_test <-glmQLFTest(br2_lib_vs_br_fit)
br2_lib_vs_br_res <- topTags(br2_lib_vs_br_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
br2_lib_vs_br_res_lfc2 <- br2_lib_vs_br_res %>% 
  filter(br2_lib_vs_br_res$logFC > 2)
View(br2_lib_vs_br_res_lfc2)
View(br2_lib_vs_br_res)
write.csv(br2_lib_vs_br_res_lfc2, file = "br2_lib_vs_br_res_lfc2.csv")

br2_lib_vs_ki_de <- estimateGLMCommonDisp(br2_lib_vs_ki_dge_filtered_norm, biopanning_design)
br2_lib_vs_ki_de <- estimateGLMTrendedDisp(br2_lib_vs_ki_de, biopanning_design)
br2_lib_vs_ki_de <- estimateGLMTagwiseDisp(br2_lib_vs_ki_de, biopanning_design)
br2_lib_vs_ki_fit <- glmQLFit(br2_lib_vs_ki_de, biopanning_design)
br2_lib_vs_ki_test <-glmQLFTest(br2_lib_vs_ki_fit)
br2_lib_vs_ki_res <- topTags(br2_lib_vs_ki_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
br2_lib_vs_ki_res_lfc2 <- br2_lib_vs_ki_res %>% 
  filter(br2_lib_vs_ki_res$logFC > 2)
View(br2_lib_vs_ki_res_lfc2)

br2_lib_vs_lu_de <- estimateGLMCommonDisp(br2_lib_vs_lu_dge_filtered_norm, biopanning_design)
br2_lib_vs_lu_de <- estimateGLMTrendedDisp(br2_lib_vs_lu_de, biopanning_design)
br2_lib_vs_lu_de <- estimateGLMTagwiseDisp(br2_lib_vs_lu_de, biopanning_design)
br2_lib_vs_lu_fit <- glmQLFit(br2_lib_vs_lu_de, biopanning_design)
br2_lib_vs_lu_test <-glmQLFTest(br2_lib_vs_lu_fit)
br2_lib_vs_lu_res <- topTags(br2_lib_vs_lu_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
br2_lib_vs_lu_res_lfc2 <- br2_lib_vs_lu_res %>% 
  filter(br2_lib_vs_lu_res$logFC > 2)
View(br2_lib_vs_lu_res_lfc2)
View(br2_lib_vs_lu_res)
write.csv(br2_lib_vs_lu_res_lfc2, file = "br2_lib_vs_lu_res_lfc2.csv")

br2_lib_vs_mus_de <- estimateGLMCommonDisp(br2_lib_vs_mus_dge_filtered_norm, biopanning_design)
br2_lib_vs_mus_de <- estimateGLMTrendedDisp(br2_lib_vs_mus_de, biopanning_design)
br2_lib_vs_mus_de <- estimateGLMTagwiseDisp(br2_lib_vs_mus_de, biopanning_design)
br2_lib_vs_mus_fit <- glmQLFit(br2_lib_vs_mus_de, biopanning_design)
br2_lib_vs_mus_test <-glmQLFTest(br2_lib_vs_mus_fit)
br2_lib_vs_mus_res <- topTags(br2_lib_vs_mus_test, sort.by = "logFC",p.value = 0.05, n = 200000)$table
br2_lib_vs_mus_res_lfc2 <- br2_lib_vs_mus_res %>% 
  filter(br2_lib_vs_mus_res$logFC > 2)
View(br2_lib_vs_mus_res_lfc2)
write.csv(br2_lib_vs_mus_res)



####DO not collapse ------ END of 2nd round BRAIN----

####round 2 BRAIN
###----
#create datasets

br2_br_vs_lu <- biopanning_all[,c(1,21,27,33,22,28,34)]
br2_br_vs_li <- biopanning_all[,c(1,21,27,33, 23,29,35)]
br2_br_vs_ki <- biopanning_all[,c(1,21,27,33, 24,30,36)]
br2_br_vs_mus <- biopanning_all[,c(1,21,27,33,25,31,37)]

#create groups
group_biopanning <- factor(c(1,1,1,2,2,2))

#create model matrix
biopanning_design <- model.matrix(~group_biopanning)

#make dge lists

br2_br_vs_lu_dge <- DGEList(counts = br2_br_vs_lu[,2:7], genes = br2_br_vs_lu[,1], group = group_biopanning)
br2_br_vs_li_dge <- DGEList(counts = br2_br_vs_li[,2:7], genes = br2_br_vs_li[,1], group = group_biopanning)
br2_br_vs_ki_dge <- DGEList(counts = br2_br_vs_ki[,2:7], genes = br2_br_vs_ki[,1], group = group_biopanning)
br2_br_vs_mus_dge <- DGEList(counts = br2_br_vs_mus[,2:7], genes = br2_br_vs_mus[,1], group = group_biopanning)

#filter out low abundance peptides

keep_br2_br_vs_lu_dge <- rowSums((br2_br_vs_lu_dge$counts)>5)>= 2
br2_br_vs_lu_dge_filtered <- br2_br_vs_lu_dge[keep_br2_br_vs_lu_dge, , keep.lib.sizes = FALSE]

keep_br2_br_vs_li_dge <- rowSums((br2_br_vs_li_dge$counts)>5)>= 2
br2_br_vs_li_dge_filtered <- br2_br_vs_li_dge[keep_br2_br_vs_li_dge, , keep.lib.sizes = FALSE]

keep_br2_br_vs_ki_dge <- rowSums((br2_br_vs_ki_dge$counts)>5)>= 2
br2_br_vs_ki_dge_filtered <- br2_br_vs_ki_dge[keep_br2_br_vs_ki_dge, , keep.lib.sizes = FALSE]

keep_br2_br_vs_mus_dge <- rowSums((br2_br_vs_mus_dge$counts)>5)>= 2
br2_br_vs_mus_dge_filtered <- br2_br_vs_mus_dge[keep_br2_br_vs_mus_dge, , keep.lib.sizes = FALSE]

nrow(br2_br_vs_mus_dge_filtered)

#normalize samples

br2_br_vs_lu_dge_filtered_norm <- calcNormFactors(br2_br_vs_lu_dge_filtered)
br2_br_vs_li_dge_filtered_norm <- calcNormFactors(br2_br_vs_li_dge_filtered)
br2_br_vs_ki_dge_filtered_norm <- calcNormFactors(br2_br_vs_ki_dge_filtered)
br2_br_vs_mus_dge_filtered_norm <- calcNormFactors(br2_br_vs_mus_dge_filtered)

#estimate dispersions, fit, perform stat tests and save results

###IT IS NECESSARY TO ADD logFC calculations

br2_br_vs_lu_de <- estimateGLMCommonDisp(br2_br_vs_lu_dge_filtered_norm, biopanning_design)
br2_br_vs_lu_de <- estimateGLMTrendedDisp(br2_br_vs_lu_de, biopanning_design)
br2_br_vs_lu_de <- estimateGLMTagwiseDisp(br2_br_vs_lu_de, biopanning_design)
br2_br_vs_lu_fit <- glmQLFit(br2_br_vs_lu_de, biopanning_design)
br2_br_vs_lu_test <-glmQLFTest(br2_br_vs_lu_fit)
br2_br_vs_lu_res <- topTags(br2_br_vs_lu_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
br2_br_vs_lu_res_lfc2 <- br2_br_vs_lu_res %>% 
  filter(br2_br_vs_lu_res$logFC > 2)
br2_br_vs_lu_res_lfc2_minus <- br2_br_vs_lu_res %>% 
  filter(br2_br_vs_lu_res$logFC < -2)
View(br2_br_vs_lu_res_lfc2)
View(br2_br_vs_lu_res_lfc2_minus)

br2_br_vs_li_de <- estimateGLMCommonDisp(br2_br_vs_li_dge_filtered_norm, biopanning_design)
br2_br_vs_li_de <- estimateGLMTrendedDisp(br2_br_vs_li_de, biopanning_design)
br2_br_vs_li_de <- estimateGLMTagwiseDisp(br2_br_vs_li_de, biopanning_design)
br2_br_vs_li_fit <- glmQLFit(br2_br_vs_li_de, biopanning_design)
br2_br_vs_li_test <-glmQLFTest(br2_br_vs_li_fit)
br2_br_vs_li_res <- topTags(br2_br_vs_li_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
br2_br_vs_li_res_lfc2 <- br2_br_vs_li_res %>% 
  filter(br2_br_vs_li_res$logFC > 2)
br2_br_vs_li_res_lfc2_minus <- br2_br_vs_li_res %>% 
  filter(br2_br_vs_li_res$logFC < -2)
View(br2_br_vs_li_res_lfc2)
View(br2_br_vs_li_res_lfc2_minus)

br2_br_vs_ki_de <- estimateGLMCommonDisp(br2_br_vs_ki_dge_filtered_norm, biopanning_design)
br2_br_vs_ki_de <- estimateGLMTrendedDisp(br2_br_vs_ki_de, biopanning_design)
br2_br_vs_ki_de <- estimateGLMTagwiseDisp(br2_br_vs_ki_de, biopanning_design)
br2_br_vs_ki_fit <- glmQLFit(br2_br_vs_ki_de, biopanning_design)
br2_br_vs_ki_test <-glmQLFTest(br2_br_vs_ki_fit)
br2_br_vs_ki_res <- topTags(br2_br_vs_ki_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
br2_br_vs_ki_res_lfc2 <- br2_br_vs_ki_res %>% 
  filter(br2_br_vs_ki_res$logFC > 2)
br2_br_vs_ki_res_lfc2_minus <- br2_br_vs_ki_res %>% 
  filter(br2_br_vs_ki_res$logFC < -2)
View(br2_br_vs_ki_res_lfc2)
View(br2_br_vs_ki_res_lfc2_minus)

br2_br_vs_mus_de <- estimateGLMCommonDisp(br2_br_vs_mus_dge_filtered_norm, biopanning_design)
br2_br_vs_mus_de <- estimateGLMTrendedDisp(br2_br_vs_mus_de, biopanning_design)
br2_br_vs_mus_de <- estimateGLMTagwiseDisp(br2_br_vs_mus_de, biopanning_design)
br2_br_vs_mus_fit <- glmQLFit(br2_br_vs_mus_de, biopanning_design)
br2_br_vs_mus_test <-glmQLFTest(br2_br_vs_mus_fit)
br2_br_vs_mus_res <- topTags(br2_br_vs_mus_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
br2_br_vs_mus_res_lfc2 <- br2_br_vs_mus_res %>% 
  filter(br2_br_vs_mus_res$logFC > 2)
br2_br_vs_mus_res_lfc2_minus <- br2_br_vs_mus_res %>% 
  filter(br2_br_vs_mus_res$logFC < -2)
View(br2_br_vs_mus_res_lfc2)
View(br2_br_vs_mus_res_lfc2_minus)

#create combined file

colnames(br2_br_vs_lu_res_lfc2_minus) <- c("pep","logFC_lu", "logCPM_lu", "F_lu", "PValue_lu", "FDR_lu") 
colnames(br2_br_vs_li_res_lfc2_minus) <- c("pep","logFC_li", "logCPM_li", "F_li", "PValue_li", "FDR_li") 
colnames(br2_br_vs_ki_res_lfc2_minus) <- c("pep","logFC_ki", "logCPM_ki", "F_ki", "PValue_ki", "FDR_ki") 
colnames(br2_br_vs_mus_res_lfc2_minus) <- c("pep","logFC_mu", "logCPM_mu", "F_mu", "PValue_mu", "FDR_mu") 

br2_joined_1 <- full_join(br2_br_vs_lu_res_lfc2_minus, br2_br_vs_li_res_lfc2_minus, by = "pep")
br2_joined_2 <- full_join(br2_joined_1, br2_br_vs_ki_res_lfc2_minus, by = "pep")
br2_joined_all <- full_join(br2_joined_2, br2_br_vs_mus_res_lfc2_minus, by = "pep")
View(br2_joined_all)
differential_binding_seq_brain <- br2_joined_all[,c(1, 2,3,6,7,8,11,12,13,16,17,18,21)]
write.csv(differential_binding_seq_brain, file = "Round2_differential_binding_brain_seq.csv")


area_br2_br_vs_lu <- nrow(subset(br2_joined_all, br2_joined_all$logFC_lu != 0))
area_br2_br_vs_li <- nrow(subset(br2_joined_all, br2_joined_all$logFC_li != 0))
area_br2_br_vs_ki <- nrow(subset(br2_joined_all, br2_joined_all$logFC_ki != 0))
area_br2_br_vs_mu <- nrow(subset(br2_joined_all, br2_joined_all$logFC_mu != 0))


area_br2_br_vs_lu_ki <- nrow(subset(br2_joined_all, br2_joined_all$logFC_lu != 0 & br2_joined_all$logCPM_ki != 0))
area_br2_br_vs_lu_li <- nrow(subset(br2_joined_all, br2_joined_all$logFC_lu != 0 & br2_joined_all$logCPM_li != 0))
area_br2_br_vs_lu_mu <- nrow(subset(br2_joined_all, br2_joined_all$logFC_lu != 0 & br2_joined_all$logCPM_mu != 0))
area_br2_br_vs_ki_li <- nrow(subset(br2_joined_all, br2_joined_all$logFC_ki != 0 & br2_joined_all$logCPM_li != 0))
area_br2_br_vs_ki_mu <- nrow(subset(br2_joined_all, br2_joined_all$logFC_ki != 0 & br2_joined_all$logCPM_mu != 0))
area_br2_br_vs_li_mu <- nrow(subset(br2_joined_all, br2_joined_all$logFC_li != 0 & br2_joined_all$logCPM_mu != 0))


area_br2_br_vs_lu_mu_li <- nrow(subset(br2_joined_all, br2_joined_all$logFC_lu != 0 & br2_joined_all$logCPM_mu != 0
                                       & br2_joined_all$logCPM_li != 0))
area_br2_br_vs_lu_ki_li <- nrow(subset(br2_joined_all, br2_joined_all$logFC_lu != 0 & br2_joined_all$logCPM_ki != 0
                                       & br2_joined_all$logCPM_li != 0))
area_br2_br_vs_ki_mu_lu <- nrow(subset(br2_joined_all, br2_joined_all$logFC_ki != 0 & br2_joined_all$logCPM_mu != 0
                                       & br2_joined_all$logCPM_lu != 0))
area_br2_br_vs_li_ki_mu <- nrow(subset(br2_joined_all, br2_joined_all$logFC_li != 0 & br2_joined_all$logCPM_ki != 0
                                       & br2_joined_all$logCPM_mu != 0))

area_br2_br_vs_lu_mu_ki_li <- nrow(subset(br2_joined_all, br2_joined_all$logFC_lu != 0 & br2_joined_all$logCPM_mu != 0
                                          & br2_joined_all$logCPM_ki != 0 & br2_joined_all$logCPM_li))

aaa <- subset(br2_joined_all, br2_joined_all$logFC_ki != 0 & br2_joined_all$logCPM_li != 0
              & br2_joined_all$logCPM_mu != 0)

View(aaa)
aaa_2 <- aaa[,c(1, 2,3,6,7,8,11,12,13,16,17,18,21)]
View(aaa_2)
write.csv(aaa_2, file = "Differential_binding_br2_most_promising_leads.csv")


###UpSetR analysis 

upset_input_br2 <- c("Lung"=area_br2_br_vs_lu, "Liver"=area_br2_br_vs_li, "Kidney"=area_br2_br_vs_ki, "Muscle" = area_br2_br_vs_mu,
                   "Lung&Kidney"=area_br2_br_vs_lu_ki, "Lung&Liver"=area_br2_br_vs_lu_li, "Lung&Muscle"=area_br2_br_vs_lu_mu,
                   "Kidney&Liver" = area_br2_br_vs_ki_li,"Kidney&Muscle" = area_br2_br_vs_ki_mu,
                   "Liver&Muscle" = area_br2_br_vs_li_mu,"Lung&Muscle&Kidney" = area_br2_br_vs_ki_mu_lu,
                   "Lung&Muscle&Liver" = area_br2_br_vs_lu_mu_li, "Lung&Kidney&Liver" = area_br2_br_vs_lu_ki_li,
                   "Liver&Kidney&Muscle" = area_br2_br_vs_li_ki_mu,
                   "Lung&Muscle&Kidney&Liver"=area_br2_br_vs_lu_mu_ki_li)

upset(fromExpression(upset_input_br2), order.by = "freq")

##### ROUND 2 LUNG
lu2_lu_vs_br <- biopanning_all[,c(1,40,46,52,39,45,51)]
lu2_lu_vs_li <- biopanning_all[,c(1,40,46,52,41,47,53)]
lu2_lu_vs_ki <- biopanning_all[,c(1,40,46,52,42,48,54)]
lu2_lu_vs_mu <- biopanning_all[,c(1,40,46,52,43,49,55)]

#create groups
group_biopanning <- factor(c(1,1,1,2,2,2))

#create model matrix
biopanning_design <- model.matrix(~group_biopanning)

#make dge lists
lu2_lu_vs_br_dge <- DGEList(counts = lu2_lu_vs_br[,2:7], genes = lu2_lu_vs_br[,1], group = group_biopanning)
lu2_lu_vs_li_dge <- DGEList(counts = lu2_lu_vs_li[,2:7], genes = lu2_lu_vs_li[,1], group = group_biopanning)
lu2_lu_vs_ki_dge <- DGEList(counts = lu2_lu_vs_ki[,2:7], genes = lu2_lu_vs_ki[,1], group = group_biopanning)
lu2_lu_vs_mu_dge <- DGEList(counts = lu2_lu_vs_mu[,2:7], genes = lu2_lu_vs_mu[,1], group = group_biopanning)



#filter out low abundance peptides

keep_lu2_lu_vs_br_dge <- rowSums((lu2_lu_vs_br_dge$counts)>5)>= 2
lu2_lu_vs_br_dge_filtered <- lu2_lu_vs_br_dge[keep_lu2_lu_vs_br_dge, , keep.lib.sizes = FALSE]

keep_lu2_lu_vs_li_dge <- rowSums((lu2_lu_vs_li_dge$counts)>5)>= 2
lu2_lu_vs_li_dge_filtered <- lu2_lu_vs_li_dge[keep_lu2_lu_vs_li_dge, , keep.lib.sizes = FALSE]

keep_lu2_lu_vs_ki_dge <- rowSums((lu2_lu_vs_ki_dge$counts)>5)>= 2
lu2_lu_vs_ki_dge_filtered <- lu2_lu_vs_ki_dge[keep_lu2_lu_vs_ki_dge, , keep.lib.sizes = FALSE]

keep_lu2_lu_vs_mu_dge <- rowSums((lu2_lu_vs_mu_dge$counts)>5)>= 2
lu2_lu_vs_mu_dge_filtered <- lu2_lu_vs_mu_dge[keep_lu2_lu_vs_mu_dge, , keep.lib.sizes = FALSE]

nrow(lu2_lu_vs_ki_dge_filtered)

#normalize samples
lu2_lu_vs_br_dge_filtered_norm <- calcNormFactors(lu2_lu_vs_br_dge_filtered)
lu2_lu_vs_li_dge_filtered_norm <- calcNormFactors(lu2_lu_vs_li_dge_filtered)
lu2_lu_vs_ki_dge_filtered_norm <- calcNormFactors(lu2_lu_vs_ki_dge_filtered)
lu2_lu_vs_mu_dge_filtered_norm <- calcNormFactors(lu2_lu_vs_mu_dge_filtered)

#estimate dispersions, fit, perform stat tests and save results

###IT IS NECESSARY TO ADD logFC calculations

lu2_lu_vs_br_de <- estimateGLMCommonDisp(lu2_lu_vs_br_dge_filtered_norm, biopanning_design)
lu2_lu_vs_br_de <- estimateGLMTrendedDisp(lu2_lu_vs_br_de, biopanning_design)
lu2_lu_vs_br_de <- estimateGLMTagwiseDisp(lu2_lu_vs_br_de, biopanning_design)
lu2_lu_vs_br_fit <- glmQLFit(lu2_lu_vs_br_de, biopanning_design)
lu2_lu_vs_br_test <-glmQLFTest(lu2_lu_vs_br_fit)
lu2_lu_vs_br_res <- topTags(lu2_lu_vs_br_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
lu2_lu_vs_br_res_lfc2 <- lu2_lu_vs_br_res %>% 
  filter(lu2_lu_vs_br_res$logFC > 2)
lu2_lu_vs_br_res_lfc2_minus <- lu2_lu_vs_br_res %>% 
  filter(lu2_lu_vs_br_res$logFC < -2)
View(lu2_lu_vs_br_res_lfc2)
View(lu2_lu_vs_br_res_lfc2_minus)

lu2_lu_vs_li_de <- estimateGLMCommonDisp(lu2_lu_vs_li_dge_filtered_norm, biopanning_design)
lu2_lu_vs_li_de <- estimateGLMTrendedDisp(lu2_lu_vs_li_de, biopanning_design)
lu2_lu_vs_li_de <- estimateGLMTagwiseDisp(lu2_lu_vs_li_de, biopanning_design)
lu2_lu_vs_li_fit <- glmQLFit(lu2_lu_vs_li_de, biopanning_design)
lu2_lu_vs_li_test <-glmQLFTest(lu2_lu_vs_li_fit)
lu2_lu_vs_li_res <- topTags(lu2_lu_vs_li_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
lu2_lu_vs_li_res_lfc2 <- lu2_lu_vs_li_res %>% 
  filter(lu2_lu_vs_li_res$logFC > 2)
lu2_lu_vs_li_res_lfc2_minus <- lu2_lu_vs_li_res %>% 
  filter(lu2_lu_vs_li_res$logFC < -2)
View(lu2_lu_vs_li_res_lfc2)
View(lu2_lu_vs_li_res_lfc2_minus)

lu2_lu_vs_ki_de <- estimateGLMCommonDisp(lu2_lu_vs_ki_dge_filtered_norm, biopanning_design)
lu2_lu_vs_ki_de <- estimateGLMTrendedDisp(lu2_lu_vs_ki_de, biopanning_design)
lu2_lu_vs_ki_de <- estimateGLMTagwiseDisp(lu2_lu_vs_ki_de, biopanning_design)
lu2_lu_vs_ki_fit <- glmQLFit(lu2_lu_vs_ki_de, biopanning_design)
lu2_lu_vs_ki_test <-glmQLFTest(lu2_lu_vs_ki_fit)
lu2_lu_vs_ki_res <- topTags(lu2_lu_vs_ki_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
lu2_lu_vs_ki_res_lfc2 <- lu2_lu_vs_ki_res %>% 
  filter(lu2_lu_vs_ki_res$logFC > 2)
lu2_lu_vs_ki_res_lfc2_minus <- lu2_lu_vs_ki_res %>% 
  filter(lu2_lu_vs_ki_res$logFC < -2)
View(lu2_lu_vs_ki_res_lfc2)
View(lu2_lu_vs_ki_res_lfc2_minus)

lu2_lu_vs_mu_de <- estimateGLMCommonDisp(lu2_lu_vs_mu_dge_filtered_norm, biopanning_design)
lu2_lu_vs_mu_de <- estimateGLMTrendedDisp(lu2_lu_vs_mu_de, biopanning_design)
lu2_lu_vs_mu_de <- estimateGLMTagwiseDisp(lu2_lu_vs_mu_de, biopanning_design)
lu2_lu_vs_mu_fit <- glmQLFit(lu2_lu_vs_mu_de, biopanning_design)
lu2_lu_vs_mu_test <-glmQLFTest(lu2_lu_vs_mu_fit)
lu2_lu_vs_mu_res <- topTags(lu2_lu_vs_mu_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
lu2_lu_vs_mu_res_lfc2 <- lu2_lu_vs_mu_res %>% 
  filter(lu2_lu_vs_mu_res$logFC > 2)
lu2_lu_vs_mu_res_lfc2_minus <- lu2_lu_vs_mu_res %>% 
  filter(lu2_lu_vs_mu_res$logFC < -2)
View(lu2_lu_vs_mu_res_lfc2)
View(lu2_lu_vs_mu_res_lfc2_minus)



#create combined file

colnames(lu2_lu_vs_br_res_lfc2_minus) <- c("pep","logFC_br", "logCPM_br", "F_br", "PValue_br", "FDR_br")
colnames(lu2_lu_vs_li_res_lfc2_minus) <- c("pep","logFC_li", "logCPM_li", "F_li", "PValue_li", "FDR_li")
colnames(lu2_lu_vs_ki_res_lfc2_minus) <- c("pep","logFC_ki", "logCPM_ki", "F_ki", "PValue_ki", "FDR_ki")
colnames(lu2_lu_vs_mu_res_lfc2_minus) <- c("pep","logFC_mu", "logCPM_mu", "F_mu", "PValue_mu", "FDR_mu")

lu2_lu_vs_mu_res_lfc2_minus

lu2_joined_1 <- full_join(lu2_lu_vs_br_res_lfc2_minus, lu2_lu_vs_li_res_lfc2_minus,by = "pep")
lu2_joined_2 <- full_join(lu2_joined_1, lu2_lu_vs_ki_res_lfc2_minus, by = "pep")
lu2_joined_all <- full_join(lu2_joined_2, lu2_lu_vs_mu_res_lfc2_minus, by = "pep")
View(lu2_joined_all)
differential_binding_seq_lung <- lu2_joined_all[,c(1, 2,3,6,7,8,11,12,13,16,17,18,21)]
write.csv(differential_binding_seq_lung, file = "Round2_differential_binding_lung_seq.csv")


#create input for setupr

area_lu2_lu_vs_br <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_br != 0))
area_lu2_lu_vs_li <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_li != 0))
area_lu2_lu_vs_ki <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_ki != 0))
area_lu2_lu_vs_mu <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_mu != 0))

area_lu2_lu_vs_br_ki <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_br != 0 & lu2_joined_all$logCPM_ki != 0))
area_lu2_lu_vs_br_li <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_br != 0 & lu2_joined_all$logCPM_li != 0))
area_lu2_lu_vs_br_mu <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_br != 0 & lu2_joined_all$logCPM_mu != 0))
area_lu2_lu_vs_ki_li <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_ki != 0 & lu2_joined_all$logCPM_li != 0))
area_lu2_lu_vs_ki_mu <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_ki != 0 & lu2_joined_all$logCPM_mu != 0))
area_lu2_lu_vs_li_mu <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_li != 0 & lu2_joined_all$logCPM_mu != 0))



area_lu2_lu_vs_br_ki_li <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_br != 0 & lu2_joined_all$logCPM_ki != 0
                                       & lu2_joined_all$logCPM_li != 0))
area_lu2_lu_vs_br_li_mu <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_br != 0 & lu2_joined_all$logCPM_li != 0
                                       & lu2_joined_all$logCPM_mu != 0))
area_lu2_lu_vs_br_mu_ki <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_br != 0 & lu2_joined_all$logCPM_mu != 0
                                       & lu2_joined_all$logCPM_ki != 0))
area_lu2_lu_vs_ki_li_mu <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_ki != 0 & lu2_joined_all$logCPM_li != 0
                                       & lu2_joined_all$logCPM_mu != 0))


area_lu2_lu_vs_br_mu_ki_li <- nrow(subset(lu2_joined_all, lu2_joined_all$logFC_br != 0 & lu2_joined_all$logCPM_mu != 0
                                          & lu2_joined_all$logCPM_ki != 0 & lu2_joined_all$logCPM_li))

aaa_lu2 <- subset(lu2_joined_all, lu2_joined_all$logFC_mu != 0 & lu2_joined_all$logCPM_br != 0
              & lu2_joined_all$logCPM_li != 0)

View(aaa_lu2)
aaa_lu2_2 <- aaa_lu2[,c(1, 2,3,6,7,8,11,12,13,16,17,18,21)]
View(aaa_lu2)
write.csv(aaa_lu2, file = "Differential_binding_lu2_most_promising_leads.csv")

###UpSetR analysis 

upset_input_lu2 <- c("Brain"=area_lu2_lu_vs_br, "Liver"=area_lu2_lu_vs_li, "Kidney"=area_lu2_lu_vs_ki, "Muscle" = area_lu2_lu_vs_mu,
                     "Brain&Kidney"=area_lu2_lu_vs_br_ki, "Brain&Liver"=area_lu2_lu_vs_br_li, "Brain&Muscle"=area_lu2_lu_vs_br_mu,
                     "Kidney&Liver" = area_lu2_lu_vs_ki_li,"Kidney&Muscle" = area_lu2_lu_vs_ki_mu,
                     "Liver&Muscle" = area_lu2_lu_vs_li_mu,
                     "Brain&Muscle&Kidney" = area_lu2_lu_vs_br_mu_ki,
                     "Brain&Muscle&Liver" = area_lu2_lu_vs_br_li_mu, 
                     "Brain&Kidney&Liver" = area_lu2_lu_vs_br_ki_li,
                     "Kidney&Muscle&Liver" = area_lu2_lu_vs_ki_li_mu,
                     "Brain&Muscle&Kidney&Liver"=area_lu2_lu_vs_br_mu_ki_li)

upset(fromExpression(upset_input_lu2), order.by = "freq")

##ROUND 1 BRAIN

#create datasets

br1_br_vs_lu <- biopanning_all[,c(1, 3,9,15,4,10,16)]
br1_br_vs_li <- biopanning_all[,c(1, 3,9,15,5,11,17)]
br1_br_vs_ki <- biopanning_all[,c(1, 3,9,15,6,12,18)]
br1_br_vs_mu <- biopanning_all[,c(1, 3,9,15,7,13,19)]

br1_br_vs_lu

#create groups
group_biopanning <- factor(c(1,1,1,2,2,2))

#create model matrix
biopanning_design <- model.matrix(~group_biopanning)

#make dge lists

br1_br_vs_lu_dge <- DGEList(counts = br1_br_vs_lu[,2:7], genes = br1_br_vs_lu[,1], group = group_biopanning)
br1_br_vs_li_dge <- DGEList(counts = br1_br_vs_li[,2:7], genes = br1_br_vs_li[,1], group = group_biopanning)
br1_br_vs_ki_dge <- DGEList(counts = br1_br_vs_ki[,2:7], genes = br1_br_vs_ki[,1], group = group_biopanning)
br1_br_vs_mu_dge <- DGEList(counts = br1_br_vs_mu[,2:7], genes = br1_br_vs_mu[,1], group = group_biopanning)

#filter out low abundance peptides

keep_br1_br_vs_lu_dge <- rowSums((br1_br_vs_lu_dge$counts)>5)>= 2
br1_br_vs_lu_dge_filtered <- br1_br_vs_lu_dge[keep_br1_br_vs_lu_dge, , keep.lib.sizes = FALSE]

keep_br1_br_vs_li_dge <- rowSums((br1_br_vs_li_dge$counts)>5)>= 2
br1_br_vs_li_dge_filtered <- br1_br_vs_li_dge[keep_br1_br_vs_li_dge, , keep.lib.sizes = FALSE]

keep_br1_br_vs_ki_dge <- rowSums((br1_br_vs_ki_dge$counts)>5)>= 2
br1_br_vs_ki_dge_filtered <- br1_br_vs_ki_dge[keep_br1_br_vs_ki_dge, , keep.lib.sizes = FALSE]

keep_br1_br_vs_mu_dge <- rowSums((br1_br_vs_mu_dge$counts)>5)>= 2
br1_br_vs_mu_dge_filtered <- br1_br_vs_mu_dge[keep_br1_br_vs_mu_dge, , keep.lib.sizes = FALSE]

nrow(br1_br_vs_mu_dge_filtered)

#normalize samples

br1_br_vs_lu_dge_filtered_norm <- calcNormFactors(br1_br_vs_lu_dge_filtered)
br1_br_vs_li_dge_filtered_norm <- calcNormFactors(br1_br_vs_li_dge_filtered)
br1_br_vs_ki_dge_filtered_norm <- calcNormFactors(br1_br_vs_ki_dge_filtered)
br1_br_vs_mu_dge_filtered_norm <- calcNormFactors(br1_br_vs_mu_dge_filtered)

#estimate dispersions, fit, perform stat tests and save results

###IT IS NECESSARY TO ADD logFC calculations

br1_br_vs_lu_de <- estimateGLMCommonDisp(br1_br_vs_lu_dge_filtered_norm, biopanning_design)
br1_br_vs_lu_de <- estimateGLMTrendedDisp(br1_br_vs_lu_de, biopanning_design)
br1_br_vs_lu_de <- estimateGLMTagwiseDisp(br1_br_vs_lu_de, biopanning_design)
br1_br_vs_lu_fit <- glmQLFit(br1_br_vs_lu_de, biopanning_design)
br1_br_vs_lu_test <-glmQLFTest(br1_br_vs_lu_fit)
br1_br_vs_lu_res <- topTags(br1_br_vs_lu_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
br1_br_vs_lu_res_lfc2 <- br1_br_vs_lu_res %>% 
  filter(br1_br_vs_lu_res$logFC > 2)
br1_br_vs_lu_res_lfc2_minus <- br1_br_vs_lu_res %>% 
  filter(br1_br_vs_lu_res$logFC < -2)
View(br1_br_vs_lu_res_lfc2)
View(br1_br_vs_lu_res_lfc2_minus)

br1_br_vs_li_de <- estimateGLMCommonDisp(br1_br_vs_li_dge_filtered_norm, biopanning_design)
br1_br_vs_li_de <- estimateGLMTrendedDisp(br1_br_vs_li_de, biopanning_design)
br1_br_vs_li_de <- estimateGLMTagwiseDisp(br1_br_vs_li_de, biopanning_design)
br1_br_vs_li_fit <- glmQLFit(br1_br_vs_li_de, biopanning_design)
br1_br_vs_li_test <-glmQLFTest(br1_br_vs_li_fit)
br1_br_vs_li_res <- topTags(br1_br_vs_li_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
br1_br_vs_li_res_lfc2 <- br1_br_vs_li_res %>% 
  filter(br1_br_vs_li_res$logFC > 2)
br1_br_vs_li_res_lfc2_minus <- br1_br_vs_li_res %>% 
  filter(br1_br_vs_li_res$logFC < -2)
View(br1_br_vs_li_res_lfc2)
View(br1_br_vs_li_res_lfc2_minus)

br1_br_vs_ki_de <- estimateGLMCommonDisp(br1_br_vs_ki_dge_filtered_norm, biopanning_design)
br1_br_vs_ki_de <- estimateGLMTrendedDisp(br1_br_vs_ki_de, biopanning_design)
br1_br_vs_ki_de <- estimateGLMTagwiseDisp(br1_br_vs_ki_de, biopanning_design)
br1_br_vs_ki_fit <- glmQLFit(br1_br_vs_ki_de, biopanning_design)
br1_br_vs_ki_test <-glmQLFTest(br1_br_vs_ki_fit)
br1_br_vs_ki_res <- topTags(br1_br_vs_ki_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
br1_br_vs_ki_res_lfc2 <- br1_br_vs_ki_res %>% 
  filter(br1_br_vs_ki_res$logFC > 2)
br1_br_vs_ki_res_lfc2_minus <- br1_br_vs_ki_res %>% 
  filter(br1_br_vs_ki_res$logFC < -2)
View(br1_br_vs_ki_res_lfc2)
View(br1_br_vs_ki_res_lfc2_minus)

br1_br_vs_mu_de <- estimateGLMCommonDisp(br1_br_vs_mu_dge_filtered_norm, biopanning_design)
br1_br_vs_mu_de <- estimateGLMTrendedDisp(br1_br_vs_mu_de, biopanning_design)
br1_br_vs_mu_de <- estimateGLMTagwiseDisp(br1_br_vs_mu_de, biopanning_design)
br1_br_vs_mu_fit <- glmQLFit(br1_br_vs_mu_de, biopanning_design)
br1_br_vs_mu_test <-glmQLFTest(br1_br_vs_mu_fit)
br1_br_vs_mu_res <- topTags(br1_br_vs_mu_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
br1_br_vs_mu_res_lfc2 <- br1_br_vs_mu_res %>% 
  filter(br1_br_vs_mu_res$logFC > 2)
br1_br_vs_mu_res_lfc2_minus <- br1_br_vs_mu_res %>% 
  filter(br1_br_vs_mu_res$logFC < -2)
View(br1_br_vs_mu_res_lfc2)
View(br1_br_vs_mu_res_lfc2_minus)

#create combined file

colnames(br1_br_vs_lu_res_lfc2_minus) <- c("pep","logFC_lu", "logCPM_lu", "F_lu", "PValue_lu", "FDR_lu") 
colnames(br1_br_vs_li_res_lfc2_minus) <- c("pep","logFC_li", "logCPM_li", "F_li", "PValue_li", "FDR_li") 
colnames(br1_br_vs_ki_res_lfc2_minus) <- c("pep","logFC_ki", "logCPM_ki", "F_ki", "PValue_ki", "FDR_ki") 
colnames(br1_br_vs_mu_res_lfc2_minus) <- c("pep","logFC_mu", "logCPM_mu", "F_mu", "PValue_mu", "FDR_mu") 

br1_joined_1 <- full_join(br1_br_vs_lu_res_lfc2_minus, br1_br_vs_li_res_lfc2_minus, by = "pep")
br1_joined_2 <- full_join(br1_joined_1, br1_br_vs_ki_res_lfc2_minus, by = "pep")
br1_joined_all <- full_join(br1_joined_2, br1_br_vs_mu_res_lfc2_minus, by = "pep")
View(br1_joined_all)

area_br1_br_vs_lu <- nrow(subset(br1_joined_all, br1_joined_all$logFC_lu != 0))
area_br1_br_vs_li <- nrow(subset(br1_joined_all, br1_joined_all$logFC_li != 0))
area_br1_br_vs_ki <- nrow(subset(br1_joined_all, br1_joined_all$logFC_ki != 0))
area_br1_br_vs_mu <- nrow(subset(br1_joined_all, br1_joined_all$logFC_mu != 0))


area_br1_br_vs_lu_ki <- nrow(subset(br1_joined_all, br1_joined_all$logFC_lu != 0 & br1_joined_all$logCPM_ki != 0))
area_br1_br_vs_lu_li <- nrow(subset(br1_joined_all, br1_joined_all$logFC_lu != 0 & br1_joined_all$logCPM_li != 0))
area_br1_br_vs_lu_mu <- nrow(subset(br1_joined_all, br1_joined_all$logFC_lu != 0 & br1_joined_all$logCPM_mu != 0))
area_br1_br_vs_ki_li <- nrow(subset(br1_joined_all, br1_joined_all$logFC_ki != 0 & br1_joined_all$logCPM_li != 0))
area_br1_br_vs_ki_mu <- nrow(subset(br1_joined_all, br1_joined_all$logFC_ki != 0 & br1_joined_all$logCPM_mu != 0))
area_br1_br_vs_li_mu <- nrow(subset(br1_joined_all, br1_joined_all$logFC_li != 0 & br1_joined_all$logCPM_mu != 0))

area_br1_br_vs_lu_mu_ki <- nrow(subset(br1_joined_all, br1_joined_all$logFC_lu != 0 & br1_joined_all$logCPM_mu != 0
                                       & br1_joined_all$logCPM_ki != 0))
area_br1_br_vs_lu_mu_li <- nrow(subset(br1_joined_all, br1_joined_all$logFC_lu != 0 & br1_joined_all$logCPM_mu != 0
                                       & br1_joined_all$logCPM_li != 0))
area_br1_br_vs_lu_ki_li <- nrow(subset(br1_joined_all, br1_joined_all$logFC_lu != 0 & br1_joined_all$logCPM_ki != 0
                                       & br1_joined_all$logCPM_li != 0))
area_br1_br_vs_ki_mu_lu <- nrow(subset(br1_joined_all, br1_joined_all$logFC_ki != 0 & br1_joined_all$logCPM_mu != 0
                                       & br1_joined_all$logCPM_lu != 0))

area_br1_br_vs_lu_mu_ki_li <- nrow(subset(br1_joined_all, br1_joined_all$logFC_lu != 0 & br1_joined_all$logCPM_mu != 0
                                          & br1_joined_all$logCPM_ki != 0 & br1_joined_all$logCPM_li))


###UpSetR analysis 

upset_input_br1 <- c("Lung"=area_br1_br_vs_lu, "Liver"=area_br1_br_vs_li, "Kidney"=area_br1_br_vs_ki, "Muscle" = area_br1_br_vs_mu,
                     "Lung&Liver"=area_br1_br_vs_lu_li, "Lung&Liver"=area_br1_br_vs_lu_li, "Lung&Muscle"=area_br1_br_vs_lu_mu,
                     "Kidney&Liver" = area_br1_br_vs_ki_li,"Kidney&Muscle" = area_br1_br_vs_ki_mu,
                     "Liver&Muscle" = area_br1_br_vs_li_mu,"Lung&Muscle&Kidney" = area_br1_br_vs_lu_mu_ki,
                     "Lung&Muscle&Liver" = area_br1_br_vs_lu_mu_li, "Lung&Kidney&Liver" = area_br1_br_vs_lu_ki_li,
                     "Kidney&Muscle&Lung" = area_br1_br_vs_ki_mu_lu,
                     "Lung&Muscle&Kidney&Liver"=area_br1_br_vs_lu_mu_ki_li)

upset(fromExpression(upset_input_br1), order.by = "freq")

##### ROUND 1 LUNG

lu1_lu_vs_br <- biopanning_all[,c(1,4,10,16, 3,9,15)]
lu1_lu_vs_li <- biopanning_all[,c(1,4,10,16, 5,11,17)]
lu1_lu_vs_ki <- biopanning_all[,c(1,4,10,16, 6,12,18)]
lu1_lu_vs_mu <- biopanning_all[,c(1,4,10,16, 7,13,19)]

#create groups
group_biopanning <- factor(c(1,1,1,2,2,2))

#create model matrix
biopanning_design <- model.matrix(~group_biopanning)

#make dge lists
lu1_lu_vs_br_dge <- DGEList(counts = lu1_lu_vs_br[,2:7], genes = lu1_lu_vs_br[,1], group = group_biopanning)
lu1_lu_vs_li_dge <- DGEList(counts = lu1_lu_vs_li[,2:7], genes = lu1_lu_vs_li[,1], group = group_biopanning)
lu1_lu_vs_ki_dge <- DGEList(counts = lu1_lu_vs_ki[,2:7], genes = lu1_lu_vs_ki[,1], group = group_biopanning)
lu1_lu_vs_mu_dge <- DGEList(counts = lu1_lu_vs_mu[,2:7], genes = lu1_lu_vs_mu[,1], group = group_biopanning)



#filter out low abundance peptides

keep_lu1_lu_vs_br_dge <- rowSums((lu1_lu_vs_br_dge$counts)>5)>= 2
lu1_lu_vs_br_dge_filtered <- lu1_lu_vs_br_dge[keep_lu1_lu_vs_br_dge, , keep.lib.sizes = FALSE]

keep_lu1_lu_vs_li_dge <- rowSums((lu1_lu_vs_li_dge$counts)>5)>= 2
lu1_lu_vs_li_dge_filtered <- lu1_lu_vs_li_dge[keep_lu1_lu_vs_li_dge, , keep.lib.sizes = FALSE]

keep_lu1_lu_vs_ki_dge <- rowSums((lu1_lu_vs_ki_dge$counts)>5)>= 2
lu1_lu_vs_ki_dge_filtered <- lu1_lu_vs_ki_dge[keep_lu1_lu_vs_ki_dge, , keep.lib.sizes = FALSE]

keep_lu1_lu_vs_mu_dge <- rowSums((lu1_lu_vs_mu_dge$counts)>5)>= 2
lu1_lu_vs_mu_dge_filtered <- lu1_lu_vs_mu_dge[keep_lu1_lu_vs_mu_dge, , keep.lib.sizes = FALSE]

nrow(lu1_lu_vs_ki_dge_filtered)

#normalize samples
lu1_lu_vs_br_dge_filtered_norm <- calcNormFactors(lu1_lu_vs_br_dge_filtered)
lu1_lu_vs_li_dge_filtered_norm <- calcNormFactors(lu1_lu_vs_li_dge_filtered)
lu1_lu_vs_ki_dge_filtered_norm <- calcNormFactors(lu1_lu_vs_ki_dge_filtered)
lu1_lu_vs_mu_dge_filtered_norm <- calcNormFactors(lu1_lu_vs_mu_dge_filtered)

#estimate dispersions, fit, perform stat tests and save results

###IT IS NECESSARY TO ADD logFC calculations

lu1_lu_vs_br_de <- estimateGLMCommonDisp(lu1_lu_vs_br_dge_filtered_norm, biopanning_design)
lu1_lu_vs_br_de <- estimateGLMTrendedDisp(lu1_lu_vs_br_de, biopanning_design)
lu1_lu_vs_br_de <- estimateGLMTagwiseDisp(lu1_lu_vs_br_de, biopanning_design)
lu1_lu_vs_br_fit <- glmQLFit(lu1_lu_vs_br_de, biopanning_design)
lu1_lu_vs_br_test <-glmQLFTest(lu1_lu_vs_br_fit)
lu1_lu_vs_br_res <- topTags(lu1_lu_vs_br_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
lu1_lu_vs_br_res_lfc2 <- lu1_lu_vs_br_res %>% 
  filter(lu1_lu_vs_br_res$logFC > 2)
lu1_lu_vs_br_res_lfc2_minus <- lu1_lu_vs_br_res %>% 
  filter(lu1_lu_vs_br_res$logFC < -2)
View(lu1_lu_vs_br_res_lfc2)
View(lu1_lu_vs_br_res_lfc2_minus)

lu1_lu_vs_li_de <- estimateGLMCommonDisp(lu1_lu_vs_li_dge_filtered_norm, biopanning_design)
lu1_lu_vs_li_de <- estimateGLMTrendedDisp(lu1_lu_vs_li_de, biopanning_design)
lu1_lu_vs_li_de <- estimateGLMTagwiseDisp(lu1_lu_vs_li_de, biopanning_design)
lu1_lu_vs_li_fit <- glmQLFit(lu1_lu_vs_li_de, biopanning_design)
lu1_lu_vs_li_test <-glmQLFTest(lu1_lu_vs_li_fit)
lu1_lu_vs_li_res <- topTags(lu1_lu_vs_li_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
lu1_lu_vs_li_res_lfc2 <- lu1_lu_vs_li_res %>% 
  filter(lu1_lu_vs_li_res$logFC > 2)
lu1_lu_vs_li_res_lfc2_minus <- lu1_lu_vs_li_res %>% 
  filter(lu1_lu_vs_li_res$logFC < -2)
View(lu1_lu_vs_li_res_lfc2)
View(lu1_lu_vs_li_res_lfc2_minus)

lu1_lu_vs_ki_de <- estimateGLMCommonDisp(lu1_lu_vs_ki_dge_filtered_norm, biopanning_design)
lu1_lu_vs_ki_de <- estimateGLMTrendedDisp(lu1_lu_vs_ki_de, biopanning_design)
lu1_lu_vs_ki_de <- estimateGLMTagwiseDisp(lu1_lu_vs_ki_de, biopanning_design)
lu1_lu_vs_ki_fit <- glmQLFit(lu1_lu_vs_ki_de, biopanning_design)
lu1_lu_vs_ki_test <-glmQLFTest(lu1_lu_vs_ki_fit)
lu1_lu_vs_ki_res <- topTags(lu1_lu_vs_ki_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
lu1_lu_vs_ki_res_lfc2 <- lu1_lu_vs_ki_res %>% 
  filter(lu1_lu_vs_ki_res$logFC > 2)
lu1_lu_vs_ki_res_lfc2_minus <- lu1_lu_vs_ki_res %>% 
  filter(lu1_lu_vs_ki_res$logFC < -2)
View(lu1_lu_vs_ki_res_lfc2)
View(lu1_lu_vs_ki_res_lfc2_minus)

lu1_lu_vs_mu_de <- estimateGLMCommonDisp(lu1_lu_vs_mu_dge_filtered_norm, biopanning_design)
lu1_lu_vs_mu_de <- estimateGLMTrendedDisp(lu1_lu_vs_mu_de, biopanning_design)
lu1_lu_vs_mu_de <- estimateGLMTagwiseDisp(lu1_lu_vs_mu_de, biopanning_design)
lu1_lu_vs_mu_fit <- glmQLFit(lu1_lu_vs_mu_de, biopanning_design)
lu1_lu_vs_mu_test <-glmQLFTest(lu1_lu_vs_mu_fit)
lu1_lu_vs_mu_res <- topTags(lu1_lu_vs_mu_test, sort.by = "logFC", n = 200000, p.value = 0.05)$table
lu1_lu_vs_mu_res_lfc2 <- lu1_lu_vs_mu_res %>% 
  filter(lu1_lu_vs_mu_res$logFC > 2)
lu1_lu_vs_mu_res_lfc2_minus <- lu1_lu_vs_mu_res %>% 
  filter(lu1_lu_vs_mu_res$logFC < -2)
View(lu1_lu_vs_mu_res_lfc2)
View(lu1_lu_vs_mu_res_lfc2_minus)



#create combined file

colnames(lu1_lu_vs_br_res_lfc2_minus) <- c("pep","logFC_br", "logCPM_br", "F_br", "PValue_br", "FDR_br")
colnames(lu1_lu_vs_li_res_lfc2_minus) <- c("pep","logFC_li", "logCPM_li", "F_li", "PValue_li", "FDR_li")
colnames(lu1_lu_vs_ki_res_lfc2_minus) <- c("pep","logFC_ki", "logCPM_ki", "F_ki", "PValue_ki", "FDR_ki")
colnames(lu1_lu_vs_mu_res_lfc2_minus) <- c("pep","logFC_mu", "logCPM_mu", "F_mu", "PValue_mu", "FDR_mu")

lu1_lu_vs_mu_res_lfc2_minus

lu1_joined_1 <- full_join(lu1_lu_vs_br_res_lfc2_minus, lu1_lu_vs_li_res_lfc2_minus,by = "pep")
lu1_joined_2 <- full_join(lu1_joined_1, lu1_lu_vs_ki_res_lfc2_minus, by = "pep")
lu1_joined_all <- full_join(lu1_joined_2, lu1_lu_vs_mu_res_lfc2_minus, by = "pep")
View(lu1_joined_all)


#create input for setupr

area_lu1_lu_vs_br <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_br != 0))
area_lu1_lu_vs_li <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_li != 0))
area_lu1_lu_vs_ki <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_ki != 0))
area_lu1_lu_vs_mu <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_mu != 0))

area_lu1_lu_vs_br_ki <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_br != 0 & lu1_joined_all$logCPM_ki != 0))
area_lu1_lu_vs_br_li <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_br != 0 & lu1_joined_all$logCPM_li != 0))
area_lu1_lu_vs_br_mu <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_br != 0 & lu1_joined_all$logCPM_mu != 0))
area_lu1_lu_vs_ki_li <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_ki != 0 & lu1_joined_all$logCPM_li != 0))
area_lu1_lu_vs_ki_mu <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_ki != 0 & lu1_joined_all$logCPM_mu != 0))
area_lu1_lu_vs_li_mu <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_li != 0 & lu1_joined_all$logCPM_mu != 0))



area_lu1_lu_vs_br_ki_li <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_br != 0 & lu1_joined_all$logCPM_ki != 0
                                       & lu1_joined_all$logCPM_li != 0))
area_lu1_lu_vs_br_li_mu <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_br != 0 & lu1_joined_all$logCPM_li != 0
                                       & lu1_joined_all$logCPM_mu != 0))
area_lu1_lu_vs_br_mu_ki <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_br != 0 & lu1_joined_all$logCPM_mu != 0
                                       & lu1_joined_all$logCPM_ki != 0))
area_lu1_lu_vs_ki_li_mu <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_ki != 0 & lu1_joined_all$logCPM_li != 0
                                       & lu1_joined_all$logCPM_mu != 0))


area_lu1_lu_vs_br_mu_ki_li <- nrow(subset(lu1_joined_all, lu1_joined_all$logFC_br != 0 & lu1_joined_all$logCPM_mu != 0
                                          & lu1_joined_all$logCPM_ki != 0 & lu1_joined_all$logCPM_li))


###UpSetR analysis 

upset_input_lu1 <- c("Brain"=area_lu1_lu_vs_br, "Liver"=area_lu1_lu_vs_li, "Kidney"=area_lu1_lu_vs_ki, "Muscle" = area_lu1_lu_vs_mu,
                     "Brain&Liver"=area_lu1_lu_vs_br_li, "Brain&Liver"=area_lu1_lu_vs_br_li, "Brain&Muscle"=area_lu1_lu_vs_br_mu,
                     "Kidney&Liver" = area_lu1_lu_vs_ki_li,"Kidney&Muscle" = area_lu1_lu_vs_ki_mu,
                     "Liver&Muscle" = area_lu1_lu_vs_li_mu,"Brain&Muscle&Kidney" = area_lu1_lu_vs_br_mu_ki,
                     "Brain&Muscle&Liver" = area_lu1_lu_vs_br_li_mu, "Brain&Kidney&Liver" = area_lu1_lu_vs_br_ki_li,
                     "Kidney&Muscle&Liver" = area_lu1_lu_vs_ki_li_mu,
                     "Brain&Muscle&Kidney&Liver"=area_lu1_lu_vs_br_mu_ki_li)

upset(fromExpression(upset_input_lu1), order.by = "freq")

