# BiocManager::install(c("affy", "limma","arrayQualityMetrics","affycoretools"))
# install.packages("rio")
library(affy)
library(limma)
library(rio)
install.packages("tidyverse")
library(tidyverse)
library(arrayQualityMetrics)
library(affycoretools)


datadir <- "data"
resultsdir <- "results"
plotsdir <- "plots"
scriptsdir <- "src"
accession <- "E-MTAB-1862"

# data preprocessing ----
# construct affy batch
affy_batch <- ReadAffy(filenames=list.celfiles(path="data/E-MTAB-1862/", full.names=TRUE))
#str(affy_batch)
#colnames(affy_batch)
#[1] "BN1.CEL"  "BN10.CEL" "BN11.CEL" "BN2.CEL" etc
#rownames(affy_batch)

# load metadata
#rm(pdata)
pdata <- rio::import("data/E-MTAB-1862/E-MTAB-1862_metadata.sdrf.txt")
colnames(pdata)

# check if all samples from pdata are represented in the CEL files
rownames(pdata) <- pdata$`Array Data File`
sample_names_abatch <- sampleNames(affy_batch)
sample_names_pdata <- pdata$`Array Data File`

length(intersect(sample_names_abatch, sample_names_pdata)) == length(sample_names_abatch)

# combine pdata and CEL data and save as RDS
pdata <- pdata[sample_names_abatch,] #????
pData(affy_batch) <- pdata

saveRDS(affy_batch, file="data/abatch_E-MTAB-1862.RDS")

#AFFYBATCH DATA ARE IN THE LINEAR SPACE SO some paramter log 2 needs to be set to 2
#THINK OF THE LINEAR MODEL YOU WANT TO USE, LOOK IT UP IN THE LIMMA VIGENTTE
#READ LITERATURE ON THE DESEASE WE ARE DEALING WITH

# generate quality control QC ----
#BiocManager::install('arrayQualityMetrics')
#library(arrayQualityMetrics)
abatch <- rio::import(here::here(datadir, "abatch_E-MTAB-1862.RDS"))

QC <- arrayQualityMetrics(abatch,
                          do.logtransform= TRUE,
                          outdir = here::here(resultsdir, "QC_abatch"))

outliers <- c("BN5.CEL","PM2.CEL") #add IDs here
selected <- !(sampleNames(abatch) %in% outliers)

abatch_outlierremoved <- abatch[,selected]
abatch_outlierremoved %>%
  rio::export(here::here(datadir, accession, paste0(accession, "_abatch_outlierremoved.RDS")))

# normalization of the data with rma() ----

abatch_outlierremoved <- rio::import(here::here(datadir, accession, paste0(accession,"_abatch_outlierremoved.RDS")))

eset <- rma(abatch_outlierremoved)

eset %>%
  rio::export(here::here(datadir, accession, paste0(accession, "_eset.RDS")))


eset <- rio::import(here::here(datadir, accession, paste0(accession, "_eset.RDS")))
                    
QC <- arrayQualityMetrics(eset,
                          do.logtransform = FALSE,
                          outdir = here::here(resultsdir, "QC_eset"))

# remove outliers ----
outliers <- c()

outliers<- QC$modules %>%
  purrr::map(function(QC){
    QC@outliers@which %>% names
  }) %>%
  unlist %>%
  unique
outliers

selected <- !(sampleNames(eset) %in% outliers)

eset_outlierremoved <- eset[,selected]

eset_outlierremoved %>%
  rio::export(here::here(datadir, accession, paste0(accession, "_eset_outlierremoved.RDS")))


# cannot run the function annotateEset???
# eset_outlierremoved <- annotateEset(eset_outlierremoved,
#                                    annotation(eset_outlierremoved),
#                                    type = "core")

# manual feature annotation ----
database <- hgu133plus2.db::hgu133plus2.db
featureAnnotation <- AnnotationDbi::select(database,
                                           columns = c("PROBEID", "SYMBOL", "ENTREZID", "ENSEMBL", "GENENAME"),
                                           keys = AnnotationDbi::keys(database, keytype = "PROBEID"),
                                           keytype = "PROBEID") %>%
  group_by(PROBEID) %>%
  summarize_all(function(x){paste(unique(x), collapse = ", ")}) %>%
  mutate(rowname = PROBEID) %>%
  column_to_rownames()

fData(eset_outlierremoved) <- featureAnnotation[featureNames(eset_outlierremoved),]


# export annotated expression set
eset_outlierremoved %>%
  rio::export(here::here(datadir, accession, paste0(accession, "_eset_outlierremoved.RDS")))

# clean data from NAs
selected_features <- !is.na(fData(eset_outlierremoved)$SYMBOL)

eset_limma <- eset_outlierremoved[selected_features,]
#View(eset_limma)

eset_limma$Label_is_benign <- as.factor(eset_limma$Label_is_benign)

design <- model.matrix(~0 + Label_is_benign,
                       data = eset_limma)

# this could be more experiments, but I have one 1 comparison
contrasts <- list(bening_vs_malignant = "Label_is_benign1 - Label_is_benign0")

contrast_matrix <- makeContrasts(contrasts = contrasts,
                                 levels = design)

colnames(contrast_matrix) <- names(contrasts)


# fit the linear model ----
fit <- lmFit(eset_limma, design) %>%
  contrasts.fit(contrast_matrix) %>%
  eBayes

coefs <- colnames(contrast_matrix) %>%
  set_names()


results <- coefs %>%
  purrr::map(function(coef, fit){
    topTable(fit,
             coef = coefs,
             number = Inf)  
  },
  fit = fit)
#View(results)

#export the results
results %>%
  rio::export(here::here(resultsdir, 'limma_result.xlsx'))

#generate volcano plots ----
#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

# The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.
EnhancedVolcano(results$bening_vs_malignant,
                lab = rownames(results$bening_vs_malignant),
                x = 'logFC',
                y = 'P.Value')

# term enrichment ----
BiocManager::install('clusterProfiler')
library(clusterProfiler)
BiocManager::install('msigdbr')
library(msigdbr)

msigdb <- msigdbr(species = "Homo sapiens")

term2gene <- msigdb %>%
  dplyr::select(gs_id, gene_symbol)

term2name <- msigdb %>%
  dplyr::select(gs_id, gs_name)

#Term enrichment analysis
TEA_result <- results %>%
  purrr::map(function(results, term2gene, term2name){
    gene <- results %>%
      dplyr::filter(adj.P.Val < 0.05) %>% #nut pvalues unter 0.05
      dplyr::select(SYMBOL) %>%
      deframe #wandelt einen datafraem in eien vector, den wir in die enricher function
    enricher(gene,
             pvalueCutoff = 1,
             qvalueCutoff = 1,
             TERM2GENE = term2gene, 
             TERM2NAME = term2name)
  },
  term2gene = term2gene, 
  term2name = term2name)

TEA_result %>%
  purrr::map(data.frame) %>%
  rio::export(here::here(resultsdir, "TEA_result.xlsx"))

# visualization ----
# https://yulab-smu.top/biomedical-knowledge-mining-book/index.html

library(enrichplot)

barplot <- TEA_result %>%
  purrr::map(function(TEA_result){
    barplot(TEA_result)
  })

# BAR plot
barplot %>%
  purrr::map2(names(.),
              function(plot, name){
                filename <- here::here(plotsdir, paste0(name, "_barplot.pdf"))
                ggsave(filename = filename,
                       plot = plot)
              })

# DOT plot
dotplot <- TEA_result %>%
  purrr::map(function(TEA_result){
    dotplot(TEA_result)
  })

dotplot %>%
  purrr::map2(names(.),
              function(plot, name){
                filename <- here::here(plotsdir, paste0(name, "_dotplot.pdf"))
                ggsave(filename = filename,
                       plot = plot)
              })
