install.packages("renv")
renv::restore()

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(enrichR)
library(msigdbr)
library(clusterProfiler)
library(apeglm)
library(EnhancedVolcano)
library(DEGreport)
library(RColorBrewer)
library(pathview)
library(gage)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(gageData)


rse <- readRDS("EwS.rds")

dds <- DESeqDataSet(rse, design = ~ condition)

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]

dds <- DESeq(dds)

result <- results(dds)

#Requirement 1 
resultsNames(dds)

resultLFC <- lfcShrink(dds, coef = "condition_shEF1_vs_shCTR", type = "apeglm")


plotCounts(dds, gene=which.min(resultLFC$padj), intgroups = "conditions")

#Requirement 2

plotMA(resultLFC, ylim=c(-2,2))

#Requirement 3
resultordered = result[order(result$pvalue),]

filteredresult = subset(resultordered, padj < 0.1)

write.csv(as.data.frame(filteredresult), file= "shEF1_vs_shCTR.csv")

#Requirement 4

EnhancedVolcano(result, lab = rownames(result), x = "log2FoldChange", y="padj")

#Requirement 5

dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)

result_tb <- result %>% data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

padj.cutoff <- 0.05
lfc.cutoff <- 0.58

sigOE <- result_tb %>% 
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)


normalized_counts <- normalized_counts %>% 
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()



top10_sigOE_genes <- result_tb %>%
  arrange(padj) %>%
  pull(gene) %>%
  head(n=10)


bot10_sigOE_genes <- result_tb %>%
  arrange(padj) %>%
  pull(gene) %>%
  tail(n=10)

top10_sigOE_norm <- normalized_counts %>% 
  filter(gene %in% top10_sigOE_genes)

bot10_sigOE_norm <- normalized_counts %>% 
  filter(gene %in% bot10_sigOE_genes)

heat_colors <- brewer.pal(6,"YlOrRd")

top_norm_OEsig <- normalized_counts[,c(1:8)] %>%
  filter(gene %in% top10_sigOE_norm$gene) %>%
  data.frame() %>% 
  column_to_rownames(var="gene")

bot_norm_OEsig <- normalized_counts[,c(1:8)] %>%
  filter(gene %in% bot10_sigOE_norm$gene) %>%
  data.frame() %>% 
  column_to_rownames(var="gene")

heat_norm_sig <- rbind(top_norm_OEsig, bot_norm_OEsig)

pheatmap(heat_norm_sig, color = heat_colors,
         cluster_rows = T, show_rownames = T,
         border_color = NA, fontsize = 10,
         scale = "row", fontsize_row = 10, 
         height = 20) 


#Requirement 6 (Left unfulfilled due to Issue with the mapID function)

kg.hsa <- kegg.gsets(species = "hsa")
kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]

go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]


KEGGresult_top <- sigOE %>% 
  arrange(padj) %>%
  head(10)

KEGGresult_bot <- sigOE %>% 
  arrange(padj) %>%
  tail(10)


KEGGresult_combined <- rbind(KEGGresult_top, KEGGresult_bot)

KEGGresult_combined$refseq <- gsub("\\..*","", KEGGresult_combined$gene)

Testresult <- result

Testresult$refseq <- gsub("\\..*","", row.names(Testresult))

Testresult$symbol <- mapIds(org.Hs.eg.db, keys = Testresult$refseq,
                                     column = "SYMBOL", keytype = "ENSEMBL", multival = "first")

#Everything above works correctly. 
#
#
#NA occurs here. 
KEGGresult_combined$symbol <- mapIds(org.Hs.eg.db, keys = KEGGresult_combined$refseq,
                        column = "SYMBOL", keytype = "ENSEMBL", multival = "first")
#NA occurs here. 
KEGGresult_combined$entrez <- mapIds(org.Hs.eg.db, keys = KEGGresult_combined$refseq,
                            column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
#NA occurs here. 
KEGGresult_combined$name <- mapIds(org.Hs.eg.db, keys = KEGGresult_combined$refseq,
                            column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")

KEGGresult.fc <- KEGGresult_combined$log2FoldChange

names(KEGGresult.fc) <- KEGGresult_combined$entrez

#All these results return NA after analysis. 
fc.kegg.sigmet.p <- gage(KEGGresult.fc, gsets = kegg.sigmet.gs)
fc.kegg.dise.p <- gage(KEGGresult.fc, gsets = kegg.dise.gs)
fc.go.bp.p <- gage(KEGGresult.fc, gsets = go.bp.gs)
fc.go.mf.p <- gage(KEGGresult.fc, gsets = go.mf.gs)
fc.go.cc.p <- gage(KEGGresult.fc, gsets = go.cc.gs)

fc.kegg.sigmet.p.up <- as.data.frame(fc.kegg.sigmet.p$greater)
fc.kegg.dise.p.up <- as.data.frame(fc.kegg.dise.p$greater)

fc.kegg.sigmet.p.down <- as.data.frame(fc.kegg.sigmet.p$less)
fc.kegg.dise.p.down <- as.data.frame(fc.kegg.dise.p$less)

fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
fc.go.mf.p.up <- as.data.frame(fc.go.mf.p$greater)
fc.go.cc.p.up <- as.data.frame(fc.go.cc.p$greater)

fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
fc.go.mf.p.up <- as.data.frame(fc.go.mf.p$greater)
fc.go.cc.p.up <- as.data.frame(fc.go.cc.p$greater)

#Visualize the pathway, I assume from KEGG that hsa05202 is the pathway we want to investigate since it uses the EWSR1 gene.
fc.kegg.sigmet.p.up[grepl("hsa05202", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]

