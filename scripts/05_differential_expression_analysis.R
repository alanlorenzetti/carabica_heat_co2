# alorenzetti 20211019

# description ####
# this script will take
# the salmon quant output
# files and perform
# differential expression
# analysis of carabica
# heat and co2 response

# setting wd
#setwd("~/gdrive/carabica_heat_co2/")

# creating results and plots
# directories
for(i in c("results", "plots")){
  if(!dir.exists(i)){dir.create(i)}
}

# setting thresholds
qthr = 0.01
lfcthr = 1



# finding files
PRJNA609253 = paste0(list.dirs(path = "salmon_quant/PRJNA609253", recursive = F), "/",
                     "quant.sf")
PRJNA630692 = paste0(list.dirs(path = "salmon_quant/PRJNA630692", recursive = F), "/",
                     "quant.sf")
PRJNA630692 = PRJNA630692[grepl(pattern = "/CA_", x = PRJNA630692)]

# importing quant data
# ~12k transcripts missing in tx2gene object
quant = list()

quant$PRJNA609253 = tximport(files = PRJNA609253,
                             type = "salmon",
                             tx2gene = tx2gene)

quant$PRJNA630692 = tximport(files = PRJNA630692,
                             type = "salmon",
                             tx2gene = tx2gene)

# defining colData
sample = c(str_replace(PRJNA609253, ".*-(.*?)/quant.*", "PRJNA609253_\\1"),
           str_replace(PRJNA630692, ".*_(.*)/quant.*", "PRJNA630692_\\1"))
temperature = c(str_replace(PRJNA609253, ".*-(.*?)_.*", "\\1"),
                str_replace(PRJNA630692, ".*ppm_(.*?)_.*", "\\1"))
co2 = c(rep("Unknown", length(PRJNA609253)),
        str_replace(PRJNA630692, ".*_(.*ppm?)_..C.*", "\\1"))
cultivar = c(str_replace(PRJNA609253, ".*-..T_(.*?)_.*", "\\1"),
             str_replace(PRJNA630692, ".*/(.*?)_...ppm.*", "\\1TU"))
replicate = c(str_replace(PRJNA609253, ".*_([0-9])/quant.*", "\\1"),
              str_replace(PRJNA630692, ".*_(.*)/quant.*", "\\1"))
response = c(rep("short", length(PRJNA609253)),
             rep("long", length(PRJNA630692)))

colData = list()

colData$full = tibble(sample = sample,
                      temperature = temperature,
                      co2 = co2,
                      cultivar = cultivar,
                      replicate = replicate,
                      response = response) %>% 
  mutate(temperature = case_when(temperature == "OpT" ~ "23C",
                                 temperature == "WaT" ~ "30C",
                                 TRUE ~ temperature),
         cultivar = case_when(cultivar == "CA" ~ "Catuai",
                              cultivar == "AC" ~ "Acaua",
                              cultivar == "CATU" ~ "Icatu",
                              TRUE ~ cultivar),
         replicate = case_when(replicate == "1" ~ "A",
                               replicate == "2" ~ "B",
                               replicate == "3" ~ "C",
                               replicate == "4" ~ "D",
                               replicate == "5" ~ "E",
                               TRUE ~ replicate),
         replicate = str_replace(replicate, "[0-9]{1,2}", ""))

colData$PRJNA609253 = colData$full %>% 
  filter(response == "short") %>% 
  mutate(group = paste0(cultivar, "_", temperature))

colData$PRJNA630692 = colData$full %>% 
  filter(response == "long") %>% 
  mutate(group = paste0(co2, "_", temperature))

# dds object from tximport
#
# the function below takes care
# of implementing the method called
# "original counts and offset"
# proposed by tximport paper
dds = list()

dds$PRJNA609253 = DESeqDataSetFromTximport(txi = quant$PRJNA609253,
                                           colData = colData$PRJNA609253,
                                           design = ~ group)

dds$PRJNA630692 = DESeqDataSetFromTximport(txi = quant$PRJNA630692,
                                           colData = colData$PRJNA630692,
                                           design = ~ group)

# filtering low count entries
keep = list()
keep$PRJNA609253 = rowSums(counts(dds$PRJNA609253)) >= 10
dds$PRJNA609253 = dds$PRJNA609253[keep$PRJNA609253,]

keep$PRJNA630692 = rowSums(counts(dds$PRJNA630692)) >= 10
dds$PRJNA630692 = dds$PRJNA630692[keep$PRJNA630692,]

# running DESeq2
dds$PRJNA609253 = DESeq(dds$PRJNA609253)
dds$PRJNA630692 = DESeq(dds$PRJNA630692)

# setting function to plot volcano
volcanoPlot=function(parsedRes = parsedRes, title = title){
  volcano = parsedRes %>% 
    ggplot(aes(x = log2FoldChange,
               y = -log10(padj),
               color = status,
               label = label)) +
    geom_point(show.legend = F,
               alpha = 0.75,
               shape = 16,
               stroke = 0) +
    geom_label_repel(show.legend = F,
                     color = "black",
                     alpha = 0.75,
                     max.overlaps = 30) +
    xlab("Log<sub>2</sub> (Fold Change)") +
    ylab("-Log<sub>10</sub> (Adjusted <i>p</i>-value)") +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown()) +
    scale_color_manual(values = c("Upregulated" = tab10$red,
                                  "Downregulated" = tab10$blue,
                                  "Not significant" = "darkgrey")) +
    ggtitle(title)
  
  return(volcano)
}

# setting function to generate results
genRes=function(obj = obj, var = var,
                numerator = numerator, denominator = denominator,
                functionalAnnot = functionalAnnot){
  
  title = paste(numerator, "vs.", denominator)
  
  results = results(object = obj,
                    contrast = c(var, numerator, denominator),
                    alpha = qthr,
                    lfcThreshold = lfcthr,
                    altHypothesis="greaterAbs") %>% 
    as_tibble(rownames = "gene")
  
  minpadj = results$padj[!is.na(results$padj) & results$padj != 0] %>% min()
  minpadj = minpadj * 0.01
  
  resultsParsed = results %>% 
    mutate(status = case_when(log2FoldChange >= lfcthr & padj < qthr ~ "Upregulated",
                              log2FoldChange <= -lfcthr & padj < qthr ~ "Downregulated",
                              TRUE ~ "Not significant"),
           padj = case_when(padj == 0 ~ minpadj,
                            TRUE ~ padj))
  
  topup = resultsParsed %>%
    filter(status == "Upregulated") %>% 
    arrange(desc(log2FoldChange)) %>% 
    head(10) %>% 
    pull(gene)
  
  topdown = resultsParsed %>%
    filter(status == "Downregulated") %>% 
    arrange(log2FoldChange) %>% 
    head(10) %>% 
    pull(gene)
  
  resultsParsed = resultsParsed %>% 
    mutate(label = case_when(gene %in% c(topup, topdown) ~ gene,
                             TRUE ~ NA_character_)) %>% 
    left_join(x = .,
              y = functionalAnnot,
              by = c("gene" = "gene_id"))
  
  sigGenes = resultsParsed %>% 
    filter(status == "Upregulated" | status == "Downregulated") %>% 
    pull(gene)
  
  p = volcanoPlot(parsedRes = resultsParsed, title = title)
  
  enrichmentAnalysis = runTopGO(parsedDE = resultsParsed,
                                thr = qthr,
                                title = title)
  
  resultsList = list(results,
                     resultsParsed,
                     p,
                     sigGenes,
                     enrichmentAnalysis$BP,
                     enrichmentAnalysis$MF)
  
  names(resultsList) = c("original",
                         "parsed",
                         "volcano",
                         "sigGenes",
                         "enrichment_BP",
                         "enrichment_MF")
    
  return(resultsList)
}

# getting results for PRJNA609253
res = list()

comps = list(c("Acaua_30C", "Acaua_23C"),
             c("Catuai_30C", "Catuai_23C"),
             c("Catuai_30C", "Acaua_30C"),
             c("Catuai_23C", "Acaua_23C"))

for(i in 1:length(comps)){
  contrastName = paste0(comps[[i]][1], "_vs_", comps[[i]][2])
  res[["PRJNA609253"]][[contrastName]] = genRes(obj = dds$PRJNA609253,
                                                var = "group",
                                                numerator = comps[[i]][1],
                                                denominator = comps[[i]][2],
                                                functionalAnnot = ob)
}

# getting results for PRJNA630692
comps = list(c("380ppm_37C", "380ppm_25C"),
             c("380ppm_42C", "380ppm_25C"),
             c("700ppm_37C", "700ppm_25C"),
             c("700ppm_42C", "700ppm_25C"),
             c("700ppm_25C", "380ppm_25C"),
             c("700ppm_37C", "380ppm_37C"),
             c("700ppm_42C", "380ppm_42C"),
             c("700ppm_42C", "380ppm_25C"),
             c("700ppm_37C", "380ppm_25C"))

for(i in 1:length(comps)){
  contrastName = paste0(comps[[i]][1], "_vs_", comps[[i]][2])
  res[["PRJNA630692"]][[contrastName]] = genRes(obj = dds$PRJNA630692,
                                                var = "group",
                                                numerator = comps[[i]][1],
                                                denominator = comps[[i]][2],
                                                functionalAnnot = ob)
}

# function to write deseq2 results
writeResults = function(resultsObj = resultsObj, var = var, contrastName = contrastName){
  
  filenamexlsx = paste0("results/deseq2-", var, "-", contrastName, ".xlsx")
  filenamevolcano = paste0("plots/volcano-", var, "-", contrastName, ".png")
  
  outputxlsx = list("original" = resultsObj[[var]][[contrastName]][["original"]],
                    "parsed" = resultsObj[[var]][[contrastName]][["parsed"]],
                    "go_enrichment_BP" = resultsObj[[var]][[contrastName]][["enrichment_BP"]],
                    "go_enrichment_MF" = resultsObj[[var]][[contrastName]][["enrichment_MF"]])
  
  write.xlsx(outputxlsx,
             file = filenamexlsx,
             overwrite = T)
  
  ggsave(plot = resultsObj[[var]][[contrastName]][["volcano"]],
         filename = filenamevolcano,
         device = "png",
         width = 10,
         height = 7.5,
         dpi = 300,
         units = "in")
}

# writing results
for(i in names(res)){
  for(j in names(res[[i]])){
    writeResults(resultsObj = res,
                 var = i,
                 contrastName = j)
  }
}

# writing quant datasets to Rdata object
save(quant,
     file ="results/quant.RData")

# writing additional files
for(i in names(quant)){
  # writing sample info
  write_tsv(x = colData[[i]],
            file = paste0("results/sampleInfo_", i, ".tsv"))
  
  # writing userfriendly abundance data
  # generated by tximport function
  countsdf = assay(dds[[i]]) %>% as_tibble(rownames = "gene")
  colnames(countsdf) = c("gene", colData[[i]][["sample"]])
  
  write_tsv(x = countsdf,
            file = paste0("results/counts_", i, ".tsv"))
}

# venn diagrams ####
# creating venn diagram for PRJNA609253
vennplot = venn(combinations = list("Acauã 30C vs. Acauã 23C" = res$PRJNA609253$Acaua_30C_vs_Acaua_23C$sigGenes,
                                    "Catuaí 30C vs. Catuaí 23C" = res$PRJNA609253$Catuai_30C_vs_Catuai_23C$sigGenes)) %>%
  plot()

ggsave(plot = vennplot,
       filename = "plots/vennplot-PRJNA609253_temperature_cultivar.png",
       device = "png",
       width = 7,
       height = 7,
       dpi = 300,
       units = "in")

# creating venn diagram for PRJNA630692
# temperature
vennplot = venn(combinations = list("380ppm 37C vs. 380ppm 25C" = res$PRJNA630692$`380ppm_37C_vs_380ppm_25C`$sigGenes,
                                    "380ppm 42C vs. 380ppm 25C" = res$PRJNA630692$`380ppm_42C_vs_380ppm_25C`$sigGenes,
                                    "700ppm 37C vs. 700ppm 25C" = res$PRJNA630692$`700ppm_37C_vs_700ppm_25C`$sigGenes,
                                    "700ppm 42C vs. 700ppm 25C" = res$PRJNA630692$`700ppm_42C_vs_700ppm_25C`$sigGenes)) %>%
  plot()

ggsave(plot = vennplot,
       filename = "plots/vennplot-PRJNA630692_temperature.png",
       device = "png",
       width = 10,
       height = 8,
       dpi = 300,
       units = "in")

# co2
vennplot = venn(combinations = list("700ppm 25C vs. 380ppm 25C" = res$PRJNA630692$`700ppm_25C_vs_380ppm_25C`$sigGenes,
                                    "700ppm 37C vs. 380ppm 37C" = res$PRJNA630692$`700ppm_37C_vs_380ppm_37C`$sigGenes,
                                    "700ppm 42C vs. 380ppm 42C" = res$PRJNA630692$`700ppm_42C_vs_380ppm_42C`$sigGenes)) %>%
  plot()

ggsave(plot = vennplot,
       filename = "plots/vennplot-PRJNA630692_co2.png",
       device = "png",
       width = 8,
       height = 8,
       dpi = 300,
       units = "in")

# co2 and temperature
vennplot = venn(combinations = list("700ppm 37C vs. 380ppm 25C" = res$PRJNA630692$`700ppm_37C_vs_380ppm_25C`$sigGenes,
                                    "700ppm 42C vs. 380ppm 25C" = res$PRJNA630692$`700ppm_42C_vs_380ppm_25C`$sigGenes)) %>%
  plot()

ggsave(plot = vennplot,
       filename = "plots/vennplot-PRJNA630692_temperature_co2.png",
       device = "png",
       width = 8,
       height = 8,
       dpi = 300,
       units = "in")

# exploring data with PCA
# vsd = vst(dds$PRJNA630692, blind=FALSE)
# pcaData = plotPCA(vsd, intgroup=c("temperature", "cultivar", "co2"), returnData=T)
# percentVar = round(100 * attr(pcaData, "percentVar"))
# ggplot(pcaData, aes(PC1, PC2, color=temperature, shape=cultivar, size=co2)) +
#   geom_point() +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#   coord_fixed() +
#   theme_bw()

