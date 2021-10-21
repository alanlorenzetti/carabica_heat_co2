# alorenzetti 20211019

# description ####
# this script will take
# the salmon quant output
# files and perform
# differential expression
# analysis of carabica
# heat and co2 response

# setting wd
setwd("~/gdrive/carabica_heat_co2/")

# loading libs ####
library(pacman)

packs = c("tidyverse",
          "DESeq2",
          "tximport",
          "rtracklayer",
          "ggtext",
          "ggthemes",
          "ggrepel",
          "openxlsx",
          "eulerr")

p_load(char = packs)

# setting ggplot theme
theme_set(theme_bw())

# getting colors from tab10 color scheme
tab10 = list()
tab10$blue = ggthemes::ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[1]
tab10$red = ggthemes::ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[3]

# creating results and plots
# directories
for(i in c("results", "plots")){
  if(!dir.exists(i)){dir.create(i)}
}

# setting thresholds
qthr = 0.01
lfcthr = 1

# getting started ####
# creating a tx2gene object
annot = rtracklayer::import(con = "data/Carabica.gff") %>% 
  as.data.frame() %>% 
  as_tibble()

tx2gene = annot %>% 
  filter(type == "mRNA") %>% 
  select(TXNAME = Name, GENEID = gene)

# finding files
PRJNA609253 = paste0(list.dirs(path = "salmon_quant/PRJNA609253", recursive = F), "/",
                     "quant.sf")
PRJNA630692 = paste0(list.dirs(path = "salmon_quant/PRJNA630692", recursive = F), "/",
                     "quant.sf")
PRJNA630692 = PRJNA630692[grepl(pattern = "/CA_", x = PRJNA630692)]

files = c(PRJNA609253,
          PRJNA630692)

# importing quant data
# ~12k transcripts missing in tx2gene object
quant = tximport(files = files,
                 type = "salmon",
                 tx2gene = tx2gene)

# defining colData
sample = c(str_replace(PRJNA609253, ".*-(.*?)/quant.*", "PRJNA609253_\\1"),
           str_replace(PRJNA630692, ".*_(.*)/quant.*", "PRJNA630692_\\1"))
temperature = c(str_replace(PRJNA609253, ".*-(.*?)_.*", "\\1"),
                str_replace(PRJNA630692, ".*ppm_(.*?)_.*", "\\1"))
co2 = c(rep("380ppm", length(PRJNA609253)),
        str_replace(PRJNA630692, ".*_(.*ppm?)_..C.*", "\\1"))
cultivar = c(str_replace(PRJNA609253, ".*-..T_(.*?)_.*", "\\1"),
             str_replace(PRJNA630692, ".*/(.*?)_...ppm.*", "\\1TU"))
replicate = c(str_replace(PRJNA609253, ".*_([0-9])/quant.*", "\\1"),
              str_replace(PRJNA630692, ".*_(.*)/quant.*", "\\1"))

colData = tibble(sample = sample,
                 temperature = temperature,
                 co2 = co2,
                 cultivar = cultivar,
                 replicate = replicate) %>% 
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
         replicate = str_replace(replicate, "[0-9]{1,2}", ""),
         group = paste0(cultivar, "_", temperature, "_", co2))

# dds object from tximport
# according to
#
# the function below takes care
# of implementing the method called
# "original counts and offset"
# proposed by tximport paper
dds = DESeqDataSetFromTximport(txi = quant,
                               colData = colData,
                               design = ~ co2 + temperature)

# filtering low count entries
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]

# running DESeq2
dds = DESeq(dds)

# exploring data with PCA
vsd = vst(dds, blind=FALSE)
pcaData = plotPCA(vsd, intgroup=c("temperature", "cultivar", "co2"), returnData=T)
percentVar = round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=temperature, shape=cultivar, size=co2)) +
  geom_point() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()

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
                     alpha = 0.75) +
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
genRes=function(obj = obj, var = var, numerator = numerator, denominator = denominator){
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
                             TRUE ~ NA_character_))
  
  sigGenes = resultsParsed %>% 
    filter(status == "Upregulated" | status == "Downregulated") %>% 
    pull(gene)
  
  p = volcanoPlot(parsedRes = resultsParsed, title = title)
  
  resultsList = list(results,
                     resultsParsed,
                     p,
                     sigGenes)
  
  names(resultsList) = c("original",
                         "parsed",
                         "volcano",
                         "sigGenes")
    
  return(resultsList)
}

# getting results for temperature
res = list()

comps = list(c("25C", "23C"),
             c("30C", "25C"),
             c("37C", "30C"),
             c("42C", "37C"),
             c("30C", "23C"),
             c("37C", "25C"),
             c("42C", "30C"))

for(i in 1:length(comps)){
  contrastName = paste0(comps[[i]][1], "_vs_", comps[[i]][2])
  res[["temperature"]][[contrastName]] = genRes(obj = dds,
                                               var = "temperature",
                                               numerator = comps[[i]][1],
                                               denominator = comps[[i]][2])
}

# getting results for co2
comps = list(c("700ppm", "380ppm"))

for(i in 1:length(comps)){
  contrastName = paste0(comps[[i]][1], "_vs_", comps[[i]][2])
  res[["co2"]][[contrastName]] = genRes(obj = dds,
                                        var = "co2",
                                        numerator = comps[[i]][1],
                                        denominator = comps[[i]][2])
}

# function to write deseq2 results
writeResults = function(resultsObj = resultsObj, var = var, contrastName = contrastName){
  
  filenamexlsx = paste0("results/deseq2-", var, "-", contrastName, ".xlsx")
  filenamevolcano = paste0("plots/deseq2-", var, "-", contrastName, ".png")
  
  outputxlsx = list("original" = resultsObj[[var]][[contrastName]][["original"]],
                    "parsed" = resultsObj[[var]][[contrastName]][["parsed"]])
  
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

# writing quant dataset to Rdata object
save(quant,
     file = "results/quant.RData")

# writing sample info
write_tsv(x = colData,
          file = "results/sampleInfo.tsv")

# writing userfriendly abundance data
# generated by tximport function
countsdf = assay(dds) %>% as_tibble(rownames = "gene")
colnames(countsdf) = c("gene", colData$sample)

write_tsv(x = countsdf,
          file = "results/counts.tsv")

# creating venn diagram of consecutive
# temperature changes
# drawing venn diagrams
vennplot = venn(combinations = list("25C_vs_23C" = res$temperature$`25C_vs_23C`$sigGenes,
                                    "30C_vs_25C" = res$temperature$`30C_vs_25C`$sigGenes,
                                    "37C_vs_30C" = res$temperature$`37C_vs_30C`$sigGenes,
                                    "42C_vs_37C" = res$temperature$`42C_vs_37C`$sigGenes)) %>%
  plot()

ggsave(plot = vennplot,
       filename = "plots/vennplot_temperature.png",
       device = "png",
       width = 7,
       height = 7,
       dpi = 300,
       units = "in")
