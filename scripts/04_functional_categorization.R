# alorenzetti 20211116

# description ####
# this script will take as input
# the funcitonal categorization performed
# using OmicsBox for the Coffea arabica transcriptome
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/713/225/GCF_003713225.1_Cara_1.0/GCF_003713225.1_Cara_1.0_rna.fna.gz

# getting started ####
# loading file
omicsbox = read_tsv(file = "data/omicsbox_table-final2.txt") %>% 
  dplyr::select(-`...1`)

# creating a tx2gene object
annot = rtracklayer::import(con = "data/Carabica.gff") %>% 
  as.data.frame() %>% 
  as_tibble()

# considering most of entries (except rRNA and partial mRNA)
# in RefSeq transcriptome file
tx2gene = annot %>%
  filter(gbkey == "mRNA" | gbkey == "ncRNA" | gbkey == "misc_RNA") %>%
  dplyr::select(TXNAME = Name, GENEID = gene) %>%
  drop_na() %>%
  unique()

# getting NCBI gene products
ncbiproducts = annot %>%
  filter(gbkey == "mRNA" | gbkey == "ncRNA" | gbkey == "misc_RNA") %>%
  dplyr::select(TXNAME = Name,
                GENEID = gene,
                ncbi_product = product) %>%
  mutate(ncbi_product = str_replace(ncbi_product, ", transcript var.*$", "")) %>% 
  drop_na() %>%
  unique() %>% 
  group_by(GENEID) %>% 
  summarise(ncbi_product = paste0(c(ncbi_product) %>% unique(), collapse = ";"))

# parsing omicsbox data
ob = omicsbox %>% 
  dplyr::select(SeqName,
                Description,
                `GO IDs`,
                `GO Names`,
                `Enzyme Codes`,
                `Enzyme Names`,
                `InterPro IDs`)

names(ob) = c("SeqName",
              "description",
              "go_id",
              "go_names",
              "enzyme_id",
              "enzyme_names",
              "interpro_id")

# declaring functions to collapse
# redundant values of columns below
collapseNaRm = function(x){
  colnarm = x %>%
    .[!is.na(.)] %>%
    unique() %>%
    paste0(collapse = "; ")
  
  if(colnarm == ""){colnarm = NA}
  return(colnarm)
}

remRed = function(x){
  nonRed = x %>%
    str_split(string = .,
              pattern = "; ",
              simplify = T) %>%
    as.character() %>% 
    .[!is.na(.)] %>%
    unique() %>% 
    paste0(collapse = ";")
  
  if(nonRed == ""){nonRed = NA}
  return(nonRed)
}

ob = ob %>% 
  mutate(description = case_when(str_detect(description, "---NA---") ~ NA_character_,
                                 TRUE ~ description),
         interpro_id = case_when(str_detect(interpro_id, "no IPS match") ~ NA_character_,
                                 TRUE ~ interpro_id),
         description = str_replace(description, " isoform.*$", "")) %>% 
  left_join(x = .,
            y = tx2gene,
            by = c("SeqName" = "TXNAME")) %>% 
  filter(!is.na(GENEID)) %>% 
  group_by(GENEID) %>% 
  summarise(across(.cols = everything(),
                   .fns = ~ collapseNaRm(.x))) %>% 
  ungroup() %>% 
  dplyr::select(-SeqName) %>%
  rowwise() %>% 
  mutate(across(.cols = -GENEID,
                .fns = ~ remRed(.x))) %>% 
  dplyr::rename(gene_id = "GENEID",
                omicsbox_product = "description")

# adding NCBI gene products
ob = left_join(x = ob,
               y = ncbiproducts,
               by = c("gene_id" = "GENEID")) %>% 
  dplyr::select(gene_id,
                omicsbox_product,
                ncbi_product,
                go_id,
                go_names,
                enzyme_id,
                enzyme_names,
                interpro_id)

# creating a table to store GO terms per gene
# and a list GO terms per genes
# to be used by topGO
GO = list()

GO$fullTable = omicsbox %>%
  dplyr::select(SeqName, `GO IDs`) %>%
  left_join(x = .,
            y = tx2gene,
            by = c("SeqName" = "TXNAME")) %>% 
  separate_rows(`GO IDs`, sep = "; ") %>%
  distinct() %>% 
  mutate(go_class = case_when(str_detect(string = `GO IDs`, pattern = "^P") ~ "biological_process",
                              str_detect(string = `GO IDs`, pattern = "^F") ~ "molecular_function"),
         transcript_id = SeqName,
         go_id = str_replace(`GO IDs`, "^.:", "")) %>% 
  dplyr::select(transcript_id, go_id, go_class, gene_id = GENEID)

GO$MF$tib = GO$fullTable %>% 
  filter(go_class == "molecular_function") %>% 
  group_by(gene_id) %>% 
  summarise(go_id = paste0(c(go_id) %>% unique(), collapse = ";")) %>% 
  drop_na()
  
GO$BP$tib = GO$fullTable %>% 
  filter(go_class == "biological_process") %>% 
  group_by(gene_id) %>% 
  summarise(go_id = paste0(c(go_id) %>% unique(), collapse = ";")) %>% 
  drop_na()

GO$MF$list = str_split(string = GO$MF$tib$go_id, pattern = ";")
names(GO$MF$list) = GO$MF$tib$gene_id

GO$BP$list = str_split(string = GO$BP$tib$go_id, pattern = ";")
names(GO$BP$list) = GO$BP$tib$gene_id

# creating function to retrieve enriched terms table
# to be applied on the next script
# testing gene enrichment analysis
runTopGO = function(parsedDE = parsedDE,
                    thr = thr,
                    title = title){
  
  title = title %>% 
    str_replace(" vs. ", "_vs_")
  
  enrich = list()
  final = list()
  
  nms = parsedDE$gene
  qvl = parsedDE %>% 
    mutate(padj = case_when(is.na(padj) ~ 1,
                            TRUE ~ padj)) %>% 
    pull(padj)
  
  names(qvl) = nms
  
  topDiffGenes = function(allScore) {
    return(allScore < thr)
  }
  
  for(i in c("BP", "MF")){
    topgoobj = new("topGOdata",
                   description = title,
                   ontology = i,
                   allGenes = qvl,
                   geneSel = topDiffGenes,
                   nodeSize = 5,
                   gene2GO = GO[[i]][["list"]],
                   annotationFun = annFUN.gene2GO)
    
    fisher = new("classicCount",
                 testStatistic = GOFisherTest,
                 name = "Fisher test")
    
    ks = new("classicScore",
             testStatistic =  GOKSTest,
             name = "KS tests")
    
    enrich$fisher = getSigGroups(topgoobj, fisher)
    enrich$ks = getSigGroups(topgoobj, ks)
    
    final[[i]] = GenTable(object = topgoobj,
                          fisher = enrich$fisher,
                          KS = enrich$ks,
                          orderBy = "fisher",
                          topNodes = 100)
    
    printGraph(topgoobj,
               enrich$fisher,
               firstSigNodes = 10,
               fn.prefix = paste0("plots/", title, "_enrichFisher_", i),
               useInfo = "all",
               pdfSW = TRUE)
  }
  
  return(final)
}
