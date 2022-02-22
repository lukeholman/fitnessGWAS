# Needs the following packages: clusterProfiler, fgsea, tidyverse

# Function to run GSEA implemented in the fgsea package
# df needs to have column called FBID


# Get the PANTHER 'GO slim' (a custom sub-set of the full GO ontology; increases power by not running too many tests):
panther_slim <- readLines("http://data.pantherdb.org/PANTHER16.0/ontology/PANTHERGOslim.obo")
panther_slim <- panther_slim[substr(panther_slim, 1,6) == "id: GO"] %>%
  str_remove_all("id[:] ")

GO_and_KEGG_gsea <- function(df, column, min.size = 50, keep.all = FALSE){

  df <- as.data.frame(df)
  names(df)[names(df) == column] <- "focal_metric"

  df <- df %>% arrange(focal_metric)

  p <- 0.05; if(keep.all) p <- 1 # Set the significance threshold

  # Set up the geneList object in the form needed by the fgsea function - named, ranked vector of pheromone sensitivity per gene
  geneList <- df$focal_metric
  names(geneList) <- df$FBID
  gene_universe <- names(geneList)

  # Internal function to run GO enrichment
  GO_enrichment <- function(geneList, ontol){

    pathways <- tbl(db, "GO") %>%
      filter(go_id %in% panther_slim) %>%  # reduce to only terms in the PANTHER slim
      filter(FBID %in% gene_universe,
             ontology == ontol) %>%
      select(-ontology) %>% collect(n=Inf)

    pathways <- with(pathways, split(FBID, go_id))

    result <- fgsea::fgsea(pathways, geneList, nperm = 10000,
                           minSize = min.size, maxSize = 500) %>%
      filter(pval <= p)

    collapse_pathways <- fgsea::collapsePathways(result, pathways, geneList)
    pathways_to_keep <- c(collapse_pathways[[1]], names(collapse_pathways[[2]]))
    result <- result %>% filter(pathway %in% pathways_to_keep)

    result <- result %>%
      left_join(tbl(db, "go_meanings") %>%
                  select(-ontology) %>% collect(n=Inf), by = c("pathway" = "GO"))

    if(nrow(result) == 0) return(NULL)
    if(ontol == "BP") Test_type <- "GO: Biological process"
    if(ontol == "MF") Test_type <- "GO: Molecular function"
    if(ontol == "CC") Test_type <- "GO: Cellular component"
    data.frame(Test_type = Test_type,
               result %>% arrange(padj),
               stringsAsFactors = FALSE)
  }

  # Internal function to run KEGG enrichment
  kegg_enrichment <- function(geneList){

    pathways <- tbl(db, "KEGG") %>%
      filter(FBID %in% gene_universe) %>%
      collect(n = Inf)
    pathways <- with(pathways, split(FBID, kegg_id))
    result <- fgsea::fgsea(pathways, geneList, nperm = 10000, minSize = min.size, maxSize = 500) %>%
      filter(pval <= p)

    dmel_kegg <- clusterProfiler::download_KEGG("dme")$KEGGPATHID2NAME

    result <- result %>%
      left_join(dmel_kegg %>% rename(term = to) %>% mutate(from = gsub("dme", "", from)), by = c("pathway" = "from")) %>%
      mutate(pathway = str_replace_all(pathway, "dme", "KEGG:"))

    if(nrow(result) == 0) return(NULL)
    data.frame(Test_type = "KEGG",
               result %>% arrange(pval),
               stringsAsFactors = FALSE)
  }

  rbind(GO_enrichment(geneList, "BP"),
        GO_enrichment(geneList, "MF"),
        GO_enrichment(geneList, "CC"),
        kegg_enrichment(geneList)) %>%
    as_tibble() %>%
    mutate(genes = map_chr(leadingEdge, ~ tbl(db, "genes") %>% # Add column of gene names
                             filter(FBID %in% .x) %>%
                             pull(gene_name) %>%
                             paste0(collapse = "; ")),
           leadingEdge = map_chr(leadingEdge, ~ paste0(.x, collapse = "; "))) %>%
    mutate(pval = round(pval, 5), # hard rounding
           padj = round(padj, 5),
           ES = round(ES, 2),
           NES = round(NES, 2))
}
