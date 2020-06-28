# REMOVED POTENTIAL IDENTIFYING INFORMATION FOR GITHUB

# set options and load libraries

options(stringsAsFactors = F)
set.seed(1)
cores <- 16

library(rio)
library(BiocParallel)

setwd("")
spearman <- import("results/spearman.tsv")

# module fun

module.fun <- function(summary.stats) {
  
  modules <- unique(module_info$module)
  
  res_module <- lapply(modules, function(module) {
    
    mgs <- module_info$mgs[which(module_info$module == module)]
    
    tryCatch({
      
      as.numeric(unlist(ks.test(summary.stats$p.value[which(summary.stats$mgs %in% mgs)], summary.stats$p.value[which(!summary.stats$mgs %in% mgs)])$p.value))
      
    }, error = function(e) {
      
      NA
      
    })
    
  })
  
  data.frame(module = modules, p.value = unlist(res_module))
  
}

spearman.module <- bplapply(unique(spearman$compid), function(compid) {
  
  res <- module.fun(spearman[spearman$compid == compid, ])
  data.frame(compid = compid, res)
  
}, BPPARAM = MulticoreParam(cores))
spearman.module <- do.call(rbind, spearman.module)
spearman.module$fdr <- p.adjust(spearman.module$p.value, method = "fdr")
export(spearman.module, "results/spearman.module.tsv")

# phylum enrichment analysis

phylum.fun <- function(summary.stats) {
  
  phylae <- unique(tax_info$phylum)
  
  res_phylum <- lapply(phylae, function(phylum) {
    
    mgs <- tax_info$MGS[which(tax_info$phylum == phylum)]
    
    tryCatch({
      
      as.numeric(unlist(ks.test(summary.stats$p.value[which(summary.stats$mgs %in% mgs)], summary.stats$p.value[which(!summary.stats$mgs %in% mgs)])$p.value))

    }, error = function(e) {
      
      NA
      
    })
    
  })
  
  data.frame(phylum = phylae, p.value = unlist(res_phylum))
  
}

spearman.phylum <- bplapply(unique(spearman$compid), function(compid) {
  
  res <- phylum.fun(spearman[spearman$compid == compid, ])
  data.frame(compid = compid, res)
  
}, BPPARAM = MulticoreParam(cores))
spearman.phylum <- do.call(rbind, spearman.phylum)
spearman.phylum$fdr <- p.adjust(spearman.phylum$p.value, method = "fdr")
export(spearman.phylum, "results/spearman.phylum.tsv")

# sub pathway enrichment

pathway.fun <- function(summary.stats) {
  
  pathways <- unique(met_info$SUB_PATHWAY)
  
  res_pathway <- lapply(pathways, function(pathway) {
    
    compids <- met_info$COMP_ID[which(met_info$SUB_PATHWAY == pathway)]
    
    tryCatch({
      
      as.numeric(unlist(ks.test(summary.stats$p.value[which(summary.stats$compid %in% compids)], summary.stats$p.value[which(!summary.stats$compid %in% compids)])$p.value))

    }, error = function(e) {
      
      NA
      
    })
    
  })
  
  data.frame(pathway = pathways, p.value = unlist(res_pathway))
  
}

spearman.pathway <- bplapply(unique(spearman$mgs), function(mgs) {
  
  res <- pathway.fun(spearman[spearman$mgs == mgs, ])
  data.frame(mgs = mgs, res)
  
}, BPPARAM = MulticoreParam(cores))
spearman.pathway <- do.call(rbind, spearman.pathway)
spearman.pathway$fdr <- p.adjust(spearman.pathway$p.value, method = "fdr")
export(spearman.pathway, "results/spearman.pathway.tsv")
