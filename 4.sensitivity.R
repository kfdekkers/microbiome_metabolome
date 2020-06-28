# REMOVED POTENTIAL IDENTIFYING INFORMATION FOR GITHUB

# set options and load libraries

options(stringsAsFactors = F)
set.seed(1)
cores <- 16

library(rio)
library(BiocParallel)
library(caret)
library(ppcor)

setwd("")

load("data/data.RData")
spearman <- import("results/spearman.tsv")

# remove id

met <- met[, 2:ncol(met)]
tax <- tax[, 2:ncol(tax)]

# remove features with zero variance

met <- predict(preProcess(met, method = "zv"), met)
tax <- predict(preProcess(tax, method = "zv"), tax)

drugs <- met_info$COMP_ID[grepl("Drug", met_info$SUB_PATHWAY)]
drugs <- met_info$COMP_ID[met_info$BIOCHEMICAL %in% c("losarthan", "metformin", "doxycycline", "fluconazole", "N-acetyl sulfapyridine", "quinine")]
drugs <- met[, colnames(met) %in% drugs]
drugs <- apply(drugs, 2, function(x) {
  
  ifelse(x > min(x), "Yes", "No")

})
sens <- data.frame(pheno[, c("bmi", "bp", "smoke")], drugs)
colnames(sens) <- gsub("X", "", colnames(sens))

# correlation function

cor.fun <- function(y, x, z) {
  
  data <- data.frame(x, z)
  complete <- complete.cases(cbind(y, data))
  data <- data[complete, ]
  y <- y[complete]
  data <- model.matrix(~ ., data)
  
  fit <- tryCatch({
    
    data.frame(pcor.test(y, data[, 2], data[, 3:ncol(data)]), warning = NA)
    
  }, warning = function(w) {
    
    data.frame(pcor.test(y, data[, 2], data[, 3:ncol(data)]), warning = w$message)
    
  })
  
  fit[, c("estimate", "p.value", "n", "warning")]
  
}

# model functions

sensitivity.fun <- function() {
  
  met <- apply(met, 2, function(x) rank(x, "keep"))
  tax <- apply(tax, 2, function(x) rank(x, "keep"))
  
  pairs <- paste(spearman$compid, spearman$mgs, sep = "_")[spearman$fdr < 0.05]

  res <- bplapply(colnames(sens), function(x) {
    
    res <- bplapply(pairs, function(pair) {
      
      pair <- unlist(strsplit(pair, split = "_"))
      
      # if (pair[1] == x) {
      #   
      #   data.frame(estimate = NA, p.value = NA, n = NA, warning = NA)
      #   
      # } else {
      #   
      #   cor.fun(met[, pair[1]], tax[, pair[2]], data.frame(pheno[, c("age", "gender", "ethnicity")], sens[, x]))
      #   
      # }
      
      cor.fun(met[, pair[1]], tax[, pair[2]], data.frame(pheno[, c("age", "gender", "ethnicity")], sens[, x]))
      
    })
    
    res <- do.call(rbind, res)
    res$fdr <- p.adjust(res$p.value, method = "fdr")
    res <- res[, c("estimate", "p.value", "fdr", "n", "warning")]
    colnames(res) <- paste(x, colnames(res), sep = "_")
    res
    
  }, BPPARAM = MulticoreParam(cores))
    
  res <- do.call(cbind, res)
  data.frame(compid = spearman$compid[spearman$fdr < 0.05], mgs = spearman$mgs[spearman$fdr < 0.05], res)
    
}

sensitivity <- sensitivity.fun()
colnames(sensitivity) <- gsub("X", "", colnames(sensitivity))
export(sensitivity, "results/sensitivity.tsv")

