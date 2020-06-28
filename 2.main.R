# REMOVED POTENTIAL IDENTIFYING INFORMATION FOR GITHUB

# set options and load libraries

options(stringsAsFactors = F)
set.seed(1)
cores <- 16
filter <- 20

library(rio)
library(BiocParallel)
library(caret)
library(ppcor)

setwd("")

# regression functions

lm.fun <- function(y, x, z) {
  
  data <- data.frame(x, z)
  complete <- complete.cases(cbind(y, data))
  data <- data[complete, ]
  y <- y[complete, ]
  fit <- lm(as.matrix(y) ~ ., data = data)
  coef <- lapply(summary(fit), function(x) x$coefficients[2, ])
  coef <- as.data.frame(do.call(rbind, coef))
  data.frame(estimate = coef[, 1], se = coef[, 2], p.value = coef[, 4], n = nrow(data))
  
}

cor.fun <- function(y, x, z) {

  data <- data.frame(x, z)
  complete <- complete.cases(cbind(y, data))
  data <- data[complete, ]
  y <- y[complete, ]
  data <- model.matrix(~ ., data)

  fit <- apply(y, 2, function(ycol) {

    tryCatch({

      data.frame(pcor.test(ycol, data[, 2], data[, 3:ncol(data)]), warning = NA)

    }, warning = function(w) {

      data.frame(pcor.test(ycol, data[, 2], data[, 3:ncol(data)]), warning = w$message)

    })

  })

  fit <- do.call(rbind, fit)
  fit[, c("estimate", "p.value", "n", "warning")]

}

# model functions

spearman.fun <- function() {
  
  met <- apply(met, 2, function(x) rank(x, "keep"))
  tax <- apply(tax, 2, function(x) rank(x, "keep"))
  
  res <- bplapply(colnames(tax), function(x) {
    
    coef <- cor.fun(met, tax[, x], pheno[, c("age", "gender", "ethnicity")])
    
    data.frame(compid = colnames(met), mgs = x, coef)
    
  }, BPPARAM = MulticoreParam(cores))
  
  res <- do.call(rbind, res)
  res$fdr <- p.adjust(res$p.value, method = "fdr")
  res[, c("compid", "mgs", "estimate", "p.value", "fdr", "n", "warning")]
  
}

lmlog1p.fun <- function() {
  
  met <- apply(met, 2, log1p)
  tax <- apply(tax, 2, log1p)
  
  res <- bplapply(colnames(tax), function(x) {
    
    coef <- lm.fun(met, tax[, x], pheno[, c("age", "gender", "ethnicity")])
    
    data.frame(compid = colnames(met), mgs = x, coef)
    
  }, BPPARAM = MulticoreParam(cores))
  
  res <- do.call(rbind, res)
  res$fdr <- p.adjust(res$p.value, method = "fdr")
  res[, c("compid", "mgs", "estimate", "se", "p.value", "fdr", "n")]
  
}

# import data

load("data/data.RData")

# remove id

met <- met[, 2:ncol(met)]
tax <- tax[, 2:ncol(tax)]

# remove features with zero variance

met <- predict(preProcess(met, method = "zv"), met)
tax <- predict(preProcess(tax, method = "zv"), tax)

# filter

met_filter <- apply(met, 2, function(x) sum(x > min(x)) > filter)
met <- met[, met_filter]
tax_filter <- apply(tax, 2, function(x) sum(x > min(x)) > filter)
tax <- tax[, tax_filter]

# run models

spearman <- spearman.fun()
export(spearman, "results/spearman.tsv")
lmlog1p <- lmlog1p.fun()
export(lmlog1p, "results/lm.tsv")
