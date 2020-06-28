# REMOVED POTENTIAL IDENTIFYING INFORMATION FOR GITHUB

options(stringsAsFactors = F)
library(rio)
setwd("")

# load

keys <- list.files("", full.names = T, recursive = T)
key <- lapply(keys, function(key) {
  
  key <- import(key)
  key[, c("SCAPIS_RID", "Subject ID", "Label ID")]
  
})
key <- do.call(rbind, key)

pheno <- import("")
tax <- import("")
met_samples1 <- import("", which = 2, n_max = 12)
met_samples2 <- import("", which = 2, n_max = 12)
met1 <- import("", which = 2, skip = 12)
met2 <- import("", which = 2, skip = 12)
load("")
module_descriptions <- import("")

# key

key$upugut <- paste0("upugut_", key$"Label ID")

# pheno

pheno$SCAPIS_RID <- pheno$rid

# met

met_samples1 <- met_samples1[, 13:ncol(met_samples1)]
rownames(met_samples1) <- met_samples1[, 1]
met_samples1 <- met_samples1[, -1]
met_samples1 <- as.data.frame(t(met_samples1))
colnames(met_samples1)[colnames(met_samples1) == "Subject ID HMDB"] <- "Subject ID"
met_samples1$SCAPIS_RID <- key$SCAPIS_RID[match(met_samples1$"Subject ID", key$"Subject ID")]

met_info1 <- as.data.frame(met1[, 1:13])
colnames(met_info1)[colnames(met_info1) == "Subject ID HMDB"] <- "HMDB"
met_info1$SUPER_PATHWAY[is.na(met_info1$"SUPER PATHWAY")] <- "Uncharacterized Molecules"

met1 <- met1[, 14:ncol(met1)]
met1 <- as.data.frame(t(met1))
colnames(met1) <- met_info1$"COMP ID"

met_samples2 <- met_samples2[, 13:ncol(met_samples2)]
rownames(met_samples2) <- met_samples2[, 1]
met_samples2 <- met_samples2[, -1]
met_samples2 <- as.data.frame(t(met_samples2))
colnames(met_samples2)[colnames(met_samples2) == "Group HMDB"] <- "Label ID"
met_samples2$SCAPIS_RID <- key$SCAPIS_RID[match(met_samples2$"Label ID", key$"Label ID")]

met_info2 <- as.data.frame(met2[, 1:13])
colnames(met_info2)[colnames(met_info2) == "Group HMDB"] <- "HMDB"
met_info2$SUPER_PATHWAY[is.na(met_info2$SUPER_PATHWAY)] <- "Uncharacterized Molecules"

met2 <- met2[, 14:ncol(met2)]
met2 <- as.data.frame(t(met2))
colnames(met2) <- met_info2$"COMP_ID"

# combine

id <- c(met_info1$"COMP ID", met_info2$COMP_ID)
id <- id[duplicated(id)]
met_info <- met_info2[match(id, met_info2$COMP_ID), ]
met <- rbind(met1[, match(id, colnames(met1))], met2[, match(id, colnames(met2))])

# normalize

met <- apply(met, 2, function(x) {
  
  x <- x  / median(x, na.rm = T)
  x[is.na(x)] <- min(x, na.rm = T)
  x
  
})

met <- data.frame(SCAPIS_RID = c(met_samples1$SCAPIS_RID, met_samples2$SCAPIS_RID), met)
colnames(met) <- gsub("X", "", colnames(met))

# remove outliers 

outliers <- met_samples1$SCAPIS_RID[met_samples1$"SAMPLE NAME" %in% outliers]
met <- met[!met$SCAPIS_RID %in% outliers, ]

# tax

tax_info <- tax[, 1:11]
tax_info$phylum[tax_info$phylum == "unclassified"] <- "Unclassified"
tax <- tax[, 12:ncol(tax)]
tax <- t(tax)
colnames(tax) <- tax_info$MGS
tax <- data.frame(SCAPIS_RID = key$SCAPIS_RID[match(rownames(tax), key$upugut)], tax)

# match data 

key <- key[which(key$SCAPIS_RID %in% met$SCAPIS_RID & key$SCAPIS_RID %in% tax$SCAPIS_RID & key$SCAPIS_RID %in% pheno$SCAPIS_RID), ]
key <- key[!duplicated(key$SCAPIS_RID), ]
met <- met[match(key$SCAPIS_RID, met$SCAPIS_RID), ]
tax <- tax[match(key$SCAPIS_RID, tax$SCAPIS_RID), ]
pheno <- pheno[match(key$SCAPIS_RID, pheno$SCAPIS_RID), ]

# phenotypes

age <- pheno$agev1
gender <- ifelse(pheno$Gender == 1, "Female", "Male")
smoke <- factor(ifelse(pheno$smokestatus == 3, "Yes", "No"))
bp <- pheno$bbps
bmi <- pheno$bmi
ethnicity <- ifelse(pheno$q005 == 0, "Other", "Swedish")

pheno <- data.frame(SCAPIS_RID = pheno$SCAPIS_RID, age, gender, ethnicity, bmi, bp, smoke)

# module_info

module_info <- lapply(1:length(HGMGS.keggModule2MGS), function(i) {
  
  if (length(HGMGS.keggModule2MGS[[i]]) > 0) data.frame(module = names(HGMGS.keggModule2MGS)[i], mgs = HGMGS.keggModule2MGS[[i]])
  
})
module_info <- do.call(rbind, module_info)
module_info <- merge(module_info, module_descriptions[, 1:4], by.x = "module", by.y = "Module")

save(tax, tax_info, met, met_info, pheno, module_info, file = "data/data.RData")
