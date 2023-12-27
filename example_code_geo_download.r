
library(AnnoProbe)
library(GEOquery)
library(limma)
library(umap)
library(tidyverse)
library(readr)
library(dplyr)
library(R.utils)
library(shiny)
library(devtools)


#gpls = getGPLList()[1]

#for (i in 1:length(gpls[,1])) {
#  gpl = gpls[i,]
#  print(gpl)
#  full_table <- idmap(platform)
#  filename = paste("/Users/4472414/Documents/Current_Manuscripts/GSE_Download_Tools/GPLS/", gpl, ".txt", sep="")
#  write.table(full_table, file=filename, col.names=NA, sep="\t");
#}

GSE_number = "GSE50532"

outputpath = "/Users/4472414/Documents/Current_Manuscripts/GSE_Download_Tools/Microarray_Data_Download/"


gset <- getGEO(GSE_number, GSEMatrix = TRUE, getGPL = TRUE)
gse <- getGEO(GSE_number, GSEMatrix = FALSE, getGPL = TRUE)
platform <- names(GPLList(gse))
print(platform)

if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
expr <- exprs(gset)
df <- as.data.frame(expr)
df <- arrange(df)
df$ID <- rownames(df)
df <- df %>% relocate(ID)
df$ID <- trimws(df$ID)
df$ID <- as.character(df$ID)
?idmap
full_table <- idmap(platform)
# full_table <- idmap(platform, type = 'pipe')
colnames(full_table)[1] = "ID"
colnames(full_table)[2] = "Symbol"
print(head(full_table))
full_table <- arrange(full_table)

# convert id to symbol for expr matrix
converted_expr_matrix <- inner_join(df, full_table, by = "ID")
converted_expr_matrix <- relocate(converted_expr_matrix, "Symbol", .before = 1)
converted_expr_matrix$ID <- NULL
converted_expr_matrix <- converted_expr_matrix[complete.cases(converted_expr_matrix$Symbol), ]
converted_expr_matrix$Symbol <- trimws(converted_expr_matrix$Symbol)
colnames(converted_expr_matrix)[1] <- "Gene"

# Remove Expression with NA
converted_expr_matrix <- converted_expr_matrix %>%
  drop_na()

# Check that expression data is numeric
isChar <- unname(which(sapply(converted_expr_matrix, function(x) is.character(x))))
isChar <-  isChar[-1]
if (length(isChar) > 0) {
  converted_expr_matrix[isChar] <- sapply(converted_expr_matrix[isChar],as.numeric)
}

# Remove Duplicate genes
if (TRUE %in% duplicated(converted_expr_matrix[,1])) {
  converted_expr_matrix <- converted_expr_matrix %>%
    group_by(Gene) %>%
    summarise_all(max)
}
converted_expr_matrix <- as.data.frame(converted_expr_matrix)

geo_accession = GSE_number;
meta <- pData(gset)
meta <- relocate(meta, geo_accession, .before = title)

# cleaning the meta data
# Select columns that start with "characteristic"
char_cols <- grep("^characteristic", colnames(meta), value = TRUE)
extr_meta <- meta[, c("geo_accession", char_cols)]

original_df <- extr_meta

# Initialize an empty data frame to store results
new_df <- data.frame(sample = character(0), title = character(0), value = character(0), stringsAsFactors = FALSE)

# Iterate through rows and columns to extract values and create new rows
for (i in seq_along(original_df$geo_accession)) {
  row_values <- as.character(original_df[i, -1])  # Convert to character
  
  for (j in seq_along(row_values)) {
    split_value <- strsplit(row_values[j], ":")[[1]]
    title <- split_value[1]
    value <- paste(split_value[-1], collapse = ":")
    
    new_row <- data.frame(geo_accession = original_df$geo_accession[i], title = title, value = value, stringsAsFactors = FALSE)
    new_df <- rbind(new_df, new_row)
  }
}

final_extr_meta <- tidyr::pivot_wider(new_df, names_from = title, values_from = value)
meta_merged <- dplyr::full_join(final_extr_meta, meta, by = "geo_accession")

outputmatrix = paste(outputpath, "/", geo_accession, ".matrix.txt", sep="")
outputmeta = paste(outputpath, "/", geo_accession, ".meta.txt", sep="")
write.table(converted_expr_matrix, file=outputmatrix,col.names=NA, sep="\t")
write.table(meta_merged, file=outputmeta, col.names=NA, sep="\t")
