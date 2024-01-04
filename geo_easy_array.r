
rm(list = ls())

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

# Set the following three parameters
GSE_number = "GSE15907" # human example
# GSE_number = "GSE66356" # mouse example

#organism = "Human" # or "Mouse"
organism = "Mouse" # or "Human"

# Output Folder path
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

platform
full_table <- read_delim(url(paste("https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/GPL_Files/", platform, ".txt", sep="")), delim = '\t', col_names = T, )
#rownames(full_table) <- full_table[,1]
full_table[,1] <- NULL

head(full_table)
#full_table <- idmap(platform)
# full_table <- idmap(platform, type = 'pipe')
colnames(full_table)[1] = "ID"
colnames(full_table)[2] = "Symbol"
print(head(full_table))
full_table <- arrange(full_table)
full_table = data.frame(full_table)
full_table$ID <- as.character(full_table$ID)
df$ID <- as.character(df$ID)
head(full_table)
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

#colnames(converted_expr_matrix)[1] <- "Symbol"
#converted_expr_matrix <- converted_expr_matrix %>%
#  group_by(Symbol) %>%
#  summarise_all(max)

converted_expr_matrix <- as.data.frame(converted_expr_matrix)

head(converted_expr_matrix)

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

meta_merged %>%
  mutate_if(is.character, str_trim)

# check whether the matrix has been logged
hist(converted_expr_matrix[,2])
hist(log(converted_expr_matrix[,2]))

if (quantile(converted_expr_matrix[,2])[4] < 12) {
  converted_expr_matrix_exp = 2^converted_expr_matrix[,-1]
  new_exp = cbind(converted_expr_matrix[,1], converted_expr_matrix_exp);
  colnames(new_exp)[1] = "ID"
  converted_expr_matrix = new_exp
}

outputmatrix = paste(outputpath, "/", geo_accession, ".matrix.txt", sep="")
outputmeta = paste(outputpath, "/", geo_accession, ".meta.txt", sep="")
outputnote = paste(outputpath, "/", geo_accession, ".note.txt", sep="")
write.table(converted_expr_matrix, file=outputmatrix,row.names=F, sep="\t")
write.table(meta_merged, file=outputmeta, row.names=FALSE, sep="\t")
write.table(rbind(meta$treatment_protocol_ch1[1], meta$organism_ch1[1]), file=outputnote)

# generation of the EASY app
# It'll be a good idea to double check the files now before proceeding with generating the apps
UPPERCASE_organism = toupper(organism)
if (UPPERCASE_organism == "HUMAN") {

  if (meta$organism_ch1[1] != "Homo sapiens") {
    print("Please double check whether samples are Human in the meta file...");
    break;
  }
  output_dir <- file.path(outputpath, paste(geo_accession, "_Human_EASY_App", sep=""))
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    #print("Dir already exists!")
  }
  R_script_output_dir <- file.path(outputpath, paste(geo_accession, "_Human_EASY_App/R", sep=""))

  if (!dir.exists(R_script_output_dir)){
    dir.create(R_script_output_dir)
  } else {
    #print("Dir already exists!")
  }


  outputmatrix = paste(output_dir, "/", geo_accession, ".matrix.txt", sep="")
  outputmeta = paste(output_dir, "/", geo_accession, ".meta.txt", sep="")
  outputnote = paste(output_dir, "/", geo_accession, ".note.txt", sep="")
  write.table(converted_expr_matrix, file=outputmatrix,row.names=F, sep="\t")
  write.table(meta_merged, file=outputmeta, row.names=FALSE, sep="\t")
  write.table(rbind(meta$treatment_protocol_ch1[1], meta$organism_ch1[1]), file=outputnote)

  hs_LINC1000_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/e065189a65f861024ac10f683f1595d44f479ffa/EASYAppFiles/Template_Human_EASY_App/LINCS_L1000_gsNsym_HS_v2.zip";
  hs_linc1000_rdata_url = "https://github.com/shawlab-moffitt/GEO_DataDownload_Example/raw/main/EASYAppFiles/Template_Human_EASY_App/LINCS_L1000_gs_HS_v2.RData";
  hs_app_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Human_EASY_App/app.R";
  hs_msigdb_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Human_EASY_App/msigdb_gsNcat_HS_v2.txt"
  hs_msigdb_zip_url = "https://github.com/shawlab-moffitt/GEO_DataDownload_Example/raw/main/EASYAppFiles/Template_Human_EASY_App/msigdb_gsNsym_HS_v2.zip"
  hs_msigdb_rdata_url = "https://github.com/shawlab-moffitt/GEO_DataDownload_Example/raw/main/EASYAppFiles/Template_Human_EASY_App/msigdb_gs_HS_v2.RData"

  internalR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Human_EASY_App/R/internal.R";
  loginR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Human_EASY_App/R/login.R"
  logoutR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Human_EASY_App/R/logout.R"
  runExampleR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Human_EASY_App/R/runExample.R"


  output_file_LINC1000 = paste(output_dir, "/LINCS_L1000_gsNsym_HS_v2.zip", sep="")
  output_file_rdata = paste(output_dir, "/LINCS_L1000_gs_HS_v2.RData", sep="")
  output_file_app = paste(output_dir, "/app.R", sep="")
  output_file_msigdb = paste(output_dir, "/msigdb_gsNcat_HS_v2.txt", sep="")
  output_file_msigdbzip = paste(output_dir, "/msigdb_gsNsym_HS_v2.zip", sep="")
  output_file_msigdbrdata = paste(output_dir, "/msigdb_gs_HS_v2.RData", sep="")
  output_file_internalR = paste(R_script_output_dir, "/internal.R", sep="")
  output_file_loginR = paste(R_script_output_dir, "/login.R", sep="")
  output_file_logoutR = paste(R_script_output_dir, "/logout.R", sep="")
  output_file_runExampleR = paste(R_script_output_dir, "/runExample.R", sep="")


  download.file(hs_LINC1000_url, output_file_LINC1000)
  download.file(hs_linc1000_rdata_url, output_file_rdata)
  download.file(hs_app_url, output_file_app);
  download.file(hs_msigdb_url, output_file_msigdb);
  download.file(hs_msigdb_zip_url, output_file_msigdbzip);
  download.file(hs_msigdb_rdata_url, output_file_msigdbrdata);

  download.file(internalR_url, output_file_internalR);
  download.file(loginR_url, output_file_loginR);
  download.file(logoutR_url, output_file_logoutR);
  download.file(runExampleR_url, output_file_runExampleR);

  # replace the expr_and meta file
  tx  <- readLines(output_file_app)
  tx2 <- gsub(pattern = "expr_file <- ''", replace = paste("expr_file <- '", geo_accession, ".expr.matrix.txt", "'", sep=""), x = tx)
  tx3 <- gsub(pattern = "meta_file <- ''", replace = paste("meta_file <- '", geo_accession, ".meta.txt", "'", sep=""), x = tx2)
  tx4 <- gsub(pattern = "ProjectName <- ''", replace = paste("ProjectName <- '", geo_accession, " Study", "'", sep=""), x = tx3)
  writeLines(tx4, con=output_file_app)


}
# write mouse samples
if (UPPERCASE_organism == "MOUSE") {

  if (meta$organism_ch1[1] != "Mus musculus") {
    print("Please double check whether samples are Human in the meta file...");
    break;
  }
  output_dir <- file.path(outputpath, paste(geo_accession, "_Mouse_EASY_App", sep=""))
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  } else {
    #print("Dir already exists!")
  }
  R_script_output_dir <- file.path(outputpath, paste(geo_accession, "_Mouse_EASY_App/R", sep=""))

  if (!dir.exists(R_script_output_dir)){
    dir.create(R_script_output_dir)
  } else {
    #print("Dir already exists!")
  }


  outputmatrix = paste(output_dir, "/", geo_accession, ".expr.matrix.txt", sep="")
  outputmeta = paste(output_dir, "/", geo_accession, ".meta.txt", sep="")
  outputnote = paste(output_dir, "/", geo_accession, ".note.txt", sep="")
  write.table(converted_expr_matrix, file=outputmatrix,row.names=F, sep="\t")
  write.table(meta_merged, file=outputmeta, row.names=FALSE, sep="\t")
  write.table(rbind(meta$treatment_protocol_ch1[1], meta$organism_ch1[1]), file=outputnote)

  CellMarker_GS_MM_rdata_url = "https://github.com/shawlab-moffitt/GEO_DataDownload_Example/raw/main/EASYAppFiles/Template_Mouse_EASY_App/CellMarker_GS_MM.RData";
  CellMarker_gsNsym_MM_tsv_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Mouse_EASY_App/CellMarker_gsNsym_MM.tsv";
  mm_app_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Mouse_EASY_App/app.R";
  mm_msigdb_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Mouse_EASY_App/msigdb_gsNcat_MM.tsv"
  mm_msigdb_zip_url = "https://github.com/shawlab-moffitt/GEO_DataDownload_Example/raw/main/EASYAppFiles/Template_Mouse_EASY_App/msigdb_gsNsym_MM.zip"
  mm_msigdb_rdata_url = "https://github.com/shawlab-moffitt/GEO_DataDownload_Example/raw/main/EASYAppFiles/Template_Mouse_EASY_App/msigdb_gs_MM.RData"

  internalR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Mouse_EASY_App/R/internal.R";
  loginR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Mouse_EASY_App/R/login.R"
  logoutR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Mouse_EASY_App/R/logout.R"
  runExampleR_url = "https://raw.githubusercontent.com/shawlab-moffitt/GEO_DataDownload_Example/main/EASYAppFiles/Template_Mouse_EASY_App/R/runExample.R"


  output_file_CellMarker_rdata = paste(output_dir, "/CellMarker_GS_MM.RData", sep="")
  output_file_CellMarker_tsv = paste(output_dir, "/CellMarker_gsNsym_MM.tsv", sep="")
  output_file_app = paste(output_dir, "/app.R", sep="")
  output_file_msigdb = paste(output_dir, "/msigdb_gsNcat_MM.tsv", sep="")
  output_file_msigdbzip = paste(output_dir, "/msigdb_gsNsym_MM.zip", sep="")
  output_file_msigdbrdata = paste(output_dir, "/msigdb_gs_MM.RData", sep="")
  output_file_internalR = paste(R_script_output_dir, "/internal.R", sep="")
  output_file_loginR = paste(R_script_output_dir, "/login.R", sep="")
  output_file_logoutR = paste(R_script_output_dir, "/logout.R", sep="")
  output_file_runExampleR = paste(R_script_output_dir, "/runExample.R", sep="")


  download.file(CellMarker_GS_MM_rdata_url, output_file_CellMarker_rdata)
  download.file(CellMarker_gsNsym_MM_tsv_url, output_file_CellMarker_tsv)
  download.file(mm_app_url, output_file_app);
  download.file(mm_msigdb_url, output_file_msigdb);
  download.file(mm_msigdb_zip_url, output_file_msigdbzip);
  download.file(mm_msigdb_rdata_url, output_file_msigdbrdata);

  download.file(internalR_url, output_file_internalR);
  download.file(loginR_url, output_file_loginR);
  download.file(logoutR_url, output_file_logoutR);
  download.file(runExampleR_url, output_file_runExampleR);

  # replace the expr_and meta file
  tx  <- readLines(output_file_app)
  tx2 <- gsub(pattern = "expr_file <- ''", replace = paste("expr_file <- '", geo_accession, ".matrix.txt", "'", sep=""), x = tx)
  tx3 <- gsub(pattern = "meta_file <- ''", replace = paste("meta_file <- '", geo_accession, ".meta.txt", "'", sep=""), x = tx2)
  tx4 <- gsub(pattern = "ProjectName <- ''", replace = paste("ProjectName <- '", geo_accession, " Study", "'", sep=""), x = tx3)
  writeLines(tx4, con=output_file_app)


}
