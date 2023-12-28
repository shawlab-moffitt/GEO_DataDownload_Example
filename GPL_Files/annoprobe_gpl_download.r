install.packages("magrittr")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("GEOquery")

#library(devtools)
#install_github("jmzeng1314/AnnoProbe")
library(AnnoProbe)
gpls = getGPLList()[1]

for (i in 1:length(gpls[,1])) {
  print(i)
  gpl = gpls[i,]
  print(gpl)
  tryCatch(
    
    expr = {
      full_table = "";
      full_table <- idmap(gpl);
      filename = paste("/Users/4472414/Documents/Current_Manuscripts/GSE_Download_Tools/GPLs/", gpl, ".txt", sep="");
      write.table(full_table, file=filename, col.names=NA, sep="\t");
      
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    }
  )
  
}
a
getwd()
