library(AnnoProbe)
gpls = getGPLList()[1]

for (i in 1:length(gpls[,1])) {
  gpl = gpls[i,]
  print(gpl)
  full_table <- idmap(platform)
  filename = paste("/Users/4472414/Documents/Current_Manuscripts/GSE_Download_Tools/GPLS/", gpl, ".txt", sep="")
  write.table(full_table, file=filename, col.names=NA, sep="\t");
}
