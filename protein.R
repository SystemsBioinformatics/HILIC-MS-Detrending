ProteinConc <- function(file, outputdir){
  mydir <- getwd()
  data <- read.csv(file, sep = "\t", header = TRUE, dec = ".")
  names(data) <- tolower(names(data))
  setwd(outputdir)
  write.table(data, file = "biomass.txt", sep = '\t', row.names = F)
  setwd(mydir)
}
do <- ProteinConc