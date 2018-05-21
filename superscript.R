if(!"tcltk" %in% rownames(installed.packages())){install.packages('tcltk')}
if(!"mgcv" %in% rownames(installed.packages())){install.packages('mgcv')}
if(!"plyr" %in% rownames(installed.packages())){install.packages('plyr')}

library("tcltk")

parent <- paste(as.character(tkchooseDirectory()), collapse = " ")
setwd(parent)

dir.create(paste(as.character(parent), '/output', sep=""), showWarnings = F)
outputdir <- paste(as.character(parent), '/output', sep="")
dirs <- list.dirs(recursive = T)


file <- list.files(path = dirs, full.names = T)
scripts <- file[grepl('R$', file)]
scripts <- scripts[!grepl('superscript', scripts)]
flux <- scripts[grepl('getFlux', scripts)]
scripts <- scripts[!grepl('getFlux', scripts)]


sapply(scripts, function(s){
  setwd(parent)
  mydir <- gsub('^.', parent, dirname(s))
  source(s)
  setwd(mydir)
  file <- list.files()
  file <- file[grepl('.txt$', file)]
  sapply(file, function(f){
    setwd(mydir)
    do(f, outputdir)
  })
})
print('data detrended')
setwd(parent)
if(length(flux) == 0){print('no getFlux.R found')
}else{
  source(flux)
  dir.create(paste(as.character(parent), '/rates output', sep=""), showWarnings = F)
  rateoutputdir <- paste(as.character(parent), '/rates output', sep="")
  getFlux(dir = outputdir, outputdir = rateoutputdir)
  print('flux calculated')
  setwd(parent)
}

