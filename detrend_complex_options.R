require("mgcv")

#Function to calculate new concentrations from model
getconc <- function(yvalue, model){
  return((yvalue - coef(model)[[1]]) / coef(model)[[2]])
}

detrend <- function(file, outputdir, calibrations = 2, date.format = "%d/%m/%Y %H:%M", do.detrend = F, ...){
  #reading the data
  mydir <- getwd()
  data <- read.csv(file, sep = "\t", header = TRUE, dec = ".", comment.char = "#")
  names(data) <- tolower(names(data))
  #####added#10-12-14#######
  data$type <- tolower(data$type)
  data <- data[!is.na(data[,1]),]
  data <- data[!data[,1]=="",]
  try(data <- data[,-c(which(TRUE == grepl('^x$', colnames(data))):ncol(data))], silent = T)
  #######stop#########
  data$date.acquired <- as.POSIXct(as.character(data$date.acquired), format = date.format)
  #####added#10-12-14#######
  if(length(unique(data$date.acquired))==1){
    data$date.acquired <- 1:length(data$date.acquired)
    print(paste('date.format', as.character(file), ' != ',date.format))
  }
  ######stop########
  setwd(outputdir)
  pdf(paste('detrending ', gsub(".txt", ".pdf", file)))
  
  #creating factors for GAM
  
  if(dim(data)[2] < 9){
    d <- data[,c('type', 'concentration')]
  }else {
    d <- data[,c('type','concentration', names(data[,9:dim(data)[2]]))]
  }
  
  all.factors <- ifelse(!is.na(d$concentration), gsub(" ", "_",do.call('paste', d[-which(names(d)=='type')])), gsub(" ", "_",do.call('paste', d)) )
  all.factors <- as.factor(all.factors)
  
  
  #getting the data of calibration points
  calcurvedata <- data[data$type == 'calibration.point',]
  if(nrow(calcurvedata)==0){
    print(paste('no calibration points in file ','\'', file, '\'', sep = ""))
    stop()
  }
    
  calcurves <- as.list(1:calibrations)
  r <- c(as.integer(rownames(calcurvedata)), 1)
  ####ordering calibration curves####  
  rows <- lapply(calcurves, function(x){
    s <- r[1]
    i <- 2
    while(r[i]==r[i-1]+1){
      s <- c(s, r[i])
      i <- i+1
    }
    r <<- r[-c(1:(i-1))]
    return(s)
  })
  calcurves <- lapply(rows, function(i){
    i <- as.data.frame(data[as.character(i),])
    return(i)
  })
  missing <- rownames(calcurvedata)
  missing <- lapply(calcurves, function(i){
    missing <<- setdiff(missing, rownames(i))
    return(missing)
  })
  if(length(missing[[calibrations]])>0){
    print(paste('In ',file,' calibrations are incorrectly ordered, or incorrect number of calibration curves'))
  }
  
 
  
  blanks <- data[data$type == 'blank',]
  controls <- data[data$type == 'control',]
  
  #making the basic plots to investigate the data
  plot(data$date.acquired, data$istd.area, pch = 18, ylab = "Internal Standard Area", xlab ="Date Acquired", main = "Internal controls")
  points(controls$date.acquired, controls$istd.area, pch = 18, col = 'sienna')
  points(blanks$date.acquired, blanks$istd.area, pch = 18, col = 'green')
  points(data[data$type == 'calibration.point',]$date.acquired, data[data$type == 'calibration.point',]$istd.area, pch = 18, col = 'darkblue')
  legend('bottomleft', c('Samples','Blanks', 'Controls', 'Calibration'), pch = 18, col = c('black', 'green', 'sienna', 'darkblue'))
  
  plot(data$date.acquired, data$area, pch = 18, ylab = "Area", xlab ="Date Acquired", main = "Samples")
  points(controls$date.acquired, controls$area, pch = 18, col = 'sienna')
  points(blanks$date.acquired, blanks$area, pch = 18, col = 'green')
  points(data[data$type == 'calibration.point',]$date.acquired, data[data$type == 'calibration.point',]$area, pch = 18, col = 'darkblue')
  legend('bottomleft', c('Samples','Blanks', 'Controls', 'Calibration'), pch = 18, col = c('black', 'green', 'sienna', 'darkblue'))
  
  plot(data$date.acquired, data$area.ratio, pch = 18, ylab = "Area Ratio", xlab ="Date Acquired", main = "area ratio")
  points(controls$date.acquired, controls$area.ratio, pch = 18, col = 'sienna')
  points(blanks$date.acquired, blanks$area.ratio, pch = 18, col = 'green')  
  points(data[data$type == 'calibration.point',]$date.acquired, data[data$type == 'calibration.point',]$area.ratio, pch = 18, col = 'darkblue')
  legend('bottomleft', c('Samples','Blanks', 'Controls', 'Calibration'), pch = 18, col = c('black', 'green', 'sienna', 'darkblue'))
  
  #making the calibration curves
  models <- lapply(calcurves, function(l){
    mdls <- lm(l$area.ratio ~ l$concentration)
    return(mdls)
  })
  joint.model <- lm(data[!is.na(data$concentration),]$area.ratio ~ data[!is.na(data$concentration),]$concentration)
  
  #plotting the calibration curves
  y <- 1:calibrations
  plot(data[!is.na(data$concentration),]$area.ratio ~ data[!is.na(data$concentration),]$concentration, pch = 20, xlab = 'Concentration', ylab = 'area ratio', main = 'calibration curves')
  abline(joint.model)
  sapply(y, function(i){
    points(calcurves[[i]]$concentration, calcurves[[i]]$area.ratio, pch = 20, col = rainbow(calibrations)[i])
    abline(models[[i]], col = rainbow(calibrations)[i])
  })
  legend('bottomright', c(paste('calibration', as.character(y)), 'joint'), pch = 20, col = c(rainbow(calibrations), 'black'))
  
  #log transforming the data
  if(do.detrend == T){
    logdata <- data.frame(area = log(data$area, 10), istd.area = log(data$istd.area, 10))
    logdata$area.ratio <- logdata$area - logdata$istd.area
    
    gam.fit <- try(gam(logdata$area.ratio ~ s(as.numeric(date.acquired), ...) + all.factors, data = data), silent = T)
    if(inherits(gam.fit, 'try-error')){
      writeLines(paste('for', file,'GAM Model has more coefficients than data. \nAdd proper k argument to detrend()'))
      dev.off()
      setwd(mydir)
      stop()
    }
    detrended <- data[,1:4]; detrended$concentration <- data$concentration
    detrended$area.ratio <- 10^(predict(gam.fit, type = "terms")[, "all.factors"]) + 10^(residuals(gam.fit))
  
  
    #getting the detrended data of calibration points
    calcurvedata <- detrended[detrended$type == 'calibration.point',]
    calcurves <- lapply(rows, function(i){
      i <- as.data.frame(calcurvedata[as.character(i),])
      return(i)
    })
  
  #plotting the detrended data
    plot(data$date.acquired, detrended$area.ratio, pch = 18, ylab = "Area Ratio", xlab ="Date Acquired", main = "detrended area ratio")
    points(data[data$type=='control',]$date.acquired, detrended[data$type=='control',]$area.ratio, pch = 18, col = 'sienna')
    points(data[data$type=='blank',]$date.acquired, detrended[data$type=='blank',]$area.ratio, pch = 18, col = 'green')  
    points(data[data$type == 'calibration.point',]$date.acquired, detrended[data$type == 'calibration.point',]$area, pch = 18, col = 'darkblue')
    legend('bottomleft', c('Internal controls','Blanks', 'Controls', 'Calibration'), pch = 18, col = c('black', 'green', 'sienna', 'darkblue'))
  
  #plotting the detrended calibration curve
    detrended.model <- lm(detrended[!is.na(detrended$concentration),]$area.ratio ~ detrended[!is.na(data$concentration),]$concentration)
    plot(area.ratio ~ concentration, data = detrended[!is.na(detrended$concentration),], pch = 18)
    sapply(y, function(i){
      points(calcurves[[i]]$concentration, calcurves[[i]]$area.ratio, pch = 20, col = rainbow(calibrations)[i])
    })
    abline(detrended.model)
    legend('bottomright', c(paste('calibration', as.character(y))), pch = 20, col = c(rainbow(calibrations)))
  }
  
  #organizing results
  results <- data.frame(RealConc = data$concentration)
  out <- lapply(models, function(m){
    m <- getconc(data$area.ratio, m)
    return(m)
  })
  out$joint <- getconc(data$area.ratio, joint.model)
  if(do.detrend == T){ out$detrended <- getconc(detrended$area.ratio, detrended.model) }
  out <- lapply(out, function(r){
    r * data$dilution
  })
  if(do.detrend == T){
    results <- cbind(results, out);names(results) <- c('concentration', paste('calibration', y), 'joint', 'detrended')
  }else{
    results <- cbind(results, out);names(results) <- c('concentration', paste('calibration', y), 'joint')
  }
  
  #organizing accuracy measurements
  accuracy <- apply(results[,-1], 2, function(r){
    accuracy <- abs(1 - r / results[,1]) * 100
    return(accuracy)
  })
  relative.deviation <- apply(accuracy, 2, function(r){
    sum(r[is.finite(r)], na.rm = T) 
  })
  absolute.deviation <- apply(results[,-1], 2, function(r){
    sum(abs(r - results[,1]), na.rm = T)
  })
  
  #print accuracy and results
  sink(paste('accuracy', file))
  cat('absolute deviation:\n')
  print(absolute.deviation)
  cat("\n")
  cat('relative deviation:\n')
  print(relative.deviation)
  cat("\n Complete results: \n")
  print(results)
  sink()
  
  #print a usable outputfile
  output <- data.frame(
    name = data$name
  )
  if(do.detrend == T){
    output$concentration <- results$detrended
  }else{
    output$concentration <- results$joint
  }
  output <- cbind(output, data[,c('type', names(data[,9:ncol(data)]))])
  write.table(output, file = paste('processed', file), sep = "\t", row.names = F)
  dev.off()
  setwd(mydir)
}

do <- detrend
