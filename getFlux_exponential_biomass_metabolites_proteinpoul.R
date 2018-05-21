require("plyr")

##correct for background trends
correctMedium <- function(data, medium){
  #Make and plot a linear line through background samples
  correction <- lm(medium[['concentration']] ~ as.numeric(medium[['time']]))
  plot(data[['concentration']] ~ as.numeric(data[['time']]), ylab = 'concentration', xlab = 'time', pch = 20, main = 'before correction')
  points(medium[['concentration']] ~ as.numeric(medium[['time']]), ylab = 'concentration', xlab = 'time', pch = 20, col = 'red')
  abline(correction)
  legend('bottomright', c('samples', 'medium'), pch = 20, col = c('black', 'red'))
  #Correct all data, and recalculate and plot a linear line through background samples
  dc <- apply(data, 1, function(r){
    r['concentration'] <- as.numeric(r['concentration']) - correction$coefficients[[2]] * as.numeric(r['time'])
  })
  mc <- apply(medium, 1, function(r){
    r['concentration'] <- as.numeric(r['concentration']) - correction$coefficients[[2]] * as.numeric(r['time'])
  })
  data['concentration'] <- dc
  medium['concentration'] <- mc
  plot(data[['concentration']] ~ as.numeric(data[['time']]), ylab = 'concentration', xlab = 'time', pch = 20, main = 'after correction')
  points(medium[['concentration']] ~ as.numeric(medium[['time']]), ylab = 'concentration', xlab = 'time', pch = 20, col = 'red')
  correction <- lm(medium[['concentration']] ~ as.numeric(medium[['time']]))
  abline(correction)
  legend('bottomright', c('samples', 'medium'), pch = 20, col = c('black', 'red'))
  return(data)
}



getFlux <- function(dir = getwd(), files = list.files(dir), inputfiles = files[grepl('processed', files)], biomassfile = files[grepl('biomass', files, ignore.case = T)], metabolites = NULL, outputdir = dir, medium = 'medium', outputname = 'rates.txt', exp = 1, background.correction = F, pool.biomass.tzero = F){
  setwd(dir)
  #inputfiles <- files[grepl('processed', files)]
  biomass <- read.csv(biomassfile, sep = "\t", header = T, stringsAsFactors = F, comment.char = "#")
  biomass <- biomass[!biomass[,1]==""|is.na(biomass[,1]),]
  
  # determine the 
  if(is.null(metabolites)){
    metabolites <- gsub(".txt", "", unlist(lapply((strsplit(inputfiles, ' ')), tail, 1)))
  }
  all.m <- unique(metabolites)
  
  #organize all the data
  metabolite.data <- as.list(1:length(all.m))
  all.data <- as.list(1:length(inputfiles))
  all.data <- lapply(all.data, function(l){
    data <- read.csv(inputfiles[l], sep = "\t", header = TRUE, dec = ".", na.strings = "-----", stringsAsFactors = F, comment.char = "#")
    l <- data[!(is.na(data$time)|data$time == 'NA'),]
  })
  all.data <- lapply(all.data, function(l){
    keep <- intersect(names(l), names(biomass))
    l <- l[,keep]
  })
  log <- sapply(all.m, function(m){
    grepl(m, metabolites)
  })
  if(length(all.m)==1){
    metabolite.data <- lapply(metabolite.data, function(m){
      m <- do.call('rbind', all.data)
      return(m)})
  } else {
    metabolite.data <- lapply(metabolite.data, function(m){
      m <- do.call('rbind', all.data[log[,m]])
      return(m)
    })
  }
  
  setwd(outputdir)
  #ly is used to lapply through the metabolite.data
  ly <- as.list(1:length(metabolite.data))
  
  
  sample.data <- lapply(ly, function(i){
    l <- metabolite.data[[i]]
    l <- l[!l['type']==medium,]
    return(l)
  })
  
  ########## Getting all the factor combinations##########
  niveaus <- lapply(ly, function(i) {
    n <- intersect(names(sample.data[[i]]), names(biomass))
    n <- n[-c(which('name' == n), which('time'== n), which('concentration'==n))]
    return(n)
  })
  value <- lapply(ly, function(i){
    value <- unique(sample.data[[i]][,niveaus[[i]]])
    return(value)
  })
  
  ##### correct data for trends in medium #######
  if(background.correction == T){
    metabolite.data <- lapply(ly, function(i){
      m <- metabolite.data[[i]]
      v <- value[[i]]
      if(is.null(dim(v))){
        mset <- m
        cset <- m[m['type'] == medium,]
        
        pdf(paste('background correction ', all.m[i], '.pdf', sep = ""))
        mset <- correctMedium(mset, cset)
        dev.off()
      
        return(mset)
      } else {
        v <- v[-which(names(v) == 'type')]
        v <- unique(v)
        log <- apply(v, 1, function(r){
          apply(m, 1, function(M){
            log <- ifelse(F %in% (r %in% M), F, T)
            return(log)
          })
        })
        out <- lapply(as.list(1:dim(log)[2]), function(p){
          l <- log[,p]
          mset <- m[l,]
          cset <- mset[mset['type']== medium,]
        
          pdf(paste('background correction ', all.m[i],' ', paste(v[p,], collapse = " "), '.pdf', collapse = ""))
          out <- correctMedium(mset, cset)
          dev.off()    
          return(out)
        
        })
        out <- ldply(out, data.frame)
        return(out)
      }
    })
  }
  #############
  
  
  # Value now gives all the possible combination of factor levels
  # The issue is that value can be a list of data.frames, or a list of vectors 
  
  log <- lapply(ly, function(i){
    m <- metabolite.data[[i]]
    v <- value[[i]]
    if(is.null(dim(v))){
      sapply(v, function(V){
        apply(m, 1, function(M){
          return(ifelse(FALSE %in% (V %in% M) , F, T))
        })
      })
    }else{
      apply(v, 1, function(V){
        apply(m, 1, function(M){
          return(ifelse(FALSE %in% (V %in% M) , F, T))
        })
      })}    
  })
  # we also want to include samples that are medium at t = 0, 
  # and that adhere to the other factors included by user
  remaining.factors <- lapply(ly, function(i){
    remaining.factors <- setdiff(names(metabolite.data[[i]]), c('name', 'concentration', 'type', 'time'))
    return(remaining.factors)
  })
  
  value2 <- lapply(ly, function(i){
    v <- value[[i]]
    if(is.null(dim(v))){
      q <- data.frame(type = rep(medium, length(v)),
                      time = rep(0, length(v)))
      return(q)
    }else{
      q <- as.data.frame(v[,which(names(v) == remaining.factors[[i]])])
      q['type'] <- rep(medium, dim(v)[1])
      q['time'] <- rep(0, dim(v)[1])
      names(q) <- c(remaining.factors[[i]], 'type', 'time')
      return(q)
    }
  })
  
  log2 <- lapply(ly, function(i){
    m <- metabolite.data[[i]]
    v <- value2[[i]]
    apply(v, 1, function(V){
      apply(m, 1, function(M){
        return(ifelse(FALSE %in% (V %in% M) , F, T))
      })
    })        
  })
  
  ###########combining the logs, and selecting the data###########
  nlog <- lapply(ly, function(i){
    nlog <- log[[i]]|log2[[i]]
    return(nlog)
  })
  
  dm <- lapply(ly, function(i){
    nlog <- nlog[[i]]
    m <- metabolite.data[[i]]
    apply(nlog, 2, function(c){
      m[c,]
    })
  })
  
  ################organising data per factor combination: biomass#########
  plog <- lapply(ly, function(i){
    p <- biomass
    v <- value[[i]]
    if(is.null(dim(v))){
      sapply(v, function(V){
        apply(p, 1, function(P){
          return(ifelse(FALSE %in% (V %in% P), F, T))
        })
      })
    }else{
      apply(v, 1, function(V){
        apply(p, 1, function(P){
          return(ifelse(FALSE %in% (V %in% P) , F, T))
        })
      })
    }
  })
  # If we want to pool biomass data at time 0, we want to also add those values at t0 of same cell lines.
  if(pool.biomass.tzero == T){
    pvalue2 <- lapply(ly, function(i){
      v <- value[[i]]
      f <- remaining.factors[[i]]
      if(length(f)==0){
        v <- data.frame(
          type = v,
          time = rep(0, length(v)))
      }else{
        v <- v[,-which(names(v)==f)]
        v <- data.frame(
          type = v,
          time = rep(0, length(v)))
      }
      return(v)
    })
  
    plog2 <- lapply(ly, function(i){
      p <- biomass
      v <- pvalue2[[i]]
      apply(v, 1, function(V){
        apply(p, 1, function(P){
          P <- gsub(" ", "", P)
          return(ifelse(FALSE %in% (V %in% P) , F, T))
        })
      })
    })
    plog <- lapply(ly, function(i){
      nlog <- plog[[i]]|plog2[[i]]
      return(nlog)
    })
  }
  
  dp <- lapply(ly, function(i){
    plog <- plog[[i]]
    p <- biomass
    apply(plog, 2, function(c){
      p[c,]
    })
  })
  
  conditions <- lapply(ly, function(i){
    as.list(1:length(dp[[i]]))
  })
  
  ##getting the models
  lines <- lapply(ly, function(i){
    lapply(conditions[[i]], function(set){
      v <- value[[i]]
      if(is.null(dim(v))){
        pdf(paste('rate calculation', all.m[i], v[set], '.pdf'))
      } else {
        pdf(paste(paste(c('rate calculation', all.m[i], v[set,]), collapse = " "),'.pdf',collapse = ""))
      }
      P <- dp[[i]][[set]]
      M <- dm[[i]][[set]]
      q <- get_l(P, M, exp)
      dev.off()
      return(q)
    })
  })
  ##getting the actual parameters
  q <- lapply(ly, function(i){
    lapply(conditions[[i]], function(set){
      v <- value[[i]]
      
      P <- dp[[i]][[set]]
      M <- dm[[i]][[set]]
      q <- get_k(P, M, exp)
      
      return(q)
    })
  })
  
  
  # get p values by comparing all lines with each other within metabolite
  pvalues <- lapply(ly, function(i){
    results <- c()
    mlines <- lines[[i]]
    set <- conditions[[i]]
    lapply(set, function(subset){
      current.line <- (mlines[[subset]])
      while(length(mlines) > subset){
        subset <- subset + 1
        test.line <- mlines[[subset]]
        B1 <- summary(current.line)$coef[[2]]
        B2 <- summary(test.line)$coef[[2]]
        Sb1 <- summary(current.line)$coef[[4]]
        Sb2 <- summary(test.line)$coef[[4]]
        tvalue <- -abs(B1 - B2) / sqrt(Sb1^2 + Sb2^2)
        df <- current.line$df + test.line$df
        results <<- c(results, 2 * pt(tvalue, df))
        
      }
    })
    return(results)
  })
  
  
  #Organizing the data to get readable output
  all.combos <- sapply(ly, function(i){
    combo <- c()
    v <- value[[i]]
    set <- conditions[[i]]
    sapply(set, function(subset) {
      if(is.null(nrow(v))){
        combo <- c(combo, (paste(all.m[i], paste(v, collapse = " - "))))
      }else{
        currentv <- paste(v[subset,], collapse = " & ")
        while(nrow(v) > subset){
          combo <- c(combo, paste(all.m[i], currentv, '-', paste(v[(subset+1),], collapse=" & "), collapse = " "))
          subset <- subset + 1
        }
      }
      return(combo)
    })
  })
  p.results <- data.frame(set = unique(unlist(all.combos)), p = unlist(pvalues))
  p.results['p<0.05'] <- ifelse(p.results$p < 0.05, '*', '')
  
  n <- lapply(ly, function(i){
    lapply(conditions[[i]], function(set){
      v <- value[[i]]
      if(is.null(dim(v))){
        n <- (paste(all.m[i], v[set]))
      } else {
        n <- (paste(c(all.m[i], v[set,]), collapse = " "))
      }
      return(n)
    })
  })
  
  
  sink(outputname, type = 'output')
  dummy <- lapply(ly, function(i) {
    lapply(conditions[[i]], function(set){
      cat(n[[i]][[set]])
      cat('\n')
      cat(q[[i]][[set]])
      cat('\n')
    })
  })
  print(p.results)
  sink()
  
  
  
}

get_k <- function(p, m, exp) {
  p <- p[with(p,order(time)),]
  if(exp == 0){
    line11 <- lm(concentration ~ time, data = p)
    line12 <- nls(concentration ~ p0 * exp(mu * time), data = p, start = list(p0 = mean(p$concentration), mu = 0))
    if(sum(resid(line11)^2) <= sum(resid(line12)^2)){
      line1 <- line11
      ip <- as.numeric(m[,'time']) * (line1[[1]][1] + 0.5 * line1[[1]][2] * as.numeric(m[,'time']) )
    } else {
      line1 <- line12
      ip <- ((summary(line1)[[10]][[1]] * exp(summary(line1)[[10]][[2]] * as.numeric(m$time)))/summary(line1)[[10]][[2]])
    }
  }else if(exp == 1){
    line1 <- lm(concentration ~ time, data = p)
    ip <- as.numeric(m[,'time']) * (line1[[1]][1] + 0.5 * line1[[1]][2] * as.numeric(m[,'time']) )
  }else if(exp == 2){
    line1 <- nls(concentration ~ p0 * exp(mu * time), data = p, start = list(p0 = mean(p$concentration), mu = 0))
    ip <- ((summary(line1)[[10]][[1]] * exp(summary(line1)[[10]][[2]] * as.numeric(m$time)))/summary(line1)[[10]][[2]])
  }
  line2 <- lm(m[,'concentration'] ~ ip)
  k <- line2$coefficients[[2]]
  error <- coefficients(summary(line2))['ip' ,'Std. Error']
  return(c(k, error, line2$df))
}


get_l <- function(p, m, exp) {
  p <- p[with(p,order(time)),]
  plot(concentration ~ time, data = m, pch = 20, xlab = 'time', ylab = 'metabolite concentration', main = 'yield')
  plot(concentration ~ time, data = p, pch = 20, xlab = 'time', ylab = 'biomass', main = 'growth over time')
  if(exp == 0){
    line11 <- lm(concentration ~ time, data = p)
    line12 <- nls(concentration ~ p0 * exp(mu * time), data = p, start = list(p0 = mean(p$concentration), mu = 0))
    lines(p$time, predict(line11)); lines(p$time, predict(line12))
    if(sum(resid(line11)^2) <= sum(resid(line12)^2)){
      line1 <- line11
      ip <- as.numeric(m[,'time']) * (line1[[1]][1] + 0.5 * line1[[1]][2] * as.numeric(m[,'time']) )
      exp <- 1
    } else {
      line1 <- line12
      ip <- ((summary(line1)[[10]][[1]] * exp(summary(line1)[[10]][[2]] * as.numeric(m$time)))/summary(line1)[[10]][[2]])
      exp <- 2
    }
    lines(p$time, predict(line1), col = 'red')
    legend('bottomright', c('linear', 'exponential'), col = c(if(exp == 1){ c('red', 'black')}else{ c('black', 'red')}), pch = 20)
  }else if(exp == 1){
    line1 <- lm(concentration ~ time, data = p)
    lines(p$time, predict(line1))
    ip <- as.numeric(m[,'time']) * (line1[[1]][1] + 0.5 * line1[[1]][2] * as.numeric(m[,'time']) )
  }else if(exp == 2){
    line1 <-nls(concentration ~ p0 * exp(mu * time), data = p, start = list(p0 = mean(p$concentration), mu = 0))
    lines(p$time, predict(line1))
    ip <- ((summary(line1)[[10]][[1]] * exp(summary(line1)[[10]][[2]] * as.numeric(m$time)))/summary(line1)[[10]][[2]])
  }
  plot(ip, m[,'concentration'], xlab = 'integrated biomass', ylab = 'metabolite concentration', pch = 20)
  line2 <- lm(m[,'concentration'] ~ ip)
  abline(line2)
  
  return(line2)
}

