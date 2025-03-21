#####################################################################################
##################################### FUNCTIONS #####################################
#####################################################################################
### 0. Pre-requisits ###
initialize.analysis = function(){
  freq <- "5m"
  risk.level <- as.character(5) # alpha
  library(tidyverse)
  library(panelvar)
  #library(devtools)
  #install_github("GabauerDavid/ConnectednessApproach")
  library(ConnectednessApproach)
  library(zoo)
  library(dataPreparation)
  library(lubridate)
  library(moments)
  library(latex2exp)
  library(xtable)
  library(doParallel)
  library(quantreg)
  library(data.table)
  library(latex2exp)
  library(DescTools)
  return(list(freq=freq,risk.level=risk.level))
}
### 1. Systemic risk network functions ###
dCoVaR.network = function(freq, risk.level){
  
  ## 1.1. List of tickers
  inRiskFolder <- "./data/"
  tickers <- setNames(read.csv(paste0(inRiskFolder,"fin_tickers.txt"), sep = "\t", header = FALSE), c("ticker","name"))
  rownames(tickers) <- tickers$ticker
  # Supplement ticker data with actual company details (including Refinitiv ID data)
  inSentimentFolder <- "./data/"
  companies <- read.csv(paste0(inSentimentFolder,"MI.Companies.BASIC.04030.txt"), sep = "\t")
  for(t in tickers$ticker){
    cmpnys <- companies%>%dplyr::filter(Ticker%in%t)%>%arrange(Ticker)
    if(nrow(cmpnys)>1){
      for(i in 1:nrow(cmpnys)) cmpnys[i,"similarity"] <- adist(cmpnys$Name[i], tickers[t,"name"])
      cmpnys <- cmpnys%>%arrange(similarity)%>%dplyr::select(-similarity)
      if(t=="PNC") cmpnys <- cmpnys[-1,]
    }
    tickers[t,colnames(cmpnys)] <- cmpnys[1,]
  }
  rm(cmpnys,companies)
  rownames(tickers) <- as.character(tickers$OrgPermID)
  
  ## 1.2. List of days for which we have data
  files <- list.files(paste0(inRiskFolder,"risk/"))
  if(!file.exists("./data/days")){
    days <- numeric(0)
    for(file in files)if(grepl(freq,file)){
      print(file)
      data <- read.csv(paste0(inRiskFolder,"risk/",file))
      days <- c(days,unique(data$day))
    }
    save(days,file ="./data/days")
  }else{
    load("./data/days")
  }
  
  ## 1.3. Get/construct bank-level daily tail risk spillover measures
  if(!file.exists("./data/network")){
    network <- array(dim = c(nrow(tickers),nrow(tickers),length(days)), dimnames = list(tickers$ticker, tickers$ticker, days))
    for(file in files)if(grepl(freq,file)){
      print(file)
      data <- read.csv(paste0(inRiskFolder,"risk/",file))
      for(d in unique(data$day)){
        print(d)
        dta <- (data%>%dplyr::filter(day==d,from%in%tickers$ticker,to%in%tickers$ticker))[,c("from","to",paste0("DCoVaR.",risk.level),paste0("rDCoVaR.",risk.level))]
        dnetwork <- matrix(nrow = nrow(tickers), ncol = nrow(tickers), dimnames = list(tickers$ticker, tickers$ticker))
        dnetwork[cbind(dta$from, dta$to)] <- dta[,3]
        dnetwork[cbind(dta$to, dta$from)] <- dta[,4]
        network[,,d] <- dnetwork
      }
    }
    save(network, file = "./data/network")
  }else{
    load("./data/network")
  }
  
  ## 1.4. Get/construct bank-level daily systemic risk measures
  if(!file.exists("./data/sysrisk")){
    files <- list.files(paste0(inRiskFolder,"risk.systemic/"))
    for(file in files)if(grepl(freq,file)){
      print(file)
      data <- read.csv(paste0(inRiskFolder,"risk.systemic/",file))
      data <- (data%>%
                 dplyr::filter(from%in%tickers$ticker)%>%
                 dplyr::mutate(day=as.Date(day,"%d.%m.%Y"),ticker=from))[,c("day","ticker",paste0("DCoVaR.",risk.level),paste0("rDCoVaR.",risk.level))]
      colnames(data) <- c("day","ticker","from","to")
      if(exists("sysrisk")) sysrisk <- rbind(sysrisk,data) else sysrisk <- data
    }
    sysrisk$net <- sysrisk$from - sysrisk$to
    sysnetrisks <- data.frame(array(dim = c(length(days),nrow(tickers)), dimnames = list(format(as.Date(days,"%d.%m.%Y"),"%Y-%m-%d"), tickers$ticker)))
    for(tick in sort(unique(sysrisk$ticker))){
      print(tick)
      temp <- sysrisk%>%dplyr::filter(ticker==tick,day%in%rownames(sysnetrisks))
      sysnetrisks[format(temp$day,"%Y-%m-%d"),tick] <- temp$net
    }
    save(sysrisk, file = "./data/sysrisk")
    save(sysnetrisks, file = "./data/sysnetrisks")
  }else{
    load("./data/sysrisk")
    load("./data/sysnetrisks")
  }
  
  # 1.5. Return results
  return(list(tickers=tickers,days=days,network=network,sysrisk=sysrisk,sysnetrisks=sysnetrisks))
}
dCoVaR.network.exports = function(network, sk, export.label = "systemic_risk"){
  
  exportfilename <- paste0("./results/",export.label,".png")
  if(file.exists(exportfilename)) return()
  
  tickers <- dimnames(network)[[1]]
  days <- dimnames(network)[[3]]
  FROM <- data.frame(row.names = days)
  TO <- FROM
  NET <- FROM
  for(tick in tickers){
    aaa <- network[,tick,]
    aaa <- aaa[-which(rownames(aaa)==tick),]
    FROM[,tick] <- colSums(aaa, na.rm = T)
    aaa <- network[tick,,]
    aaa <- aaa[-which(rownames(aaa)==tick),]
    TO[,tick] <- colSums(aaa, na.rm = T)
    NET[,tick] <- TO[,tick] - FROM[,tick]
  }
  rm(aaa)
  
  ## Make plot: time dynamics of network
  if(TRUE){
    # Plot data
    days <- as.Date(days, "%d.%m.%Y")
    days[length(days)] <- as.Date("2024-01-01")
    years <- unique(format(days, "%Y"))
    year_positions <- sapply(years, function(y) min(which(format(days, "%Y") == y)))
    TOTAL <- rowSums(FROM)/100 # == rowSums(TO)
    TOTALROLL <- rollmean(TOTAL, 11, align = "center", na.pad = TRUE)
    NETROLL <- NET
    for(tick in tickers){
      NETROLL[,tick] <- rollmean(Winsorize(NET[,tick],val = quantile(NET[,tick], probs = c(0.025, 0.975), na.rm = T))/100, 11, align = "center", na.pad = TRUE)
    }
    NETPLOT <- NETROLL[,rownames(sk$sample.properties)[seq(1,nrow(sk$sample.properties),9)]]
    # Initiate file export
    png(file=exportfilename, width=8,height=4,units="in",res=960)
    par(mfrow=c(2,1),mar=c(2,2,1.5,1.5), tck = -0.02)
    # TOTAL Connectedness
    plot(days, TOTAL, xaxt = 'n', type = 'l', cex.lab=0.5, main = "TOTAL", xlab = "", ylab = "", ylim=c(min(TOTALROLL,na.rm=T),max(TOTALROLL,na.rm=T)), col = "lightgray")
    axis(1, at = days[year_positions], labels = FALSE) #axis(1, at = days[year_positions], labels = years, cex.axis = 0.6)
    text(days[year_positions], par("usr")[3]*1.15, labels = years, srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
    abline(h=0,lty=3)
    lines(days, TOTALROLL, col = "red")
    # NET Conectedness
    cols <- rgb(col2rgb(1:ncol(NETPLOT))[1,], col2rgb(1:ncol(NETPLOT))[2,], col2rgb(1:ncol(NETPLOT))[3,], alpha = 255*0.5, maxColorValue = 255)
    matplot(days, NETPLOT, xaxt = 'n',
            type = 'l', cex.lab=0.5,      # Type of plot ('l' for lines)
            lty = 1,         # Line type
            col = cols,       # Different colors for each line
            main = "NET",
            xlab = '',  # X-axis label
            ylab = '', # Y-axis label
    )
    abline(h=0,lty=3)
    axis(1, at = days[year_positions], labels = FALSE) #axis(1, at = days[year_positions], labels = years, cex.axis = 0.6)
    text(days[year_positions], par("usr")[3]-2, labels = years, srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
    legend("top", legend = colnames(NETPLOT), col = cols, lty = 1, lwd =2, cex = 0.33, horiz = TRUE)
    # Finish file export
    par(mfrow=c(1,1))
    dev.off()
  }
}
WinsorizeNet = function(net, probs = c(0.025, 0.975)){
  if(all(class(net)=="array",length(dim(net))==2)){
    for(i in dimnames(net)[[2]])if(sum(!is.na(net[,i])>1)){
      net[,i] <- Winsorize(net[,i],val = quantile(net[,i], probs = probs, na.rm = T))
    }
  }else if(all(class(net)=="array",length(dim(net))==3)){
    for(i in dimnames(net)[[1]])for(j in dimnames(net)[[2]])if(sum(!is.na(net[i,j,])>1)){
      net[i,j,] <- Winsorize(net[i,j,],val = quantile(net[i,j,], probs = probs, na.rm = T))
    }
  }else if(class(net)=="data.frame"){
    for(tick in unique(net$ticker)){
      target <- net$ticker==tick
      net[target,"from"] <- Winsorize(net[target,"from"],quantile(net[target,"from"], probs = probs, na.rm = T))
      net[target,"to"] <- Winsorize(net[target,"to"],quantile(net[target,"to"], probs = probs, na.rm = T))
    }
  }
  return(net)
}
### 2. Sentiment network functions ###
sentiment.indices = function(days, tickers, freq, freq.return = c("d", "h", "m")){
  
  ## Get/construct bank-level daily sentiment indices
  if("d"%in%freq.return){
    if(!file.exists("./data/sentiments")){
      # Get daily sentiment data for each ticker (from monthly files)
      inSentimentFolder <- "C:/Users/User/Documents/My Research/1. Working papers & idea pool/WP.9. Sentiment networks (w Caraiani & Pochea)/CMPNY_AMER/"
      files <- list.files(paste0(inSentimentFolder,"WDAI_UDAI/monthly/"))
      files <- files[-seq(1,which(grepl(".2002",files))[1]-1)]
      sentiments <- array(dim = c(length(days),nrow(tickers)), dimnames = list(as.character(days), tickers$ticker))
      for(zip_file in files){
        print(zip_file)
        # system(paste("unrar x", rar_file, "./data/"))
        fname <- paste0(inSentimentFolder,"WDAI_UDAI/monthly/",zip_file)
        file <- try(unzip(fname, list = TRUE))
        if(class(file)=="try-error") next
        unzip(fname, exdir = "./data/")
        data <- read.csv(paste0("./data/",file["Name"]), sep = "\t")%>%
          filter(assetCode%in%tickers$OrgPermID,dataType=="News_Social")%>%
          dplyr::select(id,assetCode,ticker,windowTimestamp,sentiment)
        data[,"day"] <- format(as.POSIXct(data$windowTimestamp, format = "%Y-%m-%d", tz = "UTC"), "%d.%m.%Y")
        if(sum(data[,"day"]%in%days)>0){
          data <- data%>%filter(day%in%days)
          sentiments[cbind(as.character(data$day), tickers[as.character(data$assetCode),"ticker"])] <- data$sentiment
        }
        invisible(file.remove(paste0("./data/",file["Name"])))
      }
      # Save sentiments file
      save(sentiments, file = "./data/sentiments")
    }else{
      load("./data/sentiments")
    }
    # Fill in zero-sentiment (no-news) days
    for(t in colnames(sentiments))if(sum(!is.na(sentiments[,t]))>0){
      fill <- which(!is.na(sentiments[,t]))
      to <- fill[length(fill)]
      from <- fill[1]
      fill <- which(is.na(sentiments[,t]))
      fill <- fill[fill>=from&fill<=to]
      sentiments[fill,t] <- 0
    }
    # Remove tickers with no data
    sentiments <- sentiments[,colSums(!is.na(sentiments))>0]
    # Remove days with no data
    sentiments <- sentiments[rowSums(!is.na(sentiments))>0,]
  }
  
  ## Get/construct bank-level hourly sentiment indices
  if("h"%in%freq.return){
    # Get list of raw-data files
    inSentimentFolder <- "C:/Users/User/Documents/My Research/1. Working papers & idea pool/WP.9. Sentiment networks (w Caraiani & Pochea)/CMPNY_AMER/WDAI_UHOU/1998.01-2023.06/"
    files <- list.files(inSentimentFolder)
    if(!is_empty(files))files <- files[-seq(1,which(grepl(".2002",files))[1]-1)]
    # Initialize hsentiments dataframe
    if(!file.exists("./data/hsentiments")){
      # Get data
      date_hour_combinations <- expand.grid(hour = hours(0:23), date = seq(
        as.Date(paste(year(min(as.Date(days, "%d.%m.%Y"))), month(min(as.Date(days, "%d.%m.%Y"))), "1", sep = "-")), 
        as.Date(paste(year(max(as.Date(days, "%d.%m.%Y"))), month(max(as.Date(days, "%d.%m.%Y"))), "1", sep = "-")) + months(1) - days(1), "1 day"))
      d1 <- as.POSIXct(paste(date_hour_combinations$date, date_hour_combinations$hour), tz = "UTC", format = "%Y-%m-%d %HH %MM %SS")
      d2 <- as.POSIXct(paste(date_hour_combinations$date, date_hour_combinations$hour), tz = "UTC", format = "%Y-%m-%d 0S")
      d1[is.na(d1)] <- d2[!is.na(d2)]
      hsentiments <- array(dim = c(length(d1),nrow(tickers)), dimnames = list(format(d1, "%Y-%m-%d %H:%M"), tickers$ticker))
      start.at <- 1
    }else{
      load("./data/hsentiments")
      # Last available data point
      not.na.rows <- which(rowSums(!is.na(hsentiments))>0)
      last.obs.date <- rownames(hsentiments)[not.na.rows[length(not.na.rows)]]
      if(!is_empty(files)){
        if(!is_empty(files))start.at <- which(grepl(paste0(".",year(last.obs.date), month(last.obs.date),"."),files))
        start.at <- which(grepl(paste0(".",year(last.obs.date), month(last.obs.date),"."),files))
      }else{
        start.at <- 1
      }
    }
    
    # Supplement data if needed: hourly sentiment data for each ticker (from monthly files)
    if(all(!is_empty(files),length(start.at)>0,start.at<=length(files))){
      for(i in start.at:length(files)){
        zip_file <- files[i]
        print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] Processing ",zip_file,".."))
        # system(paste("unrar x", rar_file, "./data/"))
        fname <- paste0(inSentimentFolder,zip_file)
        file <- try(unzip(fname, list = TRUE))
        if(class(file)=="try-error") next
        unzip(fname, exdir = "./data/")
        data <- read.csv(paste0("./data/",file["Name"]), sep = "\t")%>%
          filter(assetCode%in%tickers$OrgPermID,dataType=="News_Social")%>%
          dplyr::select(id,assetCode,ticker,windowTimestamp,sentiment)
        hsentiments[cbind(format(as.POSIXct(data$windowTimestamp, tz = "UTC", format = "%Y-%m-%dT%H:%M:%OSZ"), "%Y-%m-%d %H:%M"),
                          tickers[as.character(data$assetCode),"ticker"])] <- data$sentiment
        invisible(file.remove(paste0("./data/",file["Name"])))
        if(i%%12==0) save(hsentiments, file = "./data/hsentiments")
      }
      # Save hsentiments file
      save(hsentiments, file = "./data/hsentiments")
    }
    # Fill in zero-sentiment (no-news) days
    for(t in colnames(hsentiments))if(sum(!is.na(hsentiments[,t]))>0){
      fill <- which(!is.na(hsentiments[,t]))
      to <- fill[length(fill)]
      from <- fill[1]
      fill <- which(is.na(hsentiments[,t]))
      fill <- fill[fill>=from&fill<=to]
      hsentiments[fill,t] <- 0
    }
    # Remove tickers with no data
    hsentiments <- hsentiments[,colSums(!is.na(hsentiments))>0]
    # Remove days with no data
    hsentiments <- hsentiments[rowSums(!is.na(hsentiments))>0,]
  }
  
  ## Get/construct bank-level minutely sentiment indices
  if("m"%in%freq.return){
    # Get list of raw-data files
    inSentimentFolder <- "C:/Users/User/Documents/My Research/1. Working papers & idea pool/WP.9. Sentiment networks (w Caraiani & Pochea)/CMPNY_AMER/W01M_U01M/monthly/"
    files <- list.files(inSentimentFolder)
    if(!is_empty(files))files <- files[-seq(1,which(grepl(".2002",files))[1]-1)]
    outfile <- paste0("./data/m",str_replace_all(freq,"m",""),"sentiments")
    # Initialize msentiments dataframe
    if(all(!file.exists(outfile),!is_empty(files))){
      date_minute_combinations <- expand.grid(minute = minutes(seq(0,(24*60-1),as.numeric(str_replace_all(freq,"m","")))), date = seq(
        as.Date(paste(year(min(as.Date(days, "%d.%m.%Y"))), month(min(as.Date(days, "%d.%m.%Y"))), "1", sep = "-")), 
        as.Date(paste(year(max(as.Date(days, "%d.%m.%Y"))), month(max(as.Date(days, "%d.%m.%Y"))), "1", sep = "-")) + months(1) - days(1), "1 day"))
      d1 <- format(date_minute_combinations$date + date_minute_combinations$minute, "%Y-%m-%d %H:%M")
      msentiments <- array(dim = c(length(d1),nrow(tickers)), dimnames = list(d1, tickers$ticker))
      rm(d1,date_minute_combinations)
      start.at <- 1
    }else{
      load(outfile)
      # Last available data point
      not.na.rows <- which(rowSums(!is.na(msentiments))>0)
      last.obs.date <- rownames(msentiments)[not.na.rows[length(not.na.rows)]]
      start.at <- which(grepl(paste0(".",year(last.obs.date), month(last.obs.date),"."),files))
      rm(not.na.rows,last.obs.date)
    }
    
    # Supplement data if needed: hourly sentiment data for each ticker (from monthly files)
    if(all(!is_empty(files),length(start.at)>0,start.at<=length(files))){
      for(i in start.at:length(files)){
        zip_file <- files[i]
        print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] Processing ",zip_file,".."))
        fname <- paste0(inSentimentFolder,zip_file)
        file <- try(unzip(fname, list = TRUE))
        if(class(file)=="try-error") next
        unzip(fname, exdir = "./data/")
        data <- read.csv(paste0("./data/",file["Name"]), sep = "\t")%>%
          filter(assetCode%in%tickers$OrgPermID,dataType=="News_Social", !is.na(sentiment))%>%
          dplyr::select(assetCode,ticker,windowTimestamp,buzz,sentiment)
        invisible(file.remove(paste0("./data/",file["Name"])))
        data$rtime <- format(floor_date(as.POSIXct(data$windowTimestamp, tz = "UTC", format = "%Y-%m-%dT%H:%M:%OSZ"), unit=paste0(str_replace_all(freq,"m","")," mins")), "%Y-%m-%d %H:%M")
        data <- data%>%group_by(rtime,assetCode)%>%mutate(weighted_sentiment = weighted.mean(sentiment, buzz))
        # Fill in month with zeros
        DaysInMonth <- as.Date(paste(year(data$rtime[1]),month(data$rtime[1]),1:31,sep="-"))
        DaysInMonth <- DaysInMonth[!is.na(DaysInMonth)]
        TimesInMonth <- expand.grid(minute = minutes(seq(0,(24*60-1),as.numeric(str_replace_all(freq,"m","")))), date = DaysInMonth)
        TimesInMonth <- format(TimesInMonth$date + TimesInMonth$minute, "%Y-%m-%d %H:%M")
        FillPositions <- expand.grid(times = TimesInMonth, ticker = unique(tickers[as.character(data$assetCode),"ticker"]))
        msentiments[cbind(FillPositions$times,FillPositions$ticker)] <- 0
        rm(DaysInMonth,TimesInMonth,FillPositions)
        # Write nonzero sentiment data
        msentiments[cbind(data$rtime,tickers[as.character(data$assetCode),"ticker"])] <- data$weighted_sentiment
        if(i%%12==0){
          print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] Saving ",outfile," file (intermediate).."))
          save(msentiments, file = outfile)
        }
      }
      # Save msentiments file
      print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] Saving ",outfile," file (final!).."))
      save(msentiments, file = outfile)
    }
  }
  
  # Return results
  res <- list()
  if("d"%in%freq.return) res[["sentiments"]] <- sentiments
  if("h"%in%freq.return) res[["hsentiments"]] <- hsentiments
  if("m"%in%freq.return) res[["msentiments"]] <- msentiments
  return(res)
  
}
add.sentiment2sysrisk = function(sysrisk,sentiments){
  
  dys <- as.Date(rownames(sentiments),"%d.%m.%Y")
  sysrisk$sentiment <- NA
  for(cn in colnames(sentiments)){
    print(cn)
    sysrisk[sysrisk$day%in%dys & sysrisk$ticker%in%cn,"sentiment"] <- sentiments[dys%in%sysrisk$day[sysrisk$ticker==cn],cn]
  }
  sysrisk$Ds <- (sysrisk$sentiment<0)*1
  
  tkers <- unique(sysrisk$ticker)
  for(i in 1:length(tkers))
    sysrisk[sysrisk$ticker==tkers[i],"idno"] <- i
  
  dys <- unique(sysrisk$day)
  for(i in 1:length(dys))
    sysrisk[sysrisk$day==dys[i],"dayno"] <- i
  
  #plot(sysrisk$sentiment,sysrisk$from,main="FROM",xlab="",ylab="")
  summary(lm(from~I(sentiment*Ds)+I(sentiment*(1-Ds)),data=sysrisk))
  
  #plot(sysrisk$sentiment,sysrisk$to,main="FROM",xlab="",ylab="")
  summary(lm(to~I(sentiment*Ds)+I(sentiment*(1-Ds)),data=sysrisk))
  
  return(sysrisk)
} 
sentiment.sample.check = function(ssample, from = "2002-01-01 00:00", to = "2023-07-01 00:00", plt = FALSE){
  
  fromday <- which(rownames(ssample)==from)
  today <- (which(rownames(ssample)==to)-1)
  ssample <- data.frame(ssample[ifelse(is_empty(fromday),1,fromday):ifelse(is_empty(today),nrow(ssample),today),])
  ssample <- ssample%>%mutate(date=as.Date(rownames(ssample), "%Y-%m-%d %H:%M"))
  aggssample <- ssample%>%group_by(date)%>%summarise(across(everything(), ~mean(!is.na(.))))
  ssample <- ssample%>%dplyr::select(-date)
  
  # By company
  sample.properties <- data.frame()
  for(c in colnames(ssample)){
    sample.properties[c,"n"] <- nrow(ssample)
    sample.properties[c,"na"] <- round(100*sum(is.na(ssample[,c]))/nrow(ssample),2)
    sample.properties[c,"zero"] <- round(100*sum(ssample[!is.na(ssample[,c]),c]==0)/nrow(ssample),2)
    sample.properties[c,"valid"] <- round(100 - sample.properties[c,"na"] - sample.properties[c,"zero"],2)
    s <- !is.na(ssample[,c]) #s <- !is.na(ssample[,c])&ssample[,c]!=0
    sample.properties[c,"P(>0)"] <- round(100*mean(ssample[s,c]>0),2)
    sample.properties[c,"mean"] <- round(mean(ssample[s,c]),4)
    sample.properties[c,"sd"] <- round(sd(ssample[s,c]),4)
    sample.properties[c,"skew"] <- round(moments::skewness(ssample[s,c]),4)
    sample.properties[c,"kurt"] <- round(moments::kurtosis(ssample[s,c]),4)
    sample.properties[c,c("min","q1","median","q3","max")] <- round(as.numeric(quantile(ssample[s,c], probs = c(0,0.25,0.5,0.75,1))),4)
    sample.properties[c,"valid_days"] <- round(100*mean(aggssample[,c]>0),2)
  }
  sample.properties <- sample.properties%>%dplyr::arrange(desc(valid))
  sample.properties <- sample.properties%>%filter(valid>0)
  xtable::xtable(sample.properties)
  
  # Plot sentiment exposure vs. sentiment
  if(plt){
    png(file=paste0("results/sentiment_exposure.png"), width=8,height=8,units="in",res=600)
    par(mfrow=c(2,2), mar=c(4,4,4,4))
    mdl <- lm(sample.properties$`P(>0)` ~ sample.properties$valid)
    #summary(mdl)
    plot(sample.properties$valid,sample.properties$`P(>0)`, xlab="Sentiment minutes", ylab = "Positive sentiment")
    abline(a=mdl$coefficients[1],b=mdl$coefficients[2],col="red")
    abline(h=50,lty=3)
    
    mdl <- lm(sample.properties$mean ~ sample.properties$valid)
    #summary(mdl)
    plot(sample.properties$valid,sample.properties$mean, xlab="Sentiment minutes", ylab = "Mean sentiment")
    abline(a=mdl$coefficients[1],b=mdl$coefficients[2],col="red")
    abline(h=0,lty=3)
    
    mdl <- lm(sample.properties$`P(>0)` ~ sample.properties$valid_days)
    #summary(mdl)
    plot(sample.properties$valid_days,sample.properties$`P(>0)`, xlab="Sentiment days", ylab = "Positive sentiment")
    abline(a=mdl$coefficients[1],b=mdl$coefficients[2],col="red")
    abline(h=50,lty=3)
    
    mdl <- lm(sample.properties$mean ~ sample.properties$valid_days)
    #summary(mdl)
    plot(sample.properties$valid_days,sample.properties$mean, xlab="Sentiment days", ylab = "Mean sentiment")
    abline(a=mdl$coefficients[1],b=mdl$coefficients[2],col="red")
    abline(h=0,lty=3)
    par(mfrow=c(1,1))
    dev.off()
  }
  
  # Plot sentiment minutes by year
  if(plt){
    ssample <- ssample%>%mutate(year=year(as.Date(rownames(ssample), "%Y-%m-%d %H:%M")))
    aggssample <- ssample%>%group_by(year)%>%summarise(across(everything(), ~mean(!is.na(.))))
    aggssample <- aggssample[,colMeans(aggssample>0)>0]
    rns <- as.character(aggssample$year)
    aggssample <- aggssample%>%dplyr::select(-year)
    rownames(aggssample) <- rns
    png(file=paste0("results/sentiment_minutes.png"), width=6,height=6,units="in",res=600)
    plot(as.numeric(rownames(aggssample)),100*rowMeans(aggssample),type = 'b', xlab = "", ylab = "Valid sentiment minutes", ylim=c(0,max(100*rowMeans(aggssample))))
    abline(h=0,lty=3)
    dev.off()
  }
  
  # Number of sentiment-exposed stocks per day
  aaa <- round(100*rowSums(!is.na(ssample))/ncol(ssample),2)
  if(FALSE)if(plt){
    par(mfrow = c(1,2))
    plot(as.Date(names(aaa)), aaa, type = 'l', xlab = "", ylab = "")
    w <- weekdays(as.Date(names(aaa)))
    s <- !w%in%c("Saturday","Sunday")
    plot(as.Date(names(aaa[s])), aaa[s], type = 'l', xlab = "", ylab = "")
    par(mfrow = c(1,1))
  }
  print(paste0("There are ",length(unique(as.Date(names(aaa[aaa==0]))))," observations with no valid cross-sectional data: "))
  print(unique(as.Date(names(aaa[aaa==0]))))
  
  # Return results
  return(list(sample.properties = sample.properties, valid.tickers = rownames(sample.properties%>%filter(valid>0))))
}
sentiment.index.on.day = function(sentiments, day = as.Date("2020-03-15")){
  res <- sentiments[as.Date(rownames(sentiments)) == day,]
  return(res)
}
sentiment.network = function(sdata, slabel, nlag = 1, ON_split = TRUE){
  
  nfore <- 10 # slabel=="D"
  if(any(slabel=="H",slabel=="M")) nfore <- 12 * 60 / as.numeric(as.POSIXct(rownames(sdata)[2]) - as.POSIXct(rownames(sdata)[1]))
  
  ## Filter valid tickers
  check <- sentiment.sample.check(sdata)
  sdata <- sdata[,check$valid.tickers]
  till <- which(as.Date(rownames(sdata))=="2023-07-01")
  if(length(till)>0) sdata <- sdata[1:(till[1]-1),]
  sdata[is.na(sdata)] <- 0
  
  ## Full-sample sentiment connectedness network
  if(FALSE)if(is.null(mdls[[slabel]][["FS"]])){
    
    # Get sentiment network models
    if(!exists("mdls")){
      mdlfilename <- paste0("./results/SNETmdls",as.numeric(as.POSIXct(rownames(sdata)[2]) - as.POSIXct(rownames(sdata)[1])),slabel)
      if(file.exists(mdlfilename)) load(mdlfilename) else mdls <- list()
    }
    
    # Estimate and save model
    mdls[[slabel]][["FS"]] <- suppressWarnings(try(ConnectednessApproach(zoo(sdata, order.by = as.Date(rownames(sdata),"%Y-%m-%d %H:%M")), nlag = nlag, nfore = nfore, model = "LASSO", connectedness = "Time"), silent=TRUE))
    save(mdls, file = mdlfilename)
  }
  
  ## Daily sentiment networks from intraday data
  filename <- paste0("./data/s",slabel,"network")
  if(!file.exists(filename)){
    
    # Get sentiment network models
    if(!exists("mdls")){
      mdlfilename <- paste0("./results/SNETmdls",as.numeric(as.POSIXct(rownames(sdata)[2]) - as.POSIXct(rownames(sdata)[1])),slabel)
      if(file.exists(mdlfilename)) load(mdlfilename) else mdls <- list()
    }
    
    # Initialize daily sentiment network
    sdays <- as.Date(rownames(sdata), "%Y-%m-%d %H:%M")
    usdays <- unique(sdays)
    if(!is.null(mdls[[slabel]][["RW"]])){
      lastday <- as.Date(names(mdls[[slabel]][["RW"]])[length(names(mdls[[slabel]][["RW"]]))],"%d.%m.%Y")
      irange <- (which(usdays==lastday)+1):length(usdays)
      load(filename)
      if(ON_split)load(paste0(filename,"_on"))
    }else{
      irange <- 1:length(usdays)
      snetwork <- array(dim = c(nrow(tickers),nrow(tickers),length(days)), dimnames = list(tickers$ticker, tickers$ticker, days))
      if(ON_split)sonetwork <- array(dim = c(nrow(tickers),nrow(tickers),length(days)), dimnames = list(tickers$ticker, tickers$ticker, days))
    }
    
    # Estimate daily sentiment models
    for(i in irange){
      
      # Estimate sentiment network in day d
      d <- format(usdays[i], "%d.%m.%Y")
      wd <- weekdays(as.Date(d,"%d.%m.%Y"))
      if(any(wd=="Saturday",wd=="Sunday")) next else if(wd=="Monday") lg <- 3 else lg <-1
      
      #d <- usdays[i]
      if(is.null(mdls[[slabel]][["RW"]][[d]])){
        print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] Estimating sentiment network for ",d," (",wd,").."))
        if(ON_split&&i>lg){
          if(lg==3){
            Y <- sdata[sdays==usdays[i]|sdays==usdays[i-1]|sdays==usdays[i-2]|sdays==usdays[i-3],]
            hrs <- format(as.POSIXct(rownames(Y),tz="GMT"), "%A %H:%M")
            on.from <- which(hrs=="Friday 16:00")
            on.till <- which(hrs=="Monday 09:30")
            th.till <- which(hrs=="Monday 16:00")
          }else{
            Y <- sdata[sdays==usdays[i]|sdays==usdays[i-1],]
            hrs <- format(as.POSIXct(rownames(Y),tz="GMT"), "%H:%M")
            on.from <- which(hrs=="16:00")[1]
            on.till <- which(hrs=="09:30")[2]
            th.till <- which(hrs=="16:00")[2]
          }
          Yon <- Y[(on.from+1):on.till,]
          wc.on <- which_are_constant(Yon, verbose = F)
          if(length(wc.on)>0) Yon <- Yon[,-wc.on]
          Y <- Y[(on.till+1):th.till,]
        }else{
          Y <- sdata[sdays==usdays[i],]
        }
        wc <- which_are_constant(Y, verbose = F)
        if(length(wc)>0) Y <- Y[,-wc]
        mdls[[slabel]][["RW"]][[d]][["TH"]] <- suppressMessages(sentiment.network.estimate(Y,nlag,nfore))
        if(ON_split&&i>1) mdls[[slabel]][["RW"]][[d]][["ON"]] <- suppressMessages(sentiment.network.estimate(Yon,nlag,nfore))
      }
      
      # Retrieve daily sentiment networks
      if(all(!is.null(mdls[[slabel]][["RW"]][[d]][["TH"]]),d%in%days)){
        print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] Writing TH sentiment network for ",d," (",wd,").."))
        cn <- rownames(mdls[[slabel]][["RW"]][[d]][["TH"]]$CT) # colnames(mdls[[slabel]][["RW"]][[d]][["TH"]]$CT)
        snetwork[,,d] <- 0
        snetwork[cn,cn,d] <- round(100*mdls[[slabel]][["RW"]][[d]][["TH"]]$CT[,,1],2)
      }
      if(all(i>1,!is.null(mdls[[slabel]][["RW"]][[d]][["ON"]]),d%in%days)){
        print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] Writing ON sentiment network for ",d," (",wd,").."))
        cn <- rownames(mdls[[slabel]][["RW"]][[d]][["ON"]]$CT) # colnames(mdls[[slabel]][["RW"]][[d]][["ON"]]$CT)
        sonetwork[,,d] <- 0
        sonetwork[cn,cn,d] <- round(100*mdls[[slabel]][["RW"]][[d]][["ON"]]$CT[,,1],2)
      }
      
      if(day(usdays[i])==28){
        print("Saving progress..")
        save(mdls, file = mdlfilename)
        save(snetwork, file = filename)
        save(sonetwork, file = paste0(filename,"_on"))
      }
      
    }
    print("Saving results..")
    save(mdls, file = mdlfilename)
    save(snetwork, file = filename)
    save(sonetwork, file = paste0(filename,"_on"))
    beepr::beep(8)
    
  }else{
    load(filename)
    load(paste0(filename,"_on"))
  }
  
  if(ON_split){
    return(list(snetwork=snetwork,sonetwork=sonetwork))
  }else{
    return(list(snetwork=snetwork))
  }
  
}
sentiment.network.sample.select = function(){
  aaa <- data.frame(missing=names(snetwork[2,3,is.na(snetwork[2,3,])&!is.na(network[2,3,])]))
  aaa$date <- as.Date(aaa$missing,"%d.%m.%Y")
  aaa$year <- year(aaa$date)
  plot(table(aaa$year))
  
  aaa%>%filter(year==2018)
  
  aaa$dow <- weekdays(aaa$date)
  bbb <- aaa%>%filter(date>=as.Date("2012-01-01"),date<as.Date("2023-07-01"))
  plot(as.Date(names(snetwork[2,3,]),"%d.%m.%Y"),snetwork[2,3,], type = 'l')
  plot(table(bbb$year))
  
  return(list(sample.start=as.Date("2012-01-01"),sample.end=as.Date("2023-06-30")))
}
sentiment.network.estimate = function(Y,nlag,nfore,silent=TRUE){
  ret.mdl <- NULL
  while(all(class(Y)[1]=="matrix",ncol(Y)>1)){
    mdl <- suppressWarnings(try(ConnectednessApproach(zoo(Y, order.by = as.Date(rownames(Y),"%Y-%m-%d %H:%M")), nlag = nlag, nfore = nfore, model = "LASSO", connectedness = "Time"), silent=TRUE))
    if(class(mdl)!="try-error"){
      ret.mdl <- mdl
      break
    }else{
      # Remove the asset with:
      #   1) the lowest number of unique sentiment values in day, 
      #   2) lowest standard deviation of sentiment values, 
      #   3) first in list # last in list (the lowest number of unique sentiment values in all the sample)
      u.values <- apply(Y, 2, function(x) length(unique(x)))
      sd.values <- apply(Y, 2, sd)
      rmv <- which(u.values==min(u.values))
      rmv <- which(sd.values[rmv]==min(sd.values[rmv]))
      rmv <- which(names(rmv)[1]==colnames(Y))
      to.rmv <- colnames(Y)[rmv]
      Y <- Y[,-rmv]
      if(!silent)print(paste0("Eliminating `",to.rmv,"` and trying again. Y has ", ncol(Y), " columns: ", paste(colnames(Y),collapse = ",")))
    }
  }
  return(ret.mdl)
}
### 3. Network Autoregression functions ###
TriangPeerEffects = function(net){
  # net <- network
  tpe <- net
  days <- dimnames(net)[[3]]
  for(t in 1:length(days)){
    if(t%%floor(length(days)/20)==1) print(paste0("Processing day ",t,"/",length(days),".."))
    tpe[,,days[t]] <- TPE(tpe[,,t])
  }
  return(tpe)
}
TPE = function(net){
  tpe <- net
  for(i in 1:nrow(net)-1)for(j in 1:ncol(net)-1)if(all(j!=i,!is.na(net[i,j]))){
    pe <- suppressWarnings(sqrt(net[i,]*net[,j]))
    tpe[i,j] <- sum(pe,na.rm = T) / sum(is.finite(pe))
  }
  return(tpe)
}
NAR = function(network, sysrisk, snetwork, sonetwork, tickers, sfreq = "30", slabel = "M", to.do.mdls=list("net","pairwise")){
  
  if(FALSE){
    sfreq = "30"
    slabel = "M"
    to.do.mdls=list("net","pairwise")
    sysrisk_back <- sysrisk
  }
  
  ## Trim data
  print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] Loading data and preparing analysis.."))
  from <- "02.01.2002"
  to <- "30.06.2023"
  network <- network[,,which(dimnames(network)[[3]]==from):which(dimnames(network)[[3]]==to)]
  snetwork <- snetwork[,,which(dimnames(snetwork)[[3]]==from):which(dimnames(snetwork)[[3]]==to)]
  sonetwork <- sonetwork[,,which(dimnames(sonetwork)[[3]]==from):which(dimnames(sonetwork)[[3]]==to)]

  ## Winsorize Network CoVaR
  network <- WinsorizeNet(network)
  sysrisk <- WinsorizeNet(sysrisk)
  
  ## Construct/retrieve data panel required for estimating PVARs with `pvargmm`
  pvarpaneldata <- "./data/pvar_panel_data"
  if(any(!file.exists(pvarpaneldata),!file.exists(paste0(pvarpaneldata,"_net")))){
    
    reg.update <- !file.exists(pvarpaneldata)
    net.update <- !file.exists(paste0(pvarpaneldata,"_net"))
    
    # Construct triangular peer effect networks
    {
      if(file.exists("./data/tpe")){
        load("./data/tpe")
      }else{
        tpe <- TriangPeerEffects(network)
        save(tpe, file = "./data/tpe")
      }
      if(file.exists("./data/stpe")){
        load("./data/stpe")
      }else{
        stpe <- TriangPeerEffects(snetwork)
        save(stpe, file = "./data/stpe")
      }
      if(file.exists("./data/sontpe")){
        load("./data/sontpe")
      }else{
        sotpe <- TriangPeerEffects(sonetwork)
        save(sotpe, file = "./data/sontpe")
      }
    }
    
    # Get day/night mean sentimentindices for each day and ticker
    {
      dt <- as.data.table(msentiments)
      dt[, timestamp := as.POSIXct(rownames(msentiments), format = "%Y-%m-%d %H:%M", tz = "UTC")]
      dt[, day := as.Date(timestamp)]
      dt[, dow := weekdays(day)]
      dt[, interval := 
           fifelse(format(timestamp, "%H:%M") > "09:30" & format(timestamp, "%H:%M") <= "16:00" & dow %nin% c("Saturday","Sunday"),
                   paste0(as.Date(timestamp), "_day"),
                   fifelse(format(timestamp, "%H:%M") > "16:00" & dow=="Friday",
                           paste0(as.Date(timestamp + days(2) + hours(14)), "_night"),
                           fifelse(dow=="Saturday",
                                   paste0(as.Date(timestamp + days(2)), "_night"),
                                   fifelse(dow=="Sunday",
                                           paste0(as.Date(timestamp + days(1)), "_night"),
                                           paste0(as.Date(timestamp + hours(14)), "_night")
                                   )
                           )
                   )
           )
      ]
      categorical_cols <- setdiff(names(dt), c("timestamp", "interval", "day", "dow"))  # Exclude timestamp and interval
      result <- dt[, lapply(.SD, mean, na.rm=T), by = interval, .SDcols = categorical_cols]
      result[, `:=`(
        date = as.Date(sub("_.*", "", interval)),   # Extract the date portion
        period = sub(".*_", "", interval) # Extract the day/night portion
      )]
      result[, interval := NULL]
      result <- result%>%filter(date%in%days)
      result[, (categorical_cols) := lapply(.SD, function(x) fifelse(is.finite(x), x, 0)), .SDcols = categorical_cols]
    }
    
    # Arrange the data into appropriate panel form
    {
      Sys.setlocale("LC_TIME", "C")
      panel.data <- data.frame()
      panel.net.data <- data.frame()
      no <- 0
      days <- as.Date(names(network[1,1,]), "%d.%m.%Y")
      dow <- weekdays(days)
      dowmat <- matrix(data = 0, nrow = length(days), ncol = length(unique(dow)), dimnames = list(days,unique(dow)))
      for(d in unique(dow)) dowmat[dow==d,d] <- 1
      rm(dow)
      for(i in tickers$ticker){
        
        # Vs MKT data
        if(all(net.update,i%in%unique(sysrisk$ticker))){
          ssrisk <- sysrisk[sysrisk$ticker==i,]
          sysriskdays <- unique(ssrisk$day)
          rr <- (nrow(panel.net.data)+1):(nrow(panel.net.data)+length(days))
          panel.net.data[rr,"id"] <- i
          panel.net.data[rr,"idno"] <- which(tickers$ticker==i)
          panel.net.data[rr,"day"] <- days
          panel.net.data[rr,"dayno"] <- 1:length(days)
          panel.net.data[rr[days%in%sysriskdays],"RISKtoMKT"] <- ssrisk[ssrisk$day%in%days,"to"]
          panel.net.data[rr[days%in%sysriskdays],"RISKfromMKT"] <- ssrisk[ssrisk$day%in%days,"from"]
          panel.net.data[rr,"s1"] <- lag(result%>%filter(period=="day")%>%select(i))
          panel.net.data[rr,"so"] <- result%>%filter(period=="night")%>%select(i)
          rmv <- which(rownames(snetwork[,i,])==i)
          panel.net.data[rr,"soconn_to"] <- rowSums(t(snetwork[-rmv,i,]))
          panel.net.data[rr,"soconn_from"] <- rowSums(t(snetwork[i,-rmv,]))
          panel.net.data[rr,"sconn_to"] <- rowSums(t(sonetwork[-rmv,i,]))
          panel.net.data[rr,"sconn_from"] <- rowSums(t(sonetwork[i,-rmv,]))
          panel.net.data[rr,colnames(dowmat)] <- dowmat
        }
        
        # Pairwise data
        if(reg.update)for(j in tickers$ticker)if(all(j != i, sum(!is.na(snetwork[i,i,]))>0, sum(!is.na(snetwork[j,j,]))>0)){
          rr <- (nrow(panel.data)+1):(nrow(panel.data)+length(days))
          no <- no + 1
          print(paste(no, i, j, sep = "->"))
          # General record information
          panel.data[rr,"from"] <- i
          panel.data[rr,"to"] <- j
          panel.data[rr,"id"] <- paste(i, j, sep = "->")
          panel.data[rr,"idno"] <- no
          panel.data[rr,"day"] <- days
          panel.data[rr,"dayno"] <- 1:length(days)
          # Delta CoVaR
          panel.data[rr,"dcovar1"] <- lag(network[i,j,])
          panel.data[rr,"dcovar"] <- network[i,j,]
          # CoVaR Triangular Peer Effects
          panel.data[rr,"tpe1"] <- lag(tpe[i,j,])
          panel.data[rr,"tpe"] <- tpe[i,j,]
          # Sentiment
          panel.data[rr,"s1_i"] <- lag(result%>%filter(period=="day")%>%select(i))
          panel.data[rr,"s1_j"] <- lag(result%>%filter(period=="day")%>%select(j))
          panel.data[rr,"so_i"] <- result%>%filter(period=="night")%>%select(i)
          panel.data[rr,"so_j"] <- result%>%filter(period=="night")%>%select(j)
          # NET Sentiment Connectedness
          panel.data[rr,"s1_net_i"] <- lag(rowSums(t(snetwork[,i,]))-rowSums(t(snetwork[i,,])))
          panel.data[rr,"s1_net_j"] <- lag(rowSums(t(snetwork[,j,]))-rowSums(t(snetwork[j,,])))
          panel.data[rr,"so_net_i"] <- rowSums(t(sonetwork[,i,]))-rowSums(t(sonetwork[i,,]))
          panel.data[rr,"so_net_j"] <- rowSums(t(snetwork[,j,]))-rowSums(t(snetwork[j,,]))
          # Sentiment Connectedness
          panel.data[rr,"sconn1"] <- lag(snetwork[i,j,])
          panel.data[rr,"soconn"] <- sonetwork[i,j,]
          panel.data[rr,"sconn"] <- snetwork[i,j,]
          # Sentiment Triangular Peer Effects
          panel.data[rr,"stpe1"] <- lag(stpe[i,j,])
          panel.data[rr,"sotpe"] <- sotpe[i,j,]
          panel.data[rr,"stpe"] <- stpe[i,j,]
          # Day of Week
          panel.data[rr,colnames(dowmat)] <- dowmat
          # Save progress
          if(no%%500==0){
            print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] Saving data for nrow(panel.data) = ", nrow(panel.data), ".."))
            save(panel.data, file = paste0(pvarpaneldata,"-temp"))
            file.copy(paste0(pvarpaneldata,"-temp"), pvarpaneldata, overwrite = TRUE)
          }
        }
      }
      # Final save
      if(reg.update){
        #save(panel.data, file = pvarpaneldata)
        file.remove(paste0(pvarpaneldata,"-temp"))
      }
      if(net.update){
        #panel.net.data <- panel.net.data%>%filter(!is.na(RISKtoMKT),!is.na(RISKfromMKT),!is.na(sconn_to),!is.na(sconn_from))
        save(panel.net.data, file = paste0(pvarpaneldata,"_net"))
      }
    }
    
  }else{
    load(pvarpaneldata)
    load(paste0(pvarpaneldata,"_net"))
  }
  panel.data$Year <- year(as.Date(panel.data$day, "%d.%m.%Y"))
  panel.data$Month <- month(as.Date(panel.data$day, "%d.%m.%Y"))
  panel.data$YM <- paste(panel.data$Year,panel.data$Month, sep=".")

  ## Define data samples
  if(TRUE){
    # TICKERS
    oldrn <- rownames(tickers)
    rownames(tickers) <- tickers$ticker
    for(tick in tickers$ticker){
      tickers[tick,"s.na"] <- round(100*sum(is.na(snetwork[tick, tick,]))/length(snetwork[tick, tick,]),2)
      tickers[tick,"s.0"] <- round(100*sum(snetwork[tick, tick,]==0,na.rm=T)/length(snetwork[tick, tick,]),2)
      tickers[tick,"s.valid"] <- round(100*sum(!is.na(snetwork[tick, tick,])&snetwork[tick, tick,]!=0, na.rm = T)/length(snetwork[tick, tick,]),2)
    }
    rownames(tickers) <- oldrn
    rm(oldrn)
    tickers <- tickers%>%dplyr::arrange(desc(s.valid))
    samples.tick <- list(full = tickers$ticker[tickers$s.valid>1])
    
    # PERIODS OF TIME
    samples.time <- list(months = list(from = seq(as.Date("2002-01-01"), as.Date("2023-06-30"), "1 month"), to = seq(as.Date("2002-02-01"), as.Date("2023-07-01"), "1 month")-days(1)))
  }
  
  ## To,From connectedness model
  minconn <- 1/3
  wkdays <- c("Wednesday","Thursday","Friday","Monday","Tuesday")
  if("net"%in%to.do.mdls){
    for(aname in names(samples.tick))for(sname in names(samples.time)) for(s in 1:length(samples.time[[sname]]$from)){
      assets <- samples.tick[[aname]]
      samples <- samples.time[[sname]]
      sn <- paste0(format(samples$from[s],"%Y.%m"),"-",format(samples$to[s],"%Y.%m"))
      if(format(samples$from[s],"%Y.%m")==format(samples$to[s],"%Y.%m")) sn <- format(samples$from[s],"%Y.%m")
      mfilename <- paste0("./results/mdls_net/mdl",sfreq,slabel,".",aname,".",sname,".",sn)
      if(!file.exists(mfilename)){
        print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"][",s,"/",length(samples.time[[sname]]$from),"] Estimating remaining NET models for ",paste(aname,sname,sn, sep = "_"),".."))
        sample <- panel.net.data%>%dplyr::filter(id%in%assets,day>=samples$from[s],day<=samples$to[s])
        if(nrow(sample)<30){
          print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ..skipping due to insufficient data!"))
          next
        }
        sample <- sample%>%mutate(m_soconn_from = Monday*soconn_from, m_sconn_to = Monday*sconn_to)
        mdl <- suppressWarnings(try(pvargmm(
              dependent_vars = c("RISKtoMKT", "RISKfromMKT", "sconn_to", "sconn_from"),
              lags = 1,
              predet_vars = c("s1", "so", "soconn_to", "soconn_from"),
              transformation = "fod",
              data = sample,
              panel_identifier = c("idno", "dayno"),
              steps = c("onestep"),
              system_instruments = TRUE,
              system_constant = TRUE,
              max_instr_dependent_vars = 99,
              max_instr_predet_vars = 99,
              collapse = FALSE
            )))
        gc()
        if(class(mdl)=="try-error"){
          print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ..failed!"))
        }else{
          mdl_lite <- mdl
          mdl_lite$weighting_matrix_first_step <- NULL
          save(mdl, file = mfilename)
          save(mdl_lite, file = str_replace_all(mfilename,paste0(sfreq,slabel,"."),paste0(sfreq,slabel,"lite.")))
          print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ..`full` model success!"))
        }
        rm(mdl)
        gc()
      }
    }
  }
  
  # Estimate detailed PVAR model(s) for each sample
  if("pairwise"%in%to.do.mdls){
    for(aname in names(samples.tick))for(sname in names(samples.time)) 
      for(s in 1:length(samples.time[[sname]]$from)){
        assets <- samples.tick[[aname]]
        samples <- samples.time[[sname]]
        sn <- paste0(format(samples$from[s],"%Y.%m"),"-",format(samples$to[s],"%Y.%m"))
        if(format(samples$from[s],"%Y.%m")==format(samples$to[s],"%Y.%m")) sn <- format(samples$from[s],"%Y.%m")
        # Estimate remaining models for this sample
        mfilename <- paste0("./results/mdls/mdl",sfreq,slabel,".",aname,".",sname,".",sn,".full")
        if(!file.exists(mfilename)){
          # Define the data sample
          print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"][",s,"/",length(samples.time[[sname]]$from),"] Estimating remaining models for ",paste(aname,sname,sn, sep = "_"),".."))
          sample <- panel.data%>%dplyr::filter(from%in%assets,to%in%assets,day>=samples$from[s],day<=samples$to[s])
          for(wkd in wkdays) sample[is.na(sample[,wkd]),wkd] <- 0
          # Remove rarely connected directional pairs pairs from the sample
          if(sum(!is.na(sample$dcovar))>30 && sum(!is.na(sample$sconn))){
            sampleagg <- aggregate(cbind(dcovar, sconn) ~ id, sample, function(x) mean(x != 0))
            sampleagg <- sampleagg[sampleagg$dcovar<minconn|sampleagg$sconn<minconn,"id"]
            if(!is_empty(sampleagg)) sample <- sample%>%dplyr::filter(id%nin%sampleagg)
            rm(sampleagg)
          }else{
            sample <- sample[1,]
          }
          # Validate sample
          if(length(unique(sample$id))<20){
            print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ..omitting due to insufficient data!"))
            next
          }
          # Improve sample
          sample <- sample%>%mutate(m_soconn = Monday*soconn, m_sotpe = Monday*sotpe, m_stpe1 = Monday*stpe1,
                                    aso_i = abs(so_i), aso_j = abs(so_j),
                                    m_s1_i = Monday*s1_i, m_s1_j = Monday*s1_j, m_so_i = Monday*so_i, m_so_j = Monday*so_j)
          mainvars <- c("soconn", "sotpe", "stpe1", "so_i", "so_j", "s1_i", "s1_j")
          allvar <- mainvars

          # Full VAR with all variable
          if(TRUE){
            sample <- sample[,c("idno", "dayno", "sconn", "dcovar", allvar)]
            sample <- sample[rowSums(is.na(sample))==0,]
            aggssample <- sample%>%group_by(idno)%>%summarise(across(everything(), ~mean(IsZero(.))))
            aggssample <- aggssample%>%arrange(desc(sconn))
            sample <- sample%>%filter(idno%in%as.numeric(unlist(aggssample[aggssample$sconn<0.2,"idno"]))) # <0.33; 2018.08: 0.2
            mdl <- suppressWarnings(try(pvargmm(
              dependent_vars = c("sconn", "dcovar"),
              lags = 1,
              predet_vars = allvar,
              transformation = "fod",
              data = sample,
              panel_identifier = c("idno", "dayno"),
              steps = c("onestep"),
              system_instruments = TRUE,
              system_constant = TRUE,
              max_instr_dependent_vars = 99,
              max_instr_predet_vars = 99,
              collapse = FALSE
            )))
            if(class(mdl)=="try-error"){
              print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ..failed!"))
            }else{
              save(mdl, file = mfilename)
              print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ..`full` model success!"))
            }
          }
          # Save models
          rm(sample)
          rm(mdl)
          gc()
        }
        else{
          print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"][",s,"/",length(samples.time[[sname]]$from),"] Skipping models ",paste(aname,sname,sn, sep = "_")," because they are already estimated!"))
        }
      }
  }
}

### 4. Network/model plots ###
NAR.analyze = function(sfreq = "30", slabel = "M", aname = "full", sname = "months", mdltype = NULL, n = 12){
  
  if(all(!is.null(mdltype),!is.na(mdltype),mdltype=="net")){
    mdlfolder <- "./results/mdls_net/"
  }else{
    mdlfolder <- "./results/mdls/"
    mdltype <- NULL
  }
  also.irfs <- all(!is.null(n),is.finite(n))
  
  ## Get NAR regression coefficients
  resfilename <- paste0("./results/res",ifelse(any(is.null(mdltype),is.na(mdltype)),"",mdltype),sfreq, slabel,".",aname,".",sname,".csv")
  if(file.exists(resfilename)){
    res <- read.csv(file = resfilename, row.names = 1)
    colnames(res) <- str_replace(colnames(res),"\\.\\.","->")
    res$Date<-as.Date(res$Date)
  }else{
    mdlfiles <- list.files(mdlfolder,pattern = paste0("mdl",sfreq, slabel,".",aname,".",sname), full.names = T)
    res <- data.frame()
    rownum <- 0
    for(mfile in mdlfiles){
      
      mname <- unlist(strsplit(mfile,"\\."))
      if(mname[length(mname)]=="full")mname <- mname[-length(mname)]
      mname <- paste(mname[length(mname)-1],mname[length(mname)],sep=".")
      print(paste0("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] Retrieving model for ", mname, ".."))
      
      # Model coefficients
      ans <- try(load(mfile), silent=TRUE)
      if(any(inherits(ans, "try-error"),class(mdl)=="try-error")) next
      
      rownum <- rownum + 1
      res[rownum,"YM"] <- mname
      res[rownum,"n"] <- mdl$nof_observations
      res[rownum,"N"] <- length(unique(mdl[["Set_Vars"]][["category"]]))
      res[rownum,"T"] <- length(unique(mdl[["Set_Vars"]][["period"]]))
      
      parameters <- list()
      if(is.null(mdl$second_step)){
        for(i in 1:nrow(mdl$first_step))for(j in 1:(ncol(mdl$first_step)-1)){
          parameters[[rownames(mdl$first_step)[i]]][[colnames(mdl$first_step)[j]]] <- 
            list(coef = mdl$first_step[i,j], se = mdl$standard_error_first_step[i,j], t = mdl$first_step[i,j] / mdl$standard_error_first_step[i,j])
        }
      }else{
        for(i in 1:nrow(mdl$second_step))for(j in 1:(ncol(mdl$second_step)-1)){
          parameters[[rownames(mdl$second_step)[i]]][[colnames(mdl$second_step)[j]]] <- 
            list(coef = mdl$second_step[i,j], se = mdl$standard_error_second_step[i,j], t = mdl$second_step[i,j] / mdl$standard_error_second_step[i,j])
        }
      }
      res[rownum,names(unlist(parameters))] <- unlist(parameters)
      
      # Diagnostic tests
      test <- stability(mdl)
      res[rownum,"stable.condition"] <- ifelse(all(test$Modulus<1),"yes","no")
      test <- Andrews_Lu_MMSC(mdl)
      res[rownum,names(unlist(test))] <- unlist(test)
      test <- suppressWarnings(hansen_j_test(mdl))
      res[rownum,paste("HansenJ",names(unlist(test)[-5]))] <- unlist(test)[-5]
      
      # GIRF's && FEVD's
      if(also.irfs){
        # GIRFs and FEVDs
        mdl.irfs <- girf(mdl, n.ahead = n, ma_approx_steps = 20)
        mdl.fevds <- fevd_orthogonal(mdl,n)
        for(implusename in names(mdl.irfs))for(responsename in colnames(mdl.irfs[[implusename]])){
          # Cumulative impulse response
          res[rownum, paste0("girf.",implusename,"->",responsename)] <- sum(mdl.irfs[[implusename]][,responsename])
          # Minimum and maximum responses
          maxresi <- which(mdl.irfs[[implusename]][,responsename]==max(mdl.irfs[[implusename]][,responsename]))
          minresi <-which(mdl.irfs[[implusename]][,responsename]==min(mdl.irfs[[implusename]][,responsename]))
          res[rownum, paste0("girf.",implusename,"->",responsename,"maxres")] <- mdl.irfs[[implusename]][maxresi,responsename]
          res[rownum, paste0("girf.",implusename,"->",responsename,"maxres_i")] <- maxresi
          res[rownum, paste0("girf.",implusename,"->",responsename,"minres")] <- mdl.irfs[[implusename]][minresi,responsename]
          res[rownum, paste0("girf.",implusename,"->",responsename,"minres_i")] <- minresi
          # Forecast Error Variance Decomposition
          res[rownum, paste0("fevd.",implusename,"->",responsename)] <- round(100*sum(mdl.fevds[[responsename]][n,implusename]),2)
        }
      }
    }
    rm(mdl,mdl.irfs)
    gc()
    res <- res%>%mutate(Year = as.numeric(substr(YM,1,4)), Month = as.numeric(substr(YM,6,length(YM))), Day = 1, 
                        Date = ISOdate(Year,Month,Day)+months(1)-days(1))
    res <- res%>%arrange(Date)
    write.csv(res,resfilename)
  }
  
  ## Check for missing months
  allmonths <- format(seq(as.Date("2002-02-01"),as.Date("2023-06-01"),"1 month"),"%Y.%m")
  if(sum(!allmonths%in%format(res$YM))>0){
    print(paste0("Missing months: ",paste0(allmonths[!allmonths%in%format(res$YM)], collapse = ", ")))
  }else{
    print("There are no missing months for these options!")
  }
  
  ## Plot GIRFs and FEVDs
  endovarnames <- colnames(res)[grepl("->",colnames(res))&grepl("girf.",colnames(res))]
  endovarnames <- str_remove_all(endovarnames,"girf.")
  endovarnames <- unique(str_remove_all(str_remove_all(str_remove_all(str_remove_all(unique(unlist(strsplit(endovarnames,"->"))),"maxres_i"),"minres_i"),"maxres"),"minres"))
  if(also.irfs){
    
    varnames <- endovarnames
    defcol <- rep("red",nrow(res))
    neutralcol <- rgb(0.5, 0.5, 0.5, 0.35)
    valid.models <- res[,"stable.condition"]=="yes"
    
    for(responsename in varnames){
      
      # Plot GIRFs
      nr <- ceiling(length(varnames)/2)
      png(file=paste0("results/GIRFs_",mdltype,slabel,"_",aname,"_",responsename,".png"), width=8,height=4*nr,units="in",res=30*8*nr) # 60x8=480 pixels per row
      par(mfrow=c(nr,2), mar = c(2.5,2.5,3.5,2.5))
      for(implusename in varnames){
        
        plotdates <- res[valid.models,"Date"]
        plotcol <- paste0("girf.",implusename,"->",responsename)
        plotdata <- res[valid.models,plotcol]
        plotdatamin <- res[valid.models,paste0(plotcol,"minres")]
        plotdatamax <- res[valid.models,paste0(plotcol,"maxres")]
        
        # Baseline plot
        cols <- defcol
        cols[plotdata>=0] <- "forestgreen"
        plot(plotdates, plotdata, pch=21, cex = 0.4, lwd=0.25, col="gray", bg=cols, 
             ylim=c(min(0,as.numeric(quantile(c(plotdata,plotdatamin),0))),max(0,as.numeric(quantile(c(plotdata,plotdatamax),1)))), xlab="",ylab=""
             , main = TeX(paste0(varname2label(implusename, F), " impluse,     ", varname2label(responsename, T), " response")), cex.main = 1,
        )
        abline(h=0,lty=3)
        
        # Vertical lines
        for (i in 1:length(plotdates)){
          segments(plotdates[i], plotdatamin[i], plotdates[i], plotdatamax[i], col = neutralcol)
        }
        
        # Min and Max responses
        pchs <- rep(25,length(plotdates))
        cols <- rep(neutralcol,length(plotdates))
        pchs[res[valid.models,paste0(plotcol,"minres_i")] < res[valid.models,paste0(plotcol,"maxres_i")]] <- 22
        cols[res[valid.models,paste0(plotcol,"minres_i")] < res[valid.models,paste0(plotcol,"maxres_i")]] <- rgb(1, 1, 1, 0.33)
        points(plotdates, plotdatamin, pch=pchs, cex = 0.25, col=neutralcol, bg = cols)
        
        pchs <- rep(24,length(plotdates))
        cols <- rep(neutralcol,length(plotdates))
        pchs[res[valid.models,paste0(plotcol,"minres_i")] > res[valid.models,paste0(plotcol,"maxres_i")]] <- 22
        cols[res[valid.models,paste0(plotcol,"minres_i")] > res[valid.models,paste0(plotcol,"maxres_i")]] <- rgb(1, 1, 1, 0.33)
        points(plotdates, plotdatamax, pch=pchs, cex = 0.25, col=neutralcol, bg = cols)
        
        # Cummulative responses
        cols <- defcol
        cols[plotdata>=0] <- "forestgreen"
        points(plotdates, plotdata, pch=21, cex = 0.4, lwd=0.25, col="gray", bg=cols) # , main=TeX(paste0("$",str_replace(plotcol,"->","\\\\to "),"$"))
        
      }
      par(mfrow=c(1,1), mar = c(5, 4, 4, 2) + 0.1)
      dev.off()
      
      # Plot FEVDs
      png(file=paste0("results/FEVDs_",mdltype,slabel,"_",aname,"_",responsename,".png"), width=8,height=4*nr,units="in",res=30*8*nr)
      par(mfrow=c(nr,2), mar = c(2.5,2.5,2.5,2.5))
      for(implusename in varnames){
        
        plotcol <- paste0("fevd.",implusename,"->",responsename)
        plotdata <- res[,plotcol]
        cols <- rep(rgb(0.5, 0.5, 0.5, 0.75),nrow(res))
        for(i in 1:nrow(res)) cols[i] <- rgb(0,0,0, 0.2+(plotdata[i]/100)*0.6)
        plot(res$Date, plotdata, pch=21, cex = 0.4, col=neutralcol, bg=cols, ylim=c(min(0,plotdata),max(0,plotdata)), xlab="",ylab=""
             #, main=TeX(paste0(varname2label(implusename,F)," $\\to$ ",varname2label(responsename,T))), cex.main = 0.75
             , main=TeX(paste0("Share of ", varname2label(implusename,F)," in ",varname2label(responsename,T))), cex.main = 1
        )#, type="o")
        abline(h=0,lty=3)
        abline(h=100,lty=3)
        abline(h=mean(plotdata),lty=2,col=rgb(1,0,0,0.5))#rgb(0,0,0, 0.2+(mean(plotdata)/100)*0.6))
      }
      par(mfrow=c(1,1), mar = c(5, 4, 4, 2) + 0.1)
      dev.off()
      
    }
    
  }
  
  ## Plot results
  golden_ratio <- (1+sqrt(5))/2
  varnames <- gsub("fd_","",gsub("fod_","",gsub("lag1_","",colnames(res))))
  varnames <- gsub(".coef","",varnames[grepl(".coef",varnames)])
  varnames <- unique(unlist(strsplit(varnames,".",fixed = T)))
  fno <- length(varnames)
  # AR coefficients with error bands
  if(TRUE){
    for(dep_var in endovarnames){
      #png(file=paste0("results/coefs_",slabel,"_",aname,".png"), width=8,height=8,units="in",res=600)
      #par(mfrow=c(fno,fno))
      for(indep_var in varnames){
        cn <- paste0("fd_",dep_var,".fd_lag1_",indep_var)
        if(!paste0(cn,".coef")%in%colnames(res)) cn <- paste0("fod_",dep_var,".fod_lag1_",indep_var)
        if(!paste0(cn,".coef")%in%colnames(res)) cn <- paste0("fod_",dep_var,".fod_",indep_var)
        if(paste0(cn,".coef")%in%colnames(res)){
          png(file=paste0("results/coefs/coef_",mdltype,slabel,"_",aname,"_",dep_var,"-",indep_var,".png"), width=4*golden_ratio,height=4,units="in",res=480)
          par(mar = c(2.5,2.5,3.5,2.5))
          plot.with.errors(x = res$Date,
                           avg = res[,paste0(cn,".coef")],
                           sdev = res[,paste0(cn,".se")],
                           #paste0(indep_var,"_{t-1}"), dep_var,
                           "",TeX(paste0("Coefficient of ", varname2label(indep_var,F)," on ",varname2label(dep_var,T))),
                           only.significant.coefs = FALSE, point.size = 0.35, arrow.size = 0.2)
          par(mfrow=c(1,1), mar = c(5, 4, 4, 2) + 0.1)
          dev.off()
        }
      }
      #par(mfrow=c(1,1))
      #dev.off()
    }
  }
  # AR t stat of coefficients
  if(FALSE){
    png(file=paste0("results/tstats_",slabel,"_",aname,".png"), width=8,height=8,units="in",res=600)
    par(mfrow=c(fno,fno))
    for(dep_var in varnames)for(indep_var in varnames){
      cn <- paste0("fd_",dep_var,".fd_lag1_",indep_var)
      plot.t.stats(x = res$Date,
                   ts = res[,paste0(cn,".t")],
                   paste0(indep_var,"_{t-1}"), dep_var,
                   only.significant.coefs = FALSE)
      
    }
    par(mfrow=c(1,1))
    dev.off()
  }
  # Custom plots
  if(FALSE){
    dep_var <- "dcovar"
    indep_var <- "sconn"
    cn <- paste0("fd_",dep_var,".fd_lag1_",indep_var)
    if(paste0(cn,".coef")%in%colnames(res)){
      png(file=paste0("results/coefs_",slabel,"_",aname,"_s-r.png"), width=8,height=8,units="in",res=600)
      plot.with.errors(x = res$Date,
                       avg = res[,paste0(cn,".coef")],
                       sdev = res[,paste0(cn,".se")],
                       paste0(indep_var,"_{t-1}"), dep_var,
                       only.significant.coefs = TRUE)
      dev.off()
      png(file=paste0("results/coefs_",slabel,"_",aname,"_s-r_t.png"), width=8,height=8,units="in",res=600)
      plot.t.stats(x = res$Date,
                   ts = res[,paste0(cn,".t")],
                   paste0(indep_var,"_{t-1}"), dep_var,
                   only.significant.coefs = FALSE)
      dev.off()
    }
  }
  # IRF surfaces
  if(FALSE){
    
    mdl.irsf <- girf(mdl, n.ahead = 10, ma_approx_steps = 10)
    mdl.irsf <- bootstrap_irf(mdl, typeof_irf = "OIRF", n.ahead = 10, nof_Nstar_draws = 100, confidence.band = 0.95, mc.cores = 1)
    par(mfrow = c(2,2))
    plot(1:nrow(mdl.irsf$sconn), mdl.irsf$sconn[,"sconn"], type = "l", col = "red", xlab = "", ylab = "", main = "S->S")
    plot(1:nrow(mdl.irsf$sconn), mdl.irsf$dcovar[,"sconn"], type = "l", col = "red", xlab = "", ylab = "", main = "dCoVaR-?>S")
    plot(1:nrow(mdl.irsf$sconn), mdl.irsf$sconn[,"dcovar"], type = "l", col = "red", xlab = "", ylab = "", main = "S-?>dCoVaR")
    plot(1:nrow(mdl.irsf$sconn), mdl.irsf$dcovar[,"dcovar"], type = "l", col = "red", xlab = "", ylab = "", main = "dCoVaR->dCoVaR")
    par(mfrow = c(2,2))
    
  }
  # FVED surfaces???
  if(FALSE){
    mdl.fved <- fevd_orthogonal(mdl, n.ahead = 10)
  }
  # Table of coefficients
  if(FALSE){
    rownames(res) <- res$Date
    res <- res%>%dplyr::select(grepl(varnames,colnames(res)))
    
    res2 <- res[,grep(paste(c(paste0("lag1_",varnames),"stable.condition","HansenJ p.value"),collapse="|"),colnames(res), value=TRUE)]
    for(cn in colnames(res2))if(grepl(".coef",cn)){
      res2[,paste0(cn,".stars")] <- coef.with.stars(res2[,cn],res2[,gsub(".coef",".t",cn)],df=res$n)
    }
    res2 <- res2[,grepl(".stars",colnames(res2))|grepl("stable.condition",colnames(res2))|grepl("HansenJ",colnames(res2))]
    #res2 <- res2[res2$stable.condition=="yes"&as.numeric(res2$`HansenJ p.value`)>0.05,-c(1,2)]
    res2 <- res2[month(as.Date(rownames(res2)))==1,-c(1,2)]
    rownames(res2) <- format(as.Date(rownames(res2)),"%Y %b")
    
    print(xtable::xtable(res2, type = "latex", align = paste(rep("l",ncol(res2)+1),collapse = ""), display = rep("s",ncol(res2)+1),collapse = ""))
    
  }
  
  # Return model results
  return(res)
}
varname2label = function(varname, endog = TRUE){
  
  # NET varnames
  if(endog){
    if(varname=="RISKtoMKT") return("$\\delta^{i\\to}_d$")
    if(varname=="RISKfromMKT") return("$\\delta^{\\to i}_d$")
    if(varname=="sconn_to") return("$\\theta^{i\\to}_{TH,d}$")
    if(varname=="sconn_from") return("$\\theta^{\\to i}_{TH,d}$")
  }else{
    if(varname=="RISKtoMKT") return("$\\delta^{i\\to}_{d-1}$")
    if(varname=="RISKfromMKT") return("$\\delta^{\\to i}_{d-1}$")
    if(varname=="sconn_to") return("$\\theta^{i\\to}_{TH,d-1}$")
    if(varname=="sconn_from") return("$\\theta^{\\to i}_{TH,d-1}$")
  }
  if(varname=="so") return("$ONS^i_{d}$")
  if(varname=="s1") return("$THS^i_{d-1}$")
  if(varname=="soconn_to") return("$\\theta^{i\\to}_{ON,d}$")
  if(varname=="soconn_from") return("$\\theta^{\\to i}_{ON,d}$")
  
  # PAIRWISE varnames
  if(endog){
    if(varname=="dcovar") return("$\\delta^{i\\to j}_d$")
    if(varname=="sconn") return("$\\theta^{i\\to j}_{TH,d}$")
  }else{
    if(varname=="dcovar") return("$\\delta^{i\\to j}_{d-1}$")
    if(varname=="sconn") return("$\\theta^{i\\to j}_{TH,d-1}$")
  }
  if(varname=="soconn") return("$\\theta^{i\\to j}_{ON,d}$")
  if(varname=="sotpe") return("$\\pi^{i\\to j}_{ON,d}$")
  if(varname=="stpe1") return("$\\pi^{i\\to j}_{TH,d-1}$")
  if(varname=="so_i") return("$ONS^i_{d}$")
  if(varname=="so_j") return("$ONS^j_{d}$")
  if(varname=="s1_i") return("$THS^i_{d-1}$")
  if(varname=="s1_j") return("$THS^j_{d-1}$")
  
}
plot.with.errors = function(x, avg, sdev, impulse, response, only.significant.coefs = FALSE, point.size = 0.5, arrow.size = 0.25){
  
  
  #x <- as.Date(x)
  years <- unique(format(x, "%Y"))
  year_positions <- sapply(years, function(y) min(which(format(x, "%Y") == y)))
  
  # Validation
  keep <- is.finite(avg)&is.finite(sdev)
  
  # Color scheme
  ub <- avg+1.96*sdev
  lb <- avg-1.96*sdev
  clrs <- rep("darkgray",length(x)) # clrs <- rep("darkgray",length(x))
  clrs[ub<0] <- "red"
  clrs[lb>0] <- "forestgreen"
  significant <- numeric(length(x))
  significant[keep&ub<0] <- -1
  significant[keep&lb>0] <- 1
  if(only.significant.coefs) keep <- keep&abs(significant)==1
  
  # Plot
  plot(x, rep(0,length(x)),xaxs="i", xaxt="n", col = "White", type = "l", xlab = "", ylim = c(min(lb[keep], na.rm = T),max(ub[keep], na.rm = T)), ylab = TeX(paste0("$",impulse,"$")), main = response)
  
  axis(1, at = x[year_positions], labels = FALSE) #axis(1, at = x[year_positions], labels = years, cex.axis = 0.6)
  text(x[year_positions], par("usr")[3]*1.15, labels = years, srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
  
  abline(v=x[month(x)==1], lty = 3, col = rgb(105,105,105,0.25*255,maxColorValue = 255))
  abline(h = 0, col = "DarkGrey", lty=3)
  arrows(x[keep], lb[keep], x[keep], ub[keep], length=0.025, angle=90, code=3, col = clrs[keep], lwd = arrow.size)
  if(!only.significant.coefs)points(x[keep&significant==0], avg[keep&significant==0], pch=21, col = clrs[keep&significant==0], bg = clrs[keep&significant==0], cex = point.size)
  points(x[keep&significant==1], avg[keep&significant==1], pch=21, col = clrs[keep&significant==1], bg = clrs[keep&significant==1], cex = point.size)
  points(x[keep&significant==-1], avg[keep&significant==-1], pch=21, col = clrs[keep&significant==-1], bg = clrs[keep&significant==-1], cex = point.size)
  
}
plot.t.stats = function(x, ts, impulse, response, only.significant.coefs = FALSE){
  
  # Validation
  keep <- is.finite(ts)
  
  # Color scheme
  clrs <- rep("black",length(x))
  clrs[ts<= -1.96] <- "red"
  clrs[ts>= 1.96] <- "forestgreen"
  if(only.significant.coefs) keep <- keep&abs(ts)>=1.96
  
  # Plot
  plot(x, rep(0,length(x)),xaxs="i", col = "White", type = "l", xlab = "", ylim = c(min(ts[keep], na.rm = T),max(ts[keep], na.rm = T)), ylab = TeX(paste0("$",impulse,"$")), main = response)
  abline(h = 0, col = "DarkGrey")
  abline(h = c(-1.96,1.96),lty=3, col = "DarkGrey")
  points(x[keep], ts[keep], pch=22, col = "DarkGrey", bg = clrs[keep], cex = 1)
  
  # Vertical year delimiters
  abline(v=x[month(x)==12], lty = 3, col = "DarkGrey")
  
}
coef.with.stars = function(coefs,ts,df=Inf){
  stars <- ts
  stars[abs(ts)<qt(.95,df)] <- ""
  stars[abs(ts)>=qt(.95,df)] <- "^{*}"
  stars[abs(ts)>=qt(.975,df)] <- "^{**}"
  stars[abs(ts)>=qt(.995,df)] <- "^{***}"
  return(paste0("$",round(coefs,4),stars,"$"))
}
### 5. Sample descriptions and others ###
sample.description = function(tickers, snetwork){
  
  # Add sentiment (level) statistics
  tickers <- add.summary.stats(tickers, sentiments, "s")
  
  # Add sentiment network statistics: frequency, net properties
  tickers <- add.network.stats(tickers, sHnetwork, "sn")
  
  # Add systemic risk statistics
  tickers <- add.summary.stats(tickers, sysnetrisks, "r")
  
  # Add Network-DCoVaR statistics
  #tickers <- add.network.stats(tickers, network, "rn") # Repair network (Compute Delta VaR's and store on main diagonal!!)
  
  
  # Select/Arrange table + Export in LaTeX
  rownames(tickers) <- 1:nrow(tickers)
  tbl <- xtable(tickers%>%mutate(vld=100-s.0)%>%select(ticker,name,vld), auto = TRUE)
  print(tbl, include.rownames = T)
  
  ## (Point) Figures showing sentiment statistics vs. SysRisk statistics
  # Sentiment mean vs. Net SysRisk mean
  png(file=paste0("results/SentimentMeanVsNetSysRiskMean.png"), width=8,height=8,units="in",res=600)
  plot.relationship(tickers%>%dplyr::select(s.mu,r.mu),c(TeX("$\\hat{\\mu}_{Sentiment}$"),TeX("$\\hat{\\mu}_{SysRisk}$")))
  dev.off()
  png(file=paste0("results/SentimentVsNetSysRisk.png"), width=8,height=8,units="in",res=600)
  par(mfrow=c(2,2))
  plot.relationship(tickers%>%dplyr::select(s.mu,r.mu),c(TeX("$\\hat{\\mu}_{Sentiment}$"),TeX("$\\hat{\\mu}_{SysRisk}$")))
  plot.relationship(tickers%>%dplyr::select(s.mu,r.sd),c(TeX("$\\hat{\\mu}_{Sentiment}$"),TeX("$\\hat{\\sigma}_{SysRisk}$")))
  plot.relationship(tickers%>%dplyr::select(s.sd,r.mu),c(TeX("$\\hat{\\sigma}_{Sentiment}$"),TeX("$\\hat{\\mu}_{SysRisk}$")))
  plot.relationship(tickers%>%dplyr::select(s.sd,r.sd),c(TeX("$\\hat{\\sigma}_{Sentiment}$"),TeX("$\\hat{\\sigma}_{SysRisk}$")))
  par(mfrow=c(1,1))
  dev.off()
  
  ## Add sentiment sample properties
  mu <- colMeans(sentiments, na.rm=T)
  hist(mu, breaks = 50)
  
  mean(sentiments, na.rm=T)
  sd(sentiments, na.rm=T)
  quantile(sentiments, probs = seq(0, 1, 0.1), na.rm=T)
  moments::skewness(sentiments, na.rm=T)
  moments::kurtosis(sentiments, na.rm=T)
  
}
add.summary.stats = function(tickers, data, lbl){
  oldrn <- rownames(tickers)
  rownames(tickers) <- tickers$ticker
  new <- colSums(!is.na(data))
  tickers[names(new),paste0(lbl,".n")] <- new
  tickers[names(new),paste0(lbl,".f")] <- round(100*new/nrow(data),2)
  new <- colMeans(data, na.rm=T)
  tickers[names(new),paste0(lbl,".mu")] <- round(new,3)
  new <- apply(data, 2, sd, na.rm = T)
  tickers[names(new),paste0(lbl,".sd")] <- round(new,3)
  new <- apply(data, 2, quantile, probs = seq(0, 1, 0.1), na.rm = T)
  tickers[colnames(new),gsub("%","",paste0(lbl,".q",rownames(new)))] <- round(t(new),3)
  new <- apply(data, 2, moments::skewness, na.rm = T)
  tickers[names(new),paste0(lbl,".skew")] <- round(new,3)
  new <- apply(data, 2, moments::kurtosis, na.rm = T)
  tickers[names(new),paste0(lbl,".kurt")] <- round(new,3)
  rownames(tickers) <- oldrn
  return(tickers)
}
add.network.stats = function(tickers, data, lbl){
  oldrn <- rownames(tickers)
  rownames(tickers) <- tickers$ticker
  for(tick in tickers$ticker){
    tickers[tick,paste0(lbl,".na")] <- round(100*sum(is.na(data[tick, tick,]))/length(data[tick, tick,]),2)
    tickers[tick,paste0(lbl,".0")] <- round(100*sum(data[tick, tick,]==0,na.rm=T)/length(data[tick, tick,]),2)
    tickers[tick,paste0(lbl,".valid")] <- round(100*sum(!is.na(data[tick, tick,])&data[tick, tick,]!=0, na.rm = T)/length(data[tick, tick,]),2)
  }
  rownames(tickers) <- oldrn
  return(tickers)
}
plot.relationship = function(data,lbls){
  if(ncol(data)!=2) stop("Data does not have 2 columns!")
  g <- c(mean(data[,1], na.rm = T),mean(data[,2], na.rm = T))
  mdl <- lm(data[,2]~data[,1])
  plot(data[,1], data[,2], xlab = lbls[1], ylab = lbls[2],pch=19,cex=0.75,col="Gray")
  abline(a = mdl$coefficients[1], b = mdl$coefficients[2], col = "red")
  abline(h=0,v=0,lty=3,col="DarkGray")
  points(g[1],g[2],col="red",pch=19,cex=1.25)
  text(g[1],g[2], pos = 3, offset=1, labels=TeX(paste0("\\hat{\\beta} = ",round(mdl$coefficients[2],2)," (",round(summary(mdl)$coefficients[2,"Pr(>|t|)"],4),")")), col="red", cex=0.75)
}
### 6. Results vs. market movements
add_spx_returns = function(res){
  #install.packages("quantmod")
  library(quantmod)
  getSymbols("^GSPC", src = "yahoo", from = "1900-01-01", to = Sys.Date())
  sp500_monthly <- to.monthly(GSPC, indexAt = "lastof", OHLC = FALSE)
  plot(index(GSPC), coredata(GSPC)[,"GSPC.Adjusted"], type = "l", col = "blue", xlab = "", ylab = "", log = "y")
  sp500_log_returns <- diff(log(Cl(sp500_monthly))) * 100
  res$spx <- as.numeric(coredata(sp500_log_returns[res$Date]))
  return(res)
}
plot_with_reg = function(x,y){
  plot(x,y,xlab="",ylab="")
  abline(h=0,lty=3)
  mdl <- lm(y ~ x+I(x^2)+I(x^3))
  abline(a=coefficients(mdl)[1],b=coefficients(mdl)[2],lty=3,lwd=2,col="red")
  points(x,mdl$fitted.values, col = "blue", cex=0.5)
  print(summary(mdl))
}
#####################################################################################


#####################################################################################
#################################### CODE TO RUN ####################################
#####################################################################################
list2env(initialize.analysis(), .GlobalEnv)
list2env(dCoVaR.network(freq,risk.level), .GlobalEnv)
sk <- sentiment.sample.check(sysnetrisks)
dCoVaR.network.exports(network, sk, "systemic_risk")
list2env(sentiment.indices(days,tickers,"30m",c("d","m","h")), .GlobalEnv) #list2env(sentiment.indices(days,tickers,"15m",c("d","m","h")), .GlobalEnv)
sysrisk <- add.sentiment2sysrisk(sysrisk,sentiments)
list2env(sentiment.network(sdata=msentiments, slabel="M", nlag=1, ON_split=TRUE), .GlobalEnv)
sk <- sentiment.sample.check(msentiments)
dCoVaR.network.exports(snetwork, sk, "sentiment_network")
dCoVaR.network.exports(sonetwork, sk, "sentiment_onetwork")
rm(sk)
###########NAR(network, sysrisk, snetwork, sonetwork, tickers, to.do=list("net"))
sample.description(tickers, sHnetwork)
NAR.analyze(slabel = "M", aname = "above75", n = NULL)
resnet <- NAR.analyze(slabel = "M", aname = "full", mdltype = "net", n = 12)
resnet <- add_spx_returns(resnet)
res <- NAR.analyze(slabel = "M", aname = "full", mdltype = NULL, n = 12)
res <- add_spx_returns(res)
