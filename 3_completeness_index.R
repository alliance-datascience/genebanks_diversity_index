################################################################################
library(countrycode)
library(ggplot2)
library(factoextra)
library(glue)
library(plsdepot)
library(ggrepel)
library(MetBrewer)
################################################################################
################################################################################

completeness_func <- function(outdir,collection_name,PCDI_df,QGI_df,TQI_df){
  message(paste("reading  previous results for ",collection_name))
  #load(paste0(outdir, "/", collection_name, "_subsets_1.RData"))
  subsets <- readRDS(paste0(
    outdir,
    "/",
    collection_name,
    "/",
    collection_name,
    "_subsets_new_1.RDS"
  ))
  passport_data_w <- subsets[[1]]
  passport_data_l <- subsets[[2]]
  passport_data_h <- subsets[[3]]
  df_total <- subsets[[4]]
  df_total_nonNA <- subsets[[5]]
  df_total_NA <- subsets[[6]]
  x_total <- subsets[[7]]
  x2_total <- subsets[[8]]
  x2_non_na <- subsets[[9]]
  df_wild <- subsets[[10]]
  df_landrace <- subsets[[11]]
  df_hybrid <- subsets[[12]]
  ################################################################################
  message("Joining PDCI, GQI, and TQI to wild, landraces and hybrids subsets")
  message("This step is really slow... Take it easy")
  #Each measurement is added to the passport data used as basis
  
  ################################################################################
  #crop wild relatives
  message("1. Joining PDCI, GQI, and TQI to wild")
  if(nrow(passport_data_w)>0){
    passport_data_w$PDCI <- NA
    passport_data_w$GQI <- NA
    passport_data_w$TQI <- NA
    for(i in 1:nrow(passport_data_w)){
      #i <- 1
      x <- PCDI_df$PDCI_final[
        which(PCDI_df$ACCENUMB==passport_data_w$ACCENUMB[[i]])]
      if(length(x)>0){
        passport_data_w$PDCI[[i]] <- PCDI_df$PDCI_final[
          which(PCDI_df$ACCENUMB==passport_data_w$ACCENUMB[[i]])]
      } else {
        passport_data_w$PDCI[[i]] <- NA
      };rm(x)
      x <- QGI_df$SCORE[
        which(QGI_df$ACCENUMB==passport_data_w$ACCENUMB[[i]])]
      if(length(x)>0){
      passport_data_w$GQI[[i]] <- QGI_df$SCORE[
        which(QGI_df$ACCENUMB==passport_data_w$ACCENUMB[[i]])]
      } else {
        passport_data_w$GQI[[i]] <- NA
      };rm(x)
      x <- TQI_df$TAXA_SCORE_VALUE[
        which(TQI_df$ACCENUMB==passport_data_w$ACCENUMB[[i]])]
      if(length(x)>0){
      passport_data_w$TQI[[i]] <- TQI_df$TAXA_SCORE_VALUE[
        which(TQI_df$ACCENUMB==passport_data_w$ACCENUMB[[i]])]
      } else {
        passport_data_w$TQI[[i]] <- NA
      };rm(x)
    };rm(i)
  } else {
    warning("NO WILD DATA AVAILABLE")
  }
  ################################################################################
  #landraces
  message("2. Joining PDCI, GQI, and TQI to landraces")
  if(nrow(passport_data_l)>0){
    passport_data_l$PDCI <- NA
    passport_data_l$GQI <- NA
    passport_data_l$TQI <- NA
    for(i in 1:nrow(passport_data_l)){
      #i <- 1
      x <- PCDI_df$PDCI_final[
        which(PCDI_df$ACCENUMB==passport_data_l$ACCENUMB[[i]])]
      if(length(x)>0){
        passport_data_l$PDCI[[i]] <- PCDI_df$PDCI_final[
          which(PCDI_df$ACCENUMB==passport_data_l$ACCENUMB[[i]])]
      } else {
        passport_data_l$PDCI[[i]] <- NA
      };rm(x)
      x <- QGI_df$SCORE[
        which(QGI_df$ACCENUMB==passport_data_l$ACCENUMB[[i]])]
      
      if(length(x)>0){
        passport_data_l$GQI[[i]] <- QGI_df$SCORE[
        which(QGI_df$ACCENUMB==passport_data_l$ACCENUMB[[i]])]
      } else {
        passport_data_l$GQI[[i]] <- NA
      };rm(x)
      x <- TQI_df$TAXA_SCORE_VALUE[
        which(TQI_df$ACCENUMB==passport_data_l$ACCENUMB[[i]])]
      if(length(x)>0){
      passport_data_l$TQI[[i]] <- TQI_df$TAXA_SCORE_VALUE[
        which(TQI_df$ACCENUMB==passport_data_l$ACCENUMB[[i]])]
      } else {
        passport_data_l$TQI[[i]] <- NA
      };rm(x)
    };rm(i)
  } else {
    warning("NO LANDRACE DATA AVAILABLE")
  }
  ################################################################################
  #hybrids
  message("3. Joining PDCI GQI, and TQI to hybrids")
  if(nrow(passport_data_h)>0){
    passport_data_h$PDCI <- NA
    passport_data_h$GQI <- NA
    passport_data_h$TQI <- NA
    for(i in 1:nrow(passport_data_h)){
      #i <- 1
      x <- PCDI_df$PDCI_final[
        which(PCDI_df$ACCENUMB==passport_data_h$ACCENUMB[[i]])]
      if(length(x)>0){
        passport_data_h$PDCI[[i]] <- PCDI_df$PDCI_final[
          which(PCDI_df$ACCENUMB==passport_data_h$ACCENUMB[[i]])]
      } else {
        passport_data_h$PDCI[[i]] <- NA
      };rm(x)
      x <- QGI_df$SCORE[
        which(QGI_df$ACCENUMB==passport_data_h$ACCENUMB[[i]])]
      if(length(x)>0){
        passport_data_h$GQI[[i]] <- QGI_df$SCORE[
          which(QGI_df$ACCENUMB==passport_data_h$ACCENUMB[[i]])]
      } else {
        passport_data_h$GQI[[i]] <- NA
      };rm(x)
      x <- TQI_df$TAXA_SCORE_VALUE[
        which(TQI_df$ACCENUMB==passport_data_h$ACCENUMB[[i]])]
      if(length(x)>0){
        passport_data_h$TQI[[i]] <- TQI_df$TAXA_SCORE_VALUE[
          which(TQI_df$ACCENUMB==passport_data_h$ACCENUMB[[i]])]
      } else {
        passport_data_h$TQI[[i]] <- NA
      };rm(x)

    };rm(i)
  } else {
    warning("NO HYBRID DATA AVAILABLE")
  }
  ################################################################################
  #spps approach
  message("Summarizing results per species")
  ################################################################################
  #wild
  #obtaining species for wild
  if(nrow(passport_data_w)>0){
    spp_w <- unique(passport_data_w$NAME)
    #creating data frame for mean and sd for PDCI, GQI, and TQI
    completeness_wild_df <- as.data.frame(matrix(nrow = length(spp_w),ncol=9))
    colnames(completeness_wild_df) <- c(
      "Taxon","status","records","PDCI_mean","PDCI_sd","GQI_mean","GQI_sd","TQI_mean","TQI_sd"
    )
    
    for(i in 1:nrow(completeness_wild_df)){
      #i <- 1
      x_i <- passport_data_w[which(passport_data_w$NAME==spp_w[[i]]),]
      completeness_wild_df[i,1] <- spp_w[[i]]
      completeness_wild_df[i,2] <- "wild"
      completeness_wild_df[i,3] <- nrow(x_i)
      completeness_wild_df[i,4] <- mean(x_i$PDCI,na.rm = T)
      completeness_wild_df[i,5] <- sd(x_i$PDCI,na.rm = T)
      completeness_wild_df[i,6] <- mean(x_i$GQI,na.rm = T)
      completeness_wild_df[i,7] <- sd(x_i$GQI,na.rm = T)
      completeness_wild_df[i,8] <- mean(x_i$TQI,na.rm = T)
      completeness_wild_df[i,9] <- sd(x_i$TQI,na.rm = T)
    };rm(i)
    
    
    completeness_wild_df$PDCI_mean <- completeness_wild_df$PDCI_mean/10
    completeness_wild_df$GQI_mean <- completeness_wild_df$GQI_mean/12
    completeness_wild_df$TQI_mean <- completeness_wild_df$TQI_mean/3
    
  }
  ################################################################################
  #landrace
  if(nrow(passport_data_l)>0){
    spp_l <- unique(passport_data_l$NAME)
    #creating data frame for mean and sd for PDCI, GQI, and TQI
    completeness_land_df <- as.data.frame(matrix(nrow = length(spp_l),ncol=9))
    colnames(completeness_land_df) <- c(
      "Taxon","status","records","PDCI_mean","PDCI_sd","GQI_mean","GQI_sd","TQI_mean","TQI_sd"
    )
    for(i in 1:nrow(completeness_land_df)){
      #i <- 1
      x_i <- passport_data_l[which(passport_data_l$NAME==spp_l[[i]]),]
      completeness_land_df[i,1] <- spp_l[[i]]
      completeness_land_df[i,2] <- "landrace"
      completeness_land_df[i,3] <- nrow(x_i)
      completeness_land_df[i,4] <- mean(x_i$PDCI,na.rm = T)
      completeness_land_df[i,5] <- sd(x_i$PDCI,na.rm = T)
      completeness_land_df[i,6] <- mean(x_i$GQI,na.rm = T)
      completeness_land_df[i,7] <- sd(x_i$GQI,na.rm = T)
      completeness_land_df[i,8] <- mean(x_i$TQI,na.rm = T)
      completeness_land_df[i,9] <- sd(x_i$TQI,na.rm = T)
    };rm(i)
    
    #standarizing measurment to take values from 0 to 1
    completeness_land_df$PDCI_mean <- completeness_land_df$PDCI_mean/10
    completeness_land_df$GQI_mean <- completeness_land_df$GQI_mean/12
    completeness_land_df$TQI_mean <- completeness_land_df$TQI_mean/3
    
  }
  ################################################################################
  #hybrid
  if(nrow(passport_data_h)>0){
    spp_h <- unique(passport_data_h$NAME)
    #creating data frame for mean and sd for PDCI, GQI, and TQI
    completeness_hybd_df <- as.data.frame(matrix(nrow = length(spp_h),ncol=9))
    colnames(completeness_hybd_df) <- c(
      "Taxon","status","records","PDCI_mean","PDCI_sd","GQI_mean","GQI_sd","TQI_mean","TQI_sd"
    )
    for(i in 1:nrow(completeness_hybd_df)){
      #i <- 1
      x_i <- passport_data_h[which(passport_data_h$NAME==spp_h[[i]]),]
      completeness_hybd_df[i,1] <- spp_h[[i]]
      completeness_hybd_df[i,2] <- "hybrid"
      completeness_hybd_df[i,3] <- nrow(x_i)
      completeness_hybd_df[i,4] <- mean(x_i$PDCI,na.rm = T)
      completeness_hybd_df[i,5] <- sd(x_i$PDCI,na.rm = T)
      completeness_hybd_df[i,6] <- mean(x_i$GQI,na.rm = T)
      completeness_hybd_df[i,7] <- sd(x_i$GQI,na.rm = T)
      completeness_hybd_df[i,8] <- mean(x_i$TQI,na.rm = T)
      completeness_hybd_df[i,9] <- sd(x_i$TQI,na.rm = T)
    };rm(i)
    
    #standarizing measurment to take values from 0 to 1
    completeness_hybd_df$PDCI_mean <- completeness_hybd_df$PDCI_mean/10
    completeness_hybd_df$GQI_mean <- completeness_hybd_df$GQI_mean/12
    completeness_hybd_df$TQI_mean <- completeness_hybd_df$TQI_mean/3
  }
  ################################################################################
  ################################################################################
  #joining data in one data frame
  
  if(exists("completeness_wild_df") & 
     exists("completeness_land_df") &
     exists("completeness_hybd_df")){
    completeness_df <- rbind(completeness_wild_df,completeness_land_df)
    completeness_df <- rbind(completeness_df,completeness_hybd_df)  
  } else if(exists("completeness_wild_df") & 
            exists("completeness_land_df") &
            !exists("completeness_hybd_df")){
    completeness_df <- rbind(completeness_wild_df,completeness_land_df)
  } else if(exists("completeness_wild_df") & 
            !exists("completeness_land_df") &
            exists("completeness_hybd_df")){
    completeness_df <- rbind(completeness_wild_df,completeness_hybd_df)
  }
  
  
  #################################
  #Doing weigthed mean using records proportion in the collection (CWR)
  #w <- completeness_df$records/sum(completeness_df$records)
  if(exists("completeness_wild_df")){
    PDCI_wild_W <- weighted.mean(x = completeness_wild_df$PDCI_mean,
                                 w = completeness_wild_df$records/
                                   sum(completeness_df$records) ,
                                 na.rm = T)
    GQI_wild_W <- weighted.mean(x = completeness_wild_df$GQI_mean,
                                w = completeness_wild_df$records/
                                  sum(completeness_df$records) ,
                                na.rm = T)
    TQI_wild_W <- weighted.mean(x = completeness_wild_df$GQI_mean,
                                w = completeness_wild_df$records/
                                  sum(completeness_df$records) ,
                                na.rm = T)
  } else {
    PDCI_wild_W <- NA
    GQI_wild_W <- NA
    TQI_wild_W <- NA
  }
  ###############################
  #Doing weigthed mean using records proportion in the collection (landraces)
  if(exists("completeness_land_df")){
    PDCI_land_W <- weighted.mean(x = completeness_land_df$PDCI_mean,
                                 w = completeness_land_df$records/
                                   sum(completeness_df$records) ,
                                 na.rm = T)
    GQI_land_W <- weighted.mean(x = completeness_land_df$GQI_mean,
                                w = completeness_land_df$records/
                                  sum(completeness_df$records) ,
                                na.rm = T)
    TQI_land_W <- weighted.mean(x = completeness_land_df$TQI_mean,
                                w = completeness_land_df$records/
                                  sum(completeness_df$records) ,
                                na.rm = T)
  } else {
    PDCI_land_W <- NA
    GQI_land_W <- NA
    TQI_land_W <- NA
  }
  ###############################
  #Doing weigthed mean using records proportion in the collection (hybrids)
  if(exists("completeness_hybd_df")){
    PDCI_hybd_W <- weighted.mean(x = completeness_hybd_df$PDCI_mean,
                                 w = completeness_hybd_df$records/
                                   sum(completeness_df$records) ,
                                 na.rm = T)
    GQI_hybd_W <- weighted.mean(x = completeness_hybd_df$GQI_mean,
                                w = completeness_hybd_df$records/
                                  sum(completeness_df$records) ,
                                na.rm = T)
    TQI_hybd_W <- weighted.mean(x = completeness_hybd_df$TQI_mean,
                                w = completeness_hybd_df$records/
                                  sum(completeness_df$records) ,
                                na.rm = T)
  } else {
    PDCI_hybd_W <- NA
    GQI_hybd_W <- NA
    TQI_hybd_W <- NA
  }
  ###############################
  #Doing weigthed mean using records proportion in the collection (TOTAL)
  
  PDCI_total_W <- weighted.mean(x = completeness_df$PDCI_mean[
    completeness_df$status %in% c("wild","landrace")],
    w = completeness_df$records[
      completeness_df$status %in% c("wild","landrace")]/
      sum(completeness_df$records),
    na.rm = T)
  GQI_total_W <- weighted.mean(x = completeness_df$GQI_mean[
    completeness_df$status %in% c("wild","landrace")],
    w = completeness_df$records[
      completeness_df$status %in% c("wild","landrace")]/
      sum(completeness_df$records),
    na.rm = T)
  
  TQI_total_W <- weighted.mean(x = completeness_df$TQI_mean[
    completeness_df$status %in% c("wild","landrace")],
    w = completeness_df$records[
      completeness_df$status %in% c("wild","landrace")]/
      sum(completeness_df$records),
    na.rm = T)
  
  
  #preparing summary table
  df_completeness_final <- as.data.frame(matrix(nrow=4,ncol=4)
  )
  colnames(df_completeness_final) <- c("status","PDCI","GQI","TQI")
  df_completeness_final[,1] <- c("wild","landrace","hybrid","total")
  #CWR
  df_completeness_final[1,2] <-  PDCI_wild_W
  df_completeness_final[1,3] <-  GQI_wild_W
  df_completeness_final[1,4] <-  TQI_wild_W
  #LANDRACES
  df_completeness_final[2,2] <-  PDCI_land_W
  df_completeness_final[2,3] <-  GQI_land_W
  df_completeness_final[2,4] <-  TQI_land_W
  #HYBRIDS
  df_completeness_final[3,2] <-  PDCI_hybd_W
  df_completeness_final[3,3] <-  GQI_hybd_W
  df_completeness_final[3,4] <-  TQI_hybd_W
  #SUBSET
  df_completeness_final[4,2] <-  PDCI_total_W
  df_completeness_final[4,3] <-  GQI_total_W
  df_completeness_final[4,4] <-  TQI_total_W
  
  message("Calculating completeness index")
  
  df_completeness_final$final_score <- NA
  df_completeness_final$final_score <- 
    
    unlist(lapply(1:nrow(df_completeness_final), function(i){
      
      x <- exp(mean(log(c(df_completeness_final$PDCI[[i]],
                          df_completeness_final$GQI[[i]],
                          df_completeness_final$TQI[[i]]))))
      return(x)
    }))
  
  
    #exp(mean(log(c(df_completeness_final$PDCI,df_completeness_final$GQI,df_completeness_final$TQI))))
    #df_completeness_final$PDCI*df_completeness_final$GQI*df_completeness_final$TQI
  ################################################################################
  ################################################################################
  ################################################################################
  #Getting countries
  countries <- colnames(df_total)
  ################################################################################
  message("Calculating completeness indexes per country (wild)")
  if(exists("spp_w")){
    pdci_country_w <- as.data.frame(matrix(nrow = length(spp_w),ncol=length(countries)))
    gqi_country_w <- as.data.frame(matrix(nrow = length(spp_w),ncol=length(countries)))
    tqi_country_w <- as.data.frame(matrix(nrow = length(spp_w),ncol=length(countries)))
    
    #colnames for PDCI, GQI, and TQI
    colnames(pdci_country_w) <- countries
    row.names(pdci_country_w) <- spp_w
    
    colnames(gqi_country_w) <- countries
    row.names(gqi_country_w) <- spp_w
    
    colnames(tqi_country_w) <- countries
    row.names(tqi_country_w) <- spp_w
    
    
    #calculating PDCI, GQI, and TQI in one loop
    for(i in 1:length(countries)){
      #i <- 3
      x_i <- passport_data_w[
        which(passport_data_w$ORIGCTY==countries[[i]]),]
      if(nrow(x_i)==0){
        pdci_country_w[,i] <- NA
        gqi_country_w[,i] <- NA
        tqi_country_w[,i] <- NA
      } else {
        for(j in 1:length(spp_w)){
          #j <- 1
          x_i_j <- x_i[which(x_i$NAME==spp_w[[j]]),]
          if(nrow(x_i_j)>0){
            pdci_country_w[j,i] <- mean(x_i_j$PDCI,na.rm = T)
            gqi_country_w[j,i] <- mean(x_i_j$GQI,na.rm = T)
            tqi_country_w[j,i] <- mean(x_i_j$TQI,na.rm = T)
          } else {
            pdci_country_w[j,i] <- NA
            gqi_country_w[j,i] <- NA 
            tqi_country_w[j,i] <- NA 
          }
        };rm(j,x_i,x_i_j)
      }
    };rm(i)
    pdci_country_w <- pdci_country_w/10
    gqi_country_w <- gqi_country_w/12
    tqi_country_w <- tqi_country_w/3
  } else {
    pdci_country_w <- NA
    gqi_country_w <- NA
    tqi_country_w <- NA
  }
  ################################################################################
  message("Calculating completeness indexes per country (landraces)")
  
  if(exists("spp_l")){
    pdci_country_l <- as.data.frame(matrix(nrow = length(spp_l),ncol=length(countries)))
    gqi_country_l <- as.data.frame(matrix(nrow = length(spp_l),ncol=length(countries)))
    tqi_country_l <- as.data.frame(matrix(nrow = length(spp_l),ncol=length(countries)))
    
    
    #colnames for PDCI, GQI, and TQI
    
    colnames(pdci_country_l) <- countries
    row.names(pdci_country_l) <- spp_l
    
    colnames(gqi_country_l) <- countries
    row.names(gqi_country_l) <- spp_l
    
    colnames(tqi_country_l) <- countries
    row.names(tqi_country_l) <- spp_l
    
    for(i in 1:length(countries)){
      #i <- 3
      x_i <- passport_data_l[
        which(passport_data_l$ORIGCTY==countries[[i]]),]
      if(nrow(x_i)==0){
        pdci_country_l[,i] <- NA
        gqi_country_l[,i] <- NA
        tqi_country_l[,i] <- NA
      } else {
        for(j in 1:length(spp_l)){
          #j <- 1
          x_i_j <- x_i[which(x_i$NAME==spp_l[[j]]),]
          if(nrow(x_i_j)>0){
            pdci_country_l[j,i] <- mean(x_i_j$PDCI,na.rm = T)
            gqi_country_l[j,i] <- mean(x_i_j$GQI,na.rm = T)
            tqi_country_l[j,i] <- mean(x_i_j$TQI,na.rm = T)
            
          } else {
            pdci_country_l[j,i] <- NA
            gqi_country_l[j,i] <- NA   
            tqi_country_l[j,i] <- NA   
          }
        };rm(j,x_i,x_i_j)
      }
    };rm(i)
    pdci_country_l <- pdci_country_l/10
    gqi_country_l <- gqi_country_l/12
    tqi_country_l <- tqi_country_l/3
    
  } else {
    pdci_country_l <- NA
    gqi_country_l <- NA
    tqi_country_l <- NA
  }
  ################################################################################
  if(exists("spp_h")){
    message("Calculating completeness indexes per country (hybrids)")
    pdci_country_h <- as.data.frame(matrix(nrow = length(spp_h),ncol=length(countries)))
    gqi_country_h <- as.data.frame(matrix(nrow = length(spp_h),ncol=length(countries)))
    tqi_country_h <- as.data.frame(matrix(nrow = length(spp_h),ncol=length(countries)))
    
    #colnames for PDCI, GQI, and TQI
    
    colnames(pdci_country_h) <- countries
    row.names(pdci_country_h) <- spp_h
    
    colnames(gqi_country_h) <- countries
    row.names(gqi_country_h) <- spp_h
    
    colnames(tqi_country_h) <- countries
    row.names(tqi_country_h) <- spp_h
    
    for(i in 1:length(countries)){
      #i <- 3
      x_i <- passport_data_h[
        which(passport_data_h$ORIGCTY==countries[[i]]),]
      if(nrow(x_i)==0){
        pdci_country_h[,i] <- NA
        gqi_country_h[,i] <- NA
        tqi_country_h[,i] <- NA
      } else {
        for(j in 1:length(spp_h)){
          #j <- 1
          x_i_j <- x_i[which(x_i$NAME==spp_h[[j]]),]
          if(nrow(x_i_j)>0){
            pdci_country_h[j,i] <- mean(x_i_j$PDCI,na.rm = T)
            gqi_country_h[j,i] <- mean(x_i_j$GQI,na.rm = T)
            tqi_country_h[j,i] <- mean(x_i_j$TQI,na.rm = T)
            
          } else {
            pdci_country_h[j,i] <- NA
            gqi_country_h[j,i] <- NA  
            tqi_country_h[j,i] <- NA  
          }
        };rm(j,x_i,x_i_j)
      }
    };rm(i)
    pdci_country_h <- pdci_country_h/10
    gqi_country_h <- gqi_country_h/12
    tqi_country_h <- tqi_country_h/3
  } else {
    pdci_country_h <- NA
    gqi_country_h <- NA
    tqi_country_h <- NA
  }
  ################################################################################
  #joining data per country in one dataframe
  message("Joining PDCI, GQI, and TQI per country in one data.frame")
  completeness_country_final <- as.data.frame(matrix(nrow=length(countries),ncol=10))
  colnames(completeness_country_final) <- c("country",
                                            "PDCI_wild","PDCI_landrace","PDCI_hybrid",
                                            "GQI_wild","GQI_landrace","GQI_hybrid",
                                            "TQI_wild","TQI_landrace","TQI_hybrid"
  )
  completeness_country_final$country <- countries
  
  #adding PDCI and GQI values 
  
  ##WILD
  if(exists("spp_w")){
    #PDCI
    completeness_country_final$PDCI_wild <- unlist(lapply(1:ncol(df_wild),function(i){
      x <- weighted.mean(x = pdci_country_w[,i],
                         w=df_wild[,i]/sum(df_wild[,i],na.rm = T),na.rm = T)
      return(x)
    }))
    
    #GQI
    completeness_country_final$GQI_wild <- unlist(lapply(1:ncol(df_wild),function(i){
      x <- weighted.mean(x = gqi_country_w[,i],
                         w=df_wild[,i]/sum(df_wild[,i],na.rm = T),na.rm = T)
      return(x)
    }))
    
    
    #TQI
    completeness_country_final$TQI_wild <- 
      unlist(lapply(1:ncol(df_wild),function(i){
        x <- weighted.mean(x = tqi_country_w[,i],
                           w=df_wild[,i]/sum(df_wild[,i],na.rm = T),na.rm = T)
        return(x)
      }))
    
  } else {
    completeness_country_final$PDCI_wild <- NA
    completeness_country_final$GQI_wild <- NA
    completeness_country_final$TQI_wild <- NA
  }
  
  
  if(exists("spp_l")){
    completeness_country_final$PDCI_landrace <- unlist(lapply(1:ncol(df_landrace),function(i){
      x <- weighted.mean(x = pdci_country_l[,i],
                         w=df_landrace[,i]/sum(df_landrace[,i],na.rm = T),na.rm = T)
      
      return(x)
    }))
    
    #GQI
    completeness_country_final$GQI_landrace <- 
      unlist(lapply(1:ncol(df_landrace),function(i){
        x <- weighted.mean(x = gqi_country_l[,i],
                           w=df_landrace[,i]/sum(df_landrace[,i],na.rm = T),na.rm = T)
        return(x)
      }))
    
    #TQI
    completeness_country_final$TQI_landrace <- 
      unlist(lapply(1:ncol(df_landrace),function(i){
        x <- weighted.mean(x = tqi_country_l[,i],
                           w=df_landrace[,i]/sum(df_landrace[,i],na.rm = T),na.rm = T)
        return(x)
      }))
  } else {
    completeness_country_final$PDCI_landrace <- NA
    completeness_country_final$GQI_landrace <- NA
    completeness_country_final$TQI_landrace <- NA
  }
  
  if(exists("spp_h")){
    completeness_country_final$PDCI_hybrid <- 
      unlist(lapply(1:ncol(df_hybrid),function(i){
        x <- weighted.mean(x = pdci_country_h[,i],
                           w=df_hybrid[,i]/sum(df_hybrid[,i],na.rm = T),na.rm = T)
        return(x)
      }))
    completeness_country_final$GQI_hybrid <- 
      unlist(lapply(1:ncol(df_hybrid),function(i){
        x <- weighted.mean(x = gqi_country_h[,i],
                           w=df_hybrid[,i]/sum(df_hybrid[,i],na.rm = T),na.rm = T)#
        
        return(x)
      }))
    completeness_country_final$TQI_hybrid <- 
      unlist(lapply(1:ncol(df_hybrid),function(i){
        x <- weighted.mean(x = tqi_country_h[,i],
                           w=df_hybrid[,i]/sum(df_hybrid[,i],na.rm = T),na.rm = T)
        return(x)
      }))
  } else {
    completeness_country_final$PDCI_hybrid <- NA
    completeness_country_final$GQI_hybrid <- NA
    completeness_country_final$TQI_hybrid <- NA
  }
  
  completeness_country_final[length(countries)+1,1] <- "collection"
  
  #adding collection values
  #PDCI
  completeness_country_final$PDCI_wild[[length(countries)+1]] <- 
    df_completeness_final$PDCI[which(df_completeness_final$status=="wild")]
  completeness_country_final$PDCI_landrace[[length(countries)+1]] <- 
    df_completeness_final$PDCI[which(df_completeness_final$status=="landrace")]
  completeness_country_final$PDCI_hybrid[[length(countries)+1]] <- 
    df_completeness_final$PDCI[which(df_completeness_final$status=="hybrid")]
  #GQI
  completeness_country_final$GQI_wild[[length(countries)+1]] <- 
    df_completeness_final$GQI[which(df_completeness_final$status=="wild")]
  completeness_country_final$GQI_landrace[[length(countries)+1]] <- 
    df_completeness_final$GQI[which(df_completeness_final$status=="landrace")]
  completeness_country_final$GQI_hybrid[[length(countries)+1]] <- 
    df_completeness_final$GQI[which(df_completeness_final$status=="hybrid")]
  #TQI
  completeness_country_final$TQI_wild[[length(countries)+1]] <- 
    df_completeness_final$TQI[which(df_completeness_final$status=="wild")]
  completeness_country_final$TQI_landrace[[length(countries)+1]] <- 
    df_completeness_final$TQI[which(df_completeness_final$status=="landrace")]
  completeness_country_final$TQI_hybrid[[length(countries)+1]] <- 
    df_completeness_final$TQI[which(df_completeness_final$status=="hybrid")]
  ################################################################################
  #COMPLETENESS score per country
  message("Calculating completeness index per country")
  completeness_country_final$SCORE_wild <- 
    
    unlist(lapply(1:nrow(completeness_country_final), function(i){
      
      x <- exp(mean(log(c(completeness_country_final$PDCI_wild[[i]],
                           completeness_country_final$GQI_wild[[i]],
                           completeness_country_final$TQI_wild[[i]]))))
      return(x)
    }))

    #completeness_country_final$PDCI_wild * completeness_country_final$GQI_wild *
    #completeness_country_final$TQI_wild
  completeness_country_final$SCORE_landrace <- 
    unlist(lapply(1:nrow(completeness_country_final), function(i){
      
      x <- exp(mean(log(c(completeness_country_final$PDCI_landrace[[i]],
                          completeness_country_final$GQI_landrace[[i]],
                          completeness_country_final$TQI_landrace[[i]]))))
      return(x)
    }))
    

   # completeness_country_final$PDCI_landrace * completeness_country_final$GQI_landrace *
   # completeness_country_final$TQI_landrace
  completeness_country_final$SCORE_hybrid <- 
    
    unlist(lapply(1:nrow(completeness_country_final), function(i){
      
      x <- exp(mean(log(c(completeness_country_final$PDCI_hybrid[[i]],
                          completeness_country_final$GQI_hybrid[[i]],
                          completeness_country_final$TQI_hybrid[[i]]))))
      return(x)
    }))

  ################################################################################
  #creating regions for plots
  completeness_country_final$REGION <- NA
  completeness_country_final$REGION <- countrycode(completeness_country_final$country,
                                                   origin = 'iso3c',
                                                   destination = 'region')
  
  completeness_country_final$country[which(completeness_country_final$country== "collection")] <- "Subset"
  completeness_country_final$country[which(completeness_country_final$country=="NA1")] <- "No_country"
  completeness_country_final$REGION[which(is.na(completeness_country_final$REGION))] <- "N/A"
  completeness_country_final$REGION[which(completeness_country_final$country== "Subset")] <- "Subset"
  completeness_country_final$REGION[which(completeness_country_final$country == "No_country")] <- "Subset"
  completeness_country_final$REGION[which(completeness_country_final$country == "SCG")] <- "Europe & Central Asia"
  completeness_country_final$REGION[which(completeness_country_final$country == "YUG")] <- "Europe & Central Asia"
  completeness_country_final$REGION[which(completeness_country_final$country == "ZAR")] <- "Sub-Saharan Africa"
  completeness_country_final$REGION[which(completeness_country_final$country == "BUR")] <- "East Asia & Pacific"
  ###
  row.names(completeness_country_final) <- completeness_country_final$country
  ###
  
  ################################################################################
  #subsetting for NIPALS
  message("subsetting for NIPALS graphics")
  ################################################################################
  x_res_W <- #simp_list_nona[,c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27)]
    completeness_country_final[, c(#2, #3,
      "PDCI_wild",
      "GQI_wild",
      "TQI_wild",
      "SCORE_wild",
      "REGION")]
  
  colnames(x_res_W) <- c(
    "Passport data completeness index (PDCI)",
    "Geographical Quality Score index  (GQS)",
    "Taxonomic completeness index  (TQI)",
    "Documentation Completeness index (DCI)",
    "REGION"
  )
  
  #complete cases
  x_res_W <- x_res_W[complete.cases(x_res_W), ]
  x_res_W$REGION <- factor(x_res_W$REGION)
  ##############################################################################
  #subsetting for NIPALS
  x_res_L <- #simp_list_nona[,c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27)]
    completeness_country_final[, c(#2, #3,
      "PDCI_landrace",
      "GQI_landrace",
      "TQI_landrace",
      "SCORE_landrace",
      "REGION")]
  
  colnames(x_res_L) <- c(
    "Passport data completeness index (PDCI)",
    "Geographical Quality Score index  (GQS)",
    "Taxonomic completeness index  (TQI)",
    "Documentation Completeness index (DCI)",
    "REGION"
  )
  
  #complete cases
  x_res_L <- x_res_L[complete.cases(x_res_L), ]
  x_res_L$REGION <- factor(x_res_L$REGION)
  ################################################################################
  set.seed(1000)
  ################################################################################
  #plots
  message("Plotting crop wild relatives NIPALS")
  
  if(nrow(x_res_W)>=2){
    x_W_A = nipals(x_res_W[, -ncol(x_res_W)], comps = 3,scaled = T)
    
    x_res_W_NA <- x_res_W[!row.names(x_res_W) %in% c("No_country","Subset"),]
    # Create a scatter plot
    #i <- ggplot() + geom_point(aes(x_W_A$cor.xt[,1], x_W_A$c[,2]))
    ind_coords <- data.frame(Dim.1=x_W_A$scores[,1],Dim.2=x_W_A$scores[,2],
                             REGION=x_res_W$REGION,countries=row.names(x_W_A$score))
    #ind_coords$``
    var_coords <- as.data.frame(x_W_A$cor.xt)
    var_coords$vars <-  row.names(var_coords)
    # Add arrow
    
    r <- min((max(ind_coords[, "Dim.1"]) - min(ind_coords[, "Dim.1"])/(max(var_coords[, "t1"]) - 
                                                                         min(var_coords[, "t1"]))), (max(ind_coords[, "Dim.2"]) - min(ind_coords[, "Dim.2"])/(max(var_coords[, 
                                                                                                                                                                             "t2"]) - min(var_coords[, "t2"]))))
    
    #var_coords <- var_coords*(r * 0.7)
    PCA_WILD <- ggplot()  +
      geom_hline(yintercept = 0,linetype = "dashed")+
      geom_vline(xintercept = 0,linetype = "dashed")+
      xlab(paste0("Dim1 (",round(x_W_A$values$percentage[[1]],1),"%)"))+
      ylab(paste0("Dim2 (",round(x_W_A$values$percentage[[2]],1),"%)"))+
      ggtitle(paste0("Crop Wild Relatives documentation completeness (", nrow(x_res_W_NA), " countries)"))+
      
      geom_segment(data = var_coords, aes(x = 0, y = 0, xend = t1*(r * 0.7), yend = t2*(r * 0.7)),colour = "steelblue",
                   arrow = arrow(length = unit(0.2, "cm")))+
      
      geom_text_repel(
        data          = var_coords,
        mapping       = aes(t1*(r * 0.7), t2*(r * 0.7), label = vars),
        size = 3,    # Repel away from the left edge, not from the right.
        color = "steelblue", 
        #xlim = c(NA, Inf),
        # Do not repel from top or bottom edges.
        #ylim = c(-Inf, Inf),
        max.overlaps = 10000,#200,  
        max.iter = 10000,show.legend = NA,box.padding = 0.05,
        point.padding = 0.02, 
        nudge_x = 4.5,
        nudge_y = .1,
        segment.linetype = 6,
        segment.curvature = -1e-20,hjust         = 0.5
      ) +
      geom_point(data = ind_coords,aes(x = Dim.1,y = Dim.2,colour= REGION),size = 5)+
      scale_color_manual(values = c("Middle East & North Africa" = "#88a0dc",
                                    "Europe & Central Asia" = "#381a61",
                                    "Latin America & Caribbean" = "#7c4b73",
                                    "South Asia"  = "#ed968c",
                                    "East Asia & Pacific" = "#ab3329",
                                    "Sub-Saharan Africa"  = "#e78429",
                                    "North America"  = "#f9d14a",
                                    "Subset" = "grey50",
                                    "N/A"  = "black"
      )
      ) +
      geom_text_repel(
        data          = ind_coords,
        mapping       = aes(Dim.1, Dim.2, label = countries),
        size = 2.7,    # Repel away from the left edge, not from the right.
        #xlim = c(NA, Inf),
        # Do not repel from top or bottom edges.
        #ylim = c(-Inf, Inf),
        max.overlaps = 1000,#200,
        max.iter = 10000,show.legend = NA,box.padding = 0.05,
        point.padding = 0.02, 
        nudge_x = .1,
        nudge_y = .1,
        segment.linetype = 6,
        segment.curvature = -1e-20,
        hjust= 1
      ) +
      
      theme_minimal()
    
    #PCA_WILD
    
    ggplot2::ggsave(
      paste0(
        outdir,
        "/",
        collection_name,
        "/",
        collection_name,
        "_completeness_WILD.png"
      ),
      PCA_WILD,
      dpi = 300,
      units = "in",
      width = 10,
      height = 8
    )
  } else {
    warning("NO WILD AVAILABLE")
  }
  ###############################################################################
  message("Plotting landraces NIPALS")
  
  if(nrow(x_res_L)>=2){
    x_W_A = nipals(x_res_L[, -ncol(x_res_L)], comps = 3,scaled = T)
    
    x_res_L_NA <- x_res_L[!row.names(x_res_L) %in% c("No_country","Subset"),]
    # Create a scatter plot
    #i <- ggplot() + geom_point(aes(x_W_A$cor.xt[,1], x_W_A$c[,2]))
    ind_coords <- data.frame(Dim.1=x_W_A$scores[,1],Dim.2=x_W_A$scores[,2],
                             REGION=x_res_L$REGION,countries=row.names(x_W_A$score))
    #ind_coords$``
    var_coords <- as.data.frame(x_W_A$cor.xt)
    var_coords$vars <-  row.names(var_coords)
    # Add arrow
    
    r <- min((max(ind_coords[, "Dim.1"]) - min(ind_coords[, "Dim.1"])/(max(var_coords[, "t1"]) - 
                                                                         min(var_coords[, "t1"]))), (max(ind_coords[, "Dim.2"]) - min(ind_coords[, "Dim.2"])/(max(var_coords[, 
                                                                                                                                                                             "t2"]) - min(var_coords[, "t2"]))))
    
    #var_coords <- var_coords*(r * 0.7)
    PCA_LAND <- ggplot()  +
      geom_hline(yintercept = 0,linetype = "dashed")+
      geom_vline(xintercept = 0,linetype = "dashed")+
      xlab(paste0("Dim1 (",round(x_W_A$values$percentage[[1]],1),"%)"))+
      ylab(paste0("Dim2 (",round(x_W_A$values$percentage[[2]],1),"%)"))+
      ggtitle(paste0("Landraces documentation completeness (", nrow(x_res_L_NA), " countries)"))+
      
      geom_segment(data = var_coords, aes(x = 0, y = 0, xend = t1*(r * 0.7), yend = t2*(r * 0.7)),colour = "steelblue",
                   arrow = arrow(length = unit(0.2, "cm")))+
      
      geom_text_repel(
        data          = var_coords,
        mapping       = aes(t1*(r * 0.7), t2*(r * 0.7), label = vars),
        size = 3,    # Repel away from the left edge, not from the right.
        color = "steelblue", 
        #xlim = c(NA, Inf),
        # Do not repel from top or bottom edges.
        #ylim = c(-Inf, Inf),
        max.overlaps = 10000,#200,  
        max.iter = 10000,show.legend = NA,box.padding = 0.05,
        point.padding = 0.02, 
        nudge_x = 4.5,
        nudge_y = .1,
        segment.linetype = 6,
        segment.curvature = -1e-20,hjust         = 0.5
      ) +
      geom_point(data = ind_coords,aes(x = Dim.1,y = Dim.2,colour= REGION),size = 5)+
      scale_color_manual(values = c("Middle East & North Africa" = "#88a0dc",
                                    "Europe & Central Asia" = "#381a61",
                                    "Latin America & Caribbean" = "#7c4b73",
                                    "South Asia"  = "#ed968c",
                                    "East Asia & Pacific" = "#ab3329",
                                    "Sub-Saharan Africa"  = "#e78429",
                                    "North America"  = "#f9d14a",
                                    "Subset" = "grey50",
                                    "N/A"  = "black"
      )
      ) +
      geom_text_repel(
        data          = ind_coords,
        mapping       = aes(Dim.1, Dim.2, label = countries),
        size = 2.7,    # Repel away from the left edge, not from the right.
        #xlim = c(NA, Inf),
        # Do not repel from top or bottom edges.
        #ylim = c(-Inf, Inf),
        max.overlaps = 1000,#200,
        max.iter = 10000,show.legend = NA,box.padding = 0.05,
        point.padding = 0.02, 
        nudge_x = .1,
        nudge_y = .1,
        segment.linetype = 6,
        segment.curvature = -1e-20,
        hjust= 1
      ) +
      
      theme_minimal()
    
    #PCA_WILD
    
    ggplot2::ggsave(
      paste0(
        outdir,
        "/",
        collection_name,
        "/",
        collection_name,
        
        "_completeness_LANDRACE.png"
      ),
      PCA_LAND,
      dpi = 300,
      units = "in",
      width = 10,
      height = 8
    )
  } else {
    warning("NO WILD AVAILABLE")
  }
  ###############################################################################
  #saving all
  message("Saving file")
  write.csv(
    completeness_country_final,
    paste0(
      outdir,
      "/",
      collection_name,
      "/",
      collection_name,
      "_completeness_countries.csv"
    ),
    row.names = F,
    na = ""
  )
  
  message(
    paste(collection_name, "processed for countries completeness!")
  )
  
  return(completeness_country_final)
}
################################################################################
dir <- "D:/OneDrive - CGIAR/GERMPLASM_INDEX"
#dir <- "D:/ONEDRIVE/cgiar/OneDrive - CGIAR/GERMPLASM_INDEX"
################################################################################
#outdir
outdir <- paste0(dir, "/BEANS/RESULTS")
################################################################################
#reading dataframes with quality scores
#PDCI
PCDI_df <- read.csv("D:/OneDrive - CGIAR/GERMPLASM_INDEX/CIAT Data/ciat_pdci.csv")
#GEOGRAPHICAL QUALITY SCORE
QGI_df <- read.csv("D:/OneDrive - CGIAR/GERMPLASM_INDEX/CIAT Data/genesys_quality_score_08-2024.csv")
#TAXONOMIC QUALITY SCORE INDEX
TQI_df <- read.csv("D:/OneDrive - CGIAR/GERMPLASM_INDEX/TAXONOMIC QUALITY INDEX/ciat_v6.1.csv")


################################################################################
collection_name <- "beans"
PCDI_df_beans <- PCDI_df[which(PCDI_df$CROPNAME==collection_name),]
QGI_df_beans <- QGI_df[which(QGI_df$CROPNAME==collection_name),]
TQI_df_beans <- TQI_df[which(TQI_df$CROPNAME==collection_name),]
x1 <-completeness_func(outdir,collection_name,PCDI_df_beans,QGI_df_beans,TQI_df_beans)
################################################################################
collection_name <- "Cassava"
PCDI_df_cassava <- PCDI_df[which(PCDI_df$CROPNAME==collection_name),]
QGI_df_cassava <- QGI_df[which(QGI_df$CROPNAME==collection_name),]
TQI_df_cassava <- TQI_df[which(TQI_df$CROPNAME==collection_name),]

x2 <-completeness_func(outdir,collection_name,PCDI_df_cassava,QGI_df_cassava,TQI_df_cassava)
################################################################################
collection_name <- "forages"
PCDI_df_forages <- PCDI_df[which(PCDI_df$CROPNAME==collection_name),]
QGI_df_forages <- QGI_df[which(QGI_df$CROPNAME==collection_name),]
TQI_df_forages <- TQI_df[which(TQI_df$CROPNAME==collection_name),]
x3 <-completeness_func(outdir,collection_name,PCDI_df_forages,QGI_df_forages,TQI_df_forages)
