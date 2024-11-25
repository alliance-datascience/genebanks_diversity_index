require(data.table)
require(readxl)
library(ggplot2)
library(factoextra)
library(glue)
library(fossil)
library(plsdepot)
library(ggrepel)
require(parallel)
library(plsdepot)
library(countrycode)


genetics_ind_function <- function(outdir,
                                  collection_name,
                                  numCores,
                                  accessions_df) {
  
  message(paste0("Processing collection: ",collection_name,"... "))
  message("Loading results of 1_Taxonomic_index")
  ################################################################################
  subsets <- readRDS(paste0(outdir, "/", collection_name, "/", collection_name, "_subsets_new_1.RDS"))
  ################################################################################
  #loading summary table to get IUCN
  dist_prop <- read.csv(
    paste0(outdir, "/", collection_name,"/", collection_name,  "_summary_table.csv"))
  ################################################################################
  #polish accession numbers by removing extra spaces
  accessions_df[, 1] <- trimws(accessions_df[, 1])
  #Obtaining unique accessions
  un <- unique(accessions_df[, 1])
  ################################################################################
  #loading files required to obtain index
  passport_data_w <- subsets[[1]]
  passport_data_l <- subsets[[2]]
  passport_data_h <- subsets[[3]]
  df_wild <- subsets[[10]]
  df_landrace <- subsets[[11]]
  df_hybrid <- subsets[[12]]
  
  if (!is.data.frame(df_wild)) {
    rm(df_wild)
  }
  if (!is.data.frame(df_landrace)) {
    rm(df_landrace)
  }
  if (!is.data.frame(df_hybrid)) {
    rm(df_hybrid)
  }
  ################################################################################
  ################################################################################
  message("Creating preprocessing summary tables")
  ################################################################################
  #summary file for wild
  if (exists("df_wild")) {
    stat <- strsplit(row.names(df_wild), "_")
    wild_seq <- data.frame(
      taxon = unlist(lapply(stat, `[[`, 1)),
      status = "wild",
      seq_records = NA,
      records = rowSums(df_wild),
      prop_seq = NA
    )
  }
  ################################################################################
  #summary file for landrace
  if (exists("df_landrace")) {
    stat <- strsplit(row.names(df_landrace), "_")
    
    wild_lan <- data.frame(
      taxon = unlist(lapply(stat, `[[`, 1)),
      status = "landrace",
      seq_records = NA,
      records = rowSums(df_landrace),
      prop_seq = NA
    )
  }
  ################################################################################
  if (exists("df_hybrid")) {
    stat <- strsplit(row.names(df_hybrid), "_")
    #summary file for hybrids
    wild_hyb <- data.frame(
      taxon = unlist(lapply(stat, `[[`, 1)),
      status = "hybrid",
      seq_records = NA,
      records = rowSums(df_hybrid),
      prop_seq = NA
    )
  }
  ################################################################################
  message("Checking cores to use parallel approach")
  x_det <- parallel::detectCores()
  if (numCores > x_det) {
    stop(
      "Number of cores exceed the maximum allowed by the machine, use a coherent number of cores such as four"
    )
  }
  
  ################################################################################
  message("Matching sequenced accessions using ACCENUMB")
  ################################################################################
  #sequenced for wild
  ################################################################################
  #Assigning sequenced flag for wild records using parallel library
  
  message("Sequenced accessions for wild records")
  
  cl <- parallel::makeCluster(numCores)
  parallel::clusterExport(cl,
                          varlist = c("passport_data_w", "un"),
                          envir = environment())
  
  res <- parallel::parLapplyLB(
    cl,
    X = seq_len(nrow(passport_data_w)),
    fun = function (i) {
      if (isTRUE(passport_data_w$ACCENUMB[[i]] %in% un)) {
        x <- 1
      } else {
        x <- 0
      }
    }
  )
  
  parallel::stopCluster(cl)
  rm(cl)
  
  #assigning flags
  passport_data_w$SEQ_FLAG <- unlist(res)
  rm(res)
  #sum(passport_data_w$SEQ_FLAG)
  ################################################################################
  #Obtaining summaries for sequenced accessions
  if (nrow(passport_data_w) > 0) {
    x <- tapply(passport_data_w$SEQ_FLAG[which(passport_data_w$SEQ_FLAG == 1)], 
                passport_data_w$NAME[which(passport_data_w$SEQ_FLAG ==1)], length)
    x_total <- tapply(passport_data_w$ACCENUMB, passport_data_w$NAME, length)
    #records per country
    x_country_w <- list()
    IUCN_country_w <- list()
    for(i in 1:ncol(df_wild)){
      #i <- 1
      x_i_c <- passport_data_w[which(passport_data_w$SEQ_FLAG == 1 &
                                       passport_data_w$ORIGCTY==colnames(df_wild)[i]),]
      x2 <- tapply(x_i_c$SEQ_FLAG[which(x_i_c$SEQ_FLAG == 1)], 
                   x_i_c$NAME[which(x_i_c$SEQ_FLAG ==1)], length)
      if(length(x2)>0){
        for(j in 1:nrow(x2)){
          x2[[j]] <- x2[[j]]/x_total[names(x_total) %in% names(x2[j])]
        }
      } else {
        x2 <- NA
      }
      x_country_w[[i]] <- data.frame(COUNTRY=colnames(df_wild)[i],
                                     pRseq= mean(x2,na.rm = T),
                                     pTseq=sum(x2>0,na.rm = T)/nrow(dist_prop),#exp(mean(log(x2),na.rm = T)),
                                     SUBSET="wild")
      #get IUCN proportion per country
      x_dist_W <- dist_prop[dist_prop$IUCN %in%
                              c("Vulnerable",
                                "Endangered",
                                "Critically Endangered",
                                "Extinct in the Wild") &
                              dist_prop$status == "wild", ]
      #obtaining IUCN species
      WILD_SP_DANGERED <- paste0(x_dist_W$Taxa,"_",x_dist_W$status)
      df_wild_DANGERED <- df_wild[row.names(df_wild) %in% WILD_SP_DANGERED,]
      df_wild_DANGERED[df_wild_DANGERED>0] <- 1
      #df_wild_DANGERED[df_wild_DANGERED==0] <- NA
      x_w_sp <- colSums(df_wild_DANGERED,na.rm = T)/nrow(dist_prop)
      sp_dangered_c_w <- data.frame(COUNTRY=names(x_w_sp),IUCN=x_w_sp)
      #sp_dangered_c_w$IUCN[which(sp_dangered_c_w$IUCN==0)] <- NA
    }
    x_country_w <- do.call(rbind,x_country_w)  
    
    
  } else {
    x <- NULL
  }
  
  
  #creating a df for wild species and sequenced records
  if (length(x) > 0) {
    seq_w <- data.frame(taxa = paste0(names(x), ""), freq = x)
    rm(x)
    
    for (i in 1:nrow(wild_seq)) {
      x <- seq_w$freq[seq_w$taxa %in% wild_seq$taxon[[i]]]
      if (length(x) > 0) {
        wild_seq$seq_records[[i]]  <- x
      } else {
        wild_seq$seq_records[[i]]  <- NA
      }
    }
    rm(i)
    rm(x)
    
    #proportion of sequenced accessions of total
    wild_seq$prop_seq <- wild_seq$seq_records / wild_seq$records
  }
  ################################################################################
  #sequenced for landraces
  ################################################################################
  
  message("Sequenced accessions for landrace records")
  
  cl <- parallel::makeCluster(numCores)
  parallel::clusterExport(cl,
                          varlist = c("passport_data_l", "un"),
                          envir = environment())
  
  res <- parallel::parLapplyLB(
    cl,
    X = seq_len(nrow(passport_data_l)),
    fun = function (i) {
      if (isTRUE(passport_data_l$ACCENUMB[[i]] %in% un)) {
        x <- 1
      } else {
        x <- 0
      }
    }
  )
  
  parallel::stopCluster(cl)
  rm(cl)
  passport_data_l$SEQ_FLAG <- unlist(res)
  rm(res)
  
  #Obtaining summaries for sequenced accessions
  if (nrow(passport_data_l) > 0) {
    x <- tapply(passport_data_l$SEQ_FLAG[which(passport_data_l$SEQ_FLAG == 1)], 
                passport_data_l$NAME[which(passport_data_l$SEQ_FLAG == 1)], length)
    x_total <- tapply(passport_data_l$ACCENUMB, passport_data_l$NAME, length)
    
    #records per country
    x_country_l <- list()
    for(i in 1:ncol(df_landrace)){
      
      x_i_c <- passport_data_l[which(passport_data_l$SEQ_FLAG == 1 &
                                       passport_data_l$ORIGCTY==colnames(df_landrace)[i]),]
      x2 <- tapply(x_i_c$SEQ_FLAG[which(x_i_c$SEQ_FLAG == 1)], 
                   x_i_c$NAME[which(x_i_c$SEQ_FLAG ==1)], length)
      if(length(x2)>0){
        for(j in 1:nrow(x2)){
          x2[[j]] <- x2[[j]]/x_total[names(x_total) %in% names(x2[j])]
          
        }
      } else {
        x2 <- NA
      }
      x_country_l[[i]] <- data.frame(COUNTRY=colnames(df_landrace)[i],
                                     pRseq=mean(x2,na.rm = T),#exp(mean(log(x2),na.rm = T)),
                                     pTseq=sum(x2>0,na.rm = T)/nrow(dist_prop),
                                     SUBSET="landraces")
      #get IUCN proportion per country
      x_dist_l <- dist_prop[dist_prop$IUCN %in%
                              c("Vulnerable",
                                "Endangered",
                                "Critically Endangered",
                                "Extinct in the Wild") &
                              dist_prop$status == "landrace", ]
      
      LAND_SP_DANGERED <- paste0(x_dist_l$Taxa,"_",x_dist_l$status)
      df_landrace_DANGERED <- df_landrace[row.names(df_landrace) %in% LAND_SP_DANGERED,]
      df_landrace_DANGERED[df_landrace_DANGERED>0] <- 1
      #df_landrace_DANGERED[df_landrace_DANGERED==0] <- NA
      
      x_l_sp <- colSums(df_landrace_DANGERED,na.rm = T)/nrow(dist_prop)
      sp_dangered_c_l <- data.frame(COUNTRY=names(x_l_sp),IUCN=x_l_sp)
      #sp_dangered_c_l$IUCN[which(sp_dangered_c_l$IUCN==0)] <- NA
      
    }
    x_country_l <- do.call(rbind,x_country_l)   
  } else {
    x <- NULL
  }
  
  #creating a df for landraces species and sequenced records
  if (length(x) > 0) {
    seq_l <- data.frame(taxa = paste0(names(x), ""), freq = x)
    rm(x)
    
    stat <- strsplit(row.names(df_landrace), "_")
    
    
    
    for (i in 1:nrow(wild_lan)) {
      x <- seq_l$freq[seq_l$taxa %in% wild_lan$taxon[[i]]]
      if (length(x) > 0) {
        wild_lan$seq_records[[i]]  <- x
      } else {
        wild_lan$seq_records[[i]]  <- NA
      }
    }
    rm(i)
    rm(x)
    #proportion of sequenced accessions of total
    wild_lan$prop_seq <- wild_lan$seq_records / wild_lan$records
  }
  ################################################################################
  #sequenced for hybrids
  ################################################################################
  
  message("Sequenced accessions for hybrid records")
  
  cl <- parallel::makeCluster(numCores)
  parallel::clusterExport(cl,
                          varlist = c("passport_data_h", "un"),
                          envir = environment())
  
  res <- parallel::parLapplyLB(
    cl,
    X = seq_len(nrow(passport_data_h)),
    fun = function (i) {
      if (isTRUE(passport_data_h$ACCENUMB[[i]] %in% un)) {
        x <- 1
      } else {
        x <- 0
      }
    }
  )
  
  parallel::stopCluster(cl)
  rm(cl)
  passport_data_h$SEQ_FLAG <- unlist(res)
  rm(res)
  
  #Obtaining summaries for sequenced accessions
  if (nrow(passport_data_h) > 0) {
    x <- tapply(passport_data_h$SEQ_FLAG[which(passport_data_h$SEQ_FLAG == 1)], 
                passport_data_h$NAME[which(passport_data_h$SEQ_FLAG ==1)], length)
    x_total <- tapply(passport_data_h$ACCENUMB, passport_data_h$NAME, length)
    
    #records per country
    x_country_h <- list()
    
    for(i in 1:ncol(df_hybrid)){
      
      x_i_c <- passport_data_h[which(passport_data_h$SEQ_FLAG == 1 &
                                       passport_data_h$ORIGCTY==colnames(df_hybrid)[i]),]
      x2 <- tapply(x_i_c$SEQ_FLAG[which(x_i_c$SEQ_FLAG == 1)], 
                   x_i_c$NAME[which(x_i_c$SEQ_FLAG ==1)], length)
      if(length(x2)>0){
        for(j in 1:nrow(x2)){
          x2[[j]] <- x2[[j]]/x_total[names(x_total) %in% names(x2[j])]
        }
      } else {
        x2 <- NA
      }
      
      x_country_h[[i]] <- data.frame(COUNTRY=colnames(df_hybrid)[i],
                                     pRseq=mean(x2,na.rm = T),#exp(mean(log(x2),na.rm = T)),
                                     pTseq=sum(x2>0,na.rm = T)/nrow(dist_prop),
                                     SUBSET="hybrid")
      
      x_dist_h <- dist_prop[dist_prop$IUCN %in%
                              c("Vulnerable",
                                "Endangered",
                                "Critically Endangered",
                                "Extinct in the Wild") &
                              dist_prop$status == "hybrid", ]
      #obtaining IUCN proportion
      HYBD_SP_DANGERED <- paste0(x_dist_h$Taxa,"_",x_dist_h$status)
      df_hybd_DANGERED <- x_dist_h[row.names(df_hybrid) %in% HYBD_SP_DANGERED,]
      df_hybd_DANGERED[df_hybd_DANGERED>0] <- 1
      #df_hybd_DANGERED[df_hybd_DANGERED==0] <- NA
      
      x_h_sp <- colSums(df_hybd_DANGERED,na.rm = T)/nrow(dist_prop)
      sp_dangered_c_h <- data.frame(COUNTRY=names(x_h_sp),IUCN=x_h_sp)
      #sp_dangered_c_h$IUCN[which(sp_dangered_c_h$IUCN==0)] <- NA
      
    }
    x_country_h <- do.call(rbind,x_country_h)
  } else {
    x <- NULL
  }
  #creating a df for hybrid species and sequenced records
  if (length(x) > 0) {
    seq_h <- data.frame(taxa = paste0(names(x), ""), freq = x)
    rm(x)
    
    stat <- strsplit(row.names(df_hybrid), "_")
    
    
    for (i in 1:nrow(wild_hyb)) {
      x <- seq_h$freq[seq_h$taxa %in% wild_hyb$taxon[[i]]]
      if (length(x) > 0) {
        wild_hyb$seq_records[[i]]  <- x
      } else {
        wild_hyb$seq_records[[i]]  <- NA
      }
    }
    rm(i)
    
    #proportion of sequenced accessions of total
    wild_hyb$prop_seq <- wild_hyb$seq_records / wild_hyb$records
  }
  
  ##############################################################################
  #creating geometric mean of the proportion of records per taxa and country sequenced
  if(exists("df_wild")){
    x_country <- data.frame(country=colnames(df_wild),
                            pREC_wild=NA,pREC_landrace=NA,pREC_hybrid=NA,
                            pSeq_wild=NA,pSeq_landrace=NA,pSeq_hybrid=NA,
                            IUCN_wild=NA,IUCN_landrace=NA,IUCN_hybrid=NA)
  } else if(exists("df_landrace")){
    x_country <- data.frame(country=colnames(df_landrace),
                            pREC_wild=NA,pREC_landrace=NA,pREC_hybrid=NA,
                            pSeq_wild=NA,pSeq_landrace=NA,pSeq_hybrid=NA,
                            IUCN_wild=NA,IUCN_landrace=NA,IUCN_hybrid=NA)
  }
  
  for(i in 1:nrow(x_country)){
    if(exists("x_country_w")){
      x <- x_country_w$pRseq[which(
        x_country_w$COUNTRY==x_country$country[[i]])]
      x_T <- x_country_w$pTseq[which(
        x_country_w$COUNTRY==x_country$country[[i]])]
      x_IUCN <- sp_dangered_c_w$IUCN[which(
        sp_dangered_c_w$COUNTRY==x_country$country[[i]])]
      
      if(length(x)>0){
        x_country$pREC_wild[[i]] <- x
      }
      if(length(x_T)>0){
        x_country$pSeq_wild[[i]] <- x_T
      }      
      if(length(x_IUCN)>0){
        x_country$IUCN_wild[[i]] <- x_IUCN
      }
    }
    
    if(exists("x_country_l")){
      x <- x_country_l$pRseq[which(
        x_country_l$COUNTRY==x_country$country[[i]])]
      x_T <- x_country_l$pTseq[which(
        x_country_l$COUNTRY==x_country$country[[i]])]
      x_IUCN <- sp_dangered_c_l$IUCN[which(
        sp_dangered_c_l$COUNTRY==x_country$country[[i]])]
      if(length(x)>0){
        x_country$pREC_landrace[[i]] <- x
      }
      if(length(x_T)>0){
        x_country$pSeq_landrace[[i]] <- x_T
      }   
      if(length(x_IUCN)>0){
        x_country$IUCN_landrace[[i]] <- x_IUCN
      }
    }
    
    if(exists("x_country_h")){
      x <- x_country_h$pRseq[which(
        x_country_h$COUNTRY==x_country$country[[i]])]
      x_T <- x_country_h$pTseq[which(
        x_country_h$COUNTRY==x_country$country[[i]])]      
      x_IUCN <- sp_dangered_c_h$IUCN[which(
        sp_dangered_c_h$COUNTRY==x_country$country[[i]])]
      if(length(x)>0){
        x_country$pREC_hybrid[[i]] <- x
      }
      if(length(x_T)>0){
        x_country$pSeq_hybrid[[i]] <- x_T
      }   
      if(length(x_IUCN)>0){
        x_country$IUCN_hybrid[[i]] <- x_IUCN
      }
    }
  }
  
  #calculating usability index per country (creating empty variables!)
  x_country$usability_index_wild <- NA
  x_country$usability_index_landrace <- NA
  x_country$usability_index_collection <- NA
  #calculating the index
  for(i in 1:nrow(x_country)){
    x_country$usability_index_wild[[i]] <- mean(c(x_country$pREC_wild[[i]],
                                                  x_country$pSeq_wild[[i]],
                                                  x_country$IUCN_wild[[i]]),na.rm = T)
    x_country$usability_index_landrace[[i]] <- mean(c(x_country$pREC_landrace[[i]],
                                                      x_country$pSeq_landrace[[i]],
                                                      x_country$IUCN_landrace[[i]]),na.rm = T)
    x_country$usability_index_collection[[i]] <- mean(c(x_country$usability_index_wild[[i]],
                                                        x_country$usability_index_landrace[[i]]),na.rm = T)
  }
  ###############################################################################
  #loading IUCN data for collection
  TAX_DF <-
    read.csv(
      paste0(outdir, "/", collection_name, "/", collection_name, "_IT_table_1.csv")
    )
  ################################################################################
  message("Creating final summary table")
  if (exists("wild_seq") & exists("wild_lan")) {
    status_total <- rbind(wild_seq, wild_lan)
  } else if (exists("wild_seq") & !exists("wild_lan")) {
    status_total <- (wild_seq)
  } else if (!exists("wild_seq") & exists("wild_lan")) {
    status_total <- (wild_lan)
  }
  
  if (exists("wild_hyb")) {
    status_total <- rbind(status_total, wild_hyb)
  } else {
    status_total <- status_total
  }
  
  status_total$prop_rec_collection <- NA
  status_total$prop_rec_collection <- status_total$records / sum(status_total$records, na.rm = T)
  status_total$prop_seq_collection <- NA
  status_total$prop_seq_collection <- status_total$seq_records / sum(status_total$seq_records, na.rm = T)
  status_total$sp_status <- NA
  status_total$sp_status[which(status_total$seq_records > 0)] <- 1
  
  status_total$collection <- NA
  status_total$collection <- collection_name
  #sum(status_total$prop_seq,na.rm = T)
  # IND <- sum(status_total$prop_seq * status_total$prop_seq_collection,
  #            na.rm = T) / sum(status_total$prop_seq_collection, na.rm = T)
  # IND_a <- sum(status_total$prop_seq * status_total$prop_rec_collection,
  #            na.rm = T) / sum(status_total$prop_rec_collection, na.rm = T)
  
  #############################################################################
  #SEQUENCED ACCESION PROPORTION IN A SUBSET!
  
  # IND <- weighted.mean(status_total$prop_seq, 
  #                      status_total$prop_rec_collection,
  #                      na.rm=T)
  
  IND_W <- weighted.mean(status_total$prop_seq[which(status_total$status=="wild")], 
                         status_total$prop_rec_collection[which(status_total$status=="wild")],
                         na.rm=T)
  
  IND_L <- weighted.mean(status_total$prop_seq[which(status_total$status=="landrace")], 
                         status_total$prop_rec_collection[which(status_total$status=="landrace")],
                         na.rm=T)
  
  IND_H <- weighted.mean(status_total$prop_seq[which(status_total$status=="hybrid")], 
                         status_total$prop_rec_collection[which(status_total$status=="hybrid")],
                         na.rm=T)
  IND <- (IND_W+IND_L)/2#
  # IND_a <- weighted.mean(status_total$prop_seq, 
  #                      status_total$prop_rec_collection,
  #                      na.rm=T)
  #############################################################################
  #TAXA PROPORTION IN A SUBSET WITH AT LEAST ONE RECORD SEQUENCED!
  IND2_W <- sum(status_total$sp_status[which(status_total$status=="wild")],na.rm = T)/
    nrow(status_total[which(status_total$status=="wild"),])
  IND2_L <- sum(status_total$sp_status[which(status_total$status=="landrace")],na.rm = T)/
    nrow(status_total[which(status_total$status=="landrace"),])
  IND2_H <- sum(status_total$sp_status[which(status_total$status=="hybrid")],na.rm = T)/
    nrow(status_total[which(status_total$status=="hybrid"),])
  
  IND2 <- (IND2_W+IND2_L)/2#sum(status_total$sp_status,na.rm = T)/nrow(status_total)
  
  message("Saving results")
  
  x_div_df <- data.frame(matrix(ncol=6,nrow = 4))
  colnames(x_div_df) <- c("status","pRseq", "pTseq","IUCN","usability_index","collection")
  x_div_df[,1] <- c("wild","landrace","hybrid","total")
  #
  x_div_df[1,2] <- IND_W
  x_div_df[2,2] <- IND_L
  x_div_df[3,2] <- IND_H
  x_div_df[4,2] <- IND
  #
  x_div_df[1,3] <- IND2_W
  x_div_df[2,3] <- IND2_L
  x_div_df[3,3] <- IND2_H
  x_div_df[4,3] <- IND2
  #
  #
  x_div_df[1,4] <- TAX_DF$propTaxa_IUCN[1]
  x_div_df[2,4] <- TAX_DF$propTaxa_IUCN[2]
  x_div_df[3,4] <- TAX_DF$propTaxa_IUCN[3]
  x_div_df[4,4] <- sum(TAX_DF$propTaxa_IUCN[1:2],na.rm = T)/2
  #
  
  x_div_df[1,5] <- mean(as.numeric(x_div_df[1,c(2,3,4)]),na.rm = T)
  x_div_df[2,5] <- mean(as.numeric(x_div_df[2,c(2,3,4)]),na.rm = T)
  x_div_df[3,5] <- mean(as.numeric(x_div_df[3,c(2,3,4)]),na.rm = T)
  x_div_df[4,5] <- sum(x_div_df$usability_index[1:2],na.rm = T)/2
  
  x_div_df[,6] <- collection_name
  
  results <- list(summary_table = status_total, 
                  indexes= x_div_df,
                  countries = x_country)
  write.csv(
    status_total,
    paste0(outdir, "/", collection_name,"/", collection_name, "_4_genetics_summary_table.csv"),
    row.names = F,
    na = ""
  )
  
  write.csv(
    x_div_df,
    paste0(outdir, "/", collection_name,"/", collection_name, "_4_GI_table.csv"),
    row.names = F,
    na = ""
  )
  ##############################################################################
  ##############################################################################
  #  row.names(x_country)  <- x_country$country
  x_country[nrow(x_country)+1,] <- NA
  x_country$country[nrow(x_country)] <- "Subset"
  #row.names(x_country)  <- x_country$country
  x_country$pREC_wild[nrow(x_country)] <- x_div_df[1,2]
  x_country$pREC_landrace[nrow(x_country)] <- x_div_df[2,2]
  x_country$pSeq_wild[nrow(x_country)] <- x_div_df[1,3]
  x_country$pSeq_landrace[nrow(x_country)] <- x_div_df[2,3]
  x_country$IUCN_wild[nrow(x_country)] <- x_div_df[1,4]
  x_country$IUCN_landrace[nrow(x_country)] <- x_div_df[1,4]
  x_country$usability_index_wild[nrow(x_country)] <- x_div_df[1,5]
  x_country$usability_index_landrace[nrow(x_country)] <- x_div_df[2,5]
  ################################################################################
  x_country$REGION <- NA
  x_country$REGION <- countrycode(x_country$country,
                                  origin = 'iso3c',
                                  destination = 'region')
  
  x_country$country[which(x_country$country == "NA1")] <- "No_country"
  row.names(x_country)  <- x_country$country
  x_country$REGION[which(is.na(x_country$REGION))] <- "N/A"
  x_country$REGION[which(x_country$country == "Subset")] <- "Subset"
  x_country$REGION[which(x_country$country == "No_country")] <- "N/A"
  x_country$REGION[which(x_country$country == "SCG")] <- "Europe & Central Asia"
  x_country$REGION[which(x_country$country == "YUG")] <- "Europe & Central Asia"
  x_country$REGION[which(x_country$country == "ZAR")] <- "Sub-Saharan Africa"
  x_country$REGION[which(x_country$country == "BUR")] <- "East Asia & Pacific"
  
  
  write.csv(
    x_country,
    paste0(outdir, "/", collection_name,"/", collection_name, "_4_GI_country_table.csv"),
    row.names = F,
    na = ""
  )
  ############################################################################
  message(paste0("Processing data for plots"))
  ############################################################################
  x_res_W <- 
    x_country[, c(
      "pREC_wild",
      "pSeq_wild",
      "IUCN_wild",
      "usability_index_wild",
      "REGION")]
  
  colnames(x_res_W) <- c("Proportion of accessions sequenced (mean)",
                         "Proportion of taxa sequenced (mean)",
                         "Proportion of IUCN taxa threatened (mean)",
                         "Usability index",
                         "REGION"
  )
  
  x_res_W <- x_res_W[complete.cases(x_res_W),]
  x_res_W$REGION <- factor(x_res_W$REGION)
  
  #removing variables without variability
  if(nrow(x_res_W)>0){
    x_del <- list()
    for(i in 1:(ncol(x_res_W)-1)){
      #print(i)
      if(var(x_res_W[,i])==0){
        x_del[[i]] <- i
      } else {
        x_del[[i]] <- NA
      }
    }
    x_del <- na.omit(unlist(x_del))
    if(sum(x_del,na.rm = T)>0){
      message(paste("removed:",colnames(x_res_W)[x_del]))
      x_res_W <- x_res_W[,-x_del]
    }
  }
  #######################################
  x_res_L <- 
    x_country[, c(
      "pREC_landrace",
      "pSeq_landrace",
      "IUCN_landrace",
      "usability_index_landrace",
      "REGION")]
  
  colnames(x_res_L) <- c("Proportion of accessions sequenced (mean)",
                         "Proportion of taxa sequenced (mean)",
                         "Proportion of IUCN taxa threatened (mean)",
                         "Usability index",
                         "REGION"
  )
  
  x_res_L <- x_res_L[complete.cases(x_res_L),]
  x_res_L$REGION <- factor(x_res_L$REGION)
  
  if(nrow(x_res_L)>0){
    #removing variables without variability
    x_del <- list()
    for(i in 1:(ncol(x_res_L)-1)){
      #print(i)
      if(var(x_res_L[,i])==0){
        x_del[[i]] <- i
      } else {
        x_del[[i]] <- NA
      }
    }
    x_del <- na.omit(unlist(x_del))
    if(sum(x_del,na.rm = T)>0){
      message(paste("removed:",colnames(x_res_L)[x_del]))
      x_res_L <- x_res_L[,-x_del]
    }
  }
  #rm(i)
  ###
  ################################################################################
  #seed for NIPALS
  set.seed(1000)
  ################################################################################
  message("plotting results for wild")
  
  if(nrow(x_res_W)>=2){
    x_W_A = nipals(x_res_W[, -ncol(x_res_W)], comps = 3,scaled = T)
    x_res_W_NA <- x_res_W[!row.names(x_res_W) %in% c("No_country","Subset"),]
    # Create a scatter plot
    ind_coords <- data.frame(Dim.1=x_W_A$scores[,1],Dim.2=x_W_A$scores[,2],
                             REGION=x_res_W$REGION,countries=row.names(x_W_A$score))
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
      ggtitle(paste0("Usability index for ","Crop Wild Relatives (", nrow(x_res_W_NA), " countries)"))+
      
      geom_segment(data = var_coords, aes(x = 0, y = 0, xend = t1*(r * 0.7), yend = t2*(r * 0.7)),colour = "steelblue",
                   arrow = arrow(length = unit(0.2, "cm")))+
      
      geom_text_repel(
        data          = var_coords,
        mapping       = aes(t1*(r * 0.7), t2*(r * 0.7), label = vars),
        size = 3,    # Repel away from the left edge, not from the right.
        color = "steelblue", 
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
        "_",
        "usability_WILD.png"
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
  ################################################################################
  message("plotting results for landrace")
  if(nrow(x_res_L)>=2){
    
    REGION <- x_res_L$REGION
    if(length(unique(x_res_L[,1]))>1){
      x_W_A = nipals(x_res_L[, -ncol(x_res_L)], comps = 3,scaled = T)
      
      x_res_L_NA <- x_res_L[!row.names(x_res_L) %in% c("No_country","Subset"),]
      # Create a scatter plot
      ind_coords <- data.frame(Dim.1=x_W_A$scores[,1],Dim.2=x_W_A$scores[,2],
                               REGION=x_res_L$REGION,countries=row.names(x_W_A$score))
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
        ggtitle(paste0("Usability index for ","landraces (", nrow(x_res_L_NA), " countries)"))+
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
          max.overlaps = 1000,#200,  
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
          max.overlaps = 300,#200,
          max.iter = 10000,show.legend = NA,box.padding = 0.05,
          point.padding = 0.02, 
          nudge_x = .1,
          nudge_y = .1,
          segment.linetype = 6,
          segment.curvature = -1e-20,
          hjust= 1
        ) +
        
        
        theme_minimal()
      
      ggplot2::ggsave(
        paste0(
          outdir,
          "/",
          collection_name,
          "/",
          collection_name,
          "_",
          
          "usability_land.png"
        ),
        PCA_LAND,
        dpi = 300,
        units = "in",
        width = 10,
        height = 8
      )
    } else {
      message("NO NIPALS AVAILABLE FOR LANDRACE")
    }
  } else {
    warning("NO LANDRACES AVAILABLE")
  }
  ###############################################################################
  
  message("DONE!")
  return(results)
}


################################################################################
dir <- "D:/OneDrive - CGIAR/GERMPLASM_INDEX"
#outdir
outdir <- paste0(dir, "/BEANS/RESULTS")
################################################################################
numCores <- 10
################################################################################
collection_name <- "beans"
################################################################################
#loading acccessions sequenced
accessions_df <- as.data.frame(
  readxl::read_xlsx(
    #"D:/OneDrive - CGIAR/GERMPLASM_INDEX/CIAT Data/ListaGenotipadoYuca_frijol_forrajesSept2024.xlsx",
    "D:/OneDrive - CGIAR/GERMPLASM_INDEX/CIAT Data/ListaGenotipadoYuca_frijol_forrajesOct2024.xlsx",
    sheet = collection_name
  )
)
x1 <- genetics_ind_function(outdir, collection_name, numCores, accessions_df)
################################################################################
################################################################################
################################################################################
collection_name <- "cassava"
################################################################################
#loading acccessions sequenced
accessions_df <- as.data.frame(
  readxl::read_xlsx(
    "D:/OneDrive - CGIAR/GERMPLASM_INDEX/CIAT Data/ListaGenotipadoYuca_frijol_forrajesOct2024.xlsx",
    sheet = collection_name
  )
)
x2 <- genetics_ind_function(outdir, collection_name, numCores, accessions_df)
################################################################################
################################################################################
################################################################################
collection_name <- "forages"
################################################################################
#loading acccessions sequenced
accessions_df <- as.data.frame(
  readxl::read_xlsx(
    "D:/OneDrive - CGIAR/GERMPLASM_INDEX/CIAT Data/ListaGenotipadoYuca_frijol_forrajesOct2024.xlsx",
    sheet = collection_name
  )
)

accessions_df$Accession <- gsub("^F", "", accessions_df$Accession) 
x3 <- genetics_ind_function(outdir, collection_name, numCores, accessions_df)

# x1$indicator1_seq_prop
# x2$indicator1_seq_prop
# x3$indicator1_seq_prop


x_ind <- rbind(x1$indexes,x2$indexes)
x_ind <- rbind(x_ind,x3$indexes)


write.csv(
  x_ind,
  paste0(outdir, "/", "CIAT_4_GI_table.csv"),
  row.names = F,
  na = ""
)
