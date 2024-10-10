require(data.table)
require(readxl)
require(parallel)


genetics_ind_function <- function(outdir,
                                  collection_name,
                                  numCores,
                                  accessions_df) {
  message(paste0("Processing collection: ",collection_name,"... "))
  message("Loading results of 1_Taxonomic_index")
  ################################################################################
  subsets <- readRDS(paste0(outdir, "/", collection_name, "/", collection_name, "_subsets_new_1.RDS"))
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
    x <- tapply(passport_data_w$SEQ_FLAG[which(passport_data_w$SEQ_FLAG == 1)], passport_data_w$NAME[which(passport_data_w$SEQ_FLAG ==
                                                                                                             1)], length)
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
    x <- tapply(passport_data_l$SEQ_FLAG[which(passport_data_l$SEQ_FLAG == 1)], passport_data_l$NAME[which(passport_data_l$SEQ_FLAG ==
                                                                                                             1)], length)
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
    x <- tapply(passport_data_h$SEQ_FLAG[which(passport_data_h$SEQ_FLAG == 1)], passport_data_h$NAME[which(passport_data_h$SEQ_FLAG ==
                                                                                                             1)], length)
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
   IND <- weighted.mean(status_total$prop_seq, 
                        status_total$prop_rec_collection,
                        na.rm=T)
   
   IND_W <- weighted.mean(status_total$prop_seq[which(status_total$status=="wild")], 
                        status_total$prop_rec_collection[which(status_total$status=="wild")],
                        na.rm=T)
   
   IND_L <- weighted.mean(status_total$prop_seq[which(status_total$status=="landrace")], 
                          status_total$prop_rec_collection[which(status_total$status=="landrace")],
                          na.rm=T)
   
   IND_H <- weighted.mean(status_total$prop_seq[which(status_total$status=="hybrid")], 
                          status_total$prop_rec_collection[which(status_total$status=="hybrid")],
                          na.rm=T)
   
   # IND_a <- weighted.mean(status_total$prop_seq, 
   #                      status_total$prop_rec_collection,
   #                      na.rm=T)
   IND2 <- sum(status_total$sp_status,na.rm = T)/nrow(status_total)
   IND2_W <- sum(status_total$sp_status[which(status_total$status=="wild")],na.rm = T)/
     nrow(status_total[which(status_total$status=="wild"),])
   IND2_L <- sum(status_total$sp_status[which(status_total$status=="landrace")],na.rm = T)/
     nrow(status_total[which(status_total$status=="landrace"),])
   IND2_H <- sum(status_total$sp_status[which(status_total$status=="hybrid")],na.rm = T)/
     nrow(status_total[which(status_total$status=="hybrid"),])
  message("Saving results")
  
  x_div_df <- data.frame(matrix(ncol=4,nrow = 4))
  colnames(x_div_df) <- c("status","index_1", "index_2","collection")
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
  x_div_df[,4] <- collection_name
  results <- list(summary_table = status_total, 
                  indexes= x_div_df)
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
