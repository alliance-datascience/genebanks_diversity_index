require(data.table)
library(parallel)
library(factoextra)
library(fossil)
library(plsdepot)
################################################################################
genetics_dist_function <- function(outdir,
                                   collection_name,
                                   numCores,
                                   gen_data) {
  
  #checking if there is genetics distance matrix!
  if (is.null(gen_data)) {
    stop(paste("Please include a genetic distance matrix for ",collection_name))
  }
  
  #Checking if the genetic distance is squared!
  #getting dimensions
  dim_gen <- dim(gen_data)
  
  
  if (dim_gen[1]!=dim_gen[2]){
    stop(paste("Please include a valid genetic distance matrix (squared matrix!) for ",collection_name))
  }
  message(paste("Reading  previous results for ", collection_name))
  #load(paste0(outdir, "/", collection_name, "_subsets_1.RData"))
  subsets <- readRDS(paste0(
    outdir,
    "/",
    collection_name,
    "/",
    collection_name,
    "_subsets_new_1.RDS"
  ))
  ##############################################################################
  #calling files to join information with genetic data 
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
  
  ##############################################################################
  #ensuring that data is available. Thus, if there is no wild or landrance, or hybrids there is no count table!
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
  #checking cores to use in approach
  message("Checking cores to use parallel approach")
  x_det <- parallel::detectCores()
  if (numCores > x_det) {
    stop(
      "Number of cores exceed the maximum allowed by the machine, use a coherent number of cores such as four"
    )
  }
  
  ################################################################################
  message("Matching sequenced accessions using ACCENUMB")
  ###############################################################################

  ################################################################################
  #sequenced for wild records
  ################################################################################
  #Assigning sequenced flag for wild records using parallel library (0 and 1)
  #for lack and match of accession in genetic distance respectively
  
  message("Sequenced accessions for wild records")
  
  cl <- parallel::makeCluster(numCores)
  parallel::clusterExport(cl,
                          varlist = c("passport_data_w", "gen_data"),
                          envir = environment())
  
  res <- parallel::parLapplyLB(
    cl,
    X = seq_len(nrow(passport_data_w)),
    fun = function (i) {
      if (isTRUE(passport_data_w$ACCENUMB[[i]] %in% row.names(gen_data))) {
        x <- 1
      } else {
        x <- 0
      }
    }
  )
  
  parallel::stopCluster(cl)
  rm(cl)
  
  
  # obtaining data that match with genetic data
  acc_to_W <- data.frame(ACCESSION = passport_data_w$ACCENUMB, EVAL = unlist(res))
  acc_to_W <- acc_to_W[which(acc_to_W$EVAL == 1), ]
  
  #subsetting genetic data
  gen_data_subset_W <- gen_data[, acc_to_W$ACCESSION]
  gen_data_subset_W <- as.data.frame(gen_data_subset_W[row.names(gen_data_subset_W) %in% acc_to_W$ACCESION, ])
  
  
  #calculating genetic distance summaries only if there is information
  if(nrow(gen_data_subset_W)>0){
    #obtaining median of genetic distances per individual
    x_median_W <- list()
    for (i in 1:ncol(gen_data_subset_W)) {
      gen_data_subset_W[, i] <- as.numeric(as.character(gen_data_subset_W[, i]))
      x_median_W[[i]] <- median(gen_data_subset_W[, i], na.rm = T)
    }
    rm(i)  
    
    #filling 
    acc_to_W$GEN_DISTANCE <- NA
    acc_to_W$GEN_DISTANCE <- unlist(x_median_W)
    passport_data_w$GEN_DISTANCE <- NA
    
    for (i in 1:nrow(acc_to_W)) {
      passport_data_w$GEN_DISTANCE[which(passport_data_w$ACCENUMB == acc_to_W$ACCESION[[i]])] <- acc_to_W$GEN_DISTANCE[[i]]
    }
    rm(i)
    
    #obtaining summary per country (all species together)
    dist_per_country_W <- data.frame(COUNTRY=colnames(df_wild),DIST_GEN=NA)#,MAD_GEN=NA)
    for(i in 1:ncol(df_wild)){
      #i <- 1
      dist_per_country_W$DIST_GEN[[i]] <- median(passport_data_w$GEN_DISTANCE
                                                 [which(passport_data_w$ORIGCTY==colnames(df_wild)[[i]])],na.rm = T)
      
      dist_per_country_W$SUBSET <- "wild"
      #dist_per_country$MAD_GEN[[i]] <- mad(passport_data_l$GEN_DISTANCE
      #[which(passport_data_W$ORIGCTY==colnames(df_wild)[[i]])],na.rm = T)
    }
    rm(i)
    
  } else {
    acc_to_W <- data.frame(ACCESSION = NULL, EVAL = NULL,GEN_DISTANCE=NULL)
    dist_per_country_W <- NULL
    
    warning("NO GENETIC DATA FOR CROP WILD RELATIVES!")
  }

  rm(res)
  ################################################################################
  #preparing summary file for wild
  
  #summary file for landrace
  
  if (nrow(acc_to_W)>0){
    stat <- strsplit(row.names(df_wild), "_")
    
    wild_gen_sum <- data.frame(
      taxon = unlist(lapply(stat, `[[`, 1)),
      status = "landrace",
      records = rowSums(df_wild),
      records_seq = NA,
      prop_seq = NA,
      gen_distance = NA
    )
    
    for (i in 1:nrow(wild_gen_sum)) {
      x_i <- passport_data_w$GEN_DISTANCE[which(passport_data_w$NAME == wild_gen_sum$taxon[[i]])]
      wild_gen_sum$gen_distance[[i]] <- median(x_i,na.rm = T)#exp(mean(log(x_i), na.rm = T))
      wild_gen_sum$records_seq[[i]] <- length(na.omit(x_i))
      wild_gen_sum$prop_seq[[i]] <- wild_gen_sum$records_seq[[i]] / wild_gen_sum$records[[i]]
    }
    rm(i)
    
  
} else {
  wild_gen_sum <- data.frame(
    taxon = NA,
    status = "wild",
    records = NA,
    records_seq = NA,
    prop_seq = NA,
    gen_distance = NA
  )
}


  
  
  ################################################################################
  #sequenced for landraces
  ################################################################################
  #Assigning sequenced flag for wild records using parallel library (0 and 1)
  #for lack and match of accession in genetic distance respectively
  
  message("Sequenced accessions for landrace records")
  
  cl <- parallel::makeCluster(numCores)
  parallel::clusterExport(cl,
                          varlist = c("passport_data_l", "gen_data"),
                          envir = environment())
  
  res <- parallel::parLapplyLB(
    cl,
    X = seq_len(nrow(passport_data_l)),
    fun = function (i) {
      if (isTRUE(passport_data_l$ACCENUMB[[i]] %in% row.names(gen_data))) {
        x <- 1
      } else {
        x <- 0
      }
    }
  )
  
  parallel::stopCluster(cl)
  rm(cl)
  
  
  # obtaining data that match with genetic data
  acc_to_L <- data.frame(ACCESION = passport_data_l$ACCENUMB, EVAL = unlist(res))
  acc_to_L <- acc_to_L[which(acc_to_L$EVAL == 1), ]
  
  #subsetting genetic data
  gen_data_subset_L <- gen_data[, acc_to_L$ACCESION]
  gen_data_subset_L <- as.data.frame(gen_data_subset_L[row.names(gen_data_subset_L) %in% acc_to_L$ACCESION, ])
  
  
  #calculating genetic distance summaries only if there is information
  if(nrow(gen_data_subset_L)>0){
    #obtaining median of genetic distances per individual
    x_median_L <- list()
    for (i in 1:ncol(gen_data_subset_L)) {
      gen_data_subset_L[, i] <- as.numeric(as.character(gen_data_subset_L[, i]))
      x_median_L[[i]] <- median(gen_data_subset_L[, i], na.rm = T)
    }
    rm(i)  
    
    #filling 
    acc_to_L$GEN_DISTANCE <- NA
    acc_to_L$GEN_DISTANCE <- unlist(x_median_L)
    passport_data_l$GEN_DISTANCE <- NA
    
    for (i in 1:nrow(acc_to_L)) {
      passport_data_l$GEN_DISTANCE[which(passport_data_l$ACCENUMB == acc_to_L$ACCESION[[i]])] <- acc_to_L$GEN_DISTANCE[[i]]
    }
    rm(i)
    
    
  } else {
    acc_to_L <- data.frame(ACCESION = NA, EVAL = NA,GEN_DISTANCE=NA)
    warning("NO GENETIC DATA FOR LANDRACES!")
  }
  
  rm(res)
  
 
  #summary file for landrace
  if(nrow(gen_data_subset_L)>0){
    stat <- strsplit(row.names(df_landrace), "_")
    
    land_gen_sum <- data.frame(
      taxon = unlist(lapply(stat, `[[`, 1)),
      status = "landrace",
      records = rowSums(df_landrace),
      records_seq = NA,
      prop_seq = NA,
      gen_distance = NA
    )
    
    for (i in 1:nrow(land_gen_sum)) {
      x_i <- passport_data_l$GEN_DISTANCE[which(passport_data_l$NAME == land_gen_sum$taxon[[i]])]
      land_gen_sum$gen_distance[[i]] <- median(x_i,na.rm = T)#exp(mean(log(x_i), na.rm = T))
      land_gen_sum$records_seq[[i]] <- length(na.omit(x_i))
      land_gen_sum$prop_seq[[i]] <- land_gen_sum$records_seq[[i]] / land_gen_sum$records[[i]]
    }
    rm(i)
    #obtaining summary per country (all species together)
    dist_per_country_L <- data.frame(COUNTRY=colnames(df_landrace),DIST_GEN=NA)#,MAD_GEN=NA)
    for(i in 1:ncol(df_landrace)){
      #i <- 1
      dist_per_country_L$DIST_GEN[[i]] <- median(passport_data_l$GEN_DISTANCE
                                               [which(passport_data_l$ORIGCTY==colnames(df_landrace)[[i]])],na.rm = T)
      dist_per_country_L$SUBSET <- "landrace"
      #dist_per_country$MAD_GEN[[i]] <- mad(passport_data_l$GEN_DISTANCE
      #[which(passport_data_l$ORIGCTY==colnames(df_landrace)[[i]])],na.rm = T)
    }
    rm(i)
    
    
  } else {
    land_gen_sum <- data.frame(
      taxon = NA,
      status = "landrace",
      records = NA,
      records_seq = NA,
      prop_seq = NA,
      gen_distance = NA
    )
    
    dist_per_country_L <- NULL
  }

  
  ################################################################################
  #Joining information
  ################################################################################
  #summary files
  distance_summary <- rbind(wild_gen_sum,land_gen_sum)
  distance_summary$collection <- collection_name
  ################################################################################
  #country summary files
  dist_per_country <- rbind(dist_per_country_W,dist_per_country_L)
    ################################################################################
  #saving values per collection subset
  write.csv(
    distance_summary,
    paste0(outdir, "/", collection_name,"/", collection_name, "_4_1_Genetic_distance_summary.csv"),
    row.names = F,
    na = ""
  )
  #saving values per country  
  write.csv(
    dist_per_country,
    paste0(outdir, "/", collection_name,"/", collection_name, "_4_1_Genetic_distance_country.csv"),
    row.names = F,
    na = ""
  )
  return(distance_summary)
}

################################################################################
################################################################################
################################################################################
dir <- "D:/OneDrive - CGIAR/GERMPLASM_INDEX"
#dir <- "D:/ONEDRIVE/cgiar/OneDrive - CGIAR/GERMPLASM_INDEX"
################################################################################
#outdir
outdir <- paste0(dir, "/BEANS/RESULTS")
################################################################################
#number of cores to use in the function
numCores <- 4
################################################################################
###
#RUNS
################################################################################
#Collection name
collection_name <- "cassava"
#loading genetics distances data (This must be a distance matrix!) IF DATA IS NOT AVAILABLE USE NULL!

genetics_data_cassava <- as.matrix(
  data.table::fread(
    "D:/OneDrive - CGIAR/GERMPLASM_INDEX/CIAT Data/GENETICS/distances_glm_snp5302samples&7180SNPs_accessionName.csv"
  )
)
#Ensuring that matrix is squared!
row.names(genetics_data_cassava) <- genetics_data_cassava[, 1]
genetics_data_cassava[, 1] <- NA
genetics_data_cassava <- genetics_data_cassava[, -c(1)]
#excluding diagonal 
diag(genetics_data_cassava) <- NA
gen_data <- genetics_data_cassava

x <- genetics_dist_function(outdir=outdir,
                       collection_name=collection_name,
                       numCores=numCores,
                       gen_data=gen_data)
################################################################################
collection_name <- "beans"
#loading genetics distances data (This must be a distance matrix!) IF DATA IS NOT AVAILABLE USE NULL!
gen_data <- NULL
genetics_dist_function(outdir=outdir,
                       collection_name=collection_name,
                       numCores=numCores,
                       gen_data=gen_data)
################################################################################
collection_name <- "forages"
#loading genetics distances data (This must be a distance matrix!) IF DATA IS NOT AVAILABLE USE NULL!
gen_data <- NULL
genetics_dist_function(outdir=outdir,
                       collection_name=collection_name,
                       numCores=numCores,
                       gen_data=gen_data)
################################################################################