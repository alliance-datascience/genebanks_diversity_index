require(data.table)
require(readr)
require(stringr)
require(tidyr)
require(psych)
#install.packages("ineq")
#library(ineq)
#install.packages("igraph")
#require(igraph)
library(abdiv)
#library(vegan)
library(matrixStats)
#require(factoextra)
#require(ecodist)
require(DescTools)
require(httr)
require(WorldFlora)
#require(parallel)
################################################################################
# minMax <- function(x) {
#   (x - min(x)) / (max(x) - min(x))
# }

#require(strigdist)

#read taxonomy for Phaseolus
#taxonomy <- read.csv(paste0(dir,"/","BEANS/GRIN.csv"),na.strings = NA)
#GRIN

################################################################################
################################################################################
index_function <- function(passport_data_orig,
                           tax,
                           geo,
                           tax_geo,
                           outdir,
                           API,
                           collection_name) {
  
  message("Preprocessing: creating subfolder")
  if(!dir.exists(paste0(outdir, "/", collection_name))){
    dir.create(paste0(outdir, "/", collection_name))
  }
  
  message("PART 0: ommiting plant breeding material and weedy")
  
  if(!file.exists(paste0(outdir, "/", collection_name,"/",collection_name, "_subsets_new_1.RDS"))){
  #ommiting plant breeding material
  passport_data_orig <- passport_data_orig[!passport_data_orig$SAMPSTAT %in% c(200, 412, 413, 414, 415, 416, 420, 421, 422, 423, 500, 600), ]
  
  
  ################################################################################
  message("PART1: CREATING COUNT TABLE FOR COUNTRIES AND TAXA")
  
  #changing countries without ISO3 to NA1
  message("changing countries without ISO3 to NA1")
  passport_data_orig$ORIGCTY[which(is.na(passport_data_orig$ORIGCTY))] <- "NA1"
  ################################################################################
  #Using name from genus, species, subtaxa instead of GRIN NAME
  message("Compiling name from GENUS, SPECIES, and SUBTAXA")
  NAME <- list()
  for (i in 1:nrow(passport_data_orig)) {
    if (is.na(passport_data_orig$GENUS[[i]]) &
        is.na(passport_data_orig$SPECIES[[i]]) &
        is.na(passport_data_orig$SUBTAXA[[i]])) {
      NAME[[i]] <- NA
    } else if (is.na(passport_data_orig$SUBTAXA[[i]]) &
               !is.na(passport_data_orig$SPECIES[[i]]) &
               !is.na(passport_data_orig$GENUS[[i]])) {
      NAME[[i]] <- paste(passport_data_orig$GENUS[[i]],
                         passport_data_orig$SPECIES[[i]])
    } else if (!is.na(passport_data_orig$GENUS[[i]]) &
               is.na(passport_data_orig$SPECIES[[i]])) {
      NAME[[i]] <- paste(passport_data_orig$GENUS[[i]], "spp.")
    } else {
      NAME[[i]] <- paste(
        passport_data_orig$GENUS[[i]],
        passport_data_orig$SPECIES[[i]],
        passport_data_orig$SUBTAXA[[i]]
      )
    }
  }
  rm(i)
  
  NAME <- unlist(NAME)
  passport_data_orig$NAME <- trimws(NAME)
  #Add check to avoid NA
  message("Obtaining species names with issues")
  passport_data_orig$NAME[which(is.na(passport_data_orig$NAME))] <- "CHECK"
  tax_all <- unique(passport_data_orig$NAME)
  
  
  
  message("detecting hybrid taxa")
  # Calling str_detect() function
  hybrid_table <- data.frame(taxa = tax_all,
                             hybrid = stringr::str_detect(tax_all, " x "))
  hybrid_table_WILD <- hybrid_table$taxa[hybrid_table$hybrid == TRUE]
  SPP <- hybrid_table$taxa[hybrid_table$hybrid == FALSE]
  ################################################################################
  message("CREATING COUNT TABLE FOR COUNTRIES AND WILD SPECIES")
  ################################################################################
  #subsetting to passport data without being landraces
  passport_data_w <- passport_data_orig[passport_data_orig$SAMPSTAT %in% c(100, 110, 120, 130, 999, NA), ] #removing 200 (weedy)
  #passport_data_w2 <- passport_data_orig[!passport_data_orig$SAMPSTAT %in% c(100, 110, 120, 130, 200,999), ]
  
  passport_data_w <- passport_data_w[passport_data_w$NAME %in% SPP, ]
  
  #passport_data <- passport_data[which(passport_data$GENUS=="Phaseolus"),]
  if (nrow(passport_data_w) > 0) {
    #obtaining unique countries and species
    coun <- unique((passport_data_orig$ORIGCTY))
    spp <- unique((passport_data_w$NAME))
    
    df <- data.frame(matrix(nrow = length(spp), ncol = length(coun)))
    colnames(df) <- coun
    row.names(df) <- spp
    for (i in 1:length(spp)) {
      for (j in 1:length(coun)) {
        x <- passport_data_w[which(passport_data_w$NAME == spp[[i]] &
                                     passport_data_w$ORIGCTY == coun[[j]]), ]
        
        df[i, j] <- nrow(x)
      }
    }
    df_wild <- df
    rm(df)
    #adding wild to names
    row.names(df_wild) <- paste0(row.names(df_wild), "_", "wild")
  } else {
    df_wild <- NA
  }
  ################################################################################
  #subsetting to passport data  being landraces
  message("CREATING COUNT TABLE FOR COUNTRIES AND LANDRACES")
  passport_data_l <- passport_data_orig[passport_data_orig$SAMPSTAT %in% c(300), ]
  passport_data_l <- passport_data_l[passport_data_l$NAME %in% SPP, ]
  
  if (nrow(passport_data_l) > 0) {
    #obtaining unique countries and species
    coun <- unique((passport_data_orig$ORIGCTY))
    spp <- unique((passport_data_l$NAME))
    df <- data.frame(matrix(nrow = length(spp), ncol = length(coun)))
    colnames(df) <- coun
    row.names(df) <- spp
    for (i in 1:length(spp)) {
      for (j in 1:length(coun)) {
        x <- passport_data_l[which(passport_data_l$NAME == spp[[i]] &
                                     passport_data_l$ORIGCTY == coun[[j]]), ]
        
        df[i, j] <- nrow(x)
      }
    }
    df_landrace <- df
    rm(df)
    #adding wild to names
    row.names(df_landrace) <- paste0(row.names(df_landrace), "_", "landrace")
  } else {
    df_landrace <- NA
  }
  
  ################################################################################
  ###
  #subsetting to passport data  being landraces
  message("CREATING COUNT TABLE FOR COUNTRIES AND hybrids")
  passport_data_h <- passport_data_orig[passport_data_orig$NAME %in% hybrid_table_WILD, ]
  
  if (nrow(passport_data_h) > 0) {
    #obtaining unique countries and species
    coun <- unique((passport_data_orig$ORIGCTY))
    spp <- unique((passport_data_h$NAME))
    df <- data.frame(matrix(nrow = length(spp), ncol = length(coun)))
    colnames(df) <- coun
    row.names(df) <- spp
    for (i in 1:length(spp)) {
      for (j in 1:length(coun)) {
        x <- passport_data_h[which(passport_data_h$NAME == spp[[i]] &
                                     passport_data_h$ORIGCTY == coun[[j]]), ]
        
        df[i, j] <- nrow(x)
      }
    }
    df_hybrid <- df
    rm(df)
    #adding wild to names
    row.names(df_hybrid) <- paste0(row.names(df_hybrid), "_", "hybrid")
  } else {
    df_hybrid <- NA
  }
  
  
  ###

  
  ################################################################################
  #validating counts for wild and landraces
  # common_passport <- rbind(passport_data_l, passport_data_w)
  #
  # counts <- (sum(df_wild) + sum(df_landrace))
  # nrow(common_passport)
  # taxon_names <- unique((common_passport$NAME))
  
  ################################################################################
  #count matrix where taxa are rows and countries are cols
  message("joining wild and landraces count matrices in one")
  if (is.data.frame(df_wild) == TRUE & !is.data.frame(df_landrace)) {
    df_total <- df_wild
  } else if (!is.data.frame(df_wild) &
             is.data.frame(df_landrace) == TRUE) {
    df_total <- df_landrace
  } else if (is.data.frame(df_wild) == TRUE &
             is.data.frame(df_landrace) == TRUE) {
    df_total <- rbind(df_wild, df_landrace)
  } else {
    stop("NO RECORDS, check the data")
  }
  
  message("joining wild and landraces to hybrids in one")
  if (is.data.frame(df_total) == TRUE & !is.data.frame(df_hybrid)) {
    df_total <- df_total
  } else {
    df_total <- rbind(df_total, df_hybrid)
  }
  
  
  
  
  #sum(df_total)
  df_total <- df_total[rowSums(df_total) > 0, ]
  df_total2 <- t(df_total)
  
  write.csv(
    df_total2,
    paste0(outdir, "/", collection_name, "/", collection_name,"_count_matrix_countries_1.csv"),
    row.names = T,
    na = ""
  )
  message(paste0(
    "Saved:",
    outdir,
    "/",
    collection_name,
    "/", 
    collection_name,
    "_count_matrix_countries_1.csv"
  ))
  

  #df_total
  ################################################################################
  #obtaining distances for countries
  #removing records without countries data
  # x <- df_total[, !colnames(df_total) %in% "NA1"]
  # #ommiting countries with 0 counts
  # x <- x[, colSums(x,na.rm = T) > 0]
  # #ommiting taxa with 0 counts
  #
  # x <- x[rowSums(x,na.rm = T) > 0, ]
  # #creating presence-absence matrix
  # x2 <- x
  # x2[x2 > 0] <- 1
  # #######################
  x_total <- df_total
  x_total <- x_total[, colSums(x_total, na.rm = T) > 0]
  #ommiting taxa with 0 counts
  
  x_total <- x_total[rowSums(x_total, na.rm = T) > 0, ]
  #creating presence-absence matrix
  x2_total <- x_total
  x2_total[x2_total > 0] <- 1
  
  x2_non_na <- x2_total[, !colnames(x2_total) %in% "NA1"]
  
  
  ################################################################################
  #removing NA from counts
  df_total_nonNA <- df_total
  df_total_nonNA <- df_total_nonNA[, !colnames(df_total_nonNA) %in% "NA1"]
  #df that are NA in country
  df_total_NA <- df_total
  df_total_NA <- df_total_NA[, colnames(df_total_NA) %in% "NA1"]
  ################################################################################
  #saving subsets to use for further steps
  subsets <- list(
  passport_data_w,
  passport_data_l,
  passport_data_h,
  df_total,
  df_total_nonNA,
  df_total_NA,
  x_total,
  x2_total,
  x2_non_na,
  df_wild,
  df_landrace,
  df_hybrid
  )
  message("saving  previous results")
  #save.image(paste0(outdir, "/", collection_name, "_subsets_1.RData"))
  saveRDS(subsets,paste0(outdir, "/", collection_name, "/", collection_name, "_subsets_new_1.RDS"))
  } else {
    message("reading  previous results")
    #load(paste0(outdir, "/", collection_name, "_subsets_1.RData"))
    subsets <- readRDS(paste0(outdir, "/", collection_name, "/", collection_name,"_subsets_new_1.RDS"))
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
  }
  
  
  ################################################################################
  message("CREATING SUMMARY TABLE FOR COUNTRIES AND TAXA")
  message("USING ONLY RECORDS WITH ORIGCTY INFORMATION")
  #creating matrix to analyze countries uniqueness per taxon
  dist_prop <- as.data.frame(matrix(nrow = nrow(x2_total), ncol = 15))
  colnames(dist_prop) <- c(
    "Taxa",
    #####################
    "countries/countries_nontaxa", # countries where the taxa is available compare with the countries of other taxa different
    "countries/countries_total", #countries where taxa is available/total countries of the collection
    "ncountries",
    "ncountries_total_collection",
    ######################
    "countries/countries_nontaxa_noNA",
    "countries/countries_total_noNA",
    "ncountries_noNA",
    "ncountries_total_collection_noNA",
    ######################
    "ncounts_Country",
    "ncounts_nonnaCountry_noNA",
    "ncounts_NA",
    "prop_taxa_in_collection",
    "prop_taxa_in_collection_noNA",
    "prop_taxa_in_collectionNA"
    
  )#row.names(x2)
  row.names(dist_prop) <- row.names(x2_total)
  
  #for each species
  for (i in 1:nrow(x2_total)) {
   # i <- 48#1
    #i <- 9
    #############################################################################
    #splitting each species
    a <- x2_total[i, ]
    #splitting for the rest of counts which are not for the taxon i
    b <- x2_total[-i, ]
    #ensuring using the correct species to count NAs!
    spp <- row.names(a)
    message(paste0("processing...", spp))
    #############################################################################
    #creating a dummy data frame
    df <- data.frame(
      country = colnames(x2_total),
      a = colSums(a),
      b = colSums(b),
      c = colSums(x2_total)
    )
    # changing sums to presence and obtain countries collecting
    #############################################################################
    #splitting each species
    a_nonna <- a[, !colnames(a) %in% "NA1"]
    #splitting for the rest of counts which are not for the taxon i
    b_nonna <- b[, !colnames(b) %in% "NA1"]
    #############################################################################
    #creating a dummy data frame
    df_nonna <- data.frame(
      country = colnames(x2_total)[!colnames(x2_total) %in% "NA1"],
      a = colSums(a_nonna),
      b = colSums(b_nonna),
      c = colSums(x2_total[, !colnames(x2_total) %in% "NA1"])
    )
    # changing sums to presence and obtain countries collecting
    df_nonna$b[df_nonna$b > 0] <- 1
    df_nonna$c[df_nonna$c > 0] <- 1
    #############################################################################
    
    #taxon name
    dist_prop[i, "Taxa"] <- row.names(a)
    #############################################################################
    #message("Calculating metrics including NA (countries) records ")
    #countries of taxon i / countries where !=i
    dist_prop[i, "countries/countries_nontaxa"] <- sum(df$a) / sum(df$b)
    #countries of taxon i / all countries
    dist_prop[i, "countries/countries_total"] <- sum(df$a) / ncol(x_total)
    #taxon i countries
    dist_prop[i, "ncountries"] <- sum(df$a)
    #all countries
    dist_prop[i, "ncountries_total_collection"] <- ncol(x_total)#sum(df$c > 0)
    ############################################################################
    #############################################################################
    #message("Calculating metrics for non NA (countries) records ")
    #countries of taxon i / countries where !=i
    dist_prop[i, "countries/countries_nontaxa_noNA"] <- sum(df_nonna$a) / sum(df_nonna$b)
    #countries of taxon i / all countries
    dist_prop[i,"countries/countries_total_noNA"] <- sum(df_nonna$a) / ncol(x2_non_na)
    #taxon i countries
    dist_prop[i, "ncountries_noNA"] <- sum(df_nonna$a)
    #all countries
    dist_prop[i, "ncountries_total_collection_noNA"] <- ncol(x2_non_na)#sum(df_nonna$c > 0)
  
    #############################################################################
    #numbers of records for taxon i ALL
    dist_prop[i, "ncounts_Country"] <- sum(df_total[i, ])
    #proportion of records for taxon i/ records of collection x 100 ALL
    dist_prop[i, "prop_taxa_in_collection"] <- (sum(df_total[i, ]) / sum(df_total))
    #numbers of records for taxon i NON NA
    dist_prop[i, "ncounts_nonnaCountry_noNA"] <- sum(df_total_nonNA[i, ])  
    #proportion of records for taxon i/ records of collection x 100 NON NA
    dist_prop[i, "prop_taxa_in_collection_noNA"] <- (sum(df_total_nonNA[i, ]) / sum(df_total_nonNA))
    #numbers of records for taxon i NA
    dist_prop[i, "ncounts_NA"] <-df_total_NA[i]
    #proportion of records for taxon i/ records of collection x 100 NA
    dist_prop[i, "prop_taxa_in_collectionNA"] <- df_total_NA[i]/ sum(df_total_NA,na.rm = T)
    #message("Calculating metrics for non NA (countries) records ")
    
    rm(df, a, b, a_nonna, b_nonna, df_nonna)
  }
  rm(i)
  ################################################################################
  message("Matching taxa with GRIN taxonomy")
  #obtaining taxon status (wild or landrace)
  stat <- strsplit(row.names(dist_prop), "_")
  #status
  dist_prop$status <- unlist(lapply(stat, `[[`, 2))
  #taxon name
  dist_prop$taxon <- unlist(lapply(stat, `[[`, 1))
  
  unique_taxa <- unique(dist_prop$taxon)
  
  
  dist_prop$GRIN <- NA
  
  ################################################################################
  ################################################################################
  message("OBTAINING SPECIES MATCHING WITH GRIN TAXONOMY AND NATIVE COUNTRIES (ISO3)")
  
  accepted_grin_tax_list <- list()
  native_countries_list <- list()
  for (i in 1:length(unique_taxa)) {
    #i <- 1
    tax_i <- tax[tax$name %in% unique_taxa[[i]], ]
    tax_i <- tax_i[!tax_i$synonym_code %in% c("I"),]#, "S"), ] #no invalidos no sinonimo
    if(nrow(tax_i)>1){
      tax_i <- tax_i[!tax_i$synonym_code %in% c("S"),]#, "S"), ] 
    }
    tax_i_cur <- tax[tax$current_taxonomy_species_id==tax_i$current_taxonomy_species_id, ]
    tax_i_cur <- tax_i_cur[tax_i_cur$taxonomy_species_id==tax_i_cur$current_taxonomy_species_id,]
    tax_i_cur <- tax_i_cur$name
    #message(nrow(tax_i))
    if (nrow(tax_i) > 0) {
      dist_prop$GRIN[which(dist_prop$taxon == unique_taxa[[i]])] <- 1
      accepted_grin_tax_list[[i]] <- data.frame(
        taxa = unique_taxa[[i]],
        grin_id = tax_i$taxonomy_species_id,
        grin_id_cur = tax_i$current_taxonomy_species_id,
        grin_taxon = tax_i_cur
      )
      
      tax_geo_i <- tax_geo[which(tax_geo$taxonomy_species_id == tax_i$current_taxonomy_species_id), ]
      
      tax_geo_i <- tax_geo_i$geography_id[which(tax_geo_i$geography_status_code ==
                                                  "n")]
      geo_i <- geo[geo$geography_id %in% tax_geo_i, ]
      geo_i <- unique(geo_i$country_code)
      
      ##########################################################################
      if (length(geo_i) > 0) {
        geo_i <- geo_i
      } else {
        geo_i <- NA
      }
      native_countries_list[[i]] <- data.frame(
        taxa = unique_taxa[[i]],
        grin_id = tax_i$taxonomy_species_id,
        grin_id_cur = tax_i$current_taxonomy_species_id,
        ISO3_n = geo_i
      )
    } else {
      dist_prop$GRIN[which(dist_prop$taxon == unique_taxa[[i]])] <- 0
      accepted_grin_tax_list[[i]] <- data.frame(taxa = unique_taxa[[i]],
                                                grin_id = NA,
                                                grin_id_cur = NA,
                                                grin_taxon=NA)
      native_countries_list[[i]] <- data.frame(
        taxa = unique_taxa[[i]],
        grin_id = NA,
        grin_id_cur = NA,
        ISO3_n = NA
      )
    }
    
  };rm(i)
  
  message("OBTAINING species matching with GRIN TAXONOMY")
  accepted_grin_tax_list <- do.call(rbind, accepted_grin_tax_list)
  
  
  dist_prop$GRIN_ID <- NA
  dist_prop$GRIN_taxon <- NA
  for(i in 1:nrow(dist_prop)){
    #i <- 1
    dist_prop$GRIN_ID[[i]] <-accepted_grin_tax_list$grin_id_cur[accepted_grin_tax_list$taxa %in% dist_prop$taxon[[i]]]
    dist_prop$GRIN_taxon[[i]] <-accepted_grin_tax_list$grin_taxon[accepted_grin_tax_list$taxa %in% dist_prop$taxon[[i]]]
  };rm(i)
#  dist_prop$GRIN_ID <-    accepted_grin_tax_list$grin_id_cur
  ###############################################################################
  message("OBTAINING NATIVE COUNT RIES (ISO3) FORSPECIES MATCHING WITH GRIN TAXONOMY")
  native_countries_list <- do.call(rbind, native_countries_list)
  
  ###############################################################################
  #species in grin
  message("OBTAINING species matching with World Flora")
  #
  # if(is.data.frame(df_wild)){
  # WILD_SPP <- accepted_grin_tax_list[accepted_grin_tax_list$taxa %in% (unlist(lapply(strsplit(
  #   row.names(df_wild), "_"
  # ), `[[`, 1))), ]
  # } else {
  #   WILD_SPP <- NA
  # }
  #
  # if(is.data.frame(df_landrace)){
  # LAND_SPP <- accepted_grin_tax_list[accepted_grin_tax_list$taxa %in% (unlist(lapply(strsplit(
  #   row.names(df_landrace), "_"
  # ), `[[`, 1))), ]
  # } else {
  #   LAND_SPP <- NA
  # }
  #species in WorlFlora
  
  dist_prop$WOF_taxa <- NA 
  dist_prop$WOF <- NA 
  dist_prop$WOF_taxa_id <- NA 
  
  #WorldFlora::WFO.remember(WFO.file = "D:/OneDrive - CGIAR/GERMPLASM_INDEX/COMPRESSED_FILES/WFO/classification.csv")
  for (i in 1:nrow(dist_prop)){
    #message(i)
    #i <- 1
    #dist_prop$taxon[[i]]
    #x1 <- WorldFlora::WFO.match(dist_prop$taxon[[i]],WFO.data = WFO.data,verbose = F)
    #x1 <- WorldFlora::WFO.match("Phaseolus grayanus",WFO.data = WFO.data,verbose = F)
    #x1 <- WorldFlora::WFO.one(WorldFlora::WFO.match.fuzzyjoin(dist_prop$taxon[[i]],
    x1 <- WorldFlora::WFO.one(WorldFlora::WFO.match(dist_prop$taxon[[i]],
                        WFO.data=WFO.data,verbose = F
                        #fuzzydist.max=3
                        ),verbose = F)

    #x1 <- WFO.one(x1,verbose = F,priority = "Accepted")#x1[x1$Subseq==1,]
    #WorldFlora::WFO.match("Phaseolus nanus",WFO.data = WFO.data)
    dist_prop$WOF_taxa[[i]] <- x1$scientificName
    dist_prop$WOF_taxa_id[[i]] <- x1$taxonID
    if(is.na(x1$scientificName)){
      dist_prop$WOF[[i]] <- 0
    } else {
      dist_prop$WOF[[i]] <- 1
    }
    rm(x1)
  };rm(i)
  ################################################################################
  message("Subsetting to native countries")
  #initializing empty rows for native countries and collected n_countries
  
  dist_prop$native_countries <- NA
  dist_prop$total_native_countries <- NA
  dist_prop$non_native_countries <- NA
  dist_prop$native_countries_R <- NA
  dist_prop$non_native_countries_R <- NA
  dist_prop$na_countries_R <- NA
  for (i in 1:nrow(dist_prop)) {
    #i <- 2
    x_i <- dist_prop[i, ]
    nat_i <- native_countries_list$ISO3_n[which(native_countries_list$taxa ==
                                                  x_i$taxon)]
    if (x_i$status == "wild") {
      x_count <- df_wild[which(row.names(df_wild) == x_i$Taxa), ]
      x_count_nat <- x_count[, colnames(x_count) %in% nat_i]
      x_count_nonnat <- x_count[, !colnames(x_count) %in% nat_i]
      x_count_na <- x_count_nonnat
      x_count_nonnat <- x_count_nonnat[which(colnames(x_count_nonnat) !=
                                               "NA1")]
      x_count_na <- x_count_na[which(colnames(x_count_na) == "NA1")]
      dist_prop$native_countries[[i]] <- sum(x_count_nat > 0)
      dist_prop$total_native_countries[[i]] <- length(nat_i)
      dist_prop$non_native_countries[[i]] <- sum(x_count_nonnat > 0)
      dist_prop$native_countries_R[[i]] <- sum(x_count_nat)
      dist_prop$non_native_countries_R[[i]] <- sum(x_count_nonnat)
      dist_prop$na_countries_R[[i]] <- sum(x_count_na)
    } else if (x_i$status == "landrace") {
      x_count <- df_landrace[which(row.names(df_landrace) == x_i$Taxa), ]
      x_count_nat <- x_count[, colnames(x_count) %in% nat_i]
      x_count_nonnat <- x_count[, !colnames(x_count) %in% nat_i]
      x_count_na <- x_count_nonnat
      x_count_nonnat <- x_count_nonnat[which(colnames(x_count_nonnat) !=
                                               "NA1")]
      x_count_na <- x_count_na[which(colnames(x_count_na) == "NA1")]
      dist_prop$native_countries[[i]] <- sum(x_count_nat > 0)
      dist_prop$total_native_countries[[i]] <- length(nat_i)
      dist_prop$non_native_countries[[i]] <- sum(x_count_nonnat > 0)
      dist_prop$native_countries_R[[i]] <- sum(x_count_nat)
      dist_prop$non_native_countries_R[[i]] <- sum(x_count_nonnat)
      dist_prop$na_countries_R[[i]] <- sum(x_count_na)
    } else if (x_i$status == "hybrid") {
      x_count <- df_hybrid[which(row.names(df_hybrid) == x_i$Taxa), ]
      x_count_nat <- x_count[, colnames(x_count) %in% nat_i]
      x_count_nonnat <- x_count[, !colnames(x_count) %in% nat_i]
      x_count_na <- x_count_nonnat
      x_count_nonnat <- x_count_nonnat[which(colnames(x_count_nonnat) !=
                                               "NA1")]
      x_count_na <- x_count_na[which(colnames(x_count_na) == "NA1")]
      dist_prop$native_countries[[i]] <- sum(x_count_nat > 0)
      dist_prop$total_native_countries[[i]] <- length(nat_i)
      dist_prop$non_native_countries[[i]] <- sum(x_count_nonnat > 0)
      dist_prop$native_countries_R[[i]] <- sum(x_count_nat)
      dist_prop$non_native_countries_R[[i]] <- sum(x_count_nonnat)
      dist_prop$na_countries_R[[i]] <- sum(x_count_na)
    }
  }
  
  dist_prop$total_records <- NA
  
  dist_prop$total_records <-   dist_prop$native_countries_R +
  dist_prop$non_native_countries_R + dist_prop$na_countries_R
  #ecogeography indexes
  message("Calculating ecogeographic passport based indexes")
  dist_prop$nCounnat_countotal <- dist_prop$native_countries / dist_prop$total_native_countries
  dist_prop$nCounnonnat_countotal <- dist_prop$native_countries / dist_prop$ncountries_total_collection
  dist_prop$nat_nonnat <- dist_prop$native_countries / dist_prop$non_native_countries
  dist_prop$nat_nonnat[which(dist_prop$nat_nonnat == Inf)] <- NA
  
  
  write.csv(
    native_countries_list,
    paste0(outdir, "/", collection_name, "/", collection_name, "_native_iso3_new_1.csv"),
    row.names = F,
    na = ""
  )
  message(paste0("Saved: ", outdir, "/", collection_name, "/", collection_name, "_native_iso3_1.csv"))
  ################################################################################
  ################################################################################
  
  # df_wild$species <- NA
  # df_wild$species <-  dist_prop$taxon[which(dist_prop$status=="wild")]
  ################################################################################
  #adding IUCN
  message("Obtaining IUCN status (preparing data)")
  
  stat <- lapply(1:nrow(dist_prop), function(i) {
    #i <- 53
    # message(i)
    x <-  unlist(strsplit(dist_prop$taxon[[i]], " "))
    x_tax_df <- data.frame(GENUS = NA,
                           SPECIES = NA,
                           SUBTAXA = NA)
    if (length(x) > 2) {
      x_tax_df$GENUS <-  x[1]
      
      ####
      #DETECTING parenthesis
      par_x <- stringr::str_detect(x, "\\(")
      if (par_x[[2]] == TRUE) {
        #x_tax_df$SPECIES <-  x[2]
        x_tax_df$SPECIES <-  stringr::str_c(x[2:length(x)], collapse = " ")
        x_tax_df$SUBTAXA <-  ""
      } else if (x[2] == "x") {
        x_tax_df$SPECIES <-  stringr::str_c(x[2:length(x)], collapse = " ")
        x_tax_df$SUBTAXA <-  ""
      } else {
        x_tax_df$SPECIES <- x[2]
        x_tax_df$SUBTAXA <-  stringr::str_c(x[3:length(x)], collapse = " ")
      }
    } else {
      x_tax_df$GENUS <-  x[1]
      x_tax_df$SPECIES <-  x[2]
      x_tax_df$SUBTAXA <-  ""
    }
    return(x_tax_df)
  })
  
  stat <- do.call(rbind, stat)
  
  stat$IUCN <- NA
  
  
  ################################################################################
  message("Obtaining IUCN status")
  
  headers = c(accept = "application/json", Authorization = API)
  
  for (i in 1:nrow(stat)) {
    #i <- 1
    #message(i)
    #adding parameters
    params = list(
      genus_name = stat$GENUS[[i]],
      species_name = stat$SPECIES[[i]],
      infra_name = stat$SUBTAXA[[i]]
      
    )
    #obtaining data associated to the scientific name
    res <- httr::GET(
      url = "https://api.iucnredlist.org/api/v4/taxa/scientific_name",
      accept_json(),
      httr::add_headers(.headers = headers),
      query = params
    )
    #obtaining the taxon information
    x <- httr::content(res)
    #obtaining data (if x is an error the length is 1)
    if (length(x) > 1) {
      if (length(x$assessments) > 0) {
        #obtaining the most recent assesment to obtain IUCN status
        res2 <- httr::GET(
          url = paste0(
            "https://api.iucnredlist.org/api/v4/assessment/",
            x$assessments[[length(x$assessments)]]$assessment_id
          ),
          httr::add_headers(.headers = headers)
        )
        
        x2 <- httr::content(res2)
        #x2
        stat$IUCN[[i]] <- x2$red_list_category$description$en
      } else {
        stat$IUCN[[i]] <- NA
      }
    } else {
      stat$IUCN[[i]] <- NA
    }
    Sys.sleep(0.01)
  };rm(i)  # params = list(
  #   genus_name = "Phaseolus",
  #   species_name = "vulgaris",
  #   infra_name = ""
  #
  # )
  
  
  #adding info to summary data
  dist_prop$IUCN <- NA
  
  for (i in 1:nrow(dist_prop)) {
    dist_prop$IUCN[[i]] <- stat$IUCN[[i]]
    
  }
  rm(i)
  
  #changing order
  dist_prop$Taxa <- dist_prop$taxon
  dist_prop$taxon <- NULL
  
 dist_prop <- 
  dist_prop[,c(1,16,2:15,17,18,19,21,22,20,33,23:32)]
  write.csv(
    dist_prop,
    paste0(outdir, "/", collection_name,"/", collection_name,  "_summary_table.csv"),
    row.names = F,
    na = ""
  )
  message(paste0("Saved: ", outdir, "/", collection_name, "/", collection_name, "_summary_table.csv"))
  
  
  ################################################################################
  
  message("CALCULATING GINI-SIMPSON FOR COUNTRIES")
  
  #taxonomy index is 1- Gini (inequality) * 1- alpha (diversity) * Piellou (eveness)
#  "
#if the index is high then the collection has a high diversity, good balance of records,
#and good eveness (similar number of records)
#"
  # #taxonomy index all
  # I1_all <- (1 - ineq::ineq(dist_prop$ncountries, type = "Gini")) *
  #   abdiv::simpson(dist_prop$ncountries)# *
  # #abdiv::pielou_e(dist_prop$ncountries)
  # #taxonomy index crop wild
  # I1_W <- (1 - ineq::ineq(dist_prop$ncountries[which(dist_prop$status == "wild")], type =
  #                           "Gini")) *
  #   abdiv::simpson(dist_prop$ncountries[which(dist_prop$status == "wild")]) #*
  # #abdiv::pielou_e(dist_prop$ncountries[which(dist_prop$status=="wild")])
  # #taxonomy index landrace
  # I1_L  <- (1 - ineq::ineq(dist_prop$ncountries[which(dist_prop$status == "landrace")], type =
  #                            "Gini")) *
  #   abdiv::simpson(dist_prop$ncountries[which(dist_prop$status == "landrace")]) #*
  # #abdiv::pielou_e(dist_prop$ncountries[which(dist_prop$status=="landrace")])
  #
  # write.csv(dist_prop,
  #           paste0(outdir, "/", collection_name, "_diss_country.csv"),row.names = F,na = "")
  # message(paste0("Saved:",outdir, "/", collection_name, "_diss_country.csv"))
  
  ##############################################################################
  #proportion of species with one record wild
  PROPSP1R_W <- nrow(dist_prop[which(dist_prop$ncounts_Country == 1 &
                                       dist_prop$status == "wild"), ]) / nrow(dist_prop) #+#
  #proportion of species with one record landraces
  PROPSP1R_L <- nrow(dist_prop[which(dist_prop$ncounts_Country == 1 &
                                       dist_prop$status == "landrace"), ]) / nrow(dist_prop) #+#
  #proportion of species with one record hybrids
  PROPSP1R_H <- nrow(dist_prop[which(dist_prop$ncounts_Country == 1 &
                                       dist_prop$status == "hybrid"), ]) / nrow(dist_prop) #+#
  ##############################################################################
  #proportion of species in wild and landrace
  PROPSP_W <- nrow(dist_prop[which(dist_prop$status == "wild"), ]) / nrow(dist_prop) #+#
  #proportion of species with one record landraces
  PROPSP_L <- nrow(dist_prop[which(dist_prop$status == "landrace"), ]) /
    nrow(dist_prop) #+#
  #proportion of species with one record hybrids
  PROPSP_H <- nrow(dist_prop[which(dist_prop$status == "hybrid"), ]) / nrow(dist_prop) #+#
  ##############################################################################
  #proportion of records in wild and landrace
  PROP_R_W <- sum(dist_prop$ncounts_Country[which(dist_prop$status == "wild")]) /
    sum(dist_prop$ncounts_Country) #+#
  #proportion of species with one record landraces
  PROP_R_L <- sum(dist_prop$ncounts_Country[which(dist_prop$status == "landrace")]) /
    sum(dist_prop$ncounts_Country) #+#
  #proportion of species with one record hybrid
  PROP_R_H <- sum(dist_prop$ncounts_Country[which(dist_prop$status == "hybrid")]) /
    sum(dist_prop$ncounts_Country) #+#
  ##############################################################################
  #proportion of IUCN (wild)
  PROP_W_IUCN <-
    nrow(dist_prop[dist_prop$IUCN %in%
                     c("Vulnerable",
                       "Endangered",
                       "Critically Endangered",
                       "Extinct in the Wild") & dist_prop$status == "wild", ]) / nrow(dist_prop)
  #proportion of IUCN (landraces)
  PROP_L_IUCN <-
    nrow(dist_prop[dist_prop$IUCN %in%
                     c("Vulnerable",
                       "Endangered",
                       "Critically Endangered",
                       "Extinct in the Wild") &
                     dist_prop$status == "landrace", ]) / nrow(dist_prop)
  #proportion of IUCN (hybrids)
  PROP_H_IUCN <-
    nrow(dist_prop[dist_prop$IUCN %in%
                     c("Vulnerable",
                       "Endangered",
                       "Critically Endangered",
                       "Extinct in the Wild") &
                     dist_prop$status == "hybrid", ]) / nrow(dist_prop)
  ##############################################################################
  #GINI SIMPSON
  #WILD
  if (length(dist_prop$ncounts_Country[which(dist_prop$status == "wild")]) >
      0) {
    GINI_SIMPSON_W <- abdiv::simpson(dist_prop$ncounts_Country[which(dist_prop$status ==
                                                                       "wild")]) #*
  } else {
    GINI_SIMPSON_W <- 0
  }
  #LANDRACE
  if (length(dist_prop$ncounts_Country[which(dist_prop$status == "landrace")]) >
      0) {
    GINI_SIMPSON_L <- abdiv::simpson(dist_prop$ncounts_Country[which(dist_prop$status ==
                                                                       "landrace")]) #*
  } else {
    GINI_SIMPSON_L <- 0
  }
  #HYBRID
  if (length(dist_prop$ncounts_Country[which(dist_prop$status == "hybrid")]) >
      0) {
    GINI_SIMPSON_H <- abdiv::simpson(dist_prop$ncounts_Country[which(dist_prop$status ==
                                                                       "hybrid")]) #*
  } else {
    GINI_SIMPSON_H <- 0
  }
  ##############################################################################
  #PROP RECORDS IN w AND L
  #WILD
  PROP_R_W <- sum(dist_prop$ncounts_Country[which(dist_prop$status == "wild")]) /
    sum(dist_prop$ncounts_Country)
  #LANDRACE
  PROP_L_W <- sum(dist_prop$ncounts_Country[which(dist_prop$status == "landrace")]) /
    sum(dist_prop$ncounts_Country)
  #HYBRIDS
  PROP_L_H <- sum(dist_prop$ncounts_Country[which(dist_prop$status == "hybrid")]) /
    sum(dist_prop$ncounts_Country)
  ##############################################################################
  #N SPP IN taxonomy
  #wild
  SGRIN_W <- sum(dist_prop$GRIN[which(dist_prop$status == "wild")]) / nrow(dist_prop)
  #landrace
  SGRIN_L <- sum(dist_prop$GRIN[which(dist_prop$status == "landrace")]) /
    nrow(dist_prop)
  #hybrid
  SGRIN_H <- sum(dist_prop$GRIN[which(dist_prop$status == "hybrid")]) /
    nrow(dist_prop)
  ##############################################################################
  #N SPP IN taxonomy
  #wild
  SWOF_W <- sum(dist_prop$WOF[which(dist_prop$status == "wild")]) / nrow(dist_prop)
  #landrace
  SWOF_L <- sum(dist_prop$WOF[which(dist_prop$status == "landrace")]) /
    nrow(dist_prop)
  #hybrid
  SWOF_H <- sum(dist_prop$WOF[which(dist_prop$status == "hybrid")]) /
    nrow(dist_prop)
  ##############################################################################
  #wild
  
  x <- dist_prop$ncounts_Country[which(dist_prop$status == "wild")]
  if (length(x) > 0) {
    ATK_W <- 1 - DescTools::Atkinson(x)
  } else {
    ATK_W <- NA
  }
  rm(x)
  
  #landrace
  x <- dist_prop$ncounts_Country[which(dist_prop$status == "landrace")]
  if (length(x) > 0) {
    ATK_L <- 1 - DescTools::Atkinson(x)
  } else {
    ATK_L <- NA
  }
  rm(x)
  
  #hybrids
  x <- dist_prop$ncounts_Country[which(dist_prop$status == "hybrid")]
  if (length(x) > 0) {
    ATK_H <- 1 - DescTools::Atkinson(x)
  } else {
    ATK_H <- NA
  }
  rm(x)
  ##############################################################################
  # #wild
  # IT1_W  <- (1 - PROPSP1R_W) * GINI_SIMPSON_W * SGRIN_W
  # #landrace
  # IT1_L  <-   (1 - PROPSP1R_L) * GINI_SIMPSON_L * SGRIN_L
  # #hybrid
  # IT1_H  <-   (1 - PROPSP1R_H) * GINI_SIMPSON_H * SGRIN_H
  # 
  ##############################################################################
  #IT1_FINAL_MEAN <- (IT1_W + IT1_L)/2
  #IT_FINAL_MIN <- min(IT_W,IT_L)
  
  ##############################################################################
  # IT2_W  <- ATK_W * GINI_SIMPSON_W * SGRIN_W
  # IT2_L  <-   ATK_L * GINI_SIMPSON_L * SGRIN_L
  # #IT2_FINAL_MEAN <- (IT2_W + IT2_L)/2
  ##############################################################################
  # #parameter D
  # D_W <- 0.5
  # D_L <- 0.5
  # 
  # if (PROP_R_W == D_W) {
  #   PROPSTATUS_W <- 0
  # } else {
  #   PROPSTATUS_W <- 1 - PROP_R_W
  # }
  # 
  # 
  # if (PROP_R_L == D_L) {
  #   PROPSTATUS_L <- 0
  # } else {
  #   PROPSTATUS_L <- 1 - PROP_R_L
  # }
  # 
  # IT3_W  <- ((((1 - PROPSP1R_W) * GINI_SIMPSON_W) + SGRIN_W) - PROPSTATUS_W) /
  #   2
  # IT3_L  <- ((((1 - PROPSP1R_L) * GINI_SIMPSON_L) + SGRIN_L) - PROPSTATUS_L) /
  #   2
  # if (IT3_L < 0) {
  #   IT3_L <- 0
  # }
  # 
  # if (IT3_W < 0) {
  #   IT3_W <- 0
  # }
  # 
  # 
  # IT3_FINAL_MEAN <- (IT3_W + IT3_L) / 2
  ##############################################################################
  #NEW SCORE 
  IT3_W  <- (GINI_SIMPSON_W +(PROPSP1R_W*ATK_W) + SGRIN_W)/3
  IT3_L  <- (GINI_SIMPSON_L +(PROPSP1R_L*ATK_L) + SGRIN_L)/3
  IT3_H  <- (GINI_SIMPSON_H +(PROPSP1R_H*ATK_H) + SGRIN_H)/3
  
  IT3_FINAL_MEAN <- (IT3_W + IT3_L) / 2
  ##############################################################################
  TAX_DF <- data.frame(matrix(nrow = 3, ncol = 14))
  colnames(TAX_DF) <- c(
    "type",
    "nTaxa",
    "nRecords",
    "nRecords_NA",
    "Records/Taxa",
    "propTaxa_GRIN",
    "propTaxa_WOF",
    "propRecords",
    "propTaxa_collection",
    "propSingleton",
    "Gini_Simpson",
    #"NA",#"PROPStatus_index3",
    "propTaxa_IUCN",
    "1_Atkinson",
    #"Index_1",
    #"Index_2",
    "Composition_index"
  )
  
  #COLLECTION TYPES
  TAX_DF[, 1] <- c("wild", "landraces", "hybrid")
  
  
  #WILD
  TAX_DF[1, 2] <- nrow(dist_prop[which(dist_prop$status == "wild"), ])
  TAX_DF[1, 3] <-  sum(dist_prop$ncounts_Country[which(dist_prop$status ==
                                                         "wild")])
  TAX_DF[1, 4] <-  sum(dist_prop$ncounts_NA[which(dist_prop$status ==
                                                                   "wild")])
  TAX_DF[1, 5] <-  TAX_DF[1, 3]/TAX_DF[1, 2] 
  TAX_DF[1, 6] <-  SGRIN_W
  TAX_DF[1, 7] <-  SWOF_W
  TAX_DF[1, 8] <-  PROP_R_W
  TAX_DF[1, 9] <-  PROPSP_W
  TAX_DF[1, 10] <-  PROPSP1R_W
  TAX_DF[1, 11] <-  GINI_SIMPSON_W
  #TAX_DF[1, 11] <- #NA#PROPSTATUS_W
  TAX_DF[1, 12] <- PROP_W_IUCN
  TAX_DF[1, 13] <- ATK_W
  #TAX_DF[1, 13] <- IT1_W
  #TAX_DF[1, 14] <-  IT2_W
  TAX_DF[1, 14] <-  IT3_W
  
  #LANDRACE
  TAX_DF[2, 2] <- nrow(dist_prop[which(dist_prop$status == "landrace"), ])
  TAX_DF[2, 3] <-  sum(dist_prop$ncounts_Country[which(dist_prop$status ==
                                                         "landrace")])
  TAX_DF[2, 4] <-  sum(dist_prop$ncounts_NA[which(dist_prop$status ==
                                                                   "landrace")])
  TAX_DF[2, 5] <-  TAX_DF[2, 3]/TAX_DF[2, 2] 
  TAX_DF[2, 6] <-  SGRIN_L
  TAX_DF[2, 7] <-  SWOF_L
  TAX_DF[2, 8] <-  PROP_R_L
  TAX_DF[2, 9] <-  PROPSP_L
  TAX_DF[2, 10] <-  PROPSP1R_L
  TAX_DF[2, 11] <-  GINI_SIMPSON_L
  #TAX_DF[2, 11] <- NA#PROPSTATUS_L
  TAX_DF[2, 12] <-  PROP_L_IUCN
  TAX_DF[2, 13] <-  ATK_L
  #TAX_DF[2, 13] <- IT1_L
  #TAX_DF[2, 14] <-  IT2_L
  TAX_DF[2, 14] <-  IT3_L
  
  
  # HYBRIDS
  TAX_DF[3, 2] <- nrow(dist_prop[which(dist_prop$status == "hybrid"), ])
  TAX_DF[3, 3] <-  sum(dist_prop$ncounts_Country[which(dist_prop$status ==
                                                         "hybrid")])
  TAX_DF[3, 4] <-  sum(dist_prop$ncounts_NA[which(dist_prop$status ==
                                                                   "hybrid")])
  TAX_DF[3, 5] <-  TAX_DF[3, 3]/TAX_DF[3, 2] 
  TAX_DF[3, 6] <-  SGRIN_H
  TAX_DF[3, 7] <-  SWOF_H
  TAX_DF[3, 8] <-  PROP_R_H
  TAX_DF[3, 9] <-  PROPSP_H
  TAX_DF[3, 10] <-  PROPSP1R_H
  TAX_DF[3, 11] <-  GINI_SIMPSON_H
  #TAX_DF[3, 11] <- NA
  TAX_DF[3, 12] <-  NA
  TAX_DF[3, 13] <-  ATK_H
  #TAX_DF[3, 13] <- NA
  #TAX_DF[3, 14] <-  NA
  TAX_DF[3, 14] <-  NA
  
  
  
  write.csv(
    TAX_DF,
    paste0(outdir, "/", collection_name, "/", collection_name, "_IT_table_1.csv"),
    row.names = F,
    na = ""
  )
  message(paste0("Saved: ", outdir, "/", collection_name, "/", collection_name, "_IT_table_1.csv"))
  
  
  message("DONE!")
  #return(x_Gini)
  return(TAX_DF)
}

################################################################################
################################################################################
dir <- "D:/OneDrive - CGIAR/GERMPLASM_INDEX"
#dir <- "D:/ONEDRIVE/cgiar/OneDrive - CGIAR/GERMPLASM_INDEX"
################################################################################
#outdir
outdir <- paste0(dir, "/BEANS/RESULTS")
################################################################################

inDir <- "D:/OneDrive - CGIAR/GERMPLASM_INDEX/COMPRESSED_FILES/A"
#inDir <- "D:/ONEDRIVE/cgiar/OneDrive - CGIAR/GERMPLASM_INDEX/COMPRESSED_FILES/A"
tax_table <- "taxonomy_species.txt"
geo_table <- "geography.txt"
taxonomy_geography_map_table <- "taxonomy_geography_map.txt"
#the logic is use tax table to obtain current grin tax id, then join with taxonomy_geography_map_table to
#obtain the geography_id, first filter using the geography_status_code
#finally use the geo_table to get the ISO3 and is_valid is_valid
################################################################################
#reading taxonomy table 
if(!exists("tax")){
  tax <- data.table::fread(paste0(inDir, "/", tax_table))  
  #obtaining taxon status (wild or landrace)
  stat <- strsplit(tax$name, " ")
  #status
  tax$GENUS <- trimws(unlist(lapply(stat, `[[`, 1)))
  tax <- as.data.frame(tax)
}

################################################################################
#geography data 
if(!exists("geo")){
  geo <- as.data.frame(data.table::fread(paste0(inDir, "/", geo_table)))
  }

#table with geography and taxonomy data
if(!exists("tax_geo")){
tax_geo <- as.data.frame(data.table::fread(paste0(inDir, "/", taxonomy_geography_map_table)))
}
################################################################################
#read PDCI data
#passport_data_original <- read.csv(paste0(dir, "/", "CIAT Data/ciat_pdci.csv"), na.strings = NA)
if(!exists("passport_data_original")){
  passport_data_original <- read.csv(paste0(dir, "/", "CIAT Data/genesys-accessions-COL003.csv"),
                                     na.strings = NA)  
}

#subsetting to beans
# passport_data_orig <- passport_data_orig[which(passport_data_orig$GENUS ==
################################################################################
#API to obtain IUCN status
API <- "mtfcmg5AttZ2kaJSWdsq9MkEjfhW41gjeSvm"
################################################################################
#Loading WorldFlora Online
if(!exists("WFO.data")){
WorldFlora::WFO.remember("D:/OneDrive - CGIAR/GERMPLASM_INDEX/COMPRESSED_FILES/WFO_Backbone.zip")
}
  #WorldFlora::WFO.remember(WFO.file = "D:/OneDrive - CGIAR/GERMPLASM_INDEX/COMPRESSED_FILES/WFO/classification.csv")

################################################################################
#CIAT codes= "beans"   "forages" "cassava"
passport_data_orig <- passport_data_original[which(passport_data_original$CROPCODE ==
                                                     "cassava"), ]
collection_name = "cassava"
x1 <- index_function(
  passport_data_orig = passport_data_orig,
  tax = tax,
  geo = geo,
  tax_geo = tax_geo,
  outdir = outdir,
  API = API,
  collection_name = collection_name
)

################################################################################
#CIAT codes= "beans"   "forages" "cassava"
passport_data_orig <- passport_data_original[which(passport_data_original$CROPCODE ==
                                                     "beans"), ]
collection_name = "beans"
x2 <- index_function(
  passport_data_orig = passport_data_orig,
  tax = tax,
  geo = geo,
  tax_geo = tax_geo,
  outdir = outdir,
  API = API,
  collection_name = collection_name
)
################################################################################
#CIAT codes= "beans"   "forages" "cassava"
passport_data_orig <- passport_data_original[which(passport_data_original$CROPCODE ==
                                                     "forages"), ]
collection_name = "forages"
x3 <- index_function(
  passport_data_orig = passport_data_orig,
  tax = tax,
  geo = geo,
  tax_geo = tax_geo,
  outdir = outdir,
  API = API,
  collection_name = collection_name
)


X_TOTAL <- rbind(x1, x2)
X_TOTAL <- rbind(X_TOTAL, x3)
X_TOTAL$COLLECTION <- c(rep("cassava",3), rep("beans",3),rep("forages",3))
X_TOTAL <- X_TOTAL[,c(17,1:16)]
#save.image("~/AAA.RData")
#load("~/AAA.RData")
write.csv(X_TOTAL,
          paste0(outdir, "/", "IT_CIAT.csv"),
          row.names = F,
          na = "")
