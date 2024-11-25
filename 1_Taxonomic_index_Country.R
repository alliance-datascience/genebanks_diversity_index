  ################################################################################
  library(countrycode)
  library(ggplot2)
  library(factoextra)
  library(glue)
  library(fossil)
  library(plsdepot)
  library(ggrepel)
  library(MetBrewer)
  ################################################################################
  # Define a function for decimal scaling 
  decimal_scale <- function(x) { 
    # Find the maximum absolute value of x 
    max_abs <- max(abs(x),na.rm = T) 
    # Find the smallest power of 10 that is equal to or larger than max_abs 
    power <- ceiling(log10(max_abs))
    # Divide x by 10^power 
    x <- x/ (10^power)
    return(x)
  }
  
  ################################################################################
  sigmoid_scale <- function(x){
    x <- 1/(1+exp((-(x - mean(x,na.rm = T)))/sd(x,na.rm = T)))
    return(x)
  }
  ################################################################################
  
  composition_country <- function(outdir, collection_name,plot_type) {
    
    if (is.null(plot_type) | !plot_type %in% c("Margalef","Simpson")) {
      stop("Please use a valid option: Margalef,Simpson")
    }
    
    message("Reading  previous results")
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
    ##############################################################################
    #reading summary table
    dist_prop <-
      read.csv(paste0(
        outdir,
        "/",
        collection_name,
        "/",
        collection_name,
        "_summary_table.csv"
      ))
    ##############################################################################
    #reading collection composition index
    TAX_DF <- read.csv(paste0(
      outdir,
      "/",
      collection_name,
      "/",
      collection_name,
      "_IT_table_1.csv"
    ))
    ##############################################################################
    #subsetting data in CWR, landraces and hybrids
    distw <- dist_prop[which(dist_prop$status == "wild"), ]
    distl <- dist_prop[which(dist_prop$status == "landrace"), ]
    disth <- dist_prop[which(dist_prop$status == "hybrid"), ]
    ################################################################################
    simp_list <- list()
    
    for (i in 1:ncol(df_wild)) {
      #message(paste(i,"/",colnames(df_wild)[i]))
      #i <- 89
      #i <- 3
      ############################################################################
      #Obtaining grin proportion taxa (wild)
      if(nrow(distw)>0){
        distw_copy <- NA
        distw$coun <- NA
        distw$coun <- df_wild[, i]
        f1_w <- length(distw$coun[which(distw$coun==1)])
        f2_w <- length(distw$coun[which(distw$coun==2)])
        n_w <- sum(df_wild[, i])
        distw$coun[which(distw$coun > 0)] <- 1
        distw$coun <- distw$coun * distw$GRIN
        G_W_i <- sum(distw$coun) / nrow(dist_prop)
      } else {
        G_W_i <- NA
        f1_w <- NA
        f2_w <- NA
        n_w <- NA
      }    
      #Obtaining grin proportion taxa (landraces)
      if(nrow(distl)>0){
        distl$coun <- NA
        distl$coun <- df_landrace[, i]
        f1_l <- length(distl$coun[which(distl$coun==1)])
        f2_l <- length(distl$coun[which(distl$coun==2)])
        n_l <- sum(df_landrace[, i])
        distl$coun[which(distl$coun > 0)] <- 1
        distl$coun <- distl$coun * distl$GRIN
        G_L_i <- sum(distl$coun) / nrow(dist_prop) 
      } else {
        G_L_i <- NA
        f1_l <- NA
        f2_l <- NA
        n_l <- NA
      }
      #Obtaining grin proportion taxa (hybrids)
      if(nrow(disth)>0){
        disth$coun <- NA
        disth$coun <- df_hybrid[, i]
        f1_h <- length(disth$coun[which(disth$coun==1)])
        f2_h <- length(disth$coun[which(disth$coun==2)])
        n_h <- sum(df_hybrid[, i])
        disth$coun[which(disth$coun > 0)] <- 1
        disth$coun <- disth$coun * disth$GRIN
        G_H_i <- sum(disth$coun) / nrow(dist_prop)
        disth$coun <- NULL
      } else {
        G_H_i <- NA
        f1_h <- NA
        f2_h <- NA
        n_h <- NA
      }
      
      #Obtaining coverage (Chao 2012)
      cov_w <-  1-(f1_w/n_w)* (((n_w-1)*f1_w)/ (((n_w-1)*f1_w)+(2*f2_w)))
      cov_l <-  1-(f1_l/n_l)* (((n_l-1)*f1_l)/ (((n_l-1)*f1_l)+(2*f2_l)))
      cov_h <-  1-(f1_h/n_h)* (((n_h-1)*f1_h)/ (((n_h-1)*f1_h)+(2*f2_h)))
      ############################################################################
      #summing species
      ############################################################################
      #"using wild data" (Gini Simpson, 1- Atkinson, singletons)
      if(nrow(distw)>0){
        SIMP_WILD = abdiv::dominance(df_wild[, i])
        ATK_WILD = 1-DescTools::Atkinson(df_wild[, i],parameter = 0.5)
        SP1_WILD = length(df_wild[, i][which(df_wild[, i] == 1)])
        PROP_1_SP1_WILD = 1 - (length(df_wild[, i][which(df_wild[, i] ==
                                                           1)]) /
                                 (nrow(dist_prop)))
        SPP_WILD = length(df_wild[, i][which(df_wild[, i] > 0)])
        PROP_SPP_WILD = length(df_wild[, i][which(df_wild[, i] > 0)]) /
          (nrow(dist_prop))
        REC_WILD = sum(df_wild[, i][which(df_wild[, i] > 0)],na.rm = T)
        PROP_R_WILD = sum(df_wild[, i][which(df_wild[, i] > 0)]) / sum(df_wild, df_landrace, df_hybrid,na.rm = T)
        ENS_GINI_WILD = 1 / abdiv::dominance(df_wild[, i])
        MARGALEF_WILD = abdiv::margalef(df_wild[,i])
      } else {
        SIMP_WILD = NA
        ATK_WILD = NA
        SP1_WILD = NA
        PROP_1_SP1_WILD =NA
        SPP_WILD = NA
        PROP_SPP_WILD = NA
        REC_WILD = NA
        PROP_R_WILD = NA
        ENS_GINI_WILD = NA
        MARGALEF_WILD = NA
      }
      ############################################################################
      #"using landraces data" (Gini Simpson, 1- Atkinson, singletons)
      if(nrow(distl)>0){
        SIMP_LAND = abdiv::dominance(df_landrace[, i])
        ATK_LAND = 1-DescTools::Atkinson(df_landrace[, i],parameter = 0.5)
        SP1_LAND = length(df_landrace[, i][which(df_landrace[, i] == 1)])
        PROP_1_SP1_LAND = 1 - (length(df_landrace[, i][which(df_landrace[, i] ==
                                                               1)]) /
                                 (nrow(dist_prop)))
        SPP_LAND = length(df_landrace[, i][which(df_landrace[, i] > 0)])
        PROP_SPP_LAND = length(df_landrace[, i][which(df_landrace[, i] > 0)]) /
          (nrow(dist_prop))
        REC_LAND = sum(df_landrace[, i][which(df_landrace[, i] > 0)],na.rm = T)
        PROP_R_LAND = sum(df_landrace[, i][which(df_landrace[, i] > 0)]) / sum(df_wild, df_landrace, df_hybrid,na.rm = T)
        ENS_GINI_LAND = 1 / abdiv::dominance(df_landrace[, i])
        MARGALEF_LAND = abdiv::margalef(df_landrace[,i])
        
      } else {
        SIMP_LAND = NA
        ATK_LAND = NA
        SP1_LAND = NA
        PROP_1_SP1_LAND =NA
        SPP_LAND = NA
        PROP_SPP_LAND = NA
        REC_LAND = NA
        PROP_R_LAND = NA
        ENS_GINI_LAND = NA
        MARGALEF_LAND = NA
      }
      ############################################################################
      #"using hybrids data" (Gini Simpson, 1- Atkinson, singletons)
      if(nrow(disth)>0){
        SIMP_HYBD = abdiv::dominance(df_hybrid[, i])
        ATK_HYBD = 1-DescTools::Atkinson(df_hybrid[, i],parameter = 0.5)
        SP1_HYBD = length(df_hybrid[, i][which(df_hybrid[, i] == 1)])
        PROP_1_SP1_HYBD = 1 - (length(df_hybrid[, i][which(df_hybrid[, i] ==
                                                             1)]) /
                                 (nrow(dist_prop)))
        SPP_HYBD = length(df_hybrid[, i][which(df_hybrid[, i] > 0)])
        PROP_SPP_HYBD = length(df_hybrid[, i][which(df_hybrid[, i] > 0)]) /
          (nrow(dist_prop))
        REC_HYBD = sum(df_hybrid[, i][which(df_hybrid[, i] > 0)],na.rm = T)
        PROP_R_HYBD = sum(df_hybrid[, i][which(df_hybrid[, i] > 0)]) / sum(df_wild, df_landrace, df_hybrid,na.rm = T)
        ENS_GINI_HYBD = 1 / abdiv::dominance(df_hybrid[, i])
        MARGALEF_HYBD = abdiv::margalef(df_hybrid[,i])
      } else {
        SIMP_HYBD = NA
        ATK_HYBD = NA
        SP1_HYBD= NA
        PROP_1_SP1_HYBD =NA
        SPP_HYBD = NA
        PROP_SPP_HYBD = NA
        REC_HYBD = NA
        PROP_R_HYBD = NA
        ENS_GINI_HYBD = NA
        MARGALEF_HYBD = NA
      }
      ############################################################################
      #summarizing
      df <- data.frame(
        COUNTRY = colnames(df_total)[i],
        SIMP_WILD = SIMP_WILD,
        SIMP_LAND = SIMP_LAND,
        
        SIMP_HYBD = SIMP_HYBD,
        ATK_WILD = ATK_WILD,
        ATK_LAND = ATK_LAND,
        ATK_HYBD = ATK_HYBD,
        SP1_WILD = SP1_WILD,
        SP1_LAND = SP1_LAND,
        SP1_HYBD = SP1_HYBD,
        PROP_1_SP1_WILD = PROP_1_SP1_WILD,
        PROP_1_SP1_LAND = PROP_1_SP1_LAND,
        PROP_1_SP1_HYBD = PROP_1_SP1_HYBD,
        COVERAGE_WILD=cov_w,
        COVERAGE_LAND=cov_l,
        COVERAGE_HYBD=cov_h,
        SPP_WILD = SPP_WILD,
        SPP_LAND = SPP_LAND,
        SPP_HYBD = SPP_HYBD,
        PROP_SPP_WILD = PROP_SPP_WILD,
        PROP_SPP_LAND = PROP_SPP_LAND,
        PROP_SPP_HYBD = PROP_SPP_HYBD,
        REC_WILD = REC_WILD,
        REC_LAND = REC_LAND,
        REC_HYBD = REC_HYBD,
        PROP_R_WILD = PROP_R_WILD,
        PROP_R_LAND = PROP_R_LAND,
        PROP_R_HYBD = PROP_R_HYBD,
        GRIN_TAXA_PROP_W = G_W_i,
        GRIN_TAXA_PROP_L = G_L_i,
        GRIN_TAXA_PROP_H = G_H_i,
        ENS_GINI_WILD = ENS_GINI_WILD,
        ENS_GINI_LAND = ENS_GINI_LAND,
        ENS_GINI_HYBD = ENS_GINI_HYBD,
        MARGALEF_WILD = MARGALEF_WILD,
        MARGALEF_LAND = MARGALEF_LAND,
        MARGALEF_HYBD = MARGALEF_HYBD
      )
      
      simp_list[[i]] <- df
      
      
    }
    ################################################################################
    # list_w <- list()
    # for(i in 1:ncol(df_wild)){
    #   if(sum(df_wild[i,])>0){
    #     list_w[[i]] <- df_wild[,i][which(df_wild[,i]>0)]
    #   } else {
    #     list_w[[i]] <- NA
    #   }
    # };rm(i)
    
    #specpool2vect(t(df_wild),index="chao")
    
    
    
    
    #   t()
    # names(list_w) <- colnames(df_wild)
    # list_w <- list_w[!is.na(list_w)]
    # #out2 <- iNEXT(list_w, datatype="incidence_freq",q = c(0,1,2))
    # out2 <- iNEXT(list_w, datatype="abundance",q=2)#, 
    # out2$DataInfo
    # out2$iNextEst
    # out3 <- estimateD(list_w, datatype="abundance", base="coverage",q = 1) 
    # #          iNEXT="coverage",q=1)
    # ggiNEXT(out2, type=3)#, facet.var="site")
    #ggiNEXT(x, type=1, se=TRUE, facet.var="none", color.var="site", grey=FALSE)  
    
    ###############################################################################
    simp_list <- do.call(rbind, simp_list)
    simp_list$COUNTRY[which(simp_list$COUNTRY == "NA1")] <- "No_country"
    ################################################################################
    # calculating alternative to Simpson (ACE)  
    if(is.data.frame(df_wild)){
      ACE_w_list <- list()
      for(i in 1:ncol(df_wild)){
        if(sum(df_wild[,i])>0){
          ACE_w_list[[i]] <- ACE(df_wild[,i],taxa.row = F)  
        } else {
          ACE_w_list[[i]] <- NA
        }
      };rm(i)
    } else {
      ACE_w_list <- NA
    }
    if(is.data.frame(df_landrace)){
      ACE_l_list <- list()
      for(i in 1:ncol(df_landrace)){
        if(sum(df_landrace[,i])>0){
          ACE_l_list[[i]] <- ACE(df_landrace[,i],taxa.row = F)  
        } else {
          ACE_l_list[[i]] <- NA
        }
      };rm(i)
    } else {
      ACE_l_list <- NA
    }
    if(is.data.frame(df_hybrid)){
      ACE_h_list <- list()
      for(i in 1:ncol(df_hybrid)){
        if(sum(df_hybrid[,i])>0){
          ACE_h_list[[i]] <- ACE(df_hybrid[,i],taxa.row = F)  
        } else {
          ACE_h_list[[i]] <- NA
        }
      };rm(i) 
    } else {
      ACE_h_list <- NA
    }
    simp_list$ACE_WILD <- unlist(ACE_w_list)
    simp_list$ACE_LAND <- unlist(ACE_l_list)
    simp_list$ACE_HYBD <- unlist(ACE_h_list)
    ################################################################################
    #calculating index 1 considering Gini Simpson, 1-Atkinson, and Singleton
    # simp_list$COMP_INDEX_W <- (
    #   simp_list$SIMP_WILD + (simp_list$ATK_WILD * simp_list$PROP_1_SP1_WILD) +
    #     simp_list$GRIN_TAXA_PROP_W
    # ) / 3
    # simp_list$COMP_INDEX_L <- (
    #   simp_list$SIMP_LAND + (simp_list$ATK_LAND * simp_list$PROP_1_SP1_LAND) +
    #     simp_list$GRIN_TAXA_PROP_L
    # ) / 3
    # simp_list$COMP_INDEX_H <- (
    #   simp_list$SIMP_HYBD + (simp_list$ATK_HYBD * simp_list$PROP_1_SP1_HYBD) +
    #     simp_list$GRIN_TAXA_PROP_H
    # ) / 3
    ################################################################################
    message("Adding collection data")
    simp_list[(nrow(simp_list) + 1), "COUNTRY"] <- "Subset"
    
    simp_list[(nrow(simp_list)), "SIMP_WILD"] <- TAX_DF$Simpson[1]
    simp_list[(nrow(simp_list)), "SIMP_LAND"] <- TAX_DF$Simpson[2]
    simp_list[(nrow(simp_list)), "SIMP_HYBD"]  <- TAX_DF$Simpson[3]
    
    simp_list[(nrow(simp_list)), "ATK_WILD"]  <- TAX_DF$X1_Atkinson[1]
    simp_list[(nrow(simp_list)), "ATK_LAND"]   <- TAX_DF$X1_Atkinson[2]
    simp_list[(nrow(simp_list)), "ATK_HYBD"]   <- TAX_DF$X1_Atkinson[3]
    
    simp_list[(nrow(simp_list)), "SP1_WILD"]  <- TAX_DF$nRecords[1]
    simp_list[(nrow(simp_list)), "SP1_LAND"]   <- TAX_DF$nRecords[2]
    simp_list[(nrow(simp_list)), "SP1_HYBD"]   <- TAX_DF$nRecords[3]
    
    simp_list[(nrow(simp_list)), "PROP_R_WILD"]  <- TAX_DF$propRecords[1]
    simp_list[(nrow(simp_list)), "PROP_R_LAND"]   <- TAX_DF$propRecords[2]
    simp_list[(nrow(simp_list)), "PROP_R_HYBD"]   <- TAX_DF$propRecords[3]
    
    simp_list[(nrow(simp_list)), "PROP_1_SP1_WILD"]   <- TAX_DF$propSingleton[1]
    simp_list[(nrow(simp_list)), "PROP_1_SP1_LAND"]   <- TAX_DF$propSingleton[2]
    simp_list[(nrow(simp_list)), "PROP_1_SP1_HYBD"]    <- TAX_DF$propSingleton[3]
    
    simp_list[(nrow(simp_list)), "GRIN_TAXA_PROP_W"]   <- TAX_DF$propTaxa_GRIN[1]
    simp_list[(nrow(simp_list)), "GRIN_TAXA_PROP_L"]   <- TAX_DF$propTaxa_GRIN[2]
    simp_list[(nrow(simp_list)), "GRIN_TAXA_PROP_H"]    <- TAX_DF$propTaxa_GRIN[3]
    
    # simp_list[(nrow(simp_list)), "COMP_INDEX_W"]    <- NA#TAX_DF$Composition_index[1]
    # simp_list[(nrow(simp_list)), "COMP_INDEX_L"]    <- NA#TAX_DF$Composition_index[2]
    # simp_list[(nrow(simp_list)), "COMP_INDEX_H"]    <- NA#TAX_DF$Composition_index[3]
    # 
    simp_list[(nrow(simp_list)), "PROP_SPP_WILD"]    <- TAX_DF$propTaxa_collection[1]
    simp_list[(nrow(simp_list)), "PROP_SPP_LAND"]    <- TAX_DF$propTaxa_collection[2]
    simp_list[(nrow(simp_list)), "PROP_SPP_HYBD"]    <- TAX_DF$propTaxa_collection[3]
    
    simp_list[(nrow(simp_list)), "SPP_WILD"]    <- TAX_DF$nTaxa[1]
    simp_list[(nrow(simp_list)), "SPP_LAND"]    <- TAX_DF$nTaxa[2]
    simp_list[(nrow(simp_list)), "SPP_HYBD"]    <- TAX_DF$nTaxa[3]
    
    simp_list[(nrow(simp_list)), "REC_WILD"]    <- TAX_DF$nRecords[1]
    simp_list[(nrow(simp_list)), "REC_LAND"]    <- TAX_DF$nRecords[2]
    simp_list[(nrow(simp_list)), "REC_HYBD"]    <- TAX_DF$nRecords[3]
    
    ##############################################################################
    message("Calculating 1/D as Effective number of species (Hill q=2)")
    #1/D (wild)
    if(nrow(distw)>0){
      simp_list[(nrow(simp_list)), "ENS_GINI_WILD"]    <- 1/abdiv::dominance(rowSums(df_wild,na.rm = T))  
    } else {
      simp_list[(nrow(simp_list)), "ENS_GINI_WILD"]    <- NA
    }
    #1/D (landrace)
    if(nrow(distl)>0){
      simp_list[(nrow(simp_list)), "ENS_GINI_LAND"]    <- 1/abdiv::dominance(rowSums(df_landrace,na.rm = T))  
    } else {
      simp_list[(nrow(simp_list)), "ENS_GINI_LAND"]    <- NA
    }
    #1/D (hybrids)
    if(nrow(disth)>0){
      simp_list[(nrow(simp_list)), "ENS_GINI_HYBD"]    <- 1/abdiv::dominance(rowSums(df_hybrid,na.rm = T))  
    } else {
      simp_list[(nrow(simp_list)), "ENS_GINI_HYBD"]    <- NA
    }
    ################################################################################
    message("Calculating ACE for collection")
    
    #ACE (wild)
    if(nrow(distw)>0){
      simp_list[(nrow(simp_list)), "ACE_WILD"]    <- ACE(rowSums(df_wild),taxa.row = F)  
    } else {
      simp_list[(nrow(simp_list)), "ACE_WILD"]    <- NA
    }
    #ACE (landrace)
    if(nrow(distl)>0){
      simp_list[(nrow(simp_list)), "ACE_LAND"]    <- ACE(rowSums(df_landrace),taxa.row = F)  
    } else {
      simp_list[(nrow(simp_list)), "ACE_LAND"]    <- NA
    }
    #ACE (hybrid)
    if(nrow(disth)>0){
      simp_list[(nrow(simp_list)), "ACE_HYBD"]    <- ACE(rowSums(df_hybrid),taxa.row = F)   
    } else {
      simp_list[(nrow(simp_list)), "ACE_HYBD"]    <- NA
    }
    ################################################################################
    message("Calculating coverage for collection")  
    
    #1/D (wild)
    if(nrow(distw)>0){
      simp_list[(nrow(simp_list)), "ENS_GINI_WILD"]    <- 1/abdiv::dominance(rowSums(df_wild,na.rm = T))  
    } else {
      simp_list[(nrow(simp_list)), "ENS_GINI_WILD"]    <- NA
    }
    #1/D (landrace)
    if(nrow(distl)>0){
      simp_list[(nrow(simp_list)), "ENS_GINI_LAND"]    <- 1/abdiv::dominance(rowSums(df_landrace,na.rm = T))  
    } else {
      simp_list[(nrow(simp_list)), "ENS_GINI_LAND"]    <- NA
    }
    #1/D (hybrids)
    if(nrow(disth)>0){
      simp_list[(nrow(simp_list)), "ENS_GINI_HYBD"]    <- 1/abdiv::dominance(rowSums(df_hybrid,na.rm = T))  
    } else {
      simp_list[(nrow(simp_list)), "ENS_GINI_HYBD"]    <- NA
    }
    ################################################################################
    message("Calculating coverage for collection")
    
    #ACE (wild)
    if(nrow(distw)>0){
      f1_w <- length(rowSums(df_wild)[which(rowSums(df_wild)==1)])
      f2_w <- length(rowSums(df_wild)[which(rowSums(df_wild)==2)])
      n_w <- sum(df_wild,na.rm = T)
      cov_w <-  1-(f1_w/n_w)* (((n_w-1)*f1_w)/ (((n_w-1)*f1_w)+(2*f2_w)))
      simp_list[(nrow(simp_list)), "COVERAGE_WILD"]    <- cov_w
    } else {
      simp_list[(nrow(simp_list)), "COVERAGE_WILD"]    <- NA
    }
    #ACE (landrace)
    if(nrow(distl)>0){
      f1_l<- length(rowSums(df_landrace)[which(rowSums(df_landrace)==1)])
      f2_l <- length(rowSums(df_landrace)[which(rowSums(df_landrace)==2)])
      n_l <- sum(df_landrace,na.rm = T)
      cov_l <-  1-(f1_l/n_l)* (((n_l-1)*f1_l)/ (((n_l-1)*f1_l)+(2*f2_l)))
      
      simp_list[(nrow(simp_list)), "COVERAGE_LAND"]    <- cov_l
    } else {
      simp_list[(nrow(simp_list)), "COVERAGE_LAND"]    <- NA
    }
    #ACE (hybrid)
    if(nrow(disth)>0){
      f1_h <- length(rowSums(df_hybrid)[which(rowSums(df_hybrid)==1)])
      f2_h <- length(rowSums(df_hybrid)[which(rowSums(df_hybrid)==2)])
      n_h <- sum(df_hybrid,na.rm = T)
      cov_h <-  1-(f1_h/n_h)* (((n_h-1)*f1_h)/ (((n_h-1)*f1_h)+(2*f2_h)))
      
      simp_list[(nrow(simp_list)), "COVERAGE_HYBD"]    <- cov_h
    } else {
      simp_list[(nrow(simp_list)), "COVERAGE_HYBD"]    <- NA
    }
    ################################################################################
    message("Calculating Margalef index for collection")  
    
    #1/D (wild)
    if(nrow(distw)>0){
      simp_list[(nrow(simp_list)), "MARGALEF_WILD"]    <- abdiv::margalef(rowSums(df_wild,na.rm = T))  
    } else {
      simp_list[(nrow(simp_list)), "MARGALEF_WILD"]    <- NA
    }
    #1/D (landrace)
    if(nrow(distl)>0){
      simp_list[(nrow(simp_list)), "MARGALEF_LAND"]    <- abdiv::margalef(rowSums(df_landrace,na.rm = T))  
    } else {
      simp_list[(nrow(simp_list)), "MARGALEF_LAND"]    <- NA
    }
    #1/D (hybrids)
    if(nrow(disth)>0){
      simp_list[(nrow(simp_list)), "MARGALEF_HYBD"]    <- abdiv::margalef(rowSums(df_hybrid,na.rm = T))  
    } else {
      simp_list[(nrow(simp_list)), "MARGALEF_HYBD"]    <- NA
    }
    
    ################################################################################
    message("Normalizing ACE")
    #https://medium.com/@noorfatimaafzalbutt/a-comprehensive-guide-to-normalization-in-machine-learning-afead759b062
    simp_list$ACE_WILD_S <- 
      decimal_scale(simp_list$ACE_WILD)
      #(simp_list$ACE_WILD - min(simp_list$ACE_WILD,na.rm = T))/
      #(max(simp_list$ACE_WILD,na.rm = T)-min(simp_list$ACE_WILD,na.rm = T))

    simp_list$ACE_LAND_S <- 
      decimal_scale(simp_list$ACE_LAND)
      #(simp_list$ACE_LAND - min(simp_list$ACE_LAND,na.rm = T))/
      #(max(simp_list$ACE_LAND,na.rm = T)-min(simp_list$ACE_LAND,na.rm = T))

    
    simp_list$ACE_HYBD_S <- 
      decimal_scale(simp_list$ACE_HYBD)
      #(simp_list$ACE_HYBD - min(simp_list$ACE_HYBD,na.rm = T))/
      #(max(simp_list$ACE_HYBD,na.rm = T)-min(simp_list$ACE_HYBD,na.rm = T))
    ################################################################################
    message("Normalizing ENS (Simpson)")
    #https://abagen.readthedocs.io/en/stable/user_guide/normalization.html#usage-norm-zscore
    
    
    
    simp_list$ENS_GINI_WILD_S <- 
    sigmoid_scale(simp_list$ENS_GINI_WILD)
      

  
   #(median(simp_list$ENS_GINI_WILD,na.rm = T)-simp_list$ENS_GINI_WILD)/
    #mad(simp_list$ENS_GINI_WILD,na.rm = T)
      #decimal_scale(simp_list$ENS_GINI_WILD)
     #(simp_list$ENS_GINI_WILD - min(simp_list$ENS_GINI_WILD,na.rm = T))/
     #(max(simp_list$ENS_GINI_WILD,na.rm = T)-min(simp_list$ENS_GINI_WILD,na.rm = T))
    
    simp_list$ENS_GINI_LAND_S <- 
      sigmoid_scale(simp_list$ENS_GINI_LAND)
      #decimal_scale(simp_list$ENS_GINI_LAND)
      #(median(simp_list$ENS_GINI_LAND,na.rm = T)-simp_list$ENS_GINI_LAND)/
      #mad(simp_list$ENS_GINI_LAND,na.rm = T)
      #(simp_list$ENS_GINI_LAND - min(simp_list$ENS_GINI_LAND,na.rm = T))/
      #(max(simp_list$ENS_GINI_LAND,na.rm = T)-min(simp_list$ENS_GINI_LAND,na.rm = T))
    
    simp_list$ENS_GINI_HYBD_S <- 
      sigmoid_scale(simp_list$ENS_GINI_HYBD)
     #  (median(simp_list$ENS_GINI_HYBD,na.rm = T)-simp_list$ENS_GINI_HYBD)/
     #  mad(simp_list$ENS_GINI_HYBD,na.rm = T)
     # # decimal_scale(simp_list$ENS_GINI_HYBD)
    
      #(simp_list$ENS_GINI_HYBD - min(simp_list$ENS_GINI_HYBD,na.rm = T))/
      #(max(simp_list$ENS_GINI_HYBD,na.rm = T)-min(simp_list$ENS_GINI_HYBD,na.rm = T))
    
    ################################################################################
    message("Normalizing MARGALEF")
    
    simp_list$MARGALEF_WILD_S <- 
      sigmoid_scale(simp_list$MARGALEF_WILD)
    
      #decimal_scale(simp_list$MARGALEF_WILD)
      #(simp_list$MARGALEF_WILD - min(simp_list$MARGALEF_WILD,na.rm = T))/
      #(max(simp_list$MARGALEF_WILD,na.rm = T)-min(simp_list$MARGALEF_WILD,na.rm = T))
    
    simp_list$MARGALEF_LAND_S <- 
      sigmoid_scale(simp_list$MARGALEF_LAND)
    
      #decimal_scale(simp_list$MARGALEF_LAND)
      #(simp_list$MARGALEF_LAND - min(simp_list$MARGALEF_LAND,na.rm = T))/
      #(max(simp_list$MARGALEF_LAND,na.rm = T)-min(simp_list$MARGALEF_LAND,na.rm = T))
    
      simp_list$MARGALEF_HYBD_S <- 
        sigmoid_scale(simp_list$MARGALEF_HYBD)
      #decimal_scale(simp_list$MARGALEF_HYBD)
      #(simp_list$MARGALEF_HYBD - min(simp_list$MARGALEF_HYBD,na.rm = T))/
      #(max(simp_list$MARGALEF_HYBD,na.rm = T)-min(simp_list$MARGALEF_HYBD,na.rm = T))
    
    ################################################################################
    message("Calculating index as ACE * ATKINSON")
    simp_list$COMP_INDEX_W_ACE <- (
      #simp_list$ENS_GINI_WILD_S + (simp_list$ATK_WILD * simp_list$PROP_1_SP1_WILD) +
      #simp_list$GRIN_TAXA_PROP_W
      #) / 3
      
      simp_list$ACE_WILD_S * simp_list$ATK_WILD 
    )
    
    simp_list$COMP_INDEX_L_ACE <- (
      #simp_list$ENS_GINI_LAND_S + (simp_list$ATK_LAND * simp_list$PROP_1_SP1_LAND) +
      #simp_list$GRIN_TAXA_PROP_L
      #) / 3
      simp_list$ACE_LAND_S * simp_list$ATK_LAND 
    )
    
    simp_list$COMP_INDEX_H_ACE <- (
      #simp_list$ENS_GINI_HYBD_S+ (simp_list$ATK_HYBD * simp_list$PROP_1_SP1_HYBD) +
      #simp_list$GRIN_TAXA_PROP_H
      #) / 3
      simp_list$ACE_HYBD_S * simp_list$ATK_HYBD
    )
    ################################################################################
    message("Calculating index as ENS * ATKINSON")
    simp_list$COMP_INDEX_W_ENS <- (
      #simp_list$ENS_GINI_WILD_S + (simp_list$ATK_WILD * simp_list$PROP_1_SP1_WILD) +
      #simp_list$GRIN_TAXA_PROP_W
      #) / 3
      
      simp_list$ENS_GINI_WILD_S * simp_list$ATK_WILD 
    )
    
    simp_list$COMP_INDEX_L_ENS <- (
      #simp_list$ENS_GINI_LAND_S + (simp_list$ATK_LAND * simp_list$PROP_1_SP1_LAND) +
      #simp_list$GRIN_TAXA_PROP_L
      #) / 3
      simp_list$ENS_GINI_LAND_S * simp_list$ATK_LAND 
    )
    
    simp_list$COMP_INDEX_H_ENS <- (
      #simp_list$ENS_GINI_HYBD_S+ (simp_list$ATK_HYBD * simp_list$PROP_1_SP1_HYBD) +
      #simp_list$GRIN_TAXA_PROP_H
      #) / 3
      simp_list$ENS_GINI_HYBD_S * simp_list$ATK_HYBD
    )
    ################################################################################
    message("Calculating index as MARGALEF * ATKINSON")
    simp_list$COMP_INDEX_W_MARGALEF <- (
      #simp_list$ENS_GINI_WILD_S + (simp_list$ATK_WILD * simp_list$PROP_1_SP1_WILD) +
      #simp_list$GRIN_TAXA_PROP_W
      #) / 3
      
      simp_list$MARGALEF_WILD_S * simp_list$ATK_WILD 
    )
    
    simp_list$COMP_INDEX_L_MARGALEF <- (
      #simp_list$ENS_GINI_LAND_S + (simp_list$ATK_LAND * simp_list$PROP_1_SP1_LAND) +
      #simp_list$GRIN_TAXA_PROP_L
      #) / 3
      simp_list$MARGALEF_LAND_S * simp_list$ATK_LAND 
    )
    
    simp_list$COMP_INDEX_H__MARGALEF <- (
      #simp_list$ENS_GINI_HYBD_S+ (simp_list$ATK_HYBD * simp_list$PROP_1_SP1_HYBD) +
      #simp_list$GRIN_TAXA_PROP_H
      #) / 3
      simp_list$MARGALEF_HYBD_S * simp_list$ATK_HYBD
    )
    
    ################################################################################
    row.names(simp_list)  <- simp_list$COUNTRY
    
    ################################################################################
    simp_list$REGION <- NA
    simp_list$REGION <- countrycode(simp_list$COUNTRY,
                                    origin = 'iso3c',
                                    destination = 'region')
    
    simp_list$REGION[which(is.na(simp_list$REGION))] <- "N/A"
    simp_list$REGION[which(simp_list$COUNTRY == "Subset")] <- "Subset"
    simp_list$REGION[which(simp_list$COUNTRY == "No_country")] <- "N/A"
    simp_list$REGION[which(simp_list$COUNTRY == "SCG")] <- "Europe & Central Asia"
    simp_list$REGION[which(simp_list$COUNTRY == "YUG")] <- "Europe & Central Asia"
    simp_list$REGION[which(simp_list$COUNTRY == "ZAR")] <- "Sub-Saharan Africa"
    simp_list$REGION[which(simp_list$COUNTRY == "BUR")] <- "East Asia & Pacific"
    
    #row.names(simp_list)
    #simp_list_NA <- simp_list[which(simp_list$COUNTRY=="NA1"),]
    #simp_list_nona <- simp_list[which(simp_list$COUNTRY!="NA1"),]
    
    # mean(weighted.mean(simp_list_nona$COMP_INDEX_W,simp_list_nona$PROP_R_WILD),
    #      weighted.mean(simp_list_nona$COMP_INDEX_L,simp_list_nona$PROP_R_LAND)
    #      )
    # #cor.test(simp_list$COMP_INDEX_W,simp_list$SIMP_WILD)
    # mean(mean(simp_list_nona$COMP_INDEX_W,na.rm = T),
    #      mean(simp_list_nona$COMP_INDEX_L,na.rm = T))
    # #row.names(simp_list_nona) <- simp_list_nona$COUNTRY
    ############################################################################
    message(paste0("Plotting for: ",plot_type))
    ############################################################################
    if(plot_type=="Margalef"){
    x_res_W <- #simp_list_nona[,c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27)]
      simp_list[, c(#2, #3,
        "ATK_WILD",
        #"ENS_GINI_WILD_S",
        "MARGALEF_WILD_S",
        #                "COVERAGE_WILD",
        #"PROP_SPP_WILD",
        #"PROP_R_WILD",
        #"GRIN_TAXA_PROP_W",
        #                "ACE_WILD_S",
        #"COMP_INDEX_W_ENS",
        #                "COMP_INDEX_H_ACE",
        "COMP_INDEX_W_MARGALEF",
        "REGION")]
    ##############################################################################
    colnames(x_res_W) <- c("1-Atkinson",
                           #"Effective number of species (normalized)",
                           "Margalef index (normalized)",
                           #"Proportion of CWR taxa observed",
                           #"Proportion of CWR records for the collection",
                           #"Proportion of taxa available in GRIN"
                           "Composition index",
                           "REGION"
                           )
    
    } else if(plot_type=="Simpson"){
      x_res_W <- #simp_list_nona[,c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27)]
        simp_list[, c(#2, #3,
          "ATK_WILD",
          "ENS_GINI_WILD_S",
          #"MARGALEF_WILD_S",
          #                "COVERAGE_WILD",
          #"PROP_SPP_WILD",
          #"PROP_R_WILD",
          #"GRIN_TAXA_PROP_W",
          #                "ACE_WILD_S",
          #"COMP_INDEX_W_ENS",
          #                "COMP_INDEX_H_ACE",
          "COMP_INDEX_W_MARGALEF",
          "REGION")]
      ##############################################################################
      colnames(x_res_W) <- c("1-Atkinson",
                             "Effective number of species (normalized)",
                            # "Margalef index (normalized)",
                             #"Proportion of CWR taxa observed",
                             #"Proportion of CWR records for the collection",
                             #"Proportion of taxa available in GRIN"
                             "Composition index",
                             "REGION"
      )
    }

    x_res_W <- x_res_W[complete.cases(x_res_W),]
    x_res_W$REGION <- factor(x_res_W$REGION)
    ##############################################################################
    if(plot_type=="Margalef"){
    x_res_L <- #simp_list_nona[,c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27)]
      simp_list[, c(#2, #3,
        "ATK_LAND",
        #"ENS_GINI_LAND_S",
        "MARGALEF_LAND_S",
        #                "COVERAGE_WILD",
        #"PROP_SPP_WILD",
        #"PROP_R_WILD",
        #"GRIN_TAXA_PROP_W",
        #                "ACE_WILD_S",
        #"COMP_INDEX_L_ENS",
        #                "COMP_INDEX_H_ACE",
        "COMP_INDEX_L_MARGALEF",
        "REGION")]
    ##############################################################################
    colnames(x_res_L) <- c("1-Atkinson",
                           #"Effective number of species (normalized)",
                           "Margalef index (normalized)",
                           #"Proportion of CWR taxa observed",
                           #"Proportion of CWR records for the collection",
                           #"Proportion of taxa available in GRIN"
                           "Composition index",
                           "REGION"
    )
    } else if(plot_type=="Simpson"){
      x_res_L <- #simp_list_nona[,c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27)]
        simp_list[, c(#2, #3,
          "ATK_LAND",
          "ENS_GINI_LAND_S",
          #"MARGALEF_LAND_S",
          #                "COVERAGE_WILD",
          #"PROP_SPP_WILD",
          #"PROP_R_WILD",
          #"GRIN_TAXA_PROP_W",
          #                "ACE_WILD_S",
          #"COMP_INDEX_L_ENS",
          #                "COMP_INDEX_H_ACE",
          "COMP_INDEX_L_MARGALEF",
          "REGION")]
      ##############################################################################
      colnames(x_res_L) <- c("1-Atkinson",
                             "Effective number of species (normalized)",
                             #"Margalef index (normalized)",
                             #"Proportion of CWR taxa observed",
                             #"Proportion of CWR records for the collection",
                             #"Proportion of taxa available in GRIN"
                             "Composition index",
                             "REGION"
      )
    }
    x_res_L <- x_res_L[complete.cases(x_res_L), ]
    x_res_L$REGION <- factor(x_res_L$REGION)
    
    ################################################################################
    #seed for NIPALS
    set.seed(1000)
    ################################################################################
    message("plotting results for wild")
    
    if(nrow(x_res_W)>=2){
      x_W_A = nipals(x_res_W[, -ncol(x_res_W)], comps = 3,scaled = T)
      # x_W <- FactoMineR::PCA(X = x_res_W[, -ncol(x_res_W)],
      #                        graph = F,
      #                        scale.unit = T,ncp = 3)
      # 
      
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
        ggtitle(paste0("[",plot_type," index] ","Crop Wild Relatives composition (", nrow(x_res_W_NA), " countries)"))+
        
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
      # plot(x_W_A,)
      # plot(x_W_A, what="observations", show.names=TRUE, main="Missing data")
      # 
      # #eigen values (summary)
      # x_W$eig <- x_W_A$values
      # colnames(x_W$eig) <- c("eigenvalue","percentage of variance","cumulative percentage of variance")
      # row.names(x_W$eig) <- paste("comp",1:nrow(x_W_A$values))
      # 
      # #variables results
      # x_W$var$coord <- x_W_A$cor.xt
      # colnames(x_W$var$coord) <-paste0("Dim.",1:3)
      # x_W$var$cor  <- x_W_A$loadings
      # colnames(x_W$var$coord) <-paste0("Dim.",1:3)
      # x_W$var$cos2 <- NA
      # x_W$var$contrib <- NA
      # #ind results
      # x_W$ind$coord <- x_W_A$scores
      # colnames(x_W$ind$coord) <-paste0("Dim.",1:3)
      # x_W$ind$cos2 <- x_W_A$cos 
      # colnames(x_W$ind$cos2 ) <-paste0("Dim.",1:3)
      # x_W$ind$contrib <- x_W_A$contrib
      # colnames(x_W$ind$contrib) <-paste0("Dim.",1:3)
      # x_W$svd$vs <- NA
      # x_W$svd$U <- NA
      # x_W$svd$V <- NA
      # 
      # x_W$call$centre <- NA
      # x_W$call$ecart.type <- NA
      # REGION <- x_res_W$REGION
      # PCA_WILD <- fviz_pca(
      #   x_W,
      #   title = paste0(
      #     "Crop Wild Relatives composition (",
      #     nrow(x_res_W) - 2,
      #     " countries)"
      #   ),
      #   #addEllipses = TRUE, ellipse.type = "convex",
      #   
      #   repel = TRUE#,
      #   #geom.var = c("text", "point"),
      #   #col.var = "black"
      # ) +
      #   geom_point(aes(fill = REGION), shape = 21, size = 6) +
      #   scale_fill_brewer(palette = "Set2") +
      #   theme_minimal()
      
      
      #PCA_WILD
      
      ggplot2::ggsave(
        paste0(
          outdir,
          "/",
          collection_name,
          "/",
          collection_name,
          "_",
          plot_type,
         
          "_composition_WILD.png"
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
      #x_L <- FactoMineR::PCA(x_res_L[, -ncol(x_res_L)], graph = F)
      
      
      REGION <- x_res_L$REGION
      # PCA_LAND <- fviz_pca(
      #   x_L,
      #   title = paste0("Landraces composition (", nrow(x_res_L) -
      #                    2, " countries)"),
      #   #addEllipses = TRUE, ellipse.type = "convex",
      #   
      #   repel = TRUE#,
      #   #geom.var = c("text", "point"),
      #   #col.var = "black"
      # ) +
      #   geom_point(aes(fill = REGION), shape = 21, size = 6) +
      #   scale_fill_brewer(palette = "Set2") +
      #   theme_minimal()
      # 
      if(length(unique(x_res_L[,1]))>1){
      x_W_A = nipals(x_res_L[, -ncol(x_res_L)], comps = 3,scaled = T)
      # x_W <- FactoMineR::PCA(X = x_res_W[, -ncol(x_res_W)],
      #                        graph = F,
      #                        scale.unit = T,ncp = 3)
      # 
      
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
        ggtitle(paste0("[",plot_type," index] ","Landraces composition (", nrow(x_res_L_NA), " countries)"))+
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
          #xlim = c(NA, Inf),
          # Do not repel from top or bottom edges.
          #ylim = c(-Inf, Inf),
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
          plot_type,
          "_composition_LANDRACE.png"
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
    #saving all
    write.csv(
      simp_list,
      paste0(
        outdir,
        "/",
        collection_name,
        "/",
        collection_name,
        "_composition_countries.csv"
      ),
      row.names = F,
      na = ""
    )
    
    message(
      paste(collection_name, "processed for countries composition!")
    )
    
    return(simp_list)
  }
  ################################################################################
  # plot(simp_list$ACE_WILD_S,simp_list$COMP_INDEX_W_S)
  # cor.test(simp_list$ACE_WILD_S,simp_list$COMP_INDEX_W_S)
  ################################################################################
  dir <- "D:/OneDrive - CGIAR/GERMPLASM_INDEX"
  #dir <- "D:/ONEDRIVE/cgiar/OneDrive - CGIAR/GERMPLASM_INDEX"
  ################################################################################
  #outdir
  outdir <- paste0(dir, "/BEANS/RESULTS")
  ################################################################################
  collection_name <- "cassava"
  x1 <- composition_country(outdir, collection_name,plot_type="Margalef")
  x1 <- composition_country(outdir, collection_name,plot_type="Simpson")
  
  x1$collection <-  collection_name
  ################################################################################
  collection_name <- "beans"
  x2 <- composition_country(outdir, collection_name,plot_type="Margalef")
  x2 <- composition_country(outdir, collection_name,plot_type="Simpson")
  x2$collection <-  collection_name
  ################################################################################
  collection_name <- "forages"
  x3 <- composition_country(outdir, collection_name,plot_type="Margalef")
  x3 <- composition_country(outdir, collection_name,plot_type="Simpson")
  x3$collection <-  collection_name
  ################################################################################