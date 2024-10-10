################################################################################
library(countrycode)
library(ggplot2)
library(factoextra)
library(glue)
################################################################################
################################################################################

composition_country <- function(outdir, collection_name) {
  message("reading  previous results")
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
    #i <- 1
    if(nrow(distw)>0){
      distw$coun <- NA
      distw$coun <- df_wild[, i]
      distw$coun[which(distw$coun > 0)] <- 1
      distw$coun <- distw$coun * distw$GRIN
      G_W_i <- sum(distw$coun) / nrow(dist_prop)
    } else {
      G_W_i <- NA
    }    
    if(nrow(distl)>0){
      distl$coun <- NA
      distl$coun <- df_landrace[, i]
      distl$coun[which(distl$coun > 0)] <- 1
      distl$coun <- distl$coun * distl$GRIN
      G_L_i <- sum(distl$coun) / nrow(dist_prop) 
    } else {
      G_L_i <- NA
    }
    if(nrow(disth)>0){
    disth$coun <- NA
    disth$coun <- df_hybrid[, i]
    disth$coun[which(disth$coun > 0)] <- 1
    disth$coun <- disth$coun * disth$GRIN
    G_H_i <- sum(disth$coun) / nrow(dist_prop)
    disth$coun <- NULL
    } else {
      G_H_i <- NA
    }
    
    #i <- 5
    ############################################################################
    #summing species
    ############################################################################
    #"using wild data"
    if(nrow(distw)>0){
      SIMP_WILD = abdiv::simpson(df_wild[, i])
      ATK_WILD = 1 - DescTools::Atkinson(df_wild[, i])
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
    }
    ############################################################################
    #"using landraces data"
    if(nrow(distl)>0){
      SIMP_LAND = abdiv::simpson(df_landrace[, i])
      ATK_LAND = 1 - DescTools::Atkinson(df_landrace[, i])
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
    }
    ############################################################################
    #"using hybrids data"
    if(nrow(disth)>0){
      SIMP_HYBD = abdiv::simpson(df_hybrid[, i])
      ATK_HYBD = 1 - DescTools::Atkinson(df_hybrid[, i])
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
      ENS_GINI_HYBD = ENS_GINI_HYBD
    )
    
    simp_list[[i]] <- df
    
    
  }
  
  
  
  simp_list <- do.call(rbind, simp_list)
  simp_list$COUNTRY[which(simp_list$COUNTRY == "NA1")] <- "No_country"
  ################################################################################
  simp_list$COMP_INDEX_W <- (
    simp_list$SIMP_WILD + (simp_list$ATK_WILD * simp_list$PROP_1_SP1_WILD) +
      simp_list$GRIN_TAXA_PROP_W
  ) / 3
  simp_list$COMP_INDEX_L <- (
    simp_list$SIMP_LAND + (simp_list$ATK_LAND * simp_list$PROP_1_SP1_LAND) +
      simp_list$GRIN_TAXA_PROP_L
  ) / 3
  simp_list$COMP_INDEX_H <- (
    simp_list$SIMP_HYBD + (simp_list$ATK_HYBD * simp_list$PROP_1_SP1_HYBD) +
      simp_list$GRIN_TAXA_PROP_H
  ) / 3
  ################################################################################
  simp_list[(nrow(simp_list) + 1), "COUNTRY"] <- "COLLECTION"
  
  simp_list[(nrow(simp_list)), "SIMP_WILD"] <- TAX_DF$Gini_Simpson[1]
  simp_list[(nrow(simp_list)), "SIMP_LAND"] <- TAX_DF$Gini_Simpson[2]
  simp_list[(nrow(simp_list)), "SIMP_HYBD"]  <- TAX_DF$Gini_Simpson[3]
  
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
  
  simp_list[(nrow(simp_list)), "COMP_INDEX_W"]    <- TAX_DF$Composition_index[1]
  simp_list[(nrow(simp_list)), "COMP_INDEX_L"]    <- TAX_DF$Composition_index[2]
  simp_list[(nrow(simp_list)), "COMP_INDEX_H"]    <- TAX_DF$Composition_index[3]
  
  simp_list[(nrow(simp_list)), "PROP_SPP_WILD"]    <- TAX_DF$propTaxa_collection[1]
  simp_list[(nrow(simp_list)), "PROP_SPP_LAND"]    <- TAX_DF$propTaxa_collection[2]
  simp_list[(nrow(simp_list)), "PROP_SPP_HYBD"]    <- TAX_DF$propTaxa_collection[3]
  
  simp_list[(nrow(simp_list)), "SPP_WILD"]    <- TAX_DF$nTaxa[1]
  simp_list[(nrow(simp_list)), "SPP_LAND"]    <- TAX_DF$nTaxa[2]
  simp_list[(nrow(simp_list)), "SPP_HYBD"]    <- TAX_DF$nTaxa[3]
  
  simp_list[(nrow(simp_list)), "REC_WILD"]    <- TAX_DF$nRecords[1]
  simp_list[(nrow(simp_list)), "REC_LAND"]    <- TAX_DF$nRecords[2]
  simp_list[(nrow(simp_list)), "REC_HYBD"]    <- TAX_DF$nRecords[3]
  
  if(nrow(distw)>0){
    simp_list[(nrow(simp_list)), "ENS_GINI_WILD"]    <- 1/abdiv::dominance(rowSums(df_wild,na.rm = T))  
  } else {
    simp_list[(nrow(simp_list)), "ENS_GINI_WILD"]    <- NA
  }
  
  if(nrow(distl)>0){
    simp_list[(nrow(simp_list)), "ENS_GINI_LAND"]    <- 1/abdiv::dominance(rowSums(df_landrace,na.rm = T))  
  } else {
    simp_list[(nrow(simp_list)), "ENS_GINI_LAND"]    <- NA
  }
  
  if(nrow(disth)>0){
    simp_list[(nrow(simp_list)), "ENS_GINI_HYBD"]    <- 1/abdiv::dominance(rowSums(df_hybrid,na.rm = T))  
  } else {
    simp_list[(nrow(simp_list)), "ENS_GINI_HYBD"]    <- NA
  }
################################################################################
  #https://medium.com/@noorfatimaafzalbutt/a-comprehensive-guide-to-normalization-in-machine-learning-afead759b062
  simp_list$ENS_GINI_WILD_S <- (simp_list$ENS_GINI_WILD - min(simp_list$ENS_GINI_WILD,na.rm = T))/
    (max(simp_list$ENS_GINI_WILD,na.rm = T)-min(simp_list$ENS_GINI_WILD,na.rm = T))
    #simp_list$ENS_GINI_WILD/abs(max(simp_list$ENS_GINI_WILD,na.rm = T))
    #(simp_list$ENS_GINI_WILD - median(simp_list$ENS_GINI_WILD,na.rm = T))/IQR(simp_list$ENS_GINI_WILD,na.rm = T)
  simp_list$ENS_GINI_LAND_S <- (simp_list$ENS_GINI_LAND - min(simp_list$ENS_GINI_LAND,na.rm = T))/
    (max(simp_list$ENS_GINI_LAND,na.rm = T)-min(simp_list$ENS_GINI_LAND,na.rm = T))
    #simp_list$ENS_GINI_LAND/abs(max(simp_list$ENS_GINI_LAND,na.rm = T))
    #(simp_list$ENS_GINI_LAND - median(simp_list$ENS_GINI_LAND,na.rm = T))/IQR(simp_list$ENS_GINI_LAND,na.rm = T)
  simp_list$ENS_GINI_HYBD_S <- (simp_list$ENS_GINI_HYBD - min(simp_list$ENS_GINI_HYBD,na.rm = T))/
    (max(simp_list$ENS_GINI_HYBD,na.rm = T)-min(simp_list$ENS_GINI_HYBD,na.rm = T))
    #simp_list$ENS_GINI_HYBD/abs(max(simp_list$ENS_GINI_HYBD,na.rm = T))
    #(simp_list$ENS_GINI_HYBD - median(simp_list$ENS_GINI_HYBD,na.rm = T))/IQR(simp_list$ENS_GINI_HYBD,na.rm = T)
  ################################################################################
  simp_list$COMP_INDEX_W_S <- (
    simp_list$ENS_GINI_WILD_S + (simp_list$ATK_WILD * simp_list$PROP_1_SP1_WILD) +
      simp_list$GRIN_TAXA_PROP_W
  ) / 3
  simp_list$COMP_INDEX_L_S <- (
    simp_list$ENS_GINI_LAND_S + (simp_list$ATK_LAND * simp_list$PROP_1_SP1_LAND) +
      simp_list$GRIN_TAXA_PROP_L
  ) / 3
  simp_list$COMP_INDEX_H_S <- (
    simp_list$ENS_GINI_HYBD_S+ (simp_list$ATK_HYBD * simp_list$PROP_1_SP1_HYBD) +
      simp_list$GRIN_TAXA_PROP_H
  ) / 3
  ################################################################################
  
  row.names(simp_list)  <- simp_list$COUNTRY
  
  ################################################################################
  simp_list$REGION <- NA
  simp_list$REGION <- countrycode(simp_list$COUNTRY,
                                  origin = 'iso3c',
                                  destination = 'region')
  
  simp_list$REGION[which(is.na(simp_list$REGION))] <- "N/A"
  simp_list$REGION[which(simp_list$COUNTRY == "COLLECTION")] <- "Collection"
  simp_list$REGION[which(simp_list$COUNTRY == "No_country")] <- "Collection"
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
  
  x_res_W <- #simp_list_nona[,c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27)]
    simp_list[, c(#2, #3,
                  5, #6,
                  #8,#9,
                  11, #12,
                  #14,#15,
                  #17, #18,
                  #20,#21,
                  #23, #24,
                  26, #27
                  #29, 
                  #32,
                  35,
                  38,
                  41)]
  
  x_res_W <- x_res_W[complete.cases(x_res_W), ]
  x_res_W$REGION <- factor(x_res_W$REGION)
  
  x_res_L <- #simp_list_nona[,c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27)]
    simp_list[, c(#3, 
                  6, 
                  #9,
                  12, #15,
                  #18, #21,
                  #24, 
                  27, 
                  #30, 
                  #33,
                  36,
                  39,
                  41)]
  x_res_L <- x_res_L[complete.cases(x_res_L), ]
  x_res_L$REGION <- factor(x_res_L$REGION)
  ################################################################################
  if(nrow(x_res_W)>=2){
  x_W <- FactoMineR::PCA(X = x_res_W[, -ncol(x_res_W)],
                         graph = F,
                         scale.unit = T)
  #quali.sup = "REGION",graph = F)
  #[,-ncol(x_res_W)],graph = F)
  REGION <- x_res_W$REGION
  PCA_WILD <- fviz_pca(
    x_W,
    title = paste0(
      "Crop Wild Relatives composition (",
      nrow(x_res_W) - 2,
      " countries)"
    ),
    #addEllipses = TRUE, ellipse.type = "convex",
    
    repel = TRUE#,
    #geom.var = c("text", "point"),
    #col.var = "black"
  ) +
    geom_point(aes(fill = REGION), shape = 21, size = 6) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal()
  
  
  #PCA_WILD
  
  ggplot2::ggsave(
    paste0(
      outdir,
      "/",
      collection_name,
      "/",
      collection_name,
      "_COMP_W.pdf"
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
  if(nrow(x_res_L)>=2){
  x_L <- FactoMineR::PCA(x_res_L[, -ncol(x_res_L)], graph = F)
  
  
  REGION <- x_res_L$REGION
  PCA_LAND <- fviz_pca(
    x_L,
    title = paste0("Landraces composition (", nrow(x_res_L) -
                     2, " countries)"),
    #addEllipses = TRUE, ellipse.type = "convex",
    
    repel = TRUE#,
    #geom.var = c("text", "point"),
    #col.var = "black"
  ) +
    geom_point(aes(fill = REGION), shape = 21, size = 6) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal()
  
  ggplot2::ggsave(
    paste0(
      outdir,
      "/",
      collection_name,
      "/",
      collection_name,
      "_COMP_L.pdf"
    ),
    PCA_LAND,
    dpi = 300,
    units = "in",
    width = 10,
    height = 8
  )
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
################################################################################
dir <- "D:/OneDrive - CGIAR/GERMPLASM_INDEX"
#dir <- "D:/ONEDRIVE/cgiar/OneDrive - CGIAR/GERMPLASM_INDEX"
################################################################################
#outdir
outdir <- paste0(dir, "/BEANS/RESULTS")
################################################################################
collection_name <- "beans"
x1 <- composition_country(outdir, collection_name)
x1$collection <-  collection_name
################################################################################
collection_name <- "cassava"
x2 <- composition_country(outdir, collection_name)
x2$collection <-  collection_name
################################################################################
collection_name <- "forages"
x3 <- composition_country(outdir, collection_name)
x3$collection <-  collection_name
################################################################################