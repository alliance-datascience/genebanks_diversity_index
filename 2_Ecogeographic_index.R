require(data.table)
require(readr)
require(stringr)
require(tidyr)
require(psych)
install.packages("ineq")
library(ineq)
install.packages("igraph")
require(igraph)
library(abdiv)
install.packages("vegan")
library(vegan)
library(matrixStats)
require(factoextra)
require(ecodist)
require(DescTools)
require(rgdal)
require(sp)
require(sf)
require(terra)


################################################################################
################################################################################
dir <- "D:/OneDrive - CGIAR/GERMPLASM_INDEX"
#dir <-  "D:/OneDrive - CGIAR/GERMPLASM_INDEX"
################################################################################
#outdir
outdir <- paste0(dir, "/BEANS/RESULTS")
################################################################################
collection_name="beans"

passport_data_original <- read.csv("D:/OneDrive - CGIAR/Genebanks/data/genesys_quality_score_08-2024.csv", na.strings = NA)
passport_data_original <- passport_data_original[which(passport_data_original$CROPCODE ==
                                                     collection_name & passport_data_original$INSTCODE == 'COL003'), ]

passport_data_orig<- passport_data_original[which(passport_data_original$quality_score %in% c('High', 'Moderate')),]


###################################################################################

passport_data_orig_complete<- passport_data_orig[complete.cases(passport_data_orig[,c('DECLONGITUDE', 'DECLATITUDE')]),]
qs<-passport_data_original[,c('SPECIES', 'quality_score','SAMPSTAT')] %>% group_by(SPECIES,SAMPSTAT, quality_score) %>% summarise(n = n())




p = vect(passport_data_orig_complete, geom = c("DECLONGITUDE", "DECLATITUDE"), crs = "EPSG:4326")
v = vect("D:/OneDrive - CGIAR/GERMPLASM_INDEX/CIAT Data/ecoregr/wwf_terr_ecos.shp")

vls <-terra::extract(v, p)
fnl <- cbind(as.data.frame(p), vls[,c('ECO_NAME', 'ECO_ID', 'BIOME')])





conteos<-fnl[,c('SAMPSTAT','SPECIES', 'ECO_NAME', 'GADM_GID_0')] %>% group_by(SAMPSTAT, SPECIES,  GADM_GID_0 ,ECO_NAME) %>% summarise(n = n())
x<-v['ECO_NAME']
length(unique(x$ECO_NAME))
length(unique(conteos$ECO_NAME))
conteos$status<-ifelse(conteos$SAMPSTAT %in% c(100, 120), 'wild', ifelse(conteos$SAMPSTAT == 300, 'landrace', NA))
conteos$SPECIES <- paste0('Phaseolus ', conteos$SPECIES)


#loading summary table
dist_prop <- read.csv(
  paste0(outdir, "/", collection_name,"/" ,collection_name,  "_summary_table.csv"))
message(paste0("Loaded: ",outdir, "/", collection_name,"/" ,collection_name,  "_summary_table.csv"))


conteos_plus_prop <- merge(conteos, dist_prop[,c('Taxa','status','prop_taxa_in_collection_noNA')], by.x = c('SPECIES', 'status'), by.y = c('Taxa', 'status'), all.x = T)
write.csv(conteos_plus_prop, paste0(outdir, '/beans/ecoregions.csv'), row.names = F )
write.csv(fnl, paste0(outdir, '/ecoregions_high_medium_quality_bean_col003.csv'), row.names = F )
write.csv(qs, paste0(outdir, '/percentage_quality_bean_col003.csv'), row.names = F )

###################################################################################

gadm<-vect("//alliancedfs.alliance.cgiar.org/gap_analysis_landraces/runs/input_data/shapefiles/gadm36_shp/gadm36.shp")

f<- terra::merge(gadm, passport_data_orig_complete[,c('DECLONGITUDE', 'DECLATITUDE')])

################################################################################


# ecogeographic_index_function <- function(passport_data_orig,
#                                          outdir,
#                                          dist_prop,
#                                          
#                                          collection_name){
  
  ##loading native areas
  message(paste0("Loading: ",outdir, "/", collection_name,"/" ,collection_name, "_native_iso3_1.csv"))
  
  native_list <- read.csv(paste0(outdir, "/", collection_name,"/" ,collection_name,  "_native_iso3_new_1.csv"))
  x<- native_list[,c('taxa', 'ISO3_n')] %>% group_by(taxa) %>% mutate(native = paste0(ISO3_n, collapse = ", ")) 
  write.csv(x, paste0(outdir, '/native_c_beans.csv'), row.names = F )
  
  #loading count matrix  
  df_total2 <- read.csv(paste0(outdir, "/", collection_name,"/" ,collection_name, "_count_matrix_countries_1.csv"),row.names = 1)
  message(paste0("Loading:",outdir, "/", collection_name,"/" ,collection_name,"_count_matrix_countries_1.csv"))
  
  
  #loading summary table
  dist_prop <- read.csv(
            paste0(outdir, "/", collection_name,"/" ,collection_name,  "_summary_table.csv"))
  message(paste0("Loaded: ",outdir, "/", collection_name,"/" ,collection_name,  "_summary_table.csv"))
  
  
  #loading Aichi results 
  aichi<- read.csv(paste0(dir, "/GRIN_GBIF_EQUI/Species_Gap_Analysis_Wild.csv"))
  message(paste0("Loaded: ",dir, "/GRIN_GBIF_EQUI/Species_Gap_Analysis_Wild.csv"))
  
  aichi<- aichi[which(aichi$current_GRIN_id %in% dist_prop$GRIN_ID ),]
  
  #Calculating weights per coll 
  aichi$w1<- (aichi$ERS_ex)/100
  aichi$w2<-  (aichi$ERS_in)/100
  aichi$w3<-  (aichi$w1 + aichi$w2)/2
  aichi$w4<- (aichi$FCS_ex)/100
    
  
  dist_prop_1<-merge(dist_prop, aichi[,c(2,44:47)], by.x = 'GRIN_ID', by.y = 'current_GRIN_id', all.x = T)
  
  dist_prop_1<- dist_prop_1[which(dist_prop_1$status == 'wild'),]
  dist_prop_1<- dist_prop_1[!is.na(dist_prop_1$w1),]
  
  
  dist_prop_1$w1[is.na(dist_prop_1$w1)]<-0
  dist_prop_1$w2[is.na(dist_prop_1$w2)]<-0
  dist_prop_1$w3[is.na(dist_prop_1$w3)]<-0
  dist_prop_1$w4[is.na(dist_prop_1$w4)]<-0
  
  dist_prop_1$index_1_1 <- dist_prop_1$nCounnat_countotal * dist_prop_1$w1
  dist_prop_1$index_1_2 <- dist_prop_1$nCounnat_countotal  * dist_prop_1$w2
  dist_prop_1$index_1_3 <- dist_prop_1$nCounnat_countotal * dist_prop_1$w3
  dist_prop_1$index_1_4 <- dist_prop_1$nCounnat_countotal * dist_prop_1$w4
  
  dist_prop_1$index_2_1 <- dist_prop_1$nCounnonnat_countotal * dist_prop_1$w1
  dist_prop_1$index_2_2 <- dist_prop_1$nCounnonnat_countotal * dist_prop_1$w2
  dist_prop_1$index_2_3 <- dist_prop_1$nCounnonnat_countotal* dist_prop_1$w3
  dist_prop_1$index_2_4 <-dist_prop_1$nCounnonnat_countotal * dist_prop_1$w4
  
  
  
  
  #loading Landrace results 
  landraces<- read.csv('C:/Users/mvdiaz/Downloads/41477_2022_1144_MOESM6_ESM.csv')
  
  landraces$Scientific_name <- c('Hordeum vulgare', 'Eleusine coracana', 'Zea mays', 'Pennisetum glaucum', 'Oryza glaberrima', 'Oryza sativa',
                                 'Sorghum bicolor', 'Triticum aestivum', 'Triticum durum', 'Cicer arietinum', 'Phaseolus vulgaris', 'Vigna unguiculata',
                                 'Vicia faba', 'Lathyrus sativus', 'Arachis hypogaea', 'Lens culinaris', 'Pisum sativum', 'Cajanus cajan', 'Musa sp.', 'Artocarpus altilis',
                                 'Manihot esculenta', 'Solanum tuberosum', 'Ipomoea batatas', 'Colocasia esculenta', 'Dioscorea alata')

  dist_prop_2<- merge(dist_prop, landraces[,c(7,4)], by.x = 'Taxa', by.y = 'Scientific_name', all.x = T)
  
  dist_prop_2<- dist_prop_2[which(dist_prop_2$status == 'landrace'),]
  dist_prop_2<- dist_prop_2[!is.na(dist_prop_2[,ncol(dist_prop_2)]),]
  
  
  
  names(dist_prop_2)[ncol(dist_prop_2)]<- 'index_1_land'
  
  dist_prop_2$index_1_land<-  ((100-dist_prop_2$index_1_land)/100) * dist_prop_2$nat_nonnat
  
  dist_prop_2$index_1_land[is.na(dist_prop_2$index_1_land)]<-0
  


  
  
  
  ####################
  
#}
  #SHANNON
  # N <- sum(dist_prop$ncounts_Country[which(dist_prop$status=="wild")])
  # SHANNON <- -sum((dist_prop$ncounts_Country[which(dist_prop$status=="wild")]/N)*
  #  log((dist_prop$ncounts_Country[which(dist_prop$status=="wild")]/N)))
  #   SHANNON_E <-  SHANNON/
  #   log(length(dist_prop$ncounts_Country[which(dist_prop$status=="wild")]))
  
  
  #(1-DescTools::Gini(dist_prop$ncounts_Country,unbiased = T))
  #abdiv::pielou_e(dist_prop$ncounts_Country)
  ##############################################################################  
  
  #   ################################################################################
  #   message("PART 2: CREATING COUNT TABLE WITH TAXA AND COLLECTION SITE")
  #   
  #   ################################################################################
  #   
  #   ################################################################################
  #   #converting all NA COLLSITE to NA1 to avoid issues in counts
  #   message("converting all NA COLLSITE to NA1 to avoid issues in counts")
  #   passport_data_orig$COLLSITE[which(is.na(passport_data_orig$COLLSITE))] <- "NA1"
  #   passport_data_w$COLLSITE[which(is.na(passport_data_w$COLLSITE))] <- "NA1"
  #   passport_data_l$COLLSITE[which(is.na(passport_data_l$COLLSITE))] <- "NA1"
  #   
  #   coll_site_total <- unique(passport_data_orig$COLLSITE)
  #   spp <- unique((passport_data_w$NAME))
  #   #####################################
  #   #Wild species
  #   message("CREATING COUNT TABLE FOR COLLECTING SITES AND WILD TAXA")
  #   if(nrow(passport_data_w)>0){
  #   #coll_site_w<- unique(passport_data_w$COLLSITE)
  #   df_coll_w <- data.frame(matrix(nrow = length(spp), ncol = length(coll_site_total)))
  #   colnames(df_coll_w) <- coll_site_total
  #   row.names(df_coll_w) <- spp
  #   for (i in 1:length(spp)) {
  #     for (j in 1:length(coll_site_total)) {
  #       x <- passport_data_w[which(passport_data_w$NAME == spp[[i]] &
  #                                    passport_data_w$COLLSITE == coll_site_total[[j]]), ]
  #       
  #       df_coll_w[i, j] <- nrow(x)
  #       rm(x)
  #     }
  #     rm(j)
  #   }
  #   rm(i)
  #   
  #   #adding wild label
  #   row.names(df_coll_w) <- paste0(row.names(df_coll_w), "_", "wild")
  #   #rm(x)
  #   } else {
  #     df_coll_w <- NA
  #   }
  #   ################################################################################
  #   #####################################
  #   message("CREATING COUNT TABLE FOR COLLECTING SITES AND LANDRACE TAXA")
  #   
  #   #Landraces
  #   if(nrow(passport_data_l)>0){
  #   spp <- unique((passport_data_l$NAME))
  #   #coll_site_w<- unique(passport_data_w$COLLSITE)
  #   df_coll_l <- data.frame(matrix(nrow = length(spp), ncol = length(coll_site_total)))
  #   colnames(df_coll_l) <- coll_site_total
  #   row.names(df_coll_l) <- spp
  #   for (i in 1:length(spp)) {
  #     for (j in 1:length(coll_site_total)) {
  #       x <- passport_data_l[which(passport_data_l$NAME == spp[[i]] &
  #                                    passport_data_l$COLLSITE == coll_site_total[[j]]), ]
  #       
  #       df_coll_l[i, j] <- nrow(x)
  #       rm(x)
  #     }
  #     rm(j)
  #   }
  #   rm(i)
  #   
  #   #Adding landrace label
  #   row.names(df_coll_l) <- paste0(row.names(df_coll_l), "_", "landrace")
  #   } else {
  #     df_coll_l <- NA
  #   }
  #   ################################################################################
  #   #joining data
  #   df_coll_total <- rbind(df_coll_w, df_coll_l)
  #   #validating counts
  #   coll_un <- colSums(df_coll_total)
  #   coll_un <- data.frame(collsite = names(coll_un), freq = as.numeric(coll_un))
  #   coll_un <- coll_un[order(coll_un$freq, decreasing = T), ]
  #   # sum(coll_un$freq)
  #   # sum(df_coll_w)
  #   # sum(df_coll_l)
  #   # sum(df_coll_total)
  #   ################################################################################
  #   write.csv(df_coll_total,
  #             paste0(outdir, "/", collection_name, "_count_matrix_collsites.csv"),row.names = F,na = "")
  #   message(paste0("Saved:",outdir, "/", collection_name, "_count_matrix_collsites.csv"))
  #   
  #   ################################################################################
  #   #obtaining distances for countries
  #   
  #   message("CREATING SUMMATY TABLE FOR COLLECTING SITES AND TAXA")
  #   
  #   #removing records without countries data
  #   x <- df_coll_total[, !colnames(df_coll_total) %in% "NA1"]
  #   #ommiting countries with 0 counts
  #   x <- x[, colSums(x,na.rm = T) > 0]
  #   #ommiting taxa with 0 counts
  #   
  #   x <- x[rowSums(x,na.rm = T) > 0, ]
  #   #creating presence-absence matrix
  #   x2 <- x
  #   x2[x2 > 0] <- 1
  #   
  #   
  #   ################################################################################
  #   #creating matrix to analyze countries uniqueness per taxon
  #   dist_prop_coll <- as.data.frame(matrix(nrow = nrow(x2), ncol = 7))
  #   colnames(dist_prop_coll) <- c(
  #     "Sp",
  #     "sites/sites_total-sp",
  #     "sites/sites_total",
  #     "ncollsites",
  #     "ncollsites_total",
  #     #"ncollsites_comp",
  #     "ncounts",
  #     "nprop"
  #   )#row.names(x2)
  #   row.names(dist_prop_coll) <- row.names(x2)
  #   
  #   #for each species
  #   for (i in 1:nrow(x2)) {
  #     #i <- 54
  #     #splitting each species
  #     a <- x2[i, ]
  #     #splitting for the rest of counts which are not for the taxon i
  #     b <- x2[-i, ]
  #     
  #     #creating a dummy data frame
  #     df <- data.frame(
  #       country = colnames(x2),
  #       a = colSums(a),
  #       b = colSums(b),
  #       c = colSums(x2)
  #     )
  #     # changing sums to presence and obtain countries collecting
  #     df$b[df$b > 0] <- 1
  #     df$c[df$c > 0] <- 1
  #     #taxon name
  #     dist_prop_coll[i, 1] <- row.names(a)
  #     #countries of taxon i / countries where !=i
  #     dist_prop_coll[i, 2] <- sum(df$a) / sum(df$b)
  #     #countries of taxon i / all countries
  #     dist_prop_coll[i, 3] <- sum(df$a) / sum(df$c)
  #     #taxon i countries
  #     dist_prop_coll[i, 4] <- sum(df$a)
  #     #all countries
  #     dist_prop_coll[i, 5] <- sum(df$c)
  #     #numbers of records for taxon i
  #     dist_prop_coll[i, 6] <- sum(x[i, ])
  #     #proportion of records for taxon i/ records of collection x 100
  #     dist_prop_coll[i, 7] <- (sum(x[i, ]) / sum(x)) * 100
  #     rm(df, a, b)
  #   }
  #   rm(i)
  #   
  #   #obtaining taxon status (wild or landrace)
  #   stat <- strsplit(row.names(dist_prop_coll), "_")
  #   #status
  #   dist_prop_coll$status <- unlist(lapply(stat, `[[`, 2))
  #   #taxon name
  #   dist_prop_coll$taxon <- unlist(lapply(stat, `[[`, 1))
  #   ################################################################################
  #   
  #   write.csv(dist_prop_coll,
  #             paste0(outdir, "/", collection_name, "_diss_collsite.csv"),row.names = F,na = "")
  #   message(paste0("Saved:",outdir, "/", collection_name, "_diss_collsite.csv"))
  #   
  #   ################################################################################
  #   #taxonomy index is 1- Gini (inequality) * 1- alpha (diversity) * Piellou (eveness)
  #   
  #   message("CALCULATING GINI SIMPSON FOR COLLECITNG SITES")
  #   
  #   "
  # if the index is high then the collection has a high diversity, good balance of records,
  # and good eveness (similar number of records)
  # "
  #   #taxonomy index all
  #   I1C_all <- (1 - ineq::ineq(dist_prop_coll$ncountries, type = "Gini")) *
  #     abdiv::simpson(dist_prop_coll$ncountries) #*
  #   #abdiv::pielou_e(dist_prop_coll$ncountries)
  #   #taxonomy index crop wild
  #   I1C_W <- (1 - ineq::ineq(dist_prop_coll$ncountries[which(dist_prop_coll$status ==
  #                                                              "wild")], type = "Gini")) *
  #     abdiv::simpson(dist_prop_coll$ncountries[which(dist_prop_coll$status ==
  #                                                      "wild")]) #*
  #   #abdiv::pielou_e(dist_prop_coll$ncountries[which(dist_prop_coll$status=="wild")])
  #   #taxonomy index landrace
  #   I1C_L  <- (1 - ineq::ineq(dist_prop_coll$ncountries[which(dist_prop_coll$status ==
  #                                                               "landrace")], type = "Gini")) *
  #     abdiv::simpson(dist_prop_coll$ncountries[which(dist_prop_coll$status ==
  #                                                      "landrace")]) #*
  #   # abdiv::pielou_e(dist_prop_coll$ncountries[which(dist_prop_coll$status=="landrace")])
  #   ################################################################################
  #   #taxonomy index no subsets
  #   
  #   INDEX_1 <- (1 - ineq::ineq(rowSums(df_total), type = "Gini")) *
  #     abdiv::simpson(rowSums(df_total)) #*
  #   #abdiv::pielou_e(dist_prop_coll$ncountries)
  #   #taxonomy index crop wild
  #   INDEX_W <- (1 - ineq::ineq(rowSums(df_wild), type = "Gini")) *
  #     abdiv::simpson(rowSums(df_wild)) #*
  #   #abdiv::pielou_e(dist_prop_coll$ncountries[which(dist_prop_coll$status=="wild")])
  #   #taxonomy index landrace
  #   INDEX_L  <- (1 - ineq::ineq(rowSums(df_landrace), type = "Gini")) *
  #     abdiv::simpson(rowSums(df_landrace)) #*
  #   # abdiv::pielou_e(dist_prop_coll$ncountries[which(dist_prop_coll$status=="landrace")])
  #   
  #   ################################################################################
  #   #ALL
  #   message("PART 3: observing passport based geographic indexes")
  # 
  #   nat_nonnat  <- mean(dist_prop$nat_nonnat,na.rm = T)
  #   nCounnat_countotal <- mean(dist_prop$nCounnat_countotal,na.rm = T)
  #   nCounnonnat_countotal <- mean(dist_prop$nCounnonnat_countotal,na.rm = T)
  #   ################################################################################
  #   message("PART 4: CREATING COLLECTION INDEXES SUMMATY TABLE!")
  #   
  #   ################################################################################
  #   x_Gini <- data.frame(
  #     COLLECTION = collection_name,
  #     #GRIN
  #     NTAXA = length(accepted_grin_tax_list$grin_id_cur),
  #     NTAXA_GRIN = length(na.omit(accepted_grin_tax_list$grin_id_cur)) /
  #       length(accepted_grin_tax_list$grin_id_cur),
  #     PROP_GRIN_WILD = length(na.omit(WILD_SPP$grin_id_cur)) / nrow(WILD_SPP),
  #     PROP_GRIN_LAND = length(na.omit(LAND_SPP$grin_id_cur)) / nrow(LAND_SPP),
  #     NTAXA_WILD = length(row.names(df_wild[rowSums(df_wild) > 0])),
  #     NTAXA_LAND = length(row.names(df_landrace[rowSums(df_landrace) >
  #                                                 0])),
  #     PROP_NTAXA_WILD = length(row.names(df_wild[rowSums(df_wild) >
  #                                                  0])) / length(accepted_grin_tax_list$grin_id_cur),
  #     PROP_NTAXA_LANDRACE = length(row.names(df_landrace[rowSums(df_landrace) >
  #                                                          0])) / length(accepted_grin_tax_list$grin_id_cur),
  #     NRECORDS_WILD = sum(rowSums(df_wild)),
  #     NRECORDS_LAND = sum(rowSums(df_landrace)),
  #     PROP_NR_WILD = sum(rowSums(df_wild)) / sum(rowSums(df_total)),
  #     PROP_NR_LAND = sum(rowSums(df_landrace)) / sum(rowSums(df_total)),
  #     #ALL
  #     ALL = INDEX_1,
  #     ALL_GINI = (1 - ineq::ineq(rowSums(df_total), type = "Gini")),
  #     ALL_SIMP = abdiv::simpson(rowSums(df_total)),
  #     ALL_PIE = abdiv::pielou_e(rowSums(df_total)),
  #     #WILD
  #     WILD = INDEX_W,
  #     WILD_GINI = (1 - ineq::ineq(rowSums(df_wild), type = "Gini")),
  #     WILD_SIMP = abdiv::simpson(rowSums(df_wild)),
  #     WILD_PIE = abdiv::pielou_e(rowSums(df_wild)),
  #     #LANDRACE
  #     LAND = INDEX_L,
  #     LAND_GINI = (1 - ineq::ineq(rowSums(df_landrace), type = "Gini")),
  #     LAND_SIMP = abdiv::simpson(rowSums(df_landrace)),
  #     LAND_PIE = abdiv::pielou_e(rowSums(df_landrace)),
  #     #COUNTRIES
  #     COUNTRY_all = I1_all,
  #     COUNTRY_W = I1_W,
  #     COUNTRY_L = I1_L,
  #     COLL_all = I1C_all,
  #     COLL_W = I1C_W,
  #     COLL_L = I1C_L,
  #     nat_nonnat=nat_nonnat,
  #     nCounnat_countotal=nCounnat_countotal,
  #     nCounnonnat_countotal=nCounnonnat_countotal
  #   )
  #   
  #   message("saved: indexes.csv")
  #   
  #   write.csv(x_Gini, paste0(outdir, "/", collection_name, "_indexes.csv"),row.names = F,na = "")
  #   #
  
  
  ##############################################################################
  ##############################################################################
  
  # #bray  <- vegdist(x, "bray",na.rm = T,binary = F)
  # # x2<- x2[,colSums(x2)>1]
  # # x2 <- x2[rowSums(x2)>0,]
  #  jac_x  <- vegdist(x2, "jaccard",na.rm = T,binary = F)
  # # jac_x  <- (cor(t(x)))
  # # jac_x[jac_x<0] <- NA
  # # jac_x  <-  1- jac_x
  # #
  # jac_x2 <- as.matrix(jac_x)
  #  diag(jac_x2) <- NA
  #  hh <-hclust(bray,method="ward.D")
  #  plot(hh)
  # #
  # jacc_dist_metric <- data.frame(taxon= rownames(jac_x2),Median=rowMedians(jac_x2,na.rm = T))
  #
  #
  #  x_trans <- x#/colSums(x,na.rm = T)
  #
  # for(i in 1:ncol(x)){
  #   x_trans[,i] <- x_trans[,i]/sum(x_trans[,i],na.rm = T)
  # };rm(i)
  
  
  #euc_x  <- vegdist(t(x_trans), "bray",na.rm = T,binary = F)
  #euc_x2 <- as.matrix(euc_x)
  
  
  
  #x2 <- proxy::dist(x2, by_rows = TRUE, method = "Jaccard")
  
