################################################################################
library(countrycode)
library(ggplot2)
library(factoextra)
library(glue)
################################################################################
################################################################################

completeness_func <- function(collection_name,PCDI_df,QGI_df){
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
################################################################################
message("Joining PDCI and GQI to wild, landraces and hybrids subsets")
################################################################################
message("Joining PDCI and GQI to wild")

passport_data_w$PDCI <- NA
passport_data_w$GQI <- NA
for(i in 1:nrow(passport_data_w)){
  #i <- 1
  passport_data_w$PDCI[[i]] <- PCDI_df$PDCI_final[
    which(PCDI_df$ACCENUMB==passport_data_w$ACCENUMB[[i]])]
  passport_data_w$GQI[[i]] <- QGI_df$SCORE[
    which(QGI_df$ACCENUMB==passport_data_w$ACCENUMB[[i]])]
};rm(i)
################################################################################
message("Joining PDCI and GQI to landraces")
passport_data_l$PDCI <- NA
passport_data_l$GQI <- NA
for(i in 1:nrow(passport_data_l)){
  #i <- 1
  passport_data_l$PDCI[[i]] <- PCDI_df$PDCI_final[
    which(PCDI_df$ACCENUMB==passport_data_l$ACCENUMB[[i]])]
  passport_data_l$GQI[[i]] <- QGI_df$SCORE[
    which(QGI_df$ACCENUMB==passport_data_l$ACCENUMB[[i]])]
};rm(i)
################################################################################
message("Joining PDCI and GQI to hybrids")
passport_data_h$PDCI <- NA
passport_data_h$GQI <- NA
for(i in 1:nrow(passport_data_h)){
  #i <- 1
  passport_data_h$PDCI[[i]] <- PCDI_df$PDCI_final[
    which(PCDI_df$ACCENUMB==passport_data_h$ACCENUMB[[i]])]
  passport_data_h$GQI[[i]] <- QGI_df$SCORE[
    which(QGI_df$ACCENUMB==passport_data_h$ACCENUMB[[i]])]
};rm(i)

################################################################################
#spps approach
message("Summarizing results per species")
################################################################################
#wild
spp_w <- unique(passport_data_w$NAME)

completeness_wild_df <- as.data.frame(matrix(nrow = length(spp_w),ncol=7))
colnames(completeness_wild_df) <- c(
  "Taxon","status","records","PDCI_mean","PDCI_sd","GQI_mean","GQI_sd"
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
};rm(i)
################################################################################
#landrace
spp_l <- unique(passport_data_l$NAME)

completeness_land_df <- as.data.frame(matrix(nrow = length(spp_l),ncol=7))
colnames(completeness_land_df) <- c(
  "Taxon","status","records","PDCI_mean","PDCI_sd","GQI_mean","GQI_sd"
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
};rm(i)
################################################################################
#hybrid
spp_h <- unique(passport_data_h$NAME)

completeness_hybd_df <- as.data.frame(matrix(nrow = length(spp_h),ncol=7))
colnames(completeness_hybd_df) <- c(
  "Taxon","status","records","PDCI_mean","PDCI_sd","GQI_mean","GQI_sd"
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
};rm(i)

################################################################################
completeness_wild_df$PDCI_mean <- completeness_wild_df$PDCI_mean/10
completeness_land_df$PDCI_mean <- completeness_land_df$PDCI_mean/10
completeness_hybd_df$PDCI_mean <- completeness_hybd_df$PDCI_mean/10

completeness_wild_df$GQI_mean <- completeness_wild_df$GQI_mean/12
completeness_land_df$GQI_mean <- completeness_land_df$GQI_mean/12
completeness_hybd_df$GQI_mean <- completeness_hybd_df$GQI_mean/12
################################################################################

completeness_df <- rbind(completeness_wild_df,completeness_land_df)
completeness_df <- rbind(completeness_df,completeness_hybd_df)

#################################
#w <- completeness_df$records/sum(completeness_df$records)
PDCI_wild_W <- weighted.mean(x = completeness_wild_df$PDCI_mean,
                                             w = completeness_wild_df$records/
                             sum(completeness_df$records) ,
                                             na.rm = T)
GQI_wild_W <- weighted.mean(x = completeness_wild_df$GQI_mean,
                           w = completeness_wild_df$records/
                             sum(completeness_df$records) ,
                           na.rm = T)

###############################
PDCI_land_W <- weighted.mean(x = completeness_land_df$PDCI_mean,
                             w = completeness_land_df$records/
                               sum(completeness_df$records) ,
                             na.rm = T)
GQI_land_W <- weighted.mean(x = completeness_land_df$GQI_mean,
                            w = completeness_land_df$records/
                              sum(completeness_df$records) ,
                            na.rm = T)
###############################
PDCI_hybd_W <- weighted.mean(x = completeness_hybd_df$PDCI_mean,
                             w = completeness_hybd_df$records/
                               sum(completeness_df$records) ,
                             na.rm = T)
GQI_hybd_W <- weighted.mean(x = completeness_hybd_df$GQI_mean,
                            w = completeness_hybd_df$records/
                              sum(completeness_df$records) ,
                            na.rm = T)
###############################
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
df_completeness_final <- as.data.frame(matrix(nrow=4,ncol=3)
)
colnames(df_completeness_final) <- c("status","PDCI","GQI")
df_completeness_final[,1] <- c("wild","landrace","hybrid","total")
df_completeness_final[1,2] <-  PDCI_wild_W
df_completeness_final[1,3] <-  GQI_wild_W
df_completeness_final[2,2] <-  PDCI_land_W
df_completeness_final[2,3] <-  GQI_land_W
df_completeness_final[3,2] <-  PDCI_hybd_W
df_completeness_final[3,3] <-  GQI_hybd_W
df_completeness_final[4,2] <-  PDCI_total_W
df_completeness_final[4,3] <-  GQI_total_W


df_completeness_final$final_score <- NA
df_completeness_final$final_score <- 
df_completeness_final$PDCI*df_completeness_final$GQI
################################################################################
################################################################################
################################################################################
#Getting countries
countries <- colnames(df_total)
################################################################################
message("Calculating quality per country (wild)")

pdci_country_w <- as.data.frame(matrix(nrow = length(spp_w),ncol=length(countries)))
gqi_country_w <- as.data.frame(matrix(nrow = length(spp_w),ncol=length(countries)))

colnames(pdci_country_w) <- countries
row.names(pdci_country_w) <- spp_w

colnames(gqi_country_w) <- countries
row.names(gqi_country_w) <- spp_w

 for(i in 1:length(countries)){
   #i <- 3
  x_i <- passport_data_w[
    which(passport_data_w$ORIGCTY==countries[[i]]),]
  if(nrow(x_i)==0){
    pdci_country_w[,i] <- NA
    gqi_country_w[,i] <- NA
  } else {
    for(j in 1:length(spp_w)){
      #j <- 1
      x_i_j <- x_i[which(x_i$NAME==spp_w[[j]]),]
      if(nrow(x_i_j)>0){
      pdci_country_w[j,i] <- mean(x_i_j$PDCI,na.rm = T)
      gqi_country_w[j,i] <- mean(x_i_j$GQI,na.rm = T)
      } else {
        pdci_country_w[j,i] <- NA
        gqi_country_w[j,i] <- NA       
      }
    };rm(j,x_i,x_i_j)
  }
};rm(i)
pdci_country_w <- pdci_country_w/10
gqi_country_w <- gqi_country_w/12

################################################################################
message("Calculating quality per country (landraces)")
pdci_country_l <- as.data.frame(matrix(nrow = length(spp_l),ncol=length(countries)))
gqi_country_l <- as.data.frame(matrix(nrow = length(spp_l),ncol=length(countries)))

colnames(pdci_country_l) <- countries
row.names(pdci_country_l) <- spp_l

colnames(gqi_country_l) <- countries
row.names(gqi_country_l) <- spp_l

for(i in 1:length(countries)){
  #i <- 3
  x_i <- passport_data_l[
    which(passport_data_l$ORIGCTY==countries[[i]]),]
  if(nrow(x_i)==0){
    pdci_country_l[,i] <- NA
    gqi_country_l[,i] <- NA
  } else {
    for(j in 1:length(spp_l)){
      #j <- 1
      x_i_j <- x_i[which(x_i$NAME==spp_l[[j]]),]
      if(nrow(x_i_j)>0){
        pdci_country_l[j,i] <- mean(x_i_j$PDCI,na.rm = T)
        gqi_country_l[j,i] <- mean(x_i_j$GQI,na.rm = T)
      } else {
        pdci_country_l[j,i] <- NA
        gqi_country_l[j,i] <- NA       
      }
    };rm(j,x_i,x_i_j)
  }
};rm(i)
pdci_country_l <- pdci_country_l/10
gqi_country_l <- gqi_country_l/12
################################################################################'
message("Calculating quality per country (hybrids)")
pdci_country_h <- as.data.frame(matrix(nrow = length(spp_h),ncol=length(countries)))
gqi_country_h <- as.data.frame(matrix(nrow = length(spp_h),ncol=length(countries)))

colnames(pdci_country_h) <- countries
row.names(pdci_country_h) <- spp_h

colnames(gqi_country_h) <- countries
row.names(gqi_country_h) <- spp_h

for(i in 1:length(countries)){
  #i <- 3
  x_i <- passport_data_h[
    which(passport_data_h$ORIGCTY==countries[[i]]),]
  if(nrow(x_i)==0){
    pdci_country_h[,i] <- NA
    gqi_country_h[,i] <- NA
  } else {
    for(j in 1:length(spp_h)){
      #j <- 1
      x_i_j <- x_i[which(x_i$NAME==spp_h[[j]]),]
      if(nrow(x_i_j)>0){
        pdci_country_h[j,i] <- mean(x_i_j$PDCI,na.rm = T)
        gqi_country_h[j,i] <- mean(x_i_j$GQI,na.rm = T)
      } else {
        pdci_country_h[j,i] <- NA
        gqi_country_h[j,i] <- NA       
      }
    };rm(j,x_i,x_i_j)
  }
};rm(i)
pdci_country_h <- pdci_country_h/10
gqi_country_h <- gqi_country_h/12
################################################################################
completeness_country_final <- as.data.frame(matrix(nrow=length(countries),ncol=7))
colnames(completeness_country_final) <- c("country",
                                          "PDCI_wild","PDCI_landrace","PDCI_hybrid",
                                          "GQI_wild","GQI_landrace","GQI_hybrid"
                                          )
completeness_country_final$country <- countries

#adding PDCI and GQI values 
completeness_country_final$PDCI_wild <- colMeans(pdci_country_w,na.rm = T)
completeness_country_final$PDCI_landrace <- colMeans(pdci_country_l,na.rm = T)
completeness_country_final$PDCI_hybrid <- colMeans(pdci_country_h,na.rm = T)
completeness_country_final$GQI_wild <- colMeans(gqi_country_w,na.rm = T)
completeness_country_final$GQI_landrace <- colMeans(gqi_country_l,na.rm = T)
completeness_country_final$GQI_hybrid <- colMeans(gqi_country_h,na.rm = T)

completeness_country_final[length(countries)+1,1] <- "collection"

#adding collection values
completeness_country_final$PDCI_wild[[length(countries)+1]] <- 
  df_completeness_final$PDCI[which(df_completeness_final$status=="wild")]
completeness_country_final$PDCI_landrace[[length(countries)+1]] <- 
  df_completeness_final$PDCI[which(df_completeness_final$status=="landrace")]
completeness_country_final$PDCI_hybrid[[length(countries)+1]] <- 
  df_completeness_final$PDCI[which(df_completeness_final$status=="hybrid")]
completeness_country_final$GQI_wild[[length(countries)+1]] <- 
  df_completeness_final$GQI[which(df_completeness_final$status=="wild")]
completeness_country_final$GQI_landrace[[length(countries)+1]] <- 
  df_completeness_final$GQI[which(df_completeness_final$status=="landrace")]
completeness_country_final$GQI_hybrid[[length(countries)+1]] <- 
  df_completeness_final$GQI[which(df_completeness_final$status=="hybrid")]

################################################################################
#score
completeness_country_final$SCORE_wild <- 
  completeness_country_final$PDCI_wild * completeness_country_final$GQI_wild

completeness_country_final$SCORE_landrace <- 
  completeness_country_final$PDCI_landrace * completeness_country_final$GQI_landrace

completeness_country_final$SCORE_hybrid <- 
  completeness_country_final$PDCI_hybrid * completeness_country_final$GQI_hybrid
################################################################################
completeness_country_final$REGION <- NA
completeness_country_final$REGION <- countrycode(completeness_country_final$country,
                                origin = 'iso3c',
                                destination = 'region')

completeness_country_final$country[which(completeness_country_final$country=="NA1")] <- "No_country"
completeness_country_final$REGION[which(is.na(completeness_country_final$REGION))] <- "N/A"
completeness_country_final$REGION[which(completeness_country_final$country== "COLLECTION")] <- "Collection"
completeness_country_final$REGION[which(completeness_country_final$country == "No_country")] <- "Collection"
completeness_country_final$REGION[which(completeness_country_final$country == "SCG")] <- "Europe & Central Asia"
completeness_country_final$REGION[which(completeness_country_final$country == "YUG")] <- "Europe & Central Asia"
completeness_country_final$REGION[which(completeness_country_final$country == "ZAR")] <- "Sub-Saharan Africa"
completeness_country_final$REGION[which(completeness_country_final$country == "BUR")] <- "East Asia & Pacific"
###
row.names(completeness_country_final) <- completeness_country_final$country
###

################################################################################
x_res_W <- #simp_list_nona[,c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27)]
  completeness_country_final[, c(
    2,
    5,
    8,
    11
    )]

x_res_W <- x_res_W[complete.cases(x_res_W), ]
x_res_W$REGION <- factor(x_res_W$REGION)

x_res_L <- #simp_list_nona[,c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27)]
  completeness_country_final[, c(
    3,
    6,
    9,
    11
  )]

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
      "Crop Wild Relatives completeness (",
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
      "_completeness_W.pdf"
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
    title = paste0("Landraces completeness (", nrow(x_res_L) -
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
      "_completeness_L.pdf"
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

PCDI_df <- read.csv("D:/OneDrive - CGIAR/GERMPLASM_INDEX/CIAT Data/ciat_pdci.csv")
QGI_df <- read.csv("D:/OneDrive - CGIAR/GERMPLASM_INDEX/CIAT Data/genesys_quality_score_08-2024.csv")

collection_name <- "beans"
#
