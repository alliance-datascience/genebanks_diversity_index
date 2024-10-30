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
conteos_plus_prop <- merge(conteos, dist_prop[,c('Taxa','status','prop_taxa_in_collection_noNA')], by.x = c('SPECIES', 'status'), by.y = c('Taxa', 'status'), all.x = T)



write.csv(conteos_plus_prop, paste0(outdir, '/beans/ecoregions.csv'), row.names = F )
write.csv(fnl, paste0(outdir, '/ecoregions_high_medium_quality_bean_col003.csv'), row.names = F )
write.csv(qs, paste0(outdir, '/percentage_quality_bean_col003.csv'), row.names = F )


################################################################################################3

##loading native areas
message(paste0("Loading: ",outdir, "/", collection_name,"/" ,collection_name, "_native_iso3_1.csv"))

native_list <- read.csv(paste0(outdir, "/", collection_name,"/" ,collection_name,  "_native_iso3_new_1.csv"))
x<- native_list[,c('taxa', 'ISO3_n')] %>% group_by(taxa) %>% mutate(native = paste0(ISO3_n, collapse = ", ")) 
write.csv(x, paste0(outdir, '/native_c_beans.csv'), row.names = F )

species_list <- paste0(fnl$GENUS, ' ',fnl$SPECIES); species_list<-unique(species_list)
sampstat_list<-c(100, 120, 300)


df<-data.frame(species = NA,sampstat = NA,n_ecor = NA, native_countries = NA, ecor_in_native =  NA)
i = 0

for(species in species_list){
  
  for(sampstat in sampstat_list){
    
    i = i+1
    
    native_species <- native_list$ISO3_n[which(native_list$taxa == species)]
    
    data_species_sampstat<- fnl[which(paste0(fnl$GENUS, ' ',fnl$SPECIES) == species & fnl$SAMPSTAT == sampstat), ]
    ecor_native_countries<- data_species_sampstat$ECO_NAME[which(data_species_sampstat$GADM_GID_0 %in% native_species )]
    ecor_native_countries<-unique(ecor_native_countries)
    
    df[i,]<- data.frame(species = species,sampstat = sampstat,n_ecor = length(unique(data_species_sampstat$ECO_NAME[!is.na(data_species_sampstat$ECO_NAME)])), native_countries = paste0(native_species, collapse = ", "), ecor_in_native =  length(ecor_native_countries))
    
    
  }
  
}

df<- df[-which(df$n_ecor == 0),]
df$p_ecor_in_native<- df$ecor_in_native/df$n_ecor

write.csv(df, paste0(outdir, '/beans/ecoregions_counts_native_countries.csv'), row.names = F )


########################################################

species_list <- paste0(fnl$GENUS, ' ',fnl$SPECIES); species_list<-unique(species_list)
sampstat_list<-c(100, 120, 300)
country_list <-unique(fnl$GADM_GID_0); country = country_list[1]



df<-data.frame(country  = NA, is_native = NA, n_ecor_country = NA, species = NA,sampstat = NA, n_ecor_sp_in_country = NA )
i = 0

for(country in country_list){
  
  for(species in species_list){
    
    for(sampstat in sampstat_list){
      i= i+1
      
      native_species <- native_list$ISO3_n[which(native_list$taxa == species)]
      
      is_native <- ifelse(country %in% native_species, TRUE, FALSE)
      
      ecor_country<- fnl$ECO_NAME[which(fnl$GADM_GID_0 %in% country )]; n_ecor_country<-length(unique(ecor_country[!is.na(ecor_country)]))
      
      data_country<- fnl[which(fnl$GADM_GID_0 %in% country ),]
      n_ecor_sp_in_country<- data_country$ECO_NAME[which(paste0(data_country$GENUS, ' ',data_country$SPECIES) == species & data_country$SAMPSTAT == sampstat)]
      n_ecor_sp_in_country<- length(unique(n_ecor_sp_in_country))
      
      
      df[i,]<- data.frame(country  = country, is_native = is_native, n_ecor_country = n_ecor_country, species = species,sampstat = sampstat, n_ecor_sp_in_country = n_ecor_sp_in_country )
      
    }
    
  }
  
  
  
  
}

df$p_ecor_sp_in_country <- df$n_ecor_sp_in_country/df$n_ecor_country
write.csv(df, paste0(outdir, '/beans/ecoregions_counts_per_countries.csv'), row.names = F )


##################################################################################33


countries<- unique(df$country)
sampstat_list<-c(100, 120, 300)
df_country <- data.frame(country = NA,sampstat = NA ,p_ecor_sp = NA )
i = 0

for(country in countries){
  for(sampstat in sampstat_list){
    
    i = i+1
    
    
    col<-df[df$country == country,]
    col_samp<-col[col$sampstat == sampstat,]
    col_samp<- col_samp[-which(grepl('X|x', col_samp$species)),]
    col_samp<- col_samp[-which(col_samp$n_ecor_sp_in_country == 0),]
    
  
    if(nrow(col_samp) == 0){
      
      df_country[i,] <- data.frame(country = country,sampstat = sampstat ,p_ecor_sp = paste0('Not sampstat = ', sampstat, ' species present in the country'))
      
    } else{
      
      df_country[i,] <- data.frame(country = country,sampstat = sampstat ,p_ecor_sp = mean(col_samp$p_ecor_sp_in_country))
      
    }
    
  }
}


write.csv(df_country, paste0(outdir, '/beans/p_ecoregions_per_countries_per_species.csv'), row.names = F )






































#loading summary table
dist_prop <- read.csv(
  paste0(outdir, "/", collection_name,"/" ,collection_name,  "_summary_table.csv"))
message(paste0("Loaded: ",outdir, "/", collection_name,"/" ,collection_name,  "_summary_table.csv"))

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
  