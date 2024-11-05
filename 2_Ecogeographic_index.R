require(data.table)
require(readr)
require(stringr)
require(tidyr)
require(psych)
#install.packages("ineq")
library(ineq)
#install.packages("igraph")
require(igraph)
library(abdiv)
#install.packages("vegan")
library(vegan)
library(matrixStats)
require(factoextra)
require(ecodist)
require(DescTools)
require(rgdal)
require(sp)
require(sf)
require(terra)
require(dplyr)


################################################################################
################################################################################
dir <- "D:/OneDrive - CGIAR/GERMPLASM_INDEX"
################################################################################
#outdir
outdir <- paste0(dir, "/BEANS/RESULTS")
################################################################################
collection_name="forages" #beans, forages, cassava


passport_data_original <- read.csv("D:/OneDrive - CGIAR/Genebanks/data/genesys_quality_score_09-2024.csv", na.strings = NA)
passport_data_original <- passport_data_original[which(passport_data_original$CROPCODE ==
                                                     collection_name & passport_data_original$INSTCODE == 'COL003'), ]

passport_data_orig<- passport_data_original[which(passport_data_original$quality_score_fixed %in% c('High', 'Moderate')),]


###################################################################################

passport_data_orig_complete<- passport_data_orig[complete.cases(passport_data_orig[,c('DECLONGITUDE', 'DECLATITUDE')]),]
qs<-passport_data_original[,c('SPECIES', 'quality_score','SAMPSTAT')] %>% group_by(SPECIES,SAMPSTAT, quality_score) %>% summarise(n = n())


p = vect(passport_data_orig_complete, geom = c("DECLONGITUDE", "DECLATITUDE"), crs = "EPSG:4326")


if(!exists('v')){
  
  v<- vect("D:/OneDrive - CGIAR/GERMPLASM_INDEX/CIAT Data/ecoregr/wwf_terr_ecos.shp")
}
 

vls <-terra::extract(v, p)
fnl <- cbind(as.data.frame(p), vls[,c('ECO_NAME', 'ECO_ID', 'BIOME')])

fnl$sampstat_cat <-ifelse(fnl$SAMPSTAT < 200 , 'Wild', ifelse(fnl$SAMPSTAT == 300, 'Landrace', 'Other'))

if(length(which(fnl$sampstat_cat == 'Other'))!=0){
  
  fnl<-fnl[-which(fnl$sampstat_cat == 'Other'),]
  
  
}
fnl$name_sc<- paste0(fnl$GENUS, " ", fnl$SPECIES)
  
conteos<-fnl[,c('name_sc','sampstat_cat' ,'ECO_NAME', 'GADM_GID_0')] %>% group_by(name_sc, sampstat_cat,  GADM_GID_0 ,ECO_NAME) %>% summarise(n = n())


write.csv(conteos, paste0(outdir, '/',collection_name,'/Eco/ecoregions_summary.csv'), row.names = F )
write.csv(fnl[,c("ACCENUMB","GENUS", "SPECIES" ,"SPAUTHOR" ,"SUBTAXA", "SUBTAUTHOR" ,"GRIN_AUTHOR","SAMPSTAT", "ORIGCTY", "COLLSITE" ,"ELEVATION", "quality_score_fixed", "ECO_NAME", "ECO_ID", "BIOME", "sampstat_cat", "name_sc"  )], 
          paste0(outdir, '/',collection_name,'/Eco/ecoregions_high_medium_quality_col003.csv'), row.names = F )
write.csv(qs, paste0(outdir,'/',collection_name, '/Eco/percentage_quality_col003.csv'), row.names = F )


################################################################################################3

#ecore_world<-st_read("C:/Users/mvdiaz/Downloads/ECO_COUN.shp")

# Obtener las ecorregiones por país
#ecorregiones_por_pais <- ecore_world %>%
#  select(COUNTRY, ECO_NAME) %>%  # admin es la columna de país en world; ECO_NAME es la columna de nombre de ecorregión
#  group_by(COUNTRY) 

#ecorregiones_por_pais<-ecorregiones_por_pais %>% st_drop_geometry()  

# Descargar los límites de los países
#world <- ne_countries(scale = "medium", returnclass = "sf")

# Create a dataframe with country names and ISO3 codes
#country_iso3_df <- world %>%
#  select(country = admin, iso3 = iso_a3) %>%  # Select columns for country name and ISO3 code
#  st_drop_geometry()  # Remove geometry to keep it as a regular dataframe

#ecorregiones_por_pais<-merge(ecorregiones_por_pais, country_iso3_df, by.x = "COUNTRY", by.y = "country")

#rm(ecore_world)
#rm(world)

#write.csv(ecorregiones_por_pais, "C:/Users/mvdiaz/Downloads/ecor_por_pais.csv", row.names = F)
ecorregiones_por_pais<- read.csv("C:/Users/mvdiaz/Downloads/ecor_por_pais.csv")
################################################################################################3

##loading native areas

native_list <- read.csv(paste0(outdir, "/", collection_name,"/" ,collection_name,  "_native_iso3_new_1.csv"))
x<- native_list[,c('taxa', 'ISO3_n')] %>% group_by(taxa) %>% mutate(native = paste0(ISO3_n, collapse = ", ")) 
x<- x[,-2]
x<-x[!duplicated(x),]
write.csv(x, paste0(outdir, '/',collection_name,'/native_countries_summary.csv'), row.names = F )

##############################################################################################3

species_list<-unique(fnl$name_sc)
sampstat_list<-unique(fnl$sampstat_cat)

if(length(which(grepl('X|x', species_list))) != 0){
  
  species_list<- species_list[-which(grepl('X|x', species_list))]
  
}else{
  
  species_list<-species_list
  
}

#species= species_list[1]
#sampstat = sampstat_list[1]

df<-data.frame(species = NA,sampstat = NA, native_countries = NA, n_accesions = NA, n_ecor = NA,n_ecor_in_native =  NA)
i = 0

for(species in species_list){
  
  for(sampstat in sampstat_list){
    
    i = i+1
    
    
    native_species <- native_list$ISO3_n[which(native_list$taxa == species)]
    
    if(length(native_species) == 1 ){
      
      if(native_species == ""){
        
        df[i,]<- data.frame(species = species,sampstat = sampstat, native_countries = 'No native country for this species',n_accesions = NA, n_ecor = NA, n_ecor_in_native = NA)

      }else{
        
        ecor_native_countries<- ecorregiones_por_pais$ECO_NAME[which(ecorregiones_por_pais$iso3 %in% native_species)]
        ecor_native_countries<- unique(ecor_native_countries[!is.na(ecor_native_countries)])
        
        
        #ecor_native_countries<- fnl$ECO_NAME[which(fnl$GADM_GID_0 %in% native_species )]; ecor_native_countries<-unique(ecor_native_countries[!is.na(ecor_native_countries)])
        
        ecor_in_native_countries <- conteos$ECO_NAME[which(conteos$name_sc == species & conteos$sampstat_cat == sampstat & conteos$GADM_GID_0 %in% native_species)]
        ecor_in_native_countries<-unique(ecor_in_native_countries[!is.na(ecor_in_native_countries)])
        
        n_accessions <- sum(conteos$n[which(conteos$name_sc == species & conteos$sampstat_cat == sampstat & conteos$GADM_GID_0 %in% native_species)])
        
        df[i,]<- data.frame(species = species, sampstat = sampstat, native_countries = paste0(native_species, collapse = ", "),n_accesions = n_accessions , n_ecor = length(ecor_native_countries), n_ecor_in_native =  length(ecor_in_native_countries))
        
      }
      
      
    }else{
      
      ecor_native_countries<- ecorregiones_por_pais$ECO_NAME[which(ecorregiones_por_pais$iso3 %in% native_species)]
      ecor_native_countries<- unique(ecor_native_countries[!is.na(ecor_native_countries)])
      
      
      #ecor_native_countries<- fnl$ECO_NAME[which(fnl$GADM_GID_0 %in% native_species )]; ecor_native_countries<-unique(ecor_native_countries[!is.na(ecor_native_countries)])
      
      ecor_in_native_countries <- conteos$ECO_NAME[which(conteos$name_sc == species & conteos$sampstat_cat == sampstat & conteos$GADM_GID_0 %in% native_species)]
      ecor_in_native_countries<-unique(ecor_in_native_countries[!is.na(ecor_in_native_countries)])
      
      n_accessions <- sum(conteos$n[which(conteos$name_sc == species & conteos$sampstat_cat == sampstat & conteos$GADM_GID_0 %in% native_species)])
      
      df[i,]<- data.frame(species = species, sampstat = sampstat, native_countries = paste0(native_species, collapse = ", "),n_accesions = n_accessions , n_ecor = length(ecor_native_countries), n_ecor_in_native =  length(ecor_in_native_countries))
      

    }
    
    
  }
  
}


df$p_ecor_in_native<-ifelse(!is.na(df$n_ecor), df$n_ecor_in_native/df$n_ecor, NA)


write.csv(df, paste0(outdir, '/',collection_name,'/Eco/native_countries_cov_by_species.csv'), row.names = F )


########################################################


species_list<-unique(fnl$name_sc)
sampstat_list<-unique(fnl$sampstat_cat)
country_list <-unique(ecorregiones_por_pais$iso3); #country = country_list[1]



if(length(which(grepl('X|x', species_list))) != 0){
  
  species_list<- species_list[-which(grepl('X|x', species_list))]
  
}else{
  
  species_list<-species_list
  
}


df<-data.frame(country  = NA, is_native = NA,n_accessions = NA ,n_ecor_country = NA, species = NA,sampstat = NA, n_ecor_sp_in_country = NA )
i = 0

for(country in country_list){

  for(species in species_list){

    
    for(sampstat in sampstat_list){
      

      i= i+1
      
      native_species <- native_list$ISO3_n[which(native_list$taxa == species)]
      is_native <- ifelse(country %in% native_species, TRUE, FALSE)
      
      
      ecor_country<- ecorregiones_por_pais$ECO_NAME[which(ecorregiones_por_pais$iso3 == country)]
      n_ecor_country<- length(unique(ecor_country[!is.na(ecor_country)]))
      
      
      ecor_sp_in_country <- conteos$ECO_NAME[which(conteos$name_sc == species & conteos$sampstat_cat == sampstat & conteos$GADM_GID_0 == country)]
      n_ecor_sp_in_country<-length(unique(ecor_sp_in_country[!is.na(ecor_sp_in_country)]))
      
      
      if(n_ecor_sp_in_country != 0){
        
        n_accessions <- sum(conteos$n[which(conteos$name_sc == species & conteos$sampstat_cat == sampstat & conteos$GADM_GID_0 == country)])
        
        df[i,]<- data.frame(country  = country, is_native = is_native, n_accessions = n_accessions, n_ecor_country = n_ecor_country, species = species,sampstat = sampstat, n_ecor_sp_in_country = n_ecor_sp_in_country )
        
        
      }else{
        
        
        df[i,]<- data.frame(country  = country, is_native = is_native, n_accessions = 0, n_ecor_country = n_ecor_country, species = species,sampstat = sampstat, n_ecor_sp_in_country = 0 )
        
        
      }
   
    }
    
  }
}




df$p_ecor_sp_in_country <- df$n_ecor_sp_in_country/df$n_ecor_country


df<- df[-which(df$n_ecor_sp_in_country == 0),]

write.csv(df, paste0(outdir,'/',collection_name, '/Eco/countries_ecor_coverage_(per sp).csv'), row.names = F )


##################################################################################33

df_c<-df %>% group_by(country, sampstat) %>% summarise(p_ecor_sp_in_country = mean(p_ecor_sp_in_country[which(p_ecor_sp_in_country != 0)]))

write.csv(df_c, paste0(outdir,'/',collection_name, '/Eco/countries_ecor_coverage_summary.csv'), row.names = F )











