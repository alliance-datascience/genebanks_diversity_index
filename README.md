# An ex-situ conservation assessment diversity index (ECADI) for evaluating the coverage and diversity in genebank collections

This repository consists of a set of functions to calculate a series of indicators of coverage and diversity of a germplasm collection based on:

- (i) Collection composition and taxonomy
- (ii) Ecogeography
- (iii) Documentation completeness
- (iv) Genetics diversity and genetic usability information availability
 





# Description:

This code performs the following steps:
- Preprocessing (0_preprocessing.R )
- Collection composition and taxonomy (1_Taxonomic_index_Country.R)
- (ii) Ecogeography (2_Ecogeographic_index.R)
- (iii) Documentation completeness (3_completeness_index)
- (iv) Genetics diversity and genetic usability information availability (4_Genetic_sequenced_coverage.R)
- (iv) Genetics diversity (4_1_Genetics_distance.R)

A detailed explanation of each code is provided as follows:


### Preprocessing (0_preprocessing.R )
- Inputs:




```r
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
```




# Additional information
Table 1. The Multi-crop Passport Descriptors (MCPD) V.2.1 (Alercia et al., 2015) used for the TCS. Extracted from https://www.genesys-pgr.org/documentation/basics#accedoc-tax. 

field name | Description |
------------ | ------------- |
GENUS | The genus name for taxon. An initial uppercase letter is required
SPECIES | Specific epithet portion of the scientific name in lowercase letters. The abbreviation sp. or spp. is allowed when the exact species name is unknown.
SPAUTHOR |The authority for the species name.
SUBTAXA | A subtaxon can be used to store any additional taxonomic identifier. The following abbreviations are allowed: subsp. (for subspecies); convar. (for convariety); var. (for variety); f. (for form); Group (for cultivar group).
SUBTAUTHOR | The subtaxon authority at the most detailed taxonomic level.
