# An ex-situ conservation assessment diversity index (ECADI) for evaluating the coverage and diversity in genebank collections
_______

## Techincal details report:
A document with all technical documentation is avaialable here for reviewing: 
(https://github.com/alliance-datascience/genebanks_diversity_index/blob/dev/README_FILES/Index_report_extensive.pdf)
_______

## Graphical summary:

**Figure 1**. A graphical summary of the workflow is available here: 
(https://github.com/alliance-datascience/genebanks_diversity_index/blob/dev/README_FILES/DIAGRAM_METRICS.drawio.png)
_______
## Outputs examples: 

![Figure 2.](https://github.com/alliance-datascience/genebanks_diversity_index/blob/dev/README_FILES/beans_ECADI_PLOT.png)
**Figure 2**. Radar chart with the main metrics calculated for the four main components of The ex-situ conservation assessment diversity index (ECADI) (i) Collection composition and taxonomy (Composition), (ii) Documentation completeness (Documentation), (iii) Ecogeographical representativeness (Ecogeography), and (iv) Genetics diversity and genetic usability information availability (Usability). The blue color represents the values obtained for Crop Wild Relatives and red represents the values obtained for landraces for Beans CIAT collection
_______
![Figure 2.](https://github.com/alliance-datascience/genebanks_diversity_index/blob/dev/README_FILES/beans_ECADI_WILD.png)
**Figure 3**. NIPALS graphical representation of the four main indicator to evaluate crop wild relatives taxa in CIAT bean collection. Each color represents a world region according to World bank criteria.
_______
# Main components of the ECADI idea: 

This repository consists of a set of functions to calculate a series of indicators of coverage and diversity of a germplasm collection based on:

- (i) Collection composition and taxonomy
- (ii) Ecogeography
- (iii) Documentation completeness
- (iv) Genetics diversity and genetic usability information availability

_______
_______
_______
**Dependencies:** `R (>= 4.4.1)`

**Imports:** `base, utils, methods, stats, ggplot2, DescTools, parallel, stringr, WorldFlora, plsdepot, ggrepel, ggradar, countrycode, sf, raster, factoextra, terra, data.table, glue, fossil, dplyr, vegan, abdiv, DescTools`

_______
_______
_______
# Description:

This code performs the following steps:
- Preprocessing (0_preprocessing.R)
- Collection composition and taxonomy (1_Taxonomic_index_Country.R)
- (ii) Ecogeography (2_Ecogeographic_index.R)
- (iii) Documentation completeness (3_completeness_index.R)
- (iv) Genetics diversity and genetic usability information availability (4_Genetic_sequenced_coverage.R)
- (iv[a]) Genetics diversity (4_1_Genetics_distance.R)
- (v) Ex-situ conservation assessment diversity index (5_ECADI.R)
_______

# Recommendations:
> [!IMPORTANT]
> Please follow the next recommendations and instructions: 
- (0_preprocessing.R ) and Ecogeography (2_Ecogeographic_index.R) are the slowiest steps in this pipeline. Be patient
- Please provide the germplasm collection data in MCDP format. see (https://www.genesys-pgr.org/documentation/basics#mcpd)
- The Geographical quality score, and taxonomic geographical score can be obtained from https://github.com/alliance-datascience/genebank-general/tree/dev
- Only use accessions with high and medium geographical quality score to run the Ecogeography (2_Ecogeographic_index.R) step
- If genetic data is not available please use the step Genetics diversity and genetic usability information availability (4_Genetic_sequenced_coverage.R)

> [!CAUTION]
> The ECADI must be interpreted using the four components together!
- This code was designed to calculate Hill-Simpson (Effective Number of Species), and Margalef.
- The main output for the ECADI was calculated using Hill-Simpson!
  
_______
# How to run this code? This must be the order of running:
 ## 1) MANDATORY PREPROCESSING
 * [ ] (i) 0_preprocessing.R (MANDATORY) 
 ## 2) The following codes can be run regardless the order.
 * [ ] (ii) 1_Taxonomic_index_Country.R
 * [ ] (iii) 2_Ecogeographic_index.R
 * [ ] (iv) 3_completeness_index.R
 * [ ] (v) 4_Genetic_sequenced_coverage.R
  
> [!CAUTION]
> If you have a distance matrix, please use 4_1_Genetics_distance.R. Nevertheless, this pipeline was designed to use 4_Genetic_sequenced_coverage.R!

## 3) MANDATORY SUMMARY STEP
 * [ ] (vi) 5_ECADI.R (MANDATORY)


_______
# Description of each step: 
A detailed explanation of each code is provided as follows:


> [!IMPORTANT]
> Share inputs for codes:
- A path to save outputs: outdir
- Collection name: The collection name used (e.g. "beans")

## 0. Preprocessing (0_preprocessing.R )

### Inputs


> [!IMPORTANT]
> Download GRIN information from https://npgsweb.ars-grin.gov/gringlobal/uploads/documents/taxonomy_data.cab

> GRIN Inputs:
 - taxonomy_species.txt (taxonomy)
 - geography.txt (geographical information for taxa)
 - taxonomy_geography_map.txt (geographical information for taxa and connecting table)

> CIAT germplasm data in MCDP format:
- genesys-accessions-COL003.csv

> [!IMPORTANT]
> WorldFlora Online input:
> Download this file: https://files.worldfloraonline.org/files/WFO_Backbone/_WFOCompleteBackbone/WFO_Backbone.zip

 >[!IMPORTANT]
> Add a key for IUCN API: Please obtain your API signing up here:https://api.iucnredlist.org/ (The IUCN query step is quite slow!)

### Steps done:
- Filter by Crop wild relatives (CWR), landraces and hybrids. This step discard breeding materials!
- Creating count matrices per CWR and landraces: rows are taxa and columns are countries (output 1)
- Obtain taxonomic matching using GRIN taxonomy
- Obtain taxonomic matching using WorldFlora Online taxonomy
- Obtain native areas using GRIN taxonomy information (output 2)
- Obtain IUCN conservation status
- Obtain biodiversity general metrics (Atkinson, Simpson, Margalef)  (output 3) 
- Create a summary table with taxonomy status in GRIN and WorldFlora, and IUCN status (Main output)

### Outputs:
 -  Output 1: ["/", collection_name, "/", collection_name,"_count_matrix_countries_1.csv"]
 -  Output 2: ["/", collection_name, "/", collection_name, "_native_iso3_new_1.csv"]
 -  Output 3: ["/", collection_name, "/", collection_name, "_IT_table_1.csv"]
 -  Main output:  ["/", collection_name, "/", collection_name, "_summary_table.csv"]


```r
passport_data_orig <- passport_data_original[which(passport_data_original$CROPCODE =="beans"), ]
collection_name = "beans"
x2 <- index_function(
  passport_data_orig = passport_data_orig,
  tax = tax,
  geo = geo,
  tax_geo = tax_geo,
  outdir = outdir,
  API = API,
  collection_name = collection_name)
```

## 1. Collection composition and taxonomy (1_Taxonomic_index_Country.R)

### Inputs
> [!IMPORTANT]
> This code uses the RDS calculated in the preprocessing step 
>  ["/", collection_name, "/", collection_name,"_subsets_new_1.RDS"]
> Preprocessing Inputs:
 - Summary table with taxonomy status in GRIN and WorldFlora, and IUCN status - geography.txt (geographical information for taxa)
 - taxonomy_geography_map.txt (geographical information for taxa and connecting table)

### Steps done:
- Call previous results using the RDS file from 0_preprocessing.R
- Call summary table (0_preprocessing.R)
- Obtaining taxa proportion in GRIN and WorlFlora per crop wild relatives (CWR), landraces, and hybrids subsets
- Obtaining Coverage metrics (Chao 2012). (Alternative)
- Obtaining Effective number of species (Hill-Simpson, as well known as ENS), Margalef, 1- Atkison (1-A) per collection subset and countries
- Obtaining Abundance Coverage estimator per collection subset and countries (alternative)
- Normalize ENS and Margalef using Sigmoid scaling *D is either ENS or Margalef
- Calculating composition index as D*1-A in a summary file per country and collection subset (Main output)
- Obtaining regions using ISO3 information from count matrices
- Plotting using NIPALS for landraces and CWR (Output 2A or Output 2B)
          
### Outputs:
 - Main output: ["/",   collection_name,"/",collection_name,"_composition_countries.csv"]
-  (Output 2A or Output 2B): [ "/",  collection_name,"/", collection_name, "_",plot_type,"_composition_LANDRACE.png"] or  ["/",collection_name,"/",collection_name, "_",plot_type,"_composition_WILD.png"]
>plot_type refers to what kind of biodiversity index (“simpson”,”margalef”) for ENS or Margalef.

```r
  collection_name <- "beans"
#Plotting NIPALS using Margalef 
  x2 <- composition_country(outdir = outdir,
                            collection_name = collection_name,
                            plot_type="Margalef")
#Plotting NIPALS using ENS 
  x2 <- composition_country(outdir = outdir,
                            collection_name = collection_name,
                            plot_type="Simpson")
```

## 2. Ecogeography (2_Ecogeographic_index.R )
> [!IMPORTANT]
> Download Terrestrial Ecoregions of the World https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip

## 3. Documentation completeness (3_completeness_index.R)
>[!NOTE]
>A detailed explanation of the fields used in Taxonomic Quality Score Index (TQS) are in the table 1, below.
### Inputs
> [!IMPORTANT]
> This code uses the RDS calculated in the preprocessing step 
>  ["/", collection_name, "/", collection_name,"_subsets_new_1.RDS"]
> Preprocessing Inputs:
 - Passport data for CWR, landrace, and hybrids respectively: (passport_data_w, passport_data_l, passport_data_h)
- Collection data with Passport Data Completeness Index (PDCI) [PCDI_df]
- Collection data with Geographical Quality Score Index (GQS) [QGI_df]
- Collection data with Taxonomic Quality Score Index (TQS) [TQI_df] 

### Steps done:
- Call previous results using the RDS file from 0_preprocessing.R
-   Joining PDCI, GQS, and TQS to passport data in collection subsets
- Obtaining summarizing metrics (mean, sd) for PDCI, GQS, and TQS and standardize to obtain values from 0 to 1 for the three quality data metrics
- Obtaining indicator per country as the geometric mean of species per country and collection subset (Main output) 
- Obtaining regions using ISO3 information from count matrices
- Plotting using NIPALS for landraces and CWR (Output 2A or Output 2B) 
        
### Outputs:
 - Main output: ["/",   collection_name,"/",collection_name,"_completeness_countries.csv"]
 -  (Output 2A or Output 2B) [ "/",  collection_name,"/", collection_name",_"completeness_LANDRACE.png"] or  [ "/",  collection_name,"/", collection_name",_"completeness_WILD.png"] 

```r
collection_name <- "beans"
PCDI_df_beans <- PCDI_df[which(PCDI_df$CROPNAME==collection_name),]
QGI_df_beans <- QGI_df[which(QGI_df$CROPNAME==collection_name),]
TQI_df_beans <- TQI_df[which(TQI_df$CROPNAME==collection_name),]
x1 <-completeness_func(outdir = outdir,
                       collection_name = collection_name,
                       PCDI_df = PCDI_df_beans,
                       QGI_df = QGI_df_beans,
                       TQI_df =TQI_df_beans)
```

## 4. Genetics diversity and genetic usability information availability (4_Genetic_sequenced_coverage)

### Inputs
> [!IMPORTANT]
> This code uses the RDS calculated in the preprocessing step 
>  ["/", collection_name, "/", collection_name,"_subsets_new_1.RDS"]
- numCores = Cores number used to match information
- An excel file with the accessions sequenced. The information reported is the accession number (accessions_df)
> Preprocessing Inputs:
 - Summary table with taxonomy status in GRIN and WorldFlora, and IUCN status ["/", collection_name,"/", collection_name,  "_summary_table.csv"]
- Passport data for CWR, landrace, and hybrids respectively: (passport_data_w, passport_data_l, passport_data_h)

### Steps done:
- Call previous results using the RDS file from 0_preprocessing.R
-   Checking if CWR, landraces, and hybrid datasets are available
- Adding a flag in passport data to obtain what accessions were sequenced in collection subsets
- Using information from summary table to obtain IUCN threatened species in collection subsets and country using the following categories:  "Vulnerable”, “Endangered”, “Critically Endangered”, “Extinct in the Wild". This metric is IUCN threathed/total taxa in the subset.
- Obtaining summarizing metrics per country for proportion of records sequenced (pRseq) and proportion of taxa sequenced (pTseq) per each collection subset
- Obtaining usability index as the average of pRseq, pTseq, and IUCN.
> [!Note]
> The metrics for the collection subsets is a weighted mean considering the number of records in a total collection (CWR, landraces, and hybrids).
- Save a summary table of the subset collection including taxa proportions: (Output 1) 
- A general summary table for collection subsets (Output 2)
- Obtaining regions using ISO3 information from count matrices
- Obtaining indicator per country and collection per country and collection subset (Main output) 
- Plotting using NIPALS for landraces and CWR (Output 2A or Output 2B) 
       
### Outputs:
- Output 1: ["/", collection_name,"/", collection_name, "_4_genetics_summary_table.csv]
- Output 2: ["/", collection_name,"/", collection_name, "_4_GI_table.csv"]
- Main output: ["/", collection_name,"/", collection_name, "_4_GI_country_table.csv]
-  (Output 2A or Output 2B) [ "/",  collection_name,"/", collection_name",_"usability_land.png"] or  [ "/",  collection_name,"/", collection_name",_"usability_WILD.png"]

```r
collection_name = "beans"
x1 <- genetics_ind_function(outdir=outdir, 
                            collection_name=collection_name, 
                            numCores=4,
                            accessions_df=accessions_df)
```

## 4.1. Genetics diversity (When it is available!) (4_1_Genetics_distance)

### Inputs
> [!IMPORTANT]
> This code uses the RDS calculated in the preprocessing step 
>  ["/", collection_name, "/", collection_name,"_subsets_new_1.RDS"]
- numCores = Cores number used to match information
- An excel file with the genetic distance for the overall collection!  This matrix must have the same number of columns and rows! (gen_data)
> Preprocessing Inputs:
- Passport data for CWR, landrace, and hybrids respectively: (passport_data_w, passport_data_l, passport_data_h)

### Steps done:
- Call previous results using the RDS file from 0_preprocessing.R
-   Checking if CWR, landraces, and hybrid datasets are available
- Adding a flag in passport data to obtain what accessions were used for the genetic distance calculation in collection subsets
- Sub setting genetic distances per collection subsets and calculating the genetic distance median per accession (Output 1)
- Calculating the median of genetic distances per taxon and country in each collection subset (Main Output)
       
### Outputs:
 - Output 1: ["/", collection_name,"/", collection_name, "_4_1_Genetic_distance_summary.csv"]
- Main output: ["/", collection_name,"/", collection_name, "_4_1_Genetic_distance_country.csv"]

```r
x <- genetics_dist_function(outdir=outdir,
                       collection_name=collection_name,
                       numCores=numCores,
                       gen_data=gen_data)
```
## 5. ECADI index calculation  (5_ECADI.R)

### Inputs
> [!IMPORTANT]
> This code uses the results per country calculated in steps 1 to 4.
>  eco_file: the file with the ecogeographical index calculated in the 2_Ecogeographic_index.R file

```r
collection_name <- "beans"
x1 <- ECADI_function(outdir = outdir,
                           eco_file = eco_file,
                           collection_name = collection_name
)
```



### Steps done:
- Load the composition index: ["/", collection_name,"/", collection_name, "_composition_countries.csv"]
- Load the ecogeographical representativeness index from the eco_file
- Load the documentation completeness index: ["/", collection_name,"/", collection_name, "_completeness_countries.csv"]
- Load the usability index: ["/", collection_name,"/", collection_name, "_4_GI_country_table.csv"]
- Load the genetic distance index file: ["/", collection_name,"/", collection_name, "_4_1_Genetic_distance_country.csv"]
- Create a summary file using the composition index file as template 
- Match the four indexes and countries for CWR and landraces
- Calculate the ECADI as 
```
    sum(summary_file$C1_COMP_INDEX_[]_ENS[[i]] +
          summary_file$C2_ECO_INDEX_[][[i]] +
          summary_file$C3_DC_INDEX_[][[i]] +
          summary_file$C4_USAB_INDEX_[][[i]],na.rm = T)/4
```
>[] represents either CWR or landraces

-	Calculate ECADI for a collection as the average of ECADI (landraces) and ECADI (CWR)
-	Save a CSV file (Main output) 
-	-Obtaining regions using ISO3 information from count matrices
-	-Plotting using NIPALS for landraces and CWR (Output 2A or Output 2B) 
-	Plot indexes as Spider plot (Output 3)

### Outputs:
- Main output: ["/", collection_name,"/", collection_name, "_ECADI.CSV]
- (Output 2A or Output 2B): [ "/",  collection_name, "/", collection_name,      "_",  "ECADI_WILD.png"]
- Output3: ["/",  collection_name,  "/",   collection_name,  "_",  "ECADI_PLOT.png]



_______
_______
_______
# Authors

Main: Chrystian C. Sosa, and Maria Victoria Diaz 
Other contributors: Norma Giraldo
_______
_______
_______

## Additional information
Table 1. The Multi-crop Passport Descriptors (MCPD) V.2.1 (Alercia et al., 2015) used for the TCS. Extracted from https://www.genesys-pgr.org/documentation/basics#accedoc-tax. 

field name | Description |
------------ | ------------- |
GENUS | The genus name for taxon. An initial uppercase letter is required
SPECIES | Specific epithet portion of the scientific name in lowercase letters. The abbreviation sp. or spp. is allowed when the exact species name is unknown.
SPAUTHOR |The authority for the species name.
SUBTAXA | A subtaxon can be used to store any additional taxonomic identifier. The following abbreviations are allowed: subsp. (for subspecies); convar. (for convariety); var. (for variety); f. (for form); Group (for cultivar group).
SUBTAUTHOR | The subtaxon authority at the most detailed taxonomic level.
_______
