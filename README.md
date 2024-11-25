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

## Preprocessing (0_preprocessing.R )

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

> Add a key for IUCN API: Please obtain your API signing up here:https://api.iucnredlist.org/

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



## Ecogeography (2_Ecogeographic_index.R )
> [!IMPORTANT]
> Download OLSON (2001)

> 



# Additional information
Table 1. The Multi-crop Passport Descriptors (MCPD) V.2.1 (Alercia et al., 2015) used for the TCS. Extracted from https://www.genesys-pgr.org/documentation/basics#accedoc-tax. 

field name | Description |
------------ | ------------- |
GENUS | The genus name for taxon. An initial uppercase letter is required
SPECIES | Specific epithet portion of the scientific name in lowercase letters. The abbreviation sp. or spp. is allowed when the exact species name is unknown.
SPAUTHOR |The authority for the species name.
SUBTAXA | A subtaxon can be used to store any additional taxonomic identifier. The following abbreviations are allowed: subsp. (for subspecies); convar. (for convariety); var. (for variety); f. (for form); Group (for cultivar group).
SUBTAUTHOR | The subtaxon authority at the most detailed taxonomic level.
