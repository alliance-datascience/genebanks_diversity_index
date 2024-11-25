require(data.table)
require(readxl)
library(ggplot2)
library(factoextra)
library(glue)
library(fossil)
library(plsdepot)
library(ggrepel)
require(parallel)
library(plsdepot)
library(countrycode)
library(ggradar)
library(ggvanced)
#devtools::install_github("Ringomed/ggvanced")
#devtools::install_github("ricardo-bion/ggradar")
#https://www.datanovia.com/en/blog/beautiful-radar-chart-in-r-using-fmsb-and-ggplot-packages/
ECADI_function <- function(outdir,
                           eco_file,
                           collection_name
                           ){
################################################################################  
################################################################################  
message("  ")
message(paste("RUNNING: ",collection_name))
message("  ")
##
message("Loading component results")
message("  ")
################################################################################  
################################################################################  
message("Loading results from 1_Taxonomic_index_Country")
COMP_1 <-   paste0(outdir, "/", collection_name,"/", collection_name, "_composition_countries.csv")
if(file.exists(COMP_1)){
  COMP_1 <- read.csv(
    COMP_1,
    row.names = 1)  
} else {
  rm(COMP_1)
}
################################################################################
message("Loading results from 2_ECOGEOGRAPHY")
ECO_2 <-   as.data.frame(eco_file[which(eco_file$Collection==collection_name),])
row.names(ECO_2) <- ECO_2$Country
ECO_2$Country <- NULL
################################################################################
message("Loading results from 3_completeness_index")
DCI_3 <-   paste0(outdir, "/", collection_name,"/", collection_name, "_completeness_countries.csv")
if(file.exists(DCI_3)){
  DCI_3 <- read.csv(
    DCI_3,
    row.names = 1)  
} else {
  rm(DCI_3)
}
################################################################################
message("Loading results from 4_Genetic_sequenced_coverage")
USA_COV_4 <-   paste0(outdir, "/", collection_name,"/", collection_name, "_4_GI_country_table.csv")
if(file.exists(USA_COV_4)){
  USA_COV_4 <- read.csv(
    USA_COV_4,
    row.names = 1)  
} else {
  rm(USA_COV_4)
}
################################################################################
message("Loading results from 4_1_Genetics_distance")
GEN_DIST_5 <-  paste0(outdir, "/", collection_name,"/", collection_name, "_4_1_Genetic_distance_country.csv")
if(file.exists(GEN_DIST_5)){
  GEN_DIST_5 <- read.csv(
    GEN_DIST_5,
    row.names = 1)  
} else {
  message("NOT AVAILABLE GENETIC DATA!")
  rm(GEN_DIST_5)
}
################################################################################
################################################################################
#creating a summary file
message(" ")
message("creating a summary file")
message(" ")
################################################################################
################################################################################
summary_file <- data.frame(country=row.names(COMP_1),
                           REGION =COMP_1$REGION,
                           C1_COMP_INDEX_W_ENS = COMP_1$COMP_INDEX_W_ENS,
                           C1_COMP_INDEX_L_ENS = COMP_1$COMP_INDEX_L_ENS,
                           C2_ECO_INDEX_W = NA,
                           C2_ECO_INDEX_L = NA,
                           C3_DC_INDEX_W = NA,
                           C3_DC_INDEX_L = NA,
                           C4_USAB_INDEX_W = NA,
                           C4_USAB_INDEX_L = NA,
                           C4_1_GEN_DIST_W  = NA,
                           C4_1GEN_DIST_L  = NA,
                           ECADI_W_US = NA,
                           ECADI_L_US = NA,
                           ECADI_COLL_US = NA,
                           ECADI_W_G = NA,
                           ECADI_L_G = NA,
                           ECADI_COLL_G = NA
                           )

for(i in 1:nrow(summary_file)){
  #Component 2 (ECO)
  #i <- 1
  C2_i <- ECO_2[which(row.names(ECO_2) %in% summary_file$country[[i]]),]
  if(nrow(C2_i)==0){
    C2_i[1,] <- NA
  }
  summary_file$C2_ECO_INDEX_W[[i]] <- C2_i$Wild
  summary_file$C2_ECO_INDEX_L[[i]] <- C2_i$Landrace
  rm(C2_i)
  #Component 3 (DCI)
  C3_i <- DCI_3[which(row.names(DCI_3) %in% summary_file$country[[i]]),]
  if(nrow(C3_i)==0){
    C3_i[1,] <- NA
  }
  summary_file$C3_DC_INDEX_W[[i]] <- C3_i$SCORE_wild
  summary_file$C3_DC_INDEX_L[[i]] <- C3_i$SCORE_landrace
  rm(C3_i)
  #Component 4 (USABILITY)
  C4_i <- USA_COV_4[which(row.names(USA_COV_4) %in% summary_file$country[[i]]),]
  if(nrow(C4_i)==0){
    C4_i[1,] <- NA
  }
  summary_file$C4_USAB_INDEX_W[[i]] <- C4_i$usability_index_wild
  summary_file$C4_USAB_INDEX_L[[i]] <- C4_i$usability_index_landrace
  rm(C4_i)
  #Component 4.1 (GEN_DIST)
  #C4_1_i <- GEN_DIST_5[which(row.names(GEN_DIST_5) %in% summary_file$country[[i]]),]
  summary_file$C4_1_GEN_DIST_W[[i]] <- NA
  summary_file$C4_1_GEN_DIST_W[[i]] <- NA
  ##############################################################################
  #CALCULATING ECADI
  #CROP WILD RELATIVE
  summary_file$ECADI_W_US[[i]] <- 
    sum(summary_file$C1_COMP_INDEX_W_ENS[[i]] +
          summary_file$C2_ECO_INDEX_W[[i]] +
          summary_file$C3_DC_INDEX_W[[i]] +
          summary_file$C4_USAB_INDEX_W[[i]],na.rm = T)/4
  #LANDRACES
  summary_file$ECADI_L_US[[i]] <- 
    sum(summary_file$C1_COMP_INDEX_L_ENS[[i]] +
          summary_file$C2_ECO_INDEX_L[[i]] +
          summary_file$C3_DC_INDEX_L[[i]] +
          summary_file$C4_USAB_INDEX_L[[i]],na.rm = T)/4
  
}

summary_file$ECADI_COLL_US <- (summary_file$ECADI_W_US+summary_file$ECADI_L_US)/2

row.names(summary_file) <- summary_file$country
################################################################################
write.csv(
  summary_file,
  paste0(outdir, "/", collection_name,"/", collection_name, "_ECADI.CSV"),
  row.names = F,
  na = ""
)
############################################################################
message(paste0("Processing data for plots"))
############################################################################
x_res_W_ORI <- 
  summary_file[, c(
    "C1_COMP_INDEX_W_ENS",
    "C2_ECO_INDEX_W",
    "C3_DC_INDEX_W",
    "C4_USAB_INDEX_W",
    "ECADI_W_US",
    "REGION")]
colnames(x_res_W_ORI) <- c("Collection composition index",
                       "Ecogeographical representativeness index",
                       "Documentation completeness index",
                       "Usability index",
                       "ECADI",
                       "REGION"
)

x_res_W <- x_res_W_ORI[complete.cases(x_res_W_ORI),]
x_res_W$REGION <- factor(x_res_W$REGION)

###############################################
x_res_L_ORI <- 
  summary_file[, c(
    "C1_COMP_INDEX_L_ENS",
    "C2_ECO_INDEX_L",
    "C3_DC_INDEX_L",
    "C4_USAB_INDEX_L",
    "ECADI_L_US",
    "REGION")]
colnames(x_res_L_ORI) <- c("Collection composition index",
                       "Ecogeographical representativeness index",
                       "Documentation completeness index",
                       "Usability index",
                       "ECADI",
                       "REGION"
)

x_res_L <- x_res_L_ORI[complete.cases(x_res_L_ORI),]
x_res_L$REGION <- factor(x_res_L$REGION)

################################################################################
#seed for NIPALS
set.seed(1000)
################################################################################
message("plotting results for wild")

if(nrow(x_res_W)>=2){
  x_W_A = nipals(x_res_W[, -ncol(x_res_W)], comps = 3,scaled = T)
  x_res_W_NA <- x_res_W[!row.names(x_res_W) %in% c("No_country","Subset"),]
  # Create a scatter plot
  ind_coords <- data.frame(Dim.1=x_W_A$scores[,1],Dim.2=x_W_A$scores[,2],
                           REGION=x_res_W$REGION,countries=row.names(x_W_A$score))
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
    ggtitle(paste0("ex-situ conservation assessment diversity index for ","Crop Wild Relatives (", nrow(x_res_W_NA), " countries)"))+
    
    geom_segment(data = var_coords, aes(x = 0, y = 0, xend = t1*(r * 0.7), yend = t2*(r * 0.7)),colour = "steelblue",
                 arrow = arrow(length = unit(0.2, "cm")))+
    
    geom_text_repel(
      data          = var_coords,
      mapping       = aes(t1*(r * 0.7), t2*(r * 0.7), label = vars),
      size = 3,    # Repel away from the left edge, not from the right.
      color = "steelblue", 
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
  
  ggplot2::ggsave(
    paste0(
      outdir,
      "/",
      collection_name,
      "/",
      collection_name,
      "_",
      "ECADI_WILD.png"
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
  
  REGION <- x_res_L$REGION
  if(length(unique(x_res_L[,1]))>1){
    x_W_A = nipals(x_res_L[, -ncol(x_res_L)], comps = 3,scaled = T)
    
    x_res_L_NA <- x_res_L[!row.names(x_res_L) %in% c("No_country","Subset"),]
    # Create a scatter plot
    ind_coords <- data.frame(Dim.1=x_W_A$scores[,1],Dim.2=x_W_A$scores[,2],
                             REGION=x_res_L$REGION,countries=row.names(x_W_A$score))
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
      ggtitle(paste0("ex-situ conservation assessment diversity index  for ","landraces (", nrow(x_res_L_NA), " countries)"))+
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
        
        "ECADI_LAND.png"
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
message("plotting collection metrics!")

subset_to_radar_d <- data.frame(matrix(ncol = 6,nrow=2))
colnames(subset_to_radar_d) <- colnames(x_res_W_ORI)
subset_to_radar_d[1,] <- x_res_W_ORI[which(row.names(x_res_W_ORI)=="Subset"),]
subset_to_radar_d[2,] <- x_res_L_ORI[which(row.names(x_res_L_ORI)=="Subset"),]
#row.names(subset_to_radar_d) <- c("Crop wild relatives","Landraces")
subset_to_radar_d$REGION <- NULL

subset_to_radar_d <- round(subset_to_radar_d,3)
subset_to_radar_d$Subset <- c("Crop Wild Relatives","Landraces")
subset_to_radar_d <- subset_to_radar_d[,c(6,1:5)]
row.names(subset_to_radar_d) <- subset_to_radar_d$Subset
colnames(subset_to_radar_d) <-c("Collection subset",
                                "Composition",
                                "Ecogeography",
                                "Documentation completeness",
                                "Usability",
                                "ECADI")

# plt <- ggspider(subset_to_radar_d, 
#                 scaled = F, 
#                 axis_name_font_size=3,
#                 polygon = F)
plt <- ggradar(
  subset_to_radar_d,
  values.radar = c("0", "0.5", "1"),
  grid.min = 0, grid.mid = 0.5,grid.max = 1,
  # Polygons
  group.line.width = 1,
  group.point.size = 2,
  #legend.title = paste(collection_name,"collection"),
  group.colours = c("#00AFBB", "#FC4E07"),
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  grid.label.size = 3,  # Affects the grid annotations (0%, 50%, etc.)
  axis.label.size = 2, # Afftects the names of the variables
  legend.position = "bottom"
  )
ggsave(
  filename = 
    paste0(
      outdir,
      "/",
      collection_name,
      "/",
      collection_name,
      "_",
      
      "ECADI_PLOT.png"
    ),
  plot = plt,
  width = 8,
  height = 5,
  device = "png"
)###############################################################################
message(paste("FINISHED:",collection_name))
message("DONE!")
return(summary_file)
}
################################################################################
dir <- "D:/OneDrive - CGIAR/GERMPLASM_INDEX"
#outdir
outdir <- paste0(dir, "/BEANS/RESULTS")
################################################################################

################################################################################
eco_file <- "D:/OneDrive - CGIAR/GERMPLASM_INDEX/BEANS/RESULTS/eco_summary_all_coll 1.xlsx"
eco_file <- readxl::read_xlsx(eco_file,sheet = "country_summary")
################################################################################
################################################################################
################################################################################
collection_name <- "beans"
################################################################################
x1 <- ECADI_function(outdir,
                           eco_file,
                           collection_name
)
################################################################################
collection_name <- "cassava"
################################################################################
x2 <- ECADI_function(outdir,
                     eco_file,
                     collection_name
)
################################################################################
collection_name <- "forages"
################################################################################
x3 <- ECADI_function(outdir,
                     eco_file,
                     collection_name
)
