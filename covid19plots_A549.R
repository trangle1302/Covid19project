
library(calibrate)
# Volcano plot 
res = read.csv('./Foldchange_meanintensity_CZ8780.csv', header=TRUE)
res$log2FoldChange = log2(res$fc)
res$padj = p.adjust(res$pval, method = 'BH', n = length(res$pval))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pval), labs=Gene, cex=.8))


plot(res$log2FoldChange, -log10(res$pval), type = "n")  # Create an empty plot
points(res$log2FoldChange, -log10(res$pval), pch = 19, col = "black")

subset_data <- subset(res, padj < 0.05 & abs(log2FoldChange) < 1)
points(subset_data$log2FoldChange, -log10(subset_data$pval), pch = 19, col = "red")
# Subset the data
subset_data <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
# Add points to the plot
points(subset_data$log2FoldChange, -log10(subset_data$pval), pch = 19, col = "green")
# Add text labels using textxy()
#textxy(subset_data$log2FoldChange, -log10(subset_data$pval), labs = subset_data$Gene, cex = 0.8)

######################################################
library(ggplot2)
library(tidyr)
library(circlize)

#df = read.csv('C:/Users/trang.le/Downloads/Review_Analysis_Spatial_var_filtered_for_SCV.csv', header=TRUE)
df = read.csv('./Review_spatial_var.csv', header=TRUE)

df[,"Infected"] <- sapply(df[,"Infected"], as.character)
df[,"Non.infected"] <- sapply(df[,"Non.infected"], as.character)

loc <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34)
names(loc) <- c("Nucleoplasm","Nuclear membrane","Nucleoli","Nucleoli fibrillar center","Nuclear speckles",
                "Nuclear bodies","Kinetochore","Mitotic chromosome","Endoplasmic reticulum","Golgi apparatus",
                "Vesicles","Peroxisomes","Endosomes","Lysosomes","Intermediate filaments",
                "Actin filaments","Focal adhesion sites","Microtubules","Microtubule ends","Cytokinetic bridge",
                "Midbody","Midbody ring","Cleavage furrow","Mitotic spindle","Primary cilia",
                "Centriolar satellite","Centrosome","Lipid droplets","Plasma membrane","Cell Junctions",
                "Mitochondria","Aggresome","Cytosol","Cytoplasmic bodies","Rods & Rings")

names(loc) <- c("Nucleoplasm","Nuclear Membrane","Nucleoli","Nucleoli","Nuclear speckles",
                "Nuclear bodies","Nucleoplasm","Nucleoplasm","Endoplasmic Reticulum","Golgi Apparatus",
                "Vesicles","Vesicles","Vesicles","Vesicles","Intermediate Filaments",
                "Actin Filaments","Actin Filaments","Microtubules","Microtubules","Microtubules",
                "Microtubules","Microtubules","Microtubules","Microtubules","Centrosome",
                "Centrosome","Centrosome","Vesicles","Plasma Membrane","Cell Junctions",
                "Mitochondria","Aggresome","Cytosol","Cytosol","Cytosol")

loc_order = c("Nuclear Membrane","Nucleoli", "Nucleoplasm", "Nuclear speckles","Nuclear bodies",
              "Actin Filaments", "Centrosome", "Cytosol",
              "Intermediate Filaments", "Microtubules","Mitochondria",
              "Endoplasmic Reticulum", "Golgi Apparatus", "Plasma Membrane", "Vesicles", 
              "Cell Junctions", "Aggresome")

loc_color = c( "darkorange","firebrick1", "firebrick",  "darkorange4","darksalmon",
              "deepskyblue", "deepskyblue3", "dodgerblue3", 
              "chartreuse", "chartreuse4", "darkolivegreen1",
              "goldenrod1", "gold","gold3", "khaki1", 
              "aquamarine","coral")


location_df = data.frame("Location" = loc,
                         "NonInfected"= names(loc),
                         "Infected"= names(loc),
                         "flow" = 0)
location_df = crossing(NonInfected = names(loc), Infected = names(loc), flow=0)
#rownames(location_df) = loc

location_df_nonif = location_df
location_df_if = location_df


for (i in 1:nrow(df)) {
  locations_infected = na.omit(as.numeric(strsplit(unlist(df[i,'Infected'])[[1]], "[^0-9]+")[[1]]))
  locations_noninfected = na.omit(as.numeric(strsplit(unlist(df[i,'Non.infected']), "[^0-9]+")[[1]]))
  
  for (l0 in locations_noninfected) {
    for (l1 in locations_infected) {
      current_val = location_df[(location_df$NonInfected == names(loc[loc==l0]))&
                                  (location_df$Infected == names(loc[loc==l1])), "flow"]
      location_df[(location_df$NonInfected == names(loc[loc==l0]))&
                    (location_df$Infected == names(loc[loc==l1])), "flow"] = current_val+1
    }
  }  
  #if length(locations_noninfected)>1 {}
}

### Circular plot for spatial hits
circos.clear()
circos.par(start.degree=90)#,gap.degree=c(rep(2,nrow(mat)-1),20,rep(2,ncol(mat)-1),20))
library("reshape2")
mat = dcast(location_df, NonInfected ~ Infected)
rownames(mat) = mat$NonInfected
colnames(mat) = paste0(colnames(mat),' ')
mat <- as.matrix(mat[,-1])
diag(mat) = 0

order=c(loc_order, rev(paste0(loc_order,' '))) #
#order = c(rownames(mat), rev(colnames(mat)))

col_inf = loc_color#col
names(col_inf) = paste0(loc_order,' ')
col_noninf = loc_color#col
names(col_noninf) = loc_order #rownames(mat)

#chordDiagram(mat, order=order, grid.col=c(col_green,col_red))
chordDiagram(mat, order=rev(order), grid.col=rev(c(col_noninf,rev(col_inf))), 
             transparency = 0.15,
             directional = 1,reduce = 0,
             direction.type = c("arrows", "diffHeight"), 
             annotationTrack = "grid", annotationTrackHeight = c(0.05, 0.1),
             link.arr.type = "big.arrow",preAllocateTracks=list(track.height=0.3),
             link.sort = TRUE, link.largest.ontop = TRUE)

circos.trackPlotRegion(track.index = 1,
                       panel.fun = function(x, y){
                         xlim =get.cell.meta.data("xlim")
                         ylim =get.cell.meta.data("ylim")
                         sector.name =get.cell.meta.data("sector.index")
                         circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj =c(0, 0.5))
                       }, bg.border = NA)


