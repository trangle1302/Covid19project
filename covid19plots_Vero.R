library(ggplot2)
library(tidyr)
library(circlize)

df = read.csv('C:/Users/trang.le/Downloads/Review_Analysis_Spatial_var_filtered_for_SCV.csv', header=TRUE)
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



names(loc) <- c("Nucleoplasm","Nuclear Membrane","Nucleoli","Nucleoli","Nucleoplasm",
  "Nucleoplasm","Nucleoplasm","Nucleoplasm","Endoplasmic Reticulum","Golgi Apparatus",
  "Vesicles","Vesicles","Vesicles","Vesicles","Intermediate Filaments",
  "Actin Filaments","Actin Filaments","Microtubules","Microtubules","Microtubules",
  "Microtubules","Microtubules","Microtubules","Microtubules","Centrosome",
  "Centrosome","Centrosome","Vesicles","Plasma Membrane","Plasma Membrane",
  "Mitochondria","Cytosol","Cytosol","Cytosol","Cytosol")

loc_order = c("Nuclear Membrane","Nucleoli", "Nucleoplasm", 
          "Actin Filaments", "Centrosome", "Cytosol",
          "Intermediate Filaments", "Microtubules","Mitochondria",
          "Endoplasmic Reticulum", "Golgi Apparatus", "Plasma Membrane", "Vesicles")

loc_color = c("firebrick1", "firebrick", "darkorange",
              "deepskyblue", "deepskyblue3", "dodgerblue3", 
              "chartreuse", "chartreuse4", "darkolivegreen1",
              "goldenrod1", "gold3", "khaki1", "gold")


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

chordDiagram(location_df, annotationTrack = "grid", 
             grid.col=col,
             preAllocateTracks=list(track.height=0.3),
             link.sort = TRUE, link.largest.ontop = TRUE)
circos.trackPlotRegion(track.index = 1,
                       panel.fun = function(x, y){
                         xlim =get.cell.meta.data("xlim")
                         ylim =get.cell.meta.data("ylim")
                         sector.name =get.cell.meta.data("sector.index")
                         circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj =c(0, 0.5))
                       }, bg.border = NA)

### Circular plot for spatial hits
circos.clear()
circos.par(start.degree=90)#,gap.degree=c(rep(2,nrow(mat)-1),20,rep(2,ncol(mat)-1),20))

mat = dcast(location_df, NonInfected ~ Infected)
rownames(mat) = mat$NonInfected
colnames(mat) = paste0(colnames(mat),' ')
mat <- as.matrix(mat[,-1])

diag(mat) = 0

order=c(loc_order, rev(paste0(loc_order,' ')))#c(rownames(mat), rev(colnames(mat)))


colfunc <- colorRampPalette(c("gray","red","blue","green","yellow"))
col = colfunc(13)
col_inf = loc_color#col
names(col_inf) = paste0(loc_order,' ')
col_noninf = loc_color#col
names(col_noninf) = rownames(mat)

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



'''
# Circular plot
filld <- data.frame (start = c(seq(2,36),seq(38,38+34)), end = c(seq(2,36),seq(38,38+34))+1,
                     label = c(names(loc), names(loc)))

filld$p <- rowMeans(subset(filld, select = c(start, end)))
ggplot(filld, aes(xmin = start, xmax = end, ymin = 4, ymax = 5, fill = label)) + 
  geom_rect() + 
  #geom_segment(data = subset(filld, label %in% label[duplicated(label)]),
  #             aes(x = p, y = 0, xend = p, yend = 4, colour = label),
  #             size = 2, show_guide = FALSE) +
  geom_text(aes(x = p, y = 5.5, label = label), colour = "black", size=2, angle=39) +
  coord_polar() + 
  scale_y_continuous(limits = c(0, 6))
'''


library(ggplot2)
fold_changes <- c(rnorm(20000, 0, 2))
pvalues <- runif(n=20000, min=1e-50, max=.1)
dif <- data.frame(fc =fold_changes,pv =pvalues)
dif$thershold <- ifelse(dif$fc > 1 & dif$pv < 0.01, "red", 
                        ifelse(dif$fc < -1 & dif$pv < 0.01, -1, "blue"))
ggplot(data=dif, aes(x=fc, y=-log10(pv))) +
  geom_point( size=1 ,aes(color=as.factor(thershold))) +
  theme(legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")  + theme_bw()+
  annotate("label", x =c(-8,5), y = 4.75, label = c("400","120"), col=c("red","steelblue"))+
  annotate("text", x =c(-8,5), y = 5, label = c("50 FC>4","8FC <-4"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),)
        
library(calibrate)
# Volcano plot 
res = read.csv('/home/trangle/Desktop/Covid19project/Foldchange_meanintensity_includingHPA057697.csv', header=TRUE)
res$log2FoldChange = log2(res$fc)
res$padj = p.adjust(res$pval, method = 'BH', n = length(res$pval))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.8))
