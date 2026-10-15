# load required packages
library(ComplexHeatmap)
library(circlize)
library(grid)

# set working directory
setwd(here("plots"))

# create labels containing correlation coefficients
heatmap_text <- matrix(sprintf("%.2f", module.trait.corr),
                       nrow = nrow(module.trait.corr),
                       ncol = ncol(module.trait.corr),
                       dimnames = dimnames(module.trait.corr))

# add FDR significance stars
heatmap_text[module.trait.corr.FDR < .05] <- paste0(sprintf("%.2f", module.trait.corr[module.trait.corr.FDR < .05]), "\n*")
heatmap_text[module.trait.corr.FDR < .01] <- paste0(sprintf("%.2f", module.trait.corr[module.trait.corr.FDR < .01]), "\n**")
heatmap_text[module.trait.corr.FDR < .001] <- paste0(sprintf("%.2f", module.trait.corr[module.trait.corr.FDR < .001]), "\n***")

# plot heatmap of module-trait correlations
png("fig4.10_WGCNA.png", height = 2, width = 6.27, units = "in", res = 600) 
Heatmap(t(module.trait.corr),
        name = "Pearson r",
        col = colorRamp2(c(-1, -0.3, 0, 0.3, 1), c("#2674B8", "#8FB8DA", "white", "#EA9AA5", "#D9364A")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(t(heatmap_text)[i, j],
                    x, y, gp = gpar(fontsize = 7))
  },
  
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8),
  
  column_title = "Module–Trait Relationships",
  column_title_gp = gpar(fontsize = 8),
  
  heatmap_legend_param = list(title = "Pearson r",
                              at = c(-1, -0.5, 0, 0.5, 1),
                              title_gp = gpar(fontsize = 7),
                              labels_gp = gpar(fontsize = 7)))
dev.off()  
