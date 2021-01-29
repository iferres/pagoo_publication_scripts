
#####################
## Set environment ##
#####################
#setwd(tempdir())

## Download and decompress data
tar_gz <- "cfetus_pangenome.tar.gz"
if ( !file.exists(tar_gz) ){
	download.file(url = "https://ndownloader.figshare.com/files/26144075", destfile = tar_gz)
}
untar(tarfile = tar_gz, exdir = "C_fetus")

set.seed(123)

##################
## Load dataset ##
##################

## List gff files.
gffs <- list.files(path = "C_fetus", 
                   pattern = "[.]gff$", 
                   recursive = TRUE, 
                   full.names = TRUE)

## List roary's gene_presence_absence.csv output file.
gpa <- list.files(path = "C_fetus/pangenome_cfettus/", 
                  pattern = "gene_presence_absence[.]csv", 
                  full.names = TRUE)

## Read organism's metadata.
org_meta <- read.csv("C_fetus/metadata_Cfettus_Iraola_2017.tsv", sep = "\t")
org_meta$org <- org_meta$Lane
org_meta$Lane <- NULL

## Load pagoo.
library(pagoo)

## Load pangenome and add metadata to organisms.
p <- roary_2_pagoo(gene_presence_absence_csv = gpa, gffs = gffs)
p$add_metadata("org", data = org_meta)
p$save_pangenomeRDS("C_fetus/pangenome.rds")
# p <- load_pangenomeRDS("pangenome.rds")

## Set core-level to 100%.
p$core_level <- 100
p$summary_stats
# DataFrame with 4 rows and 2 columns
#      Category    Number
#   <character> <integer>
# 1       Total      7326
# 2        Core      1117
# 3       Shell      5391
# 4       Cloud       818

###################
## Clean dataset ##
###################

## Load pipe.
library(magrittr)

## Barplot: number of CDS per organism.
pdf("C_fetus/barplot_Num_CDS_all.pdf")
p$pan_matrix %>%
  rowSums() %>%
  sort(decreasing = TRUE) %>% 
  barplot(main = "Number of CDS on each genome",
          sub = "The four genomes with abnormal number of CDS will be removed from the dataset")
dev.off()

# There're 4 organisms with an abnormal number of cds. Drop them from the
#  pangenome.
p$pan_matrix %>%
  rowSums() %>%
  sort(decreasing = TRUE) %>% 
  names() %>%
  magrittr::extract(1:4) %>%
  p$drop()

p$dropped
#           36           37           38           39 
# "16473_4#71" "16473_4#73" "16473_4#77" "16473_4#79" 

p$summary_stats
# DataFrame with 4 rows and 2 columns
#      Category    Number
#   <character> <integer>
# 1       Total      4686
# 2        Core      1164
# 3       Shell      2757
# 4       Cloud       765


##############################
## General statistics plots ##
##############################

library(ggplot2)
library(patchwork)

# 1. Pangenome curves
curves <- p$gg_curves() +                                     # Plot core- and pan-genome curves
  scale_color_manual(values = c('black', 'black')) +  # Customize line colors
  geom_point(alpha = .05, size = 4, color = 'grey') + # Add semi-transparent data points
  theme_bw(base_size = 15) +                          # Customize background theme
  theme(legend.position = 'none',                     # Remove legend
        axis.title = element_text(size = 12),                     # Customize axis title
        axis.text = element_text(size = 12))                      # Customize axis text size

# 2. Gene frequency bar plots
bars <- p$gg_barplot() +                                       # Plot gene frequency distribution
  theme_bw(base_size = 15) +                                   # Customize background color
  theme(axis.title = element_text(size = 12),                  # Customize axis label size
        axis.text=element_text(size = 12)) +                   # Customize axis text size
  geom_bar(stat = 'identity', color = 'black', fill = 'black') # Customize bar color and borders

# 3. PCA of accessory genes colored by host
pca <- p$gg_pca(colour = 'Host', size = 4) +                # Plot PCA, color by host
  theme_bw(base_size = 10) +                                # Customize background theme
  theme(legend.position = 'bottom') +                       # Customize legend position
  theme(axis.title = element_text(size = 12),               # Customize axis title
        axis.text = element_text(size = 12))                # Customize axis text size

# 4. Pie chart of core and accessory genes
pie <- p$gg_pie() +                                         # Plot pie chart
  theme_bw(base_size = 10) +                                # Customize background theme
  scale_fill_discrete(guide = guide_legend(keywidth = .75,
                                           keyheight = .75)) + # Customize fill
  scale_fill_brewer(palette = "Blues") +                    # Customize fill color
  scale_x_discrete(breaks = c(0, 25, 50, 75)) +             # Customize axis scales
  theme(legend.position = 'bottom',                         # Customize legend position
        legend.title = element_blank(),                     # Remove legend title
        legend.text = element_text(size = 10),              # Change legend text size
        legend.margin = margin(0, 0, 13, 0),                # Change legend margins
        legend.box.margin = margin(0, 0, 5, 0),             # Change box margins
        axis.title.x = element_blank(),                     # Remove X-axis title
        axis.title.y = element_blank(),                     # Remove Y-axis title
        axis.ticks = element_blank(),                       # Remove axis ticks
        axis.text.x = element_blank())                      # Remove X-axis text


# 5. Use patchwork to arrange plots using math operators
stats_plots <- (curves + bars) / (pca + pie)
ggsave(filename = "stats_plots.pdf", stats_plots)
ggsave(filename = "stats_plots.png", stats_plots)

#######################################
## Concatenated core-genome phylogeny #
#######################################

# Align translation of core genes.
library(DECIPHER)
library(parallel) 
ali <- p$core_seqs_4_phylo() %>%                      # Core genome sequences.
  mclapply(DECIPHER::AlignTranslation, mc.cores = 8)  # Align translation.

# Compute phylogeny
library(phangorn)
phy <- ali %>%                               # Select neutral clusters
  do.call(Biostrings::xscat, .) %>%          # Concatenate alignments
  setNames(p$organisms$org) %>%              # Set sequence names
  as('matrix') %>%                           # Transform to matrix
  phangorn::phyDat(type = 'DNA') %T>%        # Transform to phangorn's phyDat
  assign('dat', ., .GlobalEnv) %>%           # Assign to "dat" in .GlobalEnv
  phangorn::dist.ml() %>%                    # Compute distance
  phangorn::NJ() %>%                         # Compute NJ (initial tree)
  phangorn::pml(data = dat, k = 4) %>%       # Compute likelihood with 4 discrete
                                             # gamma distributions. 
  phangorn::optim.pml(rearrangement = "NNI", # Optimize topology with NNI,
                      optGamma = TRUE,       # optimize gamma rate parameter,
                      optInv = TRUE,         # optimize prop of variable size,
                      model ="GTR")          # and use "GTR" model.

# Create ggtree
library(ggplot2)
library(ggtree)
gtree <- phy$tree %>%             # Get phylo tree.
  midpoint() %>%                  # Root at midpoint.
  ladderize() %>%                 # Ladderize to make it prettier.
  ggtree::ggtree() %<+%           # Create ggtree object.
  as.data.frame(p$organisms)      # Add isolate metadata.



####################################################
## Identify genes that better reproduce phylogeny ##
####################################################

# Compute topological distance from each gene tree to the core genome tree.
dtopo <- ali %>%
  sapply(
    function(x) {
      x %>%
        setNames(p$organisms$org) %>%           # Set sequence names.
        as('matrix') %>%                        # Transform to matrix.
        phangorn::phyDat(type = 'DNA') %>%      # Transform to phangorn's phyDat.
        phangorn::dist.ml() %>%                 # Compute distance.
        phangorn::NJ() %>%                      # Compute NJ.
        ape::dist.topo(y = phy$tree)            # Compute topological distance.
    } 
  )

head(sort(dtopo), n = 7) # mcp4 is the closest tree
#   mcp4 group_1165       mnmC        rnr       addA group_1032      murC
#    277        279        279        281        283        283       283


# Compute a NJ tree of mcp4 (again).
mcp4 <- ali[["mcp4"]] %>%
  setNames(p$organisms$org) %>%           # Set sequence names.
  as('matrix') %>%                        # Transform to matrix.
  phangorn::phyDat(type = 'DNA') %>%      # Transform to phangorn's phyDat.
  phangorn::dist.ml() %>%                 # Compute distance.
  phangorn::NJ() %>%                      # Compute NJ tree.
  midpoint() %>%                          # Midpoint root.
  ladderize()                             # Ladderize.

# Set negative branch lentghs to 0
mcp4$edge.length[mcp4$edge.length < 0] <- 0

# Create a tanglegram to compare coregenome and mcp4 phylogenies
library(tidytree)
d1 <- fortify(as.treedata(gtree))          # To phylo and then to tibble.
d2 <- fortify(mcp4)                        # To tibble.
d2$x <- max(d2$x) - d2$x + max(d1$x) + .05 # Reverse mcp4 tree and set it apart.
pp <- gtree + geom_tree(data = d2)         # Create plot object for both trees.

# Lines
library(dplyr)
dd <- bind_rows(d1, d2) %>% # rbind both tibbles
  filter(!is.na(label))     # Remove NAs

# Plot tanglegram
tangle <- pp + geom_line(aes(x, y, group=label, color = Host), data=dd) + 
  ggtitle("Tanglegram", subtitle = "Coregenome phylogeny (left) vs mcp4 phylogeny (right)")

# Save
ggsave("C_fetus/tanglegram.pdf", plot = tangle)
ggsave("C_fetus/tanglegram.png", plot = tangle)


####################
## Infer Lineages ##
####################

library(rhierbaps)
rhb <- hierBAPS(snp.matrix = ali %>%                  # Runs the hierBAPS algorithm.
                  do.call(Biostrings::xscat, .) %>%   # Concatenate alignments.
                  setNames(p$organisms$org) %>%       # Set sequence names.
                  as('matrix') %>%                    # Transform to matrix.
                  tolower(),                          # Input matrix alignment.
                n.pops = 10,                          # Max number of subpopulations.
                max.depth = 3,                        # Max depth for hierarchical clustering.
                n.extra.rounds = 5)                   # Extra rounds to ensure convergence.


Lin <- rhb$partition.df
rownames(Lin) <- Lin$Isolate
Lin$Isolate <- NULL
Lin[] <- lapply(Lin, as.factor)

 
# > tail(Lin)
#            level 1 level 2 level 3
# 18048_2#22       1       1       1
# 18048_2#23       1       2       3
# 18048_2#24       1       4       8
# 18048_2#25       1       1       2
# 18048_2#26       1       6      11
# 18048_2#32       1       6      11


######################################################
## Principal Component Analysis on accessory genome ##
######################################################

# The first PC is able to find two groups which clearly distinct between a
# a group of bovine isolates, and another with various hosts:
pcaplot <- p$gg_pca(colour = "Host")
ggsave("C_fetus/PCA.pdf", pcaplot)
ggsave("C_fetus/PCA.png", pcaplot)

# Compute PCA and plot density of cluster loadings in PC1.
pca <- p$pan_pca()# Compute prcomp().
pdf("C_fetus/PC1_loading_density.pdf")
plot(density(pca$rotation[, 1]),          # Plot cluster variance loadings on the
     main = "PC1 Loading Density",
     sub = "Clusters with loading less than -0.05 and grater than 0.05 contribute to variance the most.")          
                                          # principal component.
abline(v = c(-0.05, 0.05), col= "red")    # Most of the clusters have low loadings,
                                          # i.e. near 0. We need those with high
                                          # absolute loading ( -0.05 > x > 0.05 ).
dev.off()

# Which clusters have high absolute loading ?
hload <- which( abs( pca$rotation[,1]) > 0.05 )
#      barS1     cirA_2       dapH     dctA_2     dpnM_2       exsA     fabG_1     fabG_2      fmt_2        glf group_1024 group_1026 
#         80        138        194        200        234        256        260        261        314        348        412        414 
#  group_103 group_1149 group_1153 group_1282 group_1350 group_1359 group_1371 group_1565  group_173 group_1760 group_1798 group_1799 
#        417        523        527        647        707        712        722        878        989       1006       1031       1032 
# group_1859 group_1923 group_2327 group_2328 group_2508 group_2880 group_3321 group_3323 group_3325 group_3336 group_3361 group_3371 
#       1074       1137       1376       1377       1451       1652       2045       2046       2047       2056       2068       2076 
# group_3372 group_3390  group_446  group_447  group_449  group_454  group_455   group_46  group_472  group_475  group_476  group_482 
#       2077       2091       2788       2789       2791       2796       2797       2799       2811       2814       2815       2822 
#  group_505  group_509  group_512  group_513  group_539  group_543  group_544  group_593  group_625  group_626  group_627  group_643 
#       2844       2848       2852       2853       2874       2878       2879       2916       2950       2951       2952       2964 
# group_7159  group_746  group_764  group_766  group_830  group_832  group_909  group_978  group_986     hldD_1     lexA_2     lexA_3 
#       3589       3656       3672       3674       3726       3728       3794       3853       3859       3912       3994       3995 
#     moaA_1     moaA_2       pctA       rfbE     trmR_2 
#       4074       4075       4208       4341       4560

###########################################
## Plot a heatmap associated with a tree ##
###########################################

# ggnewscale needed to add independent heatmap to the figure.
library(ggnewscale)

# First heatmap with just the Host name.
gh1 <- gheatmap(gtree, 
                data.frame(row.names = p$organisms$org, 
                           Host=p$organisms$Host), 
                colnames = FALSE,
                color = NULL,
                width = .1, 
                offset = 0) + 
  scale_fill_brewer(name = "Host", 
                    palette = "Accent", 
                    guide = guide_legend(order = 2, nrow = 3)) +
  new_scale_fill() 

# Second:Rhierbaps Lineages (Depth 1).
gh2 <- gheatmap(gh1, 
                Lin[, 1, drop = FALSE], 
                colnames = FALSE,
                color = NULL,
                width = .1, 
                offset = 0.006) + 
  scale_fill_brewer(name = "Depth 1", 
                    palette = "Dark2", 
                    breaks = levels(Lin[, 1]),
                    guide = guide_legend(nrow = 2, order = 3)) +
  new_scale_fill() 

# Third: Rhierbaps Lineages (Depth 2).
gh3 <- gheatmap(gh2, 
                Lin[, 2, drop = FALSE], 
                colnames = FALSE,
                color = NULL,
                width = .1, 
                offset = 0.011) + 
  scale_fill_brewer(name = "Depth 2", 
                    palette = "Set3", 
                    breaks = levels(Lin[, 2]),
                    guide = guide_legend(order = 4, nrow = 3)) +
  new_scale_fill()

# Fourth: Rhierbaps Lineages (Depth 3).
gh4 <- gheatmap(gh3, 
                Lin[, 3, drop = FALSE], 
                colnames = FALSE,
                color = NULL,
                width = .1, 
                offset = 0.016) + 
  scale_fill_discrete(name = "Depth 3",
                      breaks = levels(Lin[, 3]),
                      guide = guide_legend(order = 5, nrow = 3)) + 
  new_scale_fill()


# Order columns for prettier visualization.
ord <- hclust(dist(t(p$pan_matrix[, names(hload) ])))$order

# Last heatmap: gene counts for selected clusters.
gh5 <- gheatmap(gh4, p$pan_matrix[, names(hload)[ord] ], 
                offset = 0.024, 
                width = 3, 
                color = NULL, 
                colnames_angle = 360-90, 
                colnames_offset_y = 0, 
                font.size = 2.5, 
                hjust = 0) + 
  ylim(-5, NA) + 
  scale_fill_continuous(name = "Count", 
                        breaks = c(0, 1, 3, 5, 7), 
                        trans = "sqrt")

# Arrange legend.
gct <- gh5 + theme(legend.position="top", legend.key.size = unit(10, "points"))

# Save.
ggsave("C_fetus/Gene_content_tree.pdf", gct, height = 8)
ggsave("C_fetus/Gene_content_tree.png", gct, height = 8)



