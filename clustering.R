library(tidyverse)
library(patchwork)
library(ggdendro)
library(cowplot)

# Generate some fake data. 
# 100 genes, 4 time points

genes <- 100
stuff <- matrix(rnorm(4*genes), nrow = genes)

colnames(stuff) <- sprintf("TP%s", 1:4)
rownames(stuff) <- sprintf("Gene_%s", 1:genes)

image(t(stuff))


# Make a super simple heatmap


pdf(file = "base-method.pdf", height = 14, width =8)
heatmap(stuff, Colv = NA)
dev.off()


# If you want to get fancy and have more control:

# Cluster the data
dat_clustered <- dist(stuff, method = "euclidean")


# Hierarchical clusters and then extract the order.
dat_clust <- hclust( dat_clustered, method = "ward.D" )

plot(dat_clust)

# Now advanced plotting ------------------------------------------
# gently modified from https://stackoverflow.com/a/42379926/7023826
# Turn your matrix into a dataframe or easier manipulation
sample_names <- colnames(stuff)

# Obtain the dendrogram
dend <- as.dendrogram(dat_clust)

dend_data <- dendro_data(dend)

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data <- with(
    segment(dend_data), 
    data.frame(x = y, y = x, xend = yend, yend = xend))
# Use the dendrogram label data to position the gene labels
gene_pos_table <- with(
    dend_data$labels, 
    data.frame(y_center = x, gene = as.character(label), height = 1))

# Table to position the samples
sample_pos_table <- data.frame(sample = sample_names)|>
    mutate(x_center = (1:n()), 
           width = 1)

# Neglecting the gap parameters
heatmap_data <- stuff|> 
    as.data.frame() |>
    rownames_to_column("gene") |>
    gather(sample, expr,-gene)|>
    left_join(gene_pos_table)|>
    left_join(sample_pos_table)

# Limits for the vertical axes
gene_axis_limits <- with(
    gene_pos_table, 
    c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 
    0.1 * c(-1, 1) # extra spacing: 0.1

# Heatmap plot
plt_hmap <- ggplot(heatmap_data, 
                   aes(x = x_center, y = y_center, fill = expr, 
                       height = height, width = width)) + 
    geom_tile() +
    scale_fill_gradient2("expr", high = "darkred", low = "darkblue") +
    scale_x_continuous(breaks = sample_pos_table$x_center, 
                       labels = sample_pos_table$sample, 
                       expand = c(0, 0)) + 
    # For the y axis, alternatively set the labels as: gene_position_table$gene
    scale_y_continuous(breaks = gene_pos_table[, "y_center"], 
                       labels = rep("", nrow(gene_pos_table)),
                       limits = gene_axis_limits, 
                       expand = c(0, 0)) + 
    labs(x = "Sample", y = "") +
    theme_bw() +
    theme(# margin: top, right, bottom, and left
          plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
          panel.grid.minor = element_blank(), axis.ticks.y = element_blank())

# Dendrogram plot
plt_dendr <- ggplot(segment_data) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    scale_x_reverse(expand = c(0, 0.5)) + 
    scale_y_continuous(breaks = gene_pos_table$y_center, 
                       labels = gene_pos_table$gene, 
                       limits = gene_axis_limits, 
                       expand = c(0, 0)) + 
    labs(x = "Distance", y = "", colour = "", size = "") +
    theme_bw() + 
    theme(panel.grid.minor = element_blank())

library(patchwork)

final_fig <- plt_dendr+ plt_hmap +plot_layout(widths = c(.2,.8))


cowplot::ggsave2("ggplot-method.pdf", height = 14, width = 8.5)


