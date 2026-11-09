# Install/load required packages

library(ape)
library(ggtree)
library(ggplot2)

newick <- "(Capsicum_annuum_11531:0.0341657,(Lycium_barbarum_9331:0.0663837,(Solanum_tuberosum_5289:0.0169974,(Solanum_lycopersicum_4597:0.0115913,Solanum_pennellii_2375:0.00634721)791:0.0195411)1196:0.0399495)942:0.0341657)1918;"

tree <- read.tree(text = newick)

parse_label <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  list(
    count = as.numeric(tail(parts, 1)),
    name  = paste(head(parts, -1), collapse = "_")
  )
}

tip_parsed  <- lapply(tree$tip.label,  parse_label)
node_parsed <- lapply(tree$node.label, parse_label)

tree$tip.label  <- sapply(tip_parsed,  `[[`, "name")
tree$node.label <- sapply(node_parsed, `[[`, "name")

n_tips <- length(tree$tip.label)
anno_df <- data.frame(
  node = 1:(n_tips + length(tree$node.label)),
  dups = c(
    sapply(tip_parsed,  `[[`, "count"),
    sapply(node_parsed, `[[`, "count")
  )
)

p <- ggtree(tree) %<+% anno_df +
  
  # Clean species names
  geom_tiplab(size = 3.5, offset = 0.001) +
  
  # Duplication counts on tips, positioned after the label
  geom_tiplab(aes(label = paste0("(", dups, ")")),
              size = 3, offset = 0.031, colour = "steelblue") +
  
  # Internal node labels with duplication counts
  geom_nodelab(
    aes(label = paste0(label, " (", dups, ")")),
    size   = 3,
    hjust  = 1.1,
    vjust  = -0.3,
    colour = "firebrick"
  ) +
  
  geom_treescale() +
  xlim(0, 0.22) +
  labs(
    title    = "Gene Family Duplications",
    subtitle = "Blue = species-specific dups | Red = node-level dups"
  ) +
  theme_tree2()

p

tree <- read.tree(text = newick)
plot(tree, show.tip.label = TRUE, label.offset = 0.002)
nodelabels(tree$node.label, cex = 0.7, bg = "lightyellow")


newick <- "(Capsicum_annuum:0.0341657,(Lycium_barbarum:0.0663837,(Solanum_tuberosum:0.0169974,(Solanum_lycopersicum:0.0115913,Solanum_pennellii:0.00634721)791:0.0195411)1196:0.0399495)942:0.0341657)1918;"

tree <- read.tree(text = newick)

# Define sensitive species
sensitive <- c("Solanum_lycopersicum_4597", "Capsicum_annuum_11531")

# Assign colors to tip labels based on tolerance
tip_colors <- ifelse(tree$tip.label %in% sensitive, "red", "cornflowerblue")

# Plot
plot(tree, 
     show.tip.label = TRUE, 
     label.offset   = 0.002,
     tip.color      = tip_colors)

nodelabels(tree$node.label, cex = 0.7, bg = "lightyellow")

# Add a legend
legend("bottomleft", 
       legend = c("Sensitive", "Tolerant"),
       text.col = c("red", "cornflowerblue"),
       bty = "n",
       cex = 0.8)




######
#plot unrooted species tree


library(ape)


tr <- read.tree("~/Documents/uni/Orthofinder_8_29/final-results/WorkingDirectory/SpeciesTree_unrooted.txt")

library(ape)

tolerant <- c("Solanum_pennellii", "Solanum_tuberosum", "Lycium_barbarum")
cols <- ifelse(tr$tip.label %in% tolerant, "darkgreen", "red")

## first pass: only compute the layout to get the coordinate range
plot(tr, type = "unrooted", lab4ut = "axial", plot = FALSE)
lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
xr <- range(lp$xx); yr <- range(lp$yy)
padx <- 2.5 * diff(xr)     # room for the labels; increase if still cut off
pady <- 0.3 * diff(yr)

pdf("figures/species_tree_unrooted.pdf", width = 5, height = 2)
par(mar = c(1, 1, 1, 1))
plot(tr, type = "unrooted", lab4ut = "horizontal", rotate.tree = 90,
     tip.color = cols, font = 3, cex = 0.9, edge.width = 2,
     x.lim = xr + c(-padx, padx), y.lim = yr + c(-pady, pady))

usr <- par("usr")                                  # current plot limits
add.scale.bar(x = usr[1] + 0.03 * diff(usr[1:2]),  # bottom-left corner
              y = usr[3] + 0.06 * diff(usr[3:4]),
              length = 0.02, lwd = 2, cex = 0.8)
dev.off()
