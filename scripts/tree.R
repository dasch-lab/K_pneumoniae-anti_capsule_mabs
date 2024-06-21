#!/usr/bin/Rscript

# Set some variables
align_text_size <- 2;
legend_label_size <- 10;
border.color <- 'white';
border.size <- 0.30;

# options(
#     BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = FALSE
# )

usePackage <- function(package){
  if(eval(parse(text=paste("require(",package,")")))) {
    return (TRUE)
  }

  install.packages(package, repos='http://cran.us.r-project.org')
  return (eval(parse(text=paste("require(",package,")"))))
}

usePackage.bio <- function(package){

  library.path <- .libPaths()[1];

  if(!is.element(package, installed.packages(lib.loc=library.path)[,1]))
    BiocManager::install(package, lib=library.path, update=FALSE)

  # BiocManager::install(package, update=FALSE)
  # install.packages(package, repos='http://cran.us.r-project.org')
  return (eval(parse(text=paste("require(",package,")"))))
}

# For command line parsing
usePackage("optparse"); 

# This is the options list
options = list(
  make_option(c("-t", "--tree"), type="character", default=NULL, help="input tree file", metavar="character"),
  make_option(c("-k", "--kleborate"), type="character", default=NULL, help="input kleborate file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output file name", metavar="character")
); 

# Parsing options
parser = OptionParser(option_list=options, description='Generate the TLS tree figure');
opt = parse_args(parser);

# Check options
if (is.null(opt$tree) || is.null(opt$output) || is.null(opt$kleborate)){
  print_help(parser)
  stop("Missing argument(s)!", call.=FALSE)
}

# Install packages
usePackage("ggplot2");
usePackage("BiocManager");
usePackage("cowplot");  # For plot_gridß
usePackage("ggnewscale");
usePackage('gridExtra');
usePackage("castor");
usePackage("tibble");
usePackage("ape");
usePackage("dplyr");
usePackage("tidyr");
usePackage("ggrepel");

# Colors 
usePackage("RColorBrewer");
usePackage("wesanderson");

# Install bioconductor packages
usePackage.bio('ggtree');
usePackage.bio('treeio');

# Define order and labels
mapping <- expression(
  '591661_assembly'="ST147"['NDM-9'],
  'TLS30_assembly'="ST147"["NDM-1"]*"c.i.1",
  '591662_assembly'="ST147"['NDM-1']*"c.i.2",
  '591663_assembly'="ST13"["OXA-48"],
  'ST258_KP35_assembly'="ST258",
  'ST307_NDM5'="ST307"["NDM-5"],
  'ATCC43816'="ST493",
  '591665_assembly'="ST512"["VIM"]
);

rename_assembly <- function(nameList){

  result <- c();
  for (i in 1:length(nameList)) {

    # Check if the column_name is in the mapping
    name <- nameList[[i]];
    label <- name;
    if (name %in% names(mapping)) {
      # If present, use the corresponding value as the new column name
      label <- mapping[[name]];
      label <- unlist(label);
    }
    
    result <- c(result, label);
  }
  
  return (result);
}

# Load trees
tree.tls <- read_tree(file=opt$tree, interpret_quotes=TRUE, edge_order="pruningwise");

# Reroot tree
tree.tls  <- root(tree.tls, outgroup = "ATCC43816", edgelabel = TRUE)

# Create a vector of factors with node support and hide all nodes with support >= 50
node.label <- c(rep(NA, length(tree.tls$tip.label)), as.numeric(tree.tls$node.label));
node.label <- unlist(lapply(node.label, function(x){ return(ifelse(is.na(x) || x >= 50, 0, 1)) }));
node.support <- factor(node.label);

# Load kleborate data
data.kleborate <- read.csv(opt$kleborate, header=TRUE, sep="\t", dec=".", stringsAsFactors = FALSE);
data.kleborate <- data.kleborate %>%
  mutate(label=strain) %>%
  select(label, O_locus, K_locus, ST, O_type);
  
# Get sorted list of labels
data.tree <- data.frame(label=tree.tls$tip.label, stringsAsFactors = FALSE) %>%
  left_join(data.kleborate, by='label') %>%
  replace_na(list(ST = "unknown", O_locus = "unknown", K_locus = "unknown", O_type="unknown")) %>%
  select(label, O_locus, K_locus, ST, O_type) %>%
  mutate(label_style=rename_assembly(label));
print(data.tree);
rownames(data.tree) <- tree.tls$tip.label;




# Open a PDF for plotting; units are inches by default
cat("Rendering data..", "\n");

# Generate palettes
palette.ST.label <- unique(data.tree$ST);
palette.ST.label <- palette.ST.label[! palette.ST.label %in% c('unknown')];
palette.ST.color <- colorRampPalette(wes_palette(name="GrandBudapest1"))(length(palette.ST.label));
palette.ST.color <- c(palette.ST.color, '#CECECE');
palette.type.label <- unique(data.tree$type);
palette.type.label <- palette.type.label[! palette.type.label %in% c('unknown')];
palette.type.color <- colorRampPalette(wes_palette(name="Zissou1"))(length(palette.type.label));
palette.type.color <- c(palette.type.color, '#CECECE');
palette.klocus.label <- unique(data.tree$K_locus);
palette.klocus.label <- palette.klocus.label[! palette.klocus.label %in% c('unknown')];
palette.klocus.color <- colorRampPalette(wes_palette(name="Rushmore1"))(length(palette.klocus.label));
palette.klocus.color <- c(palette.klocus.color, '#CECECE');
palette.klocus.label <- c(palette.klocus.label, 'unknown');

# Plot phylogenetic tree
p.tree <- ggtree(tree.tls, ladderize=F) %<+% data.tree +
  scale_x_continuous() +
  geom_tiplab(
    align=TRUE, 
    size=2,
    aes(label=label_style), parse = TRUE
  ) +
  theme_tree2() + 
  theme(legend.key.size = unit(0.4, 'cm'));

# Generate a list of plots
plotList <- list();

# Append kleborate data
tree.tls.max <- max(tree.tls$edge.length);
heatmap.offset <- tree.tls.max*1.2;
heatmap.item <- tree.tls.max/3;

# ST
p.tree <- p.tree %>% gheatmap(
    data.tree %>% select(ST),
    colnames_position = "top",
    custom_column_labels = c('ST'),
    colnames_angle=90, 
    font.size = 2.5,
    offset=heatmap.offset,
    width=0.3,
    hjust=0
  ) + 
  scale_fill_manual(
    values = palette.ST.color, 
    na.translate = FALSE,
    name = "ST",
    guide = guide_legend(order = 2)
  ) + 
  new_scale_fill();


# O-type
p.tree <- p.tree %>% gheatmap(
    data.tree %>% select(O_type),
    colnames_position = "top",
    custom_column_labels = c('O-type'),
    colnames_angle=90, 
    font.size = 2.5,
    offset=heatmap.offset + (heatmap.item*2),
    width=0.3,
    hjust=0
  ) + scale_fill_manual(
    values = palette.klocus.color, 
    na.translate = FALSE,
    name = "O-type",
    guide = guide_legend(order = 4)
  ) +
  new_scale_fill();

# K-locus
p.tree <- p.tree %>% gheatmap(
    data.tree %>% select(K_locus),
    colnames_position = "top",
    custom_column_labels = c('K-locus'),
    colnames_angle=90, 
    font.size = 2.5,
    offset=heatmap.offset + (heatmap.item*4),
    width=0.3,
    hjust=0
  ) + scale_fill_manual(
    values = palette.klocus.color, 
    na.translate = FALSE,
    name = "K-locus",
    guide = guide_legend(order = 5)
  ) +
  new_scale_fill();

# Store results
p.tree <- p.tree + ggtree::vexpand(.05, 1);
plotList[['tree']] <- p.tree;

# Save the final plot
ggsave(filename=opt$output, marrangeGrob(top = NULL, grobs = plotList, nrow=1, ncol=1), units="mm", width=85, height=100, dpi=300);

# Often the most useful way to look at many warnings:
summary(warnings())