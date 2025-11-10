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
# mapping <- expression(
#   '591661_assembly'="ST147"['NDM-9'],
#   'TLS30_assembly'="ST147"["NDM-1"]*"c.i.1",
#   '591662_assembly'="ST147"['NDM-1']*"c.i.2",
#   '591663_assembly'="ST13"["OXA-48"],
#   'ST258_KP35_assembly'="ST258",
#   'ST307_NDM5'="ST307"["NDM-5"],
#   'ATCC43816'="ST493",
#   '591665_assembly'="ST512"["VIM"]
# );

mapping <- expression(
  'ERR1023711'='MSB1_8A',
  'ERR1023848'='INF211',
  '591661_assembly'="ST147 NDM-9",
  'TLS30_assembly'="ST147 NDM-1 c.i.1",
  '591662_assembly'="ST147 NDM-1 c.i.2",
  '591663_assembly'="ST13 OXA-48",
  'ST258_KP35_assembly'="ST258",
  'ST307_NDM5'="ST307 NDM-5",
  'ATCC43816'="ST493",
  '591665_assembly'="ST512 VIM"
);

origin <- list(
  '166H3'='Hopital De Bicetre',
  '183J4'='Hopital De Bicetre',
  '189D7'='Hopital De Bicetre',
  '192I2'='Hopital De Bicetre',
  '198J6'='Hopital De Bicetre',
  '240F8'='Hopital De Bicetre',
  '242J2'='Hopital De Bicetre',
  '247C4'='Hopital De Bicetre',
  '265C5'='Hopital De Bicetre',
  '333H7'='Hopital De Bicetre',
  '361J6'='Hopital De Bicetre',
  '363B4'='Hopital De Bicetre',
  '363I3'='Hopital De Bicetre',
  '368J3'='Hopital De Bicetre',
  '369D5'='Hopital De Bicetre',
  '412H6'='Hopital De Bicetre',
  '591661_assembly'='University Hospital of Pisa',
  '591662_assembly'='University Hospital of Pisa',
  '591663_assembly'='University Hospital of Pisa',
  '591665_assembly'='University Hospital of Pisa',
  'ATCC43816'='ATCC',
  'ERR1023711'='Monash University',
  'ERR1023848'='Monash University',
  'INF049'='Monash University',
  'INF281'='Monash University',
  'INF298'='Monash University',
  'INF299'='Monash University',
  'INF305'='Monash University',
  'INF310'='Monash University',
  'INF357'='Monash University',
  'ST258_KP35_assembly'='University Hospital of Pisa',
  'ST307_NDM5'='University Hospital of Pisa',
  'TLS30_assembly'='University Hospital of Pisa'
);

map_origin <- function(nameList){
  result <- c();
  for (i in 1:length(nameList)) {

    # Check if the column_name is in the mapping
    name <- nameList[[i]];
    # label <- name;
    label <- 'other';
    if (name %in% names(origin)) {
      # If present, use the corresponding value as the new column name
      label <- origin[[name]];
      label <- unlist(label);
    }

    result <- c(result, label);
  }
  
  return (result);
}

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
  mutate(label_style=rename_assembly(label)) %>%
  mutate(origin=map_origin(label));
print(data.tree)
rownames(data.tree) <- tree.tls$tip.label;


# Open a PDF for plotting; units are inches by default
cat("Rendering data..", "\n");

# Generate palettes
palette.ST.label <- unique(data.tree$ST);
palette.ST.label <- palette.ST.label[! palette.ST.label %in% c('unknown')];
palette.ST.color <- colorRampPalette(wes_palette(name="FantasticFox1"))(length(palette.ST.label));
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

color.cluster <- c(
  "Hopital De Bicetre"='#2A2D34', 
  "University Hospital of Pisa"='#009DDC',
  'Monash University'='#F26430',
  'ATCC'='#ACAD94'
);

# Plot phylogenetic tree
p.tree <- ggtree(tree.tls, ladderize=F) %<+% data.tree +
  scale_x_continuous() +
  geom_tiplab(
    align=TRUE, 
    size=2.5,
    aes(label=label_style), 
    parse = FALSE,
    offset=0.0005
  ) +
  geom_tippoint(aes(colour=factor(origin)), size=1.5) +
  scale_color_manual(
    values=color.cluster
    # labels = c("Cluster 1", "Cluster 2")
  ) +
  guides(color=guide_legend("Origin")) +
  theme_tree2() + 
  theme(legend.key.size = unit(0.4, 'cm'));

# Generate a list of plots
plotList <- list();

# Append kleborate data
tree.tls.max <- max(tree.tls$edge.length);
heatmap.offset <- tree.tls.max*0.8+0.005;
heatmap.item <- tree.tls.max/4;
heatmap.width <- 0.12;

# ST
p.tree <- p.tree + new_scale_fill();
p.tree <- p.tree %>% gheatmap(
    data.tree %>% select(ST),
    colnames_position = "top",
    custom_column_labels = c('ST'),
    colnames_angle=90, 
    font.size = 2.5,
    offset=heatmap.offset,
    # width=0.2,
    width=heatmap.width,
    hjust=0
  )+ 
  scale_fill_manual(
    values = palette.ST.color, 
    na.translate = FALSE,
    name = "ST"
  ) + theme(
    legend.text = element_text(size = 8),  # Reduce legend text size
    legend.key.size = unit(0.3, 'cm')  # Reduce legend item size
  );

# K-locus
p.tree <- p.tree + new_scale_fill();
p.tree <- p.tree %>% gheatmap(
    data.tree %>% select(K_locus),
    colnames_position = "top",
    custom_column_labels = c('KL'),
    colnames_angle=90, 
    font.size = 2.5,
    offset=heatmap.offset + (heatmap.item*1),
    # width=0.2,
    width=heatmap.width,
    hjust=0
  ) + scale_fill_manual(
    values = palette.klocus.color, 
    na.translate = FALSE,
    name = "K-locus"
  ) + theme(
    legend.text = element_text(size = 8),  # Reduce legend text size
    legend.key.size = unit(0.3, 'cm')  # Reduce legend item size
  )

# Store results
p.tree <- p.tree + ggtree::vexpand(.05, 1) + theme(plot.margin = margin(t=2, r=2, b=2, l=2, unit="mm"));
plotList[['tree']] <- p.tree;

# Save the final plot
ggsave(filename=opt$output, marrangeGrob(top = NULL, grobs = plotList, nrow=1, ncol=1), units="mm", width=120, height=130, dpi=300);

# Often the most useful way to look at many warnings:
summary(warnings())