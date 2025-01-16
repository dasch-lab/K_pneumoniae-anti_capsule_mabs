#!/usr/bin/Rscript

# Script for generating virulence heatmap
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
  return (eval(parse(text=paste("require(",package,")"))))
}

# For command line parsing
usePackage("argparse");

# create parser object
parser <- ArgumentParser()
parser$add_argument("-i", "--input", default=NULL, type="character", help="input igblast data");
parser$add_argument("-o", "--output", default=NULL, type="character", help="output pdf file");

opt <- parser$parse_args();

# Check options
if (is.null(opt$input) || is.null(opt$output)){
  # print_help(parser)
  stop("Missing argument(s)!", call.=FALSE)
}

# Get the directory of the path
output_dir <- dirname(opt$output);

# Install packages
usePackage("ggplot2");
usePackage("tidyr");
usePackage('gridExtra');
usePackage('ggnewscale');
usePackage("dplyr");
usePackage('plyr');
usePackage("cowplot");  # For plot_grid
usePackage('viridis');
usePackage("Biostrings");
usePackage("ggrepel");

usePackage.bio('ggtree');
usePackage.bio('treeio');

# Define clusters
cluster_map <- c(
  '08D18'=1,
  '08O09'=1,
  '03F18'=1,
  '08H10'=1,
  '08N23'=1,
  '08K19'=1,
  '08M13'=1,
  '03L02'=1,
  '08J10'=1,
  '10H18'=1,
  '08O03'=1,
  '05K07'=1,
  '09I10'=1,
  '05N02'=2,
  '05B17'=2,
  '08F04'=2,
  '05D08'=2,
  '05D14'=2,
  '05C11'=2,
  '05M13'=2
);

set_cluster <- function(name_list){
  result <- c();
  for(name in name_list){
    cluster_idx <- -1;
    if (name %in% names(cluster_map)) {
      cluster_idx <- cluster_map[[name]];
      cluster_idx <- unlist(cluster_idx);
    }
    result <- c(result, cluster_idx);
  }
  
  return (result);
}

parse_name <- function(name_list){
  result <- c();
  for(name in name_list){
    if(nchar(name) > 0){
      tokens <- unlist(strsplit(name, split = "_",fixed=TRUE));
      name <- head(tokens, n=1);
      name <- gsub("SBJ", "", name, fixed = TRUE);
      name <- gsub("-", "", name, fixed = TRUE);
    }
    result <- c(result, name);
  }
  
  return (result);
}

parse_germline <- function(germline_list){
  result <- c();
  for(germline in germline_list){
    if(nchar(germline) > 0){
      tokens <- unlist(strsplit(germline, split = "*",fixed=TRUE));
      germline <- head(tokens, n=1);
    }
    result <- c(result, germline);
  }
  
  return (result);
}

calculate_identity <- function(sequence_list, germline_list){
  result <- c();
  for (i in 1:length(sequence_list)) {

    sequence <- sequence_list[i];
    germline <- germline_list[i];
    
    total <- 0;
    count <- 0;
    for(j in 1:nchar(sequence)){
      char_s <- substr(sequence, j, j);
      char_g <- substr(germline, j, j);
      if(char_g == 'X'){
        next;
      }

      if (char_s == char_g){
        count <- count + 1;
      }

      total <- total + 1;
    }

    result <- c(result, 100*(count/total));
  }
  
  return (result);
}

# Function to write a dataframe to a FASTA file
write_fasta <- function(df, file_path) {
  # Open the file for writing
  con <- file(file_path, "w")

  # Write sequences to the file in FASTA format
  for (i in 1:nrow(df)) {
    fasta_name <- trimws(df$name[i]);
    cat(">", fasta_name, "\n", df$sequence[i], "\n", file = con, sep = "")
  }

  # Close the file
  close(con)
}

calculate_distance_matrix <- function(names, sequences) {

  BY.Norm.SW.SS = function(AAstring1, AAstring2){
  
    metric <- "BLOSUM62";
    # metric <- gonnet_matrix();
    gap_ppening = 10;
    gap_extension = 0.2;
    # gap_extension = 0.5;
    SW12=pairwiseAlignment(pattern = AAstring1,
                          subject = AAstring2,
                          type="local",
                          substitutionMatrix=metric,
                          gapOpening=gap_ppening,
                          gapExtension=gap_extension,
                          scoreOnly=TRUE)
    
    SW11=pairwiseAlignment(pattern = AAstring1,
                          subject = AAstring1,
                          type="local",
                          substitutionMatrix=metric,
                          gapOpening=gap_ppening,
                          gapExtension=gap_extension,
                          scoreOnly=TRUE)
    
    SW22=pairwiseAlignment(pattern = AAstring2,
                          subject = AAstring2,
                          type="local",
                          substitutionMatrix=metric,
                          gapOpening=gap_ppening,
                          gapExtension=gap_extension,
                          scoreOnly=TRUE)
    
    return(SW12/(sqrt(SW11)*sqrt(SW22)))
    
  }

  # Number of sequences
  n <- length(sequences)

  # Initialize an empty distance matrix
  distance_matrix <- matrix(NA, nrow = n, ncol = n)

  # Nested for loops to fill in the distance matrix
  for (i in 1:n) {
    for (j in i:n) {
      
      distance <- 1.0 - BY.Norm.SW.SS(sequences[[i]], sequences[[j]]);
      distance_matrix[i, j] <- distance
      distance_matrix[j, i] <- distance
    }
  }

  colnames(distance_matrix) <- names;
  rownames(distance_matrix) <- names;

  return(distance_matrix)
}

# This function extends the geom_violin to create a violin plot that is splitted by a grouping variable
GeomSplitViolin <- ggplot2::ggproto(
    "GeomSplitViolin",
    ggplot2::GeomViolin,
    draw_group = function(self,
                          data,
                          ...,
                          # add the nudge here
                          nudge = 0,
                          draw_quantiles = NULL) {
        data <- transform(data,
                          xminv = x - violinwidth * (x - xmin),
                          xmaxv = x + violinwidth * (xmax - x))
        grp <- data[1, "group"]
        newdata <- plyr::arrange(transform(data,
                                           x = if (grp %% 2 == 1) xminv else xmaxv),
                                 if (grp %% 2 == 1) y else -y)
        newdata <- rbind(newdata[1, ],
                         newdata,
                         newdata[nrow(newdata), ],
                         newdata[1, ])
        newdata[c(1, nrow(newdata)-1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

        # now nudge them apart
        newdata$x <- ifelse(newdata$group %% 2 == 1,
                            newdata$x - nudge,
                            newdata$x + nudge)

        if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {

            stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))

            quantiles <- ggplot2:::create_quantile_segment_frame(data,
                                                             draw_quantiles)
            aesthetics <- data[rep(1, nrow(quantiles)),
                               setdiff(names(data), c("x", "y")),
                               drop = FALSE]
            aesthetics$alpha <- rep(1, nrow(quantiles))
            both <- cbind(quantiles, aesthetics)
            quantile_grob <- ggplot2::GeomPath$draw_panel(both, ...)
            ggplot2:::ggname("geom_split_violin",
                             grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata, ...),
                                            quantile_grob))
        }
    else {
            ggplot2:::ggname("geom_split_violin",
                             ggplot2::GeomPolygon$draw_panel(newdata, ...))
        }
    }
)

geom_split_violin <- function(mapping = NULL,
                              data = NULL,
                              stat = "ydensity",
                              position = "identity",
                              # nudge param here
                              nudge = 0,
                              ...,
                              draw_quantiles = NULL,
                              trim = TRUE,
                              scale = "area",
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {

    ggplot2::layer(data = data,
                   mapping = mapping,
                   stat = stat,
                   geom = GeomSplitViolin,
                   position = position,
                   show.legend = show.legend,
                   inherit.aes = inherit.aes,
                   params = list(trim = trim,
                                 scale = scale,
                                 # don't forget the nudge
                                 nudge = nudge,
                                 draw_quantiles = draw_quantiles,
                                 na.rm = na.rm,
                                 ...))
}

# Parse the input data
data.input <- read.csv(opt$input, header=TRUE, sep="\t", row.names=NULL) %>%
  dplyr::mutate(name=parse_name(sequence_id)) %>%
  dplyr::mutate(v_call=parse_germline(v_call), j_call=parse_germline(j_call), d_call=parse_germline(d_call)) %>%
  dplyr::mutate(identity=calculate_identity(sequence_alignment_aa, germline_alignment_aa)) %>%
  dplyr::mutate(cdr3_len=nchar(cdr3_aa)) %>%
  dplyr::mutate(cluster=set_cluster(name)) %>%
  dplyr::mutate(type=ifelse(locus == 'IGH', 'heavy', 'light')) %>%
  dplyr::filter(cluster > 0) %>%
  dplyr::select(name, cluster, cdr3_aa, cdr3_len, identity, v_call, d_call, j_call, v_identity, j_identity, d_identity, type, sequence_alignment_aa, locus);

data.full.heavy <- data.input %>%
  dplyr::filter(type == 'heavy') %>%
  dplyr::rename(sequence = sequence_alignment_aa);
data.full.light <- data.input %>%
  dplyr::filter(type == 'light') %>%
  dplyr::rename(sequence = sequence_alignment_aa);

data.CDR3.heavy <- data.input %>%
  dplyr::filter(type == 'heavy') %>%
  dplyr::select(name, cdr3_aa) %>%
  dplyr::rename(sequence = cdr3_aa);
data.CDR3.light <- data.input %>%
  dplyr::filter(type == 'light') %>%
  dplyr::select(name, cdr3_aa) %>%
  dplyr::rename(sequence = cdr3_aa);

# Calculate the sequence distance
distance.heavy <- calculate_distance_matrix(data.full.heavy$name, data.full.heavy$sequence);
distance.light <- calculate_distance_matrix(data.full.light$name, data.full.light$sequence);

# Hierarchical clustering
cluster.heavy <- hclust(as.dist(distance.heavy));
cluster.light <- hclust(as.dist(distance.light));
order.heavy <- cluster.heavy$order;
order.light <- cluster.light$order;

# Create dendrogram objects
dendrogram.heavy <- as.dendrogram(cluster.heavy);
dendrogram.light <- as.dendrogram(cluster.light);

# Reorder the distance matrix based on the clustering result
distance.heavy <- distance.heavy[order.heavy, order.heavy];
distance.light <- distance.light[order.light, order.light];

# Convert the reordered distance matrix to a data frame
heatmap.heavy <- as.data.frame(as.table(distance.heavy));
heatmap.light <- as.data.frame(as.table(distance.light));

# Generate the frequency data for the bubble plot
data.frequency <- data.input %>%
  dplyr::group_by(type) %>%
  dplyr::mutate(total=n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(type, v_call, j_call, total, cluster) %>%
  dplyr::summarize(count=n()) %>%
  dplyr::mutate(freq=100*(count/total)) %>%
  dplyr::mutate(label=paste(v_call, j_call, sep=';')) %>%
  dplyr::mutate(cluster = as.factor(cluster));

data.violin <- data.input %>%
  dplyr::select(cluster, type, v_identity, d_identity, j_identity) %>%
  pivot_longer(cols = c(v_identity, d_identity, j_identity ), names_to = "key", values_to = "value") %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>%
  dplyr::filter(!is.na(value));

# Start generating plots
print("Rendering");

frequency.text.size <- 6;
frequency.title.size <- 8;

# Create a bubble chart with the recombinant frequency and 
germline.top.list <-  data.frequency %>%
  filter(type == "heavy") %>%
  ungroup() %>%
  top_n(3, freq) %>%
  bind_rows(
    data.frequency %>%
      filter(type == "light") %>%
      ungroup() %>%
      top_n(3, freq)
  ) %>%
  pull(label);

# Define a common size scale
bubble.min <- min(data.frequency$freq)  # Smallest size across both datasets
bubble.max <- max(data.frequency$freq)  # Largest size across both datasets
size_range <- c(3, 10)  # Visual size range for bubbles

color.cluster <- c('#FF5C8F', '#6983C9');
p.bubble.light <- ggplot(data.frequency %>% filter(type == 'light'), aes(x = v_call, y = j_call, fill = cluster, label=label)) +
  geom_point(data = subset(data.frequency, type =='light' & !label %in% germline.top.list), aes(size=freq, alpha=0.7), shape = 21) + 
  geom_point(data = subset(data.frequency, type =='light' & label %in% germline.top.list), aes(size=freq, alpha=1.0), shape = 21) + 
  geom_text_repel(
    data = subset(data.frequency, type == 'light' & label %in% germline.top.list),
    seed=42,
    box.padding = 0.5,  # Adjust padding as needed
    max.time=1,
    segment.color = 'gray85',
    segment.size = 0.1,
    nudge_y=1,
    max.iter=90000,
    min.segment.length = 0.01,
    size=2
  ) +
  xlab("V Gene") +
  ylab("J Gene") +
  scale_fill_manual(
    values=color.cluster, 
    labels = c("Cluster 1", "Cluster 2"),
    guide = "none",
  ) +
  scale_size_continuous(
    limits = c(bubble.min, bubble.max),
    range = size_range,
    name = "Frequency"
  ) +
  coord_cartesian(clip = 'off') +
  scale_alpha(guide = 'none') +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border=element_blank(),
    axis.text.y = element_text(size=6),
    axis.line=element_line(color='black'),
    axis.title.x = element_text(size=6),
    axis.title.y = element_text(size=6),
    panel.grid.major.x = element_line(color='#ECECEC', linewidth=0.3),
    axis.ticks.y=element_blank(),
    axis.ticks.x=element_blank(),
    legend.direction="horizontal",
    legend.key.size = unit(4, 'mm'),
    plot.margin = margin(t=0, r=1, b=2, l=0, unit="mm"), 
    strip.background =element_blank(),
    strip.text = element_text(size=6),
    axis.text.x = element_text(angle = 45, hjust = 1, size=6)
  );

p.bubble.heavy <- ggplot(data.frequency %>% filter(type == 'heavy'), aes(x = v_call, y = j_call, fill = cluster, label=label)) +
  geom_point(data = subset(data.frequency, type =='heavy' & !label %in% germline.top.list), aes(size=freq, alpha=0.5), shape = 21) + 
  geom_point(data = subset(data.frequency, type =='heavy' & label %in% germline.top.list), aes(size=freq, alpha=1.0), shape = 21) + 
  geom_text_repel(
    data = subset(data.frequency, type == 'heavy' & label %in% germline.top.list),
    seed=42,
    box.padding = 0.5,  # Adjust padding as needed
    max.time=1,
    segment.color = 'gray85',
    segment.size = 0.1,
    nudge_y=1,
    max.iter=90000,
    min.segment.length = 0.01,
    size=2
  ) +
  xlab("V Gene") +
  ylab("J Gene") +
  scale_fill_manual(
    values=color.cluster, 
    labels = c("Cluster 1", "Cluster 2"),
    # guide = "none",
  ) +
  scale_size_continuous(
    limits = c(bubble.min, bubble.max),
    range = size_range,
    name = "Frequency"
  ) +
  coord_cartesian(clip = 'off') +
  guides(fill=guide_legend(nrow=1, override.aes=list(shape=22, size=10)), alpha='none', size=guide_legend(nrow = 1)) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size=6),
    axis.text.y = element_text(size=6),
    axis.line=element_line(color='black'),
    axis.title.x = element_text(size=6),
    axis.title.y = element_text(size=6),
    panel.grid.major.x = element_line(color='#ECECEC', linewidth=0.3),
    axis.ticks.y=element_blank(),
    axis.ticks.x=element_blank(),
    legend.direction="horizontal",
    legend.key.size = unit(4, 'mm'),
    plot.margin = margin(t=0, r=1, b=2, l=0, unit="mm"), 
    strip.background=element_blank(),
    strip.text = element_text(size=6)
  );

# Get annotation legend
p.bubble.legend  <- get_legend(p.bubble.heavy + theme(
  legend.title = element_blank(),
  legend.direction = "horizontal",
  legend.justification="center",
  legend.box = "horizontal",   # Arrange legends horizontally
  legend.box.just = "bottom", 
  legend.title.align=0.5,
  legend.text=element_text(size=6)));

# Remove the legend before plotting
p.bubble.heavy <- p.bubble.heavy + theme(legend.position='none');
p.bubble.light <- p.bubble.light + theme(legend.position='none');

# Create a violin plot with the somatic hypermutations for VH & VL
violin.min <- min(data.violin$value);
violin.max <- max(data.violin$value);
p.violin.heavy <- ggplot(data.violin %>% filter(type == 'heavy'), aes(x=key, y=value, fill=cluster)) + 
  geom_split_violin(nudge = 0.02) + 
  geom_point(
    shape=16,
    position=position_jitterdodge(jitter.width=0.3,dodge.width=0.6)
  ) + 
  scale_x_discrete(
    limits = c("v_identity", "d_identity", "j_identity"),
    labels = c("v_identity"="V gene", "d_identity" = "D gene", "j_identity"="J gene")
  ) +
  scale_y_continuous(limits = c(violin.min, violin.max)) + 
  scale_fill_manual(
    values=color.cluster, 
    labels = c("Cluster 1", "Cluster 2"),
    guide = "none",
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border=element_blank(),
    axis.text.y = element_text(size=6),
    axis.line=element_line(color='black'),
    axis.text.x = element_text(size=6),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.x = element_line(color='#ECECEC', linewidth=0.3),
    axis.ticks.y=element_blank(),
    axis.ticks.x=element_blank(),
    legend.direction="horizontal",
    legend.key.size = unit(4, 'mm'),
    plot.margin = margin(t=0, r=1, b=2, l=0, unit="mm"), 
    strip.background =element_blank(),
    strip.text = element_text(size=6),
    legend.position='none',
  );

p.violin.light <- ggplot(data.violin %>% filter(type == 'light'), aes(x=key, y=value, fill=cluster)) + 
  geom_split_violin(nudge = 0.02) + 
  geom_point(
    shape=16,
    position=position_jitterdodge(jitter.width=0.3,dodge.width=0.6)
  ) + 
  scale_x_discrete(
    limits = c("v_identity", "j_identity"),
    labels = c("v_identity"="V gene", "j_identity"="J gene")
  ) +
  scale_y_continuous(limits = c(violin.min, violin.max)) +  # Set y-axis min and max
  scale_fill_manual(
    values=color.cluster, 
    labels = c("Cluster 1", "Cluster 2"),
    guide = "none",
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border=element_blank(),
    axis.text.y = element_text(size=6),
    axis.line=element_line(color='black'),
    axis.text.x = element_text(size=6),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.x = element_line(color='#ECECEC', linewidth=0.3),
    axis.ticks.y=element_blank(),
    axis.ticks.x=element_blank(),
    legend.direction="horizontal",
    legend.key.size = unit(4, 'mm'),
    plot.margin = margin(t=0, r=1, b=2, l=0, unit="mm"), 
    strip.background =element_blank(),
    strip.text = element_text(size=6),
    legend.position='none',
  );

# Create a heatmap with dendrogram for VH and VL
p.heatmap.heavy <- ggplot(heatmap.heavy, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile() +
  labs(fill = "Distance") +
  guides(fill = guide_colourbar(title.position = "top")) + 
  scale_fill_viridis(
      alpha = 1,
      begin = 0,
      end = 1,
      direction = -1,
      discrete = FALSE,
      na.value = "grey90",
      limits = c(0.0, 1),
      option = "D",
      n.breaks=7,
    ) +
  # coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=6),
    axis.text.y = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank()
  );

p.heatmap.light <- ggplot(heatmap.light, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile() +
  labs(fill = "Distance") +
  guides(fill = guide_colourbar(title.position = "top")) + 
  scale_fill_viridis(
      alpha = 1,
      begin = 0,
      end = 1,
      direction = -1,
      discrete = FALSE,
      na.value = "grey90",
      limits = c(0.0, 1),
      option = "D",
      n.breaks=7,
    ) +
  # coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=6),
    axis.text.y = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank()
  );

p.tree.heavy <-  ggtree(cluster.heavy, ladderize=FALSE) %<+% data.full.heavy +
    geom_tiplab(size=2, align=TRUE) + 
    geom_tippoint(aes(colour=factor(cluster)), size=1.5) +
    scale_color_manual(values=c('#FF5C8F', '#6983C9'), labels = c("Cluster 1", "Cluster 2")) +
    xlim_tree(0.2) + 
    coord_cartesian(clip = "off") + 
    theme_tree2() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=6)
    );
p.tree.light <-  ggtree(cluster.light, ladderize=FALSE) %<+% data.full.light +
    geom_tiplab(size=2, align=TRUE) + 
    geom_tippoint(aes(colour=factor(cluster)), size=1.5) +
    scale_color_manual(values=c('#FF5C8F', '#6983C9'), labels = c("Cluster 1", "Cluster 2")) +
    xlim_tree(0.5) + 
    coord_cartesian(clip = "off") +
    theme_tree2() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=6)
    );

# Get heatmap legend
p.heatmap.legend <- get_legend(p.heatmap.heavy + theme(
  legend.title=element_text(size=6),
  legend.direction = "horizontal",
  legend.justification="center",
  legend.box.just = "bottom", 
  legend.title.align=0.5,
  legend.text=element_text(size=6)
));

# Get tree legend
p.tree.legend <- get_legend(p.tree.heavy + theme(
  legend.direction = "horizontal",
  legend.justification="center",
  legend.box.just = "bottom", 
  legend.title.align=0.5,
  legend.title=element_blank(),
  legend.text=element_text(size=6)));

# Remove the legend before plotting
p.heatmap.heavy <- p.heatmap.heavy + theme(legend.position='none');
p.heatmap.light <- p.heatmap.light + theme(legend.position='none');
p.tree.heavy <- p.tree.heavy + theme(legend.position='none');
p.tree.light <- p.tree.light + theme(legend.position='none');

# Assemble the final graph
legend.y <- 0.3;
plot <- plot_grid(
  plot_grid(
    plot_grid(
      p.bubble.heavy,
      p.violin.heavy,
      labels=c('A', 'B'),
      align="v",
      ncol=1,
      vjust = -0.2,
      hjust = -0.1
    ),
    plot_grid(
      p.bubble.light,
      p.violin.light,
      align="v",
      ncol=1,
      vjust = -0.2,
      hjust = -0.1
    ),
    ncol=2,
    align='h'
  ),
  p.bubble.legend,
  plot_grid(
    plot_grid(
      p.tree.heavy,
      p.heatmap.heavy,
      nrow=1,
      rel_widths = c(0.7, 2),
      align = "h"
    ),
    plot_grid(
      p.tree.light,
      p.heatmap.light,
      nrow=1,
      rel_widths = c(0.7, 2),
      align = "h"
    ),
    vjust=legend.y,
    labels=c('C', 'D'),
    align = "h",
    ncol=2,
    nrow=1
  ),
  plot_grid(
    p.tree.legend,
    p.heatmap.legend,
    nrow=1,
    align='h'
  ),
  vjust = -0.2,
  hjust = -0.1,
  ncol=1,
  rel_heights = c(
    1.9,
    0.3,
    1,
    0.2
  )
)+ theme(plot.margin = margin(t=6, r=4, b=4, l=4, unit="mm"));

# Save the final plot
ggsave(filename=opt$output, plot=plot, units="mm", width=170, height=210, dpi=300);

# Often the most useful way to look at many warnings:
summary(warnings());