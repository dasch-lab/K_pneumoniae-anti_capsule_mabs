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

usePackage.bio('ggtree');
usePackage.bio('treeio');

color.IGH <- '#CF3B3B';
color.IGK <- '#45D176';
color.IGL <- '#25934B';
color.border <- 'white';

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

data.input <- read.csv(opt$input, header=TRUE, sep="\t", row.names=NULL) %>%
  dplyr::mutate(name=parse_name(sequence_id)) %>%
  dplyr::mutate(v_call=parse_germline(v_call), j_call=parse_germline(j_call)) %>%
  dplyr::mutate(identity=calculate_identity(sequence_alignment_aa, germline_alignment_aa)) %>%
  dplyr::mutate(cdr3_len=nchar(cdr3_aa)) %>%
  dplyr::mutate(cluster=set_cluster(name)) %>%
  dplyr::mutate(type=ifelse(locus == 'IGH', 'heavy', 'light')) %>%
  dplyr::filter(cluster > 0) %>%
  dplyr::select(name, cluster, cdr3_aa, cdr3_len, identity, v_call, j_call, type, sequence_alignment_aa, locus);

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

write_fasta(data.full.heavy, file.path(output_dir, 'antibodies_full_heavy.fasta'));
write_fasta(data.full.light, file.path(output_dir, 'antibodies_full_light.fasta'));
write_fasta(data.CDR3.heavy, file.path(output_dir, 'antibodies_CDR3H.fasta'));
write_fasta(data.CDR3.light, file.path(output_dir, 'antibodies_CDR3L.fasta'));

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

# Generate frequency plot
data.heavy.v <- data.input %>%
  dplyr::filter(type == 'heavy') %>%
  dplyr::group_by(v_call, locus) %>%
  dplyr::summarize(count=n()) %>%
  dplyr::arrange(desc(count));

data.heavy.j <- data.input %>%
  dplyr::filter(type == 'heavy') %>%
  dplyr::group_by(j_call, locus) %>%
  dplyr::summarize(count=n()) %>%
  dplyr::arrange(desc(count));

data.light.v <- data.input %>%
  dplyr::filter(type == 'light') %>%
  dplyr::group_by(v_call, locus) %>%
  dplyr::summarize(count=n()) %>%
  dplyr::arrange(desc(count));

data.light.j <- data.input %>%
  dplyr::filter(type == 'light') %>%
  dplyr::group_by(j_call, locus) %>%
  dplyr::summarize(count=n()) %>%
  dplyr::arrange(desc(count));

# Sort by count
data.heavy.v$v_call <- factor(data.heavy.v$v_call,  levels=data.heavy.v$v_call[order(data.heavy.v$count, decreasing=FALSE)]);
data.light.v$v_call <- factor(data.light.v$v_call,  levels=data.light.v$v_call[order(data.light.v$count, decreasing=FALSE)]);
data.heavy.j$j_call <- factor(data.heavy.j$j_call,  levels=data.heavy.j$j_call[order(data.heavy.j$count, decreasing=FALSE)]);
data.light.j$j_call <- factor(data.light.j$j_call,  levels=data.light.j$j_call[order(data.light.j$count, decreasing=FALSE)]);

# Generate a frequency map for heavy and light chains
data.heavy.frequency <- data.input %>%
  dplyr::filter(type == 'heavy') %>%
  dplyr::mutate(total=n()) %>%
  dplyr::group_by(v_call, j_call, total) %>%
  dplyr::summarize(count=n()) %>%
  dplyr::mutate(freq=100*(count/total));

# Append missing elements to fill the heavy table
v_call_list <- sort(unique(data.heavy.frequency %>% distinct(v_call) %>% pull()));
j_call_list <- sort(unique(data.heavy.frequency %>% distinct(j_call) %>% pull()));
data.heavy.frequency.empty <- expand.grid(count=0, v_call=v_call_list, j_call=j_call_list);
data.heavy.frequency.missing <- anti_join(data.heavy.frequency.empty, data.heavy.frequency, by=c("v_call", 'j_call'));
if(nrow(data.heavy.frequency.missing ) > 0){
  data.heavy.frequency.missing$count <- NA;
  data.heavy.frequency <- full_join(data.heavy.frequency, data.heavy.frequency.missing, by=c("v_call", 'j_call', 'count'));
}

data.light.frequency <- data.input %>%
  dplyr::filter(type == 'light') %>%
  dplyr::mutate(total=n()) %>%
  dplyr::group_by(v_call, j_call, total) %>%
  dplyr::summarize(count=n()) %>%
  dplyr::mutate(freq=100*(count/total));

# Append missing elements to fill the light table
v_call_list <- sort(unique(data.light.frequency %>% distinct(v_call) %>% pull()));
j_call_list <- sort(unique(data.light.frequency %>% distinct(j_call) %>% pull()));
data.light.frequency.empty <- expand.grid(count=0, v_call=v_call_list, j_call=j_call_list);
data.light.frequency.missing <- anti_join(data.light.frequency.empty, data.light.frequency, by=c("v_call", 'j_call'));
if(nrow(data.light.frequency.missing ) > 0){
  data.light.frequency.missing$count <- NA;
  data.light.frequency <- full_join(data.light.frequency, data.light.frequency.missing, by=c("v_call", 'j_call', 'count'));
}

# Start generating plots
print("Rendering");
frequency.text.size <- 6;
frequency.title.size <- 8;
p.heavy.j <- ggplot(data.heavy.j, aes(x = j_call, y = count)) +
  geom_bar(stat = "identity", fill=color.IGH) +
  coord_flip() +  # Flip coordinates for vertical bars
  labs(x = "IGHJ", y = "Count") +
  theme_bw() + 
  theme(
    legend.position='none',
    axis.title.x=element_text(size=frequency.title.size),
    axis.title.y=element_text(size=frequency.title.size),
    axis.text.x=element_text(size=frequency.text.size),
    axis.text.y=element_text(size=frequency.text.size)
  );

p.heavy.v <- ggplot(data.heavy.v, aes(x = v_call, y = count)) +
  geom_bar(stat = "identity", fill=color.IGH) +
  coord_flip() +  # Flip coordinates for vertical bars
  labs(x = "IGHV", y = "Count") +
  theme_bw() + 
  theme(
    legend.position='none',
    axis.title.x=element_text(size=frequency.title.size),
    axis.title.y=element_text(size=frequency.title.size),
    axis.text.x=element_text(size=frequency.text.size),
    axis.text.y=element_text(size=frequency.text.size)
  );

p.light.j <- ggplot(data.light.j, aes(x = j_call, y = count, fill=locus)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip coordinates for vertical bars
  scale_fill_manual(values=c(color.IGK, color.IGL), limits = c('IGK', 'IGL'), labels = c("IGK", "IGL")) +
  labs(x = "IGKJ/IGLJ", y = "Count") +
  theme_bw() +
  theme(
    legend.position='none',
    axis.title.x=element_text(size=frequency.title.size),
    axis.title.y=element_text(size=frequency.title.size),
    axis.text.x=element_text(size=frequency.text.size),
    axis.text.y=element_text(size=frequency.text.size)
  );

p.light.v <- ggplot(data.light.v, aes(x = v_call, y = count, fill=locus)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip coordinates for vertical bars
  scale_fill_manual(values=c(color.IGK, color.IGL), limits = c('IGK', 'IGL'), labels = c("IGK", "IGL")) +
  labs(x = "IGKV/IGLV", y = "Count") +
  theme_bw() + 
  theme(
    legend.position='none',
    axis.title.x=element_text(size=frequency.title.size),
    axis.title.y=element_text(size=frequency.title.size),
    axis.text.x=element_text(size=frequency.text.size),
    axis.text.y=element_text(size=frequency.text.size)
  );

# Plot the heatmap
p.heavy.frequency <- ggplot(data.heavy.frequency) + 
  geom_tile(mapping=aes(x=v_call, y=j_call, fill=count), colour=color.border, size=1) +
  labs(fill = "Frequency (%)", x = "IGHV", y = "IGHJ") + 
  guides(fill = guide_colourbar(title.position = "top")) + 
    scale_fill_viridis(
      alpha = 1,
      begin = 0,
      end = 1,
      direction = -1,
      discrete = FALSE,
      na.value = "grey90",
      option = "A",
      trans = "log10",
      limits = c(1, 100),
      n.breaks=7,
    ) +
    theme_bw() + 
    theme(
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_text(size=frequency.text.size),
      axis.text.x=element_text(
        angle = 45, 
        vjust = 1, 
        size = frequency.text.size, 
        hjust = 1, 
        margin = margin(t = 2, r = 0, b = 0, l = 0)
      )
    );

# Plot the heatmap
p.light.frequency <- ggplot(data.light.frequency) + 
  geom_tile(mapping=aes(x=v_call, y=j_call, fill=count), colour=color.border, size=1) +
  labs(fill = "Frequency (%)", x = "IGKV/IGLV", y = "IGKJ/IGLJ") +
  guides(fill = guide_colourbar(title.position = "top")) + 
    scale_fill_viridis(
      alpha = 1,
      begin = 0,
      end = 1,
      direction = -1,
      discrete = FALSE,
      na.value = "grey90",
      option = "A",
      trans = "log10",
      limits = c(1, 100),
      n.breaks=7,
    ) +
    theme_bw() + 
    theme(
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_text(size=frequency.text.size),
      axis.text.x=element_text(
        angle = 45, 
        vjust = 1, 
        size = frequency.text.size, 
        hjust = 1, 
        margin = margin(t = 2, r = 0, b = 0, l = 0)
      )
    );

# Get heatmap legend
p.frequency.legend <- get_legend(p.heavy.frequency + theme(
  legend.title=element_text(size=6),
  legend.direction = "horizontal",
  legend.justification="center",
  legend.box.just = "bottom", 
  legend.title.align=0.5,
  legend.text=element_text(size=6)
  ));

# Remove the legend before plotting
p.heavy.frequency <- p.heavy.frequency + theme(legend.position='none');
p.light.frequency <- p.light.frequency + theme(legend.position='none');

# Identity plot
p.identity <- ggplot(data.input, aes(x=type, y=identity, fill=type)) + 
  geom_violin() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  scale_fill_manual(values=c(color.IGH, color.IGK), limits = c('heavy', 'light'), labels = c("Heavy", "Light")) +
  theme_bw() +
  theme(
    legend.position='none'
  );

data.cdr3 <- data.input %>%
  dplyr::group_by(type, cdr3_len) %>%
  dplyr::summarize(count=n());

p.cdr3_len <- ggplot(data.cdr3) + 
  geom_density(aes(x=cdr3_len, group=type, fill=type), alpha=0.5, adjust=2,kernel = "rectangular") + 
  facet_grid(~type) +
  scale_fill_manual(values=c(color.IGH, color.IGK), limits = c('heavy', 'light'), labels = c("Heavy", "Light")) +
  theme_bw() + 
  theme(
    legend.position='none'
  );

# Create a dummy text for generating the legend
data.ab.legend <- data.frame(
  x = c(1, 2, 3),
  y = c(3, 5, 4),
  type = c('IGH', 'IGL', 'IGK')
);

# Plotting the legend
p.ab.legend <- ggplot(data.ab.legend, aes(x=x, y=y, fill=type)) +   
  geom_tile() + 
  scale_fill_manual(values = c(color.IGH, color.IGK, color.IGL), labels=c('IGH', 'IGK', 'IGL')) + 
  theme_void() + 
  theme(
      legend.title=element_blank(),
      # legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification="center",
      legend.box.just = "bottom", 
      legend.title.align=0.5,
      legend.text=element_text(size=6),
      legend.key.height= unit(0.5, 'cm'),
      legend.key.width= unit(0.5, 'cm')
    );
p.ab.legend <- get_legend(p.ab.legend);

# Create a heatmap with dendrogram
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
      # trans = "log10",
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
      # trans = "log10",
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
  # legend.title=element_blank(),
  legend.text=element_text(size=6)
));

# Get tree legend
p.tree.legend <- get_legend(p.tree.heavy + theme(
  # legend.title=element_text(size=6),
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

# 
legend.y <- 0.3;
plot <- plot_grid(
  plot_grid(
    plot_grid(
      p.heavy.v,
      p.heavy.j,
      p.heavy.frequency,
      vjust=legend.y,
      labels=c('A', NA, 'B'),
      ncol=1,
      align="v",
    rel_heights=c(0.8,0.8,1.0)
    ),
    plot_grid(
      p.light.v,
      p.light.j,
      p.light.frequency,
      vjust=legend.y,
      # labels=c('B', NULL, 'F'),
      ncol=1,
      align="v",
    rel_heights=c(0.8,0.8,1.0)
    ),
    ncol=2,
    align='h'
  ),
  plot_grid(
    p.ab.legend,
    p.frequency.legend,
    nrow=1
  ),
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
  ncol=1,
  rel_heights = c(1.9,0.2, 1,0.2)
)+ theme(plot.margin = margin(t=4, r=4, b=0, l=4, unit="mm"));

# Save the final plot
ggsave(filename=opt$output, plot=plot, units="mm", width=170, height=210, dpi=300);

# Often the most useful way to look at many warnings:
summary(warnings());