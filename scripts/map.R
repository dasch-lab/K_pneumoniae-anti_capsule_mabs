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

  return (eval(parse(text=paste("require(",package,")"))))
}

usePackage.github <- function(package, repository){

  library.path <- .libPaths()[1];

  if(!is.element(package, installed.packages(lib.loc=library.path)[,1]))
    devtools::install_github(repository);

  return (eval(parse(text=paste("require(",package,")"))))
}

# For command line parsing
usePackage("argparse");

# create parser object
parser <- ArgumentParser()
parser$add_argument("-i", "--input", default=NULL, type="character", help="input blast file");
parser$add_argument("-o", "--output", default=NULL, type="character", help="output pdf file");

opt <- parser$parse_args();

# Install packages
usePackage("ggplot2");
usePackage("tidyr");
usePackage('gridExtra');
usePackage('ggnewscale');
usePackage("dplyr");
usePackage('plyr');
usePackage('stringr');
usePackage("cowplot");  # For plot_grid
usePackage("BiocManager");
usePackage("devtools");
usePackage("ggrepel");
usePackage("devtools");
usePackage("Biostrings");

usePackage.bio('rtracklayer');
usePackage.bio('msa');

usePackage.github("thacklr", "thackl/thacklr");
usePackage.github("gggenomes", "thackl/gggenomes");
options(max.print = 1000, width = 1000) # Adjust 'width' as needed

# Define loci genes
o_locus = c('wzm','wzt','wbbM','glf','wbbN','wbbO','kfoC','gmlC','gmlB','gmlA','wbbY');
k_locus = c('galF','cpsACP','wzi','wza','wzb','wzc','wzx','wcoV','wzy','wcoU','wcoT','wcsF','wcuK','wbaZ','wcaJ','gnd','manC','manB','rmlB','rmlA','rmlD','rmlC','ugd','rfaG','gmd','wckX','wckW','wcpO','wbaP','wcaN');

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

extractToken <- function(data_list, token_id){
  result <- c();
  for(data in data_list){
    tokens <- unlist(strsplit(data, split = "|", fixed = TRUE))
    if(length(tokens) > token_id){
      label <- paste(tokens[token_id]);
    }else{
      label <- '<none>';
    }
    
    result <- c(result, label)
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

check_protein <- function(sequenceList, strandList){
  result <- c();
  for (i in 1:length(sequenceList)) {
    sequence <- sequenceList[[i]];
    strand <- strandList[[i]];
    sequence <- gsub("-", "", sequence);
    if(strand != '+'){
      sequence <- reverse_complement(sequence);
    }
    
    # Translate to protein
    protein <- translate_dna(sequence);
    
    # Check for stop codon presence
    stop_list <- which(strsplit(protein, "")[[1]] == "*");
    has_stop <- FALSE;
    if(length(stop_list) > 0){
      first_stop <- stop_list[[1]];
      if(first_stop < nchar(protein)){
        has_stop <- TRUE;
      }
    }
    
    # Append results
    result <- c(result, has_stop)
  }

  return (result);
}

complement_data <- function(sequenceList, strandList){
  result <- c();
  for (i in 1:length(sequenceList)) {
    sequence <- sequenceList[[i]];
    strand <- strandList[[i]];
    sequence <- gsub("-", "", sequence);
    if(strand != '+'){
      sequence <- reverse_complement(sequence);
    }

     # Append results
    result <- c(result, sequence)
  }

  return (result);
}

translate_protein <- function(sequenceList, orf){
  result <- c();
  for (i in 1:length(sequenceList)) {
    result <- c(result, translate_dna(sequenceList[[i]], orf));
  }

  return (result);
}

reverse_complement <- function(sequence) {
  # Convert the sequence to uppercase for consistency
  sequence <- toupper(sequence)
  
  # Reverse the sequence
  reversed_sequence <- rev(strsplit(sequence, NULL)[[1]]);
  
  # Create a mapping for nucleotide complements
  complement_mapping <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C", "N"="N", "-"="-")
  
  # Replace each nucleotide with its complement
  complemented_sequence <- sapply(reversed_sequence, function(nucleotide) {
    complement_mapping[nucleotide]
  });
  
  # Concatenate the complemented nucleotides to form the reverse complement
  reverse_complement_sequence <- paste0(complemented_sequence, collapse = "")

  return(reverse_complement_sequence)
}

translate_dna <- function(sequence, orf=NA) {
  # Convert the sequence to uppercase for consistency
  sequence <- toupper(sequence)
  
  # Define a mapping of codons to amino acids
  codon_table <- list(
    "TTT" = "F", "TTC" = "F", "TTA" = "L", "TTG" = "L",
    "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
    "ATT" = "I", "ATC" = "I", "ATA" = "I", "ATG" = "M",
    "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
    "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
    "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
    "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
    "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
    "TAT" = "Y", "TAC" = "Y", "TAA" = "*", "TAG" = "*",
    "CAT" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
    "AAT" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
    "GAT" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
    "TGT" = "C", "TGC" = "C", "TGA" = "*", "TGG" = "W",
    "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
    "AGT" = "S", "AGC" = "S", "AGA" = "R", "AGG" = "R",
    "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"
  )
  
  # Function to translate a single codon
  translate_codon <- function(codon) {
    amino_acid <- codon_table[[codon]]
    return(amino_acid)
  }
  
  # Function to translate a DNA sequence
  translate_sequence <- function(sequence) {
    codons <- strsplit(sequence, "");
    amino_acids <- sapply(codons, function(codon) {
      codon <- paste(codon, collapse = "")
      translate_codon(codon)
    })
    return(amino_acids)
  }
  
  # Function to find the reading frame with the fewest stop codons
  find_best_reading_frame <- function(frames) {
    stop_counts <- sapply(frames, function(frame) {
      sum(frame == "*")
    })
    best_frame <- which.min(stop_counts)
    return(frames[[best_frame]])
  }
  
  # Generate three reading frames
  frames <- lapply(0:2, function(offset) {
    # frame <- substring(sequence, start = offset + 1)
    frame <- substring(sequence, seq(offset + 1, nchar(sequence), 3), seq(offset + 3, nchar(sequence), 3))
    
    amino_acids <- translate_sequence(frame)
    return(amino_acids)
  })
  
  if (is.na(orf)){
    # Find the best reading frame
    best_frame <- find_best_reading_frame(frames)
  }else{
    frame <- substring(sequence, seq(0 + 1, nchar(sequence), 3), seq(0 + 3, nchar(sequence), 3));
    best_frame <- translate_sequence(frame);
  }
  
  
  # Concatenate amino acids into a protein sequence
  protein_sequence <- paste(best_frame, collapse = "")
  
  return(protein_sequence)
}

# Load dataset
gene_gap <- 70; # Gap between genes
data.blast <- read.csv(opt$input, header = TRUE, sep="\t", row.names=NULL) %>%
  dplyr::mutate(gene=extractToken(sseqid, 2)) %>%
  dplyr::mutate(coverage=100*((qend-qstart)/slen)) %>%
  dplyr::filter(gene != "") %>%
  dplyr::filter(coverage >= 90) %>%
  dplyr::group_by(assembly, gene, locus) %>%
  dplyr::slice_max(pident, with_ties = FALSE, n=1) %>%
  dplyr::ungroup() %>%
  dplyr::filter((gene %in% k_locus & locus == 'k_locus') | (gene %in% o_locus & locus == 'o_locus')) %>%
  dplyr::arrange(assembly, qstart) %>%
  dplyr::mutate(len = qend - qstart) %>%
  dplyr::mutate(start = cumsum(len) + (row_number() - 1) * gene_gap - len) %>% 
  dplyr::mutate(end=start+len) %>%
  dplyr::mutate(strand=ifelse(sstart > send, '-', '+')) %>%
  dplyr::mutate(type = 'CDS') %>%
  dplyr::mutate(has_stop = check_protein(qseq, strand)) %>%
  dplyr::rename(seq_id = assembly, identity=pident) %>%
  dplyr::mutate(assembly_order = factor(seq_id, levels=names(mapping))) %>%
  dplyr::mutate(seq_id=assembly_order) %>%
  dplyr::select(seq_id, assembly_order, identity, locus, gene, coverage, len, start, end, strand, type, has_stop, qseq) %>%
  dplyr::arrange(assembly_order) %>%
  dplyr::mutate(bin_id=assembly_order);

# Derive genome table
data.genome <- data.blast %>%
  dplyr::group_by(seq_id) %>%
  dplyr::summarise(start=min(start), end=max(end)) %>%
  dplyr::mutate(length=end-start) %>%
  dplyr::mutate(assembly_order = factor(seq_id, levels=names(mapping))) %>%
  dplyr::arrange(assembly_order) %>%
  dplyr::mutate(assembly_label = rename_assembly(assembly_order)) %>%
  dplyr::arrange(assembly_order) %>%
  dplyr::mutate(seq_id=assembly_order) %>%
  dplyr::mutate(bin_id=assembly_order);

# Select wbbO
data.alignment <- data.blast %>%
  dplyr::filter(gene == "wbbO") %>%
  dplyr::mutate(sequence=complement_data(qseq, strand)) %>%
  dplyr::select(seq_id, gene, sequence, strand);
  
# Align sequences
alignment.size <- 26;
alignment.start <- 460;
alignment.end <- alignment.start + alignment.size;
alignment.wbbO <- msa(data.alignment$sequence, type="dna", order = "input", method = "Muscle")

# Extract aligned sequences as a string
data.alignment <- data.alignment %>%
  dplyr::mutate(align=as.character(alignment.wbbO)) %>%
  dplyr::mutate(align=substr(align, alignment.start, alignment.end)) %>%
  dplyr::select(seq_id, align) %>%
  separate(align, paste0("pos_", seq(alignment.start, alignment.end)), sep=seq(1, alignment.size)) %>%
  pivot_longer(cols = -seq_id, names_to = "key", values_to = "value") %>%
  dplyr::mutate(position = as.integer(str_replace(key, "pos_", ""))) %>%
  dplyr::select(seq_id, value, position) %>%
  dplyr::mutate(assembly_label = rename_assembly(seq_id)) %>%
  dplyr::mutate(fill=ifelse(value == '-', 'error', 'empty'));

# Sort by labels
data.alignment$seq_id <- factor(data.alignment$seq_id, levels=rev(names(mapping)));

p.sequence.label <- ggplot(subset(data.alignment, position %in% c(462)), aes(x=0, y=seq_id)) +
  geom_text(
    aes(label=assembly_label), 
    parse = TRUE, 
    size=3
  )+
  theme_void() + 
  theme(
    plot.margin = margin(t=0, r=0, b=0, l=0, unit="mm")
  )
  
p.sequence <- ggplot(data.alignment, aes(x=position, y=seq_id)) +
  geom_tile(aes(fill=fill), color='white', size=1) +
  scale_discrete_manual(
    aesthetics="fill",
    values=c("error" = "#DE5458", 'empty'='#F4F4F4')
  ) +
  scale_discrete_manual(
    aesthetics="colour",
    values=c("-" = "#FFFFFF", 'A'='#46A33C', 'T'='#EF5959', 'G'='#000000', 'C'='#258BD0')
  ) +
  scale_x_continuous(breaks = round(seq(min(data.alignment$position), max(data.alignment$position), by = 3),1)) +
  geom_vline(xintercept = as.numeric(seq(alignment.start-0.5, alignment.start+alignment.size-0.5, by = 3)), color = "#C2C2C2", linetype = "dashed", size = 0.5) + 
  geom_text(aes(label=value, colour=value), size=2) +
  guides(colour='none', alpha='none', fill='none') +
  theme_minimal() + 
  theme(
    panel.background=element_rect(fill="white", colour="white"),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=6),
    legend.key.size = unit(0.7, "lines"),
    panel.grid.major.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom", 
    legend.box = "horizontal"
  );

# Select colors
locus.color <- c('k_locus'='#258BD0', 'o_locus'='#EF5959');

# Plot map representation
p.map <- gggenomes(seqs=data.genome, genes=data.blast) +
  geom_gene(aes(fill=locus, alpha=ifelse(has_stop, 0.5, 1.0)),position="pile", shape=c(3,4)) +
  geom_gene_tag(aes(label=gene), nudge_y=0.1, check_overlap = FALSE) +
  geom_bin_label(aes(label=assembly_label),data=bins(.group=vars(assembly_label)), parse = TRUE) +
  scale_fill_manual(
    values=locus.color, 
    labels = c("K locus", "O locus")
  ) +
  guides(colour='none', alpha='none', fill=guide_legend(nrow = 1, title.position = "left")) +
  theme(
    legend.position = "bottom",
    legend.title=element_blank(),
    legend.key.size = unit(4, 'mm'),
    axis.line.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank()
  )

# Flip annotations
p.map <- p.map %>% flip(c(1,3));

# Get legends
p.map.legend <- get_legend(p.map +
  theme(
    legend.title=element_blank(), 
    legend.direction = "horizontal",
    legend.justification="center",
    legend.box.just = "bottom", 
    legend.text=element_text(size=8),
    plot.margin = margin(t=0, r=, b=0, l=0, unit="mm")
  ));

# Remove the legend before plotting
p.map <- p.map + theme(legend.position='none');

plot.sequence <- plot_grid(
    p.sequence.label,
    p.sequence,
    nrow=1,
    ncol=2,
    rel_widths=c(0.2, 1),
    align='h'
  );

plot.map <- plot_grid(
  p.map,
  p.map.legend,
  nrow=2,
  ncol=1,
  rel_heights=c(5, 0.2)
);

# Generate the grid
plot <- plot_grid(
  plot.sequence,
  plot.map,
  nrow=2,
  vjust=0.3,
  labels=c('A', 'B'),
  rel_heights=c(2, 5)
) + theme(plot.margin = margin(t=4, r=4, b=0, l=4, unit="mm"));

# Save the final plot
ggsave(filename=opt$output, plot=plot, units="mm", width=170, height=180, dpi=600);

# Get output path and directory
output_directory <- dirname(opt$output);
output_name <- tools::file_path_sans_ext(basename(opt$output));
output_extension <- tools::file_ext(opt$output);

# Save the two panels separately
ggsave(filename=file.path(output_directory, paste0(output_name, '_panelA.',output_extension)), plot=plot.sequence, units="mm", width=170, height=60, dpi=600);
ggsave(filename=file.path(output_directory, paste0(output_name, '_panelB.',output_extension)), plot=plot.map, units="mm", width=170, height=100, dpi=600);

# Often the most useful way to look at many warnings:
summary(warnings());

# Print package versions
cat("R version: ", R.version.string, "\n")
cat("ggplot2 version: ", unlist(packageVersion("ggplot2")), "\n")
cat("dplyr version: ", unlist(packageVersion("dplyr")), "\n")