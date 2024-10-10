packages = c("ggalluvial","ggplot2","EnrichedHeatmap","rtracklayer","dplyr","zeallot")
invisible(
  suppressMessages(
    if(!require("BiocManager",character.only = TRUE,quietly = TRUE)){
      cat("Installing BiocManager\n",sep = "")
      install.packages("BiocManager")
    }))

cat("#####   Loading packages   #####\n")
invisible(
  suppressMessages(
    lapply(packages,function(x){
      if(!require(x,character.only = TRUE,quietly = TRUE)){
        cat("Installing package: ",x,"\n",sep = "")
        BiocManager::install(x,update = FALSE,ask = FALSE)
        library(x,character.only = TRUE,quietly = TRUE)
      }
    })))

options(stringsAsFactors = FALSE)

init_params = list() 

if(!interactive()){
  options(rgl.useNULL = TRUE)
  args = R.utils::commandArgs(trailingOnly = TRUE,asValues = TRUE)
  must_args = c("wd","chains","blocks","species","lengths")
  if(!all(must_args %in% names(args))){
    help = matrix(data = c("--wd    ","Working directory path",
                           "--chains    ","Chain alignment files (comma-separated)",
                           "--blocks    ","tBLASTn Arabidopsis block alignment file",
                           "--species    ","Name of the species (Query species first; comma-separated)",
                           "--lengths    ","Species chromosome lengths files (comma-separated)",
                           "--chr_cutoff    ","Chromosome length cutoff (default 1Mbp = 1e6)",
                           "--chain_cutoff    ","Genome chain alignment cutoff (default 1Kbp = 1e4)",
                           "--cont_chain_cutoff    ","Contiguous chain alignment cutoff (default 1Kbp = 1e4)",
                           "--cont_block_cutoff    ","Contiguous block alignment cutoff (default 1Mbp = 1e6)",
                           "--circos_cutoff    ","Circos window cutoff (default 1Mbp = 1e6)"
    )
    ,ncol = 2,byrow = TRUE)
    prmatrix(help,quote = FALSE,rowlab = rep("",nrow(help)),collab = rep("",2))
    stop(paste0("Missing command line input --> ",paste(must_args[!must_args %in% names(args)],collapse = " | ")), call. = TRUE)
  }
  init_params[["wd"]] = normalizePath(args[["wd"]])
  setwd(init_params[["wd"]])
  
init_params[["lengths"]] = normalizePath(strsplit(args[["lengths"]],split = ",")[[1]])
init_params[["chains"]] = normalizePath(strsplit(args[["chains"]],split = ",")[[1]])
init_params[["blocks"]] = normalizePath(args[["blocks"]])
init_params[["species"]] = gsub(pattern = " ",replacement = "_",strsplit(args[["species"]],split = ",")[[1]])

if("chr_cutoff" %in% names(args)){
  init_params[["chr_cutoff"]] = as.numeric(args[["chr_cutoff"]])
}else{
  init_params[["chr_cutoff"]] = 1e6
}

if("chain_cutoff" %in% names(args)){
  init_params[["chain_cutoff"]] = as.numeric(args[["chain_cutoff"]])
}else{
  init_params[["chain_cutoff"]] = 1e4
}

if("cont_chain_cutoff" %in% names(args)){
  init_params[["cont_chain_cutoff"]] = as.numeric(args[["cont_chain_cutoff"]])
}else{
  init_params[["cont_chain_cutoff"]] = 1e4
}

if("cont_block_cutoff" %in% names(args)){
  init_params[["cont_block_cutoff"]] = as.numeric(args[["cont_block_cutoff"]])
}else{
  init_params[["cont_block_cutoff"]] = 1e6
}
}else{
  suppressMessages(invisible(library(rstudioapi)))
  init_params[["wd"]] = normalizePath(rstudioapi::selectDirectory(caption = "Choose working directory:"))
  setwd(init_params[["wd"]])
  init_params[["species"]] = sapply(1:as.numeric(readline(prompt = "Number of species: ")),function(x) gsub(pattern = " ",replacement = "_",readline(prompt = "Select species 1 name: ")))
  init_params[["chains"]] = sapply(2:length(init_params$species),function(x) normalizePath(selectFile(caption = paste0("Select chain file for ",init_params$species[x]),path = init_params[["wd"]])))
  init_params[["blocks"]] = normalizePath(selectFile(caption = "Select tBLASTn block alignment file",path = init_params[["wd"]]))
  init_params[["lengths"]] = sapply(1:length(init_params$species),function(x) normalizePath(selectFile(caption = paste0("Select chromosome lengths file for ",init_params$species[x]),path = init_params[["wd"]])))
  
  cutoff_defaults = c(chr_cutoff = 1e6,
                      chain_cutoff = 1e4,
                      cont_chain_cutoff = 1e4,
                      cont_block_cutoff = 1e6,
                      circos_cutoff = 1e6)
  init_params[["chr_cutoff"]] = as.numeric(readline(prompt = "Chromosome length cutoff (default 1e6): "))
  init_params[["chain_cutoff"]] = as.numeric(readline(prompt = "Chain alignment cutoff (default 1e4): "))
  init_params[["cont_chain_cutoff"]] = as.numeric(readline(prompt = "Contiguous chain alignment cutoff (default 1e4): "))
  init_params[["cont_block_cutoff"]] = as.numeric(readline(prompt = "Contiguous block alignment cutoff (default 1e6): "))
  init_params[["circos_cutoff"]] = as.numeric(readline(prompt = "Circos window cutoff (default 1e6): "))
  
  for(i in grep(pattern = "*cutoff",x = names(init_params))){
    if(is.na(init_params[[i]])){
      init_params[[i]] = as.numeric(cutoff_defaults[names(init_params)[i]])
    }
  }
}

cat("#####   Getting lengths   #####\n")
chr_lengths = lapply(init_params$lengths,function(x){
  res = read.table(file = x,header = FALSE,sep = "\t",col.names = c("chr","length"))
  res = res[res$length >= init_params[["chr_cutoff"]],,drop = FALSE]
})
names(chr_lengths) = init_params$species

## Links
combine_chains = function(init_params,chr_lengths){
for(comp_species in 2:length(chr_lengths)){
cat(paste0("#####   Getting ",gsub(pattern = "_",replacement = " ",init_params$species[comp_species])," links   #####\n"))
chain = import.chain(init_params[["chains"]][comp_species-1])
genome_df = data.frame(chr = chr_lengths[[comp_species]]$chr,
                       start = rep(0,length(chr_lengths[[comp_species]]$chr)),
                       end = chr_lengths[[comp_species]]$length)

genome_gr = GRanges(seqnames = genome_df[, 1], ranges = IRanges(genome_df[, 2]+1, genome_df[, 3]))
genome_bins = makeWindows(genome_gr, w = init_params$chain_cutoff)
bins = liftOver(genome_bins, chain)
species2_df = as.data.frame(genome_bins)[, 1:3]
colnames(species2_df) = c("chr","start","end")
sq = seqnames(bins)
pa = PartitioningByEnd(sq)
species1_df = data.frame(
  chr_to = as.vector(unlist(sq))[start(pa)],
  start_to = sapply(start(bins), function(x) x[1]),
  end_to = sapply(end(bins), function(x) x[1])
)
species1_chr = sort(unique(species1_df$chr_to))
species2_chr = sort(unique(species2_df$chr))
l = species2_df[, 1] %in% species2_chr & species1_df[, 1] %in% species1_chr & !is.na(species1_df$start_to) & !is.na(species2_df$start)
species2_df = species2_df[l, ]
species1_df = species1_df[l, ]

species1_chromInfo = data.frame(chr = species1_chr,
                                start = rep(0,times = length(species1_chr)),
                                end = chr_lengths[[comp_species]]$length[match(species1_chr,chr_lengths[[comp_species]]$chr)])
species2_chromInfo = data.frame(chr = species2_chr,
                                start = rep(0,times = length(species2_chr)),
                                end = chr_lengths[[comp_species]]$length[match(species2_chr,chr_lengths[[comp_species]]$chr)])


norm_factor = min(sum(species1_chromInfo$end), sum(species2_chromInfo$end))

species1_chromInfo$normalized_length = species1_chromInfo$end / sum(species1_chromInfo$end) * norm_factor
species2_chromInfo$normalized_length = species2_chromInfo$end / sum(species2_chromInfo$end) * norm_factor
chromInfo = rbind(species2_chromInfo,species1_chromInfo)

if(exists("df_combined")){
  df_combined = rbind(df_combined,
        cbind(species1_df,
              species2_df,
              data.frame(species = rep(init_params[["species"]][comp_species],times = nrow(species1_df))))
  )
}else{
  df_combined = cbind(species1_df,
                      species2_df,
                      data.frame(species = rep(init_params[["species"]][comp_species],times = nrow(species1_df))))
}
}

cat("#####   Combining multi-species data   #####\n")
df_combined = df_combined %>%
  group_by(chr_to) %>%
  arrange(start_to) %>%
  mutate(group = cumsum(c(1, diff(start_to) > init_params[["cont_chain_cutoff"]])))

return(df_combined)
}

df_combined = combine_chains(init_params = init_params,chr_lengths = chr_lengths)
alluvial_df = do.call(rbind,lapply(chr_lengths[[1]]$chr,function(chr){
  df_temp = df_combined[df_combined[["chr_to"]] %in% chr,]
  res = do.call(rbind,lapply(1:max(df_temp[["group"]]),function(grp){
    cat("Processing ",chr," | ",format(round(100*grp/max(df_temp[["group"]]),digits = 1),nsmall = 1),"%\r",sep = "")
    for(sp_idx in 2:length(init_params$species)){
      sp = init_params$species[sp_idx]
      if(sp_idx == 2){
        ref_df = as.data.frame(df_temp[df_temp[["group"]] %in% grp & df_temp[["species"]] %in% names(which.max(table(df_temp[df_temp[["group"]] %in% grp,"species"]))),
                                       c("chr_to","start_to","end_to"),
                                       drop = FALSE])
        colnames(ref_df) = c("chr1","start1","end1")
        sp_df = as.data.frame(df_temp[df_temp[["group"]] %in% grp & df_temp[["species"]] %in% sp,c("chr","start","end"),drop = FALSE])
        colnames(sp_df) = c("chr2","start2","end2")
        max_sp = max(table(df_temp[df_temp[["group"]] %in% grp,"species"])) - nrow(sp_df)
        if(max_sp){
          res_df = cbind(ref_df,
                         rbind(sp_df,
                               as.data.frame(matrix(c("No match",NA,NA),
                                                    ncol = 3,
                                                    nrow = max_sp,
                                                    byrow = TRUE,
                                                    dimnames = list(NULL,colnames(sp_df))))
                         ))
        }else{
          res_df = cbind(ref_df,sp_df)
        }
      }else{
        sp_df = as.data.frame(df_temp[df_temp[["group"]] %in% grp & df_temp[["species"]] %in% sp,c("chr","start","end"),drop = FALSE])
        colnames(sp_df) = paste0(c("chr","start","end"),sp_idx)
        max_sp = max(table(df_temp[df_temp[["group"]] %in% grp,"species"])) - nrow(sp_df)
        if(max_sp){
          res_temp_df = rbind(sp_df,
                              as.data.frame(matrix(c("No match",NA,NA),
                                                   ncol = 3,
                                                   nrow = max(table(df_temp[df_temp[["group"]] %in% grp,"species"])) - nrow(sp_df),
                                                   byrow = TRUE,
                                                   dimnames = list(NULL,colnames(sp_df)))))
          res_df = cbind(res_df,res_temp_df)
        }else{
          res_df = cbind(res_df,sp_df)
        }
      }
    }
    res_chr = unique(res_df[,grep(pattern = "^chr",x = colnames(res_df))[-1],drop = FALSE])
    res_df = do.call(rbind,apply(res_chr,1,function(y){
      tmp = res_df[apply(res_df,1,function(x) all(x[grep(pattern = "^chr",x = colnames(res_df))[-1]] == y)),]
      res = tmp[1,,drop = FALSE]
      res[,grep(pattern = "^start",x = colnames(res_df))[-1]] = apply(tmp[,grep(pattern = "^start",x = colnames(res_df))[-1],drop = FALSE],2,min)
      res[,grep(pattern = "^end",x = colnames(res_df))[-1]] = apply(tmp[,grep(pattern = "^end",x = colnames(res_df))[-1],drop = FALSE],2,max)
      return(res)
    }
      ))
    return(res_df)
  }))
  cat("\n")
  return(res)
}))


################
block_colors = function(init_params,chr_lengths,alluvial_df){
blocks = read.table(file = init_params[["blocks"]], header = FALSE, sep = "\t", col.names = c("Block", "query", "chr", "start", "end", "length", "strand", "pident", "block_start", "block_end", "evalue"))
blocks$chr = factor(blocks$chr, levels = chr_lengths[[1]]$chr)
chr_color = data.frame(Block = unique(blocks$Block), color = "")
chr_color[chr_color$Block %in% c("A", "B", "C"), "color"] = "yellow"
chr_color[chr_color$Block %in% c("D", "E"), "color"] = "red"
chr_color[chr_color$Block %in% c("F", "G", "H"), "color"] = "blue"
chr_color[chr_color$Block %in% c("I", "J"), "color"] = "purple"
chr_color[chr_color$Block %in% c("KL1", "KL2", "MN"), "color"] = "orange"
chr_color[chr_color$Block %in% c("O", "P", "Q", "R"), "color"] = "green"
chr_color[chr_color$Block %in% c("S", "T", "U"), "color"] = "cyan"
chr_color[chr_color$Block %in% c("V", "W", "X"), "color"] = "pink"

contiguous_regions = data.frame(matrix(ncol = 5,nrow = 0,dimnames = list(NULL,c("chr","group","start","end","block"))))
for(block_var in unique(blocks$Block)){
  df_positions = blocks[blocks$Block %in% block_var, c("chr", "start", "end")]
  df_positions = df_positions %>%
    group_by(chr) %>%
    arrange(start) %>%
    mutate(group = cumsum(c(1, diff(start) > init_params[["cont_block_cutoff"]])))
  
  tmp = df_positions %>%
    group_by(chr, group) %>%
    summarise(start = min(start), end = max(end), block = block_var,.groups = "drop_last")
  tmp = tmp[abs(tmp$start - tmp$end) > init_params[["cont_block_cutoff"]],]
  contiguous_regions = rbind(contiguous_regions,tmp)
  rm(tmp)
}

link_block_colors = sapply(1:nrow(alluvial_df),function(x){
  inblock = contiguous_regions$chr %in% alluvial_df$chr1[x] &
    contiguous_regions$start < alluvial_df$start1[x] &
    contiguous_regions$end > alluvial_df$start1[x]
  if(any(inblock)){
    res = chr_color$color[chr_color$Block %in% contiguous_regions$block[inblock][1]]
    names(res) = contiguous_regions$block[inblock][1]
  }else{
    res = "grey"
    names(res) = "Non-block"
  }
  return(res)
})

return(list(chr_color,link_block_colors))
}

c(chr_color,link_block_colors) %<-% block_colors(init_params = init_params,chr_lengths = chr_lengths,alluvial_df = alluvial_df)
color_mapping = c(chr_color$Block,"Non-block")
names(color_mapping) = c(chr_color$color,"grey")

alluvial_df = cbind(data.frame(color = link_block_colors,block = factor(names(link_block_colors))),alluvial_df)

levels(alluvial_df$block) = color_mapping[color_mapping %in% levels(alluvial_df$block)]
color_mapping = color_mapping[color_mapping %in% levels(alluvial_df$block)]

aes_mapping = aes(
  y = end1 - start1
)
for (i in 1:length(init_params$species)){
  aes_mapping[[paste0("axis", i)]] <- sym(paste0("chr",i))
}

species_labels = sapply(init_params$species, function(x) {
  label <- gsub(pattern = "_", replacement = " ", x)
  bquote(italic(.(label)))
})

cat("#####   Printing alluvial plot with unmatched links   #####\n")
plot = ggplot(data = alluvial_df,aes_mapping) +
            geom_alluvium(aes(fill = block), width = 0.25) +
            geom_stratum(width = 0.25) +
            geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
            scale_x_discrete(limit = init_params$species,label = species_labels, expand = c(0.15, 0.05)) +
            scale_fill_manual(values = names(color_mapping), labels = color_mapping) +
            theme_minimal() +
            labs(title = "Alluvial plot of genomic blocks in selected species",
                 y = "Position") +
            theme(
              panel.background = element_rect(fill = "white", colour = "white"),
              plot.background = element_rect(fill = "white", colour = "white"),
              axis.text.x = element_text(angle = 45, hjust = 1)
            )

  suppressWarnings(ggsave(filename = paste0("Alluvial_plot_with_unmatched_",format(Sys.time(),'%Y%m%d_%H%M%S'),".png"),plot = plot,width = 3840/300, height = 3306/300, units = "in", dpi = 1200))

  cat("#####   Printing alluvial plot without unmatched links   #####\n")
  plot = ggplot(data = alluvial_df[complete.cases(alluvial_df),],aes_mapping) +
    geom_alluvium(aes(fill = block), width = 0.25) +
    geom_stratum(width = 0.25) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limit = init_params$species,label = species_labels, expand = c(0.15, 0.05)) +
    scale_fill_manual(values = names(color_mapping), labels = color_mapping) +
    theme_minimal() +
    labs(title = "Alluvial plot of genomic blocks in selected species",
         y = "Position") +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background = element_rect(fill = "white", colour = "white"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  suppressWarnings(ggsave(filename = paste0("Alluvial_plot_matched_",format(Sys.time(),'%Y%m%d_%H%M%S'),".png"),plot = plot,width = 3840/300, height = 3306/300, units = "in", dpi = 1200))
  
  cat("#####   Done   #####\n")