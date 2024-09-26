packages = c("circlize","R.utils","EnrichedHeatmap","rtracklayer","dplyr")

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
  must_args = c("wd","chain","blocks","species1","species2","gff1","gff2","lengths1","lengths2")
  if(!all(must_args %in% names(args))){
    help = matrix(data = c("--wd    ","Working directory path",
                           "--chain    ","Chain alignment file",
                           "--blocks    ","tBLASTn Arabidopsis block alignment file",
                           "--species1    ","Name of the first species (Query of alignment & blocks)",
                           "--species2    ","Name of the second species (Subject of alignment)",
                           "--gff1    ","Species1 gff/gtf/bed file",
                           "--gff2    ","Species2 gff/gtf/bed file",
                           "--lengths1    ","Species1 chromosome lengths file",
                           "--lengths2    ","Species2 chromosome lengths file",
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
  init_params[["chain"]] = normalizePath(args[["chain"]])
  init_params[["blocks"]] = normalizePath(args[["blocks"]])
  init_params[["species1"]] = gsub(pattern = " ",replacement = "_",args[["species1"]])
  init_params[["species2"]] = gsub(pattern = " ",replacement = "_",args[["species2"]])
  init_params[["gff1"]] = normalizePath(strsplit(args[["gff1"]],split = ",")[[1]])
  init_params[["gff2"]] = normalizePath(strsplit(args[["gff2"]],split = ",")[[1]])
  init_params[["lengths1"]] = normalizePath(args[["lengths1"]])
  init_params[["lengths2"]] = normalizePath(args[["lengths2"]])
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
  
  if("circos_cutoff" %in% names(args)){
    init_params[["circos_cutoff"]] = as.numeric(args[["circos_cutoff"]])
  }else{
    init_params[["circos_cutoff"]] = 1e6
  }
  
}else{
  suppressMessages(invisible(library(rstudioapi)))
  init_params[["wd"]] = normalizePath(rstudioapi::selectDirectory(caption = "Choose working directory:"))
  setwd(init_params[["wd"]])
  init_params[["chain"]] = normalizePath(selectFile(caption = "Select chain alignment file",path = init_params[["wd"]]))
  init_params[["blocks"]] = normalizePath(selectFile(caption = "Select tBLASTn block alignment file",path = init_params[["wd"]]))
  init_params[["species1"]] = gsub(pattern = " ",replacement = "_",readline(prompt = "Select species 1 name: "))
  init_params[["species2"]] = gsub(pattern = " ",replacement = "_",readline(prompt = "Select species 2 name: "))
  num_anno = as.numeric(readline(prompt = "number of annotations per species: "))
  init_params[["gff1"]] = sapply(1:num_anno,function(x) normalizePath(selectFile(caption = "Select species 1 gff file",path = init_params[["wd"]])))
  init_params[["gff2"]] = sapply(1:num_anno,function(x) normalizePath(selectFile(caption = "Select species 2 gff file",path = init_params[["wd"]])))
  init_params[["lengths1"]] = normalizePath(selectFile(caption = "Select species 1 chromosome lengths file",path = init_params[["wd"]]))
  init_params[["lengths2"]] = normalizePath(selectFile(caption = "Select species 2 chromosome lengths file",path = init_params[["wd"]]))
  
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
chr_lengths = list(species1 = read.table(file = init_params[["lengths1"]],header = FALSE,sep = "\t",col.names = c("chr","length")),
                   species2 = read.table(file = init_params[["lengths2"]],header = FALSE,sep = "\t",col.names = c("chr","length")))
chr_lengths = list(species1 = chr_lengths$species1[chr_lengths$species1$length >= init_params[["chr_cutoff"]],],
                   species2 = chr_lengths$species2[chr_lengths$species2$length >= init_params[["chr_cutoff"]],])

conv2bed = function(file,chrs){
  if(grepl(pattern = "[.]g[tf]f3?$",x = basename(file),perl = TRUE)){
    gtff = read.table(file,
                     quote = "",
                     header = FALSE,
                     sep = "\t",
                     col.names = c("chr","source","feature","start","end","score","strand","frame","attribute"))
    res = data.frame(chr = gtff$chr,
                     start = gtff$start,
                     end = gtff$end,
                     name = rep(".",times = nrow(gtff)),
                     score = gtff$score,
                     strand = gtff$strand,
                     source = gtff$source,
                     type = gtff$feature,
                     phase = gtff$frame,
                     attribute = gtff$attribute)
  }else{
    res = read.table(file,
                     quote = "",
                     header = FALSE,
                     sep = "\t",
                     col.names = c("chr","start","end","name","score","strand","source","type","phase","attribute"))
  }
  res = res[res$chr %in% chrs,,drop = FALSE]
  return(res)
}

cat("#####   Getting links   #####\n")
chain = import.chain(init_params[["chain"]])
genome_df = data.frame(chr = chr_lengths$species2$chr,
                       start = rep(0,length(chr_lengths$species2$chr)),
                       end = chr_lengths$species2$length)

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
                                end = chr_lengths$species1$length[match(species1_chr,chr_lengths$species1$chr)])
species2_chromInfo = data.frame(chr = species2_chr,
                            start = rep(0,times = length(species2_chr)),
                            end = chr_lengths$species2$length[match(species2_chr,chr_lengths$species2$chr)])


norm_factor = min(sum(species1_chromInfo$end), sum(species2_chromInfo$end))

species1_chromInfo$normalized_length = species1_chromInfo$end / sum(species1_chromInfo$end) * norm_factor
species2_chromInfo$normalized_length = species2_chromInfo$end / sum(species2_chromInfo$end) * norm_factor
chromInfo = rbind(species2_chromInfo,species1_chromInfo)


cat("#####   Getting annotations   #####\n")
genes_list = lapply(1:length(init_params[["gff1"]]),function(anno_idx){
                 res = list(species1 = conv2bed(file = init_params$gff1[[anno_idx]],
                                                chrs = chr_lengths$species1$chr),
                            species2 = conv2bed(file = init_params$gff2[[anno_idx]],
                                                chrs = chr_lengths$species2$chr))
                 return(res)
})

df_combined = cbind(species1_df,species2_df)
colnames(df_combined) = c("chr_to","start_to","end_to","chr","start","end")
df_combined <- df_combined %>%
  group_by(chr_to) %>%
  arrange(start_to) %>%
  mutate(group = cumsum(c(1, diff(start_to) > init_params[["cont_chain_cutoff"]])))

species1_df_reduced <- df_combined %>%
  group_by(chr,chr_to, group) %>%
  summarize(
    start = min(start),
    end = max(end),
    start_to = min(start_to),
    end_to = max(end_to),
    .groups = 'drop'
  )

species2_df = as.data.frame(species1_df_reduced[,c("chr","start","end")])
species1_df = as.data.frame(species1_df_reduced[,c("chr_to","start_to","end_to")])
colnames(species1_df) = c("chr","start","end")

################
blocks = read.table(file = init_params[["blocks"]], header = FALSE, sep = "\t", col.names = c("Block", "query", "chr", "start", "end", "length", "strand", "pident", "block_start", "block_end", "evalue"))
blocks$chr = factor(blocks$chr, levels = unique(chr_lengths$species1$chr))
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

link_colors = sapply(1:nrow(species1_df),function(x){
  inblock = contiguous_regions$chr %in% species1_df$chr[x] &
    contiguous_regions$start < species1_df$start[x] &
    contiguous_regions$end > species1_df$start[x]
  if(any(inblock)){
    res = chr_color$color[chr_color$Block %in% contiguous_regions$block[inblock][1]]
  }else{
    res = "grey"
  }
  names(res) = species1_df$chr[x]
  return(res)
})
################

### per block
for(block_select in unique(chr_color$Block)){
  cat("\r### Creating circos plot for block ",block_select," ###    ",sep = "")
  png(filename = paste0(init_params$wd,"/",init_params$species1,"_",init_params$species2,"_",block_select,"_block.png"),width = 2900,height = 2160,units = "px",res = 600)
  circos.par(gap.degree = c(rep(1, length(species2_chr) - 1), 5, rep(1, length(species1_chr) - 1),5),start.degree = -2.5)

circos.genomicInitialize(chromInfo, plotType = NULL, sector.width = chromInfo$normalized_length,)
suppressMessages(
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2),
              CELL_META$sector.index, cex = 0.4, niceFacing = TRUE,facing = "downward")
}, track.height = mm_h(3), cell.padding = c(0, 0, 0, 0), bg.border = NA)
)

suppressMessages(
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  }, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)
)

highlight.chromosome(species1_chr, col = "blue", track.index = 2)
highlight.chromosome(species2_chr, col = "red", track.index = 2)

for(anno_idx in 1:length(init_params$gff1)){
  circos.genomicDensity(genes_list[[anno_idx]],
                        window.size = init_params$circos_cutoff,
                        col = c("blue","red"),
                        track.height = mm_h(3),
                        track.margin = c(0,0),bg.border = NA)
}

tmp_block = contiguous_regions[contiguous_regions$block %in% block_select,,drop = FALSE]
df_block = apply(species1_df,MARGIN = 1,function(x) any(tmp_block$chr %in% x[1] & as.numeric(x[2]) > tmp_block$start & as.numeric(x[2]) < tmp_block$end))
circos.genomicLink(species2_df[df_block,,drop = FALSE], species1_df[df_block,,drop = FALSE], col = add_transparency(col = rep(chr_color[chr_color$Block %in% block_select,"color"],times = sum(df_block)),transparency = 0.9))

text(1, -1, bquote(italic(.(gsub(pattern = "_",replacement = " ",init_params$species2)))),cex = 0.5)
text(1, 1, bquote(italic(.(gsub(pattern = "_",replacement = " ",init_params$species1)))),cex = 0.5)
while(!is.null(dev.list())) dev.off()

circos.clear()
}
cat("\n")

### per chr
for(chr_select in unique(species1_df$chr)){
  cat("\r### Creating circos plot for chromosome ",chr_select," ###            ",sep = "")  
  png(filename = paste0(init_params$wd,"/",init_params$species1,"_",init_params$species2,"_",chr_select,"_chr.png"),width = 2900,height = 2160,units = "px",res = 600)
  circos.par(gap.degree = c(rep(1, length(species2_chr) - 1), 5, rep(1, length(species1_chr) - 1),5),start.degree = -2.5)
  
  circos.genomicInitialize(chromInfo, plotType = NULL, sector.width = chromInfo$normalized_length,)
  suppressMessages(
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2),
                  CELL_META$sector.index, cex = 0.4, niceFacing = TRUE,facing = "downward")
    }, track.height = mm_h(3), cell.padding = c(0, 0, 0, 0), bg.border = NA)
  )
  suppressMessages(
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    }, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)
  )
  
  highlight.chromosome(species1_chr, col = "blue", track.index = 2)
  highlight.chromosome(species2_chr, col = "red", track.index = 2)
  
  for(anno_idx in 1:length(init_params$gff1)){
    circos.genomicDensity(genes_list[[anno_idx]],
                          window.size = init_params$circos_cutoff,
                          col = c("blue","red"),
                          track.height = mm_h(3),
                          track.margin = c(0,0),bg.border = NA)
  }
  
  circos.genomicLink(species2_df[species1_df$chr %in% chr_select,], species1_df[species1_df$chr %in% chr_select,], col = add_transparency(col = link_colors[species1_df$chr %in% chr_select],transparency = 0.9))
  text(1, -1, bquote(italic(.(gsub(pattern = "_",replacement = " ",init_params$species2)))),cex = 0.5)
  text(1, 1, bquote(italic(.(gsub(pattern = "_",replacement = " ",init_params$species1)))),cex = 0.5)
  while(!is.null(dev.list())) dev.off()
  
  circos.clear()

}
cat("\n")

## All chrs
cat("### Create circos plot for all chromosomes and blocks ###\n",sep = "")

png(filename = paste0(init_params$wd,"/",init_params$species1,"_",init_params$species2,"_Allchr.png"),width = 2900,height = 2160,units = "px",res = 600)
circos.par(gap.degree = c(rep(1, length(species2_chr) - 1), 5, rep(1, length(species1_chr) - 1),5),start.degree = -2.5)

circos.genomicInitialize(chromInfo, plotType = NULL, sector.width = chromInfo$normalized_length,)
suppressMessages(
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2),
                CELL_META$sector.index, cex = 0.4, niceFacing = TRUE,facing = "downward")
  }, track.height = mm_h(3), cell.padding = c(0, 0, 0, 0), bg.border = NA)
)
suppressMessages(
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  }, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)
)

highlight.chromosome(species1_chr, col = "blue", track.index = 2)
highlight.chromosome(species2_chr, col = "red", track.index = 2)

for(anno_idx in 1:length(init_params$gff1)){
  circos.genomicDensity(genes_list[[anno_idx]],
                        window.size = init_params$circos_cutoff,
                        col = c("blue","red"),
                        track.height = mm_h(3),
                        track.margin = c(0,0),bg.border = NA)
}

for(chr_select in unique(species1_df$chr)){
  circos.genomicLink(species2_df[species1_df$chr %in% chr_select,], species1_df[species1_df$chr %in% chr_select,], col = add_transparency(col = link_colors[species1_df$chr %in% chr_select],transparency = 0.9))
}
text(1, -1, bquote(italic(.(gsub(pattern = "_",replacement = " ",init_params$species2)))),cex = 0.5)
text(1, 1, bquote(italic(.(gsub(pattern = "_",replacement = " ",init_params$species1)))),cex = 0.5)
while(!is.null(dev.list())) dev.off()

circos.clear()
unlink(paste0(init_params$wd,"/Rplots.pdf"))