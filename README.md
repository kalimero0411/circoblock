# circoblock
Create ancestral block-wise circos/alluvial plots for Brassicales-related speices

```
circoblock circos [OPTIONS] [-g Sequence.fasta | -b species1.tblastn] [--species1 species1] [--species2 species2] [-c species1.species2.final.chain] [--gff1 species1.gff] [--gff2 species2.bed] [--lengths1 species1.lengths] [--lengths2 species2.lengths] [-o outdir]

circoblock alluvial [OPTIONS] [-g Sequence.fasta | -b species1.tblastn] [--species species1,species2,species3] [-c species1.species2.final.chain,species1.species3.final.chain] [--lengths species1.lengths,species2.lengths,species3.lengths] [-o outdir]

circoblock chromosome [OPTIONS] [--genome Sequence.fasta | --tblastn species.Arabidopsis.tblastn] [--lengths species.lengths] [-o outdir]
```

## Circos Options
| Short     | Long      | Description     |
| ------------- | ------------- | -------- |
| -g         | --genome         | Genome FASTA file (required for mapping step)  |
| -b          | --tblastn         | tBLASTn output file (skip mapping step)  |
| -s1          | --species1         | Query species name  |
| -s2          | --species2         | Subject species name  |
| -c          | --chain         | Chain file (required)  |
| -g1         | --gff1         | gff/gtf/bed file for species1 required; multiple annotations comma separated in the same order)  |
| -g2         | --gff2         | gff/gtf/bed file for species2 (required; multiple annotations comma separated in the same order)  |
| -l1         | --lengths1         | lengths file for species1 (required) |
| -l2         | --lengths2         | lengths file for species2 (required) |
|          | --chr_cutoff         | Cutoff for chromosome length in bp (Default 1Mbp = 1e6)  |
|          | --chain_cutoff         | Cutoff for chromosome length in bp (Default 1Kbp = 1e4)  |
|          | --cont_chain_cutoff         | Cutoff for chromosome length in bp (Default 1Kbp = 1e4)  |
|          | --cont_block_cutoff         | Cutoff for chromosome length in bp (Default 1Mbp = 1e6)  |
|          | --circos_cutoff         | Cutoff for chromosome length in bp (Default 1Mbp = 1e6)  |
| -o          | --outdir         | Output directory (Default = circoblock_DateTime in the current directory)  |
| -t           | --threads #        | Number of CPU threads to use (Default = Detected processors or 1)  |
|          | --version         | Version  |
| -h           | --help       | Display help  |

## Alluvial Options
| Short     | Long      | Description     |
| ------------- | ------------- | -------- |
| -g         | --genome         | Genome FASTA file (required for mapping step)  |
| -b          | --tblastn         | tBLASTn output file (skip mapping step)  |
| -s          | --species         | Species names starting with the query species in order (required; comma-separated)  |
| -c          | --chains         | Chain files starting with the query species in order (required; comma-separated)  |
| -l         | --lengths         | Lengths files starting with the query species in order (required; comma-separated) |
|          | --chr_cutoff         | Cutoff for chromosome length in bp (Default 1Mbp = 1e6)  |
|          | --chain_cutoff         | Cutoff for chromosome length in bp (Default 1Kbp = 1e4)  |
|          | --cont_chain_cutoff         | Cutoff for chromosome length in bp (Default 1Kbp = 1e4)  |
|          | --cont_block_cutoff         | Cutoff for chromosome length in bp (Default 1Mbp = 1e6)  |
|          | --circos_cutoff         | Cutoff for chromosome length in bp (Default 1Mbp = 1e6)  |
| -o          | --outdir         | Output directory (Default = circoblock_DateTime in the current directory)  |
| -t           | --threads #        | Number of CPU threads to use (Default = Detected processors or 1)  |
|          | --version         | Version  |
| -h           | --help       | Display help  |

## Chromosome Options
| Short     | Long      | Description     |
| ------------- | ------------- | -------- |
| -g         | --genome         | Genome FASTA file (required for mapping step)  |
| -b          | --tblastn         | tBLASTn output file (skip mapping step)  |
| -l         | --lengths         | Lengths files starting with the query species in order |
|          | --chr_cutoff         | Cutoff for chromosome length in bp (Default 1Mbp = 1e6)  |
|          | --cont_block_cutoff         | Cutoff for chromosome length in bp (Default 1Mbp = 1e6)  |
| -o          | --outdir         | Output directory (Default = circoblock_DateTime in the current directory)  |
| -t           | --threads #        | Number of CPU threads to use (Default = Detected processors or 1)  |
|          | --version         | Version  |
| -h           | --help       | Display help  |

## Examples

### Run tblastn and then circlize
```
circoblock circos --threads 32 --genome TAIR10_genome.fa --species1 'Arabidopsis thaliana' --species2 'Some species' --gff1 TAIR10.gff3 --gff1 TAIR10_TEs.gtf --gff1 species2_genes.bed --gff2 species2_TEs.gff --lengths1 TAIR10.lengths --lengths2 species2.lengths -c chain_dir/At.sp2.final.chain --outdir path
```

### Run just the circlize
```
circoblock circos --threads 32 --tblastn At.sp2.tblastn --species1 'Arabidopsis thaliana' --species2 'Some species' --gff1 TAIR10.gff3 --gff1 TAIR10_TEs.gtf --gff1 species2_genes.bed --gff2 species2_TEs.gff --lengths1 TAIR10.lengths --lengths2 species2.lengths -c chain_dir/At.sp2.final.chain --outdir path
```

### Run tblastn and then alluvial
```
circoblock alluvial --threads 32 --genome species1_genome.fa --lengths 'species1.lengths,species2.lengths,species3.lengths' --chains 'chain_dir/sp1.sp2.final.chain,chain_dir/sp1.sp3.final.chain' --species 'species1,species2,species3' --outdir path
```

### Run tblastn and then chromosome
```
circoblock chromosome --threads 32 --genome species_genome.fa --lengths species.Chr.lengths --outdir path
```

## Comments
- Before running `circos` or `alluvial` , please map the genome using make_lastz_chains (https://github.com/hillerlab/make_lastz_chains):
```
make_chains.py [input subject genome name] [input query genome name] Subject_genome.fasta Query_genome.fasta --pd chain_dir/ -f --chaining_memory 20
```
- Length files are tab-delimited: chromosome [TAB] length
- If `--genome` was provided and the mapping set completed, the output can be used again with `--tblastn [species1].blocks.tblastn`

# Circos
- GFF files should only include lines that should be plotted. (e.g. awk -F $'\t' '$3 == "gene"{print $0}' species.gff) 

# Alluvial
- There should be an equal number (N) of species and lengths files in the same order, and N-1 chain files in the same order (the query species against every other species).

# Chromosome
- The script outputs probabilities of ancestral chromosome Ancestral Crucifer Karyotype (ACK) and Clade E Karyotype (CEK) association with the input chromosomes in text tables, as well as graphical maps, for primary (highest probability) and secondary (second-highest probability) matches.
