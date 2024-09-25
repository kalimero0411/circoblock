# circoblock
Create ancestral block-wise circos plots for Brassicales-related speices

```
circoblock [OPTIONS] [-n Name] [-g Sequence.fasta] [-s1 species1] [-s2 species2] [-c species1.species2.final.chain] [--gff1 species1.gff] [--gff2 species2.bed] [--lengths1 species1.lengths] [--lengths2 species2.lengths] [-o outdir]
```

## Options
| Short     | Long      | Description     |
| ------------- | ------------- | -------- |
| -n          | --name         | Name of the run  |
| -g         | --genome         | Genome FASTA file  |
| -b          | --tblastn         | tBLASTn output file (skip mapping step)  |
| -s1          | --species1         | Query species name  |
| -s2          | --species2         | Subject species name  |
| -c          | --chain         | Chain file (required)  |
|          | --gff1         | gff/gtf/bed file for species1  |
|          | --gff2         | gff/gtf/bed file for species2  |
|          | --lengths1         | lengths file for species1  |
|          | --lengths2         | lengths file for species2  |
|          | --chr_cutoff         | Cutoff for chromosome length in bp (Default 1Mbp = 1e6)  |
|          | --chain_cutoff         | Cutoff for chromosome length in bp (Default 1Kbp = 1e4)  |
|          | --cont_chain_cutoff         | Cutoff for chromosome length in bp (Default 1Kbp = 1e4)  |
|          | --cont_block_cutoff         | Cutoff for chromosome length in bp (Default 1Mbp = 1e6)  |
|          | --circos_cutoff         | Cutoff for chromosome length in bp (Default 1Mbp = 1e6)  |
| -o          | --outdir         | Output directory (Default = circoblock_Name)  |
| -t           | --threads #        | Number of CPU threads to use (Default = Detected processors or 1)  |
| -h           | --help       | Display help  |


## Examples

### Run tblastn and then circlize
```
circoblock --threads 32 --name At_vs_Ah --genome TAIR10_genome.fa --gff1 TAIR10.gff3 --gff1 species2.bed --lengths1 TAIR10.lengths --lengths2 species2.lengths -c chain_dir/At.sp2.final.chain --outdir At_sp2
```

### Run just the circlize
```
circoblock --threads 32 --name At_vs_Ah --genome TAIR10_genome.fa --tblastn At.sp2.tblastn --gff1 TAIR10.gff3 --gff1 species2.bed --lengths1 TAIR10.lengths --lengths2 species2.lengths -c chain_dir/At.sp2.final.chain --outdir At_sp2
```


## Comments
- Before running, please map the genome using make_lastz_chains (https://github.com/hillerlab/make_lastz_chains):
```
make_chains.py [input subject genome name] [input query genome name] Subject_genome.fasta Query_genome.fasta --pd chain_dir/ -f --chaining_memory 20
```
- Length files are tab-delimited: chromosome [TAB] length

