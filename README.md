# Comparative_genomics_JUNE
This repository details how to run the analyses and create the figures for the comparative genomics chapter. 

# FIG 2 - RE/TE abundance and diversity, length, and TE-associated BUSCOs
I identified transposable elements from the un-masked genome assemblies using EarlGrey - https://github.com/TobyBaril/EarlGrey. Earl grey generates longer, non-redundant TE consensus sequences than older software. It feeds the input genome into RepeatMasker and compares it against the Dfam database to identify well-characterised repeat families. It also runs RepeatModeler2 on any remaining unmasked genomic regions to find novel elements unique to the organism and builds initial bases consensus sequences for newly discovered families. It uses CD-HIT to cluster and collapse duplicates, and merges overlapping definitions into a single, high-fidelity, non-redundant TE library. It finally performans a final RepeatMasker pass over to resolve any lingering overlapping annotations and has a lot of 'paper-ready' tools such as TE landscape plots and GC%. 

I ran this code for each of the six blowfly species: 

```
#!/bin/bash
#SBATCH --job-name=EG_hilli
#SBATCH --account=uow03920
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=7-00:00:00
#SBATCH --output=eg_hilli_%j.out
#SBATCH --error=eg_hilli_%j.err

source /home/pm101/miniconda3/etc/profile.d/conda.sh
conda activate earlgrey
for i in 01_hilli; do
cd ${i}
earlGrey -g ${i}.fa -s 01_hilli -o ${i}_earlgrey -r Diptera -t $SLURM_CPUS_PER_TASK;
done
```

To further improve the quality of the TE library and annotation, I ran TEtrimmer on the EarlGrey outputs. TEtrimmer relies on a core module called TEstrainer to build its library. TEtrimmer fixes divergent family over-simplification, cleans multiple sequence alignments using a sliding window strategy (EarlGrey can sometimes leave behind messy, poorly conseved regions in its final consensus sequences), and pinpoints exact terminal coordinates of the TE (which can sometimes be fuzzy following Earl Grey's elongation loop). 

Run this for each species: 

```
#!/bin/bash
#SBATCH --job-name=TEtrimmer_hilli
#SBATCH --account=uow03920
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=TEtrimmer_%j.log

conda activate /nesi/nobackup/uow03920/conda_envs/tetrimmer

TEtrimmer \
  --input_file 01_hilli_combined_library.fasta \
  --genome_file ../../37_TE_PROPER/02_earlgrey/01_hilli/01_hilli.fa \
  --output_dir 01_hilli_TEtrimmer_out \
  --num_threads 16 \
  --classify_all
```

And then, run RepeatMasker to get an updated TE summary (make sure to run with the -a flag if you are wanting to calculate kimura divergence, as you need the aligned file to do this). 

```
#!/bin/bash
#SBATCH --job-name=repeatmasker_hilli
#SBATCH --account=uow03920
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=repeatmasker_hilli_%j.log

module purge 
module load RepeatMasker 

RepeatMasker \
  -pa 16 \
  -a \
  -gff \
  -lib 01_hilli_TEtrimmer_out/TEtrimmer_consensus_merged.fasta \
  ../../37_TE_PROPER/02_earlgrey/01_hilli/01_hilli.fa
  ```

The above code will produce .tbl files that contain summary information on the % of genome comprised of various transposable element families (which I used to create Fig. 2A). 

## Fig. 2B - TE length 
To calculate the TE length, I used the following R code:

```
library(dplyr)
library(ggplot2)
library(readr)

# skip the 3 header lines
rm <- read_table(
  "02_quadrimaculata.fa.out",
  skip = 3,
  col_names = FALSE
)

colnames(rm) <- c(
  "SW_score",
  "perc_div",
  "perc_del",
  "perc_ins",
  "chr",
  "start",
  "end",
  "left_query",
  "strand",
  "repeat_name",
  "class_family",
  "repeat_begin",
  "repeat_end",
  "repeat_left",
  "ID",
  "flag"
)

# inspect
glimpse(rm)

rm <- rm %>%
  mutate(length = end - start + 1)

sort(table(rm$class_family), decreasing = TRUE)[1:50]

rm <- rm %>%
  mutate(
    TE_order = case_when(
      grepl("^LTR", class_family, ignore.case = TRUE) ~ "LTR",
      grepl("^LINE", class_family, ignore.case = TRUE) ~ "LINE",
      grepl("^SINE", class_family, ignore.case = TRUE) ~ "SINE",
      grepl("^DNA", class_family, ignore.case = TRUE) ~ "DNA",
      grepl("Helitron|RC", class_family, ignore.case = TRUE) ~ "Rolling_circle",
      grepl("Unknown", class_family, ignore.case = TRUE) ~ "Unclassified",
      TRUE ~ "Other"
    )
  )

table(rm$TE_order)

rm %>%
  group_by(TE_order) %>%
  summarise(
    n = n(),
    mean_length = mean(length),
    median_length = median(length)
  ) %>%
  arrange(desc(mean_length))


#plot all species together on a figure
library(tidyverse)

# Create dataframe
df <- data.frame(
  Species = c("HIL", "QUA", "STY",
              "VIC", "SER", "MEG"),
  LTR = c(280, 653, 227, 235, 187, 437),
  LINE = c(200, 170, 169, 202, 179, 199),
  DNA = c(175, 248, 184, 189, 122, 191),
  Other = c(40, 41, 43, 44, 52, 47),
  Rolling_circle = c(234, 236, 214, 145, 202, 84),
  Unclassified = c(163, 132, 142, 152, 91, 215)
)

# Convert to long format
df_long <- df %>%
  pivot_longer(
    cols = -Species,
    names_to = "TE_Order",
    values_to = "Median_Length"
  )

# Set species order
df_long$Species <- factor(
  df_long$Species,
  levels = c("QUA", "STY", "HIL", "VIC", "SER", "MEG")
)

# Set order of TE categories
df_long$TE_Order <- factor(
  df_long$TE_Order,
  levels = c("LTR", "LINE", "DNA", "Other", "Rolling_circle", "Unclassified")
)


# Plot
length <- ggplot(df_long,
       aes(x = Species,
           y = Median_Length,
           fill = TE_Order)) +
  
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7) + 
  
  scale_fill_manual(values = c(
    "LTR" = "#9385ED",
    "LINE" = "#85EDC7",
    "DNA" = "#85DFED",
    "Other" = "#fab27f",
    "Rolling_circle" = "#E0ED85",
    "Unclassified" = "#ED85DF"
  )) +
  
  theme_classic() +
  
  labs(
    x = element_blank(),
    y = "Median RE Length (bp)",
    fill = "TE Order"
  ) +
  
  theme(
    axis.text.x = element_text(
      size = 18, colour = "black"
    ),
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title = element_text(size = 18),
    legend.position = 'right',
    axis.line = element_line(linewidth = 0.5)
  )
length

ggsave("length.png", plot = length, dpi = 400)
``` 
It calculates both mean and median TE length. I chose to present median, since mean was likely influenced by a few elements with extremely long/short length. 

## Fig. 2C - TE-associated BUSCOs.
I got the idea for this analysis from https://pmc.ncbi.nlm.nih.gov/articles/PMC12590894/. The idea is that you can use BUSCO genes as a proxy for detecting genome-wide TE-gene associations. More scripts/information can be found at https://pubmed.ncbi.nlm.nih.gov/35217860/ and https://pubmed.ncbi.nlm.nih.gov/37739812/. 

The goal is to find the intersect of the co-ordinates of BUSCO genes and co-ordinates of TEs. The co-ordinates of TEs are in the .gff file and BUSCO are in the full_table.tsv file.  

### 1. First, I ran BUSCO on the un-masked genomes: 
```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=busco_for_TE_analysis
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=24G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output busco_%j.out
#SBATCH --error busco_%j.err

module purge

module load BUSCO/5.8.2-gimkl-2022a
ml miniprot

for i in 01_hilli 02_quadrimaculata 03_stygia 04_vicina 05_sericata 06_megacephala; do 

cd ${i}

busco \
  -i ${i}.fa \
  -o ${i}_busco \
  -l diptera_odb10 \
  -m genome \
  --long \
  -c 24
cd ../
done
```
I then ran the following bash script. Note that I only ran this on single-copy, complete BUSCO genes (not duplicated, incomplete, or missing). 

```
#!/bin/bash
module load BEDTools

SPECIES=(
    01_hilli
    02_quadrimaculata
    03_stygia
    04_vicina
    05_sericata
    06_megacephala
)

echo "Species,Total_BUSCOs,BUSCOs_with_TE,Percent_with_TE"

for s in "${SPECIES[@]}"; do
    echo "Processing $s" >&2

    #build single-copy BUSCO BED from full_table.tsv
    awk -F'\t' '
    NR>1 && $2=="Complete" {
        id=$1
        if (!(id in count)) {
            chr[id]=$3
            s=$4; e=$5
            if (s=="" || e=="") next
            start[id] = (s < e) ? s : e
            end[id]   = (s < e) ? e : s
        }
        count[id]++
    }
    END {
        for (id in count)
            if (count[id]==1) {
                out = start[id] - 1
                if (out < 0) out = 0
                print chr[id], out, end[id], id
            }
    }' OFS='\t' "${s}_full_table.tsv" \
    | sort -k1,1 -k2,2n \
    > "${s}_busco.sorted.bed"

    #convert RepeatMasker GFF to sorted BED
    grep -v "^#" "${s}.fa.out.gff" \
    | awk 'BEGIN{OFS="\t"}{ print $1, $4-1, $5, $3, $6, $7 }' \
    | sort -k1,1 -k2,2n \
    > "${s}_TE.sorted.bed"

    #find the intersect
    bedtools intersect \
        -a "${s}_busco.sorted.bed" \
        -b "${s}_TE.sorted.bed" \
        -wa -wb \
        > "${s}_BUSCO_TE_overlap.bed"

    #count number of intersects
    TOTAL_BUSCO=$(wc -l < "${s}_busco.sorted.bed")

    if [ "$TOTAL_BUSCO" -eq 0 ]; then
        echo "${s},0,0,NA"
        continue
    fi

    BUSCO_WITH_TE=$(
        cut -f1,2,3 "${s}_BUSCO_TE_overlap.bed" \
        | sort -u \
        | wc -l
    )

    PERCENT=$(awk -v a="$BUSCO_WITH_TE" -v b="$TOTAL_BUSCO" \
        'BEGIN{ printf "%.2f", 100*a/b }')

    echo "${s},${TOTAL_BUSCO},${BUSCO_WITH_TE},${PERCENT}"
done
```
This will give you the number of BUSCOs analysed, how many shared coordinates with TEs, and % of the genome this took up. I also used this information to make the RE content vs genome size plots in Fig. 3. 

# Fig. 4 - Kimura divergence/TE landscape plots
This analysis uses the .fa.align files produced by including the -a flag in RepeatMasker earlier. 

First, I used the script `CalcDivergenceFromAlign.pl` script from RepeatMasker (https://github.com/Dfam-consortium/RepeatMasker/blob/master/util/calcDivergenceFromAlign.pl). I put the script in the same working directory as my .fa.align files and ran:

```perl calcDivergenceFromAlign.pl -s 01_hilli.divsum 01_hilli.fa.align```

This generates a .divsum file which can be used in createRepeatLandscape.pl, also implemented in RepeatMaskers github page (https://github.com/Dfam-consortium/RepeatMasker/blob/master/util/createRepeatLandscape.pl). You need the genome size for this. I edited this perl script so that it would generate csv files rather than html files, as I wanted to customise the plots (I've added this customised script in this repository -- it's called 'createRepeatLandscape_generate_csv.pl').

I used this R code to generate the kimura plots for each species:

```
library(tidyverse)

df <- read_csv("01_hilli_landscape.csv")

#make long
long_df <- df %>%
  pivot_longer(
    -Divergence,
    names_to = "Class",
    values_to = "GenomePercent"
  )

#i don't want to plot every single TE family, group by higher level classifications
grouped_df <- long_df %>%
  mutate(
    Group = case_when(
      str_detect(Class, "^LINE") ~ "LINEs",
      str_detect(Class, "^LTR") ~ "LTRs",
      str_detect(Class, "^DNA") ~ "DNA transposons",
      Class == "RC/Helitron" ~ "Rolling circle",
      Class == "Unknown" | Class == "Other" ~ "Unclassified",
      TRUE ~ "Unclassified"
    )
  ) %>%
  group_by(Divergence, Group) %>%
  summarise(
    GenomePercent = sum(GenomePercent),
    .groups = "drop"
  )

#group by levels
grouped_df$Group <- factor(
  grouped_df$Group,
  levels = c(
    "DNA transposons",
    "LINEs",
    "LTRs",
    "Rolling circle",
    "Unclassified"
  )
)

#colours 
fill_colours <- c(
  "DNA transposons" = "#e69f00",
  "LINEs" = "#009e73",
  "LTRs" = "#0072b2",
  "Rolling circle" = "#f0e442",
  "Unclassified" = "#cc79a7"
)


ggplot(grouped_df, aes(x = Divergence, y = GenomePercent, fill = Group)) +
  geom_col(width = 1, color = "black", linewidth = 0.05) +
  scale_x_continuous(
    limits = c(-0.5, 50.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  labs(
    x = "Kimura substitution level (CpG adjusted)",
    y = "Percent of genome",
    fill = NULL
  ) +
  scale_fill_manual(values = fill_colours) +
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.ticks = element_line(colour = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )
```
# Fig. 5 - GC content per TE family
This uses the .fa.out file that we generated from RepeatMasker. First have to edit it a bit:

### remove the first three lines from the .fa.out file, which is just header, column label, and blank spacing
```awk 'NR>3' 06_megacephala.fa.out > 06_megacephala.noheader.out```

### Turn the file into a BEDfile and rearrange the columns a little bit to fit the BED format rules
```
grep -v "^#" 06_megacephala.noheader.out | awk 'BEGIN{OFS="\t"} {
  print $5, $6-1, $7, $11, $12
}' > 06_megacephala_TE.bed
```

This one simplifies the TE classification column in the TE bed file by changing them to higher level classifications (eg LTR, LINE, DNA, etc.)
```
awk 'BEGIN{OFS="\t"} {
  print $5, $6-1, $7, $10, $11
}' 06_megacephala.noheader.out > 06_megacephala.TE.class.bed
```
Remove invalid genomic intervals where the start position is greater or equal to the end position - ensures every single region in the file has a positive, non-zero length. 
```
awk '$3>$2' 06_megacephala.TE.class.fixed.bed > 06_megacephala.TE.class.fixed.nonzero.bed
```
Use BEDTools to calculate the GC content for each of the genomic regions (TEs) specified in the BED file. 
```
bedtools nuc   -fi ../37_TE_PROPER/02_earlgrey/06_megacephala/06_megacephala.fa   -bed 06_megacephala.TE.class.fixed.nonzero.bed   > 06_megacephala_TE.gc.tsv
```
## Then, use this .tsv file to plot in R using the following code:

```
library(tidyverse)

gc <- read_tsv(
  "06_megacephala_TE.gc.tsv",
  skip = 1,   # skip the header line beginning with #1_usercol
  col_names = c(
    "chr",
    "start",
    "end",
    "family",
    "class",
    "pct_at",
    "pct_gc",
    "A",
    "C",
    "G",
    "T",
    "N",
    "other",
    "length"
  ),
  show_col_types = FALSE
)

#classify orders
gc <- gc %>%
  mutate(
    TE_order = case_when(
      class == "LTR" ~ "LTRs",
      class == "LINE" ~ "LINEs",
      class == "DNA" ~ "DNA transposon",
      grepl("RC|PLE|Helitron", class, ignore.case = TRUE) ~ "Rolling circle",
      TRUE ~ "Unclassified"
    )
  )

plot_df <- gc %>%
  mutate(
    GC_percent = pct_gc * 100
  ) %>%
  filter(!is.na(GC_percent))

# colours
my_colors <- c(
  "DNA transposon" = "#e69f00",
  "LINEs" = "#009e73",
  "LTRs" = "#0072b2",
  "Rolling circle" = "#f0e442",
  "Unclassified" = "#cc79a7"
)

plot_df$TE_order <- factor(
  plot_df$TE_order,
  levels = c(
    "DNA transposon",
    "LINEs",
    "LTRs",
    "Rolling circle",
    "Unclassified"
  )
)

#PLOT
gc_plot <- ggplot(plot_df,
                  aes(x = GC_percent,
                      fill = TE_order)) +
  
  geom_density(
    alpha = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  scale_fill_manual(values = my_colors) +
  
  scale_x_continuous(
    breaks = seq(0, 100, 10),
    limits = c(0, 100)
  ) +
  
  theme_classic() +
  
  labs(
    x = "GC (%)",
    y = "Density",
    fill = "TE class"
  ) +
  
  theme(
    legend.position = "none",
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  
  )

gc_plot
```

# Fig. 7 - Expanding/contracting gene families

First, predict proteins using homology based evidence from protein databases (Uniprot-Swissprot and BUSCO proteins from the diptera set). I used BRAKER3 to predict proteins:

```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=braker_mega
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output=braker_mega%j.out
#SBATCH --error=braker_mega%j.err

#Clean environment
module purge

#Load apptainer module
ml Apptainer
ml AUGUSTUS

for i in 13_mega; do

cd /nesi/nobackup/uow03920/06_blowfly_assembly_jan/36_orthofinder/braker/${i}


#link in the genome assembly with repetitive regions masked and also the proteins fasta we generated from ortho and busco
ln -s /nesi/nobackup/uow03920/06_blowfly_assembly_jan/11_braker/proteins.fasta ./proteins.fasta

#Run BRAKER3
apptainer exec \
  --bind /nesi/nobackup/uow03920/06_blowfly_assembly_jan/36_orthofinder/braker/${i}/augustus_config:/opt/Augustus/config \
  /nesi/nobackup/uow03920/06_blowfly_assembly_jan/11_braker/genemark/braker3.sif braker.pl \
    --threads=12 \
    --genome=${i}.fa \
    --prot_seq=proteins.fasta \
    --species=${i} \
    --gff3 \
    --AUGUSTUS_ab_initio \
    --crf

cd ../; 
done
```

You can check the quality of these predictions using BUSCO with the protein function flagged:
```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=busco_clean_proteins
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output busco_clean_%j.out
#SBATCH --error busco_clean_%j.err

module purge
module load compleasm/0.2.5-gimkl-2022a

for i in 13_mega 14_sericata; do

# ---- Input files ----
LINEAGE=/nesi/nobackup/uow03920/05_blowfly_assembly_march/29_annotation_qc/diptera_odb10
OUTDIR=/nesi/nobackup/uow03920/06_blowfly_assembly_jan/36_orthofinder/braker/${i}

mkdir -p $OUTDIR

# ---- Step 3: Run BUSCO (protein mode) on the cleaned proteome ----
compleasm.py protein \
  -p /nesi/nobackup/uow03920/06_blowfly_assembly_jan/36_orthofinder/braker/${i}_longest_no_dup.fa \
  -o $OUTDIR/busco_protein_longest \
  -l $LINEAGE \
  -t 8;
done
```

Select the longest isoforms using the following bash script:

```#!/bin/bash
# List of sample IDs
samples=(05_sericata 06_megacephala)
# Loop through each sample
for sample in "${samples[@]}"
do
    echo "Processing $sample..."

    python3 <<EOF
from Bio import SeqIO
from collections import defaultdict

input_fasta = "/nesi/nobackup/uow03920/06_blowfly_assembly_jan/36_orthofinder/braker/${sample}/braker/braker.aa_modified.faa"
output_fasta = "${sample}_longest_isoforms.fa"

# Store longest isoform per gene
longest = defaultdict(lambda: ("", 0))  # gene_id → (record, length)

for record in SeqIO.parse(input_fasta, "fasta"):
    gene_id = record.id.split(".")[0]
    if len(record.seq) > longest[gene_id][1]:
        longest[gene_id] = (record, len(record.seq))

with open(output_fasta, "w") as out:
    for rec, _ in longest.values():
        SeqIO.write(rec, out, "fasta")

print(f"Finished ${sample}")
EOF
done
```
And remove duplicate proteins:
```seqkit rmdup -s 06_cuprina_longest_isoforms_modified.faa > 06_cuprina_no_dup_longest_isoforms.fa```

## Run OrthoFinder!
I provided OrthoFinder with my own Newick Tree because orthofinder wasn't producing an accurate species tree. I believe this is because the blowflies are too closely related, so it gets a bit messy. 

Orthofinder is used to infer orthogroups (genes descended from a single gene in the last common ancestor), orthologs, rooted gene trees, and rooted species trees. It is essential for identifying gene duplication events, analyse gene family evolution (what I used it for), and compare protein sets. 

My Newick tree was called tree.nwk and was in the same folder as OrthoFinder:
```(06_megacephala_protein:30,(05_sericata_protein:19.75186453,(((01_hilli_protein:5.066546484,03_stygia_protein:5.066546484)1:4.01645853,02_quadrimaculata_protein:9.083005013)1:4.536153609,04_vicina_protein:13.61915862)1:6.132705906)1:10.24813547);
```

Run OrthoFinder. 'proteins' is a folder containing all of the protein files generated in previous steps (e.g., 01_hilli_longest_isoforms_modified.faa)

```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=orthofinder
#SBATCH --time=172:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output orthofinder_%j.out
#SBATCH --error orthofinder_%j.err

module purge
ml OrthoFinder/2.5.2 MAFFT/7.505 IQ-TREE/2.2.2.2

# Run:
orthofinder -f proteins -M msa -T iqtree -t 8 -a 4```

I mainly used OrthoFinder to generate files that could be used in CAFE3

## Gene family expansion and contraction using CAFE3
### Convert our Newick Tree to Ultramentric tree (chronos method) in R. You need to estimate total root to tip time. I used 30 MYA, as this is when all species shared a common ancestor. This part is done in R. 

```
library(ape)
# Load the rooted species tree
tree <- read.tree("SpeciesTree_rooted.txt")

# Check that it's rooted and binary
stopifnot(is.rooted(tree))
stopifnot(is.binary(tree))

# Convert to ultrametric with a relaxed clock
ultra_tree <- chronos(tree, lambda = 1, model = "correlated")

# Scale tree so root-to-tip distance = 30 Mya
tree_age <- max(node.depth.edgelength(ultra_tree))
ultra_tree$edge.length <- ultra_tree$edge.length * (30 / tree_age)

# Plot and save
plot(ultra_tree, main = "Ultrametric Tree (Root Scaled to 30 Mya)")
write.tree(ultra_tree, file = "SpeciesTree_ultrametric_30Mya.txt")  
``` 

### Filter the orthogroup gene count file (from OrthoFinder) and modify it for CAFE3
cafetutorial_clade_and_size_filter.py is available here: https://github.com/hahnlab/CAFE5/blob/master/docs/tutorial/clade_and_size_filter.py. 

```
awk 'OFS="\t" {$NF=""; print}' Orthogroups.GeneCount.tsv > tmp \
&& awk '{print "(null)""\t"$0}' tmp > cafe.input.tsv \
&& sed -i '1s/(null)/Desc/g' cafe.input.tsv \
&& rm tmp

python2 cafetutorial_clade_and_size_filter.py \
    -i cafe.input.tsv -o filtered.cafe.input.tsv -s 2
```
## Run CAFE3 using these inputs!
```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=CAFE
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output CAFE_%j.out    # save the output into a file
#SBATCH --error CAFE_%j.err     # save the error output into a file

#purge
module purge

ml GCC/12.3.0 
ml OpenBLAS

# Set variables
CAFE5_EXEC="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/32_cafe/CAFE5/build/cafe5"
CAFE_DIR="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/38_cafeagain/"
ORTHO_DIR="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/36_orthofinder/proteins/OrthoFinder/Results_May07"
TREE="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/38_cafeagain/SpeciesTree_ultrametric_17Mya.txt"
COUNTS="/nesi/nobackup/uow03920/06_blowfly_assembly_jan/38_cafeagain/cafe_input_filtered.tsv"
K_VALS=(1 2 3 4 5 6)
REPLICATES=(1 2)
PVALUE_THRESHOLD=0.05

# Step 1: Run CAFE5 with multiple k values and replicates
for k in "${K_VALS[@]}"; do
    for rep in "${REPLICATES[@]}"; do
        OUTDIR="$CAFE_DIR/k${k}_rep${rep}"
        mkdir -p "$OUTDIR"
        echo "Running CAFE5 for k=$k, replicate=$rep..."
        $CAFE5_EXEC -i "$COUNTS" -t "$TREE" -o "$OUTDIR" -k $k > "$OUTDIR/stdout_log.txt"
    done
done

# Step 2: Extract log-likelihoods, check convergence, find best k
best_k=""
best_lnl=999999

echo -e "\n--- Log-likelihoods and convergence check ---"
for k in "${K_VALS[@]}"; do
    lnls=()
    for rep in "${REPLICATES[@]}"; do
        OUTDIR="$CAFE_DIR/k${k}_rep${rep}"
        LNL=$(grep "Final -lnL:" "$OUTDIR/stdout_log.txt" | sed -E 's/.*Final -lnL: //')
        echo "k=$k, rep=$rep: lnL=$LNL"
        lnls+=("$LNL")
    done

    diff=$(echo "${lnls[0]} - ${lnls[1]}" | bc -l | tr -d '-')
    echo "k=$k: ΔlnL between replicates = $diff"

    # Use lnL from first replicate for best model selection
    cmp=$(echo "${lnls[0]} < $best_lnl" | bc)
    if [[ "$cmp" -eq 1 ]]; then
        best_k=$k
        best_lnl=${lnls[0]}
    fi
done
echo -e "\n Best k value is: $best_k (replicate 1) with lnL=$best_lnl"
BEST_RESULT_PATH="$CAFE_DIR/k${best_k}_rep1"
```

Can then use cafeplotter to generate summaries (https://github.com/moshi4/CafePlotter).

## Annotate the expanding and contraction gene families
I used EggNog and Interpro to annotate my protein files previously:

```#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=eggnog
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=15G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --output=eggnog_%j.out
#SBATCH --error=eggnog_%j.err

# Clean environment
module purge

# Load EggNOG-mapper module
module load eggnog-mapper/2.1.12-gimkl-2022a

# Loop through ecotypes for functional annotation
for i in 06_megacephala; do
  cd ${i}

  # Optional: Create softlink to protein FASTA file
 ln -s /nesi/nobackup/uow03920/06_blowfly_assembly_jan/36_orthofinder/braker/${i}/braker/braker.aa ${i}_protein.fa

  # Copy GFF3 file for functional annotation decoration
 cp /nesi/nobackup/uow03920/06_blowfly_assembly_jan/36_orthofinder/braker/${i}/braker/braker.gff3 ./${i}_braker.gff3

  # Run EggNOG-mapper
  emapper.py \
    -i ${i}_protein.fa \
    -o ${i}_eggnog \
    --data_dir /nesi/nobackup/uow03920/06_blowfly_assembly_jan/12_eggnogmapper/eggnog-mapper/data \
    --decorate_gff ${i}_braker.gff3 \
    --cpu 22
done
```
```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=interproscan_06
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=30G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --output=eggnog_%j.out    # Save standard output to file
#SBATCH --error=eggnog_%j.err     # Save error output to file

# Clean environment
module purge

# Load required modules
module load InterProScan/5.66-98.0-gimkl-2022a 
ml Perl
ml Python 
ml PCRE2/10.42-GCCcore-12.3.0 
ml GCC

for i in 05_sericata 06_megacephala; do
cd ${i}

# Run InterProScan
interproscan.sh \
  -i /nesi/nobackup/uow03920/06_blowfly_assembly_jan/36_orthofinder/braker/${i}/braker/braker.aa_modified.faa \
  -t p \
  -dp \
  -appl Pfam,SMART,PrositeProfiles,Gene3D,SUPERFAMILY,PANTHER \
  --goterms \
  --pathways \
  --iprlookup \
  -cpu 22 \
  -T temp/

cd ../
done
```
Then, I used these annotations to annotate the expanding gene families in each of the species. I used scripts that were generated by Meeran Hussain: https://github.com/meeranhussain/Comparative_genomics. This is done in R.  

## Merge the eggnog + interproscan annotations for each species
```
library(tidyverse)
library(dplyr)


samples <- c("01_hilli", "02_quadrimaculata", "03_stygia", "04_vicina", "05_sericata", "06_megacephala")

for (sample in samples) {
  # File paths
  eggnog_gff <- file.path(sample, paste0(sample, "_eggnog.emapper.decorated_modified.gff"))
  interpro_tsv <- file.path(sample, paste0(sample, "_longest_isoforms_modified.faa.tsv"))
  
  cat("Processing:", sample, "\n")
  
  #---- 1. Parse EggNOG-mapper GFF3 ----#
  eggnog <- read_tsv(
    eggnog_gff,
    comment = "##",
    col_names = FALSE
  )
  
  # Filter mRNA rows where em_* annotations are present
  eggnog_mrna <- eggnog %>% filter(X3 == "mRNA")
  
  #---- 2. Extract key fields from attribute column using regex ----#
  eggnog_clean <- eggnog_mrna %>%
    transmute(Protein_ID = str_extract(X9, "(?<=ID=)[^;]+"),
              GO_terms   = str_extract(X9, "(?<=em_GOs=)[^;]+"),
              KEGG_paths = str_extract(X9, "(?<=em_KEGG_Pathway=)[^;]+"),
              KEGG_KO    = str_extract(X9, "(?<=em_KEGG_ko=)[^;]+"),
              BRITE_class= str_extract(X9, "(?<=em_BRITE=)[^;]+"),
              EggNOG_PFAMs = str_extract(X9, "(?<=em_PFAMs=)[^;]+"),
              EggNOG_PFAMs_desc = str_extract(X9, "(?<=em_desc=)[^;]+"))
  
  # Remove rows with all NA or Protein_ID missing
  eggnog_clean <- eggnog_clean %>%
    filter(!is.na(Protein_ID)) %>%
    group_by(Protein_ID) %>%
    summarise(across(everything(), ~ toString(na.omit(.))), .groups = "drop") %>%
    mutate(
      ID_number = as.numeric(str_extract(Protein_ID, "(?<=g)\\d+"))
    ) %>%
    arrange(ID_number) %>%
    select(-ID_number)
  
  #---- 2. Parse InterProScan TSV ----#
  interpro <- read_tsv(interpro_tsv, col_names = FALSE, comment = "#")
  
  # 1. Pfam domains and descriptions (only from source = "Pfam")
  pfam_data <- interpro %>%
    filter(X4 == "Pfam") %>%
    select(Protein_ID = X1, Pfam_ID = X5, Pfam_desc = X6, Pfam_Eval = X9) %>%
    mutate(Pfam_Eval = as.numeric(Pfam_Eval)) %>%
    group_by(Protein_ID, Pfam_ID, Pfam_desc) %>%
    summarise(Min_Eval = min(Pfam_Eval, na.rm = TRUE), .groups = "drop") %>%
    group_by(Protein_ID) %>%
    summarise(
      Pfam_domains = toString(Pfam_ID),
      Pfam_descriptions = toString(Pfam_desc),
      Pfam_Evals = toString(Min_Eval),
      .groups = "drop"
    )
  
  #### Panther fetch
  panther_data <- interpro %>%
    filter(X4 == "PANTHER") %>%
    select(Protein_ID = X1, Panther_ID = X5, Panther_desc = X6, Panther_Eval = X9) %>%
    mutate(Panther_Eval = as.numeric(Panther_Eval)) %>%
    group_by(Protein_ID, Panther_ID, Panther_desc) %>%
    summarise(Min_Eval = min(Panther_Eval, na.rm = TRUE), .groups = "drop") %>%
    group_by(Protein_ID) %>%
    summarise(
      Panther_ID = toString(Panther_ID),
      Panther_descriptions = toString(Panther_desc),
      Panther_Evals = toString(Min_Eval),
      .groups = "drop"
    )
  
  # 2. InterPro IDs and descriptions (from all sources)
  ipr_data <- interpro %>%
    select(Protein_ID = X1, InterPro_ID = X12, InterPro_desc = X13) %>%
    filter(!is.na(InterPro_ID)) %>%
    group_by(Protein_ID) %>%
    summarise(
      InterPro_IDs = toString(unique(na.omit(InterPro_ID))),
      InterPro_descriptions = toString(unique(na.omit(InterPro_desc))),
      .groups = "drop"
    )
  
  # 3. GO terms (column X14), clean up GO IDs
  go_data <- interpro %>%
    select(Protein_ID = X1, GO = X14) %>%
    filter(!is.na(GO)) %>%
    separate_rows(GO, sep = ",\\s*") %>%                              # Split multiple GO term blocks
    mutate(GO = str_extract_all(GO, "GO:\\d{7}")) %>%                 # Extract only GO:#######
  unnest(GO) %>%
    group_by(Protein_ID) %>%
    summarise(
      GO_terms_interpro = toString(unique(GO)),
      .groups = "drop"
    )
  
  # 4. Merge InterProScan pieces
  interpro_clean <- reduce(
    list(pfam_data, panther_data, ipr_data, go_data),
    full_join,
    by = "Protein_ID"
  ) %>%
    mutate(
      ID_number = as.numeric(str_extract(Protein_ID, "(?<=g)\\d+"))
    ) %>%
    arrange(ID_number) %>%
    select(-ID_number)
  
  # 5. Merge with EggNOG annotations (GO terms remain separate)
  merged_annotations <- full_join(eggnog_clean, interpro_clean, by = "Protein_ID") %>%
    select(
      Protein_ID,
      EggNOG_PFAMs,
      EggNOG_PFAMs_desc,
      KEGG_paths,
      KEGG_KO,
      BRITE_class,
      Pfam_domains,
      Pfam_descriptions,
      Pfam_Evals,
      Panther_ID,
      Panther_descriptions,
      Panther_Evals,
      InterPro_IDs,
      InterPro_descriptions,
      GO_terms,               # from EggNOG
      GO_terms_interpro      # cleaned InterProScan GO terms
    )
  
  # 6. Clean NA and export
  merged_annotations <- merged_annotations %>%
    mutate(across(everything(), ~str_replace_all(., "NA|^, | ,|, NA", "")))
  
  # Reorder columns in the desired order
  merged_annotations <- merged_annotations %>%
    select(
      Protein_ID,
      EggNOG_PFAMs,
      EggNOG_PFAMs_desc,
      GO_terms,             
      KEGG_KO,
      KEGG_paths,
      BRITE_class,
      Pfam_domains,
      Pfam_descriptions,
      Pfam_Evals,
      InterPro_IDs,
      InterPro_descriptions,
      GO_terms_interpro,
      Panther_ID,
      Panther_descriptions,
      Panther_Evals,
    )
  colnames(merged_annotations) <- c("Protein_ID", "EggNOG_PFAMs_domains","EggNOG_PFAMs_descriptions", "EggNOG_GO_terms", "EggNOG_KEGG_KO", "EggNOG_KEGG_paths", "EggNOG_BRITE_class", "Interpro_Pfam_Ids",
                                    "Interpro_Pfam_descriptions", "Interpro_Pfam_Evals", "InterPro_IDs", "InterPro_descriptions", "Interpro_GO_terms", "Panther_ID", "Panther_descriptions", "Panther_Evals")
  
  # ---- 4. Write to file ----
  output_file <- file.path(sample, "merged_protein_annotations.tsv")
  write_tsv(merged_annotations, output_file)
}
```

Annotate the merged annotation file with orthogroups from orthofinder:
```
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# === File paths ===
orthogroups_file <- "Orthogroups.tsv"
og_list_dir <- "Orthogroups.tsv"
annotation_dir <- "01_annotated_genes"

# === Load OrthoFinder orthogroups file ===
orthogroups <- read_tsv(orthogroups_file)

# Create a long-format OG → gene mapping
og_gene_map <- orthogroups %>%
  pivot_longer(-Orthogroup, names_to = "Species", values_to = "Genes") %>%
  mutate(Genes = strsplit(Genes, ", ")) %>%
  unnest(Genes) %>%
  distinct(Orthogroup, Species, Genes)


# === Process each species ===
ann_files <- list.files(annotation_dir, pattern = "_annotated.tsv", full.names = TRUE)

for (ann_file in ann_files) {
  species <- str_replace(basename(ann_file), "_annotated.tsv", "")
  message("Processing species: ", species)
  
  # Load OG list
  og_ids <- unique(og_gene_map$Orthogroup)
  og_subset <- og_gene_map %>%
    filter(`Orthogroup` %in% og_ids & str_detect(Genes, species))  # Filter OG + species
  
  
  # Load merged annotation
  annotation_file <- file.path(annotation_dir, paste0(species, "_annotated.tsv"))
  ann <- read_tsv(annotation_file, show_col_types = FALSE)
  
  # Clean GeneID column name
  colnames(ann)[1] <- "GeneID"
  
  
  # Merge with OG data
  annotated <- og_subset %>%
    rename(GeneID = Genes) %>%
    left_join(ann, by = "GeneID") %>%
    select(
      Orthogroup, GeneID,
      EggNOG_PFAMs_domains,
      EggNOG_PFAMs_descriptions,
      EggNOG_GO_terms,
      EggNOG_KEGG_KO,
      EggNOG_KEGG_paths,
      EggNOG_BRITE_class,
      Interpro_Pfam_Ids,
      Interpro_Pfam_descriptions,
      Interpro_Pfam_Evals,
      InterPro_IDs,
      InterPro_descriptions,
      Interpro_GO_terms,
      Panther_ID, Panther_descriptions, Panther_Evals
    )
  
  # Output file
  out_file <- file.path("02_annotated_OGs", paste0(species, "_annotated_OGs.tsv"))
  write_tsv(annotated, out_file)
  message("[✓] Output written: ", out_file)
}
```
## split the terminal significant gene families using this python script:
```
#!/usr/bin/env python3

from collections import defaultdict
import os

# === Input settings ===
input_file = "Gamma_branch_probabilities.tab"
output_dir = "Species_specific_sig_OGs"
threshold = 0.05

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Dictionary to hold column name -> list of matching FamilyIDs
results = defaultdict(list)

# Read the input file
with open(input_file) as f:
    header = f.readline().strip().split("\t")  # first line = species
    for line in f:
        parts = line.strip().split("\t")
        family_id = parts[0]
        for i in range(1, len(parts)):
            try:
                val = float(parts[i])
                if val <= threshold:
                    results[header[i]].append(family_id)
            except ValueError:
                continue  # skip non-numeric or missing values

# Write each species' significant OGs to its own file
for species, og_list in results.items():
    if og_list:  # only if there are significant OGs
        file_path = os.path.join(output_dir, f"{species}_significant_OGs.txt")
        with open(file_path, "w") as f_out:
            f_out.write("\n".join(og_list))
        print(f"[✓] {species}: {len(og_list)} significant OGs written to {file_path}")
```
## Split these significant gene families into expanding and contracting:

```
# === Libraries ===
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# === Input files ===
gamma_file <- "Gamma_change.tab"
og_dir <- "01_species_sig_ogs"  # folder containing *_sig_ogs.txt files

# === Read CAFE Gamma_change.tab ===
gamma <- read_tsv(gamma_file)

# === Clean column names ===
colnames(gamma)[1] <- "Orthogroup"

# Clean species names
colnames(gamma) <- str_remove(colnames(gamma), "<.*?>")
colnames(gamma) <- str_replace(colnames(gamma),
                               "05_sericata_protein", "05_sericata")
colnames(gamma) <- str_replace(colnames(gamma),
                               "06_megacephala_protein", "06_megacephala")
colnames(gamma) <- str_replace(colnames(gamma),
                               "01_hilli_protein", "01_hilli")
colnames(gamma) <- str_replace(colnames(gamma),
                               "03_stygia_protein", "03_stygia")
colnames(gamma) <- str_replace(colnames(gamma),
                               "02_quadrimaculata_protein", "02_quadrimaculata")
colnames(gamma) <- str_replace(colnames(gamma),
                               "04_vicina_longest_protein", "04_vicina")

# Remove the model <?> column names
gamma <- gamma[, !colnames(gamma) %in% c("", ".1", ".2", ".3")]


# === Convert to long format ===
gamma_long <- gamma %>%
  pivot_longer(-Orthogroup, names_to = "Species", values_to = "Change") %>%
  mutate(
    Change = as.numeric(Change),
    Species = str_replace(Species, "_protein$", "")
  ) %>%
  filter(Change != 0) %>%
  mutate(Direction = ifelse(Change > 0, "expanded", "contracted"))

# === List species-specific OG files ===
og_files <- list.files(og_dir, pattern = "_sig_ogs.txt", full.names = TRUE)

# === Process each species ===
for (og_file in og_files) {
  
  # Extract species name from file
  species <- str_replace(basename(og_file), "_sig_ogs.txt", "")
  message("Processing species: ", species)
  
  # Read significant OGs
  sig_ogs <- read_lines(og_file)
  
  # Match with gamma data
  changes <- gamma_long %>%
    filter(Species == species, Orthogroup %in% sig_ogs)
  
  # Split into expanded and contracted
  expanded <- changes %>% filter(Direction == "expanded") %>% pull(Orthogroup)
  contracted <- changes %>% filter(Direction == "contracted") %>% pull(Orthogroup)
  
  # Write to files
  write_lines(expanded, file.path(og_dir, paste0(species, "_expanded_OGs.txt")))
  write_lines(contracted, file.path(og_dir, paste0(species, "_contracted_OGs.txt")))
  
  message("[✓] ", species, ": ", length(expanded), " expanded, ", length(contracted), " contracted OGs")
}
```
## Annotate the expanding and contracting gene families based on the Eggnog and Interpro annotations
```
#ANNOTATE the expanded and contracted gene families

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# === File paths ===
orthogroups_file <- "Orthogroups.tsv"
og_list_dir <- "expanded_OGs"
annotation_dir <- "../01_merge_annotations/02_annotated_OGs"
output_dir <- "/nesi/nobackup/uow03920/06_blowfly_assembly_jan/analysis/05_fetch_species_specific_gene_IDs_and_annotate_GO_terms/03_expanded_contracted_annotated"  

# === Load OrthoFinder orthogroups file ===
orthogroups <- read_tsv(orthogroups_file)

# Create a long-format OG → gene mapping
og_gene_map <- orthogroups %>%
  pivot_longer(-`Orthogroup`, names_to = "Species", values_to = "Genes") %>%
  mutate(Genes = strsplit(Genes, ", ")) %>%
  unnest(Genes)

og_gene_map <- og_gene_map %>%
  mutate(Species = str_replace(Species, "^[0-9]+_", "")) %>%  # remove 01_, 02_, etc
  mutate(Species = str_replace(Species, "_no_dup_longest_isoforms", ""))  # remove suffix


# === Process each species ===
og_files <- list.files(og_list_dir, pattern = "_expanded_OGs.txt", full.names = TRUE)

for (og_file in og_files) {
  species <- str_replace(basename(og_file), "_expanded_OGs.txt", "")
  message("Processing species: ", species)
  
  # Load OG list
  og_ids <- read_lines(og_file)
  og_subset <- og_gene_map %>%
    filter(`Orthogroup` %in% og_ids & str_detect(Genes, species))  # Filter OG + species
  
  
  # Load merged annotation
  annotation_file <- file.path(annotation_dir, paste0(species, "_annotated.tsv"))
  ann <- read_tsv(annotation_file, show_col_types = FALSE)
  
  # Clean GeneID column name
  colnames(ann)[1] <- "GeneID"
  
  
  # Merge with OG data
  annotated <- og_subset %>%
    rename(GeneID = Genes) %>%
    left_join(ann, by = "GeneID") %>%
    select(
      Orthogroup, GeneID,
      EggNOG_PFAMs_domains,
      EggNOG_GO_terms,
      EggNOG_KEGG_KO,
      EggNOG_KEGG_paths,
      EggNOG_BRITE_class,
      Interpro_Pfam_Ids,
      Interpro_Pfam_descriptions,
      InterPro_IDs,
      InterPro_descriptions,
      Interpro_GO_terms
    )
  
  # Output file
  out_file <- file.path(output_dir, paste0(species, "_annotated_expanded_OGs.tsv"))
  write_tsv(annotated, out_file)
}


##################################### Contracted ##################################################################
# === File paths ===
orthogroups_file <- "Orthogroups.tsv"
og_list_dir <- "contracted_OGs"
annotation_dir <- "../01_merge_annotations/02_annotated_OGs"
output_dir <- "/nesi/nobackup/uow03920/06_blowfly_assembly_jan/analysis/05_fetch_species_specific_gene_IDs_and_annotate_GO_terms/03_expanded_contracted_annotated"  

# === Load OrthoFinder orthogroups file ===
orthogroups <- read_tsv(orthogroups_file)

# Create a long-format OG → gene mapping
og_gene_map <- orthogroups %>%
  pivot_longer(-`Orthogroup`, names_to = "Species", values_to = "Genes") %>%
  mutate(Genes = strsplit(Genes, ", ")) %>%
  unnest(Genes)

og_gene_map <- og_gene_map %>%
  mutate(Species = str_replace(Species, "^[0-9]+_", "")) %>%  # remove 01_, 02_, etc
  mutate(Species = str_replace(Species, "_no_dup_longest_isoforms", ""))  # remove suffix


# === Process each species ===
og_files <- list.files(og_list_dir, pattern = "_contracted_OGs.txt", full.names = TRUE)

for (og_file in og_files) {
  species <- str_replace(basename(og_file), "_contracted_OGs.txt", "")
  message("Processing species: ", species)
  
  # Load OG list
  og_ids <- read_lines(og_file)
  og_subset <- og_gene_map %>%
    filter(`Orthogroup` %in% og_ids & str_detect(Genes, species))  # Filter OG + species
  
  
  # Load merged annotation
  annotation_file <- file.path(annotation_dir, paste0(species, "_annotated.tsv"))
  ann <- read_tsv(annotation_file, show_col_types = FALSE)
  
  # Clean GeneID column name
  colnames(ann)[1] <- "GeneID"
  
  
  # Merge with OG data
  annotated <- og_subset %>%
    rename(GeneID = Genes) %>%
    left_join(ann, by = "GeneID") %>%
    select(
      Orthogroup, GeneID,
      EggNOG_PFAMs_domains,
      EggNOG_GO_terms,
      EggNOG_KEGG_KO,
      EggNOG_KEGG_paths,
      EggNOG_BRITE_class,
      Interpro_Pfam_Ids,
      Interpro_Pfam_descriptions,
      InterPro_IDs,
      InterPro_descriptions,
      Interpro_GO_terms
    )
  
  # Output file
  out_file <- file.path(output_dir, paste0(species, "_annotated_contracted_OGs.tsv"))
  write_tsv(annotated, out_file)
}
```

I then manually went through these tsv files and counted how many gene families were associated with the PFAMs that are associated with transposable elements. I created a heat map of these using this script: 

```
library(tidyverse)

df <- tribble(
  ~Pfam,     ~HIL, ~QUA, ~STY, ~VIC, ~SER, ~MEG,
  "PF00075", 1,    0,    0,    0,    5,    7,
  "PF00078", 2,    0,    2,    2,    89,   36,
  "PF00665", 0,    1,    1,    0,    36,   14,
  "PF02925", 0,    0,    0,    0,    0,    0,
  "PF02992", 0,    0,    0,    0,    0,    0,
  "PF03184", 1,    0,    1,    0,    1,    0,
  "PF03221", 1,    0,    0,    0,    0,    0,
  "PF03732", 0,    0,    1,    0,    11,   5,
  "PF04687", 0,    0,    0,    0,    0,    0,
  "PF05699", 1,    0,    1,    0,    0,    1,
  "PF05840", 0,    0,    0,    0,    0,    0,
  "PF05970", 0,    0,    2,    0,    2,    3,
  "PF07727", 0,    0,    1,    0,    7,    7,
  "PF08283", 0,    0,    0,    0,    0,    0,
  "PF08284", 0,    0,    0,    0,    1,    0,
  "PF10551", 0,    0,    0,    0,    0,    0,
  "PF13358", 2,    0,    2,    0,    7,    4,
  "PF13359", 5,    0,    0,    0,    2,    0,
  "PF13456", 0,    0,    0,    0,    0,    0,
  "PF13837", 1,    0,    0,    0,    0,    0,
  "PF13976", 0,    0,    0,    0,    5,    6,
  "PF14214", 0,    0,    1,    0,    4,    3,
  "PF14223", 2,    0,    3,    0,    10,   7,
  "PF14529", 0,    0,    0,    1,    20,   13
)

# Remove rows with all zeros
df_filtered <- df %>%
  filter(rowSums(across(-Pfam)) > 0)

# Convert to long format
df_long <- df_filtered %>%
  pivot_longer(
    cols = -Pfam,
    names_to = "Species",
    values_to = "Count"
  )

df_long$Species <- factor(
  df_long$Species,
  levels = c("QUA", "STY", "HIL", "VIC", "SER", "MEG")
)

# Heatmap
ggplot(df_long, aes(x = Species, y = Pfam, fill = Count)) +
  geom_tile(color = "white") +
  
  scale_fill_gradientn(
    colours = c("#FFC7DE", "#EB71B9", "#3C0D70"),
    trans = "log1p",
    breaks = c(0, 10, 40, 80)
  ) +
  labs(
    x = NULL) +
  
  theme_minimal() +
  theme(
    text = element_text(size = 16, colour = "black"),
    axis.text.x = element_text(size = 16, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    panel.grid = element_blank()
  )

ggsave(
  "TE_gene_family_heatmap.png",
  dpi = 300
)
```














 



































