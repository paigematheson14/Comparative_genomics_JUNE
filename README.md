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

This generates a .divsum file which can be used in createRepeatLandscape.pl, also implemented in RepeatMaskers github page (https://github.com/Dfam-consortium/RepeatMasker/blob/master/util/createRepeatLandscape.pl). You need the genome size for this:

```perl createRepeatLandscape.pl -div 01_hilli.divsum -g 565000000 > 01_hilli_repeat_landscape.html```

I wanted to customise the repeat landscapes myself so ran the following R code, but made sure that it corroborated with the above. 



































