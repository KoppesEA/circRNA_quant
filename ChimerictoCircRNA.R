##Erik Koppes 06/03/2019
##Script to Identify Chimeric Reads

library(readr)
library(dplyr)
library(tidyr)
library(stringr)

##Step1 import STAR .out.junction files
ChimericColnames = c("chr_donor", "intronbase1_donor", "strand_donor",
                     "chr_acceptor", "intronbase1_acceptor", "strand_acceptor",
                     "junction_type", "repeatlength_left", "repeatlength_right",
                     "read_name", "base1_seg1", "CIGAR_seg1",
                     "base1_seg2", "CIGAR_seg2")
Con1 <- read_tsv("Control_1_S1Chimeric.out.junction", col_names = ChimericColnames) %>%
  filter(chr_donor == chr_acceptor) %>%
  filter(chr_donor != "MT") %>%
  filter(junction_type >= 0) %>%
  filter(strand_donor == strand_acceptor) %>%
  filter(abs(intronbase1_donor - intronbase1_acceptor) < 100000) %>%
  mutate(back_splice = (if_else(strand_donor == "+" & intronbase1_donor > intronbase1_acceptor, "YES",
                      if_else(strand_donor == "-" & intronbase1_donor < intronbase1_acceptor, "YES", "NO"))))
Con2 <- read_tsv("Control_2_S2Chimeric.out.junction", col_names = ChimericColnames)%>%
  filter(chr_donor == chr_acceptor) %>%
  filter(chr_donor != "MT") %>%
  filter(junction_type >= 0) %>%
  filter(strand_donor == strand_acceptor) %>%
  filter(abs(intronbase1_donor - intronbase1_acceptor) < 100000) %>%
  mutate(back_splice = (if_else(strand_donor == "+" & intronbase1_donor > intronbase1_acceptor, "YES",
                                if_else(strand_donor == "-" & intronbase1_donor < intronbase1_acceptor, "YES", "NO"))))
Con3 <- read_tsv("Control_3_S3Chimeric.out.junction", col_names = ChimericColnames)%>%
  filter(chr_donor == chr_acceptor) %>%
  filter(chr_donor != "MT") %>%
  filter(junction_type >= 0) %>%
  filter(strand_donor == strand_acceptor) %>%
  filter(abs(intronbase1_donor - intronbase1_acceptor) < 100000) %>%
  mutate(back_splice = (if_else(strand_donor == "+" & intronbase1_donor > intronbase1_acceptor, "YES",
                                if_else(strand_donor == "-" & intronbase1_donor < intronbase1_acceptor, "YES", "NO"))))


Control_cricRNA <- bind_rows(Con1, Con2, Con3) %>% filter(back_splice == "YES") %>%
  arrange(chr_donor, intronbase1_donor)

Control_BED <- Control_cricRNA %>% select(chr_donor, intronbase1_acceptor, intronbase1_donor, strand_donor)%>%
  mutate(score = 1) %>%
  mutate(intronbase1_acceptor - 1) %>%
  select(chr_donor, intronbase1_acceptor, intronbase1_donor, score, strand_donor)

  write_tsv(Control_BED, "Control_BED_Circ.BED", col_names = F)

circRNA_con_annot <- read_tsv("Rattus_norvegicus.Rnor_circ_annot.bed",
                              col_names = Annot_cols) %>%
#  str_replace_all(strand_A_chrom_B, "+", "plus")
  separate(strand_A_chrom_B, into = c("strand_A", "chrom_B"), sep = "\\|", convert = T) %>% ##doesn't seem to like + and - as non numeric
  filter(chrom_B != "NA") %>%
#  filter(!(strand_A_chrom_B %in% c("-|NA", "+|NA"))) ## removes 1176 no features from 49946 = 48770
  filter(strand_A == strand_B) ## this filter removes 46134 leaving only 2635
length(unique(circRNA_con_annot$gene_ID_B)) ##8589 genes with circRNA [w/o strand filter] or 1304 with strand filter

circRNA_con_annot_exons <- circRNA_con_annot %>% filter(feature_B == "exon") %>%
  separate(attributes_B, sep = ";", into = c("a","b","c","d","e","f","g","h", "i", "j", "k", "l", "m"))
  
Annot_cols = c("chrom_A", "chromStart_A", "chromEnd_A", "score_A", "strand_A_chrom_B",
               "chromStart_B", "chromEnd_B", "gene_ID_B", "score_B", "strand_B",
               "source_B", "feature_B", "frame_B", "attributes_B")



BedColNames_UCSC <- c("chrom", "chromStart", "chromEnd", "score", "strand", "thickStart",
                 "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
https://genome.ucsc.edu/FAQ/FAQformat.html

chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature, however, the number in position format will be represented. For example, the first 100 bases of chromosome 1 are defined as chrom=1, chromStart=0, chromEnd=100, and span the bases numbered 0-99 in our software (not 0-100), but will represent the position notation chr1:1-100. Read more here.
chromStart and chromEnd can be identical, creating a feature of length 0, commonly used for insertions. For example, use chromStart=0, chromEnd=0 to represent an insertion before the first nucleotide of a chromosome.
The 9 additional optional BED fields are:
  
  name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
shade	 	 	 	 	 	 	 	 	 
score in range  	≤ 166	167-277	278-388	389-499	500-611	612-722	723-833	834-944	≥ 945
strand - Defines the strand. Either "." (=no strand) or "+" or "-".
thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
thickEnd - The ending position at which the feature is drawn thickly (for example the stop codon in gene displays).
itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
blockCount - The number of blocks (exons) in the BED line.
blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
blockStarts - 

## need to convert to add identifier for each back-splice then export a bed file sorted
## then use either bedtool or bedops to add closest gene or exon feature; then recombine for exon level

###
ChimericList = c("Con1", "Con2", "Con3")

for (i in ChimericList) {
  print(paste(i))
}

for (i in ChimericList) {
  base <- (basename(i))
  print(dirname(i))
}
ChimericList
list.files()
