# Library Import ####


library(dplyr)
library(tidyr)
library(tibble)
library(clusterProfiler)
library(ggplot2)
library("FactoMineR")
library(factoextra)
library(stringr)
library(plotly)
library(broom)
library(DESeq2)
library(UpSetR)
library(readODS)
library(tibble)
library(openxlsx)
library(gplots)
library(RColorBrewer)

# Colour blind palette ####

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# set working directory ####

setwd("/media/rna/NERPA/Senta")

                     ##################################################################################
                     # graphic representation of raw data, trimmed, mapped and feature counts data ####
                     ##################################################################################


#data frame with general statistics of raw-, trimmed- and mapped data ####

total_sequences_rawdata <- read.delim("/media/rna/NERPA/Senta/Result_tables/total_sequences_rawdata.txt", header = TRUE, sep = "")

### 1. delete last three characters of the values in column Sample of raw data

total_sequences_rawdata$Sample <- str_sub(total_sequences_rawdata$Sample, end=-4)

### 2. reorder the data frame of raw data
total_sequences_rawdata_wide <- total_sequences_rawdata %>% 
  group_by(Sample) %>%
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = Sample, values_from = FastQC_mqc.generalstats.fastqc.total_sequences)

### 3. calculate mean values and standard devisions per column of raw data 

mean_totalseq_rawdata <- as.data.frame(colMeans(total_sequences_rawdata_wide[sapply(total_sequences_rawdata_wide, is.numeric, USE.NAMES = TRUE)]))%>%
  filter(colMeans(total_sequences_rawdata_wide[sapply(total_sequences_rawdata_wide, is.numeric)]) > 5) 

colnames(mean_totalseq_rawdata)[1] <- "mean_value_rawdata"

mean_totalseq_rawdata <- mean_totalseq_rawdata %>%
  rownames_to_column("Sample")


sd_totalseq_rawdata <- as.data.frame(sapply(total_sequences_rawdata_wide, sd)) 

colnames(sd_totalseq_rawdata)[1] <- "sd_value_rawdata"

sd_totalseq_rawdata <- sd_totalseq_rawdata %>%
  rownames_to_column("Sample") %>%
  filter(sd_value_rawdata > 3)

### 4. combine mean- and sd data frame 

mean_sd_totalseq_rawdata<- left_join(mean_totalseq_rawdata, sd_totalseq_rawdata, by="Sample")



total_sequences_trimmed <- read.delim("/media/rna/NERPA/Senta/Result_tables/total_sequences_trimmed.txt", header = TRUE, sep = "")


### 1. delete last three characters of the values in column Sample of trimmed data

total_sequences_trimmed$Sample <- str_sub(total_sequences_trimmed$Sample, end=-18)

total_sequences_trimmed$Sample <- gsub("0", "", total_sequences_trimmed$Sample)

### 2. reorder the data frame of raw data
total_sequences_trimmed_wide <- total_sequences_trimmed %>% 
  group_by(Sample) %>%
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = Sample, values_from = FastQC_mqc.generalstats.fastqc.total_sequences)

### 3. calculate mean values and standard devisions per column of trimmed data 

mean_totalseq_trimmed <- as.data.frame(colMeans(total_sequences_trimmed_wide[sapply(total_sequences_trimmed_wide, is.numeric, USE.NAMES = TRUE)]))%>%
  filter(colMeans(total_sequences_trimmed_wide[sapply(total_sequences_trimmed_wide, is.numeric)]) > 5) 

colnames(mean_totalseq_trimmed)[1] <- "mean_value_trimmed"

mean_totalseq_trimmed <- mean_totalseq_trimmed %>%
  rownames_to_column("Sample")


sd_totalseq_trimmed <- as.data.frame(sapply(total_sequences_trimmed_wide, sd)) 

colnames(sd_totalseq_trimmed)[1] <- "sd_value_trimmed"

sd_totalseq_trimmed <- sd_totalseq_trimmed %>%
  rownames_to_column("Sample") %>%
  filter(sd_value_trimmed > 3)

### 4. combine mean- and sd data frame 

mean_sd_totalseq_trimmed<- left_join(mean_totalseq_trimmed, sd_totalseq_trimmed, by="Sample")



total_sequences_mapped <- read.delim("/media/rna/NERPA/Senta/Result_tables/mapped_sequences.txt", header = TRUE, sep = "")

### 1. delete last three characters and 0 of the values in column Sample of mapped data

total_sequences_mapped$Sample <- str_sub(total_sequences_mapped$Sample, end=-14)

total_sequences_mapped$Sample <- gsub("0", "",total_sequences_mapped$Sample)

### 2. reorder the data frame of mapped data
total_sequences_mapped_wide <- total_sequences_mapped %>% 
  group_by(Sample) %>%
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = Sample, values_from = STAR_mqc.generalstats.star.uniquely_mapped)

### 3. calculate mean values and standard devisions per column of mapped data 

mean_totalseq_mapped <- as.data.frame(colMeans(total_sequences_mapped_wide[sapply(total_sequences_mapped_wide, is.numeric, USE.NAMES = TRUE)]))%>%
  filter(colMeans(total_sequences_mapped_wide[sapply(total_sequences_mapped_wide, is.numeric)]) > 5) 

colnames(mean_totalseq_mapped)[1] <- "mean_value_mapped"

mean_totalseq_mapped <- mean_totalseq_mapped %>%
  rownames_to_column("Sample")


sd_totalseq_mapped <- as.data.frame(sapply(total_sequences_mapped_wide, sd)) 

colnames(sd_totalseq_mapped)[1] <- "sd_value_mapped"

sd_totalseq_mapped <- sd_totalseq_mapped %>%
  rownames_to_column("Sample") %>%
  filter(sd_value_mapped > 3)

### 4. combine mean- and sd data frame 

mean_sd_totalseq_mapped <- left_join(mean_totalseq_mapped, sd_totalseq_mapped, by="Sample")

#combine all three data frames of raw-, mapped and trimmed data 

total_sequences_alldata <- mean_sd_totalseq_rawdata%>%
  left_join(mean_sd_totalseq_trimmed, by="Sample")%>%
  left_join(mean_sd_totalseq_mapped, by = "Sample")


#boxplot with raw-, trimmed, mapped and featureCounts data ####

### 1. formate the tables and join them 

total_sequences_rawdata <- read.delim("/media/rna/NERPA/Senta/tables/total_sequences_rawdata.txt", header = TRUE, sep = "") 
colnames(total_sequences_rawdata)[1] <- "sequence_count_rawdata"
total_sequences_rawdata<- total_sequences_rawdata %>% 
  filter(!grepl('_2', Sample))
total_sequences_rawdata$Sample <- str_sub(total_sequences_rawdata$Sample, end=-3)

total_sequences_trimmed <- read.delim("/media/rna/NERPA/Senta/tables/total_sequences_trimmed.txt", header = TRUE, sep = "")
colnames(total_sequences_trimmed)[1] <- "sequence_count_trimmed"
total_sequences_trimmed<- total_sequences_trimmed %>% 
  filter(!grepl('_2', Sample))
total_sequences_trimmed$Sample <- str_sub(total_sequences_trimmed$Sample, end=-17)
total_sequences_trimmed$Sample <- gsub("0", "", total_sequences_trimmed$Sample)

total_sequences_mapped <- read.delim("/media/rna/NERPA/Senta/tables/mapped_sequences.txt", header = TRUE, sep = "")
total_sequences_mapped$Sample <- str_sub(total_sequences_mapped$Sample, end=-13)
total_sequences_mapped$Sample <- gsub("0", "", total_sequences_mapped$Sample)
colnames(total_sequences_mapped)[2] <- "sequence_count_mapped"

total_featureCounts <- read.csv("/media/rna/NERPA/Senta/tables/feature_Counts_summary.csv", header = TRUE, sep = "")

total_sequences_alldata <- total_sequences_rawdata %>%
  left_join(total_sequences_trimmed, by="Sample") %>%
  left_join(total_sequences_mapped, by = "Sample")%>%
  left_join(total_featureCounts, by="Sample")

total_sequences_alldata$Sample <- str_sub(total_sequences_alldata$Sample, end=-2)

### 2. convert to a longer formate 

total_sequences_alldata_long <- pivot_longer(total_sequences_alldata, cols = c(sequence_count_rawdata, sequence_count_trimmed, sequence_count_mapped,feature_Counts_summary))

### 3. Create a factor to set the right order  

total_sequences_alldata_long$name <- gsub("sequence_count_", "", total_sequences_alldata_long$name)
total_sequences_alldata_long$name <- gsub("_summary", "", total_sequences_alldata_long$name)
total_sequences_alldata_long$name <- gsub("_", " ", total_sequences_alldata_long$name)
total_sequences_alldata_long$name <- gsub("C", "c", total_sequences_alldata_long$name)

total_sequences_alldata_long$name <- factor(total_sequences_alldata_long$name,
                                            levels = c('rawdata','trimmed', 'mapped', 'feature counts'),ordered = TRUE)

total_sequences_alldata_long<- total_sequences_alldata_long%>%
  rename("read_counts" = value)


### 4. Expand the colour palettes of brewer and create a boxplot 
colourCount = length(unique(total_sequences_alldata_long$Sample))
getPalette = colorRampPalette(brewer.pal(11, "Spectral"))

ggplot(total_sequences_alldata_long, aes(x=name, y=read_counts, fill=Sample)) + 
  geom_boxplot()+
  facet_grid(~Sample)+
  scale_fill_manual(values =getPalette(colourCount))+
  theme_bw() +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=24),
        strip.text = element_text(size = 24, face = "bold"),
        legend.title = element_text(size=14),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"),
        legend.position = "none",
        axis.title.x=element_blank())


### 5. Create a dataframe with columns, which shows the percentage of loss compared to raw data 

#delete the last character of the values of the column Sample
total_sequences_rawdata$Sample <- str_sub(total_sequences_rawdata$Sample, end=-2)
#convert the data frame to a wide formate
total_sequences_rawdata_wide <- total_sequences_rawdata %>% 
  group_by(Sample) %>%
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = Sample, values_from = sequence_count_rawdata)
#calculate the mean value 
mean_totalseq_rawdata <- as.data.frame(colMeans(total_sequences_rawdata_wide[sapply(total_sequences_rawdata_wide, is.numeric, USE.NAMES = TRUE)]))%>%
  filter(colMeans(total_sequences_rawdata_wide[sapply(total_sequences_rawdata_wide, is.numeric)]) > 5) 
#change column name of mean value column 
colnames(mean_totalseq_rawdata)[1] <- "mean_value_rawdata"

total_sequences_trimmed$Sample <- str_sub(total_sequences_trimmed$Sample, end=-2)

total_sequences_trimmed_wide <- total_sequences_trimmed %>% 
  group_by(Sample) %>%
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = Sample, values_from = sequence_count_trimmed)

mean_totalseq_trimmed <- as.data.frame(colMeans(total_sequences_trimmed_wide[sapply(total_sequences_trimmed_wide, is.numeric, USE.NAMES = TRUE)]))%>%
  filter(colMeans(total_sequences_trimmed_wide[sapply(total_sequences_trimmed_wide, is.numeric)]) > 5) 

colnames(mean_totalseq_trimmed)[1] <- "mean_value_trimmed"

total_sequences_mapped$Sample <- str_sub(total_sequences_mapped$Sample, end=-2)

total_sequences_mapped_wide <- total_sequences_mapped %>% 
  group_by(Sample) %>%
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = Sample, values_from = sequence_count_mapped)

mean_totalseq_mapped <- as.data.frame(colMeans(total_sequences_mapped_wide[sapply(total_sequences_mapped_wide, is.numeric, USE.NAMES = TRUE)]))%>%
  filter(colMeans(total_sequences_mapped_wide[sapply(total_sequences_mapped_wide, is.numeric)]) > 5) 

colnames(mean_totalseq_mapped)[1] <- "mean_value_mapped"

total_featureCounts$Sample <- str_sub(total_featureCounts$Sample, end=-2)

total_featureCounts_wide <- total_featureCounts %>% 
  group_by(Sample) %>%
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = Sample, values_from = feature_Counts_summary)

mean_featureCounts <- as.data.frame(colMeans(total_featureCounts_wide[sapply(total_featureCounts_wide, is.numeric, USE.NAMES = TRUE)]))%>%
  filter(colMeans(total_featureCounts_wide[sapply(total_featureCounts_wide, is.numeric)]) > 5) 

colnames(mean_featureCounts)[1] <- "mean_value_featureCounts"
#join the three tables 
mean_featureCounts <- mean_featureCounts%>%
  rownames_to_column("Sample")

mean_totalseq_rawdata <- mean_totalseq_rawdata%>%
  rownames_to_column("Sample")

mean_totalseq_trimmed <- mean_totalseq_trimmed%>%
  rownames_to_column("Sample")

mean_totalseq_mapped <- mean_totalseq_mapped%>%
  rownames_to_column("Sample")

mean_totalseq_alldata <- mean_totalseq_rawdata%>%
  left_join(mean_totalseq_trimmed, by="Sample")%>%
  left_join(mean_totalseq_mapped, by="Sample")%>%
  left_join(mean_featureCounts, by="Sample")
#calculate the percentage of sequence loss 
percentage_loss_data <- mean_totalseq_alldata%>%
  mutate(percentage_of_loss_raw_trimmed = 100-(mean_value_trimmed/mean_value_rawdata *100),
         percentage_of_loss_trimmed_mapped = 100-(mean_value_mapped/mean_value_trimmed *100),
         percentage_of_loss_mapped_featureCounts = 100-(mean_value_featureCounts/mean_value_mapped *100),
         mean_rawdata = mean2(mean_value_rawdata),
         mean_trimmed=mean2(mean_value_trimmed),
         mean_mapped = mean2(mean_value_mapped),
         mean_featureCounts = mean2(mean_value_featureCounts),
         percentage_of_loss_rawtri_alldata = 100-(mean_trimmed/mean_rawdata *100),
         percentage_of_loss_trimap_alldata = 100-(mean_mapped/mean_trimmed *100),
         percentage_of_loss_mapfC_alldata = 100-(mean_featureCounts/mean_mapped *100))
write.csv(percentage_loss_data, "/media/rna/NERPA/Senta/Result_tables/percentage_sequence_loss.csv")



                     #########################################################################
                     # Expression Analysis heat vs. control conditions over all genotypes ####     
                     #########################################################################

# Prepare feature Counts and meta data table for DESeq2 analysis ####

### 1. Load tables   
meta_data <- read.csv(here('Textdateien/meta_data.txt'), sep = '')
featureCounts <- read.table("/media/rna/NERPA/Senta/featureCounts/featureCounts.txt", header=TRUE, quote="\"")

### 2. Delete unnecessary information from the table
featureCounts2 <- featureCounts%>% select(-2, -3, -4, -5, -6)

### 3. make the Geneid column to a row 
row.names(featureCounts2) = featureCounts2$Geneid

### 4. save the new column names as a value 
class(colnames(featureCounts2)[-1])
new_names<-substr(colnames(featureCounts2)[2:65],32,37) 

### 5. replace the current colum names with the values of "new_names"
colnames(featureCounts2)[2:65] <- new_names

### 6. remove the Geneid colum 
featureCounts2 <- featureCounts2%>% select(-1)

### 7. change the data frame "featureCounts2" to a matrix 
class(featureCounts2)
featureCounts_matrix <- as.matrix(featureCounts2) 
class(featureCounts_matrix)

### 8. factorising the colums of the meta_data table 
meta_data$condition <- factor(meta_data$condition)
meta_data$genotype <- factor(meta_data$genotype)
meta_data$SampleID <- factor(meta_data$SampleID)


# Normalizing with DESeq2 ####

### 1. Creating the DeSeq2 Data Set 
Condition_Genotype<- DESeqDataSetFromMatrix(countData = featureCounts_matrix, 
                                            colData = meta_data, 
                                            design = ~genotype+condition)

### 2. Running DEseq2 (normalizing counts etc) 
Condition_Genotype <- DESeq(Condition_Genotype)
log2Data <- rlog (Condition_Genotype)
Condition_Genotype_norm <- assay(log2Data)

### 3. normalizing counts by hand 
Condition_Genotype_sizeFactor<- estimateSizeFactors(Condition_Genotype)
sizeFactors(Condition_Genotype_sizeFactor)
Condition_genotype_normCounts<-counts(Condition_Genotype_sizeFactor, normalized=TRUE)
write.csv2(Condition_Genotype_normCounts, "/media/rna/NERPA/Senta/normalisierte_Daten/normCounts.csv")


#PCA separated after genotype and condition with FactoMineR ####

log2_df <- as.matrix(assay(log2Data))
log2_df<- as.data.frame(t(log2_df)) 

Sample_ID <- rownames(log2_df)

log2_df <- log2_df %>% mutate(
  Tolerance = case_when(str_detect(Sample_ID, "RB") ~ "sensitive",
                        str_detect(Sample_ID, "CE") ~ "sensitive",
                        str_detect(Sample_ID, "ACO") ~ "sensitive",
                        str_detect(Sample_ID, "RO") ~ "sensitive",
                        str_detect(Sample_ID, "YG") ~ "tolerant",
                        str_detect(Sample_ID, "SP") ~ "tolerant",
                        str_detect(Sample_ID, "ALT") ~ "tolerant",
                        str_detect(Sample_ID, "SO") ~ "tolerant"),
  
  Condition = case_when(str_detect(Sample_ID, "_H") ~ "heat",
                        str_detect(Sample_ID, "_K") ~ "control"),
  
  Genotype = case_when(str_detect(Sample_ID, "0CE") ~ "CE",
                       str_detect(Sample_ID, "0RB") ~ "RB",
                       str_detect(Sample_ID, "0RO") ~ "RO",
                       str_detect(Sample_ID, "0SO") ~ "SO",
                       str_detect(Sample_ID, "0SP") ~ "SP",
                       str_detect(Sample_ID, "0YG") ~ "YG",
                       str_detect(Sample_ID, "ACO") ~ "ACO",
                       str_detect(Sample_ID, "ALT") ~ "ALT"))

condition_genotype_tolerance.pca <- PCA(log2_df[,c(-32920, -32919, -32918)], graph = FALSE)
  fviz_pca_ind(condition_genotype_tolerance.pca, addEllipses = FALSE, label = FALSE) +
  geom_point(aes(shape=log2_df$Condition, color=log2_df$Genotype), size=6) +
  scale_shape_manual(values = c(19,17))+
  scale_color_manual(values= c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  theme_minimal()+
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=24),
        strip.text = element_text(size = 24),
        legend.title = element_text(size=24),
        legend.text = element_text(size = 24)) +
    labs(col="Genotype",shape="Condition")

# PCA tolerance ####

fviz_pca_ind(condition_genotype_tolerance.pca, addEllipses = FALSE, labelsize=6) +
  geom_point(aes(shape=log2_df$Tolerance, color=log2_df$Tolerance), size=6) +
  scale_color_manual(values= c("#D55E00", "#56B4E9")) +
  scale_shape_manual(values = c(19, 17))+
  theme_minimal() +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=24),
        strip.text = element_text(size = 24),
        legend.title = element_text(size=24),
        legend.text = element_text(size = 24))+
    labs(col="Tolerance",shape="Tolerance")


# DEGs in heat vs. control conditions all genotypes ####

res_cond <- results(Condition_Genotype, contrast = c("condition", "heat", "control"), pAdjustMethod = "fdr")

### 1. Filtering up- and downregulated genes under heat conditions
class(res_cond)
res_cond_data <- as.data.frame(res_cond)

res_heat_up <- res_cond_data%>%
  filter(log2FoldChange>1) %>%
  filter(padj<0.05)


res_heat_down <- res_cond_data%>%
  filter(log2FoldChange < -1) %>%
  filter (padj<0.05)

### 2. Grafic representation of up- or downregulated genes 
res_cond_new <-res_cond_data%>%
  mutate(Expression = case_when(log2FoldChange<=-1 & padj<0.05 ~ "Downregulated",
                                log2FoldChange>=1 & padj < 0.05 ~ "Upregulated",
                                TRUE ~ "No differential expression"))
### 3. Compare res_cond_new table with res_heat_down and res_heat_up table 
table(res_cond_new$Expression)

res_heat_down%>%
  nrow()
res_heat_up%>%
  nrow()

### 4. Vulcano Plot 
png(filename = "/media/rna/NERPA/Senta/Figures/Vulcano.png", units = "px", width = 20000, height = 1000)
ggplot(res_cond_new, aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point(aes(color=Expression), size=3/3) +
  theme_bw()+
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log10 (p-Value)")) +
  scale_x_continuous(breaks=c(-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10)) +
  scale_color_manual(values = c("#D55E00","#999999", '#0072B2' ))+
  guides(color = guide_legend(override.aes = list(size = 8)))+
  theme(legend.text = element_text(size=20))+
  theme(legend.title = element_text(size = 24))+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24))
dev.off()


# GO enrichment of up- and downregulated genes all genotypes ####

### 1. formate the tables 
enrichment_downregulated_withoutint<- res_heat_down%>%
  rownames_to_column("Geneid")%>%
  select("Geneid")

enrichment_downregulated_withoutint$Geneid <- gsub(".v6.1", "", enrichment_downregulated_withoutint$Geneid)

enrichment_upregulated_withoutint <- res_heat_up%>%
  rownames_to_column("Geneid")%>%
  select("Geneid")

enrichment_upregulated_withoutint$Geneid <- gsub(".v6.1", "", enrichment_upregulated_withoutint$Geneid)


### 2. create the table of the enriched genes for up- and downregulated genes and GO Terms and calculate the Gene Ratio manual 

GO_enrichment_upregulated <-enricher(enrichment_upregulated_withoutint$Geneid, TERM2GENE = GO_term2gene, TERM2NAME = GO_term2name)
GO_enrichment_upregulated <- as.data.frame(GO_enrichment_upregulated)
GO_enrichment_upregulated <- GO_enrichment_upregulated%>%
  mutate(Expression="Upregulated")%>%
  separate(GeneRatio,into = c("GeneCounts", "number_of_all_genes"))


GO_enrichment_upregulated$GeneCounts <- as.numeric(GO_enrichment_upregulated$GeneCounts)
GO_enrichment_upregulated$number_of_all_genes <- as.numeric(GO_enrichment_upregulated$number_of_all_genes)

GO_enrichment_upregulated<- GO_enrichment_upregulated%>%
  mutate(GeneRatio = GeneCounts/number_of_all_genes)

GO_enrichment_upreg_top20 <- GO_enrichment_upregulated%>%
  arrange(desc(GeneCounts)) %>%
  head(20)



GO_enrichment_downregulated <- enricher(enrichment_downregulated_withoutint$Geneid, TERM2GENE = GO_term2gene, TERM2NAME = GO_term2name)
GO_enrichment_downregulated<- as.data.frame(GO_enrichment_downregulated)
GO_enrichment_downregulated <- GO_enrichment_downregulated %>%
  mutate(Expression="Downregulated") %>%
  separate(GeneRatio,into = c("GeneCounts", "number_of_all_genes"))

GO_enrichment_downregulated$GeneCounts <- as.numeric(GO_enrichment_downregulated$GeneCounts)
GO_enrichment_downregulated$number_of_all_genes <- as.numeric(GO_enrichment_downregulated$number_of_all_genes)

GO_enrichment_downregulated <- GO_enrichment_downregulated %>%
  mutate(GeneRatio = GeneCounts/number_of_all_genes)

GO_enrichment_downreg_top20 <- GO_enrichment_downregulated%>%
  arrange(desc(GeneCounts)) %>%
  head(20)

### 3. combine up- and downregulated, enriched Gene table and create a barplot with ggplot 

GO_enrichment_upreg_downreg<- bind_rows(GO_enrichment_downregulated, GO_enrichment_upregulated)

GO_downreg_upreg_top20 <- bind_rows(GO_enrichment_downreg_top20, GO_enrichment_upreg_top20)


ggplot(GO_downreg_upreg_top20, aes(x=GeneRatio, y=Description, fill =p.adjust)) +
  geom_bar(stat="identity")+
  facet_grid(~Expression)+
  scale_fill_gradient2(high= "#0072B2",low = "#D55E00", midpoint = 0.025, mid = '#F0E442') +
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text.y = element_text(size=18, face = "bold"),
        axis.text.x = element_text(size=18, face = "bold",vjust = 0.35, angle = 45),
        strip.text = element_text(size = 20, face = "bold"),
        legend.title = element_text(size=18),
        legend.text = element_text(size = 18),
        axis.title.x = element_text(vjust = -0.25)) + 
  geom_text(aes(label=GeneCounts),hjust=2, size=6)


                    #########################################################################
                    # Expression Analysis heat vs. control conditions genotype wise      ####     
                    #########################################################################



# Prepare meta data table and Running DESeq2 ####

meta_data_condgenotype <- read.csv(here('Textdateien/meta_data_condgenotype.txt'), sep = '')

### 1. factorizing the columns of the meta_data table 
meta_data_condgenotype$cond_genotype <- factor(meta_data_condgenotype$cond_genotype)
meta_data_condgenotype$SampleID <- factor(meta_data_condgenotype$SampleID)

### 2. Creating the DeSeq2 Data Set 
Condition_Genotype_DeSeq<- DESeqDataSetFromMatrix(countData = featureCounts_matrix, 
                                                  colData = meta_data_condgenotype, 
                                                  design = ~cond_genotype)
### 3. Normalizing the counts by hand 
Condition_Genotype_DeSeq <- DESeq(Condition_Genotype_DeSeq)
Condition_Genotype_DeSeq <- estimateSizeFactors(Condition_Genotype_DeSeq)
sizeFactors(Condition_Genotype_DeSeq)
Condition_Genotype_normalized<-counts(Condition_Genotype_DeSeq, normalized=TRUE)

# DEGs in heat vs. control conditions genotype wise ####

### result tables Cecile 
res_CE <- results(Condition_Genotype_DeSeq, contrast = c("cond_genotype", "heat_CE", "control_CE"), pAdjustMethod = "fdr")
res_CE_data <- as_data_frame(res_CE)
res_CE_data['Geneid'] <- rownames(res_CE)
### result tables Russet Burbank 
res_RB <- results(Condition_Genotype_DeSeq, contrast = c("cond_genotype", "heat_RB", "control_RB"), pAdjustMethod = "fdr")
res_RB_data <- as_data_frame(res_RB)
res_RB_data['Geneid'] <- rownames(res_RB)
### result tables Rosi 
res_RO <- results(Condition_Genotype_DeSeq, contrast = c("cond_genotype", "heat_RO", "control_RO"), pAdjustMethod = "fdr")
res_RO_data <- as_data_frame(res_RO)
res_RO_data['Geneid'] <- rownames(res_RO)
### result tables Solara 
res_SO <- results(Condition_Genotype_DeSeq, contrast = c("cond_genotype", "heat_SO", "control_SO"), pAdjustMethod = "fdr")
res_SO_data <- as_data_frame(res_SO)
res_SO_data['Geneid'] <- rownames(res_SO)
### result tables Spunta
res_SP <- results(Condition_Genotype_DeSeq, contrast = c("cond_genotype", "heat_SP", "control_SP"), pAdjustMethod = "fdr")
res_SP_data <- as_data_frame(res_SP)
res_SP_data['Geneid'] <- rownames(res_SP)
### result tables Yukon Gold 
res_YG <- results(Condition_Genotype_DeSeq, contrast = c("cond_genotype", "heat_YG", "control_YG"), pAdjustMethod = "fdr")
res_YG_data <- as_data_frame(res_YG)
res_YG_data['Geneid'] <- rownames(res_YG)
### result tables Acoustic 
res_ACO <- results(Condition_Genotype_DeSeq, contrast = c("cond_genotype", "heat_ACO", "control_ACO"), pAdjustMethod = "fdr")
res_ACO_data <- as_data_frame(res_ACO)
res_ACO_data['Geneid'] <- rownames(res_ACO)
### result tables Altus 
res_ALT <- results(Condition_Genotype_DeSeq, contrast = c("cond_genotype", "heat_ALT", "control_ALT"), pAdjustMethod = "fdr")
res_ALT_data <- as_data_frame(res_ALT)
res_ALT_data['Geneid'] <- rownames(res_ALT)

#Significant genes heat vs. control conditions 

#Cecile: significant gene list 
significant_heatgenes_CE<-res_CE_data%>%
  filter(padj<0.05)%>%
  mutate(genotype="1")
#Russet Brubank: significant gene list
significant_heatgenes_RB<-res_RB_data%>%
  filter(padj<0.05)%>%
  mutate(genotype="1")
#Rosi: significant gene list
significant_heatgenes_RO<-res_RO_data%>%
  filter(padj<0.05)%>%
  mutate(genotype="1")
#Solara: significant gene list
significant_heatgenes_SO<-res_SO_data%>%
  filter(padj<0.05)%>%
  mutate(genotype="1")
#Spunta: significant gene list
significant_heatgenes_SP<-res_SP_data%>%
  filter(padj<0.05)%>%
  mutate(genotype="1")
#Yukon Gold: significant gene list
significant_heatgenes_YG<-res_YG_data%>%
  filter(padj<0.05)%>%
  mutate(genotype="1")
#Acoustic: significant gene list
significant_heatgenes_ACO<-res_ACO_data%>%
  filter(padj<0.05)%>%
  mutate(genotype="1")
#Altus: significant gene list
significant_heatgenes_ALT<-res_ALT_data%>%
  filter(padj<0.05)%>%
  mutate(genotype="1")

#DEGs upregulated and downregulated genes 

#Cecile: upregulated and downregulated genes 
Cecile_upregulated<-significant_heatgenes_CE %>%
  filter(log2FoldChange>1)%>%
  mutate(genotype="1")
Cecile_downregulated <- significant_heatgenes_CE%>%
  filter(log2FoldChange < -1)%>%
  mutate(genotype="1")
#Rosi: upregulated and downregulated genes 
RO_upregulated<-significant_heatgenes_RO %>%
  filter(log2FoldChange>1)%>%
  mutate(genotype="1")
RO_downregulated <- significant_heatgenes_RO%>%
  filter(log2FoldChange < -1)%>%
  mutate(genotype="1")
#Russet Burbank: upregulated and downregulated genes 
RB_upregulated<-significant_heatgenes_RB %>%
  filter(log2FoldChange>1)%>%
  mutate(genotype="1")
RB_downregulated <- significant_heatgenes_RB %>%
  filter(log2FoldChange < -1)%>%
  mutate(genotype="1")
#Solara upregulated and downregulated
SO_upregulated<-significant_heatgenes_SO %>%
  filter(log2FoldChange>1)%>%
  mutate(genotype="1")
SO_downregulated <- significant_heatgenes_SO%>%
  filter(log2FoldChange < -1)%>%
  mutate(genotype="1")
#Spunta upregulated and downregulated 
SP_upregulated<-significant_heatgenes_SP %>%
  filter(log2FoldChange>1)%>%
  mutate(genotype="1")
SP_downregulated <- significant_heatgenes_SP%>%
  filter(log2FoldChange < -1)%>%
  mutate(genotype="1")
#Yukon Gold: upregulated and downregulated 
YG_upregulated<-significant_heatgenes_YG %>%
  filter(log2FoldChange>1)%>%
  mutate(genotype="1")
YG_downregulated <- significant_heatgenes_YG%>%
  filter(log2FoldChange < -1)%>%
  mutate(genotype="1")
#Acoustic upregulated and downregulated 
ACO_upregulated<-significant_heatgenes_ACO %>%
  filter(log2FoldChange>1)%>%
  mutate(genotype="1")
ACO_downregulated <- significant_heatgenes_ACO%>%
  filter(log2FoldChange < -1)%>%
  mutate(genotype="1")
#Altus upregulated and downregulated 
ALT_upregulated<-significant_heatgenes_ALT %>%
  filter(log2FoldChange>1)%>%
  mutate(genotype="1")
ALT_downregulated <- significant_heatgenes_ALT%>%
  filter(log2FoldChange < -1)%>%
  mutate(genotype="1")

# Create data frame with counts of upregulated or downregulated genes per genotype ####

genotypes <- c("CE", "ACO", "ALT", "RB", "RO", "SO", "SP", "YG")

upregulated <- c(nrow(Cecile_upregulated), nrow(ACO_upregulated), nrow(ALT_upregulated), nrow(RB_upregulated), nrow(RO_upregulated), nrow(SO_upregulated), nrow(SP_upregulated), nrow(YG_upregulated))

downregulated <- c(nrow(Cecile_downregulated), nrow(ACO_downregulated), nrow(ALT_downregulated), nrow(RB_downregulated), nrow(RO_downregulated), nrow(SO_downregulated), nrow(SP_downregulated), nrow(YG_downregulated))

count_DEGs_allgenotypes <- data.frame(genotypes, upregulated, downregulated)

write.csv(count_DEGs_allgenotypes,"/media/rna/NERPA/Senta/Result_tables/DEG_counts_all_genotypes.csv")

# DEG list genotype wise ####

Cecile_upregulated<- Cecile_upregulated%>%
  select(Geneid, log2FoldChange, padj) 

Cecile_downregulated <- Cecile_downregulated %>%
  select(Geneid, log2FoldChange, padj) 

Cecile_DEGs <- Cecile_upregulated%>%
  full_join(Cecile_downregulated, by= c("Geneid", "log2FoldChange", "padj")) %>%
  mutate(genotype="CE",
         expression = case_when(log2FoldChange<=-1 & padj<0.05 ~ "Downregulated",
                                log2FoldChange>=1 & padj < 0.05 ~ "Upregulated"))


ACO_upregulated<- ACO_upregulated%>%
  select(Geneid, log2FoldChange, padj) 

ACO_downregulated <- ACO_downregulated %>%
  select(Geneid, log2FoldChange, padj) 

ACO_DEGs <- ACO_upregulated%>%
  full_join(ACO_downregulated, by= c("Geneid", "log2FoldChange", "padj")) %>%
  mutate(genotype="ACO",
         expression = case_when(log2FoldChange<=-1 & padj<0.05 ~ "Downregulated",
                                log2FoldChange>=1 & padj < 0.05 ~ "Upregulated"))

ALT_upregulated<- ALT_upregulated%>%
  select(Geneid, log2FoldChange, padj) 

ALT_downregulated <- ALT_downregulated %>%
  select(Geneid, log2FoldChange, padj) 

ALT_DEGs <- ALT_upregulated%>%
  full_join(ALT_downregulated, by= c("Geneid", "log2FoldChange", "padj")) %>%
  mutate(genotype="ALT",
         expression = case_when(log2FoldChange<=-1 & padj<0.05 ~ "Downregulated",
                                log2FoldChange>=1 & padj < 0.05 ~ "Upregulated"))

RB_upregulated<- RB_upregulated%>%
  select(Geneid, log2FoldChange, padj) 

RB_downregulated <- RB_downregulated %>%
  select(Geneid, log2FoldChange, padj) 

RB_DEGs <- RB_upregulated%>%
  full_join(RB_downregulated, by= c("Geneid", "log2FoldChange", "padj")) %>%
  mutate(genotype="RB",
         expression = case_when(log2FoldChange<=-1 & padj<0.05 ~ "Downregulated",
                                log2FoldChange>=1 & padj < 0.05 ~ "Upregulated"))

RO_upregulated<- RO_upregulated%>%
  select(Geneid, log2FoldChange, padj) 

RO_downregulated <- RO_downregulated %>%
  select(Geneid, log2FoldChange, padj) 

RO_DEGs <- RO_upregulated%>%
  full_join(RO_downregulated, by= c("Geneid", "log2FoldChange", "padj")) %>%
  mutate(genotype="RO",
         expression = case_when(log2FoldChange<=-1 & padj<0.05 ~ "Downregulated",
                                log2FoldChange>=1 & padj < 0.05 ~ "Upregulated"))

SO_upregulated<- SO_upregulated%>%
  select(Geneid, log2FoldChange, padj) 

SO_downregulated <- SO_downregulated %>%
  select(Geneid, log2FoldChange, padj) 

SO_DEGs <- SO_upregulated%>%
  full_join(SO_downregulated, by= c("Geneid", "log2FoldChange", "padj")) %>%
  mutate(genotype="SO",
         expression = case_when(log2FoldChange<=-1 & padj<0.05 ~ "Downregulated",
                                log2FoldChange>=1 & padj < 0.05 ~ "Upregulated"))

SP_upregulated<- SP_upregulated%>%
  select(Geneid, log2FoldChange, padj) 

SP_downregulated <- SP_downregulated %>%
  select(Geneid, log2FoldChange, padj) 

SP_DEGs <- SP_upregulated%>%
  full_join(SP_downregulated, by= c("Geneid", "log2FoldChange", "padj")) %>%
  mutate(genotype="SP",
         expression = case_when(log2FoldChange<=-1 & padj<0.05 ~ "Downregulated",
                                log2FoldChange>=1 & padj < 0.05 ~ "Upregulated"))

YG_upregulated<- YG_upregulated%>%
  select(Geneid, log2FoldChange, padj) 

YG_downregulated <- YG_downregulated %>%
  select(Geneid, log2FoldChange, padj) 

YG_DEGs <- YG_upregulated%>%
  full_join(YG_downregulated, by= c("Geneid", "log2FoldChange", "padj")) %>%
  mutate(genotype="YG",
         expression = case_when(log2FoldChange<=-1 & padj<0.05 ~ "Downregulated",
                                log2FoldChange>=1 & padj < 0.05 ~ "Upregulated"))


genotypes_long <- rbind(ACO_DEGs, ALT_DEGs, RB_DEGs, RO_DEGs, SP_DEGs, SO_DEGs, YG_DEGs, Cecile_DEGs)
genotypes_wide <- pivot_wider(data = genotypes_long, names_from = genotype, values_from = c(log2FoldChange, padj, expression))
genotypes_wide <- genotypes_wide[c("Geneid", "log2FoldChange_ACO", "padj_ACO", "expression_ACO",
                                   "log2FoldChange_ALT", "padj_ALT", "expression_ALT",
                                   "log2FoldChange_RB", "padj_RB", "expression_RB",
                                   "log2FoldChange_RO", "padj_RO", "expression_RO",
                                   "log2FoldChange_SP", "padj_SP", "expression_SP",
                                   "log2FoldChange_SO", "padj_SO", "expression_SO",
                                   "log2FoldChange_YG", "padj_YG", "expression_YG",
                                   "log2FoldChange_CE", "padj_CE", "expression_CE")]

# combine with Potato Annotation file 

genotypes_wide$Geneid <- gsub(".v6.1", "", genotypes_wide$Geneid)

genotypes_wide <- genotypes_wide %>%
  rename("Gene" = Geneid)

Potato_Annotation_genotypes <- Potato_Annotation%>%
  select(1,2,4,6,7,8)

genotypes_wide_annotated <- merge(genotypes_wide, Potato_Annotation_genotypes, by="Gene") %>%
  rename("Gen_name" = ...22,
         "pathway" = additional.Pathway.information..2023.,
         "functional_group_v3.1" = MAPMan.functional.group.basieredn.auf.v3.1,
         "functional_category" =Functional.category..compiled.)

write.csv(genotypes_wide_annotated,"/media/rna/NERPA/Senta/Result_tables/Genotypes_DEGs.csv")

# Upset Plot upregulated genes all and most of the genotypes ####
CE_upregulated_genes<-Cecile_upregulated%>%
  select(Geneid, genotype)
RO_upregulated_genes <- RO_upregulated%>%
  select(Geneid, genotype)
RB_upregulated_genes <- RB_upregulated%>%
  select(Geneid, genotype)
SP_upregulated_genes <- SP_upregulated%>%
  select(Geneid, genotype)
SO_upregulated_genes <- SO_upregulated%>%
  select(Geneid, genotype)
ACO_upregulated_genes <- ACO_upregulated%>%
  select(Geneid, genotype)
ALT_upregulated_genes <- ALT_upregulated%>%
  select(Geneid, genotype)
YG_upregulated_genes <- YG_upregulated%>%
  select(Geneid, genotype)

upregulated_genes_all_genotypes<-CE_upregulated_genes%>%
  full_join(RO_upregulated_genes, by="Geneid")%>%
  full_join(RB_upregulated_genes, by="Geneid")%>%
  full_join(SP_upregulated_genes, by="Geneid")%>%
  full_join(SO_upregulated_genes, by="Geneid")%>%
  full_join(ACO_upregulated_genes, by="Geneid")%>%
  full_join(ALT_upregulated_genes, by="Geneid")%>%
  full_join(YG_upregulated_genes, by="Geneid")

upregulated_genes_all_genotypes<-upregulated_genes_all_genotypes%>%
  rename("Cecile"=genotype.x,
         "Rosi"=genotype.y,
         "Russet_Burbank"=genotype.x.x,
         "Spunta"=genotype.y.y,
         "Solara"=genotype.x.x.x,
         "Acoustic"=genotype.y.y.y,
         "Altus"=genotype.x.x.x.x,
         "Yukon_Gold"=genotype.y.y.y.y)

upregulated_genes_all_genotypes_upset <- upregulated_genes_all_genotypes[2:9]
upregulated_genes_all_genotypes_upset <- data.frame(lapply(upregulated_genes_all_genotypes_upset,as.numeric))
upregulated_genes_all_genotypes_upset[is.na(upregulated_genes_all_genotypes_upset)] <- 0

upset(upregulated_genes_all_genotypes_upset, sets = c("Cecile", "Rosi", "Russet_Burbank", "Acoustic", "Spunta", "Solara", "Altus", "Yukon_Gold"),
      sets.bar.color	= c('#D55E00', '#D55E00','#D55E00', '#D55E00', '#0072B2','#0072B2','#0072B2','#0072B2'),
      keep.order = TRUE,
      show.numbers = "yes",
      mb.ratio = c(0.65, 0.35),
      intersections = list(list("Cecile", "Rosi", "Russet_Burbank", "Acoustic", "Spunta", "Solara", "Altus", "Yukon_Gold"),
                           list("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank"),
                           list("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Solara"),
                           list("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Russet_Burbank", "Solara"),
                           list("Cecile", "Acoustic", "Rosi", "Spunta", "Altus", "Russet_Burbank", "Solara"),
                           list("Cecile", "Acoustic", "Rosi", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"),
                           list("Cecile", "Acoustic", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"),
                           list("Cecile", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"),
                           list("Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"),
                           list("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus"),
                           list("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Solara"),
                           list("Cecile", "Acoustic", "Rosi", "Spunta", "Russet_Burbank", "Solara"),
                           list("Cecile", "Acoustic", "Rosi", "Altus", "Russet_Burbank", "Solara"),
                           list("Cecile", "Acoustic", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"),
                           list("Cecile", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"),
                           list("Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara")),
      order.by = "freq", text.scale = 3)

# Get Geneids for all intersection with most or all genotypes upregulated ####

### 1. Gene IDs of different intersections upregulated

upregulated_allgenotypes <- get_intersect_members(Geneid_upregulated_genes, c("Cecile", "Rosi", "Russet_Burbank", "Acoustic", "Spunta", "Solara", "Altus", "Yukon_Gold"))

upregulated_CE_ACO_RO_SP_YG_ALT_RB <- get_intersect_members(Geneid_upregulated_genes, c("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank"))

upregulated_CE_ACO_RO_SP_YG_ALT_SO <- get_intersect_members(Geneid_upregulated_genes, c("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Solara"))

upregulated_CE_ACO_RO_SP_YG_RB_SO <- get_intersect_members(Geneid_upregulated_genes, c("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Russet_Burbank", "Solara"))

upregulated_CE_ACO_RO_SP_ALT_RB_SO <- get_intersect_members(Geneid_upregulated_genes, c("Cecile", "Acoustic", "Rosi", "Spunta", "Altus", "Russet_Burbank", "Solara"))

upregulated_CE_ACO_RO_YG_ALT_RB_SO <- get_intersect_members(Geneid_upregulated_genes, c("Cecile", "Acoustic", "Rosi", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"))

upregulated_CE_ACO_SP_YG_ALT_RB_SO <- get_intersect_members(Geneid_upregulated_genes, c("Cecile", "Acoustic", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"))

upregulated_CE_RO_SP_YG_ALT_RB_SO <- get_intersect_members(Geneid_upregulated_genes, c("Cecile", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"))

upregulated_ACO_RO_SP_YG_ALT_RB_SO <- get_intersect_members(Geneid_upregulated_genes, c("Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"))

upregulated_CE_ACO_RO_SP_YG_ALT <- get_intersect_members(Geneid_upregulated_genes, c("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus"))

upregulated_CE_ACO_RO_SP_YG_SO <- get_intersect_members(Geneid_upregulated_genes, c("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Solara"))

upregulated_CE_ACO_RO_SP_RB_SO <- get_intersect_members(Geneid_upregulated_genes, c("Cecile", "Acoustic", "Rosi", "Spunta", "Russet_Burbank", "Solara"))

upregulated_CE_ACO_RO_ALT_YG_SO <- get_intersect_members(Geneid_upregulated_genes, c("Cecile", "Acoustic", "Rosi", "Altus", "Russet_Burbank", "Solara"))

upregulated_CE_ACO_YG_ALT_RB_SO <- get_intersect_members(Geneid_upregulated_genes, c("Cecile", "Acoustic", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"))

upregulated_CE_SP_YG_ALT_RB_SO <- get_intersect_members(Geneid_upregulated_genes, c("Cecile", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"))

upregulated_RO_SP_YG_ALT_RB_SO <- get_intersect_members(Geneid_upregulated_genes, c("Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"))

upregulated_most_genotypes <- list("CE_RO_RB_ACO_SP_SO_ALT_YG" = upregulated_allgenotypes,
                                   "CE_ACO_RO_SP_YG_ALT_RB" = upregulated_CE_ACO_RO_SP_YG_ALT_RB,
                                   "CE_ACO_RO_SP_YG_ALT_SO" = upregulated_CE_ACO_RO_SP_YG_ALT_SO,
                                   "CE_ACO_RO_SP_YG_RB_SO" = upregulated_CE_ACO_RO_SP_YG_RB_SO,
                                   "CE_ACO_RO_SP_ALT_RB_SO" = upregulated_CE_ACO_RO_SP_ALT_RB_SO,
                                   "CE_ACO_RO_YG_ALT_RB_SO" = upregulated_CE_ACO_RO_YG_ALT_RB_SO,
                                   "CE_ACO_SP_YG_ALT_RB_SO" = upregulated_CE_ACO_SP_YG_ALT_RB_SO,
                                   "CE_RO_SP_YG_ALT_RB_SO" = upregulated_CE_RO_SP_YG_ALT_RB_SO,
                                   "ACO_RO_SP_YG_ALT_RB_SO" = upregulated_ACO_RO_SP_YG_ALT_RB_SO,
                                   "CE_ACO_RO_SP_YG_ALT" = upregulated_CE_ACO_RO_SP_YG_ALT, 
                                   "CE_ACO_RO_SP_YG_SO" = upregulated_CE_ACO_RO_SP_YG_SO,
                                   "CE_ACO_RO_SP_RB_SO" = upregulated_CE_ACO_RO_SP_RB_SO,
                                   "CE_ACO_RO_ALT_YG_SO" = upregulated_CE_ACO_RO_ALT_YG_SO,
                                   "CE_ACO_YG_ALT_RB_SO" = upregulated_CE_ACO_YG_ALT_RB_SO,
                                   "CE_SP_YG_ALT_RB_SO" = upregulated_CE_SP_YG_ALT_RB_SO,
                                   "RO_SP_YG_ALT_RB_SO" = upregulated_RO_SP_YG_ALT_RB_SO)

### 2. combine the intersection Gene IDs with the annotated Gene IDs 
Potato_Annotation <- read_ods("/media/rna/NERPA/Senta/Textdateien/S.tuberosum v6.1 vs. v3.1compiled_Sophia 8.ods", col_names = TRUE)%>%
  select(1,4,21:27)%>%
  rename("Gene"=Feature.ID)

list_upregulated_most_genotypes <- list() #create a list for annotated Gene IDs 
for (intersections in names(upregulated_most_genotypes)) {  #take names of list and save as variable "intersections"
  
  
  
  Gene_df_upregulated_most_genotypes <- as.data.frame(upregulated_most_genotypes[[intersections]]) #save as data frame 
  
  colnames(Gene_df_upregulated_most_genotypes) <- "Gene" 
  
  Gene_df_upregulated_most_genotypes$Gene <- gsub(".v6.1","",Gene_df_upregulated_most_genotypes$Gene) #add .v6.1 to intersections Gene IDs to get a match with Gene IDs in Annotation file 
  
  Gene_df_annotated_upregulated_most_genotypes <- left_join(Gene_df_upregulated_most_genotypes,Potato_Annotation,by="Gene") #join the intersection table with the Annotatione table via Gene IDs for every intersection 
  
  
  list_upregulated_most_genotypes[[intersections]] <- Gene_df_annotated_upregulated_most_genotypes #save the annotated intersection Gene IDs in the list 
}

#Save the annotated Geneids of the intersections in Excel sheets 

### 1. Create a workbook 

wb_upregulated_most_genotypes <- createWorkbook()

### 2. Loop through the list and add each dataframe as a sheet
for (i in seq_along(list_upregulated_most_genotypes)) {
  
  df_upregulated_most_genotypes <- list_upregulated_most_genotypes[[i]]
  sheet_name_upregulated_most_genotypes <- names(list_upregulated_most_genotypes)[i]
  
  addWorksheet(wb_upregulated_most_genotypes, sheetName = sheet_name_upregulated_most_genotypes)
  writeData(wb_upregulated_most_genotypes, sheet = sheet_name_upregulated_most_genotypes, x = df_upregulated_most_genotypes)
}


### 3. Save the workbook to a file
saveWorkbook(wb_upregulated_most_genotypes,
             
             "/media/rna/NERPA/Senta/Result_tables/Annotated_intersections_upregulated_most_genotypes.xlsx",
             overwrite = TRUE)

# Upset Plot downregulated genes all and most genotypes ####

CE_downregulated_genes<-Cecile_downregulated%>%
  select(Geneid, genotype)
RO_downregulated_genes <- RO_downregulated%>%
  select(Geneid, genotype)
RB_downregulated_genes <- RB_downregulated%>%
  select(Geneid, genotype)
SP_downregulated_genes <- SP_downregulated%>%
  select(Geneid, genotype)
SO_downregulated_genes <- SO_downregulated%>%
  select(Geneid, genotype)
ACO_downregulated_genes <- ACO_downregulated%>%
  select(Geneid, genotype)
ALT_downregulated_genes <- ALT_downregulated%>%
  select(Geneid, genotype)
YG_downregulated_genes <- YG_downregulated%>%
  select(Geneid, genotype)

downregulated_genes_all_genotypes<-CE_downregulated_genes%>%
  full_join(RO_downregulated_genes, by="Geneid")%>%
  full_join(RB_downregulated_genes, by="Geneid")%>%
  full_join(SP_downregulated_genes, by="Geneid")%>%
  full_join(SO_downregulated_genes, by="Geneid")%>%
  full_join(ACO_downregulated_genes, by="Geneid")%>%
  full_join(ALT_downregulated_genes, by="Geneid")%>%
  full_join(YG_downregulated_genes, by="Geneid")

downregulated_genes_all_genotypes<-downregulated_genes_all_genotypes%>%
  rename("Cecile"=genotype.x,
         "Rosi"=genotype.y,
         "Russet_Burbank"=genotype.x.x,
         "Spunta"=genotype.y.y,
         "Solara"=genotype.x.x.x,
         "Acoustic"=genotype.y.y.y,
         "Altus"=genotype.x.x.x.x,
         "Yukon_Gold"=genotype.y.y.y.y)

downregulated_genes_all_genotypes_upset <- downregulated_genes_all_genotypes[2:9]
downregulated_genes_all_genotypes_upset <- data.frame(lapply(downregulated_genes_all_genotypes_upset,as.numeric))
downregulated_genes_all_genotypes_upset[is.na(downregulated_genes_all_genotypes_upset)] <- 0

upset(downregulated_genes_all_genotypes_upset, sets = c("Cecile", "Rosi", "Russet_Burbank", "Acoustic", "Spunta", "Solara", "Altus", "Yukon_Gold"),
      sets.bar.color	= c('#D55E00', '#D55E00','#D55E00', '#D55E00', '#0072B2','#0072B2','#0072B2','#0072B2'),
      keep.order = T,
      show.numbers = "yes",
      mb.ratio = c(0.65, 0.35),
      intersections = list(list("Cecile", "Rosi", "Russet_Burbank", "Acoustic", "Spunta", "Solara", "Altus", "Yukon_Gold"),
                           list("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank"),
                           list("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Solara"),
                           list("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Russet_Burbank", "Solara"),
                           list("Cecile", "Acoustic", "Rosi", "Spunta", "Altus", "Russet_Burbank", "Solara"),
                           list("Cecile", "Acoustic", "Rosi", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"),
                           list("Cecile", "Acoustic", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"),
                           list("Cecile", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"),
                           list("Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"),
                           list("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus"),
                           list("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Solara"),
                           list("Cecile", "Acoustic", "Rosi", "Spunta", "Russet_Burbank", "Solara"),
                           list("Cecile", "Acoustic", "Rosi", "Altus", "Russet_Burbank", "Solara"),
                           list("Cecile", "Acoustic", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"),
                           list("Cecile", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"),
                           list("Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara")), 
      order.by = "freq", text.scale = 3)

# Get Geneids for all intersection with most or all genotypes downregulated ####

### 1. Gene IDs of different intersections downregulated

library(micro.gen.extra)

downregulated_allgenotypes <- get_intersect_members(Geneid_downregulated_genes, c("Cecile", "Rosi", "Russet_Burbank", "Acoustic", "Spunta", "Solara", "Altus", "Yukon_Gold"))

downregulated_CE_ACO_RO_SP_YG_ALT_RB <- get_intersect_members(Geneid_downregulated_genes, c("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank"))

downregulated_CE_ACO_RO_SP_YG_ALT_SO <- get_intersect_members(Geneid_downregulated_genes, c("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Solara"))

downregulated_CE_ACO_RO_SP_YG_RB_SO <- get_intersect_members(Geneid_downregulated_genes, c("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Russet_Burbank", "Solara"))

downregulated_CE_ACO_RO_SP_ALT_RB_SO <- get_intersect_members(Geneid_downregulated_genes, c("Cecile", "Acoustic", "Rosi", "Spunta", "Altus", "Russet_Burbank", "Solara"))

downregulated_CE_ACO_RO_YG_ALT_RB_SO <- get_intersect_members(Geneid_downregulated_genes, c("Cecile", "Acoustic", "Rosi", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"))

downregulated_CE_ACO_SP_YG_ALT_RB_SO <- get_intersect_members(Geneid_downregulated_genes, c("Cecile", "Acoustic", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"))

downregulated_CE_RO_SP_YG_ALT_RB_SO <- get_intersect_members(Geneid_downregulated_genes, c("Cecile", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"))

downregulated_ACO_RO_SP_YG_ALT_RB_SO <- get_intersect_members(Geneid_downregulated_genes, c("Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"))

downregulated_CE_ACO_RO_SP_YG_ALT <- get_intersect_members(Geneid_downregulated_genes, c("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Altus"))

downregulated_CE_ACO_RO_SP_YG_SO <- get_intersect_members(Geneid_downregulated_genes, c("Cecile", "Acoustic", "Rosi", "Spunta", "Yukon_Gold", "Solara"))

downregulated_CE_ACO_RO_SP_RB_SO <- get_intersect_members(Geneid_downregulated_genes, c("Cecile", "Acoustic", "Rosi", "Spunta", "Russet_Burbank", "Solara"))

downregulated_CE_ACO_RO_ALT_YG_SO <- get_intersect_members(Geneid_downregulated_genes, c("Cecile", "Acoustic", "Rosi", "Altus", "Russet_Burbank", "Solara"))

downregulated_CE_ACO_YG_ALT_RB_SO <- get_intersect_members(Geneid_downregulated_genes, c("Cecile", "Acoustic", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"))

downregulated_CE_SP_YG_ALT_RB_SO <- get_intersect_members(Geneid_downregulated_genes, c("Cecile", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"))

downregulated_RO_SP_YG_ALT_RB_SO <- get_intersect_members(Geneid_downregulated_genes, c("Rosi", "Spunta", "Yukon_Gold", "Altus", "Russet_Burbank", "Solara"))

downregulated_most_genotypes <- list("CE_RO_RB_ACO_SP_SO_ALT_YG" = downregulated_allgenotypes,
                                     "CE_ACO_RO_SP_YG_ALT_RB" = downregulated_CE_ACO_RO_SP_YG_ALT_RB,
                                     "CE_ACO_RO_SP_YG_ALT_SO" = downregulated_CE_ACO_RO_SP_YG_ALT_SO,
                                     "CE_ACO_RO_SP_YG_RB_SO" = downregulated_CE_ACO_RO_SP_YG_RB_SO,
                                     "CE_ACO_RO_SP_ALT_RB_SO" = downregulated_CE_ACO_RO_SP_ALT_RB_SO,
                                     "CE_ACO_RO_YG_ALT_RB_SO" = downregulated_CE_ACO_RO_YG_ALT_RB_SO,
                                     "CE_ACO_SP_YG_ALT_RB_SO" = downregulated_CE_ACO_SP_YG_ALT_RB_SO,
                                     "CE_RO_SP_YG_ALT_RB_SO" = downregulated_CE_RO_SP_YG_ALT_RB_SO,
                                     "ACO_RO_SP_YG_ALT_RB_SO" = downregulated_ACO_RO_SP_YG_ALT_RB_SO,
                                     "CE_ACO_RO_SP_YG_ALT" = downregulated_CE_ACO_RO_SP_YG_ALT, 
                                     "CE_ACO_RO_SP_YG_SO" = downregulated_CE_ACO_RO_SP_YG_SO,
                                     "CE_ACO_RO_SP_RB_SO" = downregulated_CE_ACO_RO_SP_RB_SO,
                                     "CE_ACO_RO_ALT_YG_SO" = downregulated_CE_ACO_RO_ALT_YG_SO,
                                     "CE_ACO_YG_ALT_RB_SO" = downregulated_CE_ACO_YG_ALT_RB_SO,
                                     "CE_SP_YG_ALT_RB_SO" = downregulated_CE_SP_YG_ALT_RB_SO,
                                     "RO_SP_YG_ALT_RB_SO" = downregulated_RO_SP_YG_ALT_RB_SO)


### 2. combine the intersection Gene IDs with the annotated Gene IDs 
Potato_Annotation <- read_ods("/media/rna/NERPA/Senta/Textdateien/S.tuberosum v6.1 vs. v3.1compiled_Sophia 8.ods", col_names = TRUE)%>%
  select(1,4,21:27)%>%
  rename("Gene"=Feature.ID)

list_downregulated_most_genotypes <- list() #create a list for annotated Gene IDs 
for (intersections in names(downregulated_most_genotypes)) {  #take names of list and save as variable "intersections"
  
  
  
  Gene_df_downregulated_most_genotypes <- as.data.frame(downregulated_most_genotypes[[intersections]]) #save as data frame 
  
  colnames(Gene_df_downregulated_most_genotypes) <- "Gene" 
  
  Gene_df_downregulated_most_genotypes$Gene <- gsub(".v6.1","",Gene_df_downregulated_most_genotypes$Gene) #add .v6.1 to intersections Gene IDs to get a match with Gene IDs in Annotation file 
  
  Gene_df_annotated_downregulated_most_genotypes <- left_join(Gene_df_downregulated_most_genotypes,Potato_Annotation,by="Gene") #join the intersection table with the Annotatione table via Gene IDs for every intersection 
  
  
  list_downregulated_most_genotypes[[intersections]] <- Gene_df_annotated_downregulated_most_genotypes #save the annotated intersection Gene IDs in the list 
}


# Save the annotated Geneids of the intersections in Excel sheets 

### 1. Create a workbook 

wb_downregulated_most_genotypes <- createWorkbook()

# 2. Loop through the list and add each dataframe as a sheet
for (i in seq_along(list_downregulated_most_genotypes)) {
  
  df_downregulated_most_genotypes <- list_downregulated_most_genotypes[[i]]
  sheet_name_downregulated_most_genotypes <- names(list_downregulated_most_genotypes)[i]
  
  addWorksheet(wb_downregulated_most_genotypes, sheetName = sheet_name_downregulated_most_genotypes)
  writeData(wb_downregulated_most_genotypes, sheet = sheet_name_downregulated_most_genotypes, x = df_downregulated_most_genotypes)
}


### 3. Save the workbook to a file
saveWorkbook(wb_downregulated_most_genotypes,
             
             "/media/rna/NERPA/Senta/Result_tables/Annotated_intersections_downregulated_most_genotypes.xlsx",
             overwrite = TRUE)


# boxplot of z-Score of specific genes out of the Upset plot for most genotypes ####

### 1. Creating a table with the GeneIDs per genotype and calculating the zscores

Condition_genotype_normdata <- as.data.frame(Condition_genotype_normCounts)%>%
  rownames_to_column("Geneid")

selected_genes <- data.frame(Geneid= c("Soltu.DM.02G009620.v6.1",
                                       "Soltu.DM.02G014670.v6.1",
                                       "Soltu.DM.05G024030.v6.1",
                                       "Soltu.DM.05G024040.v6.1",
                                       "Soltu.DM.08G003380.v6.1",
                                       "Soltu.DM.08G020430.v6.1",
                                       "Soltu.DM.01G015690.v6.1",
                                       "Soltu.DM.01G008160.v6.1",
                                       "Soltu.DM.01G038100.v6.1",
                                       "Soltu.DM.08G008810.v6.1"))

selected_genes_normCounts <- Condition_genotype_normdata%>%
  filter(Geneid%in%selected_genes$Geneid) %>%
  column_to_rownames("Geneid")%>%
  as.matrix()

zscores_selected_genes <- scale(t(selected_genes_normCounts))

### 2. creating a table for boxplot of zscores

##### 2.1 seperate zscore table after genotype 

zscores_selected_genes<- as.data.frame(t(zscores_selected_genes))

zscores_RB_selected_genes <- zscores_selected_genes%>%
  select(matches("0RB_"))%>%
  t()

zscores_RO_selected_genes <- zscores_selected_genes%>%
  select(matches("0RO_"))%>%
  t()

zscores_SO_selected_genes <- zscores_selected_genes%>%
  select(matches("0SO_"))%>%
  t()

zscores_SP_selected_genes <- zscores_selected_genes%>%
  select(matches("0SP_"))%>%
  t()

zscores_YG_selected_genes <- zscores_selected_genes%>%
  select(matches("YG_"))%>%
  t()

zscores_ACO_selected_genes <- zscores_selected_genes%>%
  select(matches("ACO_"))%>%
  t()

zscores_ALT_selected_genes <- zscores_selected_genes%>%
  select(matches("ALT_"))%>%
  t()
zscores_CE_selected_genes <- zscores_selected_genes%>%
  select(matches("0CE"))%>%
  t()

#### 2.2 creating a loop to add a genotype column and condition column 

list_zscores_per_genotype <- list("zscores_ALT_selected_genes" = zscores_ALT_selected_genes,
                                  "zscores_ACO_selected_genes" = zscores_ACO_selected_genes,
                                  "zscores_CE_selected_genes" = zscores_CE_selected_genes,
                                  "zscores_RB_selected_genes" = zscores_RB_selected_genes,
                                  "zscores_RO_selected_genes" = zscores_RO_selected_genes,
                                  "zscores_SO_selected_genes" = zscores_SO_selected_genes,
                                  "zscores_SP_selected_genes" = zscores_SP_selected_genes,
                                  "zscores_YG_selected_genes" = zscores_YG_selected_genes)

list_zscores_selected_genes_boxplot <- list()

for (i in names(list_zscores_per_genotype)) {
  
  zscores_boxplot_df <- as.data.frame(list_zscores_per_genotype[[i]])
  
  zscores_boxplot_df2 <- zscores_boxplot_df %>%
    mutate(genotype=rownames(zscores_boxplot_df),
           conditions=rownames(zscores_boxplot_df))
  if (startsWith(zscores_boxplot_df2$genotype[1], '0')) {
    zscores_boxplot_df2$genotype <- str_sub(zscores_boxplot_df2$genotype, 2, 3)
  } else {
    zscores_boxplot_df2$genotype <- str_sub(zscores_boxplot_df2$genotype, 1, 3)
  }
  zscores_boxplot_df2$conditions <- str_sub(zscores_boxplot_df2$conditions, -2, -2)
  
  zscores_boxplot_df2$conditions<-  gsub("H", "heat", zscores_boxplot_df2$conditions)
  zscores_boxplot_df2$conditions<-  gsub("K", "control", zscores_boxplot_df2$conditions)
  
  list_zscores_selected_genes_boxplot[[i]] <- zscores_boxplot_df2
  
}

#### 2.3 safe list components as data frame and combine the different tables 

zscores_ALT_boxplot <- as.data.frame(list_zscores_selected_genes_boxplot$zscores_ALT_selected_genes)

zscores_ACO_boxplot <- as.data.frame(list_zscores_selected_genes_boxplot$zscores_ACO_selected_genes)

zscores_CE_boxplot <- as.data.frame(list_zscores_selected_genes_boxplot$zscores_CE_selected_genes)

zscores_RO_boxplot <- as.data.frame(list_zscores_selected_genes_boxplot$zscores_RO_selected_genes)

zscores_RB_boxplot <- as.data.frame(list_zscores_selected_genes_boxplot$zscores_RB_selected_genes)

zscores_SO_boxplot <- as.data.frame(list_zscores_selected_genes_boxplot$zscores_SO_selected_genes)

zscores_SP_boxplot <- as.data.frame(list_zscores_selected_genes_boxplot$zscores_SP_selected_genes)

zscores_YG_boxplot <- as.data.frame(list_zscores_selected_genes_boxplot$zscores_YG_selected_genes)


zscores_boxplot_selected_genes <- rbind(zscores_ACO_boxplot, 
                                        zscores_ALT_boxplot, 
                                        zscores_CE_boxplot, 
                                        zscores_RB_boxplot, 
                                        zscores_RO_boxplot, 
                                        zscores_SO_boxplot, 
                                        zscores_SP_boxplot, 
                                        zscores_YG_boxplot)
rownames(zscores_boxplot_selected_genes) <- NULL

zscores_boxplot_selected_genes_long <- pivot_longer(zscores_boxplot_selected_genes, cols = c(Soltu.DM.01G008160.v6.1, 
                                                                                             Soltu.DM.01G015690.v6.1, 
                                                                                             Soltu.DM.01G038100.v6.1, 
                                                                                             Soltu.DM.02G009620.v6.1, 
                                                                                             Soltu.DM.02G014670.v6.1,
                                                                                             Soltu.DM.05G024030.v6.1,
                                                                                             Soltu.DM.05G024040.v6.1,
                                                                                             Soltu.DM.08G003380.v6.1,
                                                                                             Soltu.DM.08G008810.v6.1,
                                                                                             Soltu.DM.08G020430.v6.1))

zscores_boxplot_selected_genes_long <- zscores_boxplot_selected_genes_long%>%
  rename("ZScore" = value)

zscores_boxplot_selected_genes_long$name <- str_sub(zscores_boxplot_selected_genes_long$name, 1, 18)

ggplot(zscores_boxplot_selected_genes_long, aes(x=genotype, y=ZScore, fill=conditions)) + 
  geom_boxplot()+
  facet_wrap(~factor(name,levels = c("Soltu.DM.02G009620", "Soltu.DM.02G014670", "Soltu.DM.05G024030", "Soltu.DM.05G024040",
                                     "Soltu.DM.08G003380", "Soltu.DM.08G020430", "Soltu.DM.01G015690", "Soltu.DM.01G008160", "Soltu.DM.01G038100", "Soltu.DM.08G008810")), ncol=5)+
  scale_fill_manual(values =c("#0072B2", "#D55E00"))+
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=14),
        strip.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size=14),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"),
        axis.title.x=element_blank())


                    #########################################################################
                    #         Expression Analysis starchgenes genotype wise              ####     
                    #########################################################################

#PCA phenotypical data #### 

### 1. Creating tables with phenotypical data 
starch_content_pca <- read.csv("/media/rna/NERPA/Senta/tables/starch_content_PCA.csv", sep = ',')

sucrose_content_pca <- read.csv("/media/rna/NERPA/Senta/tables/sucrose_content_PCA.csv", sep = '')

soluble_sugar_content_pca <- read.csv("/media/rna/NERPA/Senta/tables/soluble_sugar_content_PCA.csv", sep = '')

phenotype_data <- starch_content_pca %>%
  left_join(sucrose_content_pca, by="name") %>%
  left_join(soluble_sugar_content_pca, by="name") %>%
  select(name, starch_content, sucrose_content, soluble_sugar_content, condition.x, genotype.x, tolerance.x) %>%
  rename("condition" = condition.x,
         "genotype" = genotype.x,
         "tolerance" = tolerance.x) %>%
  column_to_rownames("name")

### 2. PCA phenotype data   
phenotype.pca <- PCA(phenotype_data[,c(-4, -5, -6)], graph = FALSE)

fviz_pca_biplot(phenotype.pca, addEllipses = FALSE, arrowsize=1, labelsize=7, col.var = "#E69F00") +
  geom_point(aes(shape=phenotype_data$condition, color=phenotype_data$tolerance), size=6) +
  scale_color_manual(values= c("#D55E00", "#56B4E9")) +
  scale_shape_manual(values = c(19,17))+
  theme_minimal()+
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=20),
        strip.text = element_text(size = 20),
        legend.title = element_text(size=24),
        legend.text = element_text(size = 20))+
  labs(col="Tolerance",shape="Condition")


# filtering result tables per genotype after starchgenes ####

### 1. join every result table with Starchgene table and add expression column 

res_ACO_data$Geneid <- gsub(".v6.1", "",res_ACO_data$Geneid)

result_genotypes_list <- list(ACO = res_ACO_data,
                              ALT =res_ALT_data, 
                              CE = res_CE_data,
                              RB = res_RB_data,
                              RO = res_RO_data,
                              SO = res_SO_data,
                              SP = res_SP_data,
                              YG = res_YG_data)

starchgenes_genotypes_list <- list()

for (i in names(result_genotypes_list)) {
  
  df_result_genotypes <- as.data.frame(result_genotypes_list[[i]])
  
  
  df_result_genotypes$Geneid <- gsub(".v6.1", "", df_result_genotypes$Geneid)
  
  df_starchgenes_genotypes <- right_join(df_result_genotypes, Starchgenes, by="Geneid")
  
  df_starchgenes_genotypes <- df_starchgenes_genotypes %>%
    mutate(Expression = case_when(log2FoldChange<=-1 & padj<0.05 ~ "Down-regulated",
                                  log2FoldChange>=1 & padj < 0.05 ~ "Up-regulated",
                                  TRUE ~ "No differential expression"))
  
  starchgenes_genotypes_list[[i]] <- df_starchgenes_genotypes
}

### 2. save every genotype dataframe of the starchgenes and reorder the columns 

ALT_starchgenes <- starchgenes_genotypes_list$ALT %>%
  select(Geneid, Enzyme, Expression, log2FoldChange, padj) %>%
  mutate(genotype = "ALT")

ALT_starchgenes <- ALT_starchgenes[,c(1,2,4,5,3,6)]

ACO_starchgenes <- starchgenes_genotypes_list$ACO%>%
  select(Geneid, Enzyme, Expression, log2FoldChange, padj)%>%
  mutate(genotype = "ACO")

ACO_starchgenes <- ACO_starchgenes[,c(1,2,4,5,3,6)]

CE_starchgenes <- starchgenes_genotypes_list$CE %>%
  select(Geneid, Enzyme, Expression, log2FoldChange, padj)%>%
  mutate(genotype = "CE")

CE_starchgenes <- CE_starchgenes[,c(1,2,4,5,3,6)]

RB_starchgenes <- starchgenes_genotypes_list$RB %>%
  select(Geneid, Enzyme, Expression, log2FoldChange, padj)%>%
  mutate(genotype = "RB")

RB_starchgenes <- RB_starchgenes[,c(1,2,4,5,3,6)]

RO_starchgenes <- starchgenes_genotypes_list$RO %>%
  select(Geneid, Enzyme, Expression, log2FoldChange, padj)%>%
  mutate(genotype = "RO")

RO_starchgenes <- RO_starchgenes[,c(1,2,4,5,3,6)]

SO_starchgenes <-starchgenes_genotypes_list$SO %>%
  select(Geneid, Enzyme, Expression, log2FoldChange, padj)%>%
  mutate(genotype = "SO")

SO_starchgenes <- SO_starchgenes[,c(1,2,4,5,3,6)]

SP_starchgenes <-starchgenes_genotypes_list$SP %>%
  select(Geneid, Enzyme, Expression, log2FoldChange, padj)%>%
  mutate(genotype = "SP")

SP_starchgenes <- SP_starchgenes[,c(1,2,4,5,3,6)]

YG_starchgenes <-starchgenes_genotypes_list$YG %>%
  select(Geneid, Enzyme, Expression, log2FoldChange, padj)%>%
  mutate(genotype = "YG")

YG_starchgenes <- YG_starchgenes[,c(1,2,4,5,3,6)]

### 3. merge all tables together 

starchgenes_long <- rbind(ALT_starchgenes, ACO_starchgenes, RB_starchgenes, RO_starchgenes, SP_starchgenes, SO_starchgenes, CE_starchgenes, YG_starchgenes)
starchgenes_wide <- pivot_wider(data = starchgenes_long, names_from = genotype, values_from = c(log2FoldChange, padj, Expression))
starchgenes_wide <- starchgenes_wide[c("Geneid", "Enzyme", "log2FoldChange_ACO", "padj_ACO", "Expression_ACO",
                                       "log2FoldChange_ALT", "padj_ALT", "Expression_ALT",
                                       "log2FoldChange_RB", "padj_RB", "Expression_RB",
                                       "log2FoldChange_RO", "padj_RO", "Expression_RO",
                                       "log2FoldChange_SP", "padj_SP", "Expression_SP",
                                       "log2FoldChange_SO", "padj_SO", "Expression_SO",
                                       "log2FoldChange_YG", "padj_YG", "Expression_YG",
                                       "log2FoldChange_CE", "padj_CE", "Expression_CE")]

write.csv(starchgenes_wide, "/media/rna/NERPA/Senta/Result_tables/Starchgenes_genotypes.csv") 

write.csv(starchgenes_long, "/media/rna/NERPA/Senta/Result_tables/Starchgenes_genotypes_long.csv") 

# search for starchgenes with differential expression ####


starchgenes_long2 <- starchgenes_long%>%
  select(-padj, -log2FoldChange)


starchgenes_expression <- pivot_wider(data=starchgenes_long2, names_from = genotype, values_from = c(Expression), values_fn = list )


### 1. filter genotype starchgenes after DEGs

DEGs_ALT_starchgenes <- ALT_starchgenes%>%
  filter(Expression != 'No differential expression')
write.csv(DEGs_ALT_starchgenes, '/media/rna/NERPA/Senta/Result_tables/DEGs_ALT')

DEGs_ACO_starchgenes <- ACO_starchgenes%>%
  filter(Expression != 'No differential expression')

DEGs_RB_starchgenes <- RB_starchgenes%>%
  filter(Expression != 'No differential expression')

DEGs_RO_starchgenes <- RO_starchgenes%>%
  filter(Expression != 'No differential expression')

DEGs_SO_starchgenes <- SO_starchgenes%>%
  filter(Expression != 'No differential expression')

DEGs_SP_starchgenes <- SP_starchgenes%>%
  filter(Expression != 'No differential expression')

DEGs_CE_starchgenes <- CE_starchgenes%>%
  filter(Expression != 'No differential expression')

DEGs_YG_starchgenes <- YG_starchgenes%>%
  filter(Expression != 'No differential expression')

### 2. safe as excel sheets 

DEGs_starchgenes_allgenotypes <- list(DEGs_ALT = DEGs_ALT_starchgenes,
                                      DEGs_ACO = DEGs_ACO_starchgenes,
                                      DEGs_RB = DEGs_RB_starchgenes,
                                      DEGs_RO = DEGs_RO_starchgenes,
                                      DEGs_SO = DEGs_SO_starchgenes,
                                      DEGs_SP = DEGs_SP_starchgenes,
                                      DEGs_CE = DEGs_CE_starchgenes,
                                      DEGs_YG = DEGs_YG_starchgenes)

wb_DEGs_starchgenes <- createWorkbook()

for (i in seq_along(DEGs_starchgenes_allgenotypes)) {
  
  df_DEGs_starchgenes_allgenotypes <- DEGs_starchgenes_allgenotypes[[i]]
  sheet_name_DEGs_starchgenes_allgenotypes <- names(DEGs_starchgenes_allgenotypes)[i]
  
  addWorksheet(wb_DEGs_starchgenes, sheetName = sheet_name_DEGs_starchgenes_allgenotypes )
  writeData(wb_DEGs_starchgenes, sheet = sheet_name_DEGs_starchgenes_allgenotypes, x = df_DEGs_starchgenes_allgenotypes)
}


saveWorkbook(wb_DEGs_starchgenes,
             
             "/media/rna/NERPA/Senta/Result_tables/DEGs_starchgenes_allgenotypes.xlsx",
             overwrite = TRUE)

#creating a boxplot of starch DEGs ####
### 1. create a table of DEGs of starchgenes over all genotypes and calculate the zscore 

DEGs_starchgenes_all_genotypes<- rbind(DEGs_ACO_starchgenes, DEGs_ALT_starchgenes, DEGs_CE_starchgenes, DEGs_RB_starchgenes, DEGs_RO_starchgenes, DEGs_SP_starchgenes, DEGs_YG_starchgenes, DEGs_SO_starchgenes)

DEGs_starchgenes_all_genotypes <- DEGs_starchgenes_all_genotypes%>%
  select(Geneid)

DEGs_starchgenes_all_genotypes <- unique(DEGs_starchgenes_all_genotypes)

DEGs_starchgenes_all_genotypes$Geneid <- paste0(DEGs_starchgenes_all_genotypes$Geneid, ".v6.1")

norm_counts_starch_DEGs <- DEGs_starchgenes_all_genotypes%>%
  left_join(Condition_genotype_normdata, by="Geneid")

norm_counts_starch_DEGs$Geneid <- gsub(".v6.1", "", norm_counts_starch_DEGs$Geneid)

norm_counts_starch_DEGs <- norm_counts_starch_DEGs%>%
  column_to_rownames("Geneid")%>%
  as.matrix()

zscore_starch_DEGs <- scale(t(norm_counts_starch_DEGs)) 

### 2. creating a table for boxplot of zscores

##### 2.1 seperate zscore table after genotype 

zscore_starch_DEGs<- as.data.frame(t(zscore_starch_DEGs))

zscores_RB_starch_DEGS <- zscore_starch_DEGs%>%
  select(matches("0RB_"))%>%
  t()

zscores_RO_starch_DEGS <- zscore_starch_DEGs%>%
  select(matches("0RO_"))%>%
  t()

zscores_SO_starch_DEGS <- zscore_starch_DEGs%>%
  select(matches("0SO_"))%>%
  t()

zscores_SP_starch_DEGS <- zscore_starch_DEGs%>%
  select(matches("0SP_"))%>%
  t()

zscores_YG_starch_DEGS <- zscore_starch_DEGs%>%
  select(matches("YG_"))%>%
  t()

zscores_ACO_starch_DEGS <- zscore_starch_DEGs%>%
  select(matches("ACO_"))%>%
  t()

zscores_ALT_starch_DEGS <- zscore_starch_DEGs%>%
  select(matches("ALT_"))%>%
  t()

zscores_CE_starch_DEGs <- zscore_starch_DEGs%>%
  select(matches("0CE"))%>%
  t()

#### 2.2 creating a loop to add a genotype column and condition column 

list_zscores_starch_DEGs <- list("zscores_ALT_starch_DEGs" = zscores_ALT_starch_DEGS,
                                  "zscores_ACO_starch_DEGs" = zscores_ACO_starch_DEGS,
                                  "zscores_CE_starch_DEGs" = zscores_CE_starch_DEGs,
                                  "zscores_RB_starch_DEGs" = zscores_RB_starch_DEGS,
                                  "zscores_RO_starch_DEGs" = zscores_RO_starch_DEGS,
                                  "zscores_SO_starch_DEGs" = zscores_SO_starch_DEGS,
                                  "zscores_SP_starch_DEGs" = zscores_SP_starch_DEGS,
                                  "zscores_YG_starch_DEGs" = zscores_YG_starch_DEGS)

list_zscores_starch_DEGs_boxplot <- list()

for (i in names(list_zscores_starch_DEGs)) {
  
  zscores_starch_DEGs_boxplot_df <- as.data.frame(list_zscores_starch_DEGs[[i]])
  
  zscores_starch_DEGs_boxplot_df2 <- zscores_starch_DEGs_boxplot_df %>%
    mutate(genotype=rownames(zscores_starch_DEGs_boxplot_df),
           conditions=rownames(zscores_starch_DEGs_boxplot_df))
  if (startsWith(zscores_starch_DEGs_boxplot_df2 $genotype[1], '0')) {
    zscores_starch_DEGs_boxplot_df2 $genotype <- str_sub(zscores_starch_DEGs_boxplot_df2$genotype, 2, 3)
  } else {
    zscores_starch_DEGs_boxplot_df2$genotype <- str_sub(zscores_starch_DEGs_boxplot_df2$genotype, 1, 3)
  }
  zscores_starch_DEGs_boxplot_df2$conditions <- str_sub(zscores_starch_DEGs_boxplot_df2$conditions, -2, -2)
  
  zscores_starch_DEGs_boxplot_df2$conditions<-  gsub("H", "heat", zscores_starch_DEGs_boxplot_df2 $conditions)
  zscores_starch_DEGs_boxplot_df2$conditions<-  gsub("K", "control", zscores_starch_DEGs_boxplot_df2 $conditions)
  
  list_zscores_starch_DEGs_boxplot[[i]] <- zscores_starch_DEGs_boxplot_df2 
  
}

#### 2.3 safe list components as data frame and combine the different tables 

zscores_ALT_starch_DEGs_boxplot <- as.data.frame(list_zscores_starch_DEGs_boxplot$zscores_ALT_starch_DEGs)

zscores_ACO_starch_DEGs_boxplot <- as.data.frame(list_zscores_starch_DEGs_boxplot$zscores_ACO_starch_DEGs)

zscores_CE_starch_DEGs_boxplot <- as.data.frame(list_zscores_starch_DEGs_boxplot$zscores_CE_starch_DEGs)

zscores_RO_starch_DEGs_boxplot <- as.data.frame(list_zscores_starch_DEGs_boxplot$zscores_RO_starch_DEGs)

zscores_RB_starch_DEGs_boxplot <- as.data.frame(list_zscores_starch_DEGs_boxplot$zscores_RB_starch_DEGs)

zscores_SO_starch_DEGs_boxplot <- as.data.frame(list_zscores_starch_DEGs_boxplot$zscores_SO_starch_DEGs)

zscores_SP_starch_DEGs_boxplot <- as.data.frame(list_zscores_starch_DEGs_boxplot$zscores_SP_starch_DEGs)

zscores_YG_starch_DEGs_boxplot <- as.data.frame(list_zscores_starch_DEGs_boxplot$zscores_YG_starch_DEGs)

zscores_boxplot_starch_DEGs <- rbind(zscores_ALT_starch_DEGs_boxplot, 
                                     zscores_ACO_starch_DEGs_boxplot, 
                                     zscores_CE_starch_DEGs_boxplot, 
                                     zscores_RO_starch_DEGs_boxplot, 
                                     zscores_RB_starch_DEGs_boxplot, 
                                     zscores_SO_starch_DEGs_boxplot, 
                                     zscores_SP_starch_DEGs_boxplot, 
                                     zscores_YG_starch_DEGs_boxplot)
rownames(zscores_boxplot_starch_DEGs) <- NULL

starch_synthesis_DEGs <- read.csv("/media/rna/NERPA/Senta/tables/starch_synthesis_DEGs.csv", sep = ",")

starch_synthesis_DEGs <- starch_synthesis_DEGs%>%
  column_to_rownames("Geneid")

zscore_boxplot_starch_synthesis <- zscores_boxplot_starch_DEGs %>%
  select(rownames(starch_synthesis_DEGs), genotype, conditions) 

zscore_boxplot_starch_synthesis_long <- pivot_longer(zscore_boxplot_starch_synthesis, cols = c(Soltu.DM.08G030230,
                                                                                               Soltu.DM.12G026390,
                                                                                               Soltu.DM.02G020800,
                                                                                               Soltu.DM.03G019120,
                                                                                               Soltu.DM.07G022290,
                                                                                               Soltu.DM.09G004100,
                                                                                               Soltu.DM.07G013360,
                                                                                               Soltu.DM.07G013370,
                                                                                               Soltu.DM.10G013240,
                                                                                               Soltu.DM.12G004870,
                                                                                               Soltu.DM.08G009720,
                                                                                               Soltu.DM.09G031820,
                                                                                               Soltu.DM.05G024440))

zscore_boxplot_starch_synthesis_long <- zscore_boxplot_starch_synthesis_long%>%
  rename("ZScore" = value)

ggplot(zscore_boxplot_starch_synthesis_long, aes(x=genotype, y=ZScore, fill=conditions))+
  geom_boxplot()+
  facet_wrap(~factor(name, levels = c("Soltu.DM.08G030230",
                                      "Soltu.DM.12G026390",
                                      "Soltu.DM.02G020800",
                                      "Soltu.DM.03G019120",
                                      "Soltu.DM.07G022290",
                                      "Soltu.DM.09G004100",
                                      "Soltu.DM.07G013360",
                                      "Soltu.DM.07G013370",
                                      "Soltu.DM.10G013240",
                                      "Soltu.DM.12G004870",
                                      "Soltu.DM.08G009720",
                                      "Soltu.DM.09G031820",
                                      "Soltu.DM.05G024440")), nrow = 2)+
  scale_fill_manual(values = c("#0072B2", "#D55E00"))+
  theme_bw()+
  theme(axis.title = element_text(size=20),
              axis.text = element_text(size=14),
              strip.text = element_text(size = 14, face = "bold"),
              legend.title = element_text(size=14),
            legend.text = element_text(size = 14),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"))


starch_degradation_DEGs <- read.csv("/media/rna/NERPA/Senta/tables/starch_degradation_DEGs.csv", sep = "\t")

starch_degradation_DEGs <- starch_degradation_DEGs%>%
  column_to_rownames("Geneid")

zscore_boxplot_starch_degradation <- zscores_boxplot_starch_DEGs%>%
  select(rownames(starch_degradation_DEGs), genotype, conditions)

zscore_boxplot_starch_degradation_long <- pivot_longer(zscore_boxplot_starch_degradation, cols = c(Soltu.DM.03G024560,
                                                                                                   Soltu.DM.08G001120,
                                                                                                   Soltu.DM.04G025730,
                                                                                                   Soltu.DM.04G037250,
                                                                                                   Soltu.DM.09G011580,
                                                                                                   Soltu.DM.09G027770,
                                                                                                   Soltu.DM.05G000570,
                                                                                                   Soltu.DM.02G017070,
                                                                                                   Soltu.DM.04G033700,
                                                                                                   Soltu.DM.05G006330,
                                                                                                   Soltu.DM.07G018100,
                                                                                                   Soltu.DM.08G023420))


zscore_boxplot_starch_degradation_long <- zscore_boxplot_starch_degradation_long%>%
  rename("ZScore" = value)

ggplot(zscore_boxplot_starch_degradation_long, aes(x=genotype, y=ZScore, fill=conditions))+
  geom_boxplot()+
  facet_wrap(~factor(name, levels = c("Soltu.DM.03G024560",
                                      "Soltu.DM.08G001120",
                                      "Soltu.DM.04G025730",
                                      "Soltu.DM.04G037250",
                                      "Soltu.DM.09G011580",
                                      "Soltu.DM.09G027770",
                                      "Soltu.DM.05G000570",
                                      "Soltu.DM.02G017070",
                                      "Soltu.DM.04G033700",
                                      "Soltu.DM.05G006330",
                                      "Soltu.DM.07G018100",
                                      "Soltu.DM.08G023420")), nrow = 2)+
  scale_fill_manual(values = c("#0072B2", "#D55E00"))+
  theme_bw()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=14),
        strip.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size=14),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"))


starch_undefined_DEGs <- read.csv("/media/rna/NERPA/Senta/tables/starch_undefined_DEGs.csv", sep = "\t")

starch_undefined_DEGs <- starch_undefined_DEGs%>%
  column_to_rownames("Geneid")

zscore_boxplot_starch_undefined <- zscores_boxplot_starch_DEGs%>%
  select(rownames(starch_undefined_DEGs), genotype, conditions)

zscore_boxplot_starch_undefined_long <- pivot_longer(zscore_boxplot_starch_undefined, cols = Soltu.DM.10G004860)

zscore_boxplot_starch_undefined_long <- zscore_boxplot_starch_undefined_long%>%
  rename("ZScore" = value)

ggplot(zscore_boxplot_starch_undefined_long, aes(x=genotype, y=ZScore, fill=conditions))+
  geom_boxplot()+
  facet_wrap(~factor(name, levels = c("Soltu.DM.10G004860")), nrow = 1)+
  scale_fill_manual(values = c("#0072B2", "#D55E00"))+
  theme_bw()+
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=20),
        strip.text = element_text(size = 20, face = "bold"),
        legend.title = element_text(size=24),
        legend.text = element_text(size = 24),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold"),
        axis.title.x = element_blank())

