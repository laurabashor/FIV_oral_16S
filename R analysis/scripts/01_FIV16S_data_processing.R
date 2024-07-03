### Data processing for FIV 16S study ###

# Data from mothur run test6 2022-10-26, trimmomatic trimming
# Script modified from Dr. Zaid Abdo 2022-12-14

# Study design:
# week naming convention: Prebleed = Week -1
# Cat IDs 1-18 (Control 1-6, Placebo 7-12, cART 13-18)

################# load packages #####################

library(tidyverse)
library(phyloseq)
library(metagenomeSeq)
library(vegan)
library(openxlsx)

################# initial data processing #####################

# set working directory

############# Load mothur results
data.df <- read.table(file="01_input/final_test6.dgc.shared",header=TRUE)

#### loading taxonomy using phyloseq in different columns
#row names of taxa.df need to match column names of data.df
taxa.phs <- import_mothur(mothur_constaxonomy_file = "01_input/final_test6.dgc.0.03.cons.taxonomy")
colnames(tax_table(taxa.phs)) = c("kingdom", "phylum", "class", "order", "family",  "genus")
taxa.df <- data.frame(tax_table(taxa.phs))

meta.df <- read.csv("01_input/FIV_microbiome_design.csv") %>%
  mutate(Week = factor(Week, level = c("-1", "5", "11", "24", "")),
         dataset_ID = case_when(grepl("Cat", Cat_ID) ~ paste0(Cat_ID, "_wk", Week),
                                grepl("Positive", Cat_ID) ~ Cat_ID,
                             grepl("Negative", Cat_ID) ~ Cat_ID))
            
# set row names     
row.names(meta.df) <- meta.df[ ,1]

# save the sample names (column # 2) from the current data.df data frame to become row names
r.nm <- as.character(data.df[ , 2])
# remove first three columns (label/group/otu number) relabel row names with r.nm
data.df <- data.frame(data.df[ ,-c(1:3)], row.names=r.nm)

otu.nm <- colnames(data.df) # save otu names

# The following makes sure that the order of the rows is consistent
taxa.df <- taxa.df[row.names(taxa.df) %in% otu.nm, ]

# look at raw sequencing depth across samples
depth <- apply(data.df, 1, sum)

hist(depth) #pretty heterogeneous for depth
min(depth) #19590
max(depth) #103824
mean(depth) #61279.59
# both mock and neg ctrl seemed to get too high sequencing depth


######################## look at negative control data ################

### pull out all data negative control samples (where there are counts: >0 or !=0 is the same here)
ng.df <- data.df %>%
  filter(rownames(.) %in% c("NCon116", "NCon216", "NCon316")) %>% 
  select_if(colSums(.) > 0)

length(ng.df) #1936 taxa all together

# make a dataframe for each individual negative control sample
ng1 <- ng.df %>% filter(rownames(.) == "NCon116") %>% select_if(colSums(.) > 0) 
ng2 <- ng.df %>% filter(rownames(.) == "NCon216") %>% select_if(colSums(.) > 0)
ng3 <- ng.df %>% filter(rownames(.) == "NCon316") %>% select_if(colSums(.) > 0)

# how many taxa and sequences in each negative control sample?
length(ng1) #569 taxa
sum(ng1) #49778 sequences
length(ng1[ng1>5]) #188 taxa

length(ng2) #995 taxa
sum(ng2) #54532 sequences
length(ng2[ng2>5]) #210 taxa

length(ng3) #1060 taxa
sum(ng3) #19590
length(ng3[ng3>5]) #188 taxa


################### Zymo positive controls to process samples with mock cutoff #########################

### expected Zymo community (12% of each of the eight bacterial species)
#Listeria monocytogenes, Pseudomonas aeruginosa, Bacillus subtilis, 
#Escherichia coli, Salmonella enterica, Lactobacillus fermentum, 
#Enterococcus faecalis, Staphylococcus aureus
#also fungus: Saccharomyces cerevisiae - 2%, and Cryptococcus neoformans - 2% 

# this community includes three bacterial orders:
#(1) Bacillales (listeria, bacillus, staphylococcus), 
#(2) Enterobacterales (escherichia, salmonella)
#(3) Lactobacillales (enterococcus, lactobacillus)
#(4) Pseudomonadales (pseudomonas)

### now pull out data in Zymo mock samples #expecting eight specific taxa 
mock.df <- data.df %>%
  filter(rownames(.) %in% c("PCon116", "PCon216", "PCon316")) %>%
  select_if(colSums(.) > 0)

sum(mock.df[1,]) #82492 sequences in PCon116
sum(mock.df[2,]) #49383 sequences in PCon216
sum(mock.df[3,]) #65427 sequences in PCon316

ncol(mock.df) #1590 taxa, way too many! need to find a cutoff that will reduce 1590 to 8-10 or a few more, with the expected taxa
# if there are more than the number that should be in the Mock community (Zymo = about 8-10 taxa), we are overestimating the diversity
# some of these OTUs might have just had 1 read; we don't want any OTU columns that didn't have more than one read per sample

#### We choose a cutoff (a number that we use to subtract from all counts)
mock.cutoff <- mock.df - 102 #cutoff
mock.cutoff[mock.cutoff < 0] <- 0
mock.cutoff <- mock.cutoff %>%
  select_if(colSums(.) > 0)
ncol(mock.cutoff) 
#704 taxa if cutoff is 1, 171 taxa if cutoff is 5, 20 taxa if cutoff is 55 (look at these taxa)
#final cutoff: 12 OTUs remain if cutoff is 102

# match col names of mock with col names of taxonomy table
# when you change the cutoff, look at this again to see if you kept the expected taxa, but got rid of unknown/incorrect

m.df <- taxa.df[colnames(mock.cutoff),] # removed previous code for ordering mock.df columns, unnecessary and not functional in current version of R
g.lvls <- levels(as.factor(as.character(m.df$genus)))
g.lvls

#if we use cutoff 55 (20 taxa), we still have two unwanted taxa: "Proteocatella" and "Comamonadaceae_unclassified"
#"Proteocatella"  is in family Peptostreptococcaceae in order Eubacteriales (we should NOT have this)
#"Comamonadaceae" is a family in order Burkholderiales (we should NOT have this)

# it takes all the way to cutting off at 102 (12 OTUs, 11 taxa), we finally lose the bad taxa:
#"Bacillaceae_1_unclassified"      "Bacilli_unclassified"            "Enterobacterales_unclassified"  
#"Enterobacteriaceae_unclassified" "Enterococcaceae_unclassified"    "Escherichia/Shigella"           
#"Lactobacillales_unclassified"    "Limosilactobacillus"             "Listeria"                       
#"Pseudomonas"                     "Staphylococcus"                 

# so we have all the taxa that we expect now

#### Using what we see from the above command we apply the cutoff to the whole dataset 
# subtract that number from all counts
d.df <- data.df - 102
d.df[d.df<0] <- 0
d.df <- d.df %>%
  select_if(colSums(.) > 0)

# reset taxa dataframe with more limited dataset
t.df <- taxa.df[colnames(d.df),]

### how is seq depth now? after applying mock cutoff
rs.vt <- apply(d.df, 1, sum) # row sum (total reads per sample)
depth <- as.data.frame(rs.vt) %>%
  rename(depth = rs.vt) %>%
  arrange(depth)

median(depth$depth) #47585.5
range(depth$depth) # 9223 85797
#there is one sample that is <10,000, but overall okay

# plot sequencing depth after mock cutoff -- this includes positive and negative controls
# pdf(file="03_output/figures/sequencing_depth_after_mock_cutoff_20230303.pdf",width=4,height=4)
# hist(rs.vt, main = "Distribution of sequencing depth \n(after data processing)", 
#      xlab = "Sequencing depth", ylab = "Number of samples")
# dev.off()


######################### assess negative control samples after mock processing ########################

### pull out all data negative control samples (where there are counts: >0 or !=0 is the same here)
ng.df <- d.df %>%
  filter(rownames(.) %in% c("NCon116", "NCon216", "NCon316")) %>% 
  select_if(colSums(.) > 0)

length(ng.df) #37 taxa all together
num.ng.tax <- rowSums(ng.df != 0) # 13 in NCon1, 29 in NCon2, 8 in NCon3

# make a dataframe for each individual negative control sample
ng1 <- ng.df %>% filter(rownames(.) == "NCon116") %>% select_if(colSums(.) > 0)
ng2 <- ng.df %>% filter(rownames(.) == "NCon216") %>% select_if(colSums(.) > 0)
ng3 <- ng.df %>% filter(rownames(.) == "NCon316") %>% select_if(colSums(.) > 0)

# how many taxa in each negative control sample now?
length(ng1) #13 taxa
sum(ng1) #44686 sequences
length(ng2) #29 taxa (looking at BLAST, NCon216 had some oral/intestinal microbiota...)
sum(ng2) #46983
length(ng3) #8 taxa
sum(ng3) #13889 sequences

# ideally there should be no taxa in the negative control, or we are overestimating the diversity/counting contamination as part of the microbiome
# however, we need to report what was in our negative controls

# pull NC taxa from taxa.df
ng.taxa <- t.df[colnames(ng.df),]

# make data counts in long format
ng.long <- ng.df %>%
  filter(rownames(.) %in% c("NCon116", "NCon216", "NCon316")) %>%
  rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, 
              values_from=value) 

# combine taxa with long counts; this dataframe includes taxonomy, OTU names, and counts for all three NCs
ng.long <- ng.taxa %>%
  rownames_to_column(var = "name") %>%
  left_join(ng.long, by = "name")

# pull out just NC metadata
ng.meta <- meta.df %>%
  filter(Treatment == "negative_control")

# save negative control data post-processing
write.csv(ng.df, "03_output/ng.data.df.csv", row.names=T)
write.csv(ng.meta, "03_output/ng.meta.df.csv", row.names=T)
write.csv(ng.taxa, "03_output/ng.taxa.df.csv", row.names=T)

#write.xlsx(ng.long, "03_output/negative_control_taxa_and_counts.xlsx")

# re-set row names
d.df = d.df[row.names(meta.df),]

#### Now data is clean and well organized and we replace data.df with d.df and taxa.df with t.df 
data.df <- d.df
taxa.df <- t.df

#### write csvs
# these tables publicly available in Github repo https://github.com/laurabashor/FIV_oral_16S/
write.csv(data.df, "03_output/github/data.df.csv", row.names=T)
write.csv(meta.df, "03_output/github/meta.df.csv", row.names=T)
write.csv(taxa.df, "03_output/github/taxa.df.csv", row.names=T)



######## generate rarefaction curves ###########################

#rarefaction curve shows how many OTUs you're actually working with for each sample
#bootstrap sample of OTUs observed based on depth of coverage
#we're finding expected values, not actually "rarefying" the data
#if you see the curve flatten out, you know that sequencing more won't give you any more information

# first, prep dataset ids with just cat number and date
dataset_id <- meta.df %>% pull(dataset_ID)

# assign dataset id to rownames for plotting
row.names(meta.df) <- dataset_id

# make plotting df with dataset id row names
plot.df <- d.df
row.names(plot.df) <- dataset_id

# re-do rs.vt to get sequencing depths per sample
rs.vt = apply(data.df, 1, sum) # row sum (total reads per sample)

##### Create a .pdf file called rarefaction_curves.pdf with dimensions 10x10in (that you can change depending 
##### on the size you want and the resolution you have) to save the rarefaction curves
pdf(file="03_output/figures/rarefaction_curves.pdf", width = 10, height = 10)
# ##### rarecurve() is a function in vegan it takes the otu table x=data.df, and to have different colors per samples
# ##### we use col and set these colors to the numbers from 1 to the total number of samples "col=1:length(rs.vt)"
# ##### label of the x axes is "Number of reads" and that for the y axes is "OTU"
par(cex=0.5, cex.axis = 2, cex.lab = 2, mar=c(5,5,5,5)+0.1) #set font size for sample labels and margins for pdf plot
rarecurve(x=plot.df,col=1:length(rs.vt),
          xlab = "Number of reads", ylab="OTU")
# #### We save the file
dev.off()

# return rownames to previous
row.names(meta.df) = meta.df[ ,1]


########################## remove negative & positive control samples from dataset #########################

#### Now we are ready to remove the negative and positive control data from our dataset
# as we have used the positive controls, and we will analyze the negatives separately
data.df <- data.df %>%
  filter(!rownames(.) %in% c("NCon116", "NCon216", "NCon316")) %>%
  filter(!rownames(.) %in% c("PCon116", "PCon216", "PCon316")) %>%
  select_if(colSums(.) > 0) #remove any observations that were only in the NC/PC samples


########################## final plots & save processed data ##########################

# make sure d.df matches row names in metadata & ncon data is gone
meta.df <- meta.df %>%
  filter(!Treatment %in% c("negative_control", "positive_control")) %>%
  rename("SeqID" = id)

data.df = data.df[row.names(meta.df),]

# reset taxa dataframe with more limited dataset
taxa.df = taxa.df[colnames(data.df),]


#### Now get median sequencing depth and number of taxa
rs.vt = apply(data.df, 1, sum) # row sum (total reads per sample)
depth <- as.data.frame(rs.vt) %>%
  rename(depth = rs.vt) %>%
  arrange(depth)

median(depth$depth) #47769
range(depth$depth)
# 9223 85797 same range as before (ie lowest and highest samples are not the negative or positive controls)
# there is one sample that is <10,000, but overall okay


## what is the median #taxa in a sample? the range of #taxa?
num.tax <- rowSums(data.df != 0)
median(num.tax) #49
range(num.tax) # 25 - 105

# # sequencing depth after mock cutoff and removal of positive and negative controls
pdf(file="03_output/figures/sequencing_depth_after_mock_cutoff_20240703.pdf",width=4,height=4)
hist(rs.vt, main = "Distribution of sequencing depth \n(after data processing)",
     xlab = "Sequencing depth", ylab = "Number of samples")
dev.off()


########### write final csvs (after Mock cutoff AND removing NC/PC samples)
# write these into tables so we don't have to run the whole script again
write.csv(data.df, "03_output/data.df.csv", row.names=T)
write.csv(meta.df, "03_output/meta.df.csv", row.names=T)
write.csv(taxa.df, "03_output/taxa.df.csv", row.names=T)

# clean up; remove everything you don't need anymore
rm(taxa.phs, d.df,t.df, m.df, mock.df, g.lvls, otu.nm, r.nm, depth, ng.df, ng.meta, ng.taxa,ng.df.long,ng1,ng2,ng3, rs.vt, ncon.otu.nm, dataset_id, num.ng.tax, num.tax,plot.df)

