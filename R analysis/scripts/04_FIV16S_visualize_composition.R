### community taxa visualization for FIV 16S study ###

################# data prep #####################

# load libraries
library(tidyverse)
library(vegan)
library(phyloseq)
library(RColorBrewer)

# set working directory


# load functions to plot and look at our data (description is here and in the "functions.R" file)
source(file="02_scripts/functions_from_Zaid_Dec22.R")
source(file="02_scripts/functions_LB.R")

######## load and prep data #########
data.df <- read.csv("03_output/data.df.csv", header=T, row.names=1)
meta.df <- read.csv("03_output/meta.df.csv", header=T, row.names=1) %>%
  mutate(
         treatment = factor(Treatment, c("Control", "Placebo", "cART")),
         Cat_ID = factor(Cat_ID, c( "Cat1", "Cat2", "Cat3", "Cat4", "Cat5", "Cat6", #ctrl
                                    "Cat7", "Cat8", "Cat9", "Cat10", "Cat11", "Cat12", #placebo
                                   "Cat13", "Cat14", "Cat15", "Cat16", "Cat17", "Cat18"))) #cART

taxa.df <- read.csv("03_output/taxa.df.csv", header=T, row.names=1)

data.df <- data.df[row.names(meta.df),] %>%
  select_if(colSums(.) > 0) #remove any observations that were only in the NC/PC samples

taxa.df <- taxa.df[colnames(data.df),]

  
############# plot relative abundance #################

# code from KK
#### Filter OTU data according to KK recommendation

# transpose
counts <- as.data.frame(t(data.df))

# filter: OTU should have at least 10 reads and be in at least 5 samples 
filtered.counts <- counts[rowSums(counts>10)> 5, ] 

# save
write.csv(filtered.counts,"03_output/FIV16S_filtered_counts.csv")

#### Calculate Relative Abundance

# calculate 
ra.df <- filtered.counts
ra.df[] <- lapply(ra.df, function(x) x/sum(x))

# lil sanity check that relabun calc is working
sum(ra.df$G41PR16)
15571/68787

# save
write.csv(ra.df, "03_output/FIV16S_filtered_relabun.csv")


#### get taxa 

# clean up: reset (re-transpose) ra.df with relative abundances
data.ra.df <- as.data.frame(t(ra.df))

# just want taxa that passed filtering
taxa.ra.df <- taxa.df[colnames(data.ra.df),]

######### phyloseq container
########### We use the data.df (the otu table), taxa.df (the taxa table) and the experimental design meta.df (sample data)
data.ra.ps = phyloseq(otu_table(data.ra.df,taxa_are_rows = FALSE),tax_table(as.matrix(taxa.ra.df)),sample_data(meta.df))

#### Phylum level (the otu table and the taxa table are now on the phylum level) we call the new container data.ps.p
data.ps.p = tax_glom(data.ra.ps,"phylum")
#### We extract the new collapsed otu table and put it into a new data frame data.df.p using function otu_table in phyloseq 
data.df.p = data.frame(otu_table(data.ps.p))
#### We extract the new taxonomy table and save it in a new data frame taxa.df.p using function tax_table in phyloseq
taxa.df.p = data.frame(tax_table(data.ps.p)[,1:5])
#### We rename the columns of the data.df.f from otu00001, otu00004, ... to the phylum names
colnames(data.df.p) = taxa.df.p$phylum
#### We also extract the experimental design (which will not change) into data.df.p using function sample_data in phyloseq
meta.df.p = data.frame(sample_data(data.ps.p))


## tables & plotting from KK 

# prep for plotting/making tables
phylum_relabun <- data.df.p %>%
  rownames_to_column("SeqID") %>%
  left_join(meta.df, by = "SeqID") %>%
  select(-c("Sample", "dataset_ID", "treatment")) %>%
  pivot_longer(!c("SeqID", "Treatment", "Week", "Cat_ID"), names_to = "phylum", 
               values_to = "relabun") 

# table of predominant phyla

# calculate mean values by treatment and week
means <- phylum_relabun %>%
  group_by(Treatment, Week, phylum) %>%
  summarize(meanRelAbun = mean(relabun))

# week -1
week0 <- means %>%
  filter(Week == "-1") %>%
  slice_max(meanRelAbun, n=5)

# week 5
week4 <- means %>%
  filter(Week == "5") %>%
  slice_max(meanRelAbun, n=5)

# week 11
week12 <- means %>%
  filter(Week == "11") %>%
  slice_max(meanRelAbun, n=5)

# week 24
week25 <- means %>%
  filter(Week == "24") %>%
  slice_max(meanRelAbun, n=5)

# plot

palette_phylum <- c("#E31A1C", "#FFFF99","#B2DF8A" ,
                    "#FDBF6F", "#6A3D9A", "#A6CEE3",
                    "#1F78B4")

p <- ggplot(phylum_relabun, 
            aes(x=factor(Week) , y= relabun, fill=phylum)) +
  scale_fill_manual(values = c("#E31A1C", "#FFFF99","#B2DF8A" ,
                               "#FDBF6F", "#6A3D9A", "#A6CEE3",
                               "#1F78B4")) +
  geom_bar(stat= "identity", width =0.75) +
  facet_wrap(~Treatment,scales="free", ncol=3) +
  labs(x="Week", y= "Relative abundance", fill="Phylum") +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        strip.text = element_text(size=12),
        axis.title = element_text(size=14))


pdf("03_output/figures/barplot_relabun_phylum_by_trt_wk.pdf", width=8, height=4)
p
dev.off()



########### plot composition of study samples ############

# rename week in metadata

meta.df <- meta.df %>%
  mutate(Week = fct_recode(factor(Week), 
                           "Week -1" ="-1", "Week 5"="5", 
                           "Week 11"="11", "Week 24"="24"))


######### phyloseq container
########### We use the data.df (the otu table), taxa.df (the taxa table) and the experimental design meta.df (sample data)
data.ps = phyloseq(otu_table(data.df,taxa_are_rows = FALSE),tax_table(as.matrix(taxa.df)),sample_data(meta.df))

####### It is cleaner to look at the barplots using a level higher than the otu level
####### we can choose to collapse the otu tables (finding the counts per genus, family, order ... etc) using the 
####### tax_glom() function in phyloseq and this is why we created the phyloseq container data.ps
#### Here we apply the function tax_glom to data.ps to create a new phyloseq container that collapses the data to the 
#### family level (the otu table and the taxa table are now on the family level) we call the new container data.ps.f
data.ps.f = tax_glom(data.ps,"family")
#### We extract the new collapsed otu table and put it into a new data frame data.df.f using function otu_table in phyloseq 
data.df.f = data.frame(otu_table(data.ps.f))
#### We extract the new taxonomy table and save it in a new data frame taxa.df.f using function tax_table in phyloseq
taxa.df.f = data.frame(tax_table(data.ps.f)[,1:5])
#### We rename the columns of the data.df.f from otu00001, otu00004, ... to the family names
colnames(data.df.f) = taxa.df.f$family
#### We also extract the experimental design (which will not change) into data.df.f using function sample_data in phyloseq
meta.df.f = data.frame(sample_data(data.ps.f))

colorCount = length(unique(taxa.df.f$family))  #define number of variable colors based on number of Family (change the level accordingly to phylum/class/order)
getPalette = colorRampPalette(brewer.pal(12, "Paired")) (colorCount)  # change the palette as well as the number of colors will change according to palette.
Palette = sample(getPalette)


p.fam <- cat_plot_ftn(data = data.df.f, Trt = meta.df.f[,c("Week", "Cat_ID")],
             cutoff = 0.01, llab = "Family", palette=Palette) +
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_text(size=10, angle=45),
        axis.text.y = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.35, "cm"))

pdf("03_output/figures/barplot_families_by_trt_wk.pdf", width=8, height=10)
p.fam
dev.off()

# save palette
write.csv(Palette, "03_output/family_barplot_palette.csv")



########## plot composition of negative control samples #########

# read in negative control data
ng.df <- read.csv("03_output/ng.data.df.csv", header=T, row.names=1)
ng.taxa <- read.csv("03_output/ng.taxa.df.csv", header=T, row.names=1) 
ng.meta <- read.csv("03_output/ng.meta.df.csv", header=T, row.names=1) %>%
  mutate(id = factor(id, c("NCon116", "NCon216", "NCon316"),
                     labels = c("Negative control 1", "Negative control 2", "Negative control 3")))

# let's simplify to the family level
ng.ps = phyloseq(otu_table(ng.df, taxa_are_rows = FALSE), tax_table(as.matrix(ng.taxa)), sample_data(ng.meta))

ng.ps.f = tax_glom(ng.ps,"family")
#### We extract the new collapsed otu table and put it into a new data frame data.df.f using function otu_table in phyloseq 
ng.data.f = data.frame(otu_table(ng.ps.f))
#### We extract the new taxonomy table and save it in a new data frame taxa.df.f using function tax_table in phyloseq
ng.taxa.f = data.frame(tax_table(ng.ps.f)[,1:5])
#### We rename the columns of the data.df.f from otu00001, otu00004, ... to the family names
colnames(ng.data.f) = ng.taxa.f$family
#### We also extract the experimental design (which will not change) into data.df.f using function sample_data in phyloseq
ng.meta.f = data.frame(sample_data(ng.ps.f))

#pal <- sample(brewer.pal(11, "Paired"))
pal <- c("#E31A1C", "#FF7F00", "#B2DF8A" , "#FFFF99", "#33A02C",
         "#A6CEE3", "#1F78B4", "#FDBF6F","#CAB2D6" ,"#6A3D9A","#FB9A99")
NC_plot_ftn(ng.data.f, Trt=ng.meta.f[,c("Sample", "id")],
            cutoff = 0.01, llab = "Family", palette=pal)

pdf(file="03_output/figures/barplots_NC_families_20240208.pdf",width=9,height=5)
par(cex=1, cex.axis = 1.2, cex.lab = 1.2)
p.ng + theme(strip.text=element_blank(),
             text = element_text(size=16))
dev.off()

