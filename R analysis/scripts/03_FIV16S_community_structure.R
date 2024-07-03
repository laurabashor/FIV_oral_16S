### PERMANOVA & ordination for FIV 16S study ###

# Laura Bashor
# code modified from Drs. Zaid Abdo and Scott Carver

################# data prep #####################

# set working directory

# load libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(metagenomeSeq)

# read in processed data
data.df <- read.csv("03_output/data.df.csv", header=T, row.names=1)
meta.df <- read.csv("03_output/meta.df.csv", header=T, row.names=1 )%>%
  mutate(Week = as.factor(Week)) %>%
  rename(CatID = Sample)
taxa.df <- read.csv("03_output/taxa.df.csv", header=T, row.names=1)

# use the data.df (the otu table), taxa.df (the taxa table) and the experimental design meta.df to create a new phyloseq container
data.ps = phyloseq(otu_table(data.df,taxa_are_rows = FALSE),tax_table(as.matrix(taxa.df)),sample_data(meta.df))


################# normalization (code from Zaid Abdo) #####################

######## We use normalized data to do ordination and other analyses
##### We normalize using the Cumulative Sum Scaling (CSS) of metagenome Seq and to do so we need to create
##### a metagenomeSeq experiment using the function newMRexperiment()
##### Before that we have to do some house keeping by changing the format of the meta.df (the experimental design)
##### and the taxa.df (the taxonomic assignment) to a formate acceptable by metagenomeSeq (we use AnnotateDataFtame() function)
trt.an = AnnotatedDataFrame(meta.df)
taxa.an = AnnotatedDataFrame(taxa.df)
##### we can now create a metagenomeSeq container using newMRexperiment using the transposed raw data (t(data.df))
##### and the newly formated experimental design (trt.an) and taxonomic data (tax.an)
d.ms = newMRexperiment(t(data.df), phenoData = trt.an, featureData = taxa.an)
##### To see what is in the new metagenomeSeq container:
##### We use pData() to look at the experimental design
pData(d.ms)
##### We use fData to look at the taxonomic assignment
fData(d.ms)
##### We use MRcount() to look at the otu table (here head() shows the first 6 lines of that table)
head(MRcounts(d.ms))


##### We normalize using the CSS normalization (using this all way through)
### We use the new container in function cumNormStatFast() to find the percentile that we will use to 
### normalize based on and we add the result in variable p
p = cumNormStatFast(d.ms)
### We adjust our container to contain the normalization pecentile using cumNorm()
d.ms = cumNorm(d.ms,p)
### We find the normalization factors per sample using normFactors()
nf = normFactors(d.ms)
### We extract the normalized matrix using cumNormMat()
norm.mt = cumNormMat(d.ms)
### This matrix has rows as otus and columns as samples (all our work before used columns as otus)
## we first remove all rows (otus) that have count zero by summing per row and keeping those rows with
## sum greater than zero
norm.rs = apply(norm.mt, 1, sum)
norm.mt = norm.mt[norm.rs > 0, ]
## We then transform the normalized matrix so that rows are samples and columns are otus (for vegan to work correctly)
norm.mt = t(norm.mt)

####### We need to create a phyloseq container based on the normalized data
#### Making sure that the normalized data is in a data frame format
data.df.nrm = data.frame(norm.mt)
#### choosing only the taxa observed in the normalized data to include in the phyloseq container
taxa.mt.nrm = as.matrix(taxa.df[colnames(data.df.nrm),])
#### the experimental data didn't change so we add it as is
norm.ps = phyloseq(otu_table(data.df.nrm,taxa_are_rows = FALSE),tax_table(taxa.mt.nrm),sample_data(meta.df))

################ one PERMANOVA for all data (code from Scott Carver) #####################

# Permutational-based Multivariate Analysis of Variance
# ANOVA on multivariate/community data
# use the ‘adonis’ function in the ‘vegan’ package

# extract the normalized data frame from the phyloseq container, 
# and combine it with the metadata for a complete dataframe

df <- as.data.frame(norm.ps@otu_table) %>%
  rownames_to_column("id") %>%
  left_join(meta.df, by = c("id" = "SeqID")) %>%
  relocate(1, .after = last_col()) # move id column to the end so OTU data is first

write_csv(df, "03_output/permanova_df.csv")

# remove the final 4 columns containing metadata just to keep OTU data
dat <- df[,1:(length(df)-4)] 

# extract the Treatment and Week columns for testing with permanova
trt <- df %>% pull(Treatment)
wk <- df %>% pull(Week)
id <- df %>% pull(CatID)

permanova <- adonis(dat~wk+trt+wk*trt, strata=id, permutations=999) 

permanova[["aov.tab"]] 


#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# wk         3    5.2107 1.73691 11.5727 0.33022  0.001 ***
# trt        2    0.3742 0.18709  1.2465 0.02371  0.001 ***
# wk:trt     6    1.1893 0.19822  1.3207 0.07537  0.169    
# Residuals 60    9.0052 0.15009         0.57069           
# Total     71   15.7794                 1.00000          


################ ordination (code from Scott Carver) #####################

nmds = metaMDS(dat, k=2) 

labs = c(rep("Placebo",6),rep("Control",6),rep("cART",6))
#cols = c(rep("#7570b3",6),rep("#1b9e77",6),rep("#d95f02",6))
cols = c(rep("black",6),rep("blue",6),rep("red",6))


pdf("03_output/figures/NMDS.pdf")
ordiplot(nmds, type = "n", 
         main = "16S microbial community structure", 
         ylim = c(-1,1), xlim = c(-1,1))
orditorp(nmds, display = "sites", 
         labels = labs, col = cols,
         cex=1)
ordiellipse(nmds, display = "sites", groups=trt1, col="grey")
dev.off()


############ post-hoc tests ##################

#### week 24 ####

# extract just Week 24 OTU data for testing with permanova
wk24_dat <- df %>%
  filter(Week == 24)

#cART vs. placebo 
perm1 <- adonis(wk24_dat[-which(wk24_dat$Treatment == "Control"), 1:(length(wk24_dat)-4)] ~
                  wk24_dat[-which(wk24_dat$Treatment=="Control"), "Treatment"], permutations=999)
perm1[["aov.tab"]] 

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# wk24_dat[-which(wk24_dat$Treatment == "Control"), "Treatment"]  1   0.61477 0.61477  4.1043 0.29099  0.005 **
# Residuals                                                      10   1.49788 0.14979         0.70901          
# Total                                                          11   2.11265                 1.00000          

#cART vs. Control
perm2 <- adonis(wk24_dat[-which(wk24_dat$Treatment=="Placebo"), 1:(length(wk24_dat)-4)]~
                  wk24_dat[-which(wk24_dat$Treatment=="Placebo"), "Treatment"], permutations=999)
perm2[["aov.tab"]] 

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# wk24_dat[-which(wk24_dat$Treatment == "Placebo"), "Treatment"]  1   0.35341 0.35341  2.0371 0.16924  0.046 *
# Residuals                                                      10   1.73485 0.17349         0.83076         
# Total                                                          11   2.08826                 1.00000         

#placebo vs. Control

perm3 <- adonis(wk24_dat[-which(wk24_dat$Treatment=="cART"), 1:(length(wk24_dat)-4)]~
                  wk24_dat[-which(wk24_dat$Treatment=="cART"), "Treatment"], permutations=999)
perm3[["aov.tab"]]

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# wk24_dat[-which(wk24_dat$Treatment == "cART"), "Treatment"]  1   0.14077 0.14077 0.94496 0.08634  0.488
# Residuals                                                   10   1.48968 0.14897         0.91366       
# Total                                                       11   1.63045                 1.00000       


#### week 11 ####

# extract just Week 11 OTU data for testing with permanova
wk11_dat <- df %>%
  filter(Week == 11)

#cART vs. placebo 
perm1 <- adonis(wk11_dat[-which(wk11_dat$Treatment == "Control"), 1:(length(wk11_dat)-4)] ~
                  wk11_dat[-which(wk11_dat$Treatment=="Control"), "Treatment"], permutations=999)
perm1[["aov.tab"]] 

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# wk11_dat[-which(wk11_dat$Treatment == "Control"), "Treatment"]  1   0.14033 0.14033 0.60312 0.05688  0.923
# Residuals                                                      10   2.32677 0.23268         0.94312       
# Total                                                          11   2.46711                 1.00000       

#cART vs. Control
perm2 <- adonis(wk11_dat[-which(wk11_dat$Treatment=="Placebo"), 1:(length(wk11_dat)-4)]~
                  wk11_dat[-which(wk11_dat$Treatment=="Placebo"), "Treatment"], permutations=999)
perm2[["aov.tab"]] 

# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
# wk11_dat[-which(wk11_dat$Treatment == "Placebo"), "Treatment"]  1   0.24245 0.24245   1.073 0.0969  0.351
# Residuals                                                      10   2.25966 0.22597         0.9031       
# Total                                                          11   2.50211                 1.0000       

#placebo vs. Control

perm3 <- adonis(wk11_dat[-which(wk11_dat$Treatment=="cART"), 1:(length(wk11_dat)-4)]~
                  wk11_dat[-which(wk11_dat$Treatment=="cART"), "Treatment"], permutations=999)
perm3[["aov.tab"]]

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# wk11_dat[-which(wk11_dat$Treatment == "cART"), "Treatment"]  1   0.12016 0.12016 0.61251 0.05772  0.853
# Residuals                                                   10   1.96177 0.19618         0.94228       
# Total                                                       11   2.08193                 1.00000       


#### week 5 ####

# extract just Week 5 OTU data for testing with permanova
wk5_dat <- df %>%
  filter(Week == 5)

#cART vs. placebo 
perm1 <- adonis(wk5_dat[-which(wk5_dat$Treatment == "Control"), 1:(length(wk5_dat)-4)] ~
                  wk5_dat[-which(wk5_dat$Treatment=="Control"), "Treatment"], permutations=999)
perm1[["aov.tab"]] 

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# wk5_dat[-which(wk5_dat$Treatment == "Control"), "Treatment"]  1   0.24722 0.24722  2.2714 0.18509  0.047 *
# Residuals                                                    10   1.08842 0.10884         0.81491         
# Total                                                        11   1.33564                 1.00000         

#cART vs. Control
perm2 <- adonis(wk5_dat[-which(wk5_dat$Treatment=="Placebo"), 1:(length(wk5_dat)-4)]~
                  wk5_dat[-which(wk5_dat$Treatment=="Placebo"), "Treatment"], permutations=999)
perm2[["aov.tab"]] 

# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# wk5_dat[-which(wk5_dat$Treatment == "Placebo"), "Treatment"]  1   0.05784 0.057839   0.625 0.05882  0.806
# Residuals                                                    10   0.92542 0.092542         0.94118       
# Total                                                        11   0.98326                  1.00000       

#placebo vs. Control

perm3 <- adonis(wk5_dat[-which(wk5_dat$Treatment=="cART"), 1:(length(wk5_dat)-4)]~
                  wk5_dat[-which(wk5_dat$Treatment=="cART"), "Treatment"], permutations=999)
perm3[["aov.tab"]]

# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
# wk5_dat[-which(wk5_dat$Treatment == "cART"), "Treatment"]  1   0.15892 0.15892  1.4013 0.1229  0.226
# Residuals                                                 10   1.13411 0.11341         0.8771       
# Total                                                     11   1.29302                 1.0000       


#### week -1 ####

# extract just Week -1 OTU data for testing with permanova
wkpre_dat <- df %>%
  filter(Week == "-1")

#cART vs. placebo 
perm1 <- adonis(wkpre_dat[-which(wkpre_dat$Treatment == "Control"), 1:(length(wkpre_dat)-4)] ~
                  wkpre_dat[-which(wkpre_dat$Treatment=="Control"), "Treatment"], permutations=999)
perm1[["aov.tab"]] 

# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
# wkpre_dat[-which(wkpre_dat$Treatment == "Control"), "Treatment"]  1   0.10107 0.10107 0.69864 0.0653  0.774
# Residuals                                                        10   1.44660 0.14466         0.9347       
# Total                                                            11   1.54767                 1.0000       

#cART vs. Control
perm2 <- adonis(wkpre_dat[-which(wkpre_dat$Treatment=="Placebo"), 1:(length(wkpre_dat)-4)]~
                  wkpre_dat[-which(wkpre_dat$Treatment=="Placebo"), "Treatment"], permutations=999)
perm2[["aov.tab"]] 

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# wkpre_dat[-which(wkpre_dat$Treatment == "Placebo"), "Treatment"]  1   0.10332 0.10332 0.77944 0.07231  0.668
# Residuals                                                        10   1.32550 0.13255         0.92769       
# Total                                                            11   1.42882                 1.00000       

#placebo vs. Control

perm3 <- adonis(wkpre_dat[-which(wkpre_dat$Treatment=="cART"), 1:(length(wkpre_dat)-4)]~
                  wkpre_dat[-which(wkpre_dat$Treatment=="cART"), "Treatment"], permutations=999)
perm3[["aov.tab"]]

# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# wkpre_dat[-which(wkpre_dat$Treatment == "cART"), "Treatment"]  1   0.06500 0.065002 0.79292 0.07347  0.628
# Residuals                                                     10   0.81979 0.081979         0.92653       
# Total                                                         11   0.88479                  1.00000       

