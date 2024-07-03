### Alpha diversity analysis for FIV 16S study ###
# code modified from Dr. Zaid Abdo

################# startup #####################
library(tidyverse)
library(vegan)
library(phyloseq)
library(rstatix)
library(openxlsx)
library(ggpubr)

## set working directory

#### read in data
data.df <- read.csv("03_output/data.df.csv", header=T, row.names=1)
meta.df <- read.csv("03_output/meta.df.csv", header=T, row.names=1) %>%
  mutate(Week = as.factor(Week)) 
taxa.df <- read.csv("03_output/taxa.df.csv", header=T, row.names=1)

# expand metadata to include additional FIV study data like oral scores
meta.plus <- read.csv("01_input/dat.csv", header=T, row.names=1) %>%
  mutate(Week = as.factor(Week))

# make phyloseq object
data.ps = phyloseq(otu_table(data.df,taxa_are_rows = FALSE),tax_table(as.matrix(taxa.df)),sample_data(meta.df))

################# data prep #####################

#### We need row sums (sequencing depth) to calculate the rarefied richness below
depth = apply(data.df,1,sum)

######## Diversity analysis (old fashioned ecological approach [richness and diversity])
######## Done on the otu level!
#### We use estimate_richness() function in phyloseq to estimate diversity measures (here on non normalized data)
d.rich = estimate_richness(data.ps, measures=c("Observed", "Shannon", "InvSimpson"))

#### We use rarefy() function in vegan to normalize and calculate the expected richness based on the minimum
#### sequencing depth, and also compute the standard error se=TRUE (rarefy(data.df,min(d.rs),se=TRUE)) 
#### then we transform that calculated richness using t() (switching the rows and columns) 
rich.mt = t(rarefy(data.df,min(depth),se=TRUE))
#### We create a new data frame that incorporates columns 1, 6 and 8 from the richness estimates from phyloseq 
#### which are the observed richness, shannon diversity in InvSimpson. We add to that the first column (the expected richness)
#### from the vegan output and put it rich.df
rich.df = data.frame(cbind(d.rich,rich.mt[,1]))
#### We rename the columns to refelct what they are
names(rich.df) = c("Observed","Shannon","InvSimpson","VeganRichness")


################# tables ####################

# make a table with all diversity measures
rich_table <- rich.df %>%
  rownames_to_column("rowname") %>%
  left_join(rownames_to_column(meta.df), sto.info, by="rowname") %>%
  select(Week, Treatment,Cat_ID, Observed,Shannon,InvSimpson,VeganRichness) %>%
  rename("Observed richness" = Observed,
         "Shannon index" = Shannon, "Inverse Simpson index" = InvSimpson,
         "Normalized richness" = VeganRichness,
         "Cat" = Cat_ID) %>%
  relocate(Cat)

write.xlsx(rich_table, "03_output/FIV16S Richness Measures Table 20240208.xlsx", rowNames=T)

## Combine richness data with metadata for future statistical analysis

######## metadata for Scott July 2023 ##########
#### We will be testing richness using analysis of variance approaches
## First we combine the diversity data frame (rich.df) with the experimental design (meta.df) 
## and we add the sequencing depth to the mix (we add columns together using cbind())
meta.rich <- cbind(meta.df, rich.df, depth)

# Combine richness measures with all metadata for Scott
meta_all <- left_join(meta.plus, meta.rich, by = c("Sample", "Week", "Treatment")) %>%
  rename("SequencingDepth" = depth,
         "ObservedRichness" = Observed) 

write.csv(meta_all, "03_output/dat_20230706.csv")


# make a table with mean or median values compared across treatments
table1_mean <- meta_all %>%
  select(Sample, Week, Treatment, ObservedRichness, VeganRichness, Shannon, InvSimpson) %>%
  filter(Week %in% c("0", "4", "12", "25")) %>%
  group_by(Week, Treatment) %>%
  summarize(meanObs = mean(ObservedRichness, na.rm=T),
            meanVeg = mean(VeganRichness, na.rm=T),
            meanShan = mean(Shannon, na.rm=T),
            meanInv = mean(InvSimpson, na.rm=T)) %>% 
  gather(measure, value, meanObs:meanInv) %>%
  pivot_wider(names_from="Week", values_from="value") %>%
  relocate(measure)

table1_median <- meta_all %>%
  select(Sample, Week, Treatment, ObservedRichness, VeganRichness, Shannon, InvSimpson) %>%
  filter(Week %in% c("0", "4", "12", "25")) %>%
  group_by(Week, Treatment) %>%
  summarize(medianObs = median(ObservedRichness, na.rm=T),
            medianVeg = median(VeganRichness, na.rm=T),
            medianShan = median(Shannon, na.rm=T),
            medianInv = median(InvSimpson, na.rm=T)) %>% 
  gather(measure, value, medianObs:medianInv) %>%
  pivot_wider(names_from="Week", values_from="value") %>%
  relocate(measure)

#write.xlsx(table1_mean, "03_output/FIV16S Mean Diversity Measures Table 20230626.xlsx", rowNames=T)
#write.xlsx(table1_median, "03_output/FIV16S Median Diversity Measures Table 20230621.xlsx", rowNames=T)


####### more prep for plotting #########
#### This next step (making meta.long) is in preparation to plot. We can't plot the above measures all in the same graph
#### unless we can have all the values in one column and partition by Treatment level, and diversity measure.
#### We first create an empty container rich.plt = c()
meta.long = c()
#### We extract the column names from the rich.df data frame and put them in rich.nm
rich.nm = colnames(rich.df)
#### We have for columns and hence 4 names so we iterate from 1 to 4 (length(rich.nm))
for(i in 1:length(rich.nm)){
  #### We creat a variable rich that includes the name of a diversity measure every iteration
  #### repeated as many as the sample number (same as length of column of rich.df [length(rich.df[,1])])
  rich = rep(rich.nm[i],length(rich.df[,1]))
  #### We combine the columns of the meta.df, the new rich and the diversity measure we using cbind() function
  #### We then add the combined matrix (data frame realy) to the end of the new data frame rich.plt using rbind
  meta.long = rbind(meta.long,cbind(meta.df,rich,rich.df[,i]))
}
#### We rename the new columns to reflect what they are 
names(meta.long)[length(names(meta.long))] <- "value"
names(meta.long)[(length(names(meta.long))-1)] <- "richness_measure"

## clean up the environment
rm(d.rich, rich.df, rich, rich.mt, i, rich.nm, depth, rich_table_for_ms, meta.plus, meta.rich)


# make the meta.rich dataframe to plot from
meta.rich <- meta_all %>%
  filter(Week %in% c("0", "4", "12", "17", "25", "34")) %>%
  mutate(Treatment=factor(Treatment, levels= c("Control","Placebo", "cART")),
         Sample = factor(Sample, c(4542, 4547, 4550, 4535, 4537, 4540, #ctrl
                                   4541, 4546, 4545, 4534, 4543, 4548, #placebo
                                   4539, 4552, 4544, 4536, 4553, 4551)))
  

# make labels for weeks
wk.labs <- c("Week -1", "Week 5", "Week 11", "Week 24")
names(wk.labs) <- c("0", "4", "12", "25")

rich.labs <- c("Inverse simpson index", "Observed richness", "Shannon index", "Normalized richness")
names(rich.labs) <- c("InvSimpson", "Observed", "Shannon","VeganRichness")

################ quick visuals of the data #############
# overall increasing richness over time
ggplot(meta.rich) + 
  geom_point(aes(x= Week, y= VeganRichness, color=Treatment))

ggplot(meta.rich) + 
  geom_point(aes(x= Week, y= GingInf, color=Treatment))

ggplot(meta.rich) + 
  geom_point(aes(x= Week, y= GumLoss, color=Treatment))

# trajectory of richness for each cat
p <- ggplot(meta.rich %>% filter(!Week %in% c(17,34))) + 
  geom_point(aes(x= Week, y= VeganRichness, 
                 color=factor(Treatment, 
                              levels = c("Control", "Placebo", "cART"),
                              labels = c("Control", "Placebo", "cART")))) +
  scale_color_manual(values= c("black","blue","red"))+
  labs(y = "Normalized Richness", x= "Week", color = "Treatment")+
  facet_wrap(~Sample, ncol=6) +
  theme_bw()

pdf("03_output/figures/Richness_trajectory_by_cat.pdf", width=8, height=6)
p
dev.off()


######### 1: 0 weeks/prebleed diversity by treatment ################

week0 <- meta.rich %>% filter(Week == 0)

# check for normal distribution p> 0.05
week0 %>%
  group_by(Treatment) %>%
  shapiro_test(VeganRichness)

# check for outliers
week0 %>%
  group_by(Treatment) %>%
  identify_outliers(VeganRichness)

rich0 <- lm(VeganRichness ~ Treatment, data = week0)
anova(rich0)
summary(rich0) # overall p-value p=0.559, R-squared = 0.0746
# FIV cART vs healthy ctrl p=0.916
# FIV placebo vs healthy ctrl p=0.333

p0 <- ggplot(week0, aes(x=Treatment, y = VeganRichness)) +
  geom_boxplot()+
  theme_bw()

p0


######### 2: 4 weeks diversity by treatment (result: no differences yet) ################

week4 <- meta.rich %>% filter(Week == 4)

# check for normal distribution p> 0.05
week4 %>%
  group_by(Treatment) %>%
  shapiro_test(VeganRichness)

# check for outliers (one but not extreme)
week4 %>%
  group_by(Treatment) %>%
  identify_outliers(VeganRichness) %>%
  select(c(1:3, VeganRichness, is.outlier, is.extreme))

rich4 = lm(VeganRichness ~ Treatment, data = week4)
anova(rich4)
summary(rich4) 
# FIV cART vs healthy ctrl p=0.390
# FIV placebo vs healthy ctrl p=0.358
# overall p-value p=0.5818, R-squared = 0.0697

p4 <- ggplot(week4, aes(x=Treatment, y = VeganRichness)) +
  geom_boxplot()+
  theme_bw()

p4


######### 3: 12 weeks diversity by treatment (result: still no differences) ################

week12 <- meta.rich %>% filter(Week == 12)

# check for normal distribution p > 0.05
week12 %>%
  group_by(Treatment) %>%
  shapiro_test(VeganRichness) # uh oh! Placebo data is not normally distributed

# check for outliers -- none
week12 %>%
  group_by(Treatment) %>%
  identify_outliers(VeganRichness) %>%
  select(c(1:3, VeganRichness, is.outlier, is.extreme))

rich12 = lm(VeganRichness ~ Treatment, data = week12)
anova(rich12)
summary(rich12) 
# FIV cART vs healthy ctrl p=0.565
# FIV placebo vs healthy ctrl p=0.371
# overall p-value p=0.3399, R-squared = 0.134

p12 <- ggplot(week12, aes(x=Treatment, y = VeganRichness)) +
  geom_boxplot()+
  theme_bw()

p12


######### 4: 25 weeks (result: ) ################

week25 <- meta.rich %>% filter(Week == 25)

# check for normal distribution p > 0.05
week25 %>%
  group_by(Treatment) %>%
  shapiro_test(VeganRichness) # uh oh! control and Placebo data is not normally distributed

# check for outliers -- two, one extreme in control group
week25 %>%
  group_by(Treatment) %>%
  identify_outliers(VeganRichness) %>%
  select(c(1:3, VeganRichness, is.outlier, is.extreme))

rich25 = lm(VeganRichness ~ Treatment, data = week25)
anova(rich25)
summary(rich25) 
# FIV cART vs healthy ctrl p=0.784
# FIV placebo vs healthy ctrl p=0.345
# overall p-value p=0.6143, R-squared = 0.0629

p25 <- ggplot(week25, aes(x=Treatment, y = VeganRichness)) +
  geom_boxplot()+
  theme_bw()

p25

############## each week diversity #############
# pdf("03_output/figures/norm_richness_by_trtwk.pdf", width=5, height=4)
# ggplot(meta.rich %>% filter(!Week %in% c(17, 34)) %>%
#          mutate(Treatment=factor(Treatment, levels= c("Control","Placebo", "cART")))) +
#   geom_boxplot(aes(y = VeganRichness, x = Treatment)) +
#   geom_jitter(aes(y = VeganRichness, x = Treatment, color=Treatment), size=1.5, alpha=0.7) +
#   facet_wrap(~Week, labeller = labeller(Week = wk.labs))+
#   #scale_color_manual(values= c("#1b9e77","#7570b3","#d95f02"))+
#   scale_color_manual(values= c("black","blue","red"))+
#   theme_bw() +
#   labs(x="Treatment", y= "Normalized Richness")+
#   theme(axis.text = element_text(size=12),
#         strip.text = element_text(size=12),
#         axis.title = element_text(size=14),
#         legend.position = "none")
# dev.off()

rich1 <- ggplot(meta.rich %>% filter(!Week %in% c(17, 34)) %>%
         mutate(Treatment=factor(Treatment, 
                                 levels= c("Control","Placebo", "cART"))),
         aes(y = VeganRichness, x = Week)) +
  geom_boxplot() +
  geom_point(aes(color=Treatment), size=1.5, alpha=0.7) +
  # geom_jitter(aes(color=Treatment), size=1.5, alpha=0.7,
  #             width=0.1, height=0.1) +
  facet_wrap(~Treatment)+
  #scale_color_manual(values= c("#1b9e77","#7570b3","#d95f02"))+
  scale_color_manual(values= c("black","blue","red"))+
  scale_x_discrete(labels=c("0" = "-1", "4" = "5",
                              "12" = "11", "25" = "24")) + 
  theme_bw() +
  labs(x="Treatment", y= "Normalized Richness")+
  theme(axis.text = element_text(size=12),
        strip.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.position = "none")

#wk.labs <- c("Week -1", "Week 5", "Week 11", "Week 24")

rich1

pdf("03_output/figures/norm_richness_by_trtwk_20230630.pdf", width=6, height=4)
rich1
dev.off()

meta.long <- meta.long %>%
  filter(Week %in% c("0", "4", "12", "25"))

# pdf("03_output/figures/all_richness_by_trtwk.pdf")
# ggplot(meta.long) +
#   geom_boxplot(aes(y = value, x = Treatment, color=richness_measure)) +
#   geom_jitter(aes(y = value, x = Treatment), size=1, alpha=0.5) +
#   facet_grid(richness_measure~Week, labeller = labeller(Week = wk.labs,
#                                                         richness_measure = rich.labs), scales="free")+
#   scale_color_manual(values= c("#7570b3","#1b9e77","#d95f02","#e7298a"))+
#   theme_bw() +
#   labs(color="Richness measure", x="Treatment")+
#   theme(axis.text = element_text(size=12),
#         strip.text = element_text(size=12),
#         axis.title = element_text(size=14),
#         legend.position = "none")
# dev.off()

############## gingival inf and gum loss related to diversity #############

gingdiv <- lm(VeganRichness ~ GingInf, data = meta.rich)
summary(gingdiv) #not sig R-squared:  0.002663, p=0.6668

gumdiv <- lm(VeganRichness ~ GumLoss, data = meta.rich)
summary(gumdiv) #not sig R-squared:  0.00545, p = 0.5377

p_inf <- ggplot(meta.rich, aes(x=VeganRichness, y = GingInf)) +
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

p_inf


############## diversity over time #########

## look at the richness over study week/cat age

rich.wk <- lm(VeganRichness ~ Week, data = meta.rich)
anova(rich.wk)
summary(rich.wk) #very significant F(3,68) = 4.753, R2 = 0.173, p = 0.00455

# other richness measures?
rich.wk2 <- lm(ObservedRichness ~ Week, data = meta.rich)
summary(rich.wk2) # sig F(3,68) = 4.762, R2 = 0.1736, p = 0.004509

rich.wk3 <- lm(Shannon ~ Week, data = meta.rich)
summary(rich.wk3) # sig F(3,68) = 6.192, R2 = 0.2146, p-value: 0.0008779

rich.wk4 <- lm(InvSimpson ~ Week, data = meta.rich)
summary(rich.wk4) # borderlien F(3,68) = 2.586, R2 = 0.1024, p=0.06016


# age
rich.age <- lm(VeganRichness ~ Age, data = meta.rich)
summary(rich.age) #also significant p=0.01329

# other richness measures?
rich.age2 <- lm(ObservedRichness ~ Age, data = meta.rich)
summary(rich.age2) # sig p=0.01596

rich.age3 <- lm(Shannon ~ Age, data = meta.rich)
summary(rich.age3) # sig p-value: 0.0001973

rich.age4 <- lm(InvSimpson ~ Age, data = meta.rich)
summary(rich.age4) # sig p=0.009382



p_age <- ggplot(meta.rich, aes(x=(Age/30), 
                               y = VeganRichness))+
  geom_point(alpha=0.7) +
  geom_smooth(aes(color=Treatment), method="lm", se=F) +
  labs(y= "Normalized richness", x= "Cat age (months)") +
  scale_x_continuous(limits= c(5.5,15.5), breaks=seq(6,14,2)) +
  scale_color_manual(values= c("black","blue","red"))+
  theme_bw() + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14))

p_age

pdf("03_output/figures/cat_age_diversity.pdf", width=4, height=4)
p_age
dev.off()

ggarrange(rich1, p_age)


#### GingInf and GumLoss at Week 34 #####

dental_plot <- meta.rich %>%
  filter(Week == 34) %>%
  select(c("Sample", "Treatment", "GingInf", "GumLoss")) %>%
  pivot_longer(!c("Sample", "Treatment"),
               names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, c("GingInf", "GumLoss"), 
                          labels=c("Week 35 gingival inflammation", "Week 35 attachment loss")))

p.dental <- ggplot(dental_plot, 
                   aes(x=Treatment, y=value)) +
  geom_boxplot() +
  geom_jitter(aes(color=Treatment), size=1.5, alpha=0.7) +
  #geom_point(aes(color=Treatment)) +
  facet_wrap(~measure) +
  scale_color_manual(values= c("#1b9e77","#7570b3","#d95f02"))+
  labs(x="", y="") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        strip.text = element_text(size=12),
        axis.title = element_text(size=14))

p.dental

pdf("03_output/figures/Oral_disease_Week34.pdf", width=6, height=4)
p.dental
dev.off()
