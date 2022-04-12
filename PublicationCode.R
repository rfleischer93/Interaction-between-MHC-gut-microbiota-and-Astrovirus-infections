
## Analysis R script

## Script for processing A. jamaicensis MHC and microbiome data for the publication
## Interaction between MHC diversity and constitution, gut microbiota and Astrovirus infections in a neotropical bat
## Ramona Fleischer, Dominik Werner Schmid, Wasimuddin, Stefan Dominik Brändel, Andrea Rasche, Victor M. Corman, Christian Drosten, Marco Tschapka, Simone Sommer
## Ulm University, Institute of Evolutionary Ecology and Conservation Genomics


# load library
library(tibble)
library(tidyverse)
library(ggplot2)
library(phyloseq)
library(vegan)
library(data.table)
library(lme4)
library(ggreadClipboard())
library(tidyr)
library(ggthemes)
library(microbiome)
#library(plyr)
library(dplyr)
library(brglm)
library(effects)
library(magrittr)
library(MuMIn)
library(ggrepel)
library(car)
library(ggpubr)

# load theme file for ggplot
source("my_themes.R")



############################################################ MHC data ##############################################################
############################################################ MHC data ##############################################################
############################################################ MHC data ##############################################################

# load MHC metadata

arja_mhc <- read.csv("mhc_meta.csv")
arja_mhc # the 224 bats with joined MHC + GMB info after rarefying


# dataset overview

# Nr of individuals with / without AstV

table(arja_mhc$astro)
# NEG  POS
# 152  72


# Nr of individuals with / without AstV separated by age

table(arja_mhc$astro, arja_mhc$age)
#      Adult  Young
# NEG  91     61
# POS  40     32


# mean NrAlleles & NrST in the population

mean(arja_mhc$Nr_of_alleles_aa)
# 3.044643
sd(arja_mhc$Nr_of_alleles_aa)
# 1.123642

mean(arja_mhc$Nr_ST)
# 2.584821
sd(arja_mhc$Nr_ST)
# 0.8843137




######## MHC allele frequency ############

# make metadata into long format

mhc_long<-gather(arja_mhc, allele, freq, Arja.DRB.01:Arja.DRB.134, factor_key=TRUE) # for alleles
mhc_long$freq<-factor(mhc_long$freq)

# subset to plot allele presence
mhc_long_allele<-subset(mhc_long, freq != 0) 

# plot
allele_freq <- ggplot(mhc_long_allele,aes(x = allele, fill=freq))+
  geom_histogram(stat="count")+
  labs(x=" ",y="Allele frequency")+
  scale_fill_manual(values=c("grey45"))+ 
  geom_hline(aes(yintercept = 4), linetype = "dashed") +  # dahsed line at 4 indicating all alleles above this line are present in at least 5 individuals
  theme_bw(base_size=18)+
  theme(legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text.x = element_text(size = 8))


########## MHC ST frequency ##############

# make metadata into long format

mhc_long_st<-gather(arja_mhc, supertype, freq, ST1:ST10, factor_key=TRUE) # for supertypes
mhc_long_st$freq<-factor(mhc_long_st$freq)

# subset to plot ST presence
mhc_long_st<-subset(mhc_long_st, freq != 0)

# plot
st_freq<-ggplot(mhc_long_st,aes(x = supertype, fill=freq))+
  geom_histogram(stat="count")+
  labs(x=" ",y="Supertype frequency")+
  scale_fill_manual(values=c("grey45"))+ 
  geom_hline(aes(yintercept = 9), linetype = "dashed") +  # dahsed line at 9 indicating all ST above this line are present in at least 10 individuals
  theme_bw(base_size=18)+
  theme(legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(size = 14))+ 
  labs(y ="Supertype frequency", x=" ")


##########################################
######## Supplementary Figure S3 #########
##########################################

ggarrange(allele_freq, st_freq,labels=c("A","B"), nrow=2,ncol=1)




###########################################################################
######################## co-occurrence analysis ###########################
######################## co-occurrence analysis ###########################
######################## co-occurrence analysis ###########################
###########################################################################

library(cooccur)

############## Alleles ############## 
############## Alleles ############## 

# co-occurrence analysis of common mhc alleles and astrovirus
# we define common alleles as all mhc alleles that ocurr in at least 5 individuals

# first select all allele and virus info we want to investigate for co-occurrence associations

alleles.virus1<-subset(arja_mhc, select=c("Arja.DRB.01",
                                             "Arja.DRB.02",
                                             "Arja.DRB.03",
                                             "Arja.DRB.04",
                                             "Arja.DRB.05",
                                             "Arja.DRB.06",
                                             "Arja.DRB.07",
                                             "Arja.DRB.08",
                                             "Arja.DRB.09",
                                             "Arja.DRB.10",
                                             "Arja.DRB.11",
                                             "Arja.DRB.12",
                                             "Arja.DRB.13",
                                             "Arja.DRB.14",
                                             "Arja.DRB.15",
                                             "Arja.DRB.16",
                                             "Arja.DRB.17",
                                             "Arja.DRB.18",
                                             "Arja.DRB.19",
                                             "Arja.DRB.20",
                                             "Arja.DRB.21",
                                             "Arja.DRB.22",
                                             "Arja.DRB.23",
                                             "Arja.DRB.24",
                                             "Arja.DRB.25",
                                             "Arja.DRB.26",
                                             "Arja.DRB.27",
                                             "Arja.DRB.28",
                                             "Arja.DRB.29",
                                             "Arja.DRB.30",
                                             "Arja.DRB.31",
                                             "Arja.DRB.32",
                                             "Arja.DRB.33",
                                             "Arja.DRB.34",
                                             "Arja.DRB.35",
                                             "Arja.DRB.36",
                                             "Arja.DRB.37",
                                             "Arja.DRB.38",
                                             "Arja.DRB.39",
                                             "Arja.DRB.40",
                                             "Arja.DRB.41",
                                             "Arja.DRB.42",
                                             "Arja.DRB.43",
                                             "Arja.DRB.45",
                                             "Arja.DRB.46",
                                             "Arja.DRB.47",
                                             "Arja.DRB.48",
                                             "Arja.DRB.49",
                                             "Arja.DRB.50",
                                             "Arja.DRB.51",
                                             "Arja.DRB.52",
                                             "Arja.DRB.53",
                                             "Arja.DRB.54",
                                             "Arja.DRB.55",
                                             "Arja.DRB.58",
                                             "Arja.DRB.60",
                                             "Arja.DRB.62",
                                             "Arja.DRB.63",
                                             "astro2"))

str(alleles.virus1) # no NA-values for alleles, fine

# create matrix and transpose

alleles.virus1_long <- t(as.matrix(alleles.virus1))

# apply cooccurrence
cooccur.aj1 <- cooccur(mat = alleles.virus1_long, type = "spp_site", thresh = FALSE, spp_names = TRUE)

summary(cooccur.aj1) # short summary
plot(cooccur.aj1)    # square plot

# which alleles have a positive or negative association with astrovirus?

pair(mod = cooccur.aj1, "astro2") 
# alleles 13, 31, 42, 55 are positively associated with AstV presence
# allele 58 negatively associated with AstV presence

prob.table(cooccur.aj1) # table of expected & observed cooccurrence and p-values

#effect.sizes(cooccur.aj1, standardized = FALSE, matrix = FALSE)                    # raw effect sizes
st.effSize.alleles<-effect.sizes(cooccur.aj1, standardized = TRUE, matrix = FALSE)  # standardized effect sizes
#write.csv(st.effSize.alleles,"st.effSize.alleles.csv")

# subset only AstV-associated results
res <- cooccur.aj1$results
res1 <- subset(res, sp2_name=="astro2")


# plot the expected and observed co-occurrence in a barplot

# rename some variables for the plot
colnames(res1)[5] <- "Observed"
colnames(res1)[7] <- "Expected"

# make into long format
res_long_allele<-gather(res1, probability, value, Expected, Observed, factor_key = TRUE)
#write.csv(res_long_allele, "res_long_allele.csv") # save

# mark significant results by adding a new column "sign" to the df indicating the significance with "*" / "**" / "***"
# also change allele names from "Arja.DRB.XX" to "*XX" for plotting

res_long_allele <- read.csv("res_long_allele.csv")


#################################
####### Final Figure 2A #########
#################################

cooccurr.alleles <- ggplot(res_long_allele, aes(y=value, x=sp1_name, fill = probability, label = sign))+ 
  theme_classic(base_size = 14)+
  geom_bar(position="dodge", stat="identity") +
  geom_text(vjust = - 0.3, hjust = 0.5, size = 7) +
  scale_fill_manual(values=c("grey","deepskyblue4"))+
  labs(x="",y="AstV co-occurrence",fill="AstV co-occurrence probability") +
  theme(axis.text.x = element_text(size = 11),axis.title.y = element_text(size = 15))



############# Supertypes #############
############# Supertypes #############

# first select all ST and virus info we want to investigate for co-occurrence 

st.virus<-subset(arja_mhc, select=c("ST1","ST2","ST3","ST4","ST5","ST6","ST7","ST8","ST9","ST10","astro2"))
str(st.virus) # no NA-values for alleles, fine

#create matrix and transpose
st.virus_long <- t(as.matrix(st.virus))

# apply cooccurrence
cooccur.aj <- cooccur(mat = st.virus_long, type = "spp_site", thresh = FALSE, spp_names = TRUE)

summary(cooccur.aj) # short summary
plot(cooccur.aj)    # square plot

# which alleles have a positive or negative association with astrovirus

pair(mod = cooccur.aj, "astro2")
# ST5 is positively associated with AstV presence

prob.table(cooccur.aj) # table of expected & observed cooccurrence and p-values

#effect.sizes(cooccur.aj, standardized = FALSE, matrix = FALSE)                        # raw effect sizes
st.effSize.supertypes<-effect.sizes(cooccur.aj, standardized = TRUE, matrix = FALSE)  # standardized effect sizes
#write.csv(st.effSize.supertypes,"st.effSize.supertypes.csv")

# subset only AstV variable
res <- cooccur.aj$results
res1 <- subset(res, sp2_name=="astro2")


# plot the expected and observed co-occurrence in a barplot

# rename some variables for the plot
colnames(res1)[5] <- "Observed"
colnames(res1)[7] <- "Expected"

# make into long format
res_long_st<-gather(res1, probability, value, Expected, Observed, factor_key = TRUE)
#write.csv(res_long_st, "res_long_st.csv")

res_long_st <- read.csv("res_long_st.csv")
res_long_st$sp1_name <- factor(res_long_st$sp1_name, levels = c("astro2", "ST1","ST2","ST3","ST4","ST5","ST6","ST7","ST8","ST9","ST10")) # change order of supertypes


#################################
####### Final Figure 2B #########
#################################

cooccurr.st <- ggplot(res_long_st, aes(y=value, x=sp1_name, fill = probability, label = sign))+ 
  theme_classic(base_size = 14)+
  geom_bar(position="dodge", stat="identity") +
  geom_text(vjust = - 0.1, hjust = 0.5, size = 7) +
  scale_fill_manual(values=c("grey","deepskyblue4"))+
  labs(x="",y="AstV co-occurrence",fill="AstV co-occurrence probability")+
  theme(axis.text.x = element_text(size = 11),axis.title.y = element_text(size = 15))



#################################
######## Final Figure 2 #########
#################################

Figure2AB <- ggarrange(cooccurr.alleles, nrow=2, cooccurr.st, common.legend = TRUE, legend = "bottom", labels = c("A", "B"))





# Do individual MHC diversity estimates differ according to AstV infection status, age, sex or landscape?


############################################
######### Supplementary Table S3 ###########
############################################

# Nr_Alleles
kruskal.test(Nr_of_alleles_aa ~ astro, data = arja_mhc)
kruskal.test(Nr_of_alleles_aa ~ age, data = arja_mhc)
kruskal.test(Nr_of_alleles_aa ~ sex, data = arja_mhc)
kruskal.test(Nr_of_alleles_aa ~ landscape, data = arja_mhc)

# Nr_ST
kruskal.test(Nr_ST ~ astro, data = arja_mhc)
kruskal.test(Nr_ST ~ age, data = arja_mhc)
kruskal.test(Nr_ST ~ sex, data = arja_mhc)
kruskal.test(Nr_ST ~ landscape, data = arja_mhc)

# Mean p-distance among alleles
kruskal.test(AA_pdist_mean ~ astro, data = arja_mhc)
kruskal.test(AA_pdist_mean ~ age, data = arja_mhc)
kruskal.test(AA_pdist_mean ~ sex, data = arja_mhc)
kruskal.test(AA_pdist_mean ~ landscape, data = arja_mhc)

# Sum of p-distances among alleles
kruskal.test(AA_pdist_sum ~ astro, data = arja_mhc)
kruskal.test(AA_pdist_sum ~ age, data = arja_mhc)
kruskal.test(AA_pdist_sum ~ sex, data = arja_mhc)
kruskal.test(AA_pdist_sum ~ landscape, data = arja_mhc)

# Mean of p-distances among PSS
kruskal.test(AA_pdistPSS_mean ~ astro, data = arja_mhc)
kruskal.test(AA_pdistPSS_mean ~ age, data = arja_mhc)
kruskal.test(AA_pdistPSS_mean ~ sex, data = arja_mhc)
kruskal.test(AA_pdistPSS_mean ~ landscape, data = arja_mhc)

# Sum of p-distances among PSS
kruskal.test(AA_pdistPSS_sum ~ astro, data = arja_mhc)
kruskal.test(AA_pdistPSS_sum ~ age, data = arja_mhc)
kruskal.test(AA_pdistPSS_sum ~ sex, data = arja_mhc)
kruskal.test(AA_pdistPSS_sum ~ landscape, data = arja_mhc)

# MHC Faith's PD among alleles
kruskal.test(mhc_faiths ~ astro, data = arja_mhc)
kruskal.test(mhc_faiths ~ age, data = arja_mhc)
kruskal.test(mhc_faiths ~ sex, data = arja_mhc)
kruskal.test(mhc_faiths ~ landscape, data = arja_mhc)





# Testing the effect on AstV infection status


m1.1<-glm(astro2 ~ age + sex + mhc_faiths + Nr_ST + Nr_of_alleles_aa + AA_pdist_mean + AA_pdist_sum + AA_pdistPSS_mean + AA_pdistPSS_sum + AA_poissondist_mean, data = arja_mhc, na.action=na.fail, family=binomial)
summary(m1.1)

# many variables
# test for multicolinearity with ggally 

#install.packages("GGally")
library(GGally)

# extract variables to investigate for colinearity
exp_v_mhc<-arja_mhc %>% 
  select(Nr_ST,Nr_of_alleles_aa,AA_pdist_mean,AA_pdist_sum,AA_pdistPSS_mean,AA_pdistPSS_sum,mhc_faiths)

# check
mhc_colinearity <- ggpairs(exp_v_mhc)
mhc_colinearity # plot

corr <- ggcorr(exp_v_mhc, method = c("pairwise", "pearson")) 
corr # plot
corr$data # pearson coefficients

# A lot of MHC diversity estimates are co linear, so better include every estimate separately
# include landscape as random factor



############################################
######### Supplementary Table S4 ###########
############################################

m<-glmer(astro2 ~ age + sex + Nr_of_alleles_aa + (1|landscape), data = arja_mhc, na.action=na.fail, family=binomial)
summary(m)

m1<-glmer(astro2 ~ age + sex + Nr_ST + (1|landscape), data = arja_mhc, na.action=na.fail, family=binomial)
summary(m1)

m2<-glmer(astro2 ~ age + sex +  AA_pdist_mean + (1|landscape), data = arja_mhc, na.action=na.fail, family=binomial)
summary(m2)

m3<-glmer(astro2 ~ age + sex +  AA_pdist_sum + (1|landscape), data = arja_mhc, na.action=na.fail, family=binomial)
summary(m3)

m4<-glmer(astro2 ~ age + sex +  AA_pdistPSS_mean + (1|landscape), data = arja_mhc, na.action=na.fail, family=binomial)
summary(m4)

m5<-glmer(astro2 ~ age + sex +  AA_pdistPSS_sum + (1|landscape), data = arja_mhc, na.action=na.fail, family=binomial)
summary(m5)

m6<-glmer(astro2 ~ age + sex +  mhc_faiths + (1|landscape), data = arja_mhc, na.action=na.fail, family=binomial)
summary(m6)




##############################
### Supplementary Table S5 ###
##############################

# Testing potential associations between below/above average MHC diversity and AstV infection 
# using Pearson's Chi-squared test with Yates' continuity correction (X2)

arja_mhc_y <- subset(arja_mhc, age == "Y")
arja_mhc_a <- subset(arja_mhc, age == "A")


################ Adult

print(chisq.test(arja_mhc_a$pdist_mean_binary,arja_mhc_a$astro))
print(chisq.test(arja_mhc_a$pdist_sum_binary,arja_mhc_a$astro))
print(chisq.test(arja_mhc_a$pdistPSS_mean_binary,arja_mhc_a$astro))
print(chisq.test(arja_mhc_a$pdistPSS_sum_binary,arja_mhc_a$astro))

################ Young

print(chisq.test(arja_mhc_y$pdist_mean_binary,arja_mhc_y$astro))
print(chisq.test(arja_mhc_y$pdist_sum_binary,arja_mhc_y$astro))
print(chisq.test(arja_mhc_y$pdistPSS_mean_binary,arja_mhc_y$astro))
print(chisq.test(arja_mhc_y$pdistPSS_sum_binary,arja_mhc_y$astro))








######################################################## Microbiome data ###########################################################
######################################################## Microbiome data ###########################################################
######################################################## Microbiome data ###########################################################


####################################################################################################################################
######################################### Compare rarefied and unrarefied alpha diversity  #########################################
####################################################################################################################################

# library to calculate alpha-diversity

#install.packages("remotes")
#remotes::install_github("twbattaglia/btools")
library(remotes)
library(btools)
library(effects)


# load phyloseq datasets

###### unrarefied dataset ######

aj_unrare <- readRDS("Aj_class_unrare.RDS")


# calculate alpha diversity

aj_unrare_df<-as.data.frame(sample_data(aj_unrare))

richness<-estimate_richness(aj_unrare, measures=c("Observed","Shannon","Chao1"))
faiths_richness<-estimate_pd(aj_unrare)

aj_unrare_df$Observed<-richness$Observed
aj_unrare_df$Shannon<-richness$Shannon
aj_unrare_df$Chao1<-richness$Chao1
aj_unrare_df$Faiths<-faiths_richness$PD
names(aj_unrare_df)


####### rarefied dataset #######

aj_rare <- readRDS("Aj_class_rare11000.RDS")


# calculate alpha diversity

aj_rare_df<-as.data.frame(sample_data(aj_rare))

richness1<-estimate_richness(aj_rare, measures=c("Observed","Shannon","Chao1"))
faiths_richness1<-estimate_pd(aj_rare)

aj_rare_df$Observed_rare<-richness1$Observed
aj_rare_df$Shannon_rare<-richness1$Shannon
aj_rare_df$Chao1_rare<-richness1$Chao1
aj_rare_df$Faiths_rare<-faiths_richness1$PD
names(aj_rare_df)


# check dataframes

head(aj_unrare_df)
head(aj_rare_df)

tail(aj_unrare_df)
tail(aj_rare_df)


# both are in the same order
# now add the rarefied alpha diversity scores to the same df as the unrarefied so we can analyse them together

aj_unrare_df$Observed_rare<-aj_rare_df$Observed_rare
aj_unrare_df$Shannon_rare<-aj_rare_df$Shannon_rare
aj_unrare_df$Chao1_rare<-aj_rare_df$Chao1_rare
aj_unrare_df$Faiths_rare<-aj_rare_df$Faiths_rare

names(aj_unrare_df) # check it worked


# now compare the results from rarefied data to non-rarefied data

plot(aj_unrare_df$Observed~aj_unrare_df$Observed_rare)
plot(aj_unrare_df$Shannon~aj_unrare_df$Shannon_rare)
plot(aj_unrare_df$Chao1~aj_unrare_df$Chao1_rare)
plot(aj_unrare_df$Faiths~aj_unrare_df$Faiths_rare)


corr.obs <- ggplot(aj_unrare_df, aes(x = Observed, y = Observed_rare))+
  geom_point()+
  stat_smooth(method = "lm")+
  stat_cor()+
  theme_bw(base_size = 12)+
  labs(x= "Observed (unrarefied)", y = "Observed (rarefied)") +
  my_theme1

corr.sh <- ggplot(aj_unrare_df, aes(x = Shannon, y = Shannon_rare))+
  geom_point()+
  stat_smooth(method = "lm")+
  stat_cor()+
  theme_bw(base_size = 12)+
  labs(x= "Shannon (unrarefied)", y = "Shannon (rarefied)") +
  my_theme1

corr.ch <- ggplot(aj_unrare_df, aes(x = Chao1, y = Chao1_rare))+
  geom_point()+
  stat_smooth(method = "lm")+
  stat_cor()+
  theme_bw(base_size = 12)+
  labs(x= "Chao1 (unrarefied)", y = "Chao1 (rarefied)") +
  my_theme1

corr.f <- ggplot(aj_unrare_df, aes(x = Faiths, y = Faiths_rare))+
  geom_point()+
  stat_smooth(method = "lm")+
  stat_cor()+
  theme_bw(base_size = 12)+
  labs(x= "Faith's PD (unrarefied)", y = "Faith's PD (rarefied)") +
  my_theme1


#########################################
######## Supplementary Figure S1 ########
#########################################

ggarrange(corr.obs, corr.sh, corr.ch, corr.f, ncol=2, nrow=2)





####################################################################################################################################
######################################## Continue with UNRAREFIED alpha diversity analysis #########################################
####################################################################################################################################

aj_unrare

aj_df_final <- data.frame(sample_data(aj_unrare))


aj_df_final$ST1 <- as.factor(aj_df_final$ST1)
aj_df_final$ST2 <- as.factor(aj_df_final$ST2)
aj_df_final$ST3 <- as.factor(aj_df_final$ST3)
aj_df_final$ST4 <- as.factor(aj_df_final$ST4)
aj_df_final$ST5 <- as.factor(aj_df_final$ST5)
aj_df_final$ST6 <- as.factor(aj_df_final$ST6)
aj_df_final$ST7 <- as.factor(aj_df_final$ST7)
aj_df_final$ST8 <- as.factor(aj_df_final$ST8)
aj_df_final$ST9 <- as.factor(aj_df_final$ST9)
aj_df_final$ST10 <- as.factor(aj_df_final$ST10)


## separate astro-positive and negative bats for modelling

aj_df_final_pos<-subset(aj_df_final, astro == "POS")
aj_df_final_neg<-subset(aj_df_final, astro == "NEG")

## separate by age as well

aj_df_final_pos_a <- subset(aj_df_final_pos, age == "A")
aj_df_final_pos_y <- subset(aj_df_final_pos, age == "Y")

aj_df_final_neg_a <- subset(aj_df_final_neg, age == "A")
aj_df_final_neg_y <- subset(aj_df_final_neg, age == "Y")




####################################################################################################################################
####################################################### Astro Negatives ############################################################
####################################################################################################################################


#########################################
######## Supplementary Figure S4 ########
#########################################


aj_alpha_long<-gather(aj_df_final, alpha_measures, richness, Observed,Shannon,Chao1,Faiths, factor_key=TRUE) # for supertypes
aj_alpha_long$alpha_measures

aj_alpha_long_neg <- subset(aj_alpha_long, astro == "NEG")

age_label <- c("A" = "Adult", "Y" = "Young")


FigS4A <- ggplot(aj_alpha_long_neg, aes(x = alpha_measures, y = richness, fill=ST5, color= ST5))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(labels = c("absent","present"),values=c("#c29844","#6ab166"))+
  scale_fill_manual(labels = c("absent","present"),values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  scale_y_log10()+
  labs(x= " ", y =expression(Log[10] * " alpha diversity")) +
  theme_set(theme_bw()) + 
  my_theme4+ 
  theme(axis.title.y = element_text(size = 14))

FigS4B <- ggplot(aj_alpha_long_neg, aes(x = alpha_measures, y = richness, fill=ST6, color= ST6))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(labels = c("absent","present"),values=c("#c29844","#6ab166"))+
  scale_fill_manual(labels = c("absent","present"),values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  scale_y_log10()+
  labs(x= "Observed     Shannaon     Chao1     Faith's            Observed     Shannaon     Chao1     Faith's
        ASVs           Index            Index        PD                   ASVs           Index            Index          PD       ", y =expression(Log[10] * " alpha diversity")) +
  theme_set(theme_bw()) + 
  my_theme4+ 
  theme(axis.title.y = element_text(size = 14))


ggarrange(FigS4A,FigS4B, labels = c("A", "B"), nrow = 2, heights = c(1,1.05))



########################################
### Stats for Supplementary Table S6 ###
########################################

############## Adult AstV- #############


m1 <- lm(log(aj_df_final_neg_a$Observed) ~ aj_df_final_neg_a$SequencingDepth + aj_df_final_neg_a$ST5)
summary(m1)

m1 <- lm(aj_df_final_neg_a$Shannon ~ aj_df_final_neg_a$SequencingDepth + aj_df_final_neg_a$ST5)
summary(m1)

m1 <- lm(log(aj_df_final_neg_a$Chao1) ~ aj_df_final_neg_a$SequencingDepth + aj_df_final_neg_a$ST5)
summary(m1)

m1 <- lm(log(aj_df_final_neg_a$Faiths) ~ aj_df_final_neg_a$SequencingDepth + aj_df_final_neg_a$ST5)
summary(m1)


m1 <- lm(log(aj_df_final_neg_a$Observed) ~ aj_df_final_neg_a$SequencingDepth + aj_df_final_neg_a$ST6)
summary(m1)

m1 <- lm(aj_df_final_neg_a$Shannon ~ aj_df_final_neg_a$SequencingDepth + aj_df_final_neg_a$ST6)
summary(m1)

m1 <- lm(log(aj_df_final_neg_a$Chao1) ~ aj_df_final_neg_a$SequencingDepth + aj_df_final_neg_a$ST6)
summary(m1)

m1 <- lm(log(aj_df_final_neg_a$Faiths) ~ aj_df_final_neg_a$SequencingDepth + aj_df_final_neg_a$ST6)
summary(m1)



############# Young AstV- ############


m1 <- lm(log(aj_df_final_neg_y$Observed) ~ aj_df_final_neg_y$SequencingDepth + aj_df_final_neg_y$ST5)
summary(m1)

m1 <- lm(aj_df_final_neg_y$Shannon ~ aj_df_final_neg_y$SequencingDepth + aj_df_final_neg_y$ST5)
summary(m1)

m1 <- lm(aj_df_final_neg_y$Chao1 ~ aj_df_final_neg_y$SequencingDepth + aj_df_final_neg_y$ST5)
summary(m1)

m1 <- lm(aj_df_final_neg_y$Faiths ~ aj_df_final_neg_y$SequencingDepth + aj_df_final_neg_y$ST5)
summary(m1)


m1 <- lm(log(aj_df_final_neg_y$Observed) ~ aj_df_final_neg_y$SequencingDepth + aj_df_final_neg_y$ST6)
summary(m1)

m1 <- lm(aj_df_final_neg_y$Shannon ~ aj_df_final_neg_y$SequencingDepth + aj_df_final_neg_y$ST6)
summary(m1)

m1 <- lm(aj_df_final_neg_y$Chao1 ~ aj_df_final_neg_y$SequencingDepth + aj_df_final_neg_y$ST6)
summary(m1)

m1 <- lm(aj_df_final_neg_y$Faiths ~ aj_df_final_neg_y$SequencingDepth + aj_df_final_neg_y$ST6)
summary(m1)




########################################
### Stats for Supplementary Table S7 ###
########################################


############# Adult AstV- #############

m1 <- lm(log(aj_df_final_neg_a$Observed) ~ aj_df_final_neg_a$SequencingDepth + aj_df_final_neg_a$mhc_faiths_binary)
summary(m1)

m1 <- lm(aj_df_final_neg_a$Shannon ~ aj_df_final_neg_a$SequencingDepth + aj_df_final_neg_a$mhc_faiths_binary)
summary(m1)

m1 <- lm(log(aj_df_final_neg_a$Chao1) ~ aj_df_final_neg_a$SequencingDepth + aj_df_final_neg_a$mhc_faiths_binary)
summary(m1)

m1 <- lm(log(aj_df_final_neg_a$Faiths) ~ aj_df_final_neg_a$SequencingDepth + aj_df_final_neg_a$mhc_faiths_binary)
summary(m1)


############# Young AstV- #############

m1 <- lm(log(aj_df_final_neg_y$Observed) ~ aj_df_final_neg_y$SequencingDepth + aj_df_final_neg_y$mhc_faiths_binary)
summary(m1)

m1 <- lm(aj_df_final_neg_y$Shannon ~ aj_df_final_neg_y$SequencingDepth + aj_df_final_neg_y$mhc_faiths_binary)
summary(m1)

m1 <- lm(aj_df_final_neg_y$Chao1 ~ aj_df_final_neg_y$SequencingDepth + aj_df_final_neg_y$mhc_faiths_binary)
summary(m1)

m1 <- lm(aj_df_final_neg_y$Faiths ~ aj_df_final_neg_y$SequencingDepth + aj_df_final_neg_y$mhc_faiths_binary)
summary(m1)




####################################################################################################################################
####################################################### Astro Positives ############################################################
####################################################################################################################################


#########################################
######## Supplementary Figure S7 ########
#########################################

aj_df_long_pos <- subset(aj_df_long, astro == "POS")


FigS7A <- ggplot(aj_df_long_pos, aes(x = alpha_metrics, y = richness, fill= pdist_mean_binary, color= pdist_mean_binary))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(name=c(expression(Mean[pdist]* "     ")), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  scale_fill_manual(name=c(expression(Mean[pdist]* "     ")), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  scale_y_log10()+
  labs(x= " ", y =expression(Log[10] * " alpha diversity")) +
  theme_set(theme_bw()) + 
  my_theme4 + 
  theme(axis.title.y = element_text(size = 14))

FigS7B <- ggplot(aj_df_long_pos, aes(x = alpha_metrics, y = richness, fill= pdistPSS_mean_binary, color= pdistPSS_mean_binary))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(name=c(expression(Mean[pdistPSS])), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  scale_fill_manual(name=c(expression(Mean[pdistPSS])), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  scale_y_log10()+
  labs(x= "Observed     Shannaon     Chao1     Faith's            Observed     Shannaon     Chao1     Faith's
        ASVs           Index            Index        PD                   ASVs           Index            Index          PD       ", y =expression(Log[10] * " alpha diversity")) +
  theme_set(theme_bw()) + 
  my_theme4 + 
  theme(axis.title.y = element_text(size = 14))


ggarrange(FigS7A, FigS7B, labels = c("A", "B"), nrow = 2)





#########################################
### Stats for Supplementary Table S11 ###
#########################################


######## Adult AstV+ #######


m1 <- lm(aj_df_final_pos_a$Observed ~ aj_df_final_pos_a$SequencingDepth + aj_df_final_pos_a$pdist_mean_binary)
summary(m1)

m1 <- lm(aj_df_final_pos_a$Shannon ~ aj_df_final_pos_a$SequencingDepth + aj_df_final_pos_a$pdist_mean_binary)
summary(m1)

m1 <- lm(log(aj_df_final_pos_a$Chao1) ~ aj_df_final_pos_a$SequencingDepth + aj_df_final_pos_a$pdist_mean_binary)
summary(m1)

m1 <- lm(aj_df_final_pos_a$Faiths ~ aj_df_final_pos_a$SequencingDepth + aj_df_final_pos_a$pdist_mean_binary)
summary(m1)


m1 <- lm(aj_df_final_pos_a$Observed ~ aj_df_final_pos_a$SequencingDepth + aj_df_final_pos_a$pdistPSS_mean_binary)
summary(m1)

m1 <- lm(aj_df_final_pos_a$Shannon ~ aj_df_final_pos_a$SequencingDepth + aj_df_final_pos_a$pdistPSS_mean_binary)
summary(m1)

m1 <- lm(log(aj_df_final_pos_a$Chao1) ~ aj_df_final_pos_a$SequencingDepth + aj_df_final_pos_a$pdistPSS_mean_binary)
summary(m1)

m1 <- lm(aj_df_final_pos_a$Faiths ~ aj_df_final_pos_a$SequencingDepth + aj_df_final_pos_a$pdistPSS_mean_binary)
summary(m1)



####### Young AstV + ######

m1 <- lm(log(aj_df_final_pos_y$Observed) ~ aj_df_final_pos_y$SequencingDepth + aj_df_final_pos_y$pdist_mean_binary)
summary(m1)

m1 <- lm(aj_df_final_pos_y$Shannon ~ aj_df_final_pos_y$SequencingDepth + aj_df_final_pos_y$pdist_mean_binary)
summary(m1)

m1 <- lm(aj_df_final_pos_y$Chao1 ~ aj_df_final_pos_y$SequencingDepth + aj_df_final_pos_y$pdist_mean_binary)
summary(m1)

m1 <- lm(log(aj_df_final_pos_y$Faiths) ~ aj_df_final_pos_y$SequencingDepth + aj_df_final_pos_y$pdist_mean_binary)
summary(m1)


m1 <- lm(log(aj_df_final_pos_y$Observed) ~ aj_df_final_pos_y$SequencingDepth + aj_df_final_pos_y$pdistPSS_mean_binary)
summary(m1)

m1 <- lm(aj_df_final_pos_y$Shannon ~ aj_df_final_pos_y$SequencingDepth + aj_df_final_pos_y$pdistPSS_mean_binary)
summary(m1)

m1 <- lm(aj_df_final_pos_y$Chao1 ~ aj_df_final_pos_y$SequencingDepth + aj_df_final_pos_y$pdistPSS_mean_binary)
summary(m1)

m1 <- lm(log(aj_df_final_pos_y$Faiths) ~ aj_df_final_pos_y$SequencingDepth + aj_df_final_pos_y$pdistPSS_mean_binary)
summary(m1)




#########################################
### Stats for Supplementary Table S10 ###
#########################################


######## Adult AstV+ #######

m1 <- lm(aj_df_final_pos_a$Observed ~ aj_df_final_pos_a$SequencingDepth + aj_df_final_pos_a$ST10)
summary(m1)

m1 <- lm(aj_df_final_pos_a$Shannon ~ aj_df_final_pos_a$SequencingDepth + aj_df_final_pos_a$ST10)
summary(m1)

m1 <- lm(log(aj_df_final_pos_a$Chao1) ~ aj_df_final_pos_a$SequencingDepth + aj_df_final_pos_a$ST10)
summary(m1)

m1 <- lm(aj_df_final_pos_a$Faiths ~ aj_df_final_pos_a$SequencingDepth + aj_df_final_pos_a$ST10)
summary(m1)


####### Young AstV + ######

m1 <- lm(log(aj_df_final_pos_y$Observed) ~ aj_df_final_pos_y$SequencingDepth + aj_df_final_pos_y$ST10)
summary(m1)

m1 <- lm(aj_df_final_pos_y$Shannon ~ aj_df_final_pos_y$SequencingDepth + aj_df_final_pos_y$ST10)
summary(m1)

m1 <- lm(aj_df_final_pos_y$Chao1 ~ aj_df_final_pos_y$SequencingDepth + aj_df_final_pos_y$ST10)
summary(m1)

m1 <- lm(log(aj_df_final_pos_y$Faiths) ~ aj_df_final_pos_y$SequencingDepth + aj_df_final_pos_y$ST10)
summary(m1)





####################################################################################################################################
################################################# Continue with RAREFIED analysis ##################################################
####################################################################################################################################


####################################################################################################################################
####################################################### Astro Negatives ############################################################
####################################################################################################################################


#############################################################################################################
######################################### Alpha diversity ###################################################
#############################################################################################################

aj_rare <- readRDS("Aj_class_rare11000.RDS")


# Add alpha-diversities to phyloseq metadata


# take out phyloseq sample data

meta <- data.frame(sample_data(aj_rare))
row.names(meta) # check rownmaes

# make df with all rarefied alpha-diversities

alpha_df <- data.frame(subset(aj_rare_df, select = c("SampleID_total","Observed_rare", "Shannon_rare", "Chao1_rare", "Faiths_rare")))

# join dataframes by SampleID

meta1 <- left_join(meta, alpha_df, by="SampleID_total")
meta1$SampleID_total # check rownames

row.names(meta1) <- row.names(meta) # we must add the rownames again, otherwise phyloseq can't recognize it as sample data

# map back

sample_data(aj_rare) <- meta1
sample_data(aj_rare) # check it worked



# data frame of metadata with all necessary information

aj_df <- data.frame(sample_data(aj_rare))

# correct format
aj_df$ST1 <- as.factor(aj_df$ST1)
aj_df$ST2 <- as.factor(aj_df$ST2)
aj_df$ST3 <- as.factor(aj_df$ST3)
aj_df$ST4 <- as.factor(aj_df$ST4)
aj_df$ST5 <- as.factor(aj_df$ST5)
aj_df$ST6 <- as.factor(aj_df$ST6)
aj_df$ST7 <- as.factor(aj_df$ST7)
aj_df$ST8 <- as.factor(aj_df$ST8)
aj_df$ST9 <- as.factor(aj_df$ST9)
aj_df$ST10 <- as.factor(aj_df$ST10)


# make into long format for plotting

aj_df_long<-gather(aj_df, alpha_metrics, richness, Observed_rare,Shannon_rare,Chao1_rare,Faiths_rare, factor_key=TRUE)
aj_df_long$alpha_metrics # check it worked




#########################################
############ Final Figure 3 #############
#########################################

aj_df_long_neg <- subset(aj_df_long, astro == "NEG")

age_label <- c("A" = "Adult", "Y" = "Young")


Fig3A <- ggplot(aj_df_long_neg, aes(x = alpha_metrics, y = richness, fill=ST5, color= ST5))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75), size=1)+
  theme_bw(base_size = 12)+
  scale_color_manual(labels = c("absent","present"),values=c("#c29844","#6ab166"))+
  scale_fill_manual(labels = c("absent","present"),values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= " ", y =expression(Log[10] * " alpha diversity")) +
  theme_set(theme_bw()) + 
  my_theme4 + 
  theme(axis.title.y = element_text(size = 14))

Fig3B <- ggplot(aj_df_long_neg, aes(x = alpha_metrics, y = richness, fill=ST6, color= ST6))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75), size=1)+
  theme_bw(base_size = 12)+
  scale_color_manual(labels = c("absent","present"),values=c("#c29844","#6ab166"))+
  scale_fill_manual(labels = c("absent","present"),values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+  
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= "Observed     Shannaon     Chao1     Faith's            Observed     Shannaon     Chao1     Faith's
        ASVs           Index            Index        PD                   ASVs           Index            Index          PD       ", y =expression(Log[10] * " alpha diversity")) +
  theme_set(theme_bw()) + 
  my_theme4 + 
  theme(axis.title.y = element_text(size = 14))


Figure3 <- ggarrange(Fig3A,Fig3B, labels = c("A", "B"), nrow = 2, heights = c(1,1.05))




#########################################
######## Supplementary Figure S5 ########
#########################################


FigS5A <- ggplot(aj_df_long_neg, aes(x = alpha_metrics, y = richness, fill= pdist_mean_binary, color= pdist_mean_binary))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(name=c(expression(Mean[pdist]* "     ")), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  scale_fill_manual(name=c(expression(Mean[pdist]* "     ")), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= " ", y = " ") +
  theme_set(theme_bw()) + 
  my_theme4

FigS5B <- ggplot(aj_df_long_neg, aes(x = alpha_metrics, y = richness, fill=pdist_sum_binary, color= pdist_sum_binary))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(name=c(expression(Sum[pdist]* "      ")), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  scale_fill_manual(name=c(expression(Sum[pdist]* "      ")), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= " ", y =" ") +
  theme_set(theme_bw()) + 
  my_theme4

FigS5C <- ggplot(aj_df_long_neg, aes(x = alpha_metrics, y = richness, fill= pdistPSS_mean_binary, color= pdistPSS_mean_binary))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(name=c(expression(Mean[pdistPSS])), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  scale_fill_manual(name=c(expression(Mean[pdistPSS])), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= " ", y =expression(Log[10] * " alpha diversity")) +
  theme_set(theme_bw()) + 
  my_theme4 + 
  theme(axis.title.y = element_text(size = 14))

FigS5D <- ggplot(aj_df_long_neg, aes(x = alpha_metrics, y = richness, fill= pdistPSS_sum_binary, color= pdistPSS_sum_binary))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(name=c(expression(Sum[pdistPSS]* " ")), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  scale_fill_manual(name=c(expression(Sum[pdistPSS]* " ")), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= " ", y =" ") +
  theme_set(theme_bw()) + 
  my_theme4

FigS5E <- ggplot(aj_df_long_neg, aes(x = alpha_metrics, y = richness, fill= mhc_faiths_binary, color= mhc_faiths_binary))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(name=c(expression(MHC[FaithsPD])), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  scale_fill_manual(name=c(expression(MHC[FaithsPD])), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= "Observed     Shannaon     Chao1     Faith's            Observed     Shannaon     Chao1     Faith's
        ASVs           Index            Index        PD                   ASVs           Index            Index          PD       ", y = " ") +
  theme_set(theme_bw()) + 
  my_theme4



ggarrange(FigS5A,FigS5B,FigS5C,FigS5D,FigS5E,nrow=5)




#############################################################################################################
########################################## Beta diversity ###################################################
#############################################################################################################

# load functions 

#https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

## load R function
source("phyloseq_to_df.R")
library(BiodiversityR)



###############################################
################## Adults #####################
###############################################

# subset only astro-negative adults

aj_rare_neg <- subset_samples(aj_rare, astro == "NEG")

aj_rare_neg_a <- subset_samples(aj_rare_neg, age == "A")


# agglomorate taxa at genus level

aj_rare_neg_a_compositional <- microbiome::transform(aj_rare_neg_a, "compositional")
aj_neg_merged <- tax_glom(aj_rare_neg_a_compositional, taxrank="Genus") # agglomerate taxa to Genus level
# aj_neg_merged 
# otu_table()   OTU Table:         [ 370 taxa and 91 samples ]


df_aj_neg_a <- pssd2veg(aj_neg_merged)

pred_var_astro<- as.data.frame(df_aj_neg_a %>% dplyr::select(Nr_ST, AA_pdist_mean, ST1, ST2, ST3, ST4, ST5, ST6, ST7, ST9, ST10))
ASV_matrix_astro<- psotu2veg(aj_neg_merged)

######### access taxonomy ####
sample_data(aj_neg_merged)$SampleID_total<-make.names(sample_data(aj_neg_merged)$SampleID_total)
sample_names(aj_neg_merged) <- sample_data(aj_neg_merged)$SampleID_total

# phyloseq to data frame
Aj_phylo_na_df_astro<- phyloseq_to_df(aj_neg_merged)
Aj_phylo_na_table_astro<-phyloseq_to_df(aj_neg_merged)

#apply transformation
ASV_matrix_astro_trans <- disttransform(ASV_matrix_astro, method="chord")

###CCAs
cca_aj_astro <- cca(ASV_matrix_astro_trans ~ Nr_ST+AA_pdist_mean+ST1+ST2+ST3+ST4+ST5+ST6+ST7+ST9+ST10, scale=FALSE, data=df_aj_neg_a)
cca_aj_astro


# Testing the significance of the CCA model:
anova.cca(cca_aj_astro, permutations = how(nperm = 999)) # by anova
# p=0.098 .

# Importance of components (axis)
summary(cca_aj_astro)$concont

# Importance of estimates
envfit_m <-envfit(cca_aj_astro,pred_var_astro, na.rm=TRUE, permutations=10000)
envfit_m$vectors



######################################################################################################
############################# Repeat with subset of variables ########################################
######################################################################################################

df_aj_neg_a<- pssd2veg(aj_neg_merged)

pred_var_astro<- as.data.frame(df_aj_neg_a %>% dplyr::select(ST2, ST4, ST5, ST6, ST7, ST10))
ASV_matrix_astro<- psotu2veg(aj_neg_merged)

######### access taxonomy ####
sample_data(aj_neg_merged)$SampleID_total<-make.names(sample_data(aj_neg_merged)$SampleID_total)
sample_names(aj_neg_merged) <- sample_data(aj_neg_merged)$SampleID_total

# phyloseq to data frame
Aj_phylo_na_df_astro<- phyloseq_to_df(aj_neg_merged)
Aj_phylo_na_table_astro<-phyloseq_to_df(aj_neg_merged)

# apply transformation
ASV_matrix_astro_trans <- disttransform(ASV_matrix_astro, method="chord")

# CCA
set.seed(1) # set seed for reproducibility
cca_aj_astro <- cca(ASV_matrix_astro_trans ~ ST2 + ST4 + ST5 + ST6 + ST7 + ST10, scale=FALSE, data=df_aj_neg_a)
cca_aj_astro 
# proportion explained by constrained variables: 0.07368

# Testing the significance of the CCA model:
set.seed(1) # set seed for reproducibility
anova.cca(cca_aj_astro, permutations = how(nperm = 999)) 
# p=0.156

# Importance of components (axis)
summary(cca_aj_astro)$concont


# Importance of estimates

########################################
### Stats for Supplementary Table S8 ###
########################################

set.seed(1) # set seed for reproducibility
envfit_m <-envfit(cca_aj_astro,pred_var_astro, na.rm=TRUE, permutations=10000)
envfit_m$vectors

# Plot basic CCA
# Scaling 1: species (ASVs) scores scaled to relative eigenvalues, sites (bat individuals) are weighted averages of the species (bacterial genera)  
plot(cca_aj_astro, scaling=1, display=c('sp', 'lc', 'cn'))


# Get data in the correct shape for final ggplot

cca_aj_astro # all info is stored in this object as separate lists

# select all dimensions we want to plot
# we must also specify the correct scaling

vare_spp_sco <- scores(cca_aj_astro, scaling=1, display = "species")
vare_sam_sco <- scores(cca_aj_astro, scaling=1, display = "sites")
vare_env_sco <- scores(cca_aj_astro, scaling=1, display = "bp")[c(1,3,4,5,6),]    # MHC div estimates

vare_spp_tbl <- as_tibble(vare_spp_sco)
vare_sam_tbl <- as_tibble(vare_sam_sco)
vare_env_tbl <- as_tibble(vare_env_sco)


# the plot.cca function transforms arrow length depending on the scaling method
# thus scaling of arrows won't be correct for the environmental variables
# thus we must use a multiplier on the extracted axis scores to have the correct arrow length

rescaled_astro <- vare_env_tbl %>% 
  select(CCA1, CCA2) %>%
  as.matrix() *  20.77657

vare_spp_tbl <- mutate(vare_spp_tbl, vgntxt=rownames(vare_spp_sco), ccatype = "species")
vare_sam_tbl <- mutate(vare_sam_tbl, vgntxt=rownames(vare_sam_sco), ccatype = "sites")
vare_env_tbl <- mutate(vare_env_tbl, vgntxt=rownames(vare_env_sco), ccatype = "bp")

# save files and combine
#write.csv(vare_spp_tbl,"vare_spp_tbl_neg_a.csv")
#write.csv(vare_sam_tbl,"vare_sam_tbl_neg_a.csv")
#write.csv(vare_env_tbl, "vare_env_tbl_neg_a.csv")
#write.csv(rescaled_astro, "rescaled_astro_neg_a.csv")

# Remove X before ASVs that start with a number and import again as "vare_spp_tbl", so we can join with taxonomy info

vare_spp_tbl # taxa names that were found in the CCA

taxonomy <- tax_table(aj_rare) [, 1:6]  # we don't include bacterial species taxonomy here since we agglomerated taxa at genus level
rownames(taxonomy)

taxonomy1 <- data.frame(vgntxt = row.names(taxonomy), taxonomy)
names(taxonomy1)
head(taxonomy1$vgntxt)

# make names a column
taxonomy1$names <- rownames(taxonomy)


# we join everything, join will introduce NAs for ASVs found in the first table but not in the second, we can remove these later 
cca_taxonomy <- left_join(taxonomy1, vare_spp_tbl, by="vgntxt")

# remove NAs (these are ASVs that are not present in the merged object we creasted for the CCA)
cca_taxonomy1 <- subset(cca_taxonomy, CCA1 != "NA")
names(cca_taxonomy1)

# check everything was joined correctly
head(cca_taxonomy1)
tail(cca_taxonomy1)

vare_spp_tbl <- cca_taxonomy1
names(vare_spp_tbl)

# change order of columns 
vare_spp_tbl <- vare_spp_tbl[,c(10,11,1,12,3,4,5,6,7,8,9)]

# save file


#########################################
########### Final Figure 4A #############
#########################################

library(ggrepel)

######## Plot with manually organized df

cca_data <- read.csv("cca_data_neg_a.csv")
head(cca_data)

# for plotting microbiota as dots
data = filter(cca_data, ccatype=="species")
str(data)
data$CCA1 <- as.numeric(as.character(data$CCA1))
data$CCA2 <- as.numeric(as.character(data$CCA2))

# for plotting Diversity estimates
data1=filter(cca_data, ccatype=="bp")
str(data1)
data1$CCA1 <- as.numeric(as.character(data1$CCA1))
data1$CCA2 <- as.numeric(as.character(data1$CCA2))


# final plot
# points are taxa
# the MHC diversity estimates are shown as arrows

cca_plot <- ggplot() + 
  theme_bw(base_size = 14) +
  geom_point(data=data, aes(x = CCA1, y = CCA2, fill = vgntxt)) + 
  geom_vline(aes(xintercept = 0), linetype = "dashed", alpha = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.5) + # this will create the dashed coordinates, makes the plot more cca-like
  geom_label_repel(data=data[which(data$plot==1),], aes(x = CCA1, y = CCA2,label=Genus),
                   segment.color = "transparent",force=5, size=rel(4.5), segment.size =0, parse=TRUE, seed = 1)+
  geom_segment(data=data1, aes(x=0, y=0, xend=CCA1 , yend=CCA2), size=1.1, arrow=arrow(length = unit(0.3,"cm")), color=c("red","red","red","red","red")) +
  ylab("CCA2 (22.14 %)")+ 
  xlab("CCA1 (29.33 %)")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits=c(-10,10), breaks=seq(-10,10, by=5)) +
  scale_y_continuous(limits=c(-10,5), breaks=seq(-10,5, by=5)) + theme (legend.position = "none") 

cca_plot_1  <- cca_plot + guides(fill=FALSE)+ guides(size=FALSE)

# add labels for MHC diversity estimates
cca_plot_final <- cca_plot_1 +
  geom_label_repel(data=data1, aes(x=CCA1, y=CCA2), label=c("ST2", "ST5", "ST6", "ST7", "ST10") ,nudge_x=-.1, nudge_y=-.1, size=4, color= "red")

cca_plot_adult_neg <- cca_plot_final 





###############################################
################### Young #####################
###############################################

# subset only astro-negative youngs

aj_rare_neg <- subset_samples(aj_rare, astro == "NEG")

aj_rare_neg_y <- subset_samples(aj_rare_neg, age == "Y")


# agglomorate taxa at genus level

aj_rare_neg_y_compositional <- microbiome::transform(aj_rare_neg_y, "compositional")
aj_neg_merged_y <- tax_glom(aj_rare_neg_y_compositional, taxrank="Genus") # agglomerate taxa to Genus level
#otu_table()   OTU Table:         [ 370 taxa and 61 samples ]


df_aj_neg_y <- pssd2veg(aj_neg_merged_y)

pred_var_astro <- as.data.frame(df_aj_neg_y %>% dplyr::select(Nr_ST, AA_pdist_mean, ST1, ST2, ST3, ST4, ST5, ST6, ST7, ST9, ST10))
ASV_matrix_astro <- psotu2veg(aj_neg_merged_y)

######### access taxonomy ####
sample_data(aj_neg_merged_y)$SampleID_total<-make.names(sample_data(aj_neg_merged_y)$SampleID_total)
sample_names(aj_neg_merged_y) <- sample_data(aj_neg_merged_y)$SampleID_total

# phyloseq to data frame
Aj_phylo_na_df_astro<- phyloseq_to_df(aj_neg_merged_y)
Aj_phylo_na_table_astro<-phyloseq_to_df(aj_neg_merged_y)

#apply transformation
ASV_matrix_astro_trans <- disttransform(ASV_matrix_astro, method="chord")

###CCA
cca_aj_astro <- cca(ASV_matrix_astro_trans ~ Nr_ST+AA_pdist_mean+ST1+ST2+ST3+ST4+ST5+ST6+ST7+ST9+ST10, scale=FALSE, data= df_aj_neg_y)
cca_aj_astro


# Testing the significance of the CCA model:
anova.cca(cca_aj_astro, permutations = how(nperm = 999)) # by anova
#  0.067 .

# Importance of components (axis)
summary(cca_aj_astro)$concont

# Importance of estimates
envfit_m <-envfit(cca_aj_astro,pred_var_astro, na.rm=TRUE, permutations=10000)
envfit_m$vectors



######################################################################################################
############################# Repeat with subset of variables ########################################
######################################################################################################

df_aj<- pssd2veg(aj_neg_merged_y)

pred_var_astro<- as.data.frame(df_aj %>% dplyr::select(Nr_ST, AA_pdist_mean, ST2, ST4, ST9, ST10))
ASV_matrix_astro<- psotu2veg(aj_neg_merged_y)

######### access taxonomy ####
sample_data(aj_neg_merged_y)$SampleID_total<-make.names(sample_data(aj_neg_merged_y)$SampleID_total)
sample_names(aj_neg_merged_y) <- sample_data(aj_neg_merged_y)$SampleID_total

# phyloseq to data frame
Aj_phylo_na_df_astro<- phyloseq_to_df(aj_neg_merged_y)
Aj_phylo_na_table_astro<-phyloseq_to_df(aj_neg_merged_y)

#apply transformation
ASV_matrix_astro_trans <- disttransform(ASV_matrix_astro, method="chord")

###CCAs
set.seed(1)
cca_aj_astro <- cca(ASV_matrix_astro_trans ~ Nr_ST+AA_pdist_mean+ST2+ST4+ST9+ST10, scale=FALSE, data=df_aj)
cca_aj_astro

# Testing the significance of the CCA model:
set.seed(1)
anova.cca(cca_aj_astro, permutations = how(nperm = 999)) # by anova             
#   0.027 *


vare_spp_sco<-scores(cca_aj_astro, display="species")
vare_sam_sco<-scores(cca_aj_astro, display="sites")
vare_all_sco<- scores(cca_aj_astro)


# Importance of components (axis)
summary(cca_aj_astro)$concont


# Importance of estimates

########################################
### Stats for Supplementary Table S9 ###
########################################

set.seed(1)
envfit_m <-envfit(cca_aj_astro,pred_var_astro, na.rm=TRUE, permutations=10000)
envfit_m$vectors

# Scaling 1: species (ASVs) scores scaled to relative eigenvalues, sites (individuals) are weighted averages of the species (ASVs)  
plot(cca_aj_astro, scaling=1, display=c('sp', 'lc', 'cn'))



#### Get data for plotting in the correct shape

cca_aj_astro # all info is stored in this object as separate list

vare_spp_sco <- scores(cca_aj_astro, scaling=1, display = "species")
vare_sam_sco <- scores(cca_aj_astro, scaling=1, display = "sites")
vare_env_sco <- scores(cca_aj_astro, scaling=1, display = "bp")[1:6,]    # MHC div estimates

vare_spp_tbl <- as_tibble(vare_spp_sco)
vare_sam_tbl <- as_tibble(vare_sam_sco)
vare_env_tbl <- as_tibble(vare_env_sco)


# the plot.cca function transforms arrow length depending on the scaling method
# thus scaling of arrows won't be correct for the environmental variables
# thus we must use a multiplier on the extracted axis scores to have the correct arrow length

rescaled_astro <- vare_env_tbl %>% 
  select(CCA1, CCA2) %>%
  as.matrix() * 16.37419

vare_spp_tbl <- mutate(vare_spp_tbl, vgntxt=rownames(vare_spp_sco), ccatype = "species")
vare_sam_tbl <- mutate(vare_sam_tbl, vgntxt=rownames(vare_sam_sco), ccatype = "sites")
vare_env_tbl <- mutate(vare_env_tbl, vgntxt=rownames(vare_env_sco), ccatype = "bp")

# save files and combine
#write.csv(vare_spp_tbl,"vare_spp_tbl_neg_y.csv")
#write.csv(vare_sam_tbl,"vare_sam_tbl_neg_y.csv")
#write.csv(vare_env_tbl, "vare_env_tbl_neg_y.csv")
#write.csv(rescaled_astro, "rescaled_astro_neg_y.csv")

# Remove X before ASVs that start with a number and import again as "vare_spp_tbl", so we can join with taxonomy info

vare_spp_tbl # taxa names that were found in the CCA

taxonomy <- tax_table(aj_rare) [, 1:6]
rownames(taxonomy)

taxonomy1 <- data.frame(vgntxt = row.names(taxonomy), taxonomy)
names(taxonomy1)
head(taxonomy1$vgntxt)

# make names a column
taxonomy1$names <- rownames(taxonomy)


# we join everything, join will introduce NAs for ASVs found in the first table but not in the second, we can remove these later 
cca_taxonomy <- left_join(taxonomy1, vare_spp_tbl, by="vgntxt")

# remove NAs (these are ASVs that are not present in the merged object we creasted for the CCA)
cca_taxonomy1 <- subset(cca_taxonomy, CCA1 != "NA")

# check everything was joined correctly
head(cca_taxonomy1)
tail(cca_taxonomy1)

vare_spp_tbl <- cca_taxonomy1

# change order of columns 
vare_spp_tbl <- vare_spp_tbl[,c(10,11,1,12,3,4,5,6,7,8,9)]




#########################################
########### Final Figure 4B #############
#########################################

######## Plot with manually organized df

# for plotting we can set a threshold, so only taxa names at the outer limits of the plot are included to reduce clutter
# therefore include a labelling column in the df

cca_data <- read.csv("cca_data_neg_y.csv")
head(cca_data)

# for plotting microbiota as dots
data = filter(cca_data, ccatype=="species")
str(data)
data$CCA1 <- as.numeric(as.character(data$CCA1))
data$CCA2 <- as.numeric(as.character(data$CCA2))

# for plotting Diversity estimates
data1=filter(cca_data, ccatype=="bp")
data1$CCA1 <- as.numeric(as.character(data1$CCA1))
data1$CCA2 <- as.numeric(as.character(data1$CCA2))


# final plot
# points are taxa
# names of taxa on genus level only included if they are > 1 or < -1 on CCA1 OR CCA2 to reduce clutter

cca_plot <- ggplot() + 
  theme_bw(base_size = 14) +
  geom_point(data=data, aes(x = CCA1, y = CCA2, fill = vgntxt)) + 
  geom_vline(aes(xintercept = 0), linetype = "dashed", alpha = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.5) + # this will create the dashed coordinates, makes the plot more cca-like
  geom_label_repel(data=data[which(data$plot==1),], aes(x = CCA1, y = CCA2,label=Genus),
                   segment.color = "transparent", nudge_y = .01,force=5, size=rel(4.5), segment.size =0, parse=TRUE)+
  geom_segment(data=data1, aes(x=0, y=0, xend=CCA1 , yend=CCA2), size=1.1, arrow=arrow(length = unit(0.3,"cm")), color=c("red","red","red","red","red","red")) +
  ylab("CCA2 (27.09 %)")+ 
  xlab("CCA1 (31.95 %)")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits=c(-10,10), breaks=seq(-10,10, by=5)) +
  scale_y_continuous(limits=c(-4,8), breaks=seq(-4,8, by=2)) + theme (legend.position = "none") 

cca_plot_1  <- cca_plot + guides(fill=FALSE)+ guides(size=FALSE)

cca_plot_final <- cca_plot_1 +
  geom_label_repel(data=data1, aes(x=CCA1, y=CCA2), label=c(expression(Nr[ST]), expression(Mean[pdist]), "ST2", "ST4", "ST9", "ST10") ,nudge_x=-.1, nudge_y=-.1, size=4, color= "red")

cca_plot_young_neg <- cca_plot_final




##############################
####### Final Figure 4 #######
##############################

cca_negative <- ggarrange(cca_plot_adult_neg,cca_plot_young_neg, ncol=1, nrow=2, labels=c("A","B"))


# note that changing the plot window size changes label position slightl
# for label position as in publication
# export as .jpeg with width 880 and height 1110







####################################################################################################################################
####################################################### Astro Positives ############################################################
####################################################################################################################################


#############################################################################################################
######################################### Alpha diversity ###################################################
#############################################################################################################


#########################################
############ Final Figure 5 #############
#########################################

aj_df_long_pos <- subset(aj_df_long, astro == "POS")

age_label <- c("A" = "Adult", "Y" = "Young")


Fig5A <- ggplot(aj_df_long_pos, aes(x = alpha_metrics, y = richness, fill= pdist_mean_binary, color= pdist_mean_binary))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75), size = 1)+
  theme_bw(base_size = 12)+
  scale_color_manual(name=c(expression(Mean[pdist]* "     ")), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  scale_fill_manual(name=c(expression(Mean[pdist]* "     ")), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= " ", y =expression(Log[10] * " alpha diversity")) +
  theme_set(theme_bw()) + 
  my_theme4 + 
  theme(axis.title.y = element_text(size = 14))

Fig5B <- ggplot(aj_df_long_pos, aes(x = alpha_metrics, y = richness, fill= pdistPSS_mean_binary, color= pdistPSS_mean_binary))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75), size = 1)+
  theme_bw(base_size = 12)+
  scale_color_manual(name=c(expression(Mean[pdistPSS])), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  scale_fill_manual(name=c(expression(Mean[pdistPSS])), labels = c("Above","Below"), values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+  
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= "Observed     Shannon     Chao1     Faith's           Observed     Shannon       Chao1      Faith's
        ASVs           Index            Index        PD                 ASVs           Index            Index          PD       ", y =expression(Log[10] * " alpha diversity")) +
  theme_set(theme_bw()) + 
  my_theme4 + 
  theme(axis.title.y = element_text(size = 14))


Figure5 <- ggarrange(Fig5A,Fig5B, labels = c("A", "B"), nrow = 2, heights = c(1,1.05))
Figure5




#########################################
######## Supplementary Figure S6 ########
#########################################


FigS6A <- ggplot(aj_df_long_pos, aes(x = alpha_metrics, y = richness, fill=ST1, color= ST1))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(labels = c("absent","present"),values=c("#c29844","#6ab166"))+
  scale_fill_manual(labels = c("absent","present"),values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= " ", y = " ") +
  theme_set(theme_bw()) + 
  my_theme4+ 
  theme(axis.title.y = element_text(size = 14))

FigS6B <- ggplot(aj_df_long_pos, aes(x = alpha_metrics, y = richness, fill=ST2, color= ST2))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(labels = c("absent","present"),values=c("#c29844","#6ab166"))+
  scale_fill_manual(labels = c("absent","present"),values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= " ", y = " ") +
  theme_set(theme_bw()) + 
  my_theme4+ 
  theme(axis.title.y = element_text(size = 14))

FigS6C <- ggplot(aj_df_long_pos, aes(x = alpha_metrics, y = richness, fill=ST3, color= ST3))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+
  scale_fill_manual(labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= " ", y = " ") +
  theme_set(theme_bw()) + 
  my_theme4

FigS6D <- ggplot(aj_df_long_pos, aes(x = alpha_metrics, y = richness, fill=ST4, color= ST4))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+
  scale_fill_manual(labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= " ", y = " ") +
  theme_set(theme_bw()) + 
  my_theme4

FigS6E <- ggplot(aj_df_long_pos, aes(x = alpha_metrics, y = richness, fill=ST5, color= ST5))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+
  scale_fill_manual(labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= " ", y =expression(Log[10] * " alpha diversity")) +
  theme_set(theme_bw()) + 
  my_theme4 + 
  theme(axis.title.y = element_text(size = 14))

FigS6F <- ggplot(aj_df_long_pos, aes(x = alpha_metrics, y = richness, fill=ST6, color= ST6))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+
  scale_fill_manual(labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= " ", y = " ") +
  theme_set(theme_bw()) + 
  my_theme4

FigS6G <- ggplot(aj_df_long_pos, aes(x = alpha_metrics, y = richness, fill=ST7, color= ST7))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+
  scale_fill_manual(labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= " ", y = " ") +
  theme_set(theme_bw()) + 
  my_theme4

FigS6H <- ggplot(aj_df_long_pos, aes(x = alpha_metrics, y = richness, fill=ST9, color= ST9))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+
  scale_fill_manual(labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= " ", y = " ") +
  theme_set(theme_bw()) + 
  my_theme4

FigS6I <- ggplot(aj_df_long_pos, aes(x = alpha_metrics, y = richness, fill=ST10, color= ST10))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  theme_bw(base_size = 12)+
  scale_color_manual(labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+
  scale_fill_manual(labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+
  facet_wrap(~age, labeller = as_labeller(age_label))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.3f", as.numeric(..p.format..)))) +
  scale_y_log10()+
  labs(x= "Observed      Shannaon     Chao1     Faith's                Observed      Shannaon     Chao1     Faith's
        ASVs            Index            Index        PD                       ASVs            Index            Index          PD       ", y =" ") +
  theme_set(theme_bw()) + 
  my_theme4



ggarrange(FigS6A,FigS6B,FigS6C,FigS6D,FigS6E,FigS6F,FigS6G,FigS6H,FigS6I, nrow=9, heights=c(1,1,1,1,1,1,1,1,1.05))







#############################################################################################################
########################################## Beta diversity ###################################################
#############################################################################################################


###############################################
################## Adults #####################
###############################################

# subset only astro-positive adults

aj_rare_pos <- subset_samples(aj_rare, astro == "POS")

aj_rare_pos_a <- subset_samples(aj_rare_pos, age == "A")


# agglomerate taxa at genus level

aj_rare_pos_a_compositional <- microbiome::transform(aj_rare_pos_a, "compositional")
aj_pos_merged_a <- tax_glom(aj_rare_pos_a_compositional, taxrank ="Genus")
#otu_table()   OTU Table:         [ 370 taxa and 40 samples ]


df_aj<- pssd2veg(aj_pos_merged_a)

pred_var_astro<- as.data.frame(df_aj %>% dplyr::select(Nr_ST, AA_pdist_mean,ST1,ST2,ST3,ST4,ST5,ST6,ST7,ST9,ST10))
ASV_matrix_astro<- psotu2veg(aj_pos_merged_a)

# access taxonomy
sample_data(aj_pos_merged_a)$SampleID_total<-make.names(sample_data(aj_pos_merged_a)$SampleID_total)
sample_names(aj_pos_merged_a) <- sample_data(aj_pos_merged_a)$SampleID_total

# phyloseq to data frame
Aj_phylo_na_df_astro<- phyloseq_to_df(aj_pos_merged_a)
Aj_phylo_na_table_astro<-phyloseq_to_df(aj_pos_merged_a)

# apply transformation
ASV_matrix_astro_trans <- disttransform(ASV_matrix_astro, method="chord")

# CCA
cca_aj_astro <- cca(ASV_matrix_astro_trans ~ Nr_ST+AA_pdist_mean+ST1+ST2+ST3+ST4+ST5+ST6+ST7+ST9+ST10, scale=FALSE, data=df_aj)
cca_aj_astro


# Testing the significance of the CCA model:
anova.cca(cca_aj_astro, permutations = how(nperm = 999)) # by anova

# Importance of components (axis)
summary(cca_aj_astro)$concont

# Importance of estimates
envfit_m <-envfit(cca_aj_astro,pred_var_astro, na.rm=TRUE, permutations=10000)
envfit_m$vectors




######################################################################################################
############################# Repeat with subset of variables ########################################
######################################################################################################


df_aj<- pssd2veg(aj_pos_merged_a)

pred_var_astro<- as.data.frame(df_aj %>% dplyr::select(ST1,ST7))
ASV_matrix_astro<- psotu2veg(aj_pos_merged_a)

######### access taxonomy ####
sample_data(aj_pos_merged_a)$SampleID_total<-make.names(sample_data(aj_pos_merged_a)$SampleID_total)
sample_names(aj_pos_merged_a) <- sample_data(aj_pos_merged_a)$SampleID_total

# phyloseq to data frame
Aj_phylo_na_df_astro<- phyloseq_to_df(aj_pos_merged_a)
Aj_phylo_na_table_astro<-phyloseq_to_df(aj_pos_merged_a)

#apply transformation
ASV_matrix_astro_trans <- disttransform(ASV_matrix_astro, method="chord")

###CCAs
set.seed(1)
cca_aj_astro <- cca(ASV_matrix_astro_trans ~ ST1+ST7, scale=FALSE, data=df_aj)
cca_aj_astro


# Testing the significance of the CCA model:
set.seed(1)
anova.cca(cca_aj_astro, permutations = how(nperm = 999)) 
# 0.064

# Importance of components (axis)
summary(cca_aj_astro)$concont



# Importance of estimates

########################################
### Stats for Supplementary Table S12 ##
########################################

set.seed(1)
envfit_m <-envfit(cca_aj_astro,pred_var_astro, na.rm=TRUE, permutations=10000)
envfit_m$vectors

# Scaling 1: species (ASVs) scores scaled to relative eigenvalues, sites (individuals) are weighted averages of the species (ASVs)  
plot(cca_aj_astro, scaling=1, display=c('sp', 'lc', 'cn'))



#### Get data for plotting in the correct shape

cca_aj_astro # all info is stored in this object as separate list

vare_spp_sco <- scores(cca_aj_astro, scaling=1, display = "species")
vare_sam_sco <- scores(cca_aj_astro, scaling=1, display = "sites")
vare_env_sco <- scores(cca_aj_astro, scaling=1, display = "bp")[1:2,]    # MHC div estimates

vare_spp_tbl <- as_tibble(vare_spp_sco)
vare_sam_tbl <- as_tibble(vare_sam_sco)
vare_env_tbl <- as_tibble(vare_env_sco)


# the plot.cca function transforms arrow length depending on the scaling method
# thus scaling of arrows won't be correct for the environmental variables
# thus we must use a multiplier on the extracted axis scores to have the correct arrow length

rescaled_astro <- vare_env_tbl %>% 
  select(CCA1, CCA2) %>%
  as.matrix() * 6.748121

vare_spp_tbl <- mutate(vare_spp_tbl, vgntxt=rownames(vare_spp_sco), ccatype = "species")
vare_sam_tbl <- mutate(vare_sam_tbl, vgntxt=rownames(vare_sam_sco), ccatype = "sites")
vare_env_tbl <- mutate(vare_env_tbl, vgntxt=rownames(vare_env_sco), ccatype = "bp")

write.csv(vare_spp_tbl, "vare_spp_tbl_pos_a.csv")
write.csv(vare_sam_tbl,"vare_sam_tbl_pos_a.csv")
write.csv(vare_env_tbl, "vare_env_tbl_pos_a.csv")
write.csv(rescaled_astro, "rescaled_astro_pos_a.csv")


# Remove X before ASVs that start with a number and import again, so we can join with taxonomy info from biom file
vare_spp_tbl <- read.csv("D:\\Ramona_F\\1PhD Projects\\Bats\\R Analysis\\1New_Analysis_November2020\\beta_div_analysis\\vare_spp_tbl_pos_a.csv")
vare_spp_tbl # ASV names that were found in the CCA


taxonomy <- tax_table(aj_rare) [, 1:6]  # we don't include bacterial species taxonomy here since we agglomerated taxa at genus level
rownames(taxonomy)

taxonomy1 <- data.frame(vgntxt = row.names(taxonomy), taxonomy)
names(taxonomy1)
head(taxonomy1$vgntxt)

# make names a column
taxonomy1$names <- rownames(taxonomy)


# we join everything, join will introduce NAs for ASVs found in the first table but not in the second, we can remove these later 
cca_taxonomy <- left_join(taxonomy1, vare_spp_tbl, by="vgntxt")

# remove NAs (these are ASVs that are not present in the merged object we creasted for the CCA)
cca_taxonomy1 <- subset(cca_taxonomy, CCA1 != "NA")
names(cca_taxonomy1)

# check everything was joined correctly
head(cca_taxonomy1)
tail(cca_taxonomy1)

vare_spp_tbl <- cca_taxonomy1
names(vare_spp_tbl)

# change order of columns 
vare_spp_tbl <- vare_spp_tbl[,c(10,11,1,12,3,4,5,6,7,8,9)]

# save file


#################################################
############### Final Figure 6A #################
#################################################

#################################################
########### CCA AstV-Positive Adult #############
#################################################

library(ggrepel)

######## Plot with manually organized df

# final plot

cca_data <- read.csv("cca_data_pos_a.csv")
head(cca_data)

# for plotting bacterial genera as dots
data = filter(cca_data, ccatype=="species")
str(data)
data$CCA1 <- as.numeric(as.character(data$CCA1))
data$CCA2 <- as.numeric(as.character(data$CCA2))

# for plotting diversity estimates
data1=filter(cca_data, ccatype=="bp")
data1$CCA1 <- as.numeric(as.character(data1$CCA1))
data1$CCA2 <- as.numeric(as.character(data1$CCA2))


cca_plot <- ggplot() + 
  theme_bw(base_size = 14) +
  geom_point(data=data, aes(x = CCA1, y = CCA2, fill = vgntxt)) + 
  geom_vline(aes(xintercept = 0), linetype = "dashed", alpha = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.5) + # this will create the dashed coordinates
  geom_label_repel(data=data[which(data$plot==1),], aes(x = CCA1, y = CCA2,label=Genus),
                   segment.color = "transparent", size=rel(4.5), segment.size =0, parse=TRUE)+
  geom_segment(data=data1, aes(x=0, y=0, xend=CCA1 , yend=CCA2), size=1.1, arrow=arrow(length = unit(0.3,"cm")), color=c("red","red")) +
  ylab("CCA2 (26.18 %)")+ 
  xlab("CCA1 (73.82 %)")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits=c(-5,5), breaks=seq(-5,5, by=5)) +
  scale_y_continuous(limits=c(-5,5), breaks=seq(-5,5, by=5))+ theme (legend.position = "none") 

cca_plot_1  <- cca_plot + guides(fill=FALSE)+ guides(size=FALSE)
cca_plot_1 

cca_plot_final_adult_pos <- cca_plot_1 +
  geom_label_repel(data=data1, aes(x=CCA1, y=CCA2), label=c("ST1", "ST7") ,nudge_x=-.1, nudge_y=-.1, size=4, color = "red")





###############################################
################### Young #####################
###############################################

# subset only astro-positive youngs

aj_rare_pos_y <- subset_samples(aj_rare_pos, age == "Y")


# agglomorate taxa at genus level

aj_rare_pos_y_compositional <- microbiome::transform(aj_rare_pos_y, "compositional")
aj_pos_merged_y <- tax_glom(aj_rare_pos_y_compositional, taxrank="Genus") # agglomerate taxa to Genus level
#otu_table()   OTU Table:         [ 370 taxa and 61 samples ]


df_aj<- pssd2veg(aj_pos_merged_y)

pred_var_astro<- as.data.frame(df_aj %>% dplyr::select(Nr_ST,AA_pdist_mean, ST1, ST2, ST3,ST4,ST5, ST6, ST7, ST9, ST10))
ASV_matrix_astro<- psotu2veg(aj_pos_merged_y)


######### access taxonomy ####
sample_data(aj_pos_merged_y)$SampleID_total<-make.names(sample_data(aj_pos_merged_y)$SampleID_total)
sample_names(aj_pos_merged_y) <- sample_data(aj_pos_merged_y)$SampleID_total

# phyloseq to data frame
Aj_phylo_na_df_astro<- phyloseq_to_df(aj_pos_merged_y)
Aj_phylo_na_table_astro<-phyloseq_to_df(aj_pos_merged_y)

#apply transformation
ASV_matrix_astro_trans <- disttransform(ASV_matrix_astro, method="chord")

###CCAs
set.seed(1)
cca_aj_astro <- cca(ASV_matrix_astro_trans ~ Nr_ST+AA_pdist_mean+ST1+ST2+ST3+ST4+ST5+ST6+ST7+ST9+ST10, scale=FALSE, data=df_aj)
cca_aj_astro


# Testing the significance of the CCA model:
anova.cca(cca_aj_astro, permutations = how(nperm = 999)) # by anova   
# 0.036 *


# Importance of components (axis)
summary(cca_aj_astro)$concont

# Importance of estimates
envfit_m <-envfit(cca_aj_astro,pred_var_astro, na.rm=TRUE, permutations=10000)
envfit_m$vectors




######################################################################################################
############################# Repeat with subset of variables ########################################
######################################################################################################


df_aj<- pssd2veg(aj_pos_merged_y)

pred_var_astro<- as.data.frame(df_aj %>% dplyr::select(AA_pdist_mean, ST6, ST10))
ASV_matrix_astro<- psotu2veg(aj_pos_merged_y)

######### access taxonomy ####
sample_data(aj_pos_merged_y)$SampleID_total<-make.names(sample_data(aj_pos_merged_y)$SampleID_total)
sample_names(aj_pos_merged_y) <- sample_data(aj_pos_merged_y)$SampleID_total

# phyloseq to data frame
Aj_phylo_na_df_astro<- phyloseq_to_df(aj_pos_merged_y)
Aj_phylo_na_table_astro<-phyloseq_to_df(aj_pos_merged_y)

# apply transformation
ASV_matrix_astro_trans <- disttransform(ASV_matrix_astro, method="chord")

# CCA
set.seed(1)
cca_aj_astro <- cca(ASV_matrix_astro_trans ~  AA_pdist_mean + ST6 + ST10, scale=FALSE, data=df_aj)
cca_aj_astro


# Testing the significance of the CCA model:
set.seed(1)
anova.cca(cca_aj_astro, permutations = how(nperm = 999)) # by anova    
#  0.019 *

# Importance of components
summary(cca_aj_astro)$concont



# Importance of estimates

########################################
### Stats for Supplementary Table S13 ##
########################################

set.seed(1)
envfit_m <-envfit(cca_aj_astro,pred_var_astro, na.rm=TRUE, permutations=10000)
envfit_m$vectors


# Scaling 1: species (ASVs) scores scaled to relative eigenvalues, sites (individuals) are weighted averages of the species (ASVs)  
plot(cca_aj_astro, scaling=1, display=c('sp', 'lc', 'cn'))



#### Get data for plotting in the correct shape

cca_aj_astro # all info is stored in this object as separate list

vare_spp_sco <- scores(cca_aj_astro, scaling=1, display = "species")
vare_sam_sco <- scores(cca_aj_astro, scaling=1, display = "sites")
vare_env_sco <- scores(cca_aj_astro, scaling=1, display = "bp")[1:3,]    # MHC div estimates

vare_spp_tbl <- as_tibble(vare_spp_sco)
vare_sam_tbl <- as_tibble(vare_sam_sco)
vare_env_tbl <- as_tibble(vare_env_sco)


# scaling is not correct for the environmental variables 
# (the plot.cca function transforms this depending on the used scaling method)
# with ordiArrowMul() we can find out which arrow multiplier is used to get the arrows of the current plot

ef <- envfit(cca_aj_astro ~ ., pred_var_astro)
plot(cca_aj_astro, scaling=1, display=c('sp', 'lc', 'cn'))
ordiArrowMul(scores(ef, display="vectors")) 
# arrow multiplier is 4.418684

# now use this multiplier on the extracted axis scores
rescaled_astro <- vare_env_tbl %>% 
  select(CCA1, CCA2) %>%
  as.matrix() * 6

vare_spp_tbl <- mutate(vare_spp_tbl, vgntxt=rownames(vare_spp_sco), ccatype = "species")
vare_sam_tbl <- mutate(vare_sam_tbl, vgntxt=rownames(vare_sam_sco), ccatype = "sites")
vare_env_tbl <- mutate(vare_env_tbl, vgntxt=rownames(vare_env_sco), ccatype = "bp")

# save files and combine
#write.csv(vare_spp_tbl, "vare_spp_tbl_pos_y.csv")
#write.csv(vare_sam_tbl,"vare_sam_tbl_pos_y.csv")
#write.csv(vare_env_tbl, "vare_env_tbl_pos_y.csv")
#write.csv(rescaled_astro, "rescaled_astro_pos_y.csv")


# Remove X before ASVs that start with a number and import again as "vare_spp_tbl", so we can join with taxonomy info

vare_spp_tbl # taxa names that were found in the CCA

taxonomy <- tax_table(aj_rare) [, 1:6]
rownames(taxonomy)

# make names a column
taxonomy1$names <- rownames(taxonomy)


# we join everything, join will introduce NAs for ASVs found in the first table but not in the second, we can remove these later 
cca_taxonomy <- left_join(taxonomy1, vare_spp_tbl, by="vgntxt")

# remove NAs (these are ASVs that are not present in the merged object we creasted for the CCA)
cca_taxonomy1 <- subset(cca_taxonomy, CCA1 != "NA")
names(cca_taxonomy1)

# check everything was joined correctly
head(cca_taxonomy1)
tail(cca_taxonomy1)

vare_spp_tbl <- cca_taxonomy1

# change order of columns 
vare_spp_tbl <- vare_spp_tbl[,c(10,11,1,12,3,4,5,6,7,8,9)]

# save file



#################################################
############### Final Figure 6B #################
#################################################


######## Plot with manually organized df

# final plot

cca_data <- read.csv("cca_data_pos_y.csv")
head(cca_data)

# for plotting microbiota as dots
data = filter(cca_data, ccatype=="species")
str(data)
data$CCA1 <- as.numeric(as.character(data$CCA1))
data$CCA2 <- as.numeric(as.character(data$CCA2))

# for plotting Diversity estimates
data1=filter(cca_data, ccatype=="bp")
data1$CCA1 <- as.numeric(as.character(data1$CCA1))
data1$CCA2 <- as.numeric(as.character(data1$CCA2))


cca_plot <- ggplot() + 
  theme_bw(base_size = 14) +
  geom_point(data=data, aes(x = CCA1, y = CCA2, fill = vgntxt)) + 
  geom_vline(aes(xintercept = 0), linetype = "dashed", alpha = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", alpha = 0.5) + # this will create the dashed coordinates, makes the plot more cca-like
  geom_label_repel(data=data[which(data$plot==1),], aes(x = CCA1, y = CCA2,label=Genus),
                   segment.color = "transparent", nudge_y = .01,force=5, size=rel(4.5), segment.size =0, parse=FALSE)+
  geom_segment(data=data1, aes(x=0, y=0, xend=CCA1 , yend=CCA2), size=1.1, arrow=arrow(length = unit(0.3,"cm")), color=c("red","red","red")) +
  ylab("CCA2 (34.58 %)")+ 
  xlab("CCA1 (38.12 %)")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits=c(-6,6), breaks=seq(-6,6, by=2)) +
  scale_y_continuous(limits=c(-6,4), breaks=seq(-6,4, by=2))+ theme (legend.position = "none") 

cca_plot_1  <- cca_plot + guides(fill=FALSE)+ guides(size=FALSE)
cca_plot_1 

cca_plot_final_young_pos <- cca_plot_1 +
  geom_label_repel(data=data1, aes(x=CCA1, y=CCA2), label=c(expression(Mean[pdist]), "ST6", "ST10") ,nudge_x=-.1, nudge_y=-.1, size=4, color = "red")

cca_plot_final_young_pos



##############################
####### Final Figure 6 #######
##############################

# arrange panel A and B

Figure6 <- ggarrange(cca_plot_final_adult_pos, cca_plot_final_young_pos, labels=c("A","B"),nrow=2,ncol=1)
Figure6


# note that changing the plot window size changes label position
# for label position as in publication
# export as .jpeg with width 625 and height 925






####################################################################################################################################
############################################################## Ancom ###############################################################
####################################################################################################################################

# load R functions

library(exactRankTests)
library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
source("ancom_v2.1.R")
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5682008/
#https://github.com/FrederickHuangLin/ANCOM


# Load phyloseq object
# this is the unrarefied microbiome data filtered to taxa present in at least 30% of samples

Ancom_phylo <- readRDS("Ancom_phylo.RDS")

# Separate young and adult bats

Ancom_phylo_y <- subset_samples(Ancom_phylo, age == "Y")
#   [ 753 taxa and 93 samples ]

Ancom_phylo_a <- subset_samples(Ancom_phylo, age == "A")
#   [ 753 taxa and 131 samples ]



####################################################################################################################
################################################ Young bats: POS / NEG #############################################
####################################################################################################################

otu_data <- data.frame(otu_table(Ancom_phylo_y))

otu_data <- tibble::rownames_to_column(otu_data, "feature_id")
otu_id = otu_data$feature_id
otu_data = data.frame(otu_data[, -1], check.names = FALSE)
rownames(otu_data) = otu_id

meta_data <- sample_data(Ancom_phylo_y)

colnames(meta_data)[1] <- "Sample.ID"


# Step 1: Data preprocessing

feature_table = otu_data; sample_var = "Sample.ID"; group_var = "astro2"
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = TRUE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "astro2"; p_adj_method = "BH"; alpha = 0.05
adj_formula =  "landscape+sex"; rand_formula = NULL
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)

res$out
list <- res$out
list0.9 <- subset(list, detected_0.9 == "TRUE")
list0.8 <- subset(list, detected_0.8 == "TRUE")
list0.7 <- subset(list, detected_0.7 == "TRUE")
list0.6 <- subset(list, detected_0.6 == "TRUE")

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off_y = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off_y) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off_y["detected_0.8"], label = "W[0.8]")

fig = res$fig +  
  geom_hline(yintercept =cut_off_y["detected_0.6"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig 


ggsave("ancom_astro_y0.6.pdf")

# save "fig" and "res" results
# add taxonomic info of ASVs
# these were used for Supplementary Table S15


#####################################
######### Final Figure 7B ###########
#####################################

ancom_astro_y_tax <- read.csv("ancom_astro_y_tax.csv")


# plot with Family name
library(ggrepel)

volcano_young1 <- ggplot(ancom_astro_y_tax, aes(x = clr, y=W_fig, label = "Family"))+
  geom_point(data=ancom_astro_y_tax %>% filter(zero_ind=="No"))+
  geom_point(data=ancom_astro_y_tax %>% filter(detected_0.6=="TRUE") %>% filter(zero_ind=="No"), aes(x = clr, y=W_fig, color=Phylum, shape =astro), size=2.5)+
  geom_hline(yintercept = cut_off_y["detected_0.6"], linetype = "dashed") +
  geom_text_repel(data=ancom_astro_y_tax %>% filter(detected_0.6=="TRUE") %>% filter(zero_ind=="No"),
                  # Filter data first
                  aes(label=Family, color = Phylum), size=rel(5), force = 5)+
  scale_color_manual(guide=guide_legend(ncol=1), values = c("forestgreen","#8738de","#00569c"))+
  scale_shape_manual(guide=guide_legend(ncol=1), name= "More abundant in", labels=c("AstV-","AstV+"), values=c(17,15))+   
  ylim(0,800) +
  xlim(-2,2.5) +
  theme_bw()+
  labs(x="CLR mean difference", y="W statistic")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text = element_text(size = 17))+
  guides(fill=FALSE)+theme(legend.position = "right")

volcano_young2 <- volcano_young1 + guides(shape = guide_legend(override.aes = list(size = 3))) # increase size of the shape legend symbols
volcano_young2


############### Boxplot of ASVs identified in Ancom at 60% ###############
############### Boxplot of ASVs identified in Ancom at 60% ###############


## Add ASVs found at 60% in ancom

aj_tax <- tax_table(Ancom_phylo_y)

names <- rownames(aj_tax)
test<- as.matrix(cbind(aj_tax, names))
tax_table(Ancom_phylo_y) <- test # map back

my_ancom_asvs <- subset_taxa(Ancom_phylo_y, names=="New.ReferenceOTU40"|
                               names=="112997"|
                               names=="76821"|
                               names=="New.ReferenceOTU298"|
                               names=="New.CleanUp.ReferenceOTU7477"|
                               names=="4451477"|
                               names=="325419"|
                               names=="304779"|
                               names=="342666")


my_ancom_asvs.df <- data.frame(otu_table(my_ancom_asvs))

# add +1 to all values in the dataframe, then we are able to compare abundances later properly with log transformation
my_ancom_asvs.df_1 <- my_ancom_asvs.df + 1  # first save into new object, double check is correct
my_ancom_asvs.df <- my_ancom_asvs.df_1      # then rename

my_ancom_asvs.df  = cbind(as(my_ancom_asvs.df, "data.frame"), as(tax_table(my_ancom_asvs)[rownames(my_ancom_asvs.df), ], "matrix"))

my_ancom_asvs.df = subset(my_ancom_asvs.df, select = -c(Kingdom,Phylum, Class, Order, Genus, Species, names))
my_ancom_asvs.df = subset(my_ancom_asvs.df, select = -c(Family) )

my_ancom_asvs.df<-data.frame(t(my_ancom_asvs.df))
head(my_ancom_asvs.df)
tail(my_ancom_asvs.df)

# check if rows are ordered similar in metadata
metadata<-data.frame(sample_data(my_ancom_asvs))
head(metadata)
tail(metadata)

# both metadata and asv-df are in the same order, thus we can bind together
metadata1<-cbind(metadata, my_ancom_asvs.df)

# Change colnames of ASV columns to correct and short names
colnames(metadata1)[colnames(metadata1) %in% c("New.ReferenceOTU40", "X112997", "X76821", "New.ReferenceOTU298","New.CleanUp.ReferenceOTU7477","X4451477","X325419","X304779","X342666")] <- c("ASV40", "ASV112997","ASV76821","ASV298","ASV7477","ASV4451477","ASV325419","ASV304779","ASV342666")

#convert to long format
data_long <- gather(metadata1, ASV_ID, Abundance,  ASV40:ASV342666, factor_key=TRUE)

table(data_long$ASV_ID) # check it worked

data_long$Abundance <- as.numeric(data_long$Abundance)
#write.csv(data_long, "ancom_data_long_y.csv")


data_long <- read.csv("ancom_data_long_y.csv")



#################################################
########### Supplementary Figure S8B ############
#################################################

data_long$ASV_ID

# New facet label names for supp variable

supp.labs.y <- c("Helicobacteraceae 
ASV40", "Helicobacteraceae
ASV112997","Mycoplasmataceae
ASV76821","Mycoplasmataceae
ASV298","Mycoplasmataceae
ASV7477","Clostridiaceae
ASV4451477","Clostridiaceae
ASV325419","Clostridiaceae
ASV342666","Clostridiaceae
ASV304779")

names(supp.labs.y) <- c("ASV40", "ASV112997","ASV76821","ASV298","ASV7477","ASV4451477","ASV325419","ASV342666","ASV304779")


# change order of asvs in a new df

neworder <- c("ASV304779","ASV325419","ASV342666","ASV4451477",
              "ASV112997","ASV40","ASV298","ASV7477","ASV76821")

library(dplyr)  ## or dplyr (transform -> mutate)
data_long1 <- arrange(mutate(data_long,ASV_ID=factor(ASV_ID,levels=neworder)),ASV_ID)

all_asvs_young_log <- ggplot(data_long1, aes(x = astro, y = Abundance))+
  geom_boxplot(outlier.shape = NA, aes(fill=astro))+
  scale_fill_manual(name = "",labels= c((expression(AstV^"-")),c(expression(AstV^"+"))), values=c("dodgerblue2", "coral3"))+  
  geom_jitter(width=0.05, size = 1, alpha = 0.6)+
  facet_wrap(~ASV_ID, scales = "free", ncol = 4, labeller = labeller(ASV_ID = supp.labs.y))+
  scale_y_log10()+
  theme_bw(base_size=16)+
  stat_compare_means(method = "wilcox.test", size=4.8, aes(label = sprintf("Wilcoxon, p = %5.4f", as.numeric(..p.format..))))+
  labs(x= " ", y =expression(Log[10] * " ASV abundance")) +
  my_theme4 +
  theme(legend.position=c(0.6,0.15))+
  theme(legend.text = element_text(size = 17))+ggtitle("B")+theme(strip.text.x = element_text(face = "bold.italic"))



## Subset the most convincing ASVS

data_long <- read.csv("ancom_data_long_y.csv")


## Subset the most convincing ASVS

data_long_otu40 <- subset(data_long, ASV_ID =="ASV40")
data_long_otu112997 <- subset(data_long, ASV_ID =="ASV112997")
data_long_otu76821 <- subset(data_long, ASV_ID =="ASV76821")
data_long_otu298 <- subset(data_long, ASV_ID =="ASV298")
data_long_otu7477<- subset(data_long, ASV_ID == "ASV7477")
data_long_otu4451477<- subset(data_long, ASV_ID == "ASV4451477")
data_long_otu325419<- subset(data_long, ASV_ID == "ASV325419")
data_long_otu304779<- subset(data_long, ASV_ID == "ASV304779")
data_long_otu342666 <- subset(data_long, ASV_ID =="ASV342666")


##lets split it by ST1-10

data_long_otu76821$ST1 <- as.factor(data_long_otu76821$ST1)
data_long_otu76821$ST2 <- as.factor(data_long_otu76821$ST2)
data_long_otu76821$ST3 <- as.factor(data_long_otu76821$ST3)
data_long_otu76821$ST4 <- as.factor(data_long_otu76821$ST4)
data_long_otu76821$ST5 <- as.factor(data_long_otu76821$ST5)
data_long_otu76821$ST6 <- as.factor(data_long_otu76821$ST6)
data_long_otu76821$ST7 <- as.factor(data_long_otu76821$ST7)
data_long_otu76821$ST8 <- as.factor(data_long_otu76821$ST8)
data_long_otu76821$ST9 <- as.factor(data_long_otu76821$ST9)
data_long_otu76821$ST10 <- as.factor(data_long_otu76821$ST10)


Fig8_4 <-ggplot(data_long_otu76821, aes(x = ST2, y = Abundance, fill=ST2, color=ST2))+
  geom_boxplot(alpha= 0.7, outlier.shape = NA)+
  scale_fill_manual(name = "ST", labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+ 
  scale_color_manual(name = "ST", labels = c("Absent", "Present"),values=c("#c29844","#6ab166"))+
  geom_point(size= 1.5, position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size=4.6)+
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "italic"))+
  ggtitle("Mycoplasmataceae",subtitle = "ASV76821")+
  scale_y_log10()+
  labs(y = expression(Log[10]~"abundance")) +
  xlab("ST2")+
  theme(legend.position="right", axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(legend.title = element_text(size = 16), legend.text = element_text(size=14))  

Fig8_5 <-ggplot(data_long_otu76821, aes(x = ST4, y = Abundance, fill=ST4, color=ST4))+
  geom_boxplot(alpha= 0.7, outlier.shape = NA)+
  scale_fill_manual(name = "ST", labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+ 
  scale_color_manual(name = "ST", labels = c("Absent", "Present"),values=c("#c29844","#6ab166"))+ 
  geom_point(size= 1.5, position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size=4.6)+
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "italic"))+
  ggtitle("Mycoplasmataceae",subtitle = "ASV76821")+
  scale_y_log10()+
  ylab(" ")+ 
  xlab("ST4")+
  theme(legend.position="right", axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(legend.title = element_text(size = 16), legend.text = element_text(size=14))  




##lets split it by ST1-10

data_long_otu342666$ST1 <- as.factor(data_long_otu342666$ST1)
data_long_otu342666$ST2 <- as.factor(data_long_otu342666$ST2)
data_long_otu342666$ST3 <- as.factor(data_long_otu342666$ST3)
data_long_otu342666$ST4 <- as.factor(data_long_otu342666$ST4)
data_long_otu342666$ST5 <- as.factor(data_long_otu342666$ST5)
data_long_otu342666$ST6 <- as.factor(data_long_otu342666$ST6)
data_long_otu342666$ST7 <- as.factor(data_long_otu342666$ST7)
data_long_otu342666$ST8 <- as.factor(data_long_otu342666$ST8)
data_long_otu342666$ST9 <- as.factor(data_long_otu342666$ST9)
data_long_otu342666$ST10 <- as.factor(data_long_otu342666$ST10)


Fig8_6 <- ggplot(data_long_otu342666, aes(x = ST5, y = Abundance, fill=ST5, color=ST5))+
  geom_boxplot(alpha= 0.7, outlier.shape = NA)+
  scale_fill_manual(name = "ST", labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+ 
  scale_color_manual(name = "ST", labels = c("Absent", "Present"),values=c("#c29844","#6ab166"))+
  geom_point(size= 1.5, position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size=4.6)+
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "italic"))+
  ggtitle("Clostridiaceae",subtitle = "ASV342666")+
  scale_y_log10()+
  ylab(" ")+
  xlab("ST5")+
  theme(legend.position="right", axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(legend.title = element_text(size = 16), legend.text = element_text(size=14))  



##lets split it by ST1-10
data_long_otu7477$ST1 <- as.factor(data_long_otu7477$ST1)
data_long_otu7477$ST2 <- as.factor(data_long_otu7477$ST2)
data_long_otu7477$ST3 <- as.factor(data_long_otu7477$ST3)
data_long_otu7477$ST4 <- as.factor(data_long_otu7477$ST4)
data_long_otu7477$ST5 <- as.factor(data_long_otu7477$ST5)
data_long_otu7477$ST6 <- as.factor(data_long_otu7477$ST6)
data_long_otu7477$ST7 <- as.factor(data_long_otu7477$ST7)
data_long_otu7477$ST8 <- as.factor(data_long_otu7477$ST8)
data_long_otu7477$ST9 <- as.factor(data_long_otu7477$ST9)
data_long_otu7477$ST10 <- as.factor(data_long_otu7477$ST10)


Fig8_7 <- ggplot(data_long_otu7477, aes(x = ST2, y = Abundance, fill=ST2, color=ST2))+
  geom_boxplot(alpha= 0.7, outlier.shape = NA)+
  scale_fill_manual(name = "ST", labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+ 
  scale_color_manual(name = "ST", labels = c("Absent", "Present"),values=c("#c29844","#6ab166"))+
  geom_point(size= 1.5, position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size=4.6)+
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "italic"))+
  ggtitle("Mycoplasmataceae",subtitle = "ASV7477")+
  scale_y_log10()+
  labs(y = expression(Log[10]~"abundance")) +
  xlab("ST2")+
  theme(legend.position="right", axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(legend.title = element_text(size = 16), legend.text = element_text(size=14))  


Fig8_8 <- ggplot(data_long_otu7477, aes(x = ST4, y = Abundance, fill=ST4, color=ST4))+
  geom_boxplot(alpha= 0.7, outlier.shape = NA)+
  scale_fill_manual(name = "ST", labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+ 
  scale_color_manual(name = "ST", labels = c("Absent", "Present"),values=c("#c29844","#6ab166"))+ 
  geom_point(size= 1.5, position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size=4.6)+
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "italic"))+
  ggtitle("Mycoplasmataceae",subtitle = "ASV7477")+
  scale_y_log10()+
  ylab(" ")+
  xlab("ST4")+
  theme(legend.position="right", axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(legend.title = element_text(size = 16), legend.text = element_text(size=14))  



##lets split it by ST1-10

data_long_otu112997$ST1 <- as.factor(data_long_otu112997$ST1)
data_long_otu112997$ST2 <- as.factor(data_long_otu112997$ST2)
data_long_otu112997$ST3 <- as.factor(data_long_otu112997$ST3)
data_long_otu112997$ST4 <- as.factor(data_long_otu112997$ST4)
data_long_otu112997$ST5 <- as.factor(data_long_otu112997$ST5)
data_long_otu112997$ST6 <- as.factor(data_long_otu112997$ST6)
data_long_otu112997$ST7 <- as.factor(data_long_otu112997$ST7)
data_long_otu112997$ST8 <- as.factor(data_long_otu112997$ST8)
data_long_otu112997$ST9 <- as.factor(data_long_otu112997$ST9)
data_long_otu112997$ST10 <- as.factor(data_long_otu112997$ST10)


Fig8_9 <- ggplot(data_long_otu112997, aes(x = ST9, y = Abundance, fill=ST9, color=ST9))+
  geom_boxplot(alpha= 0.7, outlier.shape = NA)+
  scale_fill_manual(name = "ST", labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+ 
  scale_color_manual(name = "ST", labels = c("Absent", "Present"),values=c("#c29844","#6ab166"))+ 
  geom_point(size= 1.5, position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size=4.6)+
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "italic"))+
  ggtitle("") + 
  ggtitle("Helicobacteraceae",subtitle = "ASV112997")+
  scale_y_log10() +
  ylab(" ")+
  xlab("ST9")+
  theme(legend.position="right", axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(legend.title = element_text(size = 16), legend.text = element_text(size=14))  



## lets split it by ST1-10
## mycoplasma

data_long_otu298$ST1 <- as.factor(data_long_otu298$ST1)
data_long_otu298$ST2 <- as.factor(data_long_otu298$ST2)
data_long_otu298$ST3 <- as.factor(data_long_otu298$ST3)
data_long_otu298$ST4 <- as.factor(data_long_otu298$ST4)
data_long_otu298$ST5 <- as.factor(data_long_otu298$ST5)
data_long_otu298$ST6 <- as.factor(data_long_otu298$ST6)
data_long_otu298$ST7 <- as.factor(data_long_otu298$ST7)
data_long_otu298$ST8 <- as.factor(data_long_otu298$ST8)
data_long_otu298$ST9 <- as.factor(data_long_otu298$ST9)
data_long_otu298$ST10 <- as.factor(data_long_otu298$ST10)


Fig8_10 <- ggplot(data_long_otu298, aes(x = ST2, y = Abundance, fill=ST2, color=ST2))+
  geom_boxplot(alpha= 0.7, outlier.shape = NA)+
  scale_fill_manual(name = "ST", labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+ 
  scale_color_manual(name = "ST", labels = c("Absent", "Present"),values=c("#c29844","#6ab166"))+
  geom_point(size= 1.5, position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size=4.6)+
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "italic"))+
  ggtitle("Mycoplasmataceae",subtitle = "ASV298")+
  scale_y_log10()+
  labs(y = expression(Log[10]~"abundance")) +
  xlab("ST2")+
  theme(legend.position="right", axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(legend.title = element_text(size = 16), legend.text = element_text(size=14))  



## OTU40

data_long_otu40$ST1 <- as.factor(data_long_otu40$ST1)
data_long_otu40$ST2 <- as.factor(data_long_otu40$ST2)
data_long_otu40$ST3 <- as.factor(data_long_otu40$ST3)
data_long_otu40$ST4 <- as.factor(data_long_otu40$ST4)
data_long_otu40$ST5 <- as.factor(data_long_otu40$ST5)
data_long_otu40$ST6 <- as.factor(data_long_otu40$ST6)
data_long_otu40$ST7 <- as.factor(data_long_otu40$ST7)
data_long_otu40$ST8 <- as.factor(data_long_otu40$ST8)
data_long_otu40$ST9 <- as.factor(data_long_otu40$ST9)
data_long_otu40$ST10 <- as.factor(data_long_otu40$ST10)

Fig8_11 <- ggplot(data_long_otu40 , aes(x = ST4, y = Abundance, fill=ST4, color=ST4))+
  geom_boxplot(alpha= 0.7, outlier.shape = NA)+
  scale_fill_manual(name = "ST", labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+ 
  scale_color_manual(name = "ST", labels = c("Absent", "Present"),values=c("#c29844","#6ab166"))+
  geom_point(size= 1.5, position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size=4.6)+
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "italic"))+
  ggtitle("Helicobacteraceae",subtitle = "ASV40")+
  scale_y_log10()+
  ylab(" ")+
  xlab("ST4")+
  theme(legend.position="right", axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(legend.title = element_text(size = 16), legend.text = element_text(size=14))  



##lets split it by ST1-10
data_long_otu4451477$ST1 <- as.factor(data_long_otu4451477$ST1)
data_long_otu4451477$ST2 <- as.factor(data_long_otu4451477$ST2)
data_long_otu4451477$ST3 <- as.factor(data_long_otu4451477$ST3)
data_long_otu4451477$ST4 <- as.factor(data_long_otu4451477$ST4)
data_long_otu4451477$ST5 <- as.factor(data_long_otu4451477$ST5)
data_long_otu4451477$ST6 <- as.factor(data_long_otu4451477$ST6)
data_long_otu4451477$ST7 <- as.factor(data_long_otu4451477$ST7)
data_long_otu4451477$ST8 <- as.factor(data_long_otu4451477$ST8)
data_long_otu4451477$ST9 <- as.factor(data_long_otu4451477$ST9)
data_long_otu4451477$ST10 <- as.factor(data_long_otu4451477$ST10)


Fig8_12 <- ggplot(data_long_otu4451477, aes(x = ST2, y = Abundance, fill=ST2, color=ST2))+
  geom_boxplot(alpha= 0.7, outlier.shape = NA)+
  scale_fill_manual(name = "ST", labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+ 
  scale_color_manual(name = "ST", labels = c("Absent", "Present"),values=c("#c29844","#6ab166"))+ 
  geom_point(size= 1.5, position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size=4.6)+
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "italic"))+
  ggtitle("") + 
  ggtitle("Clostridiaceae",subtitle = "ASV4451477")+
  scale_y_log10()+
  ylab(" ")+
  xlab("ST2")+
  theme(legend.position="right", axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(legend.title = element_text(size = 16), legend.text = element_text(size=14))  







####################################################################################################################
################################################ Adult bats: POS / NEG #############################################
####################################################################################################################

otu_data <- data.frame(otu_table(Ancom_phylo_a))

otu_data <- tibble::rownames_to_column(otu_data, "feature_id")
otu_id = otu_data$feature_id
otu_data = data.frame(otu_data[, -1], check.names = FALSE)
rownames(otu_data) = otu_id

meta_data <- sample_data(Ancom_phylo_a)

colnames(meta_data)[1] <- "Sample.ID"


# Step 1: Data preprocessing

feature_table = otu_data; sample_var = "Sample.ID"; group_var = "astro2"
out_cut = 0.05; zero_cut = 0.90; lib_cut = 1000; neg_lb = TRUE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "astro2"; p_adj_method = "BH"; alpha = 0.05
adj_formula =  "landscape+sex"; rand_formula = NULL
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)

res$out
list <- res$out
list0.9 <- subset(list, detected_0.9 == "TRUE")
list0.8 <- subset(list, detected_0.8 == "TRUE")
list0.7 <- subset(list, detected_0.7 == "TRUE")
list0.6 <- subset(list, detected_0.6 == "TRUE")

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off_a = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off_a) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off_a["detected_0.6"], label = "W[0.6]")

fig = res$fig +  
  geom_hline(yintercept = cut_off_a["detected_0.6"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig 


# save
ggsave("ancom_astro_a0.6.pdf")

# save "fig" and "res" results
# add taxonomic info of ASVs
# these were used for Supplementary Table S14


#####################################
######### Final Figure 7A ###########
#####################################

ancom_astro_a_tax <- read.csv("ancom_astro_a_tax.csv")

# plot 
volcano_adult3 <- ggplot(ancom_astro_a_tax, aes(x = clr, y=W_fig, label = "Family"))+
  geom_point(data=ancom_astro_a_tax %>% filter(zero_ind=="No"))+
  geom_point(data=ancom_astro_a_tax %>% filter(detected_0.6=="TRUE") %>% filter(zero_ind=="No"), aes(x = clr, y=W_fig, color=Phylum, shape =astro), size=2.5)+
  geom_hline(yintercept = cut_off_a["detected_0.6"], linetype = "dashed") +
  geom_text_repel(data=ancom_astro_a_tax %>% filter(detected_0.6=="TRUE") %>% filter(zero_ind=="No") ,
                  # Filter data first
                  aes(label=Family, color = Phylum), size=rel(5), force = 5)+
  scale_color_manual(values = c("forestgreen", "#8738de","#00569c"))+
  scale_shape_manual(name= "More abundant in", labels=c("AstV-","AstV+"), values=c(17,15))+  
  ylim(0,800) +
  xlim(-2,2.5) +
  theme_bw()+
  labs(x="CLR mean difference", y="W statistic")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text = element_text(size = 17))+
  guides(fill=FALSE)+theme(legend.position = "none")

volcano_adult3 



#####################################
######### Final Figure 7 ############
#####################################

library(cowplot)

# extract the legend from one of the plots
legend_young <- get_legend(volcano_young1 + theme(legend.position="bottom")+ guides(guide=guide_legend(ncol=1)))


# arrange the plots in a single row
volcano_row <- plot_grid(volcano_adult3 + theme(legend.position="none"),
                         volcano_young2 + theme(legend.position="none"),
                         align = 'vh',
                         labels = c("A", "B"),
                         nrow = 1)

# combine
Figure7 <- plot_grid(volcano_row, legend_young, ncol = 1, rel_heights = c(1, .2))
Figure7




############### Boxplot of ASVs identified in Ancom at 60% ###############
############### Boxplot of ASVs identified in Ancom at 60% ###############


## Add ASVs found at 60% in ancom

aj_tax <- tax_table(Ancom_phylo_a)

names <- rownames(aj_tax)
test<- as.matrix(cbind(aj_tax, names))
tax_table(Ancom_phylo_a) <- test # map back

my_ancom_asvs <- subset_taxa(Ancom_phylo_a, names=="New.ReferenceOTU40"|
                               names=="247639"|
                               names=="1141646"|
                               names=="New.ReferenceOTU365"|
                               names=="712047"|
                               names=="1073436"|
                               names=="730049"|
                               names=="813945")


my_ancom_asvs.df <- data.frame(otu_table(my_ancom_asvs))

# add +1 to all values in the dataframe, then we are able to compare abundances later properly 
my_ancom_asvs.df_1 <- my_ancom_asvs.df + 1
my_ancom_asvs.df <- my_ancom_asvs.df_1 # map back

my_ancom_asvs.df  = cbind(as(my_ancom_asvs.df, "data.frame"), as(tax_table(my_ancom_asvs)[rownames(my_ancom_asvs.df), ], "matrix"))

my_ancom_asvs.df = subset(my_ancom_asvs.df, select = -c(Kingdom, Phylum, Class, Order, Family, Genus, Species, names, names.1, names.2))

my_ancom_asvs.df<-data.frame(t(my_ancom_asvs.df))
head(my_ancom_asvs.df)
tail(my_ancom_asvs.df)

# extract metadata
metadata<-data.frame(sample_data(my_ancom_asvs))
head(metadata)
tail(metadata)

metadata1<-cbind(metadata, my_ancom_asvs.df)

# Change colnames of ASV columns
colnames(metadata1)[colnames(metadata1) %in% c("New.ReferenceOTU40", "X247639","X1141646","New.ReferenceOTU365","X712047","X1073436","X730049","X813945")] <- c("ASV40", "ASV247639","ASV1141646","ASV365","ASV712047","ASV1073436","ASV730049","ASV813945")


#convert to long format
data_long <- gather(metadata1, ASV_ID, Abundance, ASV40:ASV813945, factor_key=TRUE)

table(data_long$ASV_ID)

data_long$Abundance <- as.numeric(data_long$Abundance)
#write.csv(data_long,"ancom_data_long_a.csv")

data_long <- read.csv("ancom_data_long_a.csv")


#################################################
########### Supplementary Figure S8A ############
#################################################

data_long$ASV_ID

# New facet label names for supp variable

supp.labs.a <- c("Helicobacteraceae
ASV40", "Clostridiaceae
ASV247639","Streptococcaceae
ASV1141646","Clostridiaceae
ASV365","Clostridiaceae
ASV712047","Staphylococcaceae
ASV1073436","Staphylococcaceae
ASV730049","Pseudomonadaceae
ASV813945")

names(supp.labs.a) <- c("ASV40", "ASV247639","ASV1141646","ASV365","ASV712047","ASV1073436","ASV730049","ASV813945")


# plot on log-scale

all_asvs_adult_log <- ggplot(data_long, aes(x = astro, y = Abundance))+
  geom_boxplot(outlier.shape = NA, aes(fill=astro))+
  scale_fill_manual(name = "",labels= c((expression(AstV^"-")),c(expression(AstV^"+"))), values=c("dodgerblue2", "coral3"))+  
  geom_jitter(width=0.05, size = 1, alpha = 0.6)+
  facet_wrap(~ASV_ID, scales = "free", ncol = 4, labeller = labeller(ASV_ID = supp.labs.a))+
  scale_y_log10()+
  theme_bw(base_size=16)+
  stat_compare_means(method = "wilcox.test", size=4.8, aes(label = sprintf("Wilcoxon, p = %5.4f", as.numeric(..p.format..))))+
  labs(x= " ", y =expression(Log[10] * " ASV abundance")) +
  my_theme4 +
  theme(legend.position = "none")+ggtitle("A")+theme(strip.text.x = element_text(face = "bold.italic"))



#################################################
########### Supplementary Figure S8 #############
#################################################

ggarrange(all_asvs_adult_log,all_asvs_young_log, nrow=2, ncol=1, heights=c(2,3))





## Subset the most convincing ASVS

data_long <- read.csv("ancom_data_long_a.csv")


data_long_otu40 <- subset(data_long, ASV_ID == "ASV40")
data_long_otu247639 <- subset(data_long, ASV_ID == "ASV247639")
data_long_otu1141646 <- subset(data_long, ASV_ID =="ASV1141646")
data_long_otu365 <- subset(data_long, ASV_ID =="ASV365")
data_long_otu712047 <- subset(data_long, ASV_ID =="ASV712047")
data_long_otu1073436 <- subset(data_long, ASV_ID =="ASV1073436")
data_long_otu730049 <- subset(data_long, ASV_ID == "ASV730049")
data_long_otu813945 <- subset(data_long, ASV_ID == "ASV813945")



##lets split it by ST1-10

data_long_otu730049$ST1 <- as.factor(data_long_otu730049$ST1)
data_long_otu730049$ST2 <- as.factor(data_long_otu730049$ST2)
data_long_otu730049$ST3 <- as.factor(data_long_otu730049$ST3)
data_long_otu730049$ST4 <- as.factor(data_long_otu730049$ST4)
data_long_otu730049$ST5 <- as.factor(data_long_otu730049$ST5)
data_long_otu730049$ST6 <- as.factor(data_long_otu730049$ST6)
data_long_otu730049$ST7 <- as.factor(data_long_otu730049$ST7)
data_long_otu730049$ST8 <- as.factor(data_long_otu730049$ST8)
data_long_otu730049$ST9 <- as.factor(data_long_otu730049$ST9)
data_long_otu730049$ST10 <- as.factor(data_long_otu730049$ST10)


Fig8_1  <- ggplot(data_long_otu730049, aes(x = ST1, y = Abundance, fill=ST1, color=ST1))+
  geom_boxplot(alpha= 0.7, outlier.shape = NA)+
  scale_fill_manual(name = "ST", labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+ 
  scale_color_manual(name = "ST", labels = c("Absent", "Present"),values=c("#c29844","#6ab166"))+
  geom_point(size= 1.5, position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size=4.6)+
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "italic"))+
  scale_y_log10()+
  ggtitle("Staphylococcaceae",subtitle = "ASV730049")+
  labs(y = expression(Log[10]~"abundance")) +
  xlab("ST1")+
  theme(legend.position="right", axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(legend.title = element_text(size = 16), legend.text = element_text(size=14))  



##lets split it by ST1-10

data_long_otu813945$ST1 <- as.factor(data_long_otu813945$ST1)
data_long_otu813945$ST2 <- as.factor(data_long_otu813945$ST2)
data_long_otu813945$ST3 <- as.factor(data_long_otu813945$ST3)
data_long_otu813945$ST4 <- as.factor(data_long_otu813945$ST4)
data_long_otu813945$ST5 <- as.factor(data_long_otu813945$ST5)
data_long_otu813945$ST6 <- as.factor(data_long_otu813945$ST6)
data_long_otu813945$ST7 <- as.factor(data_long_otu813945$ST7)
data_long_otu813945$ST8 <- as.factor(data_long_otu813945$ST8)
data_long_otu813945$ST9 <- as.factor(data_long_otu813945$ST9)
data_long_otu813945$ST10 <- as.factor(data_long_otu813945$ST10)


Fig8_2  <- ggplot(data_long_otu813945, aes(x = ST4, y = Abundance, fill=ST4, color=ST4))+
  geom_boxplot(alpha= 0.7, outlier.shape = NA)+
  scale_fill_manual(name = "ST", labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+ 
  scale_color_manual(name = "ST", labels = c("Absent", "Present"),values=c("#c29844","#6ab166"))+
  geom_point(size= 1.5, position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size=4.6)+
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "italic"))+
  scale_y_log10()+
  ggtitle("Pseudomonadaceae",subtitle = "ASV813945")+
  ylab(" ")+
  xlab("ST4")+
  theme(legend.position="right", axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(legend.title = element_text(size = 16), legend.text = element_text(size=14))  



##lets split it by ST1-10

data_long_otu1141646$ST1 <- as.factor(data_long_otu1141646$ST1)
data_long_otu1141646$ST2 <- as.factor(data_long_otu1141646$ST2)
data_long_otu1141646$ST3 <- as.factor(data_long_otu1141646$ST3)
data_long_otu1141646$ST4 <- as.factor(data_long_otu1141646$ST4)
data_long_otu1141646$ST5 <- as.factor(data_long_otu1141646$ST5)
data_long_otu1141646$ST6 <- as.factor(data_long_otu1141646$ST6)
data_long_otu1141646$ST7 <- as.factor(data_long_otu1141646$ST7)
data_long_otu1141646$ST8 <- as.factor(data_long_otu1141646$ST8)
data_long_otu1141646$ST9 <- as.factor(data_long_otu1141646$ST9)
data_long_otu1141646$ST10 <- as.factor(data_long_otu1141646$ST10)


Fig8_3  <- ggplot(data_long_otu1141646, aes(x = ST6, y = Abundance, fill=ST6, color=ST6))+
  geom_boxplot(alpha= 0.7, outlier.shape = NA)+
  scale_fill_manual(name = "ST", labels = c("Absent","Present"),values=c("#c29844","#6ab166"))+ 
  scale_color_manual(name = "ST", labels = c("Absent", "Present"),values=c("#c29844","#6ab166"))+
  geom_point(size= 1.5, position=position_jitterdodge(jitter.width = NULL, jitter.height = 0, dodge.width = 0.75))+
  stat_compare_means(method = "wilcox.test", aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))), size=4.6)+
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(face = "italic"))+
  scale_y_log10()+
  ggtitle("Streptococcaceae",subtitle = "ASV1141646")+
  ylab(" ") + 
  xlab("ST6")+
  theme(legend.position="right", axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(legend.title = element_text(size = 16), legend.text = element_text(size=14))  





#################################################
################ Final Figure 8 #################
#################################################

Figure8A <- ggarrange(Fig8_1,Fig8_2,Fig8_3, ncol=3, common.legend = TRUE, legend = "right", labels = c("A"))
Figure8B <- ggarrange(Fig8_4,Fig8_5,Fig8_6,Fig8_7,Fig8_8,Fig8_9,Fig8_10,Fig8_11,Fig8_12,ncol=3, nrow=3, common.legend = TRUE, legend = "right", labels = c("B"))

ggarrange(Figure8A,Figure8B, nrow=2, heights=c(1,3))
