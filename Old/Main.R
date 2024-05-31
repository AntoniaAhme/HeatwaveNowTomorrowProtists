### HNT22 mesocosms 18S rRNA metabarcoding data ###
## Antonia Ahme, 17.11.2023 ##
# Generally SRS before bargraphs + alpha diversity and CLR on raw data + euclidean PCA for ordination (Gloor 2017)

#### HOUSEKEEPING ####
require(dplyr)
require(tidyr)
require(vegan)
require(ggplot2)
require(phyloseq)
require(SRS)
require(zCompositions)
require(propr)
require(BiodiversityR)
require(ggrepel)
require(ggpubr)
require(rstatix)
require(easyCODA)
require(microbial)
require(ellipse)
require(qualpalr)
require(microViz)
require(writexl)
require(ggforce)

rm(list=ls())
options(stringsAsFactors = F)
setwd("~/AWI/RProjects/HNT22")

set.seed(22)

#### FUNCTIONS & LAYOUTS ####
treat_pal <- c("Ambient" = "skyblue2","Ambient+HW" = "darkblue", "RCP" = "goldenrod2", "RCP+HW" = "tomato2")

# Symbols for base R plotting
treat_pch <- c(16,17,15,19)

# Create the designs for plotting
bar.theme <- theme(panel.background = element_blank(),
                   panel.border = element_blank(),
                   panel.grid = element_blank(),
                   axis.title.x=element_blank(),
                   axis.text.x = element_text(size = 20, face= "bold"),
                   axis.ticks.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   plot.title = element_blank(),
                   strip.text = element_text(size = 15, face = "bold"),
                   strip.background = element_blank(), strip.placement = "outside",
                   text = element_text(size = 20, face = "bold"), legend.position = "right")

plot.theme <- theme(panel.background = element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    panel.grid = element_blank(),
                    axis.title.x = element_text(size = 15, face= "bold"),
                    axis.text.x = element_text(size = 15, face= "bold"),
                    axis.ticks.x = element_blank(),
                    axis.title.y = element_text(size = 15, face= "bold"),
                    axis.text.y = element_text(face="plain", color="black", size=15),
                    axis.ticks.y = element_blank(),
                    plot.title = element_text(size = 15, face= "bold"),
                    strip.text = element_text(size = 15, face = "bold"),
                    strip.background = element_blank(), strip.placement = "outside",
                    text = element_text(size = 15, face = "bold"), legend.position = "none")

# Standard error and mean function
se <- function(x, ...) sqrt(var(x, ...)/length(x))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = se(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#### UPLOAD DATA ####
# Import asv and tax tabs
asv <- read.delim('Data/Counts.txt')
asv <- data.frame(asv, row.names = 1)
colnames(asv) <- sub("HNT22.euk.", "", colnames(asv))
colnames(asv) <- gsub('\\.', '-', colnames(asv))
tax <- read.delim('Data/Taxonomy.txt')

# Import sample tab
sam <- read.delim('Data/Samples.txt')
sam$sample_ID <- sub("HNT22-euk-", "", sam$sample_ID)
sam$sample_ID <- as.factor(sam$sample_ID)

#### PREPARE DATA FOR PLOTTING ####
# Remove ASVs that don't occur in the dataset
asv <- asv[rowSums(asv)>0,]

# Investigate sequencing depth and remove samples where it is too low
depths <- colSums(asv) # prepare df containing the sequencing depth
plot(depths) # sequencing depth does not have any outliers or very low numbers
min(depths)
max(depths)
quantile(depths, probs = seq(0,1,.05)) # Show all sequencing depth quantiles
# Remove all samples that do not fall into the 90% quantile range
asv <- asv[,colSums(asv)<104236]
asv <- asv[,colSums(asv)>53000]
rm(depths)

# Match sam tab to new asv tab as some samples might have been removed in the last step
sam <- sam[sam$sample_ID %in% colnames(asv),]

# Check and adjust whether rownames of meta-info file match colnames of ASV counts
all(sam$sample_ID == colnames(asv))

# Remove asvs with a count of less than 10 (singletons) in replicate sample means (create rep means first)
ASV <- asv
colnames(ASV) <- sam$uni_ID[match(colnames(ASV),sam$sample_ID)]
Names <- unique(names(ASV))
ASV <- sapply(Names, function(x)  rowMeans(ASV[names(ASV) %in% x]))
ASV <- as.data.frame(ASV)
ASV <- ASV %>% filter_at(vars(1:48), any_vars(.>=10)) # number of pooled samples and cutoff-level of 10
ASV.rn <- rownames(ASV)
asv <- asv[rownames(asv) %in% ASV.rn,]
rm(ASV,Names,ASV.rn)
asv_all <- asv

tax <- tax[tax$ASV %in% rownames(asv),] # create matching taxonomy table after removal

# Remove contamination and/or unwanted groups
unique(tax$Division)
unique(tax$Supergroup)

# Remove unwanted divisions, NAs in higher taxonomic ranks, parasites and multicellular organisms
tax <- filter(tax, tax$Division!="Metazoa")
tax <- filter(tax, tax$Division!="Fungi")
tax <- filter(tax, tax$Division!="Cryptophyta:nucl")
tax <- filter(tax, tax$Division!="Chlorophyta:plas")
tax <- filter(tax, tax$Division!="Ochrophyta:plas")
tax <- filter(tax, tax$Supergroup!="Archaeplastida:plas")
tax <- filter(tax, tax$Class!="Syndiniales")
tax <- filter(tax, tax$Class!="Embryophyceae")
tax <- filter(tax, tax$Class!="Klebsormidiophyceae")
tax <- tax[!is.na(tax$Division),]
tax[is.na(tax)] <- "Other"

unique(tax$Division)
unique(tax$Class)

# Rename taxa
tax <- tax %>%
  mutate(Division = recode(Division, "Stramenopiles_X" = 'Stramenopiles')) %>%
  mutate(Division = recode(Division, "Opisthokonta_X" = 'Opisthokonta')) %>%
  mutate(Division = recode(Division, "Pseudofungi" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "Filosa-Thecofilosea" = 'Thecofilosea')) %>%
  mutate(Class = recode(Class, "Filosa-Imbricatea" = 'Imbricatea')) %>%
  mutate(Class = recode(Class, "Stramenopiles_XX" = 'Stramenopiles')) %>%
  mutate(Class = recode(Class, "Picozoa_X" = 'Picozoa')) %>%
  mutate(Class = recode(Class, "Telonemia_X" = 'Telonemia'))

unique(tax$Division)
unique(tax$Class)

tax_all <- tax

# Create taxonomy file for phyto- & mixoplankton (primary producers)
tax_chlo <- filter(tax, tax$Class=="Chlorarachniophyceae")
tax <- filter(tax, tax$Supergroup!="Cercozoa")
tax <- filter(tax, tax$Supergroup!="Apusozoa")
tax <- filter(tax, tax$Supergroup!="Amoebozoa")
tax <- filter(tax, tax$Supergroup!="Rhizaria")
tax <- filter(tax, tax$Division!="Choanoflagellida")
tax <- filter(tax, tax$Division!="Pseudofungi")
tax <- filter(tax, tax$Division!="Opalozoa")
tax <- filter(tax, tax$Division!="Picozoa")
tax <- filter(tax, tax$Division!="Telonemia")
tax <- filter(tax, tax$Division!="Katablepharidophyta")
tax <- filter(tax, tax$Division!="Sagenista")
tax <- filter(tax, tax$Division!="Apicomplexa")
tax <- filter(tax, tax$Division!="Mesomycetozoa")
tax <- filter(tax, tax$Division!="Centroheliozoa")
tax <- filter(tax, tax$Division!="Stramenopiles_X")
tax <- filter(tax, tax$Division!="Streptophyta")
tax <- filter(tax, tax$Division!="Rhodelphidia")
tax <- filter(tax, tax$Division!="Perkinsea")
tax <- filter(tax, tax$Division!="Opisthokonta_X")
tax <- filter(tax, tax$Class!="Noctilucophyceae")

tax_pp <- bind_rows(tax, tax_chlo)

tax <- tax_all

# Create taxonomy file for microzooplankton
tax_rhod <- filter(tax, tax$Class=="Rhodelphea")
tax_noct <- filter(tax, tax$Class=="Noctilucophyceae")
tax <- filter(tax, tax$Supergroup!="Archaeplastida")
tax <- filter(tax, tax$Division!="Haptophyta")
tax <- filter(tax, tax$Division!="Ochrophyta")
tax <- filter(tax, tax$Division!="Ciliophora")
tax <- filter(tax, tax$Division!="Chlorophyta")
tax <- filter(tax, tax$Division!="Rhodophyta")
tax <- filter(tax, tax$Division!="Cryptophyta")
tax <- filter(tax, tax$Division!="Dinoflagellata")
tax <- filter(tax, tax$Division!="Prasinodermophyta")
tax <- filter(tax, tax$Class!="Chlorarachniophyceae")

tax_mz <- bind_rows(tax, tax_rhod, tax_noct)

# Find missing groups and adapt tax selection above
require(sqldf)
tax_new <- bind_rows(tax_mz, tax_pp)
res <- sqldf('SELECT * FROM tax_all EXCEPT SELECT * FROM tax_new') 
print("rows from tax_new which are not in tax_all") 
print(res)

# Find double groups
newd <-  tax_new %>% group_by(ASV) %>% filter(n()>1) #
newd

# Subset asv tabs based on newly selected taxonomy
asv <- asv[rownames(asv) %in% tax_all$ASV,]
asv_pp <- asv[rownames(asv) %in% tax_pp$ASV,]
asv_mz <- asv[rownames(asv) %in% tax_mz$ASV,]

## Create a phyloseq object of all raw data for rarefaction
#tax_raw <- tax_all
#rownames(tax_raw) <- tax_raw$ASV
#tax_raw <- tax_raw[,2:9] # delete ASV column
#tax_raw <- as.matrix(tax_raw)
rownames(sam) <- sam$sample_ID
#OTU = otu_table(asv, taxa_are_rows = TRUE)
#TAX = phyloseq::tax_table(tax_raw)
sam$uni_ID <- as.factor(sam$uni_ID)
samples = sample_data(sam)
#ps_raw <- phyloseq(OTU, TAX, samples)
#rarecurve(t(otu_table(ps_raw)), step=50, cex=0.5) 
# looks fine

## Scaling with ranked subsampling (srs)
# All samples will be scaled to sample with lowest sequencing depth
depth.min <- min(colSums(asv))
asv.srs <- SRS(asv, depth.min) # running the SRS function
rownames(asv.srs) <- rownames(asv)
rm(depth.min)

# Check and adjust that rownames sam tab and colnames of asv tab match
all(rownames(sam) == colnames(asv))

# Prepare tax file for phyloseq object
rownames(tax_all) <- tax_all$ASV
tax <- tax_all[,2:9] # delete ASV column
tax <- as.matrix(tax)

# Name elements for phyloseq object with scaled data
OTU = otu_table(asv.srs, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(tax)
SAM = sample_data(sam)

# Create phyloseq object
ps <- phyloseq(OTU, TAX, SAM)

## Merge replicates
mps <- merge_samples(ps,"uni_ID")
replicates <- sample_data(ps)$uni_ID%>%sort()%>%table()%>%as.data.frame
otu_table(mps) <- otu_table(mps)/replicates$Freq
newmetdat <- sample_data(ps)%>%as_tibble()%>%distinct(uni_ID, .keep_all = TRUE)%>%arrange(uni_ID)%>%as.data.frame()
rownames(newmetdat) <- rownames(otu_table(mps))
sample_data(mps) <- newmetdat
ps_merged <- mps

# Subset for timepoints
ps_t0 <- subset_samples(ps, ju_day=="246")
ps_preHW <- subset_samples(ps, ju_day=="253")
ps_staHW <- subset_samples(ps, ju_day=="256")
ps_endHW <- subset_samples(ps, ju_day=="260")
ps_posHW <- subset_samples(ps, ju_day=="263")
ps_tfin <- subset_samples(ps, ju_day=="272")

## Normalization using CLR transformation
# Remove zeroes
czm <- cmultRepl(t(asv),  label=0, method="CZM") 

# Clr transformation
rho <- propr(czm, metric = 'rho', ivar = 'clr', symmetrize = TRUE, p=0)

# Clr data
clr <- rho@logratio

# Transpose the normalised and filtered asv tab & create dataframe for r
asv.clr <- t(clr)
asv.clr <- as.data.frame(asv.clr)

# Name elements for phyloseq object with scaled data
OTU = otu_table(asv.clr, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(tax)
SAM = sample_data(sam)

# Create clr-transformed phyloseq object
ps_clr <- phyloseq(OTU, TAX, SAM)

# Create clr transformed timepoint subsets
ps_t0_clr <- subset_samples(ps_clr, ju_day=="246")
ps_preHW_clr <- subset_samples(ps_clr, ju_day=="253")
ps_staHW_clr <- subset_samples(ps_clr, ju_day=="256")
ps_endHW_clr <- subset_samples(ps_clr, ju_day=="260")
ps_posHW_clr <- subset_samples(ps_clr, ju_day=="263")
ps_tfin_clr <- subset_samples(ps_clr, ju_day=="272")

#### BARGRAPHS GENUS LEVEL ####
## Create a table of all classes and abundances
class <- phyloseq::tax_glom(ps, "Class")
df <- plot_bar(class, fill = "Class")
df2 <- df$data
class <- df2 %>% select(Sample, Class, Abundance, Treatment, ju_day, replicate)

# Rename classes for prettier plotting
#class$Class[class$Abundance < 300] <- "Other"
class <- class %>%
  mutate(Class = recode(Class, "MAST-1" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-2" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-3" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-8" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-12" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MOCH-3" = "MOCH")) %>%
  mutate(Class = recode(Class, "Prymnesiophyceae" = "Haptophyta")) %>%
  mutate(Class = recode(Class, "Haptophyta_Clade_HAP5" = "Haptophyta"))

# Get percentages of classes (and add manually pre-processing)
class_sum <- aggregate(Abundance~Treatment+ju_day+Class, data=class, FUN=sum)
maxvalues <- aggregate(Abundance~Treatment+ju_day, data=class, FUN=sum)
maxvalues

## Classes of PR2 are spanning several taxonomic levels, rename to group
colnames(class)[2] <- "Group"

## Create color palette
class_pal <- qualpal(20, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

## Plotting with and without legend for manual post-processing
class_plot <- ggplot(class, aes(fill = Group, x = ju_day, y = Abundance)) +
  facet_wrap(~ Treatment, ncol = 2) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_fill_manual(values = class_pal$hex)

class_plot
ggsave("Output/PhytoGroups.png", class_plot, height = 8, width = 14, dpi = 320)

### Bargraph of dominant genera over time
genus <- phyloseq::tax_glom(ps, "Genus")
df <- plot_bar(genus, fill = "Genus")
df2 <- df$data
genus <- df2 %>% select(Sample, Genus, Abundance, Treatment, ju_day, replicate)

# Rename species for prettier plotting
genus$Genus[genus$Abundance < 100] <- "Other"
#species <- species %>%
#  mutate(Species = recode(Species, "Chaetoceros_diadema_1" = 'Chaetoceros_diadema'))
#species$Species <- gsub('_', ' ', species$Species)

## Create color palette
gen_pal <- qualpal(40, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

## Plotting with and without legend
genus_plot <- ggplot(genus, aes(fill = Genus, x = ju_day, y = Abundance)) +
  facet_wrap(~ Treatment, ncol = 2) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_fill_manual(values = gen_pal$hex)

genus_plot
ggsave("Output/PhytoGenera.png", genus_plot, height = 8, width = 14, dpi = 320)

### Bargraph of dominant species over time
species <- phyloseq::tax_glom(ps, "Species")
df <- plot_bar(species, fill = "Species")
df2 <- df$data
species <- df2 %>% select(Sample, Species, Abundance, Treatment, ju_day, replicate)

# Rename species for prettier plotting
species$Species[species$Abundance < 200] <- "Other"
#species <- species %>%
#  mutate(Species = recode(Species, "Chaetoceros_diadema_1" = 'Chaetoceros_diadema'))
#species$Species <- gsub('_', ' ', species$Species)

## Create color palette
spe_pal <- qualpal(36, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

## Plotting with and without legend
species_plot <- ggplot(species, aes(fill = Species, x = ju_day, y = Abundance)) +
  facet_wrap(~ Treatment, ncol = 2) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_fill_manual(values = spe_pal$hex)

species_plot
ggsave("Output/PhytoSpecies.png", species_plot, height = 8, width = 14, dpi = 320)

# looks different than microscopy results... check if species are in pr2
library("pr2database")
data("pr2")
"Chaetoceros_protuberans" %in% pr2$species
#TRUE

#### BARGRAPHS TIMEPOINTS SPECIES LEVEL REPLICATES ####
## Create a table of all classes and abundances
class <- phyloseq::tax_glom(ps_unmerged, "Class")
df <- plot_bar(class, fill = "Class")
df2 <- df$data
class <- df2 %>% select(Sample, Class, Abundance, Treatment, ju_day, replicate)

# Rename classes for prettier plotting
#class$Class[class$Abundance < 300] <- "Other"
class <- class %>%
  mutate(Class = recode(Class, "MAST-1" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-2" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-3" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-8" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-12" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MOCH-3" = "MOCH")) %>%
  mutate(Class = recode(Class, "Prymnesiophyceae" = "Haptophyta"))

# Get percentages of classes (and add manually pre-processing)
class_sum <- aggregate(Abundance~Treatment+ju_day+Class, data=class, FUN=sum)
maxvalues <- aggregate(Abundance~Treatment+ju_day, data=class, FUN=sum)
maxvalues

## Classes of PR2 are spanning several taxonomic levels, rename to group
colnames(class)[2] <- "Group"

## Create color palette
class_pal <- qualpal(21, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

## Plotting with and without legend for manual post-processing
class_plot <- ggplot(class, aes(fill = Group, x = ju_day, y = Abundance)) +
  facet_wrap(replicate ~ Treatment, ncol = 4) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_fill_manual(values = class_pal$hex)

class_plot
ggsave("Output/DNAPhytoGroups.png", class_plot, height = 10, width = 18, dpi = 320)

### Bargraph of dominant genera over time
genus <- phyloseq::tax_glom(ps_unmerged, "Genus")
df <- plot_bar(genus, fill = "Genus")
df2 <- df$data
genus <- df2 %>% select(Sample, Genus, Abundance, Treatment, ju_day, replicate)

# Rename species for prettier plotting
genus$Genus[genus$Abundance < 100] <- "Other"
#species <- species %>%
#  mutate(Species = recode(Species, "Chaetoceros_diadema_1" = 'Chaetoceros_diadema'))
#species$Species <- gsub('_', ' ', species$Species)

## Create color palette
gen_pal <- qualpal(38, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

## Plotting with and without legend
genus_plot <- ggplot(genus, aes(fill = Genus, x = ju_day, y = Abundance)) +
  facet_wrap(replicate ~ Treatment, ncol = 4) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_fill_manual(values = gen_pal$hex)

genus_plot
ggsave("Output/PhytoGenera.png", genus_plot, height = 10, width = 20, dpi = 320)

### Bargraph of dominant species over time
species <- phyloseq::tax_glom(ps_unmerged, "Species")
df <- plot_bar(species, fill = "Species")
df2 <- df$data
species <- df2 %>% select(Sample, Species, Abundance, Treatment, ju_day, replicate)

# Rename species for prettier plotting
species$Species[species$Abundance < 200] <- "Other"
#species <- species %>%
#  mutate(Species = recode(Species, "Chaetoceros_diadema_1" = 'Chaetoceros_diadema'))
#species$Species <- gsub('_', ' ', species$Species)

## Create color palette
spe_pal <- qualpal(36, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

## Plotting with and without legend
species_plot <- ggplot(species, aes(fill = Species, x = ju_day, y = Abundance)) +
  facet_wrap(replicate ~ Treatment, ncol = 4) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_fill_manual(values = spe_pal$hex)

species_plot
ggsave("Output/PhytoSpecies.png", species_plot, height = 10, width = 20, dpi = 320)

#### CALCULATE & PLOT DIVERSITY ####
ps.rich <- microbial::richness(ps,  method = c("Observed", "Evenness", "Shannon"))

# Add diversity measures to sample tab
sam$Richness <- ps.rich$Observed
sam$Evenness <- ps.rich$Evenness
sam$Shannon <- ps.rich$Shannon
sam$sample_ID <- NULL
sam$cosm <- NULL
sam$Treatment <- paste(sam$scenario, "-", sam$heatwave)
sam <- sam %>%
  mutate(Treatment = recode(Treatment, 'ambient - no' = 'Ambient')) %>%
  mutate(Treatment = recode(Treatment, 'ambient - yes' = 'Ambient+HW')) %>%
  mutate(Treatment = recode(Treatment, 'rcp - no' = 'RCP')) %>%
  mutate(Treatment = recode(Treatment, 'rcp - yes' = 'RCP+HW'))
sam$Treatment <- as.factor(sam$Treatment)

# Plot diversity
div <- sam

# Shannon
# Create summary of data for shannon indec
div_shan <- data_summary(div, varname="Shannon", 
                          groupnames=c("ju_day", "Treatment"))

shan_time <- ggplot(div_shan, aes(x=ju_day, y=Shannon, color=Treatment)) + 
  geom_point(position=position_dodge(0.05), size = 3)+
  geom_line(size = 1)+
  geom_errorbar(aes(ymin=Shannon-se, ymax=Shannon+se), size=1.1, width=.1,
                position=position_dodge(0.05)) +
  theme_classic() + theme(axis.line = element_blank(), 
                          text = element_text(size = 16, color="black"), 
                          axis.text.x = element_text(face="plain", color="black", size=16),
                          axis.text.y = element_text(face="plain", color="black", size=16),
                          panel.background = element_rect(colour = "black", size=1),
                          legend.position="none") +
  labs(x="Julian day", y=bquote("Shannon index")) + 
  scale_color_manual(values=treat_pal) +
  ggtitle("Shannon over time")

shan_time
ggsave("Output/AllShannonOverTime.png", shan_time, height = 9, width = 15, dpi = 320)

# Richness
# Create summary of data for species richness
div_rich <- data_summary(div, varname="Richness", 
                          groupnames=c("ju_day", "Treatment"))

rich_time <- ggplot(div_rich, aes(x=ju_day, y=Richness, color=Treatment)) + 
  geom_point(position=position_dodge(0.05), size = 3)+
  geom_line(size = 1)+
  geom_errorbar(aes(ymin=Richness-se, ymax=Richness+se), size=1.1, width=.1,
                position=position_dodge(0.05)) +
  theme_classic() + theme(axis.line = element_blank(), 
                          text = element_text(size = 16, color="black"), 
                          axis.text.x = element_text(face="plain", color="black", size=16),
                          axis.text.y = element_text(face="plain", color="black", size=16),
                          panel.background = element_rect(colour = "black", size=1),
                          legend.position="none") +
  labs(x="Day", y=bquote("Species Richness")) + 
  scale_color_manual(values=treat_pal) +
  ggtitle("Richness over time")

rich_time
ggsave("Output/PhytoRichnessOverTime.png", rich_time, height = 9, width = 15, dpi = 320)

# Evenness
# Create summary of data for species evenness
div_even <- data_summary(div, varname="Evenness", 
                          groupnames=c("ju_day", "Treatment"))

even_time <- ggplot(div_even, aes(x=ju_day, y=Evenness, color=Treatment)) + 
  geom_point(position=position_dodge(0.05), size = 3)+
  geom_line(size = 1)+
  geom_errorbar(aes(ymin=Evenness-se, ymax=Evenness+se), size=1.1, width=.1,
                position=position_dodge(0.05)) +
  theme_classic() + theme(axis.line = element_blank(), 
                          text = element_text(size = 16, color="black"), 
                          axis.text.x = element_text(face="plain", color="black", size=16),
                          axis.text.y = element_text(face="plain", color="black", size=16),
                          panel.background = element_rect(colour = "black", size=1),
                          legend.position="none") +
  labs(x="Day", y=bquote("Species Evenness")) + 
  scale_color_manual(values=treat_pal) +
  ggtitle("Evenness over time")

even_time
ggsave("Output/PPEvennessOverTime.png", even_time, height = 9, width = 15, dpi = 320)


div_tfin <- subset(div, ju_day == "272")

anova <- aov(Richness ~ Treatment, data = div_tfin)
summary(anova)
TukeyHSD(anova)

#### ORDINATION OF TIMEPOINTS ####
## t0 ordination
ordi_t0 <- ps_t0_clr %>%
  dist_calc("euclidean") %>%
  ord_calc() %>%
  ord_plot(shape = "Treatment", size = 4, color = "Treatment") +
  scale_color_manual(values=treat_pal) +
  ggforce::geom_mark_ellipse(aes(color = Treatment))

ordi_t0

ggsave("Output/Ordination_t0.png", ordi_t0, dpi = 300, width = 8, height = 4)

## pre_HW ordination
ordi_preHW <- ps_preHW_clr %>%
  dist_calc("euclidean") %>%
  ord_calc() %>%
  ord_plot(shape = "Treatment", size = 4, color = "Treatment") +
  scale_color_manual(values=treat_pal) +
  ggforce::geom_mark_ellipse(aes(color = Treatment))

ordi_preHW

ggsave("Output/Ordination_preHW.png", ordi_preHW, dpi = 300, width = 8, height = 4)

## sta_HW ordination
ordi_staHW <- ps_staHW_clr %>%
  dist_calc("euclidean") %>%
  ord_calc() %>%
  ord_plot(shape = "Treatment", size = 4, color = "Treatment") +
  scale_color_manual(values=treat_pal) +
  ggforce::geom_mark_ellipse(aes(color = Treatment))

ordi_staHW

ggsave("Output/Ordination_staHW.png", ordi_staHW, dpi = 300, width = 8, height = 4)

## end_HW ordination
ordi_endHW <- ps_endHW_clr %>%
  dist_calc("euclidean") %>%
  ord_calc() %>%
  ord_plot(shape = "Treatment", size = 4, color = "Treatment") +
  scale_color_manual(values=treat_pal) +
  ggforce::geom_mark_ellipse(aes(color = Treatment))

ordi_endHW

ggsave("Output/Ordination_endHW.png", ordi_endHW, dpi = 300, width = 8, height = 4)

## pos_HW ordination
ordi_posHW <- ps_posHW_clr %>%
  dist_calc("euclidean") %>%
  ord_calc() %>%
  ord_plot(shape = "Treatment", size = 4, color = "Treatment") +
  scale_color_manual(values=treat_pal) +
  ggforce::geom_mark_ellipse(aes(color = Treatment))

ordi_posHW

ggsave("Output/Ordination_posHW.png", ordi_posHW, dpi = 300, width = 8, height = 4)

## tfin ordination
ordi_tfin <- ps_tfin_clr %>%
  dist_calc("euclidean") %>%
  ord_calc() %>%
  ord_plot(shape = "Treatment", size = 4, color = "Treatment") +
  scale_color_manual(values=treat_pal) +
  ggforce::geom_mark_ellipse(aes(color = Treatment))

ordi_tfin

ggsave("Output/Ordination_tfin.png", ordi_tfin, dpi = 300, width = 8, height = 4)

#### STATISTICS ####
## Create subsets for phases and groups of interest
rcp <- subset(div, heatwave == "no")
hwn <- subset(div, scenario == "ambient")
hwt <- subset(div, scenario == "rcp")
hwn_warm <- subset(hwn, ju_day == "253" | ju_day == "256")
hwn_cool <- subset(hwn, ju_day == "260" | ju_day == "263")
hwn_rec <- subset(hwn, ju_day > "263")
hwt_warm <- subset(hwt, ju_day == "253" | ju_day == "256")
hwt_cool <- subset(hwt, ju_day == "260" | ju_day == "263")
hwt_rec <- subset(hwt, ju_day > "263")

## Investigate effect of RCP scenario alone on diversity
require(nlme)
newModel<-lme(Shannon~ju_day, random = ~1|replicate/Treatment, data = rcp, method ="ML")
summary(newModel)

#### Check for specific species & groups ####
df_tax <- as.data.frame(tax)

"Emiliania" %in% df_tax$Genus #FALSE
"Gephyrocapsa" %in% df_tax$Genus #TRUE

# Is Emiliania huxleyi in PR2?
#devtools::install_github("pr2database/pr2database")
require(pr2database)
data("pr2")
"Emiliania_huxleyi" %in% pr2$species
#TRUE

# Create a plot of all haptophytes
cocco <- subset_taxa(ps, Division == "Haptophyta")
cocco_species <- tax_glom(cocco, "Species")
taxa_names(cocco_species) <- tax_table(cocco_species)[, "Species"]
otu_table(cocco_species)[1:5, 1:5]
print(colnames(otu_table(cocco_species)))
coc <- psmelt(cocco_species)
coc$Species <- gsub('_', ' ', coc$Species)

coc <- coc %>%
  mutate(Treatment = recode(Treatment, "ambient-yes" = 'Amb+HW')) %>%
  mutate(Treatment = recode(Treatment, "ambient-no" = 'Amb')) %>%
  mutate(Treatment = recode(Treatment, "rcp-yes" = 'RCP+HW')) %>%
  mutate(Treatment = recode(Treatment, "rcp-no" = 'RCP'))

coc$Species[coc$Abundance < 300] <- "Other"

haptos <- ggplot(coc, aes(x = Treatment, y = Abundance, group = Treatment, fill=Treatment)) +
  geom_boxplot(outlier.shape  = NA, width =0.4) +
  geom_jitter(height = 0, width = .2, size =0.8) +
  labs(x = "", y = "Abundance\n") +
  facet_grid(~ Species, scales = "free")+
  RDA.theme +
  theme(strip.text = element_text(size = 15, face = "italic"))+
  scale_fill_manual(values=treat_pal)+
  ggtitle("Haptophytes")
  
haptos
ggsave("Output/Haptophyta.png",haptos, height = 7, width = 15, dpi = 320)

geph <- subset(coc, Species == "Gephyrocapsa oceanica")

geph_time <- ggplot(geph, aes(x = Treatment, y = Abundance, group = Treatment, fill=Treatment, color=Treatment)) +
  geom_boxplot(outlier.shape  = NA, width =0.4) +
  geom_point(size = 4) +
  labs(x = "", y = "Abundance\n") +
  facet_grid(~ ju_day, scales = "free")+
  RDA.theme +
  theme(strip.text = element_text(size = 15, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_fill_manual(values=treat_pal)+
  scale_color_manual(values=treat_pal)+
  ggtitle("Gephyrocapsa oceanica over time")

geph_time
ggsave("Output/GephyrocapsaOverTime.png", geph_time, height = 6, width = 10, dpi = 320)
