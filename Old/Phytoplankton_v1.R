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
treat_pal <- c("AMB" = "skyblue2","AMB+HW" = "darkblue", "FUT" = "goldenrod2", "FUT+HW" = "tomato2")

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
tax <- filter(tax, tax$Genus!="Minorisa")

tax <- rbind(tax, tax_chlo)

# Subset asv tabs based on newly selected taxonomy
asv <- asv[rownames(asv) %in% tax$ASV,]

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
rownames(tax) <- tax$ASV
tax <- tax[,2:9] # delete ASV column
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
  mutate(Treatment = recode(Treatment, 'ambient - no' = 'AMB')) %>%
  mutate(Treatment = recode(Treatment, 'ambient - yes' = 'AMB+HW')) %>%
  mutate(Treatment = recode(Treatment, 'rcp - no' = 'FUT')) %>%
  mutate(Treatment = recode(Treatment, 'rcp - yes' = 'FUT+HW'))
sam$Treatment <- as.factor(sam$Treatment)

### Plot diversity
div <- sam

# Create subsets for phases and groups of interest
rcp <- subset(div, heatwave == "no")
hwn <- subset(div, scenario == "ambient")
hwt <- subset(div, scenario == "rcp")
hwn_warm <- subset(hwn, ju_day == "253" | ju_day == "256")
hwn_cool <- subset(hwn, ju_day == "260" | ju_day == "263")
hwn_rec <- subset(hwn, ju_day > "263")
hwt_warm <- subset(hwt, ju_day == "253" | ju_day == "256")
hwt_cool <- subset(hwt, ju_day == "260" | ju_day == "263")
hwt_rec <- subset(hwt, ju_day > "263")

## Shannon
# Create summary of data for shannon index
shan_time <- ggplot(div, aes(x=ju_day, y=Shannon, color=Treatment)) + 
  geom_point(position=position_dodge(0.05), size = 3)+
  geom_smooth(method = "gam",size = 1)+
  plot.theme +
  labs(x="Julian date (d)", y=bquote("Shannon index")) + 
  scale_color_manual(values=treat_pal)

shan_time
ggsave("Output/PPShannonOverTime.png", shan_time, height = 9, width = 15, dpi = 320)


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

#### STATISTICS OF DIVERSITY PARAMETERS ####
 
## Investigate effect of RCP scenario alone on diversity


#### PACKAGES ####
## Check which packages you actually used for documentation and tidying purposes
require(NCmisc)
packages <- list.functions.in.file("~/AWI/RProjects/TopTrons/CommunityComposition.R", alphabetic = TRUE) # set to your filepath
summary(packages)

## Download citation report and citations
pacman::p_load(grateful)
cite_packages(out.format = "docx", out.dir = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/AP2-Helgoland/", pkgs = "Session", out.file = "GDA_packages")
