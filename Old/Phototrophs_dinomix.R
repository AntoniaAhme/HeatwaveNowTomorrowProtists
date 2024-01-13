### HNT22 mesocosms 18S rRNA metabarcoding data ###
### Autotrophic protists ###
## Antonia Ahme, 07.12.2023 ##

rm(list=ls())
options(stringsAsFactors = F)
setwd("~/AWI/RProjects/HNT22") # Set to your own working directory

set.seed(22)

#### HOUSEKEEPING ####
require(dplyr)
require(tidyr)
require(vegan)
require(ggplot2)
require(phyloseq)
require(SRS)
require(qualpalr)
require(writexl)
require(rstatix)
require(mgcv)
require(nlme)

#### FUNCTIONS & LAYOUTS ####
treat_pal <- c("AMB+HW" = "royalblue2", "FUT" = "goldenrod2", "FUT+HW" = "firebrick1")
treat_pal2 <- c("AMB" = "skyblue2", "AMB+HW" = "royalblue3", "FUT" = "goldenrod2", "FUT+HW" = "firebrick1")
treat_pch <- c(16,17)

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
                   text = element_text(size = 20, face = "bold"), legend.position = "right",
                   legend.text = element_text(face = "italic"))

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

## Functions for derivatives of GAM(M) models ##
# by GAVIN SIMPSON             
Deriv <- function(mod, n = 200, eps = 1e-7, newdata, term) {
  if(inherits(mod, "gamm"))
    mod <- mod$gam
  m.terms <- attr(terms(mod), "term.labels")
  if(missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x) seq(min(x), max(x), length = n))
    names(newD) <- m.terms
  } else {
    newD <- newdata
  }
  X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
  newD <- newD + eps
  X1 <- predict(mod, data.frame(newD), type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  ## number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  ## match the term with the the terms in the model
  if(!missing(term)) {
    want <- grep(term, t.labs)
    if(!identical(length(want), length(term)))
      stop("One or more 'term's not found in model!")
    t.labs <- t.labs[want]
  }
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  for(i in seq_len(nt)) {
    Xi <- Xp * 0
    want <- grep(t.labs[i], colnames(X1))
    Xi[, want] <- Xp[, want]
    df <- Xi %*% coef(mod)
    df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    lD[[i]] <- list(deriv = df, se.deriv = df.sd)
  }
  class(lD) <- "Deriv"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  lD ##return
}

confint.Deriv <- function(object, term, alpha = 0.05, ...) {
  l <- length(object) - 3
  term.labs <- names(object[seq_len(l)])
  if(missing(term)) {
    term <- term.labs
  } else { ## how many attempts to get this right!?!?
    ##term <- match(term, term.labs)
    ##term <- term[match(term, term.labs)]
    term <- term.labs[match(term, term.labs)]
  }
  if(any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  res <- vector(mode = "list", length = length(term))
  names(res) <- term
  residual.df <- df.residual(object$gamModel)
  tVal <- qt(1 - (alpha/2), residual.df)
  ##for(i in term.labs[term]) {
  for(i in term) {
    upr <- object[[i]]$deriv + tVal * object[[i]]$se.deriv
    lwr <- object[[i]]$deriv - tVal * object[[i]]$se.deriv
    res[[i]] <- list(upper = drop(upr), lower = drop(lwr))
  }
  res$alpha = alpha
  res
}

signifD <- function(x, d, upper, lower, eval = 0) {
  miss <- upper > eval & lower < eval
  incr <- decr <- x
  want <- d > eval
  incr[!want | miss] <- NA
  want <- d < eval
  decr[!want | miss] <- NA
  list(incr = incr, decr = decr)
}

plot.Deriv <- function(x, alpha = 0.05, polygon = TRUE,
                       sizer = FALSE, term,
                       eval = 0, lwd = 3,
                       col = "lightgrey", border = col,
                       ylab, xlab, main, ...) {
  l <- length(x) - 3
  ## get terms and check specified (if any) are in model
  term.labs <- names(x[seq_len(l)])
  if(missing(term)) {
    term <- term.labs
  } else {
    term <- term.labs[match(term, term.labs)]
  }
  if(any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  if(all(miss))
    stop("All terms in 'term' not found in model.")
  l <- sum(!miss)
  nplt <- n2mfrow(l)
  tVal <- qt(1 - (alpha/2), df.residual(x$gamModel))
  if(missing(ylab))
    ylab <- expression(italic(hat(f)*"'"*(x)))
  if(missing(xlab)) {
    xlab <- attr(terms(x$gamModel), "term.labels")
    names(xlab) <- xlab
  }
  if (missing(main)) {
    main <- term
    names(main) <- term
  }
  ## compute confidence interval
  CI <- confint(x, term = term)
  ## plots
  layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
  for(i in term) {
    upr <- CI[[i]]$upper
    lwr <- CI[[i]]$lower
    ylim <- range(upr, lwr)
    plot(x$eval[,i], x[[i]]$deriv, type = "n",
         ylim = ylim, ylab = ylab, xlab = xlab[i], main = main[i], ...)
    if(isTRUE(polygon)) {
      polygon(c(x$eval[,i], rev(x$eval[,i])),
              c(upr, rev(lwr)), col = col, border = border)
    } else {
      lines(x$eval[,i], upr, lty = "dashed")
      lines(x$eval[,i], lwr, lty = "dashed")
    }
    abline(h = 0, ...)
    if(isTRUE(sizer)) {
      lines(x$eval[,i], x[[i]]$deriv, lwd = 1)
      S <- signifD(x[[i]]$deriv, x[[i]]$deriv, upr, lwr,
                   eval = eval)
      lines(x$eval[,i], S$incr, lwd = lwd, col = "blue")
      lines(x$eval[,i], S$decr, lwd = lwd, col = "red")
    } else {
      lines(x$eval[,i], x[[i]]$deriv, lwd = 2)
    }
  }
  layout(1)
  invisible(x)
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
plot(depths)
min(depths)
# looks okay, no incredeibly low numbers, so keep as it is

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
# based on https://doi.org/10.1111/jeu.12691 and other literature data
tax_chlo <- filter(tax, tax$Class=="Chlorarachniophyceae")
tax <- filter(tax, tax$Supergroup!="Cercozoa")
tax <- filter(tax, tax$Supergroup!="Apusozoa")
tax <- filter(tax, tax$Supergroup!="Amoebozoa")
tax <- filter(tax, tax$Supergroup!="Rhizaria")
tax <- filter(tax, tax$Division!="Choanoflagellida")
tax <- filter(tax, tax$Division!="Ciliophora")
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
tax <- filter(tax, tax$Family!="Amphisoleniaceae")
tax <- filter(tax, tax$Genus!="Gyrodinium")
tax <- filter(tax, tax$Genus!="Islandinium")
tax <- filter(tax, tax$Genus!="Amylax")
tax <- filter(tax, tax$Genus!="Archaeperidinium")
tax <- filter(tax, tax$Genus!="Dinophysis")
tax <- filter(tax, tax$Genus!="Diplopsalis")
tax <- filter(tax, tax$Genus!="Luciella")
tax <- filter(tax, tax$Genus!="Pentapharsodinium")
tax <- filter(tax, tax$Genus!="Phalacroma")
tax <- filter(tax, tax$Genus!="Polykrikos")
tax <- filter(tax, tax$Genus!="Protoperidinium")
tax <- filter(tax, tax$Genus!="Qia")
tax <- filter(tax, tax$Genus!="Protodinium")
tax <- filter(tax, tax$Genus!="Stockeria")

tax <- rbind(tax, tax_chlo)

# Subset asv tabs based on newly selected taxonomy
asv <- asv[rownames(asv) %in% tax$ASV,]

## Create a phyloseq object of all raw data for rarefaction
tax_raw <- tax_all
rownames(tax_raw) <- tax_raw$ASV
tax_raw <- tax_raw[,2:9] # delete ASV column
tax_raw <- as.matrix(tax_raw)
rownames(sam) <- sam$sample_ID
OTU = otu_table(asv, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(tax_raw)
sam$uni_ID <- as.factor(sam$uni_ID)
samples = sample_data(sam)
ps_raw <- phyloseq(OTU, TAX, samples)
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

#### CALCULATE DIVERSITY & PLOT OVERVIEW GRAPH SHANNON ####
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
sam$time <- as.numeric(sam$time)

### Prepare diversity analyses
div <- sam

## Save diversity data as excel
write_xlsx(div, "Data/Diversity_PP.xlsx")

# Shannon
# Create summary of data for shannon indec
div_shan <- data_summary(div, varname="Shannon", 
                         groupnames=c("time", "Treatment"))

shan_time <- ggplot(div_shan, aes(x=time, y=Shannon, color=Treatment)) + 
  geom_rect(data=NULL,aes(xmin=9,xmax=17,ymin=-Inf,ymax=Inf),
            fill="lightgoldenrod") +
  geom_point(position=position_dodge(0.05), size = 3)+
  geom_line(size = 1)+
  geom_errorbar(aes(ymin=Shannon-se, ymax=Shannon+se), size=1.1, width=1,
                position=position_dodge(0.05)) +
  theme_classic() + theme(axis.line = element_blank(), 
                          text = element_text(size = 16, color="black"), 
                          axis.text.x = element_text(face="plain", color="black", size=16),
                          axis.text.y = element_text(face="plain", color="black", size=16),
                          panel.background = element_rect(colour = "black", size=1),
                          legend.position=c(.1,.2)) +
  labs(x="Incubation time (d)", y=bquote("Shannon index")) + 
  scale_color_manual(values=treat_pal2)

shan_time
ggsave("Output/AllShannonOverTime_PP.png", shan_time, height = 5, width = 10, dpi = 320)

#### CALCULATE & PLOT & ANALYSE SHANNON LRRs ####
## Create subsets for phases and groups of interest
amb <- subset(div, heatwave == "no" & scenario == "ambient")
fut <- subset(div, heatwave == "no" & scenario == "rcp")
hw_all <- subset(div, time > 7 & heatwave == "yes")
hwa <- subset(hw_all, scenario == "ambient")
hwf <- subset(hw_all, scenario == "rcp")

## Calculate log response ratios
## FUT
means_amb <- amb %>%
  group_by(time) %>%
  get_summary_stats(Shannon, type = "mean")

fut_m <- merge(fut, means_amb, by = "time")
fut_m$LRR <- log(fut_m$Shannon/fut_m$mean)

## HW-A
hwa_m <- merge(hwa, means_amb, by = "time")
hwa_m$LRR <- log(hwa_m$Shannon/hwa_m$mean)

## HW-F
means_fut <- fut %>%
  group_by(time) %>%
  get_summary_stats(Shannon, type = "mean")

hwf_m <- merge(hwf, means_fut, by = "time")
hwf_m$LRR <- log(hwf_m$Shannon/hwf_m$mean)

## Merged HW dataframe fr plotting
hw <- rbind(hwa_m, hwf_m)

### Statistics for significant in- and decreasees
# Model GAM and then use first derivatives to detect significant changes
# After Gavin Simson

### FUT
m1 <- gam(LRR ~ s(time, k=3, fx= TRUE), data = fut_m)
m1$aic # knot number to decrease aic value
summary(m1)
plot(m1, residuals = TRUE, pch = 19, cex = 0.75)

want <- seq(1, nrow(fut_m), length.out = 200)
pdat <- with(fut_m,
             data.frame(time = time[want]))
p2 <- predict(m1, newdata = pdat, type = "terms", se.fit = TRUE)
pdat <- transform(pdat, p2 = p2$fit[,1], se2 = p2$se.fit[,1])

df.res <- df.residual(m1)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat <- transform(pdat,
                  upper = p2 + (crit.t * se2),
                  lower = p2 - (crit.t * se2))

Term <- "time"
m1.d <- Deriv(m1)
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p2, d = m1.d[[Term]]$deriv,
                   +                    m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)

# Plot to see siginificant changes and manually add to own plots
ylim <- with(pdat, range(upper, lower, p2))
ylab <- expression(LRR)

plot(p2 ~ time, data = pdat, type = "n", ylab = ylab, ylim = ylim)
lines(p2 ~ time, data = pdat)
lines(upper ~ time, data = pdat, lty = "dashed")
lines(lower ~ time, data = pdat, lty = "dashed")
lines(unlist(m1.dsig$incr) ~ time, data = pdat, col = "blue", lwd = 3)
lines(unlist(m1.dsig$decr) ~ time, data = pdat, col = "red", lwd = 3)
# significant decrease from day 0 to 15

## Plot prettier for paper
# Extract gam lines
newdat <- data.frame(time=seq(0,27,0.1)) 
p2 <- predict(m1, newdata=newdat, type = "response", se.fit=TRUE)
newdat$LLRmin <- p2$fit
newdat$LLRmin_se <- p2$se.fit
newdat$LLR_upr <- p2$fit + (1.96 * p2$se.fit)
newdat$LLR_lwr <- p2$fit - (1.96 * p2$se.fit)
newdat$Treatment <- "FUT"
signi <- newdat
signi <- subset(signi, time < 15.1)
head(newdat)

# plot
fut_time <- ggplot(fut_m, aes(x=time, y=LRR, size = sig, color = Treatment)) + 
  geom_hline(yintercept=0, linetype="dashed")+
  geom_point(position=position_dodge(0.05), size = 3, alpha = 0.8) +
  geom_line(data = newdat, aes(time, LLRmin), size = 1) +
  geom_line(data = signi, aes(time, LLRmin), size = 2, color = "indianred1") +
  geom_line(data = newdat, aes(time, LLR_upr), size = .5, alpha = 0.8, linetype="dotted") +
  geom_line(data = newdat, aes(time, LLR_lwr), size = .5, alpha = 0.8, linetype="dotted") +
  plot.theme +
  labs(x="Incubation time (d)", y=bquote("LRR")) + 
  scale_color_manual(values=treat_pal) +
  coord_cartesian(ylim = c(-0.5, 0.3)) +
  scale_x_continuous(breaks = seq(0, 27, 3))
  
fut_time
ggsave("Output/PP_FUT_LRR.png", fut_time, height = 4, width = 8, dpi = 320)


### AMB-HW
m1 <- gam(LRR ~ s(time, k=5, fx= TRUE), data = hwa_m)
m1$aic
summary(m1)
plot(m1, residuals = TRUE, pch = 19, cex = 0.75)

want <- seq(1, nrow(hwa_m), length.out = 200)
pdat <- with(hwa_m,
             data.frame(time = time[want]))
p2 <- predict(m1, newdata = pdat, type = "terms", se.fit = TRUE)
pdat <- transform(pdat, p2 = p2$fit[,1], se2 = p2$se.fit[,1])
df.res <- df.residual(m1)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat <- transform(pdat,
                  upper = p2 + (crit.t * se2),
                  lower = p2 - (crit.t * se2))

Term <- "time"
m1.d <- Deriv(m1)
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p2, d = m1.d[[Term]]$deriv,
                   +                    m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)

# Plot to see sginificant changes and manually add to own plots
ylim <- with(pdat, range(upper, lower, p2))
ylab <- expression(LRR)

plot(p2 ~ time, data = pdat, type = "n", ylab = ylab, ylim = ylim)
lines(p2 ~ time, data = pdat)
lines(upper ~ time, data = pdat, lty = "dashed")
lines(lower ~ time, data = pdat, lty = "dashed")
lines(unlist(m1.dsig$incr) ~ time, data = pdat, col = "blue", lwd = 3)
lines(unlist(m1.dsig$decr) ~ time, data = pdat, col = "red", lwd = 3)

# Extract gam lines
newdat <- data.frame(time=seq(8,27,0.1)) 
p2 <- predict(m1, newdata=newdat, type = "response", se.fit=TRUE)
newdat$LLRmin <- p2$fit
newdat$LLRmin_se <- p2$se.fit
newdat$LLR_upr <- p2$fit + (1.96 * p2$se.fit)
newdat$LLR_lwr <- p2$fit - (1.96 * p2$se.fit)
newdat$Treatment <- "AMB+HW"
signi <- newdat
signi <- subset(signi, time < 15.1)
signi <- subset(signi, time > 10.9)
newdat_amb <- newdat
signi_amb1 <- signi
signi <- newdat
signi <- subset(signi, time < 22.1)
signi <- subset(signi, time > 17.9)
signi_amb2 <- signi

### FUT-HW
m1 <- gam(LRR ~ s(time, k=5, fx= TRUE), data = hwf_m)
m1$aic
summary(m1)
plot(m1, residuals = TRUE, pch = 19, cex = 0.75)

want <- seq(1, nrow(hwf_m), length.out = 200)
pdat <- with(hwf_m,
             data.frame(time = time[want]))

p2 <- predict(m1, newdata = pdat, type = "terms", se.fit = TRUE)
pdat <- transform(pdat, p2 = p2$fit[,1], se2 = p2$se.fit[,1])
df.res <- df.residual(m1)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat <- transform(pdat,
                  upper = p2 + (crit.t * se2),
                  lower = p2 - (crit.t * se2))

Term <- "time"
m1.d <- Deriv(m1)
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p2, d = m1.d[[Term]]$deriv,
                   +                    m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)

# Plot to see sginificant changes and manually add to own plots
ylim <- with(pdat, range(upper, lower, p2))
ylab <- expression(LRR)
 
plot(p2 ~ time, data = pdat, type = "n", ylab = ylab, ylim = ylim)
lines(p2 ~ time, data = pdat)
lines(upper ~ time, data = pdat, lty = "dashed")
lines(lower ~ time, data = pdat, lty = "dashed")
lines(unlist(m1.dsig$incr) ~ time, data = pdat, col = "blue", lwd = 3)
lines(unlist(m1.dsig$decr) ~ time, data = pdat, col = "red", lwd = 3)

# Extract gam lines
newdat <- data.frame(time=seq(8,27,0.1)) 
p2 <- predict(m1, newdata=newdat, type = "response", se.fit=TRUE)
newdat$LLRmin <- p2$fit
newdat$LLRmin_se <- p2$se.fit
newdat$LLR_upr <- p2$fit + (1.96 * p2$se.fit)
newdat$LLR_lwr <- p2$fit - (1.96 * p2$se.fit)
newdat$Treatment <- "FUT+HW"
signi <- newdat
signi <- subset(signi, time < 18.1)
signi <- subset(signi, time > 12.9)
newdat_fut_hw <- newdat
signi_fut_hw <- signi

## Plot both HW together
hw_time <- ggplot(hw, aes(x=time, y=LRR, color = Treatment, shape = Treatment)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_rect(data=NULL,aes(xmin=10,xmax=16,ymin=-Inf,ymax=Inf),
            fill="gold", alpha = .5)+
  geom_point(position=position_dodge(0.05), size = 3, alpha = 0.8) +
  geom_line(data = newdat_amb, aes(time, LLRmin), size = 1) +
  geom_line(data = signi_amb1, aes(time, LLRmin), size = 3, color = "royalblue1") +
  geom_line(data = signi_amb2, aes(time, LLRmin), size = 3, color = "royalblue1") +
  geom_line(data = newdat_amb, aes(time, LLR_upr), size = .5, alpha = 0.8, linetype="dotted") +
  geom_line(data = newdat_amb, aes(time, LLR_lwr), size = .5, alpha = 0.8, linetype="dotted") +
  geom_line(data = newdat_fut_hw, aes(time, LLRmin), size = 1) +
  geom_line(data = signi_fut_hw, aes(time, LLRmin), size = 3, color = "firebrick1") +
  geom_line(data = newdat_fut_hw, aes(time, LLR_upr), size = .5, alpha = 0.8, linetype="dotted") +
  geom_line(data = newdat_fut_hw, aes(time, LLR_lwr), size = .5, alpha = 0.8, linetype="dotted") +
  plot.theme +
  labs(x="Incubation time (d)", y=bquote("LRR")) + 
  scale_color_manual(values=treat_pal) +
  scale_shape_manual(values=treat_pch) +
  coord_cartesian(ylim = c(-0.5, 0.3)) +
  scale_x_continuous(breaks = seq(0, 27, 3))

hw_time
ggsave("Output/PP_HW_LRR.png", hw_time, height = 4, width = 8, dpi = 320)


#### CALCULATE & PLOT & ANALYSE RICHNESS LRRs ####
## Calculate log response ratios
## FUT
means_amb <- amb %>%
  group_by(time) %>%
  get_summary_stats(Richness, type = "mean")

fut_m <- merge(fut, means_amb, by = "time")
fut_m$LRR <- log(fut_m$Richness/fut_m$mean)

## HW-A
hwa_m <- merge(hwa, means_amb, by = "time")
hwa_m$LRR <- log(hwa_m$Richness/hwa_m$mean)

## HW-F
means_fut <- fut %>%
  group_by(time) %>%
  get_summary_stats(Richness, type = "mean")

hwf_m <- merge(hwf, means_fut, by = "time")
hwf_m$LRR <- log(hwf_m$Richness/hwf_m$mean)

## Merged HW dataframe fr plotting
hw <- rbind(hwa_m, hwf_m)

### Statistics for significant in- and decreasees
# Model GAM and then use first derivatives to detect significant changes
# After Gavin Simson

### FUT
m1 <- gam(LRR ~ s(time, k=3, fx= TRUE), data = fut_m)
m1$aic # knot number to decrease aic value
summary(m1)
plot(m1, residuals = TRUE, pch = 19, cex = 0.75)

want <- seq(1, nrow(fut_m), length.out = 200)
pdat <- with(fut_m,
             data.frame(time = time[want]))
p2 <- predict(m1, newdata = pdat, type = "terms", se.fit = TRUE)
pdat <- transform(pdat, p2 = p2$fit[,1], se2 = p2$se.fit[,1])

df.res <- df.residual(m1)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat <- transform(pdat,
                  upper = p2 + (crit.t * se2),
                  lower = p2 - (crit.t * se2))

Term <- "time"
m1.d <- Deriv(m1)
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p2, d = m1.d[[Term]]$deriv,
                   +                    m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)

# Plot to see siginificant changes and manually add to own plots
ylim <- with(pdat, range(upper, lower, p2))
ylab <- expression(LRR)

plot(p2 ~ time, data = pdat, type = "n", ylab = ylab, ylim = ylim)
lines(p2 ~ time, data = pdat)
lines(upper ~ time, data = pdat, lty = "dashed")
lines(lower ~ time, data = pdat, lty = "dashed")
lines(unlist(m1.dsig$incr) ~ time, data = pdat, col = "blue", lwd = 3)
lines(unlist(m1.dsig$decr) ~ time, data = pdat, col = "red", lwd = 3)
# significant decrease from day 0 to 15

## Plot prettier for paper
# Extract gam lines
newdat <- data.frame(time=seq(0,27,0.1)) 
p2 <- predict(m1, newdata=newdat, type = "response", se.fit=TRUE)
newdat$LLRmin <- p2$fit
newdat$LLRmin_se <- p2$se.fit
newdat$LLR_upr <- p2$fit + (1.96 * p2$se.fit)
newdat$LLR_lwr <- p2$fit - (1.96 * p2$se.fit)
newdat$Treatment <- "FUT"
signi <- newdat
signi <- subset(signi, time < 15.1)
head(newdat)

# plot
fut_time <- ggplot(fut_m, aes(x=time, y=LRR, size = sig, color = Treatment)) + 
  geom_hline(yintercept=0, linetype="dashed")+
  geom_point(position=position_dodge(0.05), size = 3, alpha = 0.8) +
  geom_line(data = newdat, aes(time, LLRmin), size = 1) +
  geom_line(data = signi, aes(time, LLRmin), size = 2, color = "indianred1") +
  geom_line(data = newdat, aes(time, LLR_upr), size = .5, alpha = 0.8, linetype="dotted") +
  geom_line(data = newdat, aes(time, LLR_lwr), size = .5, alpha = 0.8, linetype="dotted") +
  plot.theme +
  labs(x="Incubation time (d)", y=bquote("LRR")) + 
  scale_color_manual(values=treat_pal) +
  coord_cartesian(ylim = c(-0.5, 0.3)) +
  scale_x_continuous(breaks = seq(0, 27, 3))

fut_time
ggsave("Output/PP_FUT_LRR_rich.png", fut_time, height = 4, width = 8, dpi = 320)


### AMB-HW
m1 <- gam(LRR ~ s(time, k=5, fx= TRUE), data = hwa_m)
m1$aic
summary(m1)
plot(m1, residuals = TRUE, pch = 19, cex = 0.75)

want <- seq(1, nrow(hwa_m), length.out = 200)
pdat <- with(hwa_m,
             data.frame(time = time[want]))
p2 <- predict(m1, newdata = pdat, type = "terms", se.fit = TRUE)
pdat <- transform(pdat, p2 = p2$fit[,1], se2 = p2$se.fit[,1])
df.res <- df.residual(m1)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat <- transform(pdat,
                  upper = p2 + (crit.t * se2),
                  lower = p2 - (crit.t * se2))

Term <- "time"
m1.d <- Deriv(m1)
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p2, d = m1.d[[Term]]$deriv,
                   +                    m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)

# Plot to see sginificant changes and manually add to own plots
ylim <- with(pdat, range(upper, lower, p2))
ylab <- expression(LRR)

plot(p2 ~ time, data = pdat, type = "n", ylab = ylab, ylim = ylim)
lines(p2 ~ time, data = pdat)
lines(upper ~ time, data = pdat, lty = "dashed")
lines(lower ~ time, data = pdat, lty = "dashed")
lines(unlist(m1.dsig$incr) ~ time, data = pdat, col = "blue", lwd = 3)
lines(unlist(m1.dsig$decr) ~ time, data = pdat, col = "red", lwd = 3)

# Extract gam lines
newdat <- data.frame(time=seq(8,27,0.1)) 
p2 <- predict(m1, newdata=newdat, type = "response", se.fit=TRUE)
newdat$LLRmin <- p2$fit
newdat$LLRmin_se <- p2$se.fit
newdat$LLR_upr <- p2$fit + (1.96 * p2$se.fit)
newdat$LLR_lwr <- p2$fit - (1.96 * p2$se.fit)
newdat$Treatment <- "AMB+HW"


### FUT-HW
m1 <- gam(LRR ~ s(time, k=6, fx= TRUE), data = hwf_m)
m1$aic
summary(m1)
plot(m1, residuals = TRUE, pch = 19, cex = 0.75)

want <- seq(1, nrow(hwf_m), length.out = 200)
pdat <- with(hwf_m,
             data.frame(time = time[want]))

p2 <- predict(m1, newdata = pdat, type = "terms", se.fit = TRUE)
pdat <- transform(pdat, p2 = p2$fit[,1], se2 = p2$se.fit[,1])
df.res <- df.residual(m1)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat <- transform(pdat,
                  upper = p2 + (crit.t * se2),
                  lower = p2 - (crit.t * se2))

Term <- "time"
m1.d <- Deriv(m1)
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p2, d = m1.d[[Term]]$deriv,
                   +                    m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)

# Plot to see sginificant changes and manually add to own plots
ylim <- with(pdat, range(upper, lower, p2))
ylab <- expression(LRR)

plot(p2 ~ time, data = pdat, type = "n", ylab = ylab, ylim = ylim)
lines(p2 ~ time, data = pdat)
lines(upper ~ time, data = pdat, lty = "dashed")
lines(lower ~ time, data = pdat, lty = "dashed")
lines(unlist(m1.dsig$incr) ~ time, data = pdat, col = "blue", lwd = 3)
lines(unlist(m1.dsig$decr) ~ time, data = pdat, col = "red", lwd = 3)

# Extract gam lines
newdat <- data.frame(time=seq(8,27,0.1)) 
p2 <- predict(m1, newdata=newdat, type = "response", se.fit=TRUE)
newdat$LLRmin <- p2$fit
newdat$LLRmin_se <- p2$se.fit
newdat$LLR_upr <- p2$fit + (1.96 * p2$se.fit)
newdat$LLR_lwr <- p2$fit - (1.96 * p2$se.fit)
newdat$Treatment <- "FUT+HW"
signi <- newdat
signi <- subset(signi, time < 27.1)
signi <- subset(signi, time > 24.9)
newdat_fut_hw <- newdat
signi_fut_hw <- signi

## Plot both HW together
hw_time <- ggplot(hw, aes(x=time, y=LRR, color = Treatment, shape = Treatment)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_rect(data=NULL,aes(xmin=10,xmax=16,ymin=-Inf,ymax=Inf),
            fill="gold", alpha = .5)+
  geom_point(position=position_dodge(0.05), size = 3, alpha = 0.8) +
  geom_line(data = newdat_amb, aes(time, LLRmin), size = 1) +
  geom_line(data = newdat_amb, aes(time, LLR_upr), size = .5, alpha = 0.8, linetype="dotted") +
  geom_line(data = newdat_amb, aes(time, LLR_lwr), size = .5, alpha = 0.8, linetype="dotted") +
  geom_line(data = newdat_fut_hw, aes(time, LLRmin), size = 1) +
  geom_line(data = signi_fut_hw, aes(time, LLRmin), size = 3, color = "firebrick1") +
  geom_line(data = newdat_fut_hw, aes(time, LLR_upr), size = .5, alpha = 0.8, linetype="dotted") +
  geom_line(data = newdat_fut_hw, aes(time, LLR_lwr), size = .5, alpha = 0.8, linetype="dotted") +
  plot.theme +
  labs(x="Incubation time (d)", y=bquote("LRR")) + 
  scale_color_manual(values=treat_pal) +
  scale_shape_manual(values=treat_pch) +
  coord_cartesian(ylim = c(-0.5, 0.3)) +
  scale_x_continuous(breaks = seq(0, 27, 3))

hw_time
ggsave("Output/PP_HW_LRR_rich.png", hw_time, height = 4, width = 8, dpi = 320)


#### CALCULATE & PLOT & ANALYSE EVENNESS LRRs ####
## Calculate log response ratios
## FUT
means_amb <- amb %>%
  group_by(time) %>%
  get_summary_stats(Evenness, type = "mean")

fut_m <- merge(fut, means_amb, by = "time")
fut_m$LRR <- log(fut_m$Evenness/fut_m$mean)

## HW-A
hwa_m <- merge(hwa, means_amb, by = "time")
hwa_m$LRR <- log(hwa_m$Evenness/hwa_m$mean)

## HW-F
means_fut <- fut %>%
  group_by(time) %>%
  get_summary_stats(Evenness, type = "mean")

hwf_m <- merge(hwf, means_fut, by = "time")
hwf_m$LRR <- log(hwf_m$Evenness/hwf_m$mean)

## Merged HW dataframe fr plotting
hw <- rbind(hwa_m, hwf_m)

### Statistics for significant in- and decreasees
# Model GAM and then use first derivatives to detect significant changes
# After Gavin Simson

### FUT
m1 <- gam(LRR ~ s(time, k=3, fx= TRUE), data = fut_m)
m1$aic # knot number to decrease aic value
summary(m1)
plot(m1, residuals = TRUE, pch = 19, cex = 0.75)

want <- seq(1, nrow(fut_m), length.out = 200)
pdat <- with(fut_m,
             data.frame(time = time[want]))
p2 <- predict(m1, newdata = pdat, type = "terms", se.fit = TRUE)
pdat <- transform(pdat, p2 = p2$fit[,1], se2 = p2$se.fit[,1])

df.res <- df.residual(m1)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat <- transform(pdat,
                  upper = p2 + (crit.t * se2),
                  lower = p2 - (crit.t * se2))

Term <- "time"
m1.d <- Deriv(m1)
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p2, d = m1.d[[Term]]$deriv,
                   +                    m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)

# Plot to see siginificant changes and manually add to own plots
ylim <- with(pdat, range(upper, lower, p2))
ylab <- expression(LRR)

plot(p2 ~ time, data = pdat, type = "n", ylab = ylab, ylim = ylim)
lines(p2 ~ time, data = pdat)
lines(upper ~ time, data = pdat, lty = "dashed")
lines(lower ~ time, data = pdat, lty = "dashed")
lines(unlist(m1.dsig$incr) ~ time, data = pdat, col = "blue", lwd = 3)
lines(unlist(m1.dsig$decr) ~ time, data = pdat, col = "red", lwd = 3)
# significant decrease from day 0 to 15

## Plot prettier for paper
# Extract gam lines
newdat <- data.frame(time=seq(0,27,0.1)) 
p2 <- predict(m1, newdata=newdat, type = "response", se.fit=TRUE)
newdat$LLRmin <- p2$fit
newdat$LLRmin_se <- p2$se.fit
newdat$LLR_upr <- p2$fit + (1.96 * p2$se.fit)
newdat$LLR_lwr <- p2$fit - (1.96 * p2$se.fit)
newdat$Treatment <- "FUT"
signi <- newdat
signi <- subset(signi, time < 15.1)
head(newdat)

# plot
fut_time <- ggplot(fut_m, aes(x=time, y=LRR, size = sig, color = Treatment)) + 
  geom_hline(yintercept=0, linetype="dashed")+
  geom_point(position=position_dodge(0.05), size = 3, alpha = 0.8) +
  geom_line(data = newdat, aes(time, LLRmin), size = 1) +
  geom_line(data = signi, aes(time, LLRmin), size = 2, color = "indianred1") +
  geom_line(data = newdat, aes(time, LLR_upr), size = .5, alpha = 0.8, linetype="dotted") +
  geom_line(data = newdat, aes(time, LLR_lwr), size = .5, alpha = 0.8, linetype="dotted") +
  plot.theme +
  labs(x="Incubation time (d)", y=bquote("LRR")) + 
  scale_color_manual(values=treat_pal) +
  coord_cartesian(ylim = c(-0.5, 0.3)) +
  scale_x_continuous(breaks = seq(0, 27, 3))

fut_time
ggsave("Output/PP_FUT_LRR_even.png", fut_time, height = 4, width = 8, dpi = 320)


### AMB-HW
m1 <- gam(LRR ~ s(time, k=5, fx= TRUE), data = hwa_m)
m1$aic
summary(m1)
plot(m1, residuals = TRUE, pch = 19, cex = 0.75)

want <- seq(1, nrow(hwa_m), length.out = 200)
pdat <- with(hwa_m,
             data.frame(time = time[want]))
p2 <- predict(m1, newdata = pdat, type = "terms", se.fit = TRUE)
pdat <- transform(pdat, p2 = p2$fit[,1], se2 = p2$se.fit[,1])
df.res <- df.residual(m1)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat <- transform(pdat,
                  upper = p2 + (crit.t * se2),
                  lower = p2 - (crit.t * se2))

Term <- "time"
m1.d <- Deriv(m1)
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p2, d = m1.d[[Term]]$deriv,
                   +                    m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)

# Plot to see sginificant changes and manually add to own plots
ylim <- with(pdat, range(upper, lower, p2))
ylab <- expression(LRR)

plot(p2 ~ time, data = pdat, type = "n", ylab = ylab, ylim = ylim)
lines(p2 ~ time, data = pdat)
lines(upper ~ time, data = pdat, lty = "dashed")
lines(lower ~ time, data = pdat, lty = "dashed")
lines(unlist(m1.dsig$incr) ~ time, data = pdat, col = "blue", lwd = 3)
lines(unlist(m1.dsig$decr) ~ time, data = pdat, col = "red", lwd = 3)

# Extract gam lines
newdat <- data.frame(time=seq(8,27,0.1)) 
p2 <- predict(m1, newdata=newdat, type = "response", se.fit=TRUE)
newdat$LLRmin <- p2$fit
newdat$LLRmin_se <- p2$se.fit
newdat$LLR_upr <- p2$fit + (1.96 * p2$se.fit)
newdat$LLR_lwr <- p2$fit - (1.96 * p2$se.fit)
newdat$Treatment <- "AMB+HW"
signi <- newdat
signi <- subset(signi, time < 15.1)
signi <- subset(signi, time > 10.9)
newdat_amb <- newdat
signi_amb1 <- signi
signi <- newdat
signi <- subset(signi, time < 22.1)
signi <- subset(signi, time > 17.9)
signi_amb2 <- signi

### FUT-HW
m1 <- gam(LRR ~ s(time, k=5, fx= TRUE), data = hwf_m)
m1$aic
summary(m1)
plot(m1, residuals = TRUE, pch = 19, cex = 0.75)

want <- seq(1, nrow(hwf_m), length.out = 200)
pdat <- with(hwf_m,
             data.frame(time = time[want]))

p2 <- predict(m1, newdata = pdat, type = "terms", se.fit = TRUE)
pdat <- transform(pdat, p2 = p2$fit[,1], se2 = p2$se.fit[,1])
df.res <- df.residual(m1)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat <- transform(pdat,
                  upper = p2 + (crit.t * se2),
                  lower = p2 - (crit.t * se2))

Term <- "time"
m1.d <- Deriv(m1)
m1.dci <- confint(m1.d, term = Term)
m1.dsig <- signifD(pdat$p2, d = m1.d[[Term]]$deriv,
                   +                    m1.dci[[Term]]$upper, m1.dci[[Term]]$lower)

# Plot to see sginificant changes and manually add to own plots
ylim <- with(pdat, range(upper, lower, p2))
ylab <- expression(LRR)

plot(p2 ~ time, data = pdat, type = "n", ylab = ylab, ylim = ylim)
lines(p2 ~ time, data = pdat)
lines(upper ~ time, data = pdat, lty = "dashed")
lines(lower ~ time, data = pdat, lty = "dashed")
lines(unlist(m1.dsig$incr) ~ time, data = pdat, col = "blue", lwd = 3)
lines(unlist(m1.dsig$decr) ~ time, data = pdat, col = "red", lwd = 3)

# Extract gam lines
newdat <- data.frame(time=seq(8,27,0.1)) 
p2 <- predict(m1, newdata=newdat, type = "response", se.fit=TRUE)
newdat$LLRmin <- p2$fit
newdat$LLRmin_se <- p2$se.fit
newdat$LLR_upr <- p2$fit + (1.96 * p2$se.fit)
newdat$LLR_lwr <- p2$fit - (1.96 * p2$se.fit)
newdat$Treatment <- "FUT+HW"
signi <- newdat
signi <- subset(signi, time < 18.1)
signi <- subset(signi, time > 12.9)
newdat_fut_hw <- newdat
signi_fut_hw <- signi

## Plot both HW together
hw_time <- ggplot(hw, aes(x=time, y=LRR, color = Treatment, shape = Treatment)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_rect(data=NULL,aes(xmin=10,xmax=16,ymin=-Inf,ymax=Inf),
            fill="gold", alpha = .5)+
  geom_point(position=position_dodge(0.05), size = 3, alpha = 0.8) +
  geom_line(data = newdat_amb, aes(time, LLRmin), size = 1) +
  geom_line(data = signi_amb1, aes(time, LLRmin), size = 3, color = "royalblue1") +
  geom_line(data = signi_amb2, aes(time, LLRmin), size = 3, color = "royalblue1") +
  geom_line(data = newdat_amb, aes(time, LLR_upr), size = .5, alpha = 0.8, linetype="dotted") +
  geom_line(data = newdat_amb, aes(time, LLR_lwr), size = .5, alpha = 0.8, linetype="dotted") +
  geom_line(data = newdat_fut_hw, aes(time, LLRmin), size = 1) +
  geom_line(data = signi_fut_hw, aes(time, LLRmin), size = 3, color = "firebrick1") +
  geom_line(data = newdat_fut_hw, aes(time, LLR_upr), size = .5, alpha = 0.8, linetype="dotted") +
  geom_line(data = newdat_fut_hw, aes(time, LLR_lwr), size = .5, alpha = 0.8, linetype="dotted") +
  plot.theme +
  labs(x="Incubation time (d)", y=bquote("LRR")) + 
  scale_color_manual(values=treat_pal) +
  scale_shape_manual(values=treat_pch) +
  coord_cartesian(ylim = c(-0.5, 0.3)) +
  scale_x_continuous(breaks = seq(0, 27, 3))

hw_time
ggsave("Output/PP_HW_LRR_even.png", hw_time, height = 4, width = 8, dpi = 320)


#### BARGRAPHS SPECIES LEVEL ####
## Create a dataframe on species level
species <- phyloseq::tax_glom(ps_merged, "Species")
df <- plot_bar(species, fill = "Species")
df2 <- df$data
species <- df2 %>% select(Sample, Species, Abundance, Treatment, time, replicate)

# Prepare dataframe for plotting
species$Species[species$Abundance < 200] <- "Other"

species <- species %>%
  mutate(Treatment = recode(Treatment, 'RCP' = 'FUT')) %>%
  mutate(Treatment = recode(Treatment, 'RCP+HW' = 'FUT+HW')) %>%
  mutate(Treatment = recode(Treatment, 'Ambient' = 'AMB')) %>%
  mutate(Treatment = recode(Treatment, 'Ambient+HW' = 'AMB+HW'))

species$Species <- sub("_", " ", species$Species)
species$Species <- sub("X_", "", species$Species)
species$Species <- sub("XXX_", "", species$Species)
species$Species <- sub("XX", "", species$Species)

species <- species %>%
  mutate(Species = recode(Species, "Dolichomastigaceae-B sp." = 'Dolichomastigaceae indet.')) %>%
  mutate(Species = recode(Species, "MOCH-3 sp." = 'Marine ochrophyte indet.')) %>%
  mutate(Species = recode(Species, "Prasino-Clade-VIII sp." = 'Prasino-Clade indet.')) %>%
  mutate(Species = recode(Species, "Micromonas commoda_A2" = 'Micromonas commoda'))

## Create color palette
spe_pal <- qualpal(30, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

## Plotting with and without legend
species_plot <- ggplot(species, aes(fill = Species, x = time, y = Abundance)) +
  facet_wrap(~ Treatment, ncol = 2) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_fill_manual(values = spe_pal$hex)

species_plot
ggsave("Output/PPSpecies.png", species_plot, height = 8, width = 14, dpi = 320)

#### PACKAGES ####
## Check which packages you actually used for documentation and tidying purposes
#packages <- NCmisc::list.functions.in.file("~/AWI/RProjects/HNT22/Phototrophs.R", alphabetic = TRUE) # set to your filepath
#summary(packages)

## Download citation report and citations
#pacman::p_load(grateful)
#cite_packages(out.format = "docx", out.dir = "~/AWI/RProjects/HNT22/Data", pkgs = "Session", out.file = "GDA_packages")
