#!/bin/Rscript


args<-commandArgs(trailingOnly=TRUE)

library.path <- .libPaths()
library("tidyverse", lib.loc = library.path)
#Read in sample heterozygosity table
sample_het <- read.table(args[1], header = TRUE)
#Read in sample call rate table
sample_callrate <- read.table(args[2], header = TRUE)

#Calculate proportion of heterozygosity
sample_het <- sample_het %>%
  mutate(het = (N.NM. - O.HOM.)/N.NM.)

#Calculate mean and SD of heterozygosity
summary <- sample_het %>%
  summarise(mean_het = mean(het),
            sd_het = sd(het))

mean_het <- summary[1,1]
sd_het <- summary[1,2]

#Write list of samples who are > 2 SDs from mean of heterozygosity
sample_het <- sample_het %>%
  mutate(remove_het = ifelse(het > 3*sd_het + mean_het, "remove",
                             ifelse(het < mean_het - 3*sd_het, "remove", "keep")))


#Merge with callrate table
sample_stats <- sample_het %>%
  left_join(sample_callrate, by = c("FID", "IID"))

#Calculate genotyping rate (1 minus missing rate)
sample_stats <- sample_stats %>%
  mutate(callrate = F_MISS)

#Plot scatterplot
ggplot(data = sample_stats, mapping = aes(x = het, y = callrate, color = remove_het)) +
  geom_point() +
  theme_bw() +
  ggsave("het_imiss.pdf")

#Write list of samples to remove - if heterozygosity outliers or if callrate is <98%
samples_to_remove <- sample_stats %>%
  filter(remove_het == "remove" | callrate > 0.015) %>%
  select(FID, IID)

#Export as text file with just FID and IID
write.table(samples_to_remove, "samples_to_remove.txt",
            quote=F, col.names = F, row.names = F)