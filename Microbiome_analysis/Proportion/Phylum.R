
#Proportion of each species across samples was calculated for feature comparison between groups
#Mean of the count number for each species in two groups was calculated for overall compositional comparison
#Read the count table at phylum level
phylum <- read.csv("~/Microbiome_analysis/Normalized_data/level/level-2.csv")

#Sum the number of reads observed across samples for each species
sum <- rowSums(phylum[,c(2:32)])

#Make a new data frame for proportion
phylum_proportion <- phylum[,1]
phylum_proportion <- as.data.frame(phylum_proportion)


#Calculate the proportion for each phylum
for (name in 1:34){
  row <- phylum[name,]
  total <- sum[name]
  for (phy in 2:32){
    value <- row[phy]
    proportion <- value/total
    phylum_proportion[name,phy] <- proportion
  }
}

#Add pre/pro information
phylum_proportion[,33] <- phylum[,33]
colnames(Species_proportion)[33] <- "Pre/Post"

#Save the proportion table
write.csv(phylum_proportion,"~/Microbiome_analysis/Normalized_data/level/proportion_phylum2.csv", row.names = FALSE)



#Calculated the average count number for species in each group(pre/pro)
#Mean for the previous group
Pre_phylum <- colMeans(phylum_proportion[c(1:18),c(2:32)])
#Mean for the post group
Post_phylum <- colMeans(phylum_proportion[c(19:34),c(2:32)])
#Convert it into vector
Pre_phylum <- as.vector(Pre_phylum)
Post_phylum <- as.vector(Post_phylum)

#Make a new table for mean in each group
phylum_mean_table <- cbind(Pre_phylum,Post_phylum)
#Set the column names
colnames(phylum_mean_table) <-c("Pre","Post")
rownames(phylum_mean_table) <- colnames(phylum_proportion)[2:32]

#Save a table for mean proportion
write.csv(phylum_mean_table,"~/Microbiome_analysis/Normalized_data/level/mean_phylum2.csv")




#Compare the proportion at Phylum level between groups (Wilcoxon test)
#Make new table for each group
#For previous group
Pre_table <- phylum_proportion[c(1:18),c(2:32)]
rownames(Pre_table) <- phylum_proportion[1:18,1]
#For post group
Post_table <- phylum_proportion[c(19:34),c(2:32)]
rownames(Post_table) <- phylum_proportion[19:34,1]

#Remove two samples from previous table
Pre_table <- Pre_table[-c(8,9),]
p_value <- list()

#Perform wilcoxon test to identify differently abundant features
for (i in 1:31){
  Pre <- Pre_table[,i]
  Post <- Post_table[,i]
  w <- wilcox.test(Pre,Post, pared = TRUE, alternative = "two.sided")
  p <- w$p.value
  if (w$p.value <= 0.05) {
    print(paste(colnames(Pre_table)[i],p))
  }
  p_value <- append(p_value,p)
}

#Convert into a vector
p_value <- as.vector(p_value)

#Calculate adjusted p values
adjusted_p <- p.adjust(p_value, "bonferroni")

#Detect which feature is significantly different(adjusted p-value < 0.05)
for (i in 1:31){
  p <- adjusted_p[i]
  if (p <= 0.05) {
    print(paste(colnames(Pre_table)[i], p))
  }
}



