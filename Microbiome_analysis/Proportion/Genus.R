
#Proportion of each species across samples was calculated for feature comparison between groups
#Mean of the count number for each species in two groups was calculated for overall compositional comparison
#Read the count table at genus level
Genus <- read.csv("~/Microbiome_analysis/Normalized_data/level/level-6.csv")

#Sum the number of reads observed across samples for each species
sum <- rowSums(Genus[,c(2:478)])

#Make a new data frame for proportion
Genus_proportion <- Genus[,1]
Genus_proportion <- as.data.frame(Genus_proportion)


#Calculate the proportion for each Genus
for (name in 1:34){
  row <- Genus[name,]
  total <- sum[name]
  for (fam in 2:478){
    value <- row[fam]
    proportion <- value/total
    Genus_proportion[name,fam] <- proportion
  }
}

#Add pre/pro information
Genus_proportion[,479] <- Genus[,479]
colnames(Genus_proportion)[479] <- "Pre/Post"

#Save the proportion table
write.csv(Genus_proportion,"~/Microbiome_analysis/Normalized_data/level/proportion_Genus6.csv", row.names = FALSE)



#Calculated the average count number for species in each group(pre/pro)
#Mean for the previous group
Pre_Genus <- colMeans(Genus_proportion[c(1:18),c(2:478)])
#Mean for the post group
Post_Genus <- colMeans(Genus_proportion[c(19:34),c(2:478)])
#Convert it into vector
Pre_Genus <- as.vector(Pre_Genus)
Post_Genus <- as.vector(Post_Genus)

#Make a new table for mean in each group
Genus_mean_table <- cbind(Pre_Genus,Post_Genus)
#Set the column names
colnames(Genus_mean_table) <-c("Pre","Post")
rownames(Genus_mean_table) <- colnames(Genus_proportion)[2:478]

#Save a table for mean proportion
write.csv(Genus_mean_table,"~/Microbiome_analysis/Normalized_data/level/mean_Genus6.csv")




#Compare the proportion at Genus level between groups (Wilcoxon test)
#Make new table for each group
#For previous group
Pre_table <- Genus_proportion[c(1:18),c(2:478)]
rownames(Pre_table) <- Genus_proportion[1:18,1]
#For post group
Post_table <- Genus_proportion[c(19:34),c(2:478)]
rownames(Post_table) <- Genus_proportion[19:34,1]

#Remove two samples from previous table
Pre_table <- Pre_table[-c(8,9),]
p_value <- list()

#Perform wilcoxon test to identify differently abundant features
for (i in 1:477){
  Pre <- Pre_table[,i]
  Post <- Post_table[,i]
  w <- wilcox.test(Pre,Post, paired = TRUE, alternative = "two.sided")
  p <- w$p.value
  if (p <= 0.05) {
    print(paste(colnames(Pre_table)[i],p))
  }
  p_value <- append(p_value,p)
}

#Convert into a vector
p_value <- as.vector(p_value)

#Calculate adjusted p values
adjusted_p <- p.adjust(p_value, "bonferroni")

#Detect which feature is significantly different(adjusted p-value < 0.05)
for (i in 1:477){
  p <- adjusted_p[i]
  if (p <= 0.05) {
    print(paste(i,colnames(Pre_table)[i], p))
  }
}

