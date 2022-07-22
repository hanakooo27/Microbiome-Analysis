
#Proportion of each species across samples was calculated for feature comparison between groups
#Mean of the count number for each species in two groups was calculated for overall compositional comparison
#Read the count table at species level
Species <- read.csv("~/Microbiome_analysis/Normalized_data/level/level-7.csv")

#Sum the number of reads observed across samples for each species
sum <- rowSums(Species[,c(2:634)])

#Make a new data frame for proportion
Species_proportion <- Species[,1]
Species_proportion <- as.data.frame(Species_proportion)


#Calculate the proportion for each Species
for (name in 1:34){
  row <- Species[name,]
  total <- sum[name]
  for (spe in 2:634){
    value <- row[spe]
    proportion <- value/total
    Species_proportion[name,spe] <- proportion
  }
}

#Add pre/pro information
Species_proportion[,635] <- Species[,635]
colnames(Species_proportion)[635] <- "Pre/Post"

#Save the proportion table
write.csv(Species_proportion,"~/Microbiome_analysis/Normalized_data/level/proportion_Species7.csv", row.names = FALSE)



#Calculated the average count number for species in each group(pre/pro)
#Mean for the previous group
Pre_Species <- colMeans(Species_proportion[c(1:18),c(2:634)])
#Mean for the post group
Post_Species <- colMeans(Species_proportion[c(19:34),c(2:634)])
#Convert it into vector
Pre_Species <- as.vector(Pre_Species)
Post_Species <- as.vector(Post_Species)

#Make a new table for mean in each group
Species_mean_table <- cbind(Pre_Species,Post_Species)
#Set the column names
colnames(Species_mean_table) <-c("Pre","Post")
rownames(Species_mean_table) <- colnames(Species_proportion)[2:634]

#Save a table for mean proportion
write.csv(Species_mean_table,"~/Microbiome_analysis/Normalized_data/level/mean_Species7.csv")




#Compare the proportion at Species level between groups (Wilcoxon test)
#Make new table for each group
#For previous group
Pre_table <- Species_proportion[c(1:18),c(2:634)]
rownames(Pre_table) <- Species_proportion[1:18,1]
#For post group
Post_table <- Species_proportion[c(19:34),c(2:634)]
rownames(Post_table) <- Species_proportion[19:34,1]

#Remove two samples from previous table
Pre_table <- Pre_table[-c(8,9),]
p_value <- list()

#Perform wilcoxon test to identify differently abundant features
for (i in 1:633){
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
adjusted_p <- p.adjust(p_value, "bonferroni", length(p_value))

#Detect which feature is significantly different(adjusted p-value < 0.05)
for (i in 1:633){
  p <- adjusted_p[i]
  if (p <= 0.05) {
    print(paste(i,colnames(Pre_table)[i], p))
  }
}

