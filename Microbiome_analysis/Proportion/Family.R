
#Proportion of each species across samples was calculated for feature comparison between groups
#Mean of the count number for each species in two groups was calculated for overall compositional comparison
#Read the count table at family level
Family <- read.csv("~/Microbiome_analysis/Normalized_data/level/level-5.csv")

#Sum the number of reads observed across samples for each species
sum <- rowSums(Family[,c(2:251)])

#Make a new data frame for proportion
Family_proportion <- Family[,1]
Family_proportion <- as.data.frame(Family_proportion)


#Calculate the proportion for each Species Family
for (name in 1:34){
  row <- Family[name,]
  total <- sum[name]
  for (fam in 2:251){
    value <- row[fam]
    proportion <- value/total
    Family_proportion[name,fam] <- proportion
  }
}

#Add pre/pro information
Family_proportion[,252] <- Family[,252]
colnames(Family_proportion)[252] <- "Pre/Post"

#Save the proportion table
write.csv(Family_proportion,"~/Microbiome_analysis/Normalized_data/level/proportion_Family5.csv", row.names = FALSE)



#Calculated the average count number for species in each group(pre/pro)
#Mean for the previous group
Pre_Family <- colMeans(Family_proportion[c(1:18),c(2:251)])
#Mean for the post group
Post_Family <- colMeans(Family_proportion[c(19:34),c(2:251)])
#Convert it into vector
Pre_Family <- as.vector(Pre_Family)
Post_Family <- as.vector(Post_Family)

#Make a new table for mean in each group
Family_mean_table <- cbind(Pre_Family,Post_Family)
#Set the column names
colnames(Family_mean_table) <-c("Pre","Post")
rownames(Family_mean_table) <- colnames(Family_proportion)[2:251]

#Save a table for mean proportion
write.csv(Family_mean_table,"~/Microbiome_analysis/Normalized_data/level/mean_Family5.csv")




#Compare the proportion at Family level between groups (Wilcoxon test)
#Make new table for each group
#For previous group
Pre_table <- Family_proportion[c(1:18),c(2:251)]
rownames(Pre_table) <- Family_proportion[1:18,1]
#For post group
Post_table <- Family_proportion[c(19:34),c(2:251)]
rownames(Post_table) <- Family_proportion[19:34,1]

#Remove two samples from previous table
Pre_table <- Pre_table[-c(8,9),]
p_value <- list()

#Perform wilcoxon test to identify differently abundant features
for (i in 1:250){
  Pre <- Pre_table[,i]
  Post <- Post_table[,i]
  w <- wilcox.test(Pre,Post, pared = TRUE, alternative = "two.sided")
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
for (i in 1:250){
  p <- adjusted_p[i]
  if (p <= 0.05) {
    print(paste(i,colnames(Pre_table)[i], p))
  }
}

