
#Proportion of each species across samples was calculated for feature comparison between groups
#Mean of the count number for each species in two groups was calculated for overall compositional comparison
#Read the count table at class level
class <- read.csv("~/Microbiome_analysis/Normalized_data/level/level-3.csv")

#Sum the number of reads observed across samples for each species
sum <- rowSums(class[,c(2:62)])

#Make a new data frame for proportion
class_proportion <- class[,1]
class_proportion <- as.data.frame(class_proportion)


#Calculate the proportion for each class
for (name in 1:34){
  row <- class[name,]
  total <- sum[name]
  for (cla in 2:62){
    value <- row[cla]
    proportion <- value/total
    class_proportion[name,cla] <- proportion
  }
}

#Add pre/pro information
class_proportion[,63] <- class[,63]

#Save the proportion table
write.csv(class_proportion,"~/Microbiome_analysis/Normalized_data/level/proportion_class3.csv", row.names = FALSE)




#Calculated the average count number for species in each group(pre/pro)
#Mean for the previous group
Pre_class <- colMeans(class_proportion[c(1:18),c(2:62)])
#Mean for the post group
Post_class <- colMeans(class_proportion[c(19:34),c(2:62)])
#Convert it into vector
Pre_class <- as.vector(Pre_class)
Post_class <- as.vector(Post_class)

#Make a new table for mean in each group
class_mean_table <- cbind(Pre_class,Post_class)
#Set the column names
colnames(class_mean_table) <-c("Pre","Post")
rownames(class_mean_table) <- colnames(class_proportion)[2:62]

#Save a table for mean proportion
write.csv(class_mean_table,"~/Microbiome_analysis/Normalized_data/level/mean_class3.csv")




#Compare the proportion at Class level between groups (Wilcoxon test)
#Make new table for each group
#For previous group
Pre_table <- class_proportion[c(1:18),c(2:62)]
rownames(Pre_table) <- class_proportion[1:18,1]
#For post group
Post_table <- class_proportion[c(19:34),c(2:62)]
rownames(Post_table) <- class_proportion[19:34,1]

#Remove two samples from previous table
Pre_table <- Pre_table[-c(8,9),]
p_value <- list()

#Perform wilcoxon test to identify differently abundant features
for (i in 1:61){
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
for (i in 1:61){
  p <- adjusted_p[i]
  if (p <= 0.05) {
    print(paste(i,colnames(Pre_table)[i], p))
  }
}

