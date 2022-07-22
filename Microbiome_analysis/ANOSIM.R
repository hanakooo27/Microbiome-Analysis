
library(vegan)

#Read the count table
count_table <- read.csv("~/Microbiome_analysis/Normalized_data/0.1filtered_perm_Normalized_table.csv")
#Transpose the table
count_table <- t(count_table)
#Set the column name 
colnames(count_table) <- count_table[1,]
#Extract taxon names
taxo <- count_table[36,]
#Remove columns for indices and taxon names
count_table <- count_table[-c(1,36),]

#Make a list for pre/post
pre <- rep("pre", time=18)
post <- rep("post", time=16)
pre_post <- c(pre,post)
#Convert it as a data frame
pre_post <- as.data.frame(pre_post)
#Set the column and row name
colnames(pre_post) <- "pre_post"
rownames(pre_post) <- rownames(count_table)

#Make a list for age (Age was categorized into three groups:1 =<1, 2 1-2, 3 2=< years old)
#For the previous group
Age <- c(3,1,1,1,3,1,1,2,2,2,2,3,3,2,2,1,3)
Age <- as.data.frame(Age)
colnames(Age) <- "Age"
#For the post group
Age2 <- c(3,1,1,1,3,1,2,2,2,3,3,2,2,1,3)
Age2 <- as.data.frame(Age2)
colnames(Age2) <- "Age"  

#Convert all the values as numeric in the count table
count_table2 <- matrix(as.numeric(count_table), ncol = ncol(count_table))
#Set the column and row names
colnames(count_table2) <- colnames(count_table)
rownames(count_table2) <- rownames(count_table)
#
sum <- rowSums(count_table2[,c(1:633)])

#Make a new data frame for concatenating the proportion of each species
taxo_proportion <- rownames(count_table)
taxo_proportion <- as.data.frame(taxo_proportion)

#Calculate the proportion for each species
for (name in 1:34){
  row <- count_table2[name,]
  total <- sum[name]
  for (taxo in 1:633){
    value <- row[taxo]
    proportion <- value/total
    taxo_proportion[name,taxo] <- proportion
  }
}

#Set column names
colnames(taxo_proportion)[2:634] <- colnames(count_table2)

#Make new tables for previous and post groups (Here "Belle" was removed as its age was unavailable)
taxo_proportion <- cbind(pre_post,taxo_proportion)
pre_taxo_proportion <- taxo_proportion[c(1:2,4:18),]
pre_taxo_proportion <- cbind(Age,pre_taxo_proportion)
post_taxo_proportion <- taxo_proportion[c(19:20,22:34),]
post_taxo_proportion <- cbind(Age2,post_taxo_proportion)


#ANOSIM between pre-post groups
#Make community matrix
com = taxo_proportion[,2:ncol(taxo_proportion)]
#Convert it into matrix
m_com = as.matrix(com)
#Perform ANOSIM
ano = anosim(m_com, taxo_proportion$pre_post, distance = "bray", permutations = 9999)
ano

#ANOSIM among different age groups in the previous samples
#Make community matrix
com2 = pre_taxo_proportion[,4:ncol(pre_taxo_proportion)]
#Convert it into matrix
m_com2 = as.matrix(com2)
#Perform ANOSIM
ano2 = anosim(m_com2, pre_taxo_proportion$Age, distance = "bray", permutations = 9999)
ano2

#ANOSIM among different age groups in the previous samples
#Make community matrix
com3 = post_taxo_proportion[,4:ncol(post_taxo_proportion)]
#Convert it into matrix
m_com3 = as.matrix(com3)
#Perform ANOSIM
ano3 = anosim(m_com3, post_taxo_proportion$Age, distance = "bray", permutations = 9999)
ano3
