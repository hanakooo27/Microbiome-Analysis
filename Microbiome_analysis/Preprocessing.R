
# Master data--------------------------------------------------------------------------------------------------------
#Here, we will try to integrate 36 files (18 from each group) into one master data

#In current directory, 18 folders for each dog are saved; folders are named after dogs (e.g. "Asher")
#Each folder consists 2 csv files; one is the previous data, and the other is the post data
#They are named "name Previous.csv" and "name Post.csv" respectively. *name is each dogs name
#Files have following structure;
#There are nine columns(1:7 taxon levels{Kingdom,Phylm,Class,Order,Family,Genus,Species},8 read counts,9 proportion)
#                  rows(names of observed bacteria; The number depends on samples)

#Make a path where the folders are saved
path = "~/Microbiome_analysis/Dogs by name"

#Extract dog names()
dog_names <- list.files(path = path, pattern = NULL, full.names = FALSE)
dog_names

#Read one of files to extract the names of columns for the master data
asher <- read.csv(paste(path,"Asher/Asher Pre.csv",sep=""))
col_list <- colnames(asher)
#Add new columns for the master data
col_list <- c(col_list,"Dog_name","Pre/Post")
#Check the list of column names
print(col_list)

#Make a new data frame for the master data
df <- data.frame(matrix(ncol = length(col_list), nrow = 0))
#Set the column names
colnames(df) <- col_list


library(dplyr)
#Read csv files from each dog and add new columns:Dog_name and Pre/Post
for (dog in dog_names) {
  print(dog)
  df_pre <- read.csv(paste(path, "/", dog, "/",dog, " Pre.csv", sep=""))
  df_post <- read.csv(paste(path, "/", dog, "/",dog, " Post.csv", sep=""))
  df_pre$Dog_name <- dog
  df_pre$"Pre/Post" <- "Pre"
  df_post$Dog_name <- dog
  df_post$"Pre/Post" <- "Post"
  df <- rbind(df, df_pre)
  df <- rbind(df, df_post)
}

#Save the master data
write.csv(df, "~/Microbiome_analysis/master_data.csv", row.names=FALSE)
#This master data consists of following information
#11 columns(1:7 taxon levels,8 read counts,9 proportion,10 Dog_name,11 Pre/Post)
#Rows(accumulated names of observed bacteria from all samples) *overlapped names have not been removed at this point



# Make a list of full species names by combining each name from seven levels into one ------------------------------
#Read the master data
df <- read.csv("~/Microbiome_analysis/master_data.csv")

#Make a copy of "df" data frame 
df_new <- df

#Get the length of data frame (length of taxonomic names accumulated from all files)
len <- dim(df_new)[1]

#Combine names from each taxonomic levels into one
for (x in 1:len) {
  list <- list()
  for (y in 1:7) {
    z <- df_new[x,y]
    check <- z ==""
    if (check == "FALSE") {
      list <- c(list,z)
    }
  }
  df_new$name_bacteria[x] <- knitr::combine_words(list, sep="|", and =" ")
}

#Save the data with the new column, "name_bacteria"
write.csv(df_new, "~/Microbiome_analysis/master_key.csv", row.names=FALSE)



#We will remove overlapping names from the list of accumulated observed bacteria
#Read the master_key data
df_new <- read.csv("~/Microbiome_analysis/master_key.csv")

#Extract index for "pre" data
index_pre <- df_new$`Pre.Post` == "Pre"
#Make a new data frame for "pre" data
df_pre <- df_new[index_pre,]
#Extract index for "post" data
index_post <- df_new$`Pre.Post` == "Post"
#Make a new data frame for "post" data
df_post <- df_new[index_post,]

#Get list of names of bacteria from "pre" data
pre_bacteria_list <- df_pre[,c(1:7,12)]
#Get list of names of bacteria from "post" data
post_bacteria_list <- df_post[,c(1:7,12)]
#Extract unique names
pre_bacteria_list_unique <- unique(pre_bacteria_list)
#Extract unique names
post_bacteria_list_unique <- unique(post_bacteria_list)

#Make a list of names of bacteria for both "pre" and "post"
both_bacteria_list <- rbind(pre_bacteria_list_unique,post_bacteria_list_unique)
#Extract unique names
both_bacteria_list_unique <- unique(both_bacteria_list)

#Make a new data frame for "both"
df_both_final <- data.frame(both_bacteria_list_unique)
#Set column names
colnames(df_both_final) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species", "bacteria_names")

#Extract each dog information and concatenate it into "both" data
for (dog in dog_names) {
  index <- df_pre$Dog_name == dog
  df <- df_pre[index,]
  match_index <- match(df_both_final$bacteria_names,df$name_bacteria)
  df_both_final <- cbind(df_both_final,df[match_index,8:10])
}

for (dog in dog_names) {
  index <- df_post$Dog_name == dog
  df <- df_post[index,]
  match_index <- match(df_both_final$bacteria_names,df$name_bacteria)
  df_both_final <- cbind(df_both_final,df[match_index,8:10])
}

#Save the table for "both"
write.csv(df_both_final, "~/Microbiome_analysis/df_both_final.csv")
#This table consists of following information
#Columns(1:index 2:8 taxon levels,9 full name of species,10:118 {1 read counts,2 proportion,3 Dog_name} x 36samples)
#Rows(names of observed bacteria) 



# Adjust the table format to use it for the downstream analysis including QIIME2 and LEfSe--------------------------
#Read the table
df_both_final <- read.csv("~/Microbiome_analysis/df_both_final.csv")
#Adjust the index
df_both_final[,1] <- rownames(df_both_final)
#Set columns name 
colnames(df_table)[1] <- "#OTU ID"


#Make a list of multiple of 3
three <- list()
for (i in 1:36) {
  num <- i * 3
  print(num)
  three<- c(three, num)
}

#Remove column of "Dog_names"
#Make a list of unnecessary columns
list1 <- list()
for (i in three) {
  num <- i + 9
  names <- colnames(df_both_final[num])
  list1 <- c(list1,names)
}

#Remove columns of "Dog_names"
df_table <- df_both_final[,!(names(df_both_final) %in% list1)]


#Make a list of even number
even <- list()
for (i in 1:36) {
  num <- i * 2
  print(num)
  even<- c(even, num)
}

#Remove column of "proportion"
#make a list of unnecessary columns
list2 <- list()
for (i in even) {
  num <- i + 9
  names <- colnames(df_table[num])
  list2 <- c(list2,names)
}

#Remove column of "proportion" and Make a table for count data
count_table <- df_table[,!(names(df_table) %in% list2)]


#Remove column of "read number"
#make a list of names for unnecessary columns
list3 <- list()
for (i in even) {
  num <- i + 8
  names <- colnames(df_table[num])
  list3 <- c(list3,names)
}

#Remove column of "read number" and Make a table for proportion data
proportion_table <- df_table[,!(names(df_table) %in% list3)]


#Make lists for column names
pre_post <- rep(c("Pre","Post"), each=18)
names <- rep(dog_names, time=2)

#for count table
for (i in 10:45) {
  n <- i-9
  time <- pre_post[n]
  name <- names[n]
  colnames(count_table)[i] <- paste(time, name, sep="_")
  print(colnames(count_table[i]))
}

#for proportion table
for (i in 10:45) {
  n <- i-9
  time <- pre_post[n]
  name <- names[n]
  colnames(proportion_table)[i] <- paste(time, name, sep="_")
  print(colnames(proportion_table[i]))
}

#Now, we will make other two tables including prefix of each taxon level in the bacterial names
#Make new tables 
count_table_level <- count_table
proportion_table_level <- proportion_table

#Extract columns for taxon levels
taxon_levels <- count_table_level[,c(2:8)]

#Get the length of data frame (length of bacterial names)
len <- dim(count_table_level)[1]

#Make a list of prefix for each level
appendex <- c("k_",";p_",";c_",";o_",";f_",";g_",";s_")

#Add prefix of the taxon level and make a new column for the name
#For count table
for (x in 1:len) {
  list <- list()
  for (y in 1:7) {
    t <- taxon_levels[x,y]
    check <- t ==""
    a <- appendex[y]
    if (check == "FALSE") {
      list <- c(list,a,t)
    }
  }
  count_table_level$bacteria_level[x] <- knitr::combine_words(list, sep="", and =" ")
}

#Add prefix of the taxon level and make a new column for the name
#For proportion table
for (x in 1:len) {
  list <- list()
  for (y in 1:7) {
    t <- taxon_levels[x,y]
    check <- t ==""
    a <- appendex[y]
    if (check == "FALSE") {
      list <- c(list,a,t)
    }
  }
  proportion_table_level$bacteria_level[x] <- knitr::combine_words(list, sep="", and =" ")
}


#Remove unnecessary columns (Each taxon names)
count_table <- count_table[,-c(2:8)]
proportion_table <- proportion_table[,-c(2:8)]
count_table_level <- count_table_level[,-c(2:9)]
proportion_table_level <- proportion_table_level[,-c(2:9)]

#Move the column of taxon names to the end
count_table[,39] <- count_table[,2]
count_table <- count_table[,-2]
proportion_table[,39] <- proportion_table[,2]
proportion_table <- proportion_table[,-2]

#Change the name of the new column
colnames(count_OTU_table)[38] <- "Taxon"
colnames(proportion_OTU_table)[38] <- "Taxon"
colnames(count_OTU_table_level)[38] <- "Taxon"
colnames(proportion_OTU_table_level)[38] <- "Taxon"

#Replace NA to 0
count_OTU_table[is.na(count_OTU_table)] <- 0
proportion_OTU_table[is.na(proportion_OTU_table)] <- 0
count_OTU_table_level[is.na(count_OTU_table_level)] <- 0
proportion_OTU_table_level[is.na(proportion_OTU_table_level)] <- 0

#Save the tables 
write.csv(count_table,"~/Microbiome_analysis/OTU_data/count_table.csv", row.names = FALSE)
write.csv(count_table_level,"~/Microbiome_analysis/OTU_data/count_table_level.csv", row.names = FALSE)
write.csv(proportion_table,"~/Microbiome_analysis/OTU_data/proportion_table.csv", row.names = FALSE)
write.csv(proportion_table_level,"~/Microbiome_analysis/OTU_data/proportion_table_level.csv", row.names = FALSE)
#These tables consist of following information
#Columns(1:index 2:37 Samples name, 38 Taxon names)
#Rows(count number for each taxon) 



# Normalization-----------------------------------------------------------------------------------------------------
#Normalization was performed in "SRS.R" script



# After removing two samples----------------------------------------------------------------------------------------
#Count table after removing two samples was named "count_34_table.csv"
#There are rows that have no count across all samples as some species were observed only in those two samples
#We will remove those rows here
#Read table with 34 samples
table_34 <- read.csv("~/Microbiome_analysis/OTU_data/count_34_table.csv")
#Get indies of rows that have count 
non_zero <- rowSums(table_34[,2:35]) != 0
#Remove rows that have no count 
new34_table <- table_34[non_zero,]
#Save the count table with 34 samples
write.csv(new34_table, "~/Microbiome_analysis/OTU_data/count_new34_table.csv", row.names = FALSE)



# 0.1% filtering----------------------------------------------------------------------------------------------------
#Here we will perform 0.1%filtering where features that are observed less than 10 % of totall samples
#Read the count table with 34 samples
new34_table <- read.csv("~/Microbiome_analysis/OTU_data/count_new34_table.csv")
library(Matrix)
#Make a copy of "new34_table"
table <- new34_table
#Set the indices as row name of the table
rownames(table) <- new34_table[,1]
#Count the number of samples that detected each species 
#and sort the species in the order of the observed proportion across samples
abund <- sort(apply(table[,c(2:35)], 1, nnzero), decreasing = TRUE)
#Make it data frame
abund <- as.data.frame(abund)
#Get row names(indices) that have more than 10% observation across samples
taxa <- rownames(abund)[abund/34 > 0.1]
#Filter the features
filtered_table <- table[match(taxa, rownames(table)),]
#Save the 0.1% filtered count table with 34 samples
write.csv(filtered_table,"~/Microbiome_analysis/OTU_data/0.1filter_OTU_table.csv", row.names = FALSE)



# Normalization-----------------------------------------------------------------------------------------------------
#Normalization was performed in "SRS.R" script



# PERFect filtering-------------------------------------------------------------------------------------------------
#Permutation filtering using PERFect was performed in "PERFect.R" script



# Rarefaction curves------------------------------------------------------------------------------------------------
#Rarefaction curves were plotted in "Rarefaction_Curves.R" script



# Correlation-------------------------------------------------------------------------------------------------------
#Correlation was calculated in "Correlation.R" script



# ANOSIM------------------------------------------------------------------------------------------------------------
#ANOSIM was performed in "ANOSIM.R" script



# Wilcoxon singed rank test-----------------------------------------------------------------------------------------
#Proportion of each feature across samples was calculated at 6 taxonomic levels
#Scrips for each level was found in "Proportion" file
#Significantly different features detected here at each level was analyzed and plotted further in Prism.9



# PC-corr analysis--------------------------------------------------------------------------------------------------
#PC-corr was performed in "PC-corr.R" script


