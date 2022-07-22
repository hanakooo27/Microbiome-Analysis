
source("~/Microbiome_analysis/PC_corr_function.R")

#Manipulate the count table for PC-corr analysis
#Read the count table
count_table_original <- read.csv("~/Microbiome_analysis/Normalized_data/0.1filtered_perm_Normalized_table.csv")
#Set row names
rownames(count_table_original) <- count_table_original[,1]
#Remove columns for indices and taxon names
count_table <- count_table_original[,-c(1,36)]
#Transpose the table
count_table <- t(count_table)

#Make a list for pre/post
pre <- rep("pre",time=18)
post <- rep("post",time=16)
labels <- c(pre,post)
labels <- as.data.frame(labels)



#make a new table for analysis
#Read the original count matrix with taxon level prefixed names with 36 samples
count_table_level <- read.csv("~/Microbiome_analysis/OTU_data/count_table_level.csv")

#Match the indices 
match <- match(count_table_original[,1],count_table_level[,1])
#Select the indices remained in the normalized and filtered table with 34 samples
selected <- count_table_level[match,]
#Remove two samples
selected <- selected[,-c(27,28)]
#Make a new table 
level_count_table <- cbind(selected[,1],count_table_original[,2:35])
level_count_table <- cbind(level_count_table,selected[,36])
#Set column names
colnames(level_count_table)[1] <- "#OTU ID"
colnames(level_count_table)[36] <- "Taxon"
#Save the new table with taxonomic level prefixed names
write.csv(level_count_table,"./OTU_data/level_count_table.csv", row.names = FALSE)



#Read the table
taxo <- read.csv("~/Microbiome_analysis/OTU_data/level_count_table.csv")
#Make a list for species with taxonomic level prefixed names
list_names <- taxo[,36]
#Extract names at the last level of the species
feat_names <- gsub(".*\\;", "", list_names) 
feat_names2 <- sub(" ","",feat_names)
feat_names2 <- as.data.frame(feat_names2)
rownames(feat_names2) <- taxo[,1]
#Save the table for names at the last level of the species
write.csv(feat_names2,"~/Dogs by name/OTU_data/taxo_list.csv")



#Read the feature names 
feat_names2 <- read.csv("~/Microbiome_analysis/OTU_data/taxo_list.csv", row.names=1)
#Get the length of the feature names
len <- length(table(feat_names2))
#Make an empty list for 
list <- list()

#Extract names that appeared more than twice i.e. overlapping 
for (i in 1:len) {
  num <- table(feat_names2)[i] >= 2
  if (num == "TRUE") {
    print(num)
    list <- append(list,table(feat_names2)[i])
  }
}

#Convert it into data frame
list <- as.data.frame(list)
#Get the length of the overlapping names
len2 <- length(list)

#Make a list
index <- list()

#Make the list of indices that have overlapping names
for (i in 1:len2){
  name <- colnames(list)[i]
  for (n in 1:633){
    if (name == feat_names2[n,1]){
      index <- append(index,n)
    }
  }
}

#Extract full names that have overlapping feature names
for (i in index){
  print(taxo[i,c(1,36)])
}


#Manually edited the overlapping feature names in the excel (name at the second last taxonomic level was added)
#Edited table was named as "taxo_list2.csv"
#Read the table
feat_names <- read.csv("~/Microbiome_analysis/OTU_data/taxo_list2.csv")
#Extract feature names
feat_names <- feat_names[,2]
#Convert it into data frame
feat_names <- as.vector(feat_names)

#Extract sample names
samples_names <- rownames(count_table)

#Run PC_corr
PC_corr_v2(count_table,labels,feat_names,samples_names,'yes')


