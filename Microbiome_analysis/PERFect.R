library(PERFect)
library(ggplot2)
library(knitr)
library(kableExtra)


#Read the normalized table with 34 samples
normalized_34_table <- read.csv("~/Microbiome_analysis/Normalized_data/0.1filter_Normalized_table.csv")
#Get a list for the taxon names
taxon <- normalized_34_table[,"taxon"]
#Make it a data frame
taxon <- as.data.frame(taxon)
#Set row names as the indices 
rownames(taxon) <- normalized_34_table[,1]

#Transpose the matrix and make a new table
feature_table <- t(normalized_34_table)
#Remove rows of indices and taxon 
feature_table <- feature_table[-c(1,36),]
#Convert the values as numeric
feature_table <- matrix(as.numeric(feature_table), ncol = ncol(feature_table))
#Extract the name of samples
names <- colnames(normalized_34_table)[2:35]
#Set row names of the table 
rownames(feature_table) <- names
#Set column names of the table
colnames(feature_table) <- normalized_34_table[,1]
#Convert the matrix into a data frame
feature_table <- as.data.frame(feature_table)


#Filtering loss function
res_sim <- PERFect_sim(X = feature_table)

#p <- pvals_Plots(PERFect = res_sim, X = feature_table, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)
#p$plot + ggtitle("Simultanenous Filtering")


#Permutation filtering algorithm
res_perm <- PERFect_perm(X = feature_table , Order = "pvals", pvals_sim = res_sim, algorithm = "fast",
                          rollmean = FALSE, k = 20000)

#Plot permutation PERFect p-values
p <- pvals_Plots(res_perm, feature_table)
p$plot + ggtitle("Fast Algorithm")

#Get filtered table
filtered_table <- t(res_perm$filtX)

#Add a column for the indices
filtered_table <- cbind(rownames(filtered_table),filtered_table)
#Set column name
colnames(filtered_table)[1] <- "#OTU ID"

#Add taxon names to the filtered table
#Get the length of the taxon names remained after filtering
len_bacteria <- dim(filtered_table)[1]
#Make an empty list to contain the taxon names
list_names <- list()
#Match the names of species remained in the filtered table
for (i in 1:len_bacteria) {
  match <- rownames(filtered_table)[i] == rownames(taxon)
  list_names <- c(list_names,taxon[match,])
}

#Combine the list of the taxon names with the filtered table
filtered_table <- cbind(filtered_table,list_names)
#Set the column name
colnames(filtered_table)[36] <- "Taxon"

#Save the filtered table
write.csv(filtered_table,"~/Microbiome_analysis/Normalized_data/0.1filtered_perm_Normalized_table.csv",row.names = FALSE)

