#Here we checked the correlation between the number of sequencing and the number of observed species
#"correlation_feature.csv" have following format
#Columns (1 Sample ID,2 observed_features,3 feature_count)
#Rows (36 samples names)

#Read the file
correlation <- read.csv("~/Microbiome_analysis/correlation_feature.csv")
#Set row names as the sample names
rownames(correlation) <- correlation[,1]
#Remove the column of the sample names
correlation <- correlation[,-1]

library("ggpubr")

#Plot the correlation
ggscatter(correlation, x = "observed_features", y = "feature_count", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "number of observed feature", ylab = "number of reads")

