#---------------introduction to scaling with ranked subsampling (SRS)####
#....................................................................................................
#Improved normalization of species count data in ecology by scaling with ranked subsampling (SRS):...
#application to microbial communities................................................................
#....................................................................................................
#Version 1.0 from 2020-05-13.........................................................................
#....................................................................................................
#Authors: Lukas Beule* (lukas.beule@agr.uni-goettingen.de) and Petr Karlovsky* (pkarlov@gwdg.de).....
#....................................................................................................
#*Molecular Phytopathology and Mycotoxin Research, Faculty of Agricultural Sciences,.................
#University of Goettingen, Goettingen, Germany.......................................................
#....................................................................................................
#point of contact: Lukas Beule (lukas.beule@agr.uni-goettingen.de)...................................
#....................................................................................................
#Description of SRS:
#The normalization by SRS reduces the number of counts in each sample in such a way that 
#(i) the total count equals Cmin, (ii) each removed OTU is less or equally abundant that 
#any preserved OTU, and (iii) the relative frequencies of OTUs remaining in the sample 
#after normalization are as close as possible to the frequencies in the original sample. 
#The algorithm consists of two steps. In the first step, the counts for all OTUs 
#are divided by a scaling factor chosen in such a way that the sum of the scaled counts 
#(Cscaled with integer or non-integer values) equals Cmin. The relative frequencies of 
#all OTUs remain unchanged. In the second step, the non-integer count values are converted 
#into integers by an algorithm that we dub ranked subsampling. The scaled count 
#Cscaled for each OTU is split into the integer-part Cint by truncating the digits after 
#the decimal separator (Cint = floor(Cscaled)) and the fractional part Cfrac (Cfrac = Cscaled - Cint). 
#Since ΣCint ≤ Cmin, additional ∆C = Cmin - ΣCint counts have to be added to the library 
#to reach the total count of Cmin. This is achieved as follows. OTUs are ranked in the 
#descending order of their Cfrac values, which lie in the open interval (0, 1). 
#Beginning with the OTU of the highest rank, single count per OTU is added to the 
#normalized library until the total number of added counts reaches ∆C and the sum of all
#counts in the normalized library equals Cmin. For example, if ∆C = 5 and the seven top 
#Cfrac values are 0.96, 0.96, 0.88, 0.55, 0.55, 0.55, and 0.55, the following counts are 
#added: a single count for each OTU with Cfrac of 0.96; a single count for the OTU with Cfrac of
#0.88; and a single count each for two OTUs among those with Cfrac of 0.55. When the lowest Cfrag
#involved in picking ∆C counts is shared by several OTUs, the OTUs used for adding a single
#count to the library are selected in the order of their Cint values. This selection 
#minimizes the effect of normalization on the relative frequencies of OTUs. OTUs with 
#identical Cfrag as well as Cint are sampled randomly without replacement.
#---------------general comments and example data#####
#Samples should be arranged columnwise.
#Input data should not contain any categorial data such as taxonomic assignment and barcode sequences.
#An example of input data can be found below:


#Read the count table 
count_table <- read.csv("~/Microbiome_analysis/OTU_data/count_table.csv")
taxon <- count_table[,38]
table <- count_table[,-c(1,38)]

#In case you detected any errors, please contact the authors (lukas.beule@agr.uni-goettingen.de)
#---------------selection of the desired library size####
#(e.g. species counts of the library with the lowest sequencing depth):
Cmin <- min(colSums(table))
Cmin
data = table
#---------------SRS function####
SRS <- function(data, Cmin){
  if(Cmin > min(colSums(data))){
    print("ERROR: Cmin > minimum library size. Please select a Cmin that is ≤ the minimum library size of the dataset.")
  } else {
    if(Cmin < 0){
      print("ERROR: Cmin < 0. Please select a Cmin ≥ 0.")
    } else {
      if(Cmin %% 1 > 0){
        print("ERROR: Please select a Cmin without decimal places.")
      } else {
        counter = 1 #counting the loops
        for(i in seq(1, ncol(data), 1)){
          if (i == 1) {
            fixed_factor <- (data[,i]/(sum(data[,i])/Cmin))
            assign(paste(names(data)[i],sep=""),fixed_factor)
            fixed_factor_1 <- data.frame(get(names(data)[i]))
            colnames(fixed_factor_1)[i] <- names(data)[i]
          } else {
            fixed_factor <- (data[,i]/(sum(data[,i])/Cmin))
            assign(paste(names(data)[i],sep=""),fixed_factor)
            fixed_factor_1 <- cbind(fixed_factor_1, fixed_factor)
            colnames(fixed_factor_1)[{
              counter = counter + 1
            }] <- names(data)[i]
          }
        }
        
        fixed_factor_1
        
        revtrunc_fixed_factor_1 <- floor(fixed_factor_1) # floor (e.g. 1.9 will become 1)
        revtrunc_fixed_factor_1
        
        diff_counts <- Cmin-colSums(revtrunc_fixed_factor_1) #how many counts differences to the selected library size?
        diff_counts
        
        revtrunc <- function(x) { sign(x) * (x - floor(x)) } 
        revtrunc_fixed_factor <- (round(revtrunc(fixed_factor_1),10000000))
        revtrunc_fixed_factor
        
        x <- as.data.frame(revtrunc_fixed_factor)
        counter = 1
        for(i in seq(1, ncol(x), 1)){
          if (i == 1) {
            if(diff_counts[i] == 0){
              fixed_factor <- revtrunc_fixed_factor_1[,i]
              assign(paste(names(data)[i],sep=""),fixed_factor)
              fixed_factor_1 <- data.frame(get(names(data)[i]))
              colnames(fixed_factor_1)[i] <- names(data)[i]
            } #if the sum of the counts in the library = Cmin
            else {
              maxN <- function(x, N=diff_counts[i]){
                len <- length(x)
                if(N>len){
                  warning('N greater than length(x).  Setting N=length(x)')
                  N <- length(x)
                }
                sort(x,partial=len-N+1)[len-N+1]
              }
              max <- which(x[,i] == maxN(x[,i]), arr.ind = TRUE)
              max
              sum(x[,i] > unique(x[,i][max]))
              normalization_value <- diff_counts[i] - sum(x[,i] > unique(x[,i][max]))
              normalization_value
              
              lowest_level_choise <- as.data.frame(which(x[,i] == unique(maxN(x[,i]))))
              lowest_level_choise 
              length(t(lowest_level_choise)) #how many counts have to be sampled?
              
              if(sum(revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]]) == 0){
                lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise[,1]), normalization_value, replace = F)))
                y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                y[lowest_level] = 1 #set the randomly selected counts to 1
                y
                
              } #if all of the integer values of the lowest rank = 0, do random subsmpling
              else {
                sub_int <- subset(lowest_level_choise, (revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]] >= 1) == TRUE)
                sub_int_bind <- as.data.frame(cbind(sub_int,revtrunc_fixed_factor_1[,i][sub_int[,1]]))
                names(sub_int_bind)[1] <- "V1"
                names(sub_int_bind)[2] <- "V2"
                
                sub_int_bind_ordered <- sub_int_bind[order(sub_int_bind$V2, decreasing = TRUE),]  
                sub_int_bind_ordered
                
                sub_int_bind_ordered_V1 <- sub_int_bind_ordered$V1
                sub_int_bind_ordered_V2 <- sub_int_bind_ordered$V2
                
                if((length(unique(sub_int_bind_ordered_V2)) == 1 & length(sub_int_bind_ordered_V2)>as.vector(normalization_value))){
                  lowest_level <- as.numeric(as.vector(sample(as.factor(sub_int_bind_ordered_V1), normalization_value, replace = F)))
                  y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                  y[lowest_level] = 1 #set the randomly selected counts to 1
                  y
                } #if all of the integer values of the lowest rank are equal, do random subsampling
                else {
                  if(length(sub_int_bind_ordered_V1)>normalization_value){
                    
                    maxN_1 <- function(x, N=normalization_value){
                      len <- length(x)
                      if(N>len){
                        warning('N greater than length(x).  Setting N=length(x)')
                        N <- length(x)
                      }
                      sort(x,partial=len-N+1)[len-N+1]
                    }
                    max_1 <- which(as.data.frame(sub_int_bind_ordered_V2)[,1] == maxN_1(as.data.frame(sub_int_bind_ordered_V2)[,1]), arr.ind = TRUE)
                    max_1
                    sum(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1]))
                    
                    normalization_value_1 <- normalization_value - sum(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1])) # how may values are above the maxima that is the one to be normalized?
                    normalization_value_1
                    lowest_level_choise_1 <- as.data.frame(which(as.data.frame(sub_int_bind_ordered_V2)[,1] == unique(maxN_1(as.data.frame(sub_int_bind_ordered_V2)[,1]))))
                    lowest_level_choise_1
                    lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise_1[,1]), normalization_value_1, replace = F)))
                    lowest_level <- sub_int_bind_ordered_V1[lowest_level]
                    lowest_level
                    
                    lowest_level_1 <- sub_int_bind_ordered_V1[(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1]))]
                    
                    lowest_level <-c(lowest_level_1, lowest_level)
                    y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                    y[lowest_level] = 1 #set the randomly selected counts to 1
                    y
                    
                  } #if integer ranks are > normalization_value, do ranked subsampling
                  else {
                    if(length(sub_int_bind_ordered_V1)<normalization_value){
                      sub_int_zeros <- subset(lowest_level_choise, (revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]] < 1) == TRUE)
                      length(t(sub_int_zeros))
                      
                      lowest_level_2 <- as.numeric(as.vector(sample(as.factor(sub_int_zeros[,1]), (normalization_value-length(sub_int_bind_ordered_V1)), replace = F)))
                      lowest_level_2
                      lowest_level_3 <- c(sub_int_bind_ordered_V1,lowest_level_2)
                      
                      y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                      y[lowest_level_3] = 1 #set the randomly selected counts to 1
                      y
                      
                    } #if integer ranks are < normalization_value, do ranked subsampling of the zero integers
                    else {
                      y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                      y[sub_int_bind_ordered_V1] = 1 #set the randomly selected counts to 1
                      y
                    } # if integer ranks are = normalization_value, sample all
                  }
                  
                  
                }
                
              } 
              
              SRS <- revtrunc_fixed_factor_1[,i] + ceiling(x[,i] > unique(x[,i][max])) + y #sum it all u
              SRS
              sum(SRS) #verification
              assign(paste(names(data)[i],sep=""),SRS)
              fixed_factor_1 <- data.frame(get(names(data)[i]))
              colnames(fixed_factor_1)[i] <- names(data)[i]
            } #if the sum of the counts in the library > Cmin
          } #for the first libraray
          else {
            if(diff_counts[i] == 0){
              fixed_factor <- revtrunc_fixed_factor_1[,i]
              assign(paste(names(data)[i],sep=""),fixed_factor)
              fixed_factor_1 <- cbind(fixed_factor_1, fixed_factor)
              colnames(fixed_factor_1)[{
                counter = counter + 1
              }] <- names(data)[i]
            } #if the sum of the counts in the library = Cmin
            else { 
              maxN <- function(x, N=diff_counts[i]){
                len <- length(x)
                if(N>len){
                  warning('N greater than length(x).  Setting N=length(x)')
                  N <- length(x)
                }
                sort(x,partial=len-N+1)[len-N+1]
              }
              max <- which(x[,i] == maxN(x[,i]), arr.ind = TRUE)
              max
              sum(x[,i] > unique(x[,i][max]))
              normalization_value <- diff_counts[i] - sum(x[,i] > unique(x[,i][max]))
              normalization_value
              
              lowest_level_choise <- as.data.frame(which(x[,i] == unique(maxN(x[,i]))))
              lowest_level_choise 
              length(t(lowest_level_choise)) #how many counts have to be sampled?
              
              if(sum(revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]]) == 0){
                lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise[,1]), normalization_value, replace = F)))
                y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                y[lowest_level] = 1 #set the randomly selected counts to 1
                y
                
              } #if all of the integer values of the lowest rank = 0, do random subsmpling
              else {
                sub_int <- subset(lowest_level_choise, (revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]] >= 1) == TRUE)
                sub_int_bind <- as.data.frame(cbind(sub_int,revtrunc_fixed_factor_1[,i][sub_int[,1]]))
                names(sub_int_bind)[1] <- "V1"
                names(sub_int_bind)[2] <- "V2"
                
                sub_int_bind_ordered <- sub_int_bind[order(sub_int_bind$V2, decreasing = TRUE),]  
                sub_int_bind_ordered
                
                sub_int_bind_ordered_V1 <- sub_int_bind_ordered$V1
                sub_int_bind_ordered_V2 <- sub_int_bind_ordered$V2
                
                if((length(unique(sub_int_bind_ordered_V2)) == 1 & length(sub_int_bind_ordered_V2)>as.vector(normalization_value))){
                  lowest_level <- as.numeric(as.vector(sample(as.factor(sub_int_bind_ordered_V1), normalization_value, replace = F)))
                  y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                  y[lowest_level] = 1 #set the randomly selected counts to 1
                  y
                } #if all of the integer values of the lowest rank are equal, do random subsampling
                else {
                  if(length(sub_int_bind_ordered_V1)>normalization_value){
                    
                    maxN_1 <- function(x, N=normalization_value){
                      len <- length(x)
                      if(N>len){
                        warning('N greater than length(x).  Setting N=length(x)')
                        N <- length(x)
                      }
                      sort(x,partial=len-N+1)[len-N+1]
                    }
                    max_1 <- which(as.data.frame(sub_int_bind_ordered_V2)[,1] == maxN_1(as.data.frame(sub_int_bind_ordered_V2)[,1]), arr.ind = TRUE)
                    max_1
                    sum(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1]))
                    
                    normalization_value_1 <- normalization_value - sum(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1])) # how may values are above the maxima that is the one to be normalized?
                    normalization_value_1
                    lowest_level_choise_1 <- as.data.frame(which(as.data.frame(sub_int_bind_ordered_V2)[,1] == unique(maxN_1(as.data.frame(sub_int_bind_ordered_V2)[,1]))))
                    lowest_level_choise_1
                    lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise_1[,1]), normalization_value_1, replace = F)))
                    lowest_level <- sub_int_bind_ordered_V1[lowest_level]
                    lowest_level
                    
                    lowest_level_1 <- sub_int_bind_ordered_V1[(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1]))]
                    
                    lowest_level <-c(lowest_level_1, lowest_level)
                    y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                    y[lowest_level] = 1 #set the randomly selected counts to 1
                    y
                    
                  } #if integer ranks are > normalization_value, do ranked subsampling
                  else {
                    if(length(sub_int_bind_ordered_V1)<normalization_value){
                      sub_int_zeros <- subset(lowest_level_choise, (revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]] < 1) == TRUE)
                      length(t(sub_int_zeros))
                      
                      lowest_level_2 <- as.numeric(as.vector(sample(as.factor(sub_int_zeros[,1]), (normalization_value-length(sub_int_bind_ordered_V1)), replace = F)))
                      lowest_level_2
                      lowest_level_3 <- c(sub_int_bind_ordered_V1,lowest_level_2)
                      
                      y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                      y[lowest_level_3] = 1 #set the randomly selected counts to 1
                      y
                      
                    } #if integer ranks are < normalization_value, do ranked subsampling of the zero integers
                    else {
                      y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                      y[sub_int_bind_ordered_V1] = 1 #set the randomly selected counts to 1
                      y
                    } # if integer ranks are = normalization_value, sample all
                  }
                  
                  
                }
                
              }    
              
              SRS <- revtrunc_fixed_factor_1[,i] + ceiling(x[,i] > unique(x[,i][max])) + y #sum it all up
              SRS
              sum(SRS) #verification
              assign(paste(names(data)[i],sep=""),SRS)
              fixed_factor_1 <- cbind(fixed_factor_1, SRS)
              colnames(fixed_factor_1)[{
                counter = counter + 1
              }] <- names(data)[i]
            } #if the sum of the counts in the library > Cmin
          } #for all other libaries 
        }
        
        SRS_output <- fixed_factor_1
        SRS_output
      }
    }
  }
}
#---------------running the SRS function on example data####
SRS_output <- SRS(data = OTU_table, Cmin = Cmin)
SRS_output
dim(SRS_output)
#---------------confirmation of equal library size of the example data####
for(i in seq(1,ncol(SRS_output), 1)){
  if(sum(SRS_output[,i]) == Cmin){print(TRUE)} else {print(ERROR)}
}

#Make the normalized table
Normalized_table <- cbind(SRS_output,taxon)
#Add index column 
Normalized_table <- cbind(count_table[,1],Normalized_table)
#Name the index column
colnames(Normalized_table)[1] <- "#OTU ID"
#Identify the rows(Species) that have no count across all samples after the normalization
non_zero <- rowSums(Normalized_OTU_table[,2:37]) != 0
#Eliminate the rows(Species) that have no count across all samples after the normalization
Normalized_new_table <- Normalized_table[non_zero,]
#Save the normalized table
write.csv(Normalized_new_table,"~/Microbiome_analysis/Normalized_data/Normalized_new_table.csv",row.names=FALSE)





#Following normalization was performed after removing two samples

#---------------------------------------------------------------
count_table <- read.csv("~/Microbiome_analysis/OTU_data/0.1filter_table.csv")
taxon <- count_table[,36]
table <- count_table[,-c(1,36)]

#In case you detected any errors, please contact the authors (lukas.beule@agr.uni-goettingen.de)
#---------------selection of the desired library size####
#(e.g. species counts of the library with the lowest sequencing depth):
Cmin <- min(colSums(OTU_table))
Cmin
data = table
#---------------SRS function####
SRS <- function(data, Cmin){
  if(Cmin > min(colSums(data))){
    print("ERROR: Cmin > minimum library size. Please select a Cmin that is ≤ the minimum library size of the dataset.")
  } else {
    if(Cmin < 0){
      print("ERROR: Cmin < 0. Please select a Cmin ≥ 0.")
    } else {
      if(Cmin %% 1 > 0){
        print("ERROR: Please select a Cmin without decimal places.")
      } else {
        counter = 1 #counting the loops
        for(i in seq(1, ncol(data), 1)){
          if (i == 1) {
            fixed_factor <- (data[,i]/(sum(data[,i])/Cmin))
            assign(paste(names(data)[i],sep=""),fixed_factor)
            fixed_factor_1 <- data.frame(get(names(data)[i]))
            colnames(fixed_factor_1)[i] <- names(data)[i]
          } else {
            fixed_factor <- (data[,i]/(sum(data[,i])/Cmin))
            assign(paste(names(data)[i],sep=""),fixed_factor)
            fixed_factor_1 <- cbind(fixed_factor_1, fixed_factor)
            colnames(fixed_factor_1)[{
              counter = counter + 1
            }] <- names(data)[i]
          }
        }
        
        fixed_factor_1
        
        revtrunc_fixed_factor_1 <- floor(fixed_factor_1) # floor (e.g. 1.9 will become 1)
        revtrunc_fixed_factor_1
        
        diff_counts <- Cmin-colSums(revtrunc_fixed_factor_1) #how many counts differences to the selected library size?
        diff_counts
        
        revtrunc <- function(x) { sign(x) * (x - floor(x)) } 
        revtrunc_fixed_factor <- (round(revtrunc(fixed_factor_1),10000000))
        revtrunc_fixed_factor
        
        x <- as.data.frame(revtrunc_fixed_factor)
        counter = 1
        for(i in seq(1, ncol(x), 1)){
          if (i == 1) {
            if(diff_counts[i] == 0){
              fixed_factor <- revtrunc_fixed_factor_1[,i]
              assign(paste(names(data)[i],sep=""),fixed_factor)
              fixed_factor_1 <- data.frame(get(names(data)[i]))
              colnames(fixed_factor_1)[i] <- names(data)[i]
            } #if the sum of the counts in the library = Cmin
            else {
              maxN <- function(x, N=diff_counts[i]){
                len <- length(x)
                if(N>len){
                  warning('N greater than length(x).  Setting N=length(x)')
                  N <- length(x)
                }
                sort(x,partial=len-N+1)[len-N+1]
              }
              max <- which(x[,i] == maxN(x[,i]), arr.ind = TRUE)
              max
              sum(x[,i] > unique(x[,i][max]))
              normalization_value <- diff_counts[i] - sum(x[,i] > unique(x[,i][max]))
              normalization_value
              
              lowest_level_choise <- as.data.frame(which(x[,i] == unique(maxN(x[,i]))))
              lowest_level_choise 
              length(t(lowest_level_choise)) #how many counts have to be sampled?
              
              if(sum(revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]]) == 0){
                lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise[,1]), normalization_value, replace = F)))
                y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                y[lowest_level] = 1 #set the randomly selected counts to 1
                y
                
              } #if all of the integer values of the lowest rank = 0, do random subsmpling
              else {
                sub_int <- subset(lowest_level_choise, (revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]] >= 1) == TRUE)
                sub_int_bind <- as.data.frame(cbind(sub_int,revtrunc_fixed_factor_1[,i][sub_int[,1]]))
                names(sub_int_bind)[1] <- "V1"
                names(sub_int_bind)[2] <- "V2"
                
                sub_int_bind_ordered <- sub_int_bind[order(sub_int_bind$V2, decreasing = TRUE),]  
                sub_int_bind_ordered
                
                sub_int_bind_ordered_V1 <- sub_int_bind_ordered$V1
                sub_int_bind_ordered_V2 <- sub_int_bind_ordered$V2
                
                if((length(unique(sub_int_bind_ordered_V2)) == 1 & length(sub_int_bind_ordered_V2)>as.vector(normalization_value))){
                  lowest_level <- as.numeric(as.vector(sample(as.factor(sub_int_bind_ordered_V1), normalization_value, replace = F)))
                  y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                  y[lowest_level] = 1 #set the randomly selected counts to 1
                  y
                } #if all of the integer values of the lowest rank are equal, do random subsampling
                else {
                  if(length(sub_int_bind_ordered_V1)>normalization_value){
                    
                    maxN_1 <- function(x, N=normalization_value){
                      len <- length(x)
                      if(N>len){
                        warning('N greater than length(x).  Setting N=length(x)')
                        N <- length(x)
                      }
                      sort(x,partial=len-N+1)[len-N+1]
                    }
                    max_1 <- which(as.data.frame(sub_int_bind_ordered_V2)[,1] == maxN_1(as.data.frame(sub_int_bind_ordered_V2)[,1]), arr.ind = TRUE)
                    max_1
                    sum(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1]))
                    
                    normalization_value_1 <- normalization_value - sum(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1])) # how may values are above the maxima that is the one to be normalized?
                    normalization_value_1
                    lowest_level_choise_1 <- as.data.frame(which(as.data.frame(sub_int_bind_ordered_V2)[,1] == unique(maxN_1(as.data.frame(sub_int_bind_ordered_V2)[,1]))))
                    lowest_level_choise_1
                    lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise_1[,1]), normalization_value_1, replace = F)))
                    lowest_level <- sub_int_bind_ordered_V1[lowest_level]
                    lowest_level
                    
                    lowest_level_1 <- sub_int_bind_ordered_V1[(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1]))]
                    
                    lowest_level <-c(lowest_level_1, lowest_level)
                    y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                    y[lowest_level] = 1 #set the randomly selected counts to 1
                    y
                    
                  } #if integer ranks are > normalization_value, do ranked subsampling
                  else {
                    if(length(sub_int_bind_ordered_V1)<normalization_value){
                      sub_int_zeros <- subset(lowest_level_choise, (revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]] < 1) == TRUE)
                      length(t(sub_int_zeros))
                      
                      lowest_level_2 <- as.numeric(as.vector(sample(as.factor(sub_int_zeros[,1]), (normalization_value-length(sub_int_bind_ordered_V1)), replace = F)))
                      lowest_level_2
                      lowest_level_3 <- c(sub_int_bind_ordered_V1,lowest_level_2)
                      
                      y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                      y[lowest_level_3] = 1 #set the randomly selected counts to 1
                      y
                      
                    } #if integer ranks are < normalization_value, do ranked subsampling of the zero integers
                    else {
                      y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                      y[sub_int_bind_ordered_V1] = 1 #set the randomly selected counts to 1
                      y
                    } # if integer ranks are = normalization_value, sample all
                  }
                  
                  
                }
                
              } 
              
              SRS <- revtrunc_fixed_factor_1[,i] + ceiling(x[,i] > unique(x[,i][max])) + y #sum it all u
              SRS
              sum(SRS) #verification
              assign(paste(names(data)[i],sep=""),SRS)
              fixed_factor_1 <- data.frame(get(names(data)[i]))
              colnames(fixed_factor_1)[i] <- names(data)[i]
            } #if the sum of the counts in the library > Cmin
          } #for the first libraray
          else {
            if(diff_counts[i] == 0){
              fixed_factor <- revtrunc_fixed_factor_1[,i]
              assign(paste(names(data)[i],sep=""),fixed_factor)
              fixed_factor_1 <- cbind(fixed_factor_1, fixed_factor)
              colnames(fixed_factor_1)[{
                counter = counter + 1
              }] <- names(data)[i]
            } #if the sum of the counts in the library = Cmin
            else { 
              maxN <- function(x, N=diff_counts[i]){
                len <- length(x)
                if(N>len){
                  warning('N greater than length(x).  Setting N=length(x)')
                  N <- length(x)
                }
                sort(x,partial=len-N+1)[len-N+1]
              }
              max <- which(x[,i] == maxN(x[,i]), arr.ind = TRUE)
              max
              sum(x[,i] > unique(x[,i][max]))
              normalization_value <- diff_counts[i] - sum(x[,i] > unique(x[,i][max]))
              normalization_value
              
              lowest_level_choise <- as.data.frame(which(x[,i] == unique(maxN(x[,i]))))
              lowest_level_choise 
              length(t(lowest_level_choise)) #how many counts have to be sampled?
              
              if(sum(revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]]) == 0){
                lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise[,1]), normalization_value, replace = F)))
                y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                y[lowest_level] = 1 #set the randomly selected counts to 1
                y
                
              } #if all of the integer values of the lowest rank = 0, do random subsmpling
              else {
                sub_int <- subset(lowest_level_choise, (revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]] >= 1) == TRUE)
                sub_int_bind <- as.data.frame(cbind(sub_int,revtrunc_fixed_factor_1[,i][sub_int[,1]]))
                names(sub_int_bind)[1] <- "V1"
                names(sub_int_bind)[2] <- "V2"
                
                sub_int_bind_ordered <- sub_int_bind[order(sub_int_bind$V2, decreasing = TRUE),]  
                sub_int_bind_ordered
                
                sub_int_bind_ordered_V1 <- sub_int_bind_ordered$V1
                sub_int_bind_ordered_V2 <- sub_int_bind_ordered$V2
                
                if((length(unique(sub_int_bind_ordered_V2)) == 1 & length(sub_int_bind_ordered_V2)>as.vector(normalization_value))){
                  lowest_level <- as.numeric(as.vector(sample(as.factor(sub_int_bind_ordered_V1), normalization_value, replace = F)))
                  y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                  y[lowest_level] = 1 #set the randomly selected counts to 1
                  y
                } #if all of the integer values of the lowest rank are equal, do random subsampling
                else {
                  if(length(sub_int_bind_ordered_V1)>normalization_value){
                    
                    maxN_1 <- function(x, N=normalization_value){
                      len <- length(x)
                      if(N>len){
                        warning('N greater than length(x).  Setting N=length(x)')
                        N <- length(x)
                      }
                      sort(x,partial=len-N+1)[len-N+1]
                    }
                    max_1 <- which(as.data.frame(sub_int_bind_ordered_V2)[,1] == maxN_1(as.data.frame(sub_int_bind_ordered_V2)[,1]), arr.ind = TRUE)
                    max_1
                    sum(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1]))
                    
                    normalization_value_1 <- normalization_value - sum(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1])) # how may values are above the maxima that is the one to be normalized?
                    normalization_value_1
                    lowest_level_choise_1 <- as.data.frame(which(as.data.frame(sub_int_bind_ordered_V2)[,1] == unique(maxN_1(as.data.frame(sub_int_bind_ordered_V2)[,1]))))
                    lowest_level_choise_1
                    lowest_level <- as.numeric(as.vector(sample(as.factor(lowest_level_choise_1[,1]), normalization_value_1, replace = F)))
                    lowest_level <- sub_int_bind_ordered_V1[lowest_level]
                    lowest_level
                    
                    lowest_level_1 <- sub_int_bind_ordered_V1[(as.data.frame(sub_int_bind_ordered_V2)[,1] > unique(as.data.frame(sub_int_bind_ordered_V2)[,1][max_1]))]
                    
                    lowest_level <-c(lowest_level_1, lowest_level)
                    y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                    y[lowest_level] = 1 #set the randomly selected counts to 1
                    y
                    
                  } #if integer ranks are > normalization_value, do ranked subsampling
                  else {
                    if(length(sub_int_bind_ordered_V1)<normalization_value){
                      sub_int_zeros <- subset(lowest_level_choise, (revtrunc_fixed_factor_1[,i][lowest_level_choise[,1]] < 1) == TRUE)
                      length(t(sub_int_zeros))
                      
                      lowest_level_2 <- as.numeric(as.vector(sample(as.factor(sub_int_zeros[,1]), (normalization_value-length(sub_int_bind_ordered_V1)), replace = F)))
                      lowest_level_2
                      lowest_level_3 <- c(sub_int_bind_ordered_V1,lowest_level_2)
                      
                      y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                      y[lowest_level_3] = 1 #set the randomly selected counts to 1
                      y
                      
                    } #if integer ranks are < normalization_value, do ranked subsampling of the zero integers
                    else {
                      y <- as.vector(rep(0,length(x[,1]))) #create an empty vector
                      y[sub_int_bind_ordered_V1] = 1 #set the randomly selected counts to 1
                      y
                    } # if integer ranks are = normalization_value, sample all
                  }
                  
                  
                }
                
              }    
              
              SRS <- revtrunc_fixed_factor_1[,i] + ceiling(x[,i] > unique(x[,i][max])) + y #sum it all up
              SRS
              sum(SRS) #verification
              assign(paste(names(data)[i],sep=""),SRS)
              fixed_factor_1 <- cbind(fixed_factor_1, SRS)
              colnames(fixed_factor_1)[{
                counter = counter + 1
              }] <- names(data)[i]
            } #if the sum of the counts in the library > Cmin
          } #for all other libaries 
        }
        
        SRS_output <- fixed_factor_1
        SRS_output
      }
    }
  }
}
#---------------running the SRS function on example data####
SRS_output <- SRS(data = OTU_table, Cmin = Cmin)
SRS_output
dim(SRS_output)
#---------------confirmation of equal library size of the example data####
for(i in seq(1,ncol(SRS_output), 1)){
  if(sum(SRS_output[,i]) == Cmin){print(TRUE)} else {print(ERROR)}
}

#Make the normalized table
Normalized_table <- cbind(SRS_output,taxon)
#Add index column
Normalized_table <- cbind(count_table[,1],Normalized_table)
#Name the index column
colnames(Normalized_table)[1] <- "#OTU ID"
#Identify the rows(Species) that have no count across all samples after the normalization
non_zero <- rowSums(Normalized_table[,2:35]) != 0
#Eliminate the rows(Species) that have no count across all samples after the normalization
Normalized_new_table <- Normalized_table[non_zero,]
#Save the normalized table
write.csv(Normalized_new_table,"~/Microbiome_analysis/Normalized_data/0.1filter_Normalized_table.csv",row.names=FALSE)
