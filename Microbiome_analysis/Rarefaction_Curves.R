
#Read the count table
Count_table <- read.csv("~/Microbiome_analysis/OTU_data/count_table.csv")

#First arrange the table required for rarefaction curves plotting
#Extract taxon names
taxo <- Count_table[38]
#Set it as a data frame
taxo <- as.data.frame(taxo)
#Name the column
colnames(taxo)<- "taxo"
#Extract the indexies
ID <- Count_table[1]
#Set it as a data frame
ID <- as.data.frame(ID)
#Name the column
colnames(ID) <- "name"

#Make a new data frame 
df <- data.frame(matrix(ncol = 4, nrow = 0))
#Set the column names
colnames(df) <-c("name","taxo","Group","value")
#Concatenate information (indices, taxon names, names of dogs, the count of each species) in  the data framed
for (i in 2:37){
  g <- rep(colnames(Count_table[i]), time=3643)
  g <- as.data.frame(g)
  colnames(g) <- "Group"
  col <- Count_table[i]
  col <- as.data.frame(col)
  colnames(col) <- "Value"
  name <- cbind(ID,taxo,g,col)
  non_zero <- name$Value != 0
  name <- name[non_zero,]
  df <- rbind(df,name)
}


#Here we will plot the curves using the table we have just prepared
library(tidyverse)
library(tidyr)
library(magrittr) 
library(dplyr) 


collect <- function(data, group){
  data %>%
    filter(Group == group) %>%
    uncount(Value) %>%
    sample_n(n()) %>%
    mutate(observation = row_number()) %>%
    arrange(name, observation) %>%
    group_by(name) %>%
    mutate(distinct = row_number() == 1) %>%
    ungroup() %>%
    arrange(observation) %>%
    mutate(s = cumsum(distinct)) %>%
    select(observation, s) 
}

#Make the list of corours for 18 dogs
colors <- c("Red","Orange","Yellow","GreenYellow","Green","olivedrab","LightSeaGreen","LightSkyBlue","DodgerBlue","Blue",
            "cyan","Aquamarine","MediumPurple","Purple","HotPink","DarkGray","Sienna","Tan")


#Rarefaction curves for the previous samples
#Perform the collecting curves 100 times for "Reggie"
data <- map_dfr(1:100, ~collect(df,"Pre_Reggie"),.id = "iteration")  
#Plot the curve of the average of 100 collecting curves
#"Reggie" was plotted first as it has the largest scale across samples
data2 <- data %>%
  group_by(observation) %>%
  summarize(r=mean(s))
plot(data2$observation, data2$r, type = "l",ylab = "Number of obserbation", xlab = "Number of sequences",col=colors[17])

#Make the list of samples except "Reggie"
list <- 2:17
list <- append(list,19)
#Perform the collecting curves 100 times for the rest of samples and plot the average 
for (i in list){
  color <- colors[i-1]
  data <- map_dfr(1:100, ~collect(df,colnames(Count_table[i])), .id = "iteration") 
  data2 <- data %>%
    group_by(observation) %>%
    summarize(r=mean(s))
  lines(data2$observation, data2$r, type="l", col=color)
}

#Extract names of samples
name <- colnames(Count_table)[2:19]
#Add the legend
legend("bottomright", legend = name,
       lwd = 1.5, col = colors,bty="n",ncol=2,cex=0.6,pt.cex=0.7)  
title(main = "Rarefaction curves")



#Rarefaction curves for the post samples
#Make the list of names for the post samples
list_name <- colnames(Count_table)[20:37]

#Perform the collecting curves 100 times for "Digby"
data <- map_dfr(1:100, ~collect(df,"Post_Digby"),.id = "iteration") 
#Plot the curve of the average of 100 collecting curves
#"Digby" was plotted first as it has the largest scale across samples
data2 <- data %>%
  group_by(observation) %>%
  summarize(r=mean(s))
plot(data2$observation, data2$r, type = "l",ylab = "Number of obserbation", xlab = "Number of sequences",col=colors[9])

#Make the list of samples except "Digby"
list2 <- 1:8
list2 <- append(list2,10:18)
#Perform the collecting curves 100 times for the rest of samples and plot the average
for (i in list2){
  color <- colors[i]
  data <- map_dfr(1:100, ~collect(df,colnames(list_name[i])),.id = "iteration") 
  data2 <- data %>%
    group_by(observation) %>%
    summarize(r=mean(s))
  lines(data2$observation, data2$r, type="l", col=color)
}

#Add the legend
legend("bottomright", legend = list_name,
       lwd = 1.5, col = colors,bty="n",ncol=2,cex=0.6,pt.cex=0.7)  
title(main = "Rarefaction curves")



#Rarefaction curves for the post samples (Zoomed in)
#Perform the collecting curves 100 times for "Bob"
data <- map_dfr(1:100, ~collect(df,"Post_Bob"),.id = "iteration")  
#Plot the curve of the average of 100 collecting curves
#"Bob" was plotted first as it has the similar scale to the previous group and represent the post samples
data2 <- data %>%
  group_by(observation) %>%
  summarize(r=mean(s))
plot(data2$observation, data2$r, type = "l",ylab = "Number of obserbation", xlab = "Number of sequences",col=colors[4])

#Make the list of samples except "Bob"
list3 <- 1:3
list3 <- append(list3,5:18)
#Perform the collecting curves 100 times for the rest of samples and plot the average
for (i in list3){
  color <- colors[i]
  data <- map_dfr(1:100, ~collect(df,list_name[i]),.id = "iteration") 
  data2 <- data %>%
    group_by(observation) %>%
    summarize(r=mean(s))
  lines(data2$observation, data2$r, type="l", col=color)
}

#Add the legend
legend("bottomright", legend = list_name,
       lwd = 1.5, col = colors,bty="n",ncol=2,cex=0.6,pt.cex=0.7)  
title(main = "Rarefaction curves")
