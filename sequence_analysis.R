
# Replication material for 'Stratified Pathways to Italy’s “Latest-Late” Transition to Adulthood'

# Date of the last update: 2023-07-26

# Title: Stratified Pathways to Italy’s “Latest-Late” Transition to Adulthood

# Author: Badolato Luca (badolato.3@osu.edu), Department of Sociology and Institute for Population Research, The Ohio State University

# Journal: Advances in Life Course Research

# Description: Script to run the Sequence Analysis and Clustering. Takes as input "data.dta" and returns as output "data_cluster_optimal_matching.dta"

# Clean the environment and set the working directory
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required R packages
listofpackages = c("tidyverse", "ggplot2","ggthemes", "ggExtra","readr","foreign", "TraMineR","readstata13", "cluster", "WeightedCluster","factoextra","NbClust", "haven")
for (j in listofpackages){
  if(sum(installed.packages()[, 1] == j) == 0) {
    install.packages(j)
  }
  library(j, character.only = T)
}

# Load the data, output of the stata file "data_cleaning.do"
data = read.dta13("data.dta", encoding=NULL )
head(data)  

# Define the color palette used in the graphs
colpal = c("#e5fbe5","palegreen2", "chartreuse2","green3", "lightskyblue1",
           "skyblue3","steelblue","royalblue4","khaki1","yellow2","gold1","darkgoldenrod1","salmon",
           "tomato","firebrick1","red3")

# Definition of the sequences, using the "states-sequence" (STS) format,
# in which the successive states of an individual are given in consecutive
# columns. 
ages = c("age13", "age14", "age15", "age16", "age17",
         "age18","age19","age20","age21","age22", 
         "age23","age24","age25","age26","age27",
         "age28","age29","age30")

seqstatl(data, var = ages, format='STS')

data_men = data[data$gender == "Male",]
data_women = data[data$gender == "Female",]
data.SPS <- seqformat(data, var = c("age13", "age14", "age15", "age16", "age17",
                                    "age18","age19","age20","age21","age22", 
                                    "age23","age24","age25","age26","age27",
                                    "age28","age29","age30"),
                      from = "STS", to = "SPS")

head(data.SPS)

# Definition of the alphabet, the list of possible states or events. Note
# that we have 16 possible states. Note that the state "0" correspond to 
# 0000 or that the state "101" to the 0101. 
data.alphabet <- c("0", "1", "10", "11", "100", "101", "110",
                   "111","1000","1001","1010","1011","1100","1101","1110"
                   ,"1111")
                   
data.labels <- c("0000", "0001", "0010", "0011", "0100", "0101", "0110",
                 "0111","1000","1001","1010","1011","1100","1101","1110"
                 ,"1111")
                 
data.scodes <- c("0000", "0001", "0010", "0011", "0100", "0101", "0110",
                 "0111","1000","1001","1010","1011","1100","1101","1110"
                 ,"1111")

# Definition of the sequences
data.seq <- seqdef(data,var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"),
                   alphabet = data.alphabet, states = data.scodes, 
                   labels = data.labels, xtstep = 6,
                     cpal = colpal,
                   weights = data$anweight)


summary(data.seq)
print(data.seq[1:10, 1:18], ext = TRUE)
print(data.seq[1:10, 1:18], ext = TRUE, format = "SPS")

data.seq_men <- seqdef(data_men, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"),
                   alphabet = data.alphabet, states = data.scodes, 
                   labels = data.labels, xtstep = 6,
                   cpal =  colpal,
                   weights = data_men$anweight)

data.seq_women <- seqdef(data_women, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"),
                   alphabet = data.alphabet, states = data.scodes, 
                   labels = data.labels, xtstep = 6,
                   cpal =  colpal,
                   weights = data_women$anweight)

# Legend plot
seqlegend(data.seq, ncol = 4, bty = "o")

# The index plot of the first 10 sequences 
seqiplot(data.seq, with.legend = FALSE, main = "Index plot (10 first sequences)", border = NA)

# The index plot:  
seqIplot(data.seq_men, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)
seqIplot(data.seq_women, sortv = "from.start", main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

# Sequence frequency plot
seqfplot(data.seq_men, with.legend = FALSE, main = "", ylab = "Cumulated frequency", border = NA)
seqfplot(data.seq_women, with.legend = FALSE, main = "", ylab = "Cumulated frequency", border = NA)

# Sequence frequency table
seqtab(data.seq_men)
seqtab(data.seq_women)

# State distribution plot
seqdplot(data.seq_men, with.legend = FALSE, main = "", ylab = "Frequency", border = NA)
seqstatd(data.seq_men)
seqdplot(data.seq_women, with.legend = FALSE, main = "", ylab = "Frequency", border = NA)
seqstatd(data.seq_women)

# Compute the ransition rates between states
tr <- seqtrate(data.seq)
round(tr, 4)

# Cluster analysis

# Compute the optimal matching distances using substitution costs based on transition rates
# observed in the data and a 1 indel cost. 
submat <- seqsubm(data.seq, method = "TRATE")
round(submat,4)

dist.om1 <- seqdist(data.seq, method = "OM", indel = 1, sm = submat)

# Cluster algorithm (ward) and dendrogram:
clusterward <- agnes(dist.om1, diss = TRUE, method = "ward")   

fviz_dend(clusterward, k = 7, 
          cex = 0.5, 
          k_colors = c("chartreuse1", "red", "orange", "blue","steelblue4",
                       "magenta", "coral"),
          rect = TRUE, 
          rect_border =  c("chartreuse1", "red", "orange", "blue","steelblue4",
                           "magenta",  "coral"),
          rect_fill = TRUE,
          ggtheme = theme_classic(),
          show_label = FALSE,
          lwd= 0.7,
         xlab= "Sequences")

# Cutting the dendrogram to define seven clusters:
clust <- cutree(clusterward, k = 7)
clustqual <- wcClusterQuality(dist.om1, clust)
clustqual$stats

# Save the cluster results to do the analysis in Stata
data_cluster = as.data.frame(clust)
write.dta(data_cluster, "data_cluster_optimal_matching.dta")

data.seq.clusters <- cbind(data, clust)

data.seq.clu1 <- data.seq.clusters[data.seq.clusters$clust==1,]
seqstatl(data.seq.clu1, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"), format='STS')
data.seq.clu1 <- seqdef(data.seq.clu1,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu1$anweight)

data.seq.clu2 <- data.seq.clusters[data.seq.clusters$clust==2,]
seqstatl(data.seq.clu2, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"), format='STS')
data.seq.clu2 <- seqdef(data.seq.clu2,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu2$anweight)

data.seq.clu3 <- data.seq.clusters[data.seq.clusters$clust==3,]
seqstatl(data.seq.clu3, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"), format='STS')
data.seq.clu3 <- seqdef(data.seq.clu3,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu3$anweight)

data.seq.clu4 <- data.seq.clusters[data.seq.clusters$clust==4,]
seqstatl(data.seq.clu4, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"), format='STS')
data.seq.clu4 <- seqdef(data.seq.clu4,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu4$anweight)

data.seq.clu5 <- data.seq.clusters[data.seq.clusters$clust==5,]
seqstatl(data.seq.clu5, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"), format='STS')
data.seq.clu5 <- seqdef(data.seq.clu5,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu5$anweight)

data.seq.clu6 <- data.seq.clusters[data.seq.clusters$clust==6,]
seqstatl(data.seq.clu6, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"), format='STS')
data.seq.clu6 <- seqdef(data.seq.clu6,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu6$anweight)

data.seq.clu7 <- data.seq.clusters[data.seq.clusters$clust==7,]
seqstatl(data.seq.clu7, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"), format='STS')
data.seq.clu7 <- seqdef(data.seq.clu7,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu7$anweight)

# State distribution plor for each cluster
seqdplot(data.seq.clu1, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu1, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

seqdplot(data.seq.clu2, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu2, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

seqdplot(data.seq.clu3, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu3, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

seqdplot(data.seq.clu4, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu4, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

seqdplot(data.seq.clu5, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu5, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

seqdplot(data.seq.clu6, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu6, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

seqdplot(data.seq.clu7, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu7, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

# Compute the modoid sequences

# Cluster1
dist.om_cluster1 <- seqdist(data.seq.clu1, method = "OM", indel = 1, sm = submat)
medoid_cluster1 <- seqrep(data.seq.clu1, diss = dist.om_cluster1, criterion = "dist",
                           nrep = 1)
print(medoid_cluster1, format = "SPS")
data.seq.clu1 = cbind(data.seq.clu1, id = c(1:dim(data.seq.clu1)[1]))

medoid_id_1    =  data.seq.clu1[data.seq.clu1$age13 == "0000" & data.seq.clu1$age14 == "0000" & data.seq.clu1$age15 == "0000" & 
                                  data.seq.clu1$age16 == "0000" & data.seq.clu1$age17 == "0000" & data.seq.clu1$age18 == "0000" & 
                                  data.seq.clu1$age19 == "1000" & data.seq.clu1$age20 == "1000" & data.seq.clu1$age21 == "1000" & 
                                  data.seq.clu1$age22 == "1000" & data.seq.clu1$age23 == "1000" & data.seq.clu1$age24 == "1000" & 
                                  data.seq.clu1$age25 == "1000" & data.seq.clu1$age26 == "1000" & data.seq.clu1$age27 == "1000" & 
                                  data.seq.clu1$age28 == "1000" & data.seq.clu1$age29 == "1110" & data.seq.clu1$age30 == "1110" , 19]

medoid_distance_1 = dist.om_cluster1[medoid_id_1[1], - medoid_id_1[1]]
mean_distance_1 = sum(medoid_distance_1) / (dim(data.seq.clu1)[1] - 1)
max_distance_1  = max(medoid_distance_1)

# Cluster2
dist.om_cluster2 <- seqdist(data.seq.clu2, method = "OM", indel = 1, sm = submat)
medoid_cluster2 <- seqrep(data.seq.clu2, diss = dist.om_cluster2, criterion = "dist",
                          nrep = 1)
print(medoid_cluster2, format = "SPS")
data.seq.clu2 = cbind(data.seq.clu2, id = c(1:dim(data.seq.clu2)[1]))

medoid_id_2    =  data.seq.clu2[data.seq.clu2$age13 == "0000" & data.seq.clu2$age14 == "0000" & data.seq.clu2$age15 == "0000" & 
                                  data.seq.clu2$age16 == "0000" & data.seq.clu2$age17 == "0000" & data.seq.clu2$age18 == "0000" & 
                                  data.seq.clu2$age19 == "0000" & data.seq.clu2$age20 == "0000"  & data.seq.clu2$age21 == "0000"  & 
                                  data.seq.clu2$age22 == "0000" & data.seq.clu2$age23 == "0000"  & data.seq.clu2$age24 == "0000"  & 
                                  data.seq.clu2$age25 == "0000"  & data.seq.clu2$age26 == "0000"  & data.seq.clu2$age27 == "0000"  & 
                                  data.seq.clu2$age28 == "0000"  & data.seq.clu2$age29 == "0000"  & data.seq.clu2$age30 == "0000"  , 19]

medoid_distance_2 = dist.om_cluster2[medoid_id_2[1], - medoid_id_2[1]]
mean_distance_2 = sum(medoid_distance_2) / (dim(data.seq.clu2)[1] - 1)
max_distance_2  = max(medoid_distance_2)

# Cluster3
dist.om_cluster3 <- seqdist(data.seq.clu3, method = "OM", indel = 1, sm = submat)
medoid_cluster3 <- seqrep(data.seq.clu3, diss = dist.om_cluster3, criterion = "dist",
                          nrep = 1)
print(medoid_cluster3, format = "SPS")
data.seq.clu3  = cbind(data.seq.clu3, id = c(1:dim(data.seq.clu3)[1]))

medoid_id_3    =  data.seq.clu3[data.seq.clu3$age13 == "0000" & data.seq.clu3$age14 == "0000" & data.seq.clu3$age15 == "0000" & 
                                  data.seq.clu3$age16 == "0000" & data.seq.clu3$age17 == "0000" & data.seq.clu3$age18 == "0000" & 
                                  data.seq.clu3$age19 == "0000" & data.seq.clu3$age20 == "0000" & data.seq.clu3$age21 == "0000" & 
                                  data.seq.clu3$age22 == "0000" & data.seq.clu3$age23 == "0110" & data.seq.clu3$age24 == "0110"  & 
                                  data.seq.clu3$age25 == "0111" & data.seq.clu3$age26 == "0111" & data.seq.clu3$age27 == "0111" & 
                                  data.seq.clu3$age28 == "0111" & data.seq.clu3$age29 == "0111" & data.seq.clu3$age30 == "0111" , 19]

medoid_distance_3 = dist.om_cluster3[medoid_id_3[1], - medoid_id_3[1]]
mean_distance_3 = sum(medoid_distance_3) / (dim(data.seq.clu3)[1] - 1)
max_distance_3  = max(medoid_distance_3)

# Cluster4
dist.om_cluster4 = seqdist(data.seq.clu4, method = "OM", indel = 1, sm = submat)
medoid_cluster4  = seqrep(data.seq.clu4, diss = dist.om_cluster4, criterion = "dist",
                          nrep = 1)
print(medoid_cluster4, format = "SPS")
data.seq.clu4 = cbind(data.seq.clu4, id = c(1:dim(data.seq.clu4)[1]))

medoid_id_4   =  data.seq.clu4[data.seq.clu4$age13 == "0000" & data.seq.clu4$age14 == "0000" & data.seq.clu4$age15 == "0000" & 
                                  data.seq.clu4$age16 == "0000" & data.seq.clu4$age17 == "0000" & data.seq.clu4$age18 == "0000" & 
                                  data.seq.clu4$age19 == "0000" & data.seq.clu4$age20 == "0000" & data.seq.clu4$age21 == "0000" & 
                                  data.seq.clu4$age22 == "0000" & data.seq.clu4$age23 == "0000" & data.seq.clu4$age24 == "1000" & 
                                  data.seq.clu4$age25 == "1000" & data.seq.clu4$age26 == "1000" & data.seq.clu4$age27 == "1000" & 
                                  data.seq.clu4$age28 == "1000" & data.seq.clu4$age29 == "1000" & data.seq.clu4$age30 == "1110" , 19]

medoid_distance_4 = dist.om_cluster4[medoid_id_4[1], - medoid_id_4[1]]
mean_distance_4 = sum(medoid_distance_4) / (dim(data.seq.clu4)[1] - 1)
max_distance_4  = max(medoid_distance_4)

# Cluster5
dist.om_cluster5 = seqdist(data.seq.clu5, method = "OM", indel = 1, sm = submat)
medoid_cluster5  = seqrep(data.seq.clu5, diss = dist.om_cluster5, criterion = "dist",
                          nrep = 1)
print(medoid_cluster5, format = "SPS")
data.seq.clu5 = cbind(data.seq.clu5, id = c(1:dim(data.seq.clu5)[1]))

medoid_id_5    =  data.seq.clu5[data.seq.clu5$age13   == "0000" & data.seq.clu5$age14 == "0000" & data.seq.clu5$age15 == "0000" & 
                                  data.seq.clu5$age16 == "0000" & data.seq.clu5$age17 == "0000" & data.seq.clu5$age18 == "0000" & 
                                  data.seq.clu5$age19 == "0000" & data.seq.clu5$age20 == "0000" & data.seq.clu5$age21 == "0000" & 
                                  data.seq.clu5$age22 == "0000" & data.seq.clu5$age23 == "0100" & data.seq.clu5$age24 == "0100" & 
                                  data.seq.clu5$age25 == "0100" & data.seq.clu5$age26 == "0100" & data.seq.clu5$age27 == "0100" & 
                                  data.seq.clu5$age28 == "0100" & data.seq.clu5$age29 == "0100" & data.seq.clu5$age30 == "0100" , 19]

medoid_distance_5 = dist.om_cluster5[medoid_id_5[1], - medoid_id_5[1]]
mean_distance_5 = sum(medoid_distance_5) / (dim(data.seq.clu5)[1] - 1)
max_distance_5  = max(medoid_distance_5)

# Cluster6
dist.om_cluster6 = seqdist(data.seq.clu6, method = "OM", indel = 1, sm = submat)
medoid_cluster6  = seqrep(data.seq.clu6, diss = dist.om_cluster6, criterion = "dist",
                          nrep = 1)
print(medoid_cluster6, format = "SPS")
data.seq.clu6  = cbind(data.seq.clu6, id = c(1:dim(data.seq.clu6)[1]))

medoid_id_6    =  data.seq.clu6[data.seq.clu6$age13 == "0000" & data.seq.clu6$age14 == "0000" & data.seq.clu6$age15 == "0000"  & 
                                  data.seq.clu6$age16 == "0000"  & data.seq.clu6$age17 == "0000"  & data.seq.clu6$age18 == "0000"  & 
                                  data.seq.clu6$age19 == "0000"  & data.seq.clu6$age20 == "1100" & data.seq.clu6$age21 == "1100" & 
                                  data.seq.clu6$age22 == "1100" & data.seq.clu6$age23 == "1100" & data.seq.clu6$age24 == "1100" & 
                                  data.seq.clu6$age25 == "1100" & data.seq.clu6$age26 == "1100" & data.seq.clu6$age27 == "1100" & 
                                  data.seq.clu6$age28 == "1110" & data.seq.clu6$age29 == "1110" & data.seq.clu6$age30 == "1110" , 19]

medoid_distance_6 = dist.om_cluster6[medoid_id_6[1], - medoid_id_6[1]]
mean_distance_6 = sum(medoid_distance_6) / (dim(data.seq.clu6)[1] - 1)
max_distance_6  = max(medoid_distance_6)

# Cluster7
dist.om_cluster7 = seqdist(data.seq.clu7, method = "OM", indel = 1, sm = submat)
medoid_cluster7  = seqrep(data.seq.clu7, diss = dist.om_cluster7, criterion = "dist",
                          nrep = 1)
print(medoid_cluster7, format = "SPS")
data.seq.clu7  = cbind(data.seq.clu7, id = c(1:dim(data.seq.clu7)[1]))

medoid_id_7    =  data.seq.clu7[data.seq.clu7$age13 == "0000" & data.seq.clu7$age14 == "0000" & data.seq.clu7$age15 == "0000" & 
                                  data.seq.clu7$age16 == "0000" & data.seq.clu7$age17 == "0000" & data.seq.clu7$age18 == "1000" & 
                                  data.seq.clu7$age19 == "1000" & data.seq.clu7$age20 == "1000" & data.seq.clu7$age21 == "1110" & 
                                  data.seq.clu7$age22 == "1111" & data.seq.clu7$age23 == "1111" & data.seq.clu7$age24 == "1111" & 
                                  data.seq.clu7$age25 == "1111" & data.seq.clu7$age26 == "1111" & data.seq.clu7$age27 == "1111" & 
                                  data.seq.clu7$age28 == "1111" & data.seq.clu7$age29 == "1111" & data.seq.clu7$age30 == "1111" , 19]

medoid_distance_7 = dist.om_cluster7[medoid_id_7[1], - medoid_id_7[1]]
mean_distance_7 = sum(medoid_distance_7) / (dim(data.seq.clu7)[1] - 1)
max_distance_7  = max(medoid_distance_7)

# Robustness check using the Hamming distance
dist.om1 <- seqdist(data.seq, method = "HAM", indel = 1, sm = submat)

# Cluster algorithm (ward) and dendrogram:
clusterward <- agnes(dist.om1, diss = TRUE, method = "ward")   

# Cutting the dendrogram to define seven clusters:
clust <- cutree(clusterward, k = 7)
clustqual <- wcClusterQuality(dist.om1, clust)
clustqual$stats

data.seq.clusters <- cbind(data, clust)

data.seq.clu1 <- data.seq.clusters[data.seq.clusters$clust==1,]
seqstatl(data.seq.clu1, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"), format='STS')
data.seq.clu1 <- seqdef(data.seq.clu1,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu1$anweight)

data.seq.clu2 <- data.seq.clusters[data.seq.clusters$clust==2,]
seqstatl(data.seq.clu2, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"), format='STS')
data.seq.clu2 <- seqdef(data.seq.clu2,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu2$anweight)

data.seq.clu3 <- data.seq.clusters[data.seq.clusters$clust==3,]
seqstatl(data.seq.clu3, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"), format='STS')
data.seq.clu3 <- seqdef(data.seq.clu3,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu3$anweight)

data.seq.clu4 <- data.seq.clusters[data.seq.clusters$clust==4,]
seqstatl(data.seq.clu4, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"), format='STS')
data.seq.clu4 <- seqdef(data.seq.clu4,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu4$anweight)

data.seq.clu5 <- data.seq.clusters[data.seq.clusters$clust==5,]
seqstatl(data.seq.clu5, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"), format='STS')
data.seq.clu5 <- seqdef(data.seq.clu5,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu5$anweight)

data.seq.clu6 <- data.seq.clusters[data.seq.clusters$clust==6,]
seqstatl(data.seq.clu6, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"), format='STS')
data.seq.clu6 <- seqdef(data.seq.clu6,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu6$anweight)

data.seq.clu7 <- data.seq.clusters[data.seq.clusters$clust==7,]
seqstatl(data.seq.clu7, var = c("age13", "age14", "age15", "age16", "age17",
                                "age18","age19","age20","age21","age22", 
                                "age23","age24","age25","age26","age27",
                                "age28","age29","age30"), format='STS')
data.seq.clu7 <- seqdef(data.seq.clu7,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu7$anweight)

# State distribution plor for each cluster
seqdplot(data.seq.clu1, border = NA, with.legend = FALSE)
seqdplot(data.seq.clu2, border = NA, with.legend = FALSE)
seqdplot(data.seq.clu3, border = NA, with.legend = FALSE)
seqdplot(data.seq.clu4, border = NA, with.legend = FALSE)
seqdplot(data.seq.clu5, border = NA, with.legend = FALSE)
seqdplot(data.seq.clu6, border = NA, with.legend = FALSE)
seqdplot(data.seq.clu7, border = NA, with.legend = FALSE)

# Robustness check: gender-based typology

# Men
dist.om_men <- seqdist(data.seq_men, method = "OM", indel = 1, sm = submat)

# Cluster algorithm (ward) and dendrogram:
clusterward_men <- agnes(dist.om_men, diss = TRUE, method = "ward")   

fviz_dend(clusterward_men, k = 5, 
          cex = 0.5, 
          k_colors = c("chartreuse1", "red", "orange", "blue","steelblue4"),
          rect = TRUE, 
          rect_border =  c("chartreuse1", "red", "orange", "blue","steelblue4"),
          rect_fill = TRUE,
          ggtheme = theme_classic(),
          show_label = FALSE,
          lwd= 0.7,
          xlab= "Sequences")

# Cutting the dendrogram to define seven clusters:
clust <- cutree(clusterward_men, k = 5)

data.seq.clusters_men <- cbind(data_men, clust)

data.seq.clu1_men <- data.seq.clusters_men[data.seq.clusters_men$clust==1,]
data.seq.clu1_men <- seqdef(data.seq.clu1_men,var = c("age13", "age14", "age15", "age16", "age17",
                                              "age18","age19","age20","age21","age22", 
                                              "age23","age24","age25","age26","age27",
                                              "age28","age29","age30"),
                        alphabet = data.alphabet, states = data.scodes, 
                        labels = data.labels, xtstep = 6,
                        cpal = colpal,
                        weights = data.seq.clu1_men$anweight)

data.seq.clu2_men <- data.seq.clusters_men[data.seq.clusters_men$clust==2,]
data.seq.clu2_men <- seqdef(data.seq.clu2_men,var = c("age13", "age14", "age15", "age16", "age17",
                                                          "age18","age19","age20","age21","age22", 
                                                          "age23","age24","age25","age26","age27",
                                                          "age28","age29","age30"),
                              alphabet = data.alphabet, states = data.scodes, 
                              labels = data.labels, xtstep = 6,
                              cpal = colpal,
                              weights = data.seq.clu1_men$anweight)

data.seq.clu3_men <- data.seq.clusters_men[data.seq.clusters_men$clust==3,]
data.seq.clu3_men <- seqdef(data.seq.clu3_men,var = c("age13", "age14", "age15", "age16", "age17",
                                                          "age18","age19","age20","age21","age22", 
                                                          "age23","age24","age25","age26","age27",
                                                          "age28","age29","age30"),
                              alphabet = data.alphabet, states = data.scodes, 
                              labels = data.labels, xtstep = 6,
                              cpal = colpal,
                              weights = data.seq.clu1_men$anweight)

data.seq.clu4_men <- data.seq.clusters_men[data.seq.clusters_men$clust==4,]
data.seq.clu4_men <- seqdef(data.seq.clu4_men,var = c("age13", "age14", "age15", "age16", "age17",
                                                          "age18","age19","age20","age21","age22", 
                                                          "age23","age24","age25","age26","age27",
                                                          "age28","age29","age30"),
                              alphabet = data.alphabet, states = data.scodes, 
                              labels = data.labels, xtstep = 6,
                              cpal = colpal,
                              weights = data.seq.clu1_men$anweight)

data.seq.clu5_men <- data.seq.clusters_men[data.seq.clusters_men$clust==5,]
data.seq.clu5_men <- seqdef(data.seq.clu5_men,var = c("age13", "age14", "age15", "age16", "age17",
                                                          "age18","age19","age20","age21","age22", 
                                                          "age23","age24","age25","age26","age27",
                                                          "age28","age29","age30"),
                              alphabet = data.alphabet, states = data.scodes, 
                              labels = data.labels, xtstep = 6,
                              cpal = colpal,
                              weights = data.seq.clu1_men$anweight)

# State distribution plor for each cluster
seqdplot(data.seq.clu1_men, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu1_men, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

seqdplot(data.seq.clu2_men, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu2_men, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

seqdplot(data.seq.clu3_men, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu3_men, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

seqdplot(data.seq.clu4_men, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu4_men, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

seqdplot(data.seq.clu5_men, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu5_men, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

# Compute the modoid sequences for men-based trajectories
dist.om_cluster1_men <- seqdist(data.seq.clu1_men, method = "OM", indel = 1, sm = submat)
medoid_cluster1 <- seqrep(data.seq.clu1_men, diss = dist.om_cluster1_men, criterion = "dist",
                          nrep = 1)
print(medoid_cluster1, format = "SPS")

dist.om_cluster2_men <- seqdist(data.seq.clu2_men, method = "OM", indel = 1, sm = submat)
medoid_cluster2 <- seqrep(data.seq.clu2_men, diss = dist.om_cluster2_men, criterion = "dist",
                          nrep = 1)
print(medoid_cluster2, format = "SPS")

dist.om_cluster3_men <- seqdist(data.seq.clu3_men, method = "OM", indel = 1, sm = submat)
medoid_cluster3 <- seqrep(data.seq.clu3_men, diss = dist.om_cluster3_men, criterion = "dist",
                          nrep = 1)
print(medoid_cluster3, format = "SPS")

dist.om_cluster4_men <- seqdist(data.seq.clu4_men, method = "OM", indel = 1, sm = submat)
medoid_cluster4 <- seqrep(data.seq.clu4_men, diss = dist.om_cluster4_men, criterion = "dist",
                          nrep = 1)
print(medoid_cluster4, format = "SPS")

dist.om_cluster5_men <- seqdist(data.seq.clu5_men, method = "OM", indel = 1, sm = submat)
medoid_cluster5 <- seqrep(data.seq.clu5_men, diss = dist.om_cluster5_men, criterion = "dist",
                          nrep = 1)
print(medoid_cluster5, format = "SPS")

# Women
dist.om_women <- seqdist(data.seq_women, method = "OM", indel = 1, sm = submat)

# Cluster algorithm (ward) and dendrogram:
clusterward_women <- agnes(dist.om_women, diss = TRUE, method = "ward")   

fviz_dend(clusterward_women, k = 5, 
          cex = 0.5, 
          k_colors = c("chartreuse1", "red", "orange", "blue","steelblue4"),
          rect = TRUE, 
          rect_border =  c("chartreuse1", "red", "orange", "blue","steelblue4"),
          rect_fill = TRUE,
          ggtheme = theme_classic(),
          show_label = FALSE,
          lwd= 0.7,
          xlab= "Sequences")

# Cutting the dendrogram to define seven clusters:
clust <- cutree(clusterward_women, k = 5)

data.seq.clusters_women <- cbind(data_women, clust)

data.seq.clu1_women <- data.seq.clusters_women[data.seq.clusters_women$clust==1,]
data.seq.clu1_women <- seqdef(data.seq.clu1_women,var = c("age13", "age14", "age15", "age16", "age17",
                                                          "age18","age19","age20","age21","age22", 
                                                          "age23","age24","age25","age26","age27",
                                                          "age28","age29","age30"),
                              alphabet = data.alphabet, states = data.scodes, 
                              labels = data.labels, xtstep = 6,
                              cpal = colpal,
                              weights = data.seq.clu1_women$anweight)

data.seq.clu2_women <- data.seq.clusters_women[data.seq.clusters_women$clust==2,]
data.seq.clu2_women <- seqdef(data.seq.clu2_women,var = c("age13", "age14", "age15", "age16", "age17",
                                                          "age18","age19","age20","age21","age22", 
                                                          "age23","age24","age25","age26","age27",
                                                          "age28","age29","age30"),
                              alphabet = data.alphabet, states = data.scodes, 
                              labels = data.labels, xtstep = 6,
                              cpal = colpal,
                              weights = data.seq.clu1_women$anweight)

data.seq.clu3_women <- data.seq.clusters_women[data.seq.clusters_women$clust==3,]
data.seq.clu3_women <- seqdef(data.seq.clu3_women,var = c("age13", "age14", "age15", "age16", "age17",
                                                          "age18","age19","age20","age21","age22", 
                                                          "age23","age24","age25","age26","age27",
                                                          "age28","age29","age30"),
                              alphabet = data.alphabet, states = data.scodes, 
                              labels = data.labels, xtstep = 6,
                              cpal = colpal,
                              weights = data.seq.clu1_women$anweight)

data.seq.clu4_women <- data.seq.clusters_women[data.seq.clusters_women$clust==4,]
data.seq.clu4_women <- seqdef(data.seq.clu4_women,var = c("age13", "age14", "age15", "age16", "age17",
                                                          "age18","age19","age20","age21","age22", 
                                                          "age23","age24","age25","age26","age27",
                                                          "age28","age29","age30"),
                              alphabet = data.alphabet, states = data.scodes, 
                              labels = data.labels, xtstep = 6,
                              cpal = colpal,
                              weights = data.seq.clu1_women$anweight)

data.seq.clu5_women <- data.seq.clusters_women[data.seq.clusters_women$clust==5,]
data.seq.clu5_women <- seqdef(data.seq.clu5_women,var = c("age13", "age14", "age15", "age16", "age17",
                                                          "age18","age19","age20","age21","age22", 
                                                          "age23","age24","age25","age26","age27",
                                                          "age28","age29","age30"),
                              alphabet = data.alphabet, states = data.scodes, 
                              labels = data.labels, xtstep = 6,
                              cpal = colpal,
                              weights = data.seq.clu1_women$anweight)

# State distribution plor for each cluster
seqdplot(data.seq.clu1_women, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu1_women, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

seqdplot(data.seq.clu2_women, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu2_women, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

seqdplot(data.seq.clu3_women, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu3_women, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

seqdplot(data.seq.clu4_women, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu4_women, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

seqdplot(data.seq.clu5_women, border = NA, with.legend = FALSE)
seqIplot(data.seq.clu5_women, sortv = "from.start",main = "", xlab = "Age", ylab = "Sequences", with.legend = FALSE)

# Compute the medoid sequences for women-based trajectories
dist.om_cluster1_women <- seqdist(data.seq.clu1_women, method = "OM", indel = 1, sm = submat)
medoid_cluster1 <- seqrep(data.seq.clu1_women, diss = dist.om_cluster1_women, criterion = "dist",
                          nrep = 1)
print(medoid_cluster1, format = "SPS")

dist.om_cluster2_women <- seqdist(data.seq.clu2_women, method = "OM", indel = 1, sm = submat)
medoid_cluster2 <- seqrep(data.seq.clu2_women, diss = dist.om_cluster2_women, criterion = "dist",
                          nrep = 1)
print(medoid_cluster2, format = "SPS")

dist.om_cluster3_women <- seqdist(data.seq.clu3_women, method = "OM", indel = 1, sm = submat)
medoid_cluster3 <- seqrep(data.seq.clu3_women, diss = dist.om_cluster3_women, criterion = "dist",
                          nrep = 1)
print(medoid_cluster3, format = "SPS")

dist.om_cluster4_women <- seqdist(data.seq.clu4_women, method = "OM", indel = 1, sm = submat)
medoid_cluster4 <- seqrep(data.seq.clu4_women, diss = dist.om_cluster4_women, criterion = "dist",
                          nrep = 1)
print(medoid_cluster4, format = "SPS")

dist.om_cluster5_women <- seqdist(data.seq.clu5_women, method = "OM", indel = 1, sm = submat)
medoid_cluster5 <- seqrep(data.seq.clu5_women, diss = dist.om_cluster5_women, criterion = "dist",
                          nrep = 1)
print(medoid_cluster5, format = "SPS")

# Save the gender-specific clusters
data_cluster_robust <- as.data.frame(rbind(data.seq.clusters_women, data.seq.clusters_men))
write_dta(data_cluster_robust, "data_cluster_optimal_matching_gender_specific.dta")

# End
