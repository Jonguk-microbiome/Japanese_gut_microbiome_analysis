# BMC Microbiology R script

# package list
library(ggplot2)
library(reshape2)
library(superheat)  
library(phyloseq)
library(dplyr)
library(vegan)
library(car)
library(S4Vectors)

# Import biom (NIBIOHN cohort or MORINAGA cohort)
biom <- import_biom("PATH/table_even10000.biom", treefilename = "PATH/rep_set.tre") # MORINAGA cohort : table_even5000.biom

# Read phenotype metadata file (map.csv)
map <- read.csv("PATH/map.tsv", sep = "\t", row.names = 1)

# Set numeric or factor
for (k in num:num) { 
  map[,k] <- as.factor(map[,k])
}

for (k in num:num) { 
  map[,k] <- as.numeric(map[,k])
}


# Add information from metadata to biom file
sample_data(biom)<- map

# OTU table
colnames(tax_table(biom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
OTU_table <- as.matrix(biom@otu_table@.Data)

# Taxonomy genus levels table with scientific classification
genus <- tax_glom(nexis_biom, taxrank = "Genus")

rownames(genus@otu_table) <- genus@tax_table[,6]
genus@tax_table@.Data <- gsub("_",":",genus@tax_table@.Data)
genus_tax_table <- as.data.frame(genus@tax_table@.Data)

genus_tax_table_name <- tidyr::unite_(genus_tax_table, paste(colnames(genus_tax_table)[-1], collapse=""), colnames(genus_tax_table)[-1])
rownames(genus@otu_table@.Data) <- genus_tax_table_name[,2]

genus_table <- as.data.frame(genus@otu_table@.Data)
# read numbers > percentage
genus_table <- genus_table/100 # NIBIOHN: genus_table/100, MORINAGA: genus_table/50

# Dominant genus
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}
dominant_genus_table=noise.removal(t(genus_table), percent=1)

# Estimate alpha diversity
biom_alpha_div <- estimate_richness(biom, measures = c("Observed", "Shannon", "Simpson", "Fisher", "Chao1"))



# Beta diversity analysis by Bray-Curtis
genus_bray <- vegdist(genus_table, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)

# PCoA
genus_bray_pco <- dudi.pco(genus_bray, scann = FALSE)      

# envfit analysis using phenotype metadata
genus_bray_envfit <- envfit(genus_bray_pco$tab, biom@sam_data, na.rm = TRUE, permutations = 10000) 

genus_bray_envfit.sp.df <- as.data.frame(genus_bray_envfit$vectors$arrows*sqrt(genus_bray_envfit$vectors$r))
genus_bray_envfit.sp.df$parameter <- rownames(genus_bray_envfit.sp.df)
genus_bray_envfit.sp.df$pval <- genus_bray_envfit$vectors$pvals
genus_bray_envfit.sp.df.select <- extractROWS(x = genus_bray_envfit.sp.df, i = genus_bray_envfit.sp.df$pval <= 0.001)
#pickup vector sum >> arrange using vecsum
genus_bray_envfit.sp.df.select$vecsum <- sqrt(genus_bray_envfit.sp.df.select$A1^2+genus_bray_envfit.sp.df.select$A2^2)
genus_bray_envfit.sp.df.select <- genus_bray_envfit.sp.df.select[order(-genus_bray_envfit.sp.df.select$vecsum),]

# eigen value for PCoA [Bray-Curtis]
genus_bray_pco_eig1per<-genus_bray_pco$eig[1]/sum(genus_bray_pco$eig)*100
genus_bray_pco_eig2per<-genus_bray_pco$eig[2]/sum(genus_bray_pco$eig)*100

# PCoA plot [Figure 2b]
ggplot(biom@sam_data, aes(x = genus_bray_pco$li[,1], y = genus_bray_pco$li[,2]))+
  geom_point(aes(size=2)+coord_fixed()+theme_bw()+
               #circle variance display
               stat_ellipse(aes(x=genus_bray_pco$li[,1],y=genus_bray_pco$li[,2]),level = 0.5)+
               labs(x=paste("PCo 1 (", format(round(genus_bray_pco_eig1per, 2)), "%)", sep = ""), 
                    y=paste("PCo 2 (", format(round(genus_bray_pco_eig2per, 2)), "%)", sep = ""))+
               geom_segment(data=genus_bray_envfit.sp.df.select,aes(x=0,xend=A1,y=0,yend=A2),   
                            arrow = arrow(length = unit(0.5, "cm")),colour="black")+
               geom_text_repel(data=genus_bray_envfit.sp.df.select,aes(x=A1,y=A2,label=parameter),size=5)


             

# envfit analysis and Cumulative effect size of non-redundant covariates selected by stepwise redundancy analysis
env_results_vector <- data.frame(genus_bray_envfit$vectors$r, genus_bray_envfit$vectors$pvals)
env_results_factor <- data.frame(genus_bray_envfit$factors$r, genus_bray_envfit$factors$pvals)
env_sig_list <- c(rownames(env_results_vector)[env_results_vector$genus_bray_envfit.vectors.pvals<0.05],
                  rownames(env_results_factor)[env_results_factor$genus_bray_envfit.factors.pvals<0.05])

# FDR 
env_results_vector$fdr <- round(p.adjust(env_results_vector$genus_bray_envfit.vectors.pvals, 'fdr'),3)
env_results_factor$fdr <- round(p.adjust(env_results_factor$genus_bray_envfit.factors.pvals, 'fdr'),3)
# write csv file for Figure 2a (displayed by excel) and Supplementary table3
write.csv(env_results_vector_NIBIOHN, file="PATH/.csv")
write.csv(env_results_factor_NIBIOHN, file="PATH/.csv")


pheno_sig <- biom@sam_data[,env_sig_list]
#NA omit
nrow(pheno_sig)
pheno_sig_comp<-na.omit(pheno_sig) 

genus_filtering <- extractROWS(x = genus_table, i = intersect(rownames(genus_table), rownames(pheno_sig_comp)))
nrow(genus_filtering)

res0.rda <- capscale(genus_filtering~1, data=pheno_sig_comp, distance="bray") 
res1.rda <- capscale(genus_filtering~., data=pheno_sig_comp, distance="bray")
result_stepwise <- ordiR2step(res0.rda, res1.rda, direction="both", R2scope = F)
result_stepwise_anova<-result_stepwise$anova
# Figure 2c (displayed by excel)
result_stepwise_anova             


# Heatmap dominant genus - metadata extracted by ordiR2step  [Figure 3]
result_anova_factors <- as.data.frame(result_stepwise$terminfo$ordered)
rownames(result_anova_factors)
pheno_sig_ordiR2step <- pheno_sig[,rownames(result_anova_factors)]

dominant_genus_table_t <- as.data.frame(t(dominant_genus_table))
colnames(dominant_genus_table_t) <- gsub("_NA","",colnames(dominant_genus_table_t))

colnames(pheno_sig_ordiR2step)  # check vector and factor

spearman_p<-cor(NIBIOHN_genus_denoized_t, 
                pheno_sig_ordiR2step[,-c(number)],  # remove factors
                use= "pairwise.complete.obs",
                method = "spearman")

col1 <- colorRampPalette(c("#00007F", "blue", "#007FFF","cyan", "white", "yellow", "#FF7F00", "red","#7F0000"))

superheat(spearman_p, 
          dist.method="manhattan",
          # place dendrograms on columns and rows 
          row.dendrogram = T, 
          col.dendrogram = T,
          bottom.label.col = "white",
          left.label.col = "white",
          # make gridlines white for enhanced prettiness
          grid.hline.col = "white",
          grid.vline.col = "white",
          left.label.size = 5,
          bottom.label.size = 0.5,
          left.label.text.size=5,
          bottom.label.text.size = 5,
          # rotate bottom label text
          bottom.label.text.angle = 90,
          X.text = round(as.matrix(spearman_p),2),
          X.text.col="white",
          X.text.size =4,
          bottom.label.text.alignment="right",
          left.label.text.alignment="right",
          scale = F)


#
#
# Beta diversity analysis by JSD for enterotype analysis
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
KLD <- function(x,y) sum(x * log(x/y))pheno_sig_morinaga_ordiR2step
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

genus_JSD=dist.JSD(t(genus_table))

# PCoA
genus_JSD_pco <- dudi.pco(genus_JSD, scann = FALSE)      

# CH index for pam clustering
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}
nclusters=NULL
for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(genus_JSD, k)
    nclusters[k]=index.G1(genus_table, data.cluster_temp,  d = genus_JSD,
                          centrotypes = "medoids")
  }
}
#[Figure 1b]
plot(nclusters, type="h", xlab="k clusters", ylab="CH index", main="pam clustering  [genus JSD]")

genus_JSD_pam.clu <- data.frame(pam.clustering(genus_JSD, k=3))  # set k number according to CH index 
biom@sam_data$enterotype <-as.factor(genus_JSD_pam.clu[,1]) 


# envfit analysis using genus table
genus_JSD_envfit <- envfit(genus_JSD_pco$tab, dominant_genus_table, na.rm = TRUE, permutations = 10000) 

genus_JSD_envfit.sp.df <- as.data.frame(genus_JSD_envfit$vectors$arrows*sqrt(genus_JSD_envfit$vectors$r))
genus_JSD_envfit.sp.df$parameter <- rownames(genus_JSD_envfit.sp.df)
genus_JSD_envfit.sp.df$pval <- genus_JSD_envfit$vectors$pvals
genus_JSD_envfit.sp.df.select <- extractROWS(x = genus_JSD_envfit.sp.df, i = genus_JSD_envfit.sp.df$pval <= 0.001)
#pickup vector sum >> arrange using vecsum
genus_JSD_envfit.sp.df.select$vecsum <- sqrt(genus_JSD_envfit.sp.df.select$A1^2+genus_JSD_envfit.sp.df.select$A2^2)
genus_JSD_envfit.sp.df.select <- genus_JSD_envfit.sp.df.select[order(-genus_JSD_envfit.sp.df.select$vecsum),]

genus_JSD_envfit.sp.df.select$parameter<-gsub(".*D.4","D.4",genus_JSD_envfit.sp.df.select$parameter)


# eigen value for PCoA [JSD]
genus_JSD_pco_eig1per<-genus_JSD_pco$eig[1]/sum(genus_JSD_pco$eig)*100
genus_JSD_pco_eig2per<-genus_JSD_pco$eig[2]/sum(genus_JSD_pco$eig)*100

# PCoA plot [Figure 1a,b]
ggplot(biom@sam_data, aes(x = genus_JSD_pco$li[,1], y = genus_JSD_pco$li[,2]))+
  geom_point(aes(colour=enterotype), size=2)+coord_fixed()+theme_bw()+
  #circle variance display
  stat_ellipse(aes(x=genus_JSD_pco$li[,1],y=genus_JSD_pco$li[,2], colour=enterotype),level = 0.5)+
  labs(x=paste("PCo 1 (", format(round(genus_JSD_pco_eig1per, 2)), "%)", sep = ""), 
       y=paste("PCo 2 (", format(round(genus_JSD_pco_eig2per, 2)), "%)", sep = ""))+
  geom_segment(data=genus_JSD_envfit.sp.df.select[1:3,],aes(x=0,xend=A1,y=0,yend=A2),   
               arrow = arrow(length = unit(0.5, "cm")),colour="black")+
  geom_text_repel(data=genus_JSD_envfit.sp.df.select[1:3,],aes(x=A1,y=A2,label=parameter),size=5)

# PCoA plot [Figure S2a,b]
ggplot(biom@sam_data, aes(x = genus_JSD_pco$li[,1], y = genus_JSD_pco$li[,2]))+
  geom_point(aes(colour=enterotype), size=2)+coord_fixed()+theme_bw()+
  #circle variance display
  stat_ellipse(aes(x=genus_JSD_pco$li[,1],y=genus_JSD_pco$li[,2], colour=enterotype),level = 0.5)+
  labs(x=paste("PCo 1 (", format(round(genus_JSD_pco_eig1per, 2)), "%)", sep = ""), 
       y=paste("PCo 2 (", format(round(genus_JSD_pco_eig2per, 2)), "%)", sep = ""))+
  geom_segment(data=genus_JSD_envfit.sp.df.select,aes(x=0,xend=A1,y=0,yend=A2),   
               arrow = arrow(length = unit(0.5, "cm")),colour="black")+
  geom_text_repel(data=genus_JSD_envfit.sp.df.select,aes(x=A1,y=A2,label=parameter),size=5)







# using Both cohorts (NIBIOHN cohort, MORINAGA cohort)

# CH index, PCoA [Both: NIBIOHN + MORINAGA]


genus_table <- genus_table/100 # NIBIOHN: genus_table/100, MORINAGA: genus_table/50

# Dominant genus
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}
dominant_genus_table=noise.removal(t(genus_table), percent=1)




NIBIOHN_genus_melt <- melt(as.matrix(genus_table))  # NIBIOHN cohort
morinaga_genus_melt <- melt(as.matrix(genus_table))  # MORINAGA cohort

Both_genus_melt <- rbind(NIBIOHN_genus_melt,morinaga_genus_melt)
head(Both_genus_melt)
dcast(Both_genus_melt)
Both_genus <- dcast(Both_genus_melt, Var1 ~ Var2)
head(Both_genus)
rownames(Both_genus) <- Both_genus$Var1
Both_genus <- Both_genus[,-1]
dim(Both_genus)
Both_genus[is.na(Both_genus)] <- 0     # NA > 0

Both_genus_JSD=dist.JSD(t(Both_genus))


## CH index
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}
nclusters=NULL
for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(Both_genus_JSD, k)
    nclusters[k]=index.G1(Both_genus, data.cluster_temp,  d = Both_genus_JSD,
                          centrotypes = "medoids")
  }
}
#[Figure S2d]
plot(nclusters, type="h", xlab="k clusters", ylab="CH index", main="pam clustering  [genus JSD]")

#PCoA================================================================================================================
Both_genus_JSD_pco<-dudi.pco(Both_genus_JSD, scann = FALSE)      # All dimension

#PCoA [JSD]
Both_genus_JSD_pco.eig1per<-Both_genus_JSD_pco$eig[1]/sum(Both_genus_JSD_pco$eig)*100
Both_genus_JSD_pco.eig2per<-Both_genus_JSD_pco$eig[2]/sum(Both_genus_JSD_pco$eig)*100
#=====================================================================================================================
#[enterotype_plot]====================================================================================================
Both_genus_JSD.pam.clu<-data.frame(pam.clustering(Both_genus_JSD, k=2))
enterotype_Both_genus_JSD.pam.clu <- as.factor(Both_genus_JSD.pam.clu[,1])



# Create Dominant list in both cohorts
NIBIOHN_MORINAGA_dominant_list <- c(rownames(dominant_genus_table),rownames(dominant_genus_table)) # NIBIOHN cohort and MORINAGA cohort
# remove list element(s)
NIBIOHN_MORINAGA_dominant_list <- NIBIOHN_MORINAGA_dominant_list[!duplicated(NIBIOHN_MORINAGA_dominant_list, fromLast = TRUE)]


# Make table by dominant list in each cohorts
NIBIOHN_dominant_genus_both <- genus_table[,NIBIOHN_MORINAGA_dominant_list]  # NIBIOHN cohort

morinaga_dominant_genus_both <- genus_table[,NIBIOHN_MORINAGA_dominant_list]  # MORINAGA cohort

NIBIOHN_dominant_genus_both$cohort <- c("NIBIOHN")
morinaga_dominant_genus_both$cohort <- c("MORINAGA")

NIBIOHN_dominant_genus_both_melt <- melt(as.matrix(NIBIOHN_dominant_genus_both))  
morinaga_dominant_genus_both_melt <- melt(as.matrix(morinaga_dominant_genus_both))  

dominant_genus_both_melt <- rbind(NIBIOHN_dominant_genus_both_melt, morinaga_dominant_genus_both_melt)
head(dominant_genus_both_melt)
dcast(dominant_genus_both_melt)
dominant_genus_both <- dcast(dominant_genus_both_melt, Var1 ~ Var2)

rownames(dominant_genus_both) <- dominant_genus_both$Var1
dominant_genus_both <- dominant_genus_both[,-1]
dim(dominant_genus_both)
dominant_genus_both[is.na(dominant_genus_both)] <- 0     #NA >>>>>>> 0

dominant_genus_both$enterotype <-enterotype_Both_genus_JSD.pam.clu

Both_envfit <- envfit(Both_genus_JSD_pco$tab, dominant_genus_both[,-c(24,25)], na.rm = TRUE, permutations = 10000) 

#[envfit plot for ggplot-taxon]
Both_envfit.sp.df <- as.data.frame(Both_envfit$vectors$arrows*sqrt(Both_envfit$vectors$r))
Both_envfit.sp.df$parameter <- rownames(Both_envfit.sp.df)
Both_envfit.sp.df$pval <- Both_envfit$vectors$pvals
Both_envfit.sp.df.select <- extractROWS(x = Both_envfit.sp.df, i = Both_envfit.sp.df$pval<=0.001)
#pickup vector sum >> arrange using vecsum
Both_envfit.sp.df.select$vecsum <- sqrt(Both_envfit.sp.df.select$A1^2+Both_envfit.sp.df.select$A2^2)
Both_envfit.sp.df.select <- Both_envfit.sp.df.select[order(-Both_envfit.sp.df.select$vecsum),]

Both_envfit.sp.df.select$parameter <- gsub(".*D.4","D.4",Both_envfit.sp.df.select$parameter)
Both_envfit.sp.df.select$parameter <- gsub("_NA","",Both_envfit.sp.df.select$parameter)


# PCoA both cohorts [Figure 1c]
ggplot(dominant_genus_both, aes(x = Both_genus_JSD_pco$li[,1], y = -Both_genus_JSD_pco$li[,2]))+
  geom_point(aes(colour=cohort), size=2)+coord_fixed()+#theme_bw()+
  #circle variance display
  stat_ellipse(aes(x=Both_genus_JSD_pco$li[,1],y=-Both_genus_JSD_pco$li[,2], colour=cohort),level = 0.5)+
  labs(title="Enterotype   [Genus JSD]", 
       x=paste("PCo 1 (", format(round(Both_genus_JSD_pco.eig1per, 2)), "%)", sep = ""), 
       y=paste("PCo 2 (", format(round(Both_genus_JSD_pco.eig2per, 2)), "%)", sep = ""))

# PCoA both cohorts [Figure S2c]
ggplot(dominant_genus_both, aes(x = Both_genus_JSD_pco$li[,1], y = -Both_genus_JSD_pco$li[,2]))+
  geom_point(aes(colour=enterotype, shape=cohort), size=2)+coord_fixed()+#theme_bw()+
  #circle variance display
  stat_ellipse(aes(x=Both_genus_JSD_pco$li[,1],y=-Both_genus_JSD_pco$li[,2], colour=enterotype),level = 0.5)+
  labs(title="Enterotype   [Genus JSD]", 
       x=paste("PCo 1 (", format(round(Both_genus_JSD_pco.eig1per, 2)), "%)", sep = ""), 
       y=paste("PCo 2 (", format(round(Both_genus_JSD_pco.eig2per, 2)), "%)", sep = ""))+
  geom_segment(data=Both_envfit.sp.df.select,aes(x=0,xend=-A1,y=0,yend=-A2),   
               arrow = arrow(length = unit(0.5, "cm")),colour="black")+
  geom_text_repel(data=Both_envfit.sp.df.select,aes(x=A1,y=-A2,label=parameter),size=5)



# Boxplot - dominant genus in both cohorts [figure 1e]
NIBIOHN_dominant_genus_both_melt <- melt(as.matrix(log10(NIBIOHN_dominant_genus_both)))  
morinaga_dominant_genus_both_melt <- melt(as.matrix(log10(morinaga_dominant_genus_both)))  

NIBIOHN_dominant_genus_both_melt$cohort <- c("NIBIOHN")
morinaga_dominant_genus_both_melt$cohort <- c("MORINAGA")

head(NIBIOHN_dominant_genus_both_melt)
head(morinaga_dominant_genus_both_melt)

dominant_genus_both_melt <- rbind(NIBIOHN_dominant_genus_both_melt, morinaga_dominant_genus_both_melt)
dominant_genus_both_melt$Var2 <- gsub("_NA", "", dominant_genus_both_melt$Var2)

ggplot(dominant_genus_both_melt, aes_string(x="value", y="Var2", fill="cohort")) + 
geom_boxplot(size = 0.5,width = 0.5)+theme_bw()


# Boxplot - alpha diversity in both cohorts [figure 1d]
# For alpha-diversity comparison analysis between NIBIOHN and MORINAGA cohorts, 5,000 reads were randomly reselected from NIBIOHN cohort.  
# Import biom (5,000 reads were randomly reselected from NIBIOHN cohort)
biom <- import_biom("PATH/table_even5000.biom", treefilename = "PATH/rep_set.tre")
# Read phenotype metadata file (map.csv)
map <- read.csv("PATH/map.tsv", sep = "\t", row.names = 1)

# Alpha diversity  [figure 1c]
NIBIOHN_biom_alpha_div <- estimate_richness(biom, measures = c("Observed", "Shannon", "Simpson", "Fisher", "Chao1"))  # NIBIOHN cohort
morinaga_biom_alpha_div <- estimate_richness(biom, measures = c("Observed", "Shannon", "Simpson", "Fisher", "Chao1"))  # MORINAGA cohort

NIBIOHN_biom_alpha_div$cohort <- c("NIBIOHN") 
morinaga_biom_alpha_div$cohort <- c("MORINAGA")

biom_alpha_div <- rbind(NIBIOHN_biom_alpha_div, morinaga_biom_alpha_div)

ggplot(biom_alpha_div, aes_string(x="cohort", y="Observed", fill="cohort")) +  # y= "Observed", "Shannon", "Simpson", "Fisher", "Chao1"
  geom_boxplot(size = 0.5,width = 0.5)+theme_bw()









