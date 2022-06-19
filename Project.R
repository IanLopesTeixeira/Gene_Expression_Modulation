setwd("/Users/diogo/Desktop/Gene Exp. Modulation")

library(data.table)

bt <- fread("result.txt", nrows = 3)
ht <- fread("head.txt", nrows = 3)

names(bt) <- colnames(ht)

s <- t(bt)
s <- s[-c(1, 2),]

tissue <- fread("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

x <- tissue[match(names(s),unlist(tissue[,1])),"SMTSD"]
metadata <- cbind(names(s),s,x)
colnames(metadata) <- c("Sample ID", "exp", "Tissue")

age <- fread("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
x <- age[match(names(s),unlist(tissue[,1])),"SMTSD"]

new <- sapply(unlist(metadata[,"Sample ID"]), function(x) paste(unlist(strsplit(x,"-"))[1:2],collapse = "-"))
metadata <- cbind(metadata,age[match(new,unlist(age[,1])),c("SEX","AGE")])
save(metadata,file = "Project.RData")
####################################################################### 
########################## CORRER A PARTIR DAQUI ######################
#######################################################################
load("Project.RData")
View(metadata)


####################################################################### 
#######################################################################
#######################################################################
tissue <- unique(metadata$Tissue) #vetor com os tecidos 

length(which(metadata$AGE[grep("Brain",metadata$Tissue)] == "30-39"))
#Para saber o n. de 

HIF1A <- fread("result_HIF1A.txt")
HIF1A <- t(HIF1A)
HIF1A <- HIF1A[-c(1, 2),]
STAT1 <- fread("result_STAT1.txt")
STAT1 <- t(STAT1)
STAT1 <- STAT1[-c(1, 2),]
TAF1 <- fread("result_TAF1.txt")
TAF1 <- t(TAF1)
TAF1 <- TAF1[-c(1, 2),]
TRIM22 <- fread("result_TRIM22.txt")
TRIM22 <- t(TRIM22)
TRIM22 <- TRIM22[-c(1, 2),]
GABPA <- fread("result_GABPA.txt")
GABPA <- t(GABPA)
GABPA <- GABPA[-c(1, 2),]
SP1 <- fread("result_SP1.txt")
SP1 <- t(SP1)
SP1 <- SP1[-c(1, 2),]
SIN3A <- fread("result_SIN3A.txt")
SIN3A<- t(SIN3A)
SIN3A <- SIN3A[-c(1, 2),]

expression_table <- rbind("SETD2" = metadata$exp, "HIF1A" = HIF1A, "STAT1" = STAT1, "TAF1" = TAF1, "TRIM22" = TRIM22, "GABPA" = GABPA, "SP1" = SP1, "SIN3A" = SIN3A)
colnames(expression_table) <- metadata$`Sample ID`
save(expression_table,file = "expression_table.RData")
####################################################################### 
#######################################################################
#######################################################################

load("expression_table.RData")
View(expression_table)

cor_table <- matrix(NA, nrow = 7, ncol = 3, dimnames = list(rownames(expression_table)[-1], c("Correlation","P_value","P_value_adjust")))

for (n in 2:nrow(expression_table)){
  x <- cor.test(as.numeric(expression_table[1,]),as.numeric(expression_table[n,]))
  cor_table[n-1,1] <- x$estimate
  cor_table[n-1,2] <- x$p.value
}

cor_table[,3] <- p.adjust(cor_table[,2])

save(cor_table,file = "cor_table.RData")
load("cor_table.RData")

library(ggplot2)

x = which( metadata$exp >= 0  | metadata$Tissue == "Brain - Cerebellar Hemisphere")
df<-as.data.frame(metadata$exp[x])

p <- ggplot(df, aes(x="tissue", y="exp")) + geom_violin()



######## aula 16/05

mRNA <- read.csv("miRWalk_miRNA_Targets.csv")

#How many miRNAs may regulate your gene?
nrow(mRNA) #tem repetidos
length(unique(mRNA$mirnaid)) #correto
#How many were validated?
count_validated = 0
validated_index = c()
for (n in 1:length(mRNA$validated)){
  if (mRNA$validated[n] != ""){
    count_validated =  count_validated + 1
    validated_index = c(validated_index,n)
  }
}
count_validated
validated_index
validated_index_list = mRNA$mirnaid[validated_index] #todos os validados (18)
unique(mRNA$mirnaid[validated_index]) #10

#What it is their binding position?
mRNA$position[validated_index]

# Which miRNAs has the longest consecutive pairing?
#How many nucleotides and binding position?
max_consecutive = max(mRNA$longest_consecutive_pairings[validated_index])
validated_index_list[which(mRNA$longest_consecutive_pairings[validated_index] == max_consecutive)]
mRNA$position[validated_index[which(mRNA$longest_consecutive_pairings[validated_index] == max_consecutive)]]

#Step 2) Identify miRNAs that may regulate indirectly your gene

#TAF1
#SIN3A

TAF1_miRNA <- read.csv("TAF1_miRNA.csv")
SIN3A_miRNA <- read.csv("SIN3A_miRNA.csv")

#TAF1

TAF1_intersection = unique(intersect(TAF1_miRNA$mirnaid,mRNA$mirnaid))
TAF1_intersection

TAF1_miRNA_validated = TAF1_miRNA[TAF1_miRNA$validated != "",]

mRNA_validated = mRNA[validated_index,]

TAF1_intersection_validated = unique(intersect(TAF1_miRNA_validated$mirnaid,mRNA_validated$mirnaid))

SIN3A_intersection = unique(intersect(SIN3A_miRNA$mirnaid,mRNA$mirnaid))

SIN3A_miRNA_validated = SIN3A_miRNA[SIN3A_miRNA$validated != "",]
SIN3A_intersection_validated = unique(intersect(SIN3A_miRNA_validated$mirnaid,mRNA_validated$mirnaid))


eQTL <- read.csv("eQTLs.csv")
table_e <-as.matrix(table(eQTL$Tissue))
NES=table(eQTL$Tissue,eQTL$NES>0)
table_eQTL=cbind(table_e,NES[,1],NES[,2])
colnames(table_eQTL)=c("Nr of eQTLs","Nr of Positive eQTLs","Nr of Negative eQTLs")



sQTL <- read.csv("sQTLs.csv")
table_s <-as.matrix(table(sQTL$Tissue))
NES_s=table(sQTL$Tissue,sQTL$NES>0)
table_sQTL=cbind(table_s,NES_s[,1],NES_s[,2])
colnames(table_sQTL)=c("Nr of eQTLs","Nr of Positive sQTLs","Nr of Negative sQTLs")

library(data.table)
clinivar <- fread("clinvar_result.txt")
clinivar=clinivar[-1,]
clinivar=clinivar[,-16]                         
rs=clinivar[na.omit(match(unique(unlist(eQTL$SNP.Id)),unlist(clinivar$`dbSNP ID`))),]
intersect(unlist(eQTL$SNP.Id),unlist(clinivar$`dbSNP ID`))


rs2 <- clinivar[na.omit(match(unique(unlist(sQTL$SNP.Id)),unlist(clinivar$`dbSNP ID`))),]


gwas=read.table(file = 'gwas-association-downloaded_2022-06-02-ensemblMappedGenes_SETD2.tsv', sep = '\t', header = TRUE)

e_g_ano=gwas[na.omit(match(unique(unlist(eQTL$SNP.Id)),unlist(gwas$SNPS))),]
e = eQTL[na.omit(match(unique(unlist(gwas$SNPS)),unlist(eQTL$SNP.Id))),]

s_g_ano = gwas[na.omit(match(unique(unlist(sQTL$SNP.Id)),unlist(gwas$SNPS))),]
s = sQTL[na.omit(match(unique(unlist(gwas$SNPS)),unlist(sQTL$SNP.Id))),]


m = 1
#x é a exp.table com as colunas com os nomes certos
for (n in 1:ncol(x)){
  str = strsplit(colnames(x)[n]," - ")[[1]][1]
  if (str == "Brain"){
    m = m + 1
    table <- cbind(table,x[,n])
    colnames(table)[m] <- colnames(x)[n]
  }
}
table <- table[,-1]

save(table,file = "brain_table.RData")
cor_table_Brain <- matrix(NA, nrow = 7, ncol = 3, dimnames = list(rownames(expression_table)[-1], c("Correlation","P_value","P_value_adjust")))

for (n in 2:nrow(table)){
  x <- cor.test(as.numeric(table[1,]),as.numeric(table[n,]))
  cor_table_Brain[n-1,1] <- x$estimate
  cor_table_Brain[n-1,2] <- x$p.value
}

cor_table_Brain[,3] <- p.adjust(cor_table_Brain[,2])
save(cor_table_Brain,file = "cor_table_Brain.RData")


library(ggplot2)

ggplot(x, aes(x = `metadata$Tissue`, y = `metadata$exp`, fill = `metadata$Tissue`)) +
  geom_violin(trim = FALSE,draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_boxplot(width = 20) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.text = element_text(angle = 45,hjust = 1),
        axis.line = element_line(colour = "black"))


x <- as.data.frame(metadata$Tissue)
x <- cbind(x,as.numeric(metadata$exp))
x <- cbind(x,as.character(metadata$SEX))
colnames(x) <- c("Tissue","exp","sex")
###PLOT SEX####
ggplot(x, aes(x = Tissue, y = exp, fill = sex)) +
  introdataviz::geom_split_violin(trim = FALSE, alpha = 0.4,width = 2) +
  geom_boxplot(width = 0.05) + 
  scale_fill_manual(values=c("steelblue1", "palevioletred1")) + 
  labs(y = "TPM") + 
  theme_classic() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.text = element_text(angle = 45,hjust = 1),
        axis.line = element_line(colour = "black"),
        plot.margin = margin(1,0.5,0.8,1.2, "cm"),
        axis.title.x = element_blank()
  )

###PLOT TISSUE####
ggplot(x, aes(x = `Tissue`, y = `exp`, fill = `Tissue`)) +
  geom_violin(trim = FALSE, alpha = 0.4,width = 1.5) +
  geom_boxplot(width = 0.05) + 
  scale_fill_manual(values=c("chocolate1", "goldenrod2", "green2", "indianred2", "mistyrose1","firebrick1","indianred4", "yellow2","yellow2","yellow2","yellow2","yellow2","yellow2","yellow2","yellow2","yellow2","yellow2","yellow2","yellow2", "turquoise3", "skyblue", "plum","rosybrown2","thistle","wheat1","wheat2","wheat3","lightyellow4", "wheat4", "tan","lavenderblush", "mediumorchid1", "mediumorchid4", "mediumturquoise", "mediumspringgreen", "palegreen3", "palegreen2", "palegreen3", "lightslateblue", "gold1", "lightpink2","burlywood3", "lightgreen", "gray88", "mediumslateblue", "mediumpurple1", "yellow4", "yellowgreen", "wheat2", "snow3", "seagreen4", "orchid1", "indianred1", "deeppink")) + 
  labs(y = "TPM") + 
  theme_classic() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none", 
        axis.text = element_text(angle = 45,hjust = 1),
        axis.line = element_line(colour = "black"),
        plot.margin = margin(1,0.5,0.8,1.2, "cm"),
        axis.title.x = element_blank()
  )

tiff( "sQTLs_plot.tiff", units= "cm", width= 30, height=20, res=600)

ggplot(table_sQTL, aes(x=as.factor(Tissue), y=value, fill=observation)) +
  geom_bar(stat="identity",width = 0.5) +
  scale_fill_viridis(discrete=TRUE, name="") +
  labs(y = "Number of sQTLs") + 
  theme_classic() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9, 1),
        legend.direction = "vertical",
        legend.title = element_blank(),
        axis.text = element_text(angle = 45,hjust = 1),
        axis.line = element_line(colour = "black"),
        plot.margin = margin(1,0.5,0.8,1.2, "cm"),
        axis.title.x = element_blank()
  )

dev.off()

# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)



table_eQTL <- cbind(row.names(table_eQTL), table_eQTL)
colnames(table_eQTL)=c("Tissue", "Nr of eQTLs","Nr of Positive eQTLs","Nr of Negative eQTLs")
table_eQTL <- table_eQTL[-2]
table_eQTL <- table_eQTL %>% gather(key = "observation", value="value", c(2,3)) 
table_eQTL$value <- as.numeric(table_eQTL$value)

ggplot(table_eQTL, aes(x=as.factor(Tissue), y=value, fill=observation)) +
  geom_bar(stat="identity",width = 0.5) +
  scale_fill_viridis(discrete=TRUE, name="") +
  labs(y = "Number of eQTLs") + 
  theme_classic() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.8, 0.9),
        legend.direction = "vertical",
        legend.title = element_blank(),
        axis.text = element_text(angle = 45,hjust = 1),
        axis.line = element_line(colour = "black"),
        plot.margin = margin(1,0.5,0.8,1.2, "cm"),
        axis.title.x = element_blank()
  )

table_sQTL <- cbind(row.names(table_sQTL), table_sQTL)
colnames(table_sQTL)=c("Tissue", "Nr of eQTLs","Nr of Positive sQTLs","Nr of Negative sQTLs")
table_sQTL <- as.data.frame(table_sQTL)
table_sQTL <- table_sQTL[-2]
table_sQTL$`Nr of Positive sQTLs` <- as.numeric(table_sQTL$`Nr of Positive sQTLs`)
table_sQTL$`Nr of Negative sQTLs` <- as.numeric(table_sQTL$`Nr of Negative sQTLs`)
table_sQTL <- table_sQTL %>% gather(key = "observation", value="value", c(2,3)) 
table_sQTL$value <- as.numeric(table_sQTL$value)

ggplot(table_sQTL, aes(x=as.factor(Tissue), y=value, fill=observation)) +
  geom_bar(stat="identity",width = 0.5) +
  scale_fill_viridis(discrete=TRUE, name="") +
  labs(y = "Number of sQTLs") + 
  theme_classic() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9, 1),
        legend.direction = "vertical",
        legend.title = element_blank(),
        axis.text = element_text(angle = 45,hjust = 1),
        axis.line = element_line(colour = "black"),
        plot.margin = margin(1,0.5,0.8,1.2, "cm"),
        axis.title.x = element_blank()
  )

