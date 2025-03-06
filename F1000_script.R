#R protocol written and implemented by Dr Rhiannon S Morgan, Dr Marie M Phelan and Professor Richard Barrett-Jolley
#Training in the use of R and Rstudio with in-house scripts was provided by the Liverpool Shared Research Facility (LivSRF) Computational Biology Facility (Dr Eva Caamano Gutierrez & Dr Arturas Grauslys) and the Highfield NMR facility (Dr Rudi Grosman). 

#setup working directory - edit as required
setwd('my_location')

#read in publication data first column is unique sample ID second column is group
data1<-read.csv("MTBLS2396.csv",header=T,stringsAsFactors = FALSE)

#additional packages required for scripts - written and executed using R 3.6.1 (if not already installed in R)
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("ellipse")

##in-house scripts
source("NMRMetab_norm_scale.R")
source("NMRMetab_ANOVA.R")

####data normalisation and scaling
data1P= NMRMetab_norm_scale(data1, normalisation = 'PQN')
data1PP= NMRMetab_norm_scale(data1, normalisation = 'PQN', scaling ='Pareto')

####subsetting data
#chorion 1:3 /nochorion 4:6 
data1PPa=((data1PP[1:6,]))
#epp, glass,N2
data1PPb=(data1PP[c(22:24,37:42),])
#dpf
data1PPc=(data1PP[c(7:12,22:24),])
#head,tail whole
data1PPd=(data1PP[c(22:36),])
#head,tail whole (headgroups 113-181 bins only)
data1PPe=(data1PP[c(22:36),c(1,2,115:183)])

##PCAs
#chorion PCA (one example provided - change colours where indicated to swap plots)
source("NMRMetab_PCA_3Cyan.R")
NMRMetab_PCA(data1PPa,drawEllipses=T)
#N2glass PCA
source("NMRMetab_PCA_3Red.R")
NMRMetab_PCA(data1PPb,drawEllipses=T)
#dpf PCA
source("NMRMetab_PCA_3Blue.R")
NMRMetab_PCA(data1PPc,drawEllipses=T)
#head/tail PCA
source("NMRMetab_PCA_3Green.R")
NMRMetab_PCA(data1PPd,drawEllipses=T)
#head/tail PCA (headgroup 113-181 only)
source("NMRMetab_PCA_3Green.R")
NMRMetab_PCA(data1PPe,drawEllipses=T)

##Performing Q-Q plots
#q-q plots for subsetted data on peaks of interest - examples below
data1_bin181=((data1P[,183]))
qqnorm(data1_bin181, pch = 1, frame = FALSE)
qqline(data1_bin181, col = "steelblue", lwd = 2)
ggplot(data1_bin181, aes(sample = data1_bin181[,3])) +
stat_qq() +
stat_qq_line()

data1_bin140=((data1P[,142]))
qqnorm(data1_bin140, pch = 1, frame = FALSE)
qqline(data1_bin140, col = "steelblue", lwd = 2)
ggplot(data1_bin140, aes(sample = data1_bin140[,3])) +
stat_qq() +
stat_qq_line()

#glass, epp, N2, bin5
data1Pb_b=(data1P[c(22:24,37:42),10])
data1Pb1_b=(data1P[c(37:72),7])
qqnorm(data1Pb1_b, pch = 1, frame = FALSE)
qqline(data1Pb1_b, col = "steelblue", lwd = 2)
ggplot(data1Pb1_b, aes(sample = data1Pb1_b[,3])) +
stat_qq() +
stat_qq_line()

##ANOVAs of each sample condition comparisons
NMRMetab_anova(data1PPa)
NMRMetab_anova(data1PPb)
NMRMetab_anova(data1PPc)
NMRMetab_anova(data1PPd)
NMRMetab_anova(data1PPe)

####Cluster Significance for PCA####
##in-house script for scaling and normalisation
source("NMRMetab_norm_scale.R")

#install package (if not already installed in R)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ClusterSignificance")

#read library
library(ggplot2)
library(ClusterSignificance)

#Iterations
iterations=1000

#read in dataset 
data1<-read.csv("MTBLS2396.csv",header=T,stringsAsFactors = FALSE)

#data normalisation and scaling
data1PP= NMRMetab_norm_scale(data1, normalisation = 'PQN', scaling ='Pareto')

##subsetting data
#chorion 1:3 /nochorion 4:6 
data1PPa=((data1PP[1:6,]))

#perform PCA 
dat<-data1PPa[,3:ncol(data1PPa)]
groups<-as.factor(data1PPa[,2])
Labels<-as.character(data1PPa[,1])
pc<-prcomp(dat, scale = T)
pcs=c(1,2)
pcdf<-data.frame(pc1=pc$x[,pcs[1]], pc2=pc$x[,pcs[2]])
pcpMatrix<-as.matrix(pcdf)
colores<-rainbow(length(unique(groups)))
ggplot(data=pcdf, aes(x=pc1, y=pc2)) + 
  geom_point(aes(fill=groups),colour="black",size=5.5,pch=21)+
  scale_fill_manual(values=colores)+scale_y_continuous(labels=function(x)format(x,scientific=T))

# define classes
class1<-rep("Chorion", 3)
class2<-rep("No Chorion",3)
classes=cbind(t(class1),t(class2))[1,]

#Cluster significance analysis 
pe <- permute(
  seed = 35,
  mat = pcpMatrix, 
  iter = iterations, 
  classes = classes, 
  projmethod = "pcp"
)
plot(pe)
ClusterSignificance::pvalue(pe)
conf.int(pe, conf.level = 0.95)

##subsetting data
#epp, glass,N2
data1PPb=(data1PP[c(22:24,37:42),])

#read library
library(ggplot2)
library(ClusterSignificance)

#perform PCA
dat<-data1PPb[,3:ncol(data1PPb)]
groups<-as.factor(data1PPb[,2])
Labels<-as.character(data1PPb[,1])
pc<-prcomp(dat, scale = T)
pcs=c(1,2)
pcdf<-data.frame(pc1=pc$x[,pcs[1]], pc2=pc$x[,pcs[2]])
pcpMatrix<-as.matrix(pcdf)
colores<-rainbow(length(unique(groups)))
ggplot(data=pcdf, aes(x=pc1, y=pc2)) + 
  geom_point(aes(fill=groups),colour="black",size=5.5,pch=21)+
  scale_fill_manual(values=colores)+scale_y_continuous(labels=function(x)format(x,scientific=T))

#define classes
class1<-rep("Epi", 3)
class2<-rep("No Epi",3)
class3<-rep("N2",3)
classes=cbind(t(class1),t(class2),t(class3))[1,]

#perform cluster significance analysis
pe <- permute(
  seed = 35,
  mat = pcpMatrix, 
  iter = iterations, 
  classes = classes, 
  projmethod = "pcp"
)
plot(pe) 
ClusterSignificance::pvalue(pe)
conf.int(pe, conf.level = 0.95)

##subsetting data
#dpf
data1PPc=(data1PP[c(7:12,22:24),])

#read library
library(ggplot2)
library(ClusterSignificance)

#perform PCA
dat<-data1PPc[,3:ncol(data1PPc)]
groups<-as.factor(data1PPc[,2])
Labels<-as.character(data1PPc[,1])
pc<-prcomp(dat, scale = T)
pcs=c(1,2)
pcdf<-data.frame(pc1=pc$x[,pcs[1]], pc2=pc$x[,pcs[2]])
pcpMatrix<-as.matrix(pcdf)
colores<-rainbow(length(unique(groups)))
ggplot(data=pcdf, aes(x=pc1, y=pc2)) + 
  geom_point(aes(fill=groups),colour="black",size=5.5,pch=21)+
  scale_fill_manual(values=colores)+scale_y_continuous(labels=function(x)format(x,scientific=T))

#define classes
class1<-rep("1dpf",3)
class2<-rep("2dpf",3)
class3<-rep("3dpf",3)
classes=cbind(t(class1),t(class2),t(class3))[1,]
View(data1PPc)

#perform cluster significance analysis
pe <- permute(
  seed = 35,
  mat = pcpMatrix, 
  iter = iterations, 
  classes = classes, 
  projmethod = "pcp"
)
plot(pe)
ClusterSignificance::pvalue(pe)
conf.int(pe, conf.level = 0.95)

##subsetting data
#head, tail, whole
data1PPd=(data1PP[c(22:36),])

#read library
library(ggplot2)
library(ClusterSignificance)

#perform PCA
dat<-data1PPd[,3:ncol(data1PPd)]
groups<-as.factor(data1PPd[,2])
Labels<-as.character(data1PPd[,1])
pc<-prcomp(dat, scale = T)
pcs=c(1,2)
pcdf<-data.frame(pc1=pc$x[,pcs[1]], pc2=pc$x[,pcs[2]])
pcpMatrix<-as.matrix(pcdf)
colores<-rainbow(length(unique(groups)))
ggplot(data=pcdf, aes(x=pc1, y=pc2)) + 
  geom_point(aes(fill=groups),colour="black",size=5.5,pch=21)+
  scale_fill_manual(values=colores)+scale_y_continuous(labels=function(x)format(x,scientific=T))

#define classes
class1<-rep("whole", 3)
class2<-rep("head",6)
class3<-rep("tail",6)
classes=cbind(t(class1),t(class2),t(class3))[1,]

#perform cluster significance analysis
pe <- permute(
  seed = 35,
  mat = pcpMatrix, 
  iter = iterations, 
  classes = classes, 
  projmethod = "pcp"
)
plot(pe)
ClusterSignificance::pvalue(pe)
conf.int(pe, conf.level = 0.95)

##subsetting the data
#head,tail whole (headgroups 113-181 bins only)
data1PPe=(data1PP[c(22:36),c(1,2,115:183)])

#read library
library(ggplot2)
library(ClusterSignificance)

#perform PCA
dat<-data1PPe[,3:ncol(data1PPe)]
groups<-as.factor(data1PPe[,2])
Labels<-as.character(data1PPe[,1])
pc<-prcomp(dat, scale = T)
pcs=c(1,2)
pcdf<-data.frame(pc1=pc$x[,pcs[1]], pc2=pc$x[,pcs[2]])
pcpMatrix<-as.matrix(pcdf)
colores<-rainbow(length(unique(groups)))
ggplot(data=pcdf, aes(x=pc1, y=pc2)) + 
  geom_point(aes(fill=groups),colour="black",size=5.5,pch=21)+
  scale_fill_manual(values=colores)+scale_y_continuous(labels=function(x)format(x,scientific=T))

#define classes
class1<-rep("whole", 3)
class2<-rep("head",6)
class3<-rep("tail",6)
classes=cbind(t(class1),t(class2),t(class3))[1,]

#perform cluster significance analysis
pe <- permute(
#  seed = 35,
  mat = pcpMatrix, 
  iter = iterations, 
  classes = classes, 
  projmethod = "pcp"
)
plot(pe)
ClusterSignificance::pvalue(pe)
conf.int(pe, conf.level = 0.95)

