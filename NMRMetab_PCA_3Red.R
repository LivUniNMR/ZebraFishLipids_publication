#--------------------------------------------------------------------------------#
#---------- Principal Component Analysis via SVD for NMR Metabolomics data ------#
#--------------------------------------------------------------------------------#
#
# Dr Eva Caamano Gutierrez & Dr Arturas Grauslys, 2017
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
# https://creativecommons.org/licenses/by-nc-sa/4.0/
# The script's aim is to allow unexperienced R users to perform PCA in NMR metabolomics
# data prepared such as rows contain samples and columns variables
# furthermore the first two columns are data describers: row names and groups
# There are different options that the user must specify including
# (a) the principal components you want to plot (1 and 2 are default)
# (b) whether you want labels on your score plots or loading plots (set useLabelsscores/loadings=T/F)
# (c) scale=T/F if your data is not scaled you should say TRUE, otherwise FALSE
# (d) if you want other name in the legend of your plots you can change it in legendName
# (e) if you want ellipses drawn in your score plot you can set it in drawEllipses (T/F)
# In order to use this script follow this steps:
# (1) Source this file
# (2) Call the NMRMetab_PCA() function with adequate parameters
# e.g. NMRMetab_PCA(myBucketTable,useLabels=T,pcs=c(1,2),scale=T, legendName="Group")
# These options can be changed and/or expanded but any aditions should be 
# documented and dated in this preface.
#--------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#note all the packages below need to be installed prior to the use of this function
library(ggplot2)
library(ggrepel)
library(ellipse)


NMRMetab_PCA<-function(data,useLabelsScores=F,useLabelsLoadings=T,pcs=c(1,2),scale=F, legendName="Group", drawEllipses=T){
  
  path<-makeTSFolder("NMRMetab_PCA")
  dat<-data[,3:ncol(data)]
  group<-as.factor(data[,2])
  Labels<-as.character(data[,1])
  
  
  s<-insidePCA_function(dat,group,useLabelsScores,useLabelsLoadings,labels=Labels,pcs,type='scores',scale=scale,legendName=legendName, drawEllipses=drawEllipses,path=path)
  
  ggsave(plot=s,file=paste(path,"/PCAScores",".pdf",sep=""))
  
  l<-insidePCA_function(dat,group,useLabelsScores,useLabelsLoadings,labels="",pcs,type='loadings',scale=scale,legendName=legendName,path=path)
   ggsave(plot=l,file=paste(path,"/PCALoadings",".pdf",sep=""))
  
  va<- insidePCA_function(dat,group,useLabelsScores=F,useLabelsLoadings=F,pcs,type='varAcc',scale=scale,legendName=legendName,path=path)
  ggsave(plot=va,file=paste(path,"/VarAccounted",".pdf",sep=""))
  
  print("PCA function completed")
}

makeTSFolder = function(prefix){
  ts = format(Sys.time(), "%b_%d_%Y_%X")
  
  ts<-gsub(":","-",ts)
  tsDir = paste(prefix, ts, sep='-')
  if(!file.exists(file.path(getwd(),tsDir))) dir.create(file.path(getwd(),tsDir))
  return(file.path(getwd(),tsDir))
}

insidePCA_function<-function(data, groups, useLabelsScores=F,useLabelsLoadings=T, labels = "", pcs=c(1,2),
                                 type='scores', scale=T, legendName="Treatment", drawEllipses=T,path=""){
  
  # Performs PCA on a given dataset and plots a choice of results: 
  # scores, loadings, cumulative variance accounted for or biplot.
  #
  # INPUTS:
  #
  #  data - data.frame or matrix
  #   - data to analyse with variables in columns and samples in rows
  #  groups - factor
  #   - a grouping factor that determines the colours of the points in scores plots
  #  useLabels (optional) - boolean - default=FALSE
  #   - if TRUE will draw labels next to each point
  #  labels (optional) - default=""
  #   - labels to be drawn next to each point. If useLabels=TRUE and labels is empty will use rownames(data) as labels.
  #  pcs (optional) - a vector of 2 integers - default=c(1,2)
  #   - principal components to be plotted
  #  type (optional) - string - default='scores'
  #   - which plot to produce at the end. Possible: 'scores', 'loadings', 'varAcc'.
  #  scale (optional) - boolean - default=TRUE
  #   - determines if variables are to be scaled (recommended)
  #  legendName (optional) - string - default='Groups'
  #   - sets the name for the legend of colours (set by groups)
  #
  # OUTPUTS:
  # a ggplot object. If not assigned to a variable will show the plot.

#change colour & point shape here - need to have correct number for groups used:
#points 21-15 only may be used with outline/fill options
# NMRMetab_PCA_3red.R -edit and rename as needed
   colores<-c("pink","red","rosybrown")
   pointss<-c(22,21,24)
   outlines<-c("black","black","black") 

#  NMRMetab_PCA_3Green.R
#   colores<-c("springgreen","green3","darkgreen")
#   pointss<-c(25,22,24)
#   outlines<-c("black","black","black") 

#  NMRMetab_PCA_3Pink.R
# colores<-c("violet","violetred3","maroon4")
#   pointss<-c(21,22,23)
#   outlines<-c("black","black","black") 

#  NMRMetab_PCA_3Cyan.R
#   colores<-c("cyan","cyan3","darkcyan")
#   pointss<-c(21,23,24)
#   outlines<-c("black","black","black") 

#  NMRMetab_PCA_3Blue.R
#   colores<-c("steelblue","blue","darkblue")
#   pointss<-c(22,23,25)
#   outlines<-c("black","black","black")  
  
  data<-as.matrix(data)
  data<-data[,!apply(data,2,function(x) all(x==0))]
  data<-data[,!apply(data,2,function(x) any(is.na(x)))]
  
  
  if(scale){
    pc<-prcomp(data, scale = T)
  } else {
    pc <- prcomp(data)
  }
  
  
  
  if(type=='scores'){
    if(useLabelsScores & length(labels) != nrow(data)){
      print("Warning: The labels not given or given incorrectly. Using rownames.")
      labels <- rownames(data)
    }
    
    pcdf<-data.frame(pc1=pc$x[,pcs[1]], pc2=pc$x[,pcs[2]])
    
    if(path==""){
    write.csv(pc$x,paste(makeTSFolder("NMRMetab_PCA"),"/PCAScores",".csv",sep=""))
    }else{
      write.csv(pc$x,paste(path,"/PCAScores",".csv",sep=""))
    }
    
    if(useLabelsScores) pcdf$labels<-labels
    perc_accounted<-summary(pc)$importance[2,pcs]*100


    
    .e <- environment()
    p <- ggplot(data=pcdf, aes(x=pc1, y=pc2), environment=.e) + 
      geom_point(aes(fill=groups,shape=groups), colour="black",size=3) +
#      geom_point(aes(fill=groups),colour="black",size=5.5,pch=21)+
     scale_shape_manual(values=pointss)+
     scale_fill_manual(values=colores)+
     scale_color_manual(values=outlines) 
# scale_fill_manual(values=colores,name=legendName)+scale_y_continuous(labels=function(x)format(x,scientific=T))
    
    if(useLabelsScores)  p <- p + geom_text_repel(aes(label = labels))
    
    if (drawEllipses){
      # Generate an ellipse for each groups of points 
      df_ell = do.call('rbind', lapply(levels(as.factor(groups)), function(grp) {
        ellipse(cor(pcdf[which(groups==grp),'pc1'], pcdf[which(groups==grp),'pc2']), 
                scale=c(sd(pcdf[which(groups==grp),'pc1']),sd(pcdf[which(groups==grp),'pc2'])), 
                centre=c(mean(pcdf[which(groups==grp),'pc1']),mean(pcdf[which(groups==grp),'pc2'])))
        }
      ))
      
      df_ell = as.data.frame(df_ell)
      grp = levels(as.factor(groups))
      grp = do.call('c', lapply(grp, function(x) rep(x, 100)))
      df_ell$group = grp
    }
    
    p <- p+ #guides(color=FALSE)+
      #guides(colour = guide_legend(title = legendName))+
      xlab(paste("PC",pcs[1], " (", round(perc_accounted[1],2), "%)", sep=""))+
      ylab(paste("PC",pcs[2], " (", round(perc_accounted[2],2), "%)", sep=""))+
      #ggtitle("PCA scores plot")+
      theme_bw(base_size=20)+
      theme(legend.position="bottom")
    #theme(legend.position="none")
    if (drawEllipses){
      p <- p+
#CHANGE elipse style here
#        geom_path(data=df_ell, aes(x=x, y=y, colour=group ), size=1, linetype=2)+
        geom_path(data=df_ell, aes(x=x, y=y, colour=group ), size=0.5, linetype=2)+
        scale_colour_manual(values=colores, guide=F)
    }
    
    p
    
  } else if(type=='loadings'){
    
    if(useLabelsLoadings & length(labels) != nrow(pc$rotation)){
      #print("Warning: loadings labels not given or given incorrectly. Using the column names.")
      labels <- colnames(data)
    }
    
    pcdf<-data.frame(load1=pc$rotation[,pcs[1]], load2=pc$rotation[,pcs[2]], var=labels)
    if(path==""){
    write.csv(pc$rotation,paste(makeTSFolder("NMRMetab_PCA"),"/PCALoadings",".csv",sep=""))
    }else{
      write.csv(pc$rotation,paste(path,"/PCALoadings.csv",sep=""))
    }
              
    
    
    .e <- environment()
    
    p <- ggplot(data=pcdf, aes(x=load1, y=load2), environment=.e) + geom_point()+scale_y_continuous(labels=function(x)format(x,scientific=T))
    
    if(useLabelsLoadings) p <- p + geom_text_repel(aes(x=load1,y=load2),label=labels)
    
    
    p <- p+
      xlab(paste("Loadings for PC",pcs[1],sep=""))+
      ylab(paste("Loadings for PC",pcs[2],sep=""))+
      ggtitle("PCA loadings plot")+
      theme_bw(base_size=20)
    p
    
  } else if(type=='varAcc'){
    perc_accounted<-summary(pc)$importance[2,]*100
    
    perc_with_cumsum <- data.frame(pc = as.factor(1:length(perc_accounted)),
                                   perc_acc = perc_accounted,
                                   perc_cumsum = cumsum(perc_accounted))
    p<- ggplot(data = perc_with_cumsum, aes(x=pc, y=perc_cumsum))+
      geom_bar(stat='identity', col='black', fill='white')+
      geom_hline(yintercept = 95, col='red')+
      geom_hline(yintercept = 0, col='black')+
      xlab('PC')+
      ylab('% Variance')+
      ggtitle('% Variance accounted for by principle components')+
      theme_bw()
      
  }else if(type=='biplot'){
    
    p<- ggbiplot(pc,scale=1,obs.scale=T,var.scale=T,groups=groups,ellipse=T,
                 circle=F,labels=groups)+scale_colour_discrete(name="")+
      theme_bw(base_size=20)+theme(legend.position="none")
    
   
 
      }else {
    cat(sprintf("\nError: no type %s", type))
  }
  
  return(p)
}



















